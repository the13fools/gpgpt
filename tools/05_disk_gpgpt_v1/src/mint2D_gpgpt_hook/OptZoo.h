#ifndef OPTZOO_H
#define OPTZOO_H

#include <string>


#include <TinyAD/ScalarFunction.hh>
#include "Surface.h"
#include "AppState.h"

#include "UtilsMisc.h"

#include <igl/find.h>

#include "ElementVarsTinyAD.h"

/**
 * Enumeration for the different field views.
 */

template<int N>
 class OptZoo 
 {
  public:

    using ADFunc = TinyAD::ScalarFunction<N, double, Eigen::Index>; 


// the const test term tried to make all the vectors the all ones or somethgn.
    static void addConstTestTerm(ADFunc& func, AppState& appState) {

        std::cout << "add const obj" << std::endl;

    func.template add_elements<1>(TinyAD::range(appState.F.rows()), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
    {

     // Evaluate element using either double or TinyAD::Double
        using T = TINYAD_SCALAR_TYPE(element);
        using VAR = std::decay_t<decltype(element)>; 

        Eigen::Index f_idx = element.handle;
//         Eigen::VectorXi bound_face_idx = appState.bound_face_idx;

// // Exit early if on a boundary element. 
//         if (bound_face_idx(f_idx) == 1)
//         {
//             return T(0);
//         }

        ProcElement<T,VAR> e(ElementLiftType::primal); 
        e.setSelfData(appState, f_idx, element);

        if (f_idx == 0)
        {

        }


        Eigen::VectorXi bound_face_idx = appState.bound_face_idx;
      
        Eigen::VectorX<T> targ = Eigen::VectorX<T>::Ones(N);

        Eigen::VectorX<T> curr = e.self_data.dofs_curr_elem;
        
        return .001*(curr-targ).squaredNorm(); // + w_bound*delta.squaredNorm();
      



    });

 

}


// the const test term tried to make all the vectors the all ones or somethgn.
static void addUnitNormTerm(ADFunc& func, AppState& appState) {

        std::cout << "add unit norm (ginzburg-landau) obj" << std::endl;

    func.template add_elements<1>(TinyAD::range(appState.F.rows()), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
    {

        using T = TINYAD_SCALAR_TYPE(element);
        using VAR = std::decay_t<decltype(element)>; 

        Eigen::Index f_idx = element.handle;
        Eigen::VectorXi bound_face_idx = appState.bound_face_idx;

// Exit early if on a boundary element. 
        if (bound_face_idx(f_idx) == 1)
        {
            return T(0);
        }

        ProcElement<T,VAR> e(ElementLiftType::primal); 
        e.setSelfData(appState, f_idx, element);

        if (f_idx == 0)
        {

        }
        T ret = T(0);
        T targ = T(1);

        for (int i = 0; i < e.self_data.primal_norms.size(); i++)
        {
            T curr_diff = e.self_data.primal_norms[i] - targ;
            ret = ret + curr_diff*curr_diff;
        }

        return ret; 

    });

 

}



// The pinned boundary condition.  
// TODO,  add in other boundary conditions like mixed neumann and free for meshing examples.  
static void addPinnedBoundaryTerm(ADFunc& func, AppState& appState) {

    std::cout << "add boundary obj: TODO replace with hard constraint" << std::endl;

    func.template add_elements<4>(TinyAD::range(appState.F.rows()), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
    {

     // Evaluate element using either double or TinyAD::Double
        using T = TINYAD_SCALAR_TYPE(element);
        using VAR = std::decay_t<decltype(element)>; // TINYAD_VARIABLES_TYPE(element);

        // Get boundary annotations 
        Eigen::VectorXi bound_face_idx = appState.bound_face_idx;
        Eigen::Index f_idx = element.handle;

        if ((int)f_idx == 0)
        {
            // std::cout << "eval boundary obj" << std::endl;
        }

        // Dirichlet (i.e. pinned) boundary condition
        if (bound_face_idx(f_idx) == 1)
        {
            ProcElement<T,VAR> e(ElementLiftType::primal); 
            e.setSelfData(appState, f_idx, element);

            // Eigen::VectorX<T> curr = e.self_data.primals_rank1[0];
            Eigen::VectorX<T> curr = e.self_data.primals;


            Eigen::VectorX<T> targ = appState.frames_orig.row(f_idx);
            return e.w_bound*(curr-targ).squaredNorm(); // + e.w_bound*e.self_data.delta.squaredNorm();
        }


        // // Free boundary condition
        // TODO: Move this into the L2 and L4 smoothness terms/ make handling more robust.  
        // if (bound_face_idx(f_idx) == -1)
        // {
        //     T ret = e.w_bound*e.self_data.delta.squaredNorm();
        //     for(int i = 0; i < 3; i++)
        //     {
        //       int neighbor_edge_idx = e.cur_surf->data().faceNeighbors(f_idx, i);
        //       if(neighbor_edge_idx > -1)
        //       {
        //         Eigen::VectorX<T> s_n = element.variables(neighbor_edge_idx);
        //         Eigen::Vector2<T> n_i = s_n.head(2);
        //         // ret = ret + (n_i-curr).squaredNorm() * w_smooth;
        //         Eigen::Matrix2<T> nini = n_i*n_i.transpose();
        //         Eigen::Vector4<T> ninit = flatten(nini);
        //         ret = ret + (ninit-e.self_data.currcurrt).squaredNorm() * e.w_bound; // * w_smooth * w_attenuate;
        //         // ret = ret + (n_i*n_i.transpose()-currcurr).norm() * w_smooth;
        //       }
        //     }
        //     return ret;
        // }

        // Not a boundary element
        return T(0);

    });

 

}


// Deps primal vars 
static void addSmoothness_L4_Term(ADFunc& func, AppState& appState) {

    std::cout << "add L4 smoonthess obj" << std::endl;

    func.template add_elements<4>(TinyAD::range(appState.F.rows()), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
    {

     // Evaluate element using either double or TinyAD::Double
        using T = TINYAD_SCALAR_TYPE(element);
        using VAR = std::decay_t<decltype(element)>; 

        Eigen::Index f_idx = element.handle;
        Eigen::VectorXi bound_face_idx = appState.bound_face_idx;

// Exit early if on a boundary element. 
        if (bound_face_idx(f_idx) == 1)
        {
            return T(0);
        }

        ProcElement<T,VAR> e(ElementLiftType::L4_krushkal); 
        // e.setElementVars(appState, f_idx, s_curr);
        e.setSelfData(appState, f_idx, element);
        e.setNeighborData(appState, f_idx, element);


///////////////////
//// Initlaize elementwise objectives 
///////////////////

          T dirichlet_term = T(0);
          for (int i = 0; i < e.num_neighbors; i++)
          {
            dirichlet_term += (e.neighbor_data.at(i).L_4_krushkal - e.self_data.L_4_krushkal).squaredNorm(); 
          }

          appState.os->smoothness_L4(f_idx) = TinyAD::to_passive(dirichlet_term);

          T ret = T(0);//  delta_norm_term * delta_weight;
          ret = ret + e.w_attenuate * e.w_smooth * dirichlet_term;

          return ret;
    });

}




// Deps primal vars 
static void addSmoothness_L2_Term(ADFunc& func, AppState& appState) {

    std::cout << "add L2 smoonthess obj" << std::endl;

    func.template add_elements<4>(TinyAD::range(appState.F.rows()), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
    {

     // Evaluate element using either double or TinyAD::Double
        using T = TINYAD_SCALAR_TYPE(element);
        using VAR = std::decay_t<decltype(element)>; 



        // Get variable 2D vertex positions
        Eigen::Index f_idx = element.handle;
        // Eigen::VectorX<T> s_curr = element.variables(f_idx);

        ProcElement<T,VAR> e(ElementLiftType::L2_krushkal); 
        // e.setElementVars(appState, f_idx, s_curr);
        e.setSelfData(appState, f_idx, element);



        Eigen::VectorXi bound_face_idx = appState.bound_face_idx;

        if ((int)f_idx == 0)
        {
            // std::cout << "eval smoothness obj" << std::endl;
            // std::cout << w_bound << " " << w_smooth_vector << " " << w_smooth << " " << w_curl << " " << w_attenuate << std::endl;
        }


// A bit hacky but exit early if on a boundary element.  Should really do it same as in matlab mint and make boundary elements distict from the mesh. 
        if (bound_face_idx(f_idx) == 1)
        {
            return T(0);
        }

        e.setNeighborData(appState, f_idx, element);


///////////////////
//// Initlaize elementwise objectives 
///////////////////

                  
          // T primal_biharmonic_term = (a + b + c - 3*curr).squaredNorm();
          // T biharmonic_term = (aat+bbt+cct-3*currcurrt).squaredNorm();

          T primal_dirichlet_term = T(0);
          T dirichlet_term = T(0);
          for (int i = 0; i < e.num_neighbors; i++)
          {
            
            primal_dirichlet_term += (e.neighbor_data.at(i).primals - e.self_data.primals).squaredNorm();
            dirichlet_term += (e.neighbor_data.at(i).L_2_krushkal - e.self_data.L_2_krushkal).squaredNorm(); 
          }

        //   appState.os->smoothness_primal(f_idx) = TinyAD::to_passive(primal_dirichlet_term);
          appState.os->smoothness_L2(f_idx) = TinyAD::to_passive(dirichlet_term);

          // T delta_dirichlet = (a_delta+b_delta+c_delta-3*delta).squaredNorm()*delta_rescale;

          // T delta_norm_term = delta_rescale * e.self_data.delta.squaredNorm();// + delta_dirichlet;

          T delta_weight = 1.; // std::min(w_curl/100., 1./w_attenuate);

          T ret = T(0);//  delta_norm_term * delta_weight;
          if (e.w_smooth > 0)
            ret = ret + e.w_attenuate * e.w_smooth * dirichlet_term;
     

          return ret;

    });

}


// deps primal vars 
static void addCurlTerm_L2(ADFunc& func, AppState& appState) {

    std::cout << "add L2 curl obj" << std::endl;

    func.template add_elements<4>(TinyAD::range(appState.F.rows()), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
    {

     // Evaluate element using either double or TinyAD::Double
        using T = TINYAD_SCALAR_TYPE(element);
        using VAR = std::decay_t<decltype(element)>; 



        // Get variable 2D vertex positions
        Eigen::Index f_idx = element.handle;

        ProcElement<T,VAR> e(ElementLiftType::L2_tensor); 
        e.setSelfData(appState, f_idx, element);



        Eigen::VectorXi bound_face_idx = appState.bound_face_idx;

        if ((int)f_idx == 0)
        {
            // std::cout << "eval smoothness obj" << std::endl;
            // std::cout << w_bound << " " << w_smooth_vector << " " << w_smooth << " " << w_curl << " " << w_attenuate << std::endl;
        }


// A bit hacky but exit early if on a boundary element.  Should really do it same as in matlab mint and make boundary elements distict from the mesh. 
        // Eigen::VectorXi bound_face_idx = appState.bound_face_idx;
        if (bound_face_idx(f_idx) == 1)
        {
            return T(0);
        }

        e.setNeighborData(appState, f_idx, element);



///////////////////
//// Curl Term 
///////////////////



          Eigen::Vector4d ea = appState.C_sym_2.row(appState.cur_surf->data().faceEdges(f_idx, 0));
          Eigen::Vector4d eb = appState.C_sym_2.row(appState.cur_surf->data().faceEdges(f_idx, 1));
          Eigen::Vector4d ec = appState.C_sym_2.row(appState.cur_surf->data().faceEdges(f_idx, 2));

          T curl_term = pow(ea.dot(e.neighbor_data.at(0).L_2_primals ) - ea.dot(e.self_data.L_2_primals),2);
          curl_term +=  pow(eb.dot(e.neighbor_data.at(1).L_2_primals ) - eb.dot(e.self_data.L_2_primals),2);
          curl_term +=  pow(ec.dot(e.neighbor_data.at(2).L_2_primals ) - ec.dot(e.self_data.L_2_primals),2);

          appState.os->curl_L2(f_idx) = TinyAD::to_passive(curl_term);

          T w_curl_new = e.w_curl; // std::min(1e8, 1./e.w_attenuate) * e.w_curl;

          

          T ret = T(0);


          if (w_curl_new > 0)
          {
            // std::cout << w_curl_new;
            ret = ret + w_curl_new * curl_term;
            // std::cout << curl_term << " ";
          }
          

          return ret;

    });
}

    // deps primal vars 
static void addCurlTerm_L4(ADFunc& func, AppState& appState) {

    std::cout << "add L4 curl obj" << std::endl;

    func.template add_elements<4>(TinyAD::range(appState.F.rows()), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
    {

     // Evaluate element using either double or TinyAD::Double
        using T = TINYAD_SCALAR_TYPE(element);
        using VAR = std::decay_t<decltype(element)>; 



        // Get variable 2D vertex positions
        Eigen::Index f_idx = element.handle;
        // Eigen::VectorX<T> s_curr = element.variables(f_idx);

        ProcElement<T,VAR> e(ElementLiftType::L4_tensor); 
        // e.setElementVars(appState, f_idx, s_curr);
        e.setSelfData(appState, f_idx, element);



        Eigen::VectorXi bound_face_idx = appState.bound_face_idx;

        if ((int)f_idx == 0)
        {
            // std::cout << "eval smoothness obj" << std::endl;
            // std::cout << w_bound << " " << w_smooth_vector << " " << w_smooth << " " << w_curl << " " << w_attenuate << std::endl;
        }


// A bit hacky but exit early if on a boundary element.  Should really do it same as in matlab mint and make boundary elements distict from the mesh. 
        // Eigen::VectorXi bound_face_idx = appState.bound_face_idx;
        if (bound_face_idx(f_idx) == 1)
        {
            return T(0);
        }

        e.setNeighborData(appState, f_idx, element);



///////////////////
//// Curl Term 
///////////////////



          Eigen::VectorXd ea = appState.C_sym_4.row(appState.cur_surf->data().faceEdges(f_idx, 0));
          Eigen::VectorXd eb = appState.C_sym_4.row(appState.cur_surf->data().faceEdges(f_idx, 1));
          Eigen::VectorXd ec = appState.C_sym_4.row(appState.cur_surf->data().faceEdges(f_idx, 2));

          T curl_term = pow(ea.dot(e.neighbor_data.at(0).L_4_primals ) - ea.dot(e.self_data.L_4_primals),2);
          curl_term +=  pow(eb.dot(e.neighbor_data.at(1).L_4_primals ) - eb.dot(e.self_data.L_4_primals),2);
          curl_term +=  pow(ec.dot(e.neighbor_data.at(2).L_4_primals ) - ec.dot(e.self_data.L_4_primals),2);

        //   curl_term = curl_term / e.self_data.frame_norm_euclidian;
        // double norm_passive = TinyAD::to_passive(e.self_data.frame_norm_euclidian);

        //   curl_term = curl_term / norm_passive;


          appState.os->curl_L4(f_idx) = TinyAD::to_passive(curl_term);

          T w_curl_new = e.w_curl; // std::min(1e8, 1./e.w_attenuate) * e.w_curl;

          

          T ret = T(0);

          if (w_curl_new > 0)
          {
            // std::cout << w_curl_new;
            ret = ret + w_curl_new * curl_term;
            // std::cout << curl_term << " ";
          }
          

          return ret;


    });
}



 }; // namespace OptZoo

#endif // OPTZOO_H



