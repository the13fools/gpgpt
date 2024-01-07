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
    //decltype(TinyAD::scalar_function<N>(TinyAD::range(1)));


// the const test term tried to make all the vectors the all ones or somethgn.
    static void addConstTestTerm(ADFunc& func, AppState& appState) {

        std::cout << "add const obj" << std::endl;

        func.template add_elements<1>(TinyAD::range(appState.nelem), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
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
            // e.setElementVars(appState, f_idx, s_curr);
            e.setSelfData(appState, f_idx, element);
            // e.setNeighborData(appState, f_idx, element);

            if (f_idx == 0)
            {

            }


            Eigen::VectorXi bound_face_idx = appState.bound_face_idx;
        
            Eigen::VectorX<T> targ = Eigen::VectorX<T>::Ones(N);

            Eigen::VectorX<T> curr = e.self_data.dofs_curr_elem;
            
            return 1000.*(curr-targ).squaredNorm(); // + w_bound*delta.squaredNorm();
        



        });

 

    }   



// the const test term tried to make all the vectors the all ones or somethgn.
    static void addWeakOrthogonalityTerm(ADFunc& func, AppState& appState) {

        std::cout << "add const obj" << std::endl;

        func.template add_elements<1>(TinyAD::range(appState.nelem), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
        {

        // Evaluate element using either double or TinyAD::Double
            using T = TINYAD_SCALAR_TYPE(element);
            using VAR = std::decay_t<decltype(element)>; 

            Eigen::Index f_idx = element.handle;
    //         Eigen::VectorXi bound_face_idx = appState.bound_face_idx;

    // // Exit early if on a boundary element. 
            if ( f_idx < appState.cur_tet_mesh->nTets()  ) 
            {
                return T(0);
            }

            ProcElement<T,VAR> e(ElementLiftType::primal); 
            // e.setElementVars(appState, f_idx, s_curr);
            e.setSelfData(appState, f_idx, element);
            // e.setNeighborData(appState, f_idx, element);


            T ret = T(0);
            int nvecs = e.self_data.primals_rank1.size();
            for (int i = 0; i < nvecs; i++)
            {
                for (int j = i+1; j < nvecs; j++)
                {
                    T dprod = e.self_data.primals_rank1[i].transpose() * e.self_data.primals_rank1[j];
                    ret = ret + dprod*dprod;
                }
            }

            return ret * 1e-1 * e.w_attenuate;

        });

 

    }      


// the const test term tried to make all the vectors the all ones or somethgn.
static void addUnitNormTerm(ADFunc& func, AppState& appState) {

        std::cout << "add unit norm (ginzburg-landau) obj" << std::endl;

    func.template add_elements<1>(TinyAD::range(appState.nelem), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
    {

        using T = TINYAD_SCALAR_TYPE(element);
        using VAR = std::decay_t<decltype(element)>; 

        Eigen::Index f_idx = element.handle;
        Eigen::VectorXi bound_face_idx = appState.bound_face_idx;

// // Exit early if on a boundary element. 
//         if (bound_face_idx(f_idx) == 1)
//         {
//             return T(0);
//         }

        if ( f_idx < appState.cur_tet_mesh->nTets()  ) 
        {
            return T(0);
        }

        ProcElement<T,VAR> e(ElementLiftType::primal); 
        // e.setElementVars(appState, f_idx, s_curr);
        e.setSelfData(appState, f_idx, element);
        // e.setNeighborData(appState, f_idx, element);

        if (f_idx == 0)
        {

        }
        T ret = T(0);
        T targ = T(1);

        for (int i = 0; i < e.self_data.primal_norms.size(); i++)
        {
            T curr_diff = e.self_data.primal_norms[i] - targ;

            if (i == 0 )
            {
                curr_diff *= 10;
            }
            else 
            {
                // curr_diff *= 1e-1;
            }



            ret = ret + curr_diff*curr_diff;
            
        }

        return ret * 1e-3; 

    });

 

}


// the const test term tried to make all the vectors the all ones or somethgn.
static void addNormalBoundaryTerm(ADFunc& func, AppState& appState) {

        std::cout << "add unit norm (ginzburg-landau) obj" << std::endl;

    func.template add_elements<1>(TinyAD::range(appState.nelem), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
    {

        using T = TINYAD_SCALAR_TYPE(element);
        using VAR = std::decay_t<decltype(element)>; 

        Eigen::Index f_idx = element.handle;
        Eigen::VectorXi bound_face_idx = appState.bound_face_idx;

// // Exit early if on a boundary element. 
//         if (bound_face_idx(f_idx) == 1)
//         {
//             return T(0);
//         }

        if ( f_idx < appState.cur_tet_mesh->nTets()  ) 
        {
            return T(0);
        }

        ProcElement<T,VAR> e(ElementLiftType::primal); 
        // e.setElementVars(appState, f_idx, s_curr);
        e.setSelfData(appState, f_idx, element);
        // e.setNeighborData(appState, f_idx, element);


        int boundaryFace = appState.cur_tet_mesh->boundaryFace(f_idx - appState.cur_tet_mesh->nTets() );
        // Eigen::Vector3d normal = appState->cur_tet_mesh->faceNormal(boundaryFace);
        Eigen::Vector3d a = appState.V.row(appState.cur_tet_mesh->faceVertex(boundaryFace, 0));
        Eigen::Vector3d b = appState.V.row(appState.cur_tet_mesh->faceVertex(boundaryFace, 1));
        Eigen::Vector3d c = appState.V.row(appState.cur_tet_mesh->faceVertex(boundaryFace, 2));

        Eigen::Vector3<T> b1 = (a-b).normalized();
        Eigen::Vector3<T> b2 = (a-c).normalized();
        Eigen::VectorX<T> n = b1.cross(b2).normalized();
        Eigen::VectorX<T> curr = e.self_data.primals_rank1[0];

        T orth1 = b1.dot( curr );
        T orth2 = b2.dot( curr );
        T unit = 1 - curr.dot(n);
        // T unit = 1 - curr.dot(curr);


        return (orth1*orth1 + orth2*orth2 ) * e.w_bound * e.w_bound * e.w_bound * e.w_bound + unit*unit * e.w_bound;
        // return unit*unit * e.w_bound;


        // if (f_idx == 0)
        // {

        // }
        // T ret = T(0);
        // T targ = T(1);

        // for (int i = 0; i < e.self_data.primal_norms.size(); i++)
        // {
        //     T curr_diff = e.self_data.primal_norms[0] - targ;
        //     ret = ret + curr_diff*curr_diff;
        // }

        // return ret; 

    });

 

}



// The pinned boundary condition.  
// TODO,  add in other boundary conditions like mixed neumann and free for meshing examples.  
static void addPinnedBoundaryTerm(ADFunc& func, AppState& appState) {

    std::cout << "add boundary obj: TODO replace with hard constraint" << std::endl;

    func.template add_elements<4>(TinyAD::range(appState.nelem), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
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
static void addSmoothness_L2_Term(ADFunc& func, AppState& appState) {

    std::cout << "add L2 smoonthess obj" << std::endl;

    func.template add_elements<5>(TinyAD::range(appState.nelem), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
    {

     // Evaluate element using either double or TinyAD::Double
        using T = TINYAD_SCALAR_TYPE(element);
        using VAR = std::decay_t<decltype(element)>; 

        // Get variable 2D vertex positions
        Eigen::Index f_idx = element.handle;

        // ProcElement<T,VAR> e(ElementLiftType::L2_krushkal); 
        // e.setSelfData(appState, f_idx, element);

        ProcElement<T,VAR> e_sym(ElementLiftType::L2_sym_krushkal); 
        e_sym.setSelfData(appState, f_idx, element);


        Eigen::VectorXi bound_face_idx = appState.bound_face_idx;

        if ((int)f_idx == 0)
        {
            // std::cout << "eval smoothness obj" << std::endl;
            // std::cout << w_bound << " " << w_smooth_vector << " " << w_smooth << " " << w_curl << " " << w_attenuate << std::endl;
        }


// A bit hacky but exit early if on a boundary element.  Should really do it same as in matlab mint and make boundary elements distict from the mesh. 
        // Eigen::VectorXi bound_face_idx = appState.bound_face_idx;
        // if (bound_face_idx(f_idx) == 1)
        // {
        //     return T(0);
        // }
        // e.setNeighborData(appState, f_idx, element);

        e_sym.setNeighborData(appState, f_idx, element);

///////////////////
//// Initlaize elementwise objectives 
///////////////////

        //   T primal_dirichlet_term = T(0);
        //   T dirichlet_term = T(0);
        //   for (int i = 0; i < e.num_neighbors; i++)
        //   {
        //     primal_dirichlet_term += (e.neighbor_data.at(i).primals - e.self_data.primals).squaredNorm();
        //     dirichlet_term += (e.neighbor_data.at(i).L_2_krushkal - e.self_data.L_2_krushkal).squaredNorm(); 
        //   }

          T sym_dirichlet_term = T(0);

          for (int i = 0; i < e_sym.num_neighbors; i++)
          {
            Eigen::VectorX<T> curr_diff = e_sym.neighbor_data.at(i).L_2_sym_krushkal - e_sym.self_data.L_2_sym_krushkal;
            sym_dirichlet_term += curr_diff.transpose() * appState.L2_sym_tensor_weights * curr_diff;   
          }

         // This is a test to make sure that the symmetric weights give the same answer as the full tensor 
        //   appState.os->smoothness_L2(f_idx) = TinyAD::to_passive(dirichlet_term - sym_dirichlet_term);

        if ( f_idx < appState.cur_tet_mesh->nTets() )
            appState.os->smoothness_L2(f_idx) = TinyAD::to_passive(sym_dirichlet_term);


          T delta_weight = 1.; // std::min(w_curl/100., 1./w_attenuate);

          T ret = T(0);//  delta_norm_term * delta_weight;

          if (e_sym.w_smooth > 0)
            ret = ret + e_sym.w_smooth * sym_dirichlet_term;


          return ret;

    });

}


// Deps primal vars 
static void addSmoothness_L4_Term(ADFunc& func, AppState& appState) {

    std::cout << "add L4 smoonthess obj" << std::endl;

    func.template add_elements<5>(TinyAD::range(appState.nelem), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
    {

     // Evaluate element using either double or TinyAD::Double
        using T = TINYAD_SCALAR_TYPE(element);
        using VAR = std::decay_t<decltype(element)>; 

        // Get variable 2D vertex positions
        Eigen::Index f_idx = element.handle;

        // ProcElement<T,VAR> e(ElementLiftType::L4_krushkal); 
        // e.setSelfData(appState, f_idx, element);

        ProcElement<T,VAR> e_sym(ElementLiftType::L4_sym_krushkal); 
        e_sym.setSelfData(appState, f_idx, element);


        Eigen::VectorXi bound_face_idx = appState.bound_face_idx;

        if ((int)f_idx == 0)
        {
            // std::cout << "eval smoothness obj" << std::endl;
            // std::cout << w_bound << " " << w_smooth_vector << " " << w_smooth << " " << w_curl << " " << w_attenuate << std::endl;
        }


// A bit hacky but exit early if on a boundary element.  Should really do it same as in matlab mint and make boundary elements distict from the mesh. 
        // Eigen::VectorXi bound_face_idx = appState.bound_face_idx;
        // if (bound_face_idx(f_idx) == 1)
        // {
        //     return T(0);
        // }
        // e.setNeighborData(appState, f_idx, element);

        e_sym.setNeighborData(appState, f_idx, element);

///////////////////
//// Initlaize elementwise objectives 
///////////////////

        //   T primal_dirichlet_term = T(0);
        //   T dirichlet_term = T(0);
        //   for (int i = 0; i < e.num_neighbors; i++)
        //   {
        //     // primal_dirichlet_term += (e.neighbor_data.at(i).primals - e.self_data.primals).squaredNorm();
        //     dirichlet_term += (e.neighbor_data.at(i).L_4_krushkal - e.self_data.L_4_krushkal).squaredNorm(); 
        //   }

          T sym_dirichlet_term = T(0);

          for (int i = 0; i < e_sym.num_neighbors; i++)
          {
            Eigen::VectorX<T> curr_diff = e_sym.neighbor_data.at(i).L_4_sym_krushkal - e_sym.self_data.L_4_sym_krushkal;
            sym_dirichlet_term += curr_diff.transpose() * appState.L4_sym_tensor_weights * curr_diff;   
          }

         // This is a test to make sure that the symmetric weights give the same answer as the full tensor 
        //   appState.os->smoothness_L4(f_idx) = TinyAD::to_passive(dirichlet_term - sym_dirichlet_term);
          if ( f_idx < appState.cur_tet_mesh->nTets() )
            appState.os->smoothness_L4(f_idx) = TinyAD::to_passive(sym_dirichlet_term);


          T delta_weight = 1.; // std::min(w_curl/100., 1./w_attenuate);

          T ret = T(0);//  delta_norm_term * delta_weight;

          if (e_sym.w_smooth > 0)
            ret = ret + e_sym.w_smooth * sym_dirichlet_term;


          return ret;

    });

}



// k is the order of the symmetric curl term. 
static void addCurlTerms(ADFunc& func, AppState& appState) {

    int num_curls = appState.curl_orders.size();
    appState.os->curls_Lks.resize(num_curls);
    for(int term = 0; term < num_curls; term++)
    {
        int k = appState.curl_orders[term];
        std::cout << "add L" << k << " curl obj" << std::endl;
        // appState.os->curls_Lks[term].resize(appState.T.rows());
        appState.os->curls_Lks[term] = Eigen::VectorXd::Zero(appState.T.rows());
    }



    func.template add_elements<5>(TinyAD::range(appState.nelem), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
    {

        int num_curls = appState.curl_orders.size();
        
        // Evaluate element using either double or TinyAD::Double
        using T = TINYAD_SCALAR_TYPE(element);
        using VAR = std::decay_t<decltype(element)>; 

        T ret = T(0);

        if ((appState.config->w_curl) < 1e-9)
            return ret;
        
        for(int term = 0; term < num_curls; term++)
        {
            int k = appState.curl_orders[term];
            Eigen::SparseMatrix<double> weights; 
            switch(k)
            {
                case 2:
                    weights = appState.L2_curl_tensor_weights;
                    break;
                case 4:
                    weights = appState.L4_curl_tensor_weights;
                    break;
                case 6:
                    weights = appState.L6_curl_tensor_weights;
                    break;
                default:
                    std::cout << "curl order not supported" << std::endl;
                    exit(1);
            }

            // Get variable 2D vertex positions
            Eigen::Index f_idx = element.handle;
            // Eigen::VectorX<T> s_curr = element.variables(f_idx);

            ProcElement<T,VAR> e(ElementLiftType::Lk_edge_contract); 
            // e.setElementVars(appState, f_idx, s_curr);
            e.setSelfData(appState, f_idx, element);
            e.edge_contract_order = k; // set the order of the curl term.
            // std::cout << "set edge contract order to " << k << std::endl;



            Eigen::VectorXi bound_face_idx = appState.bound_face_idx;

            if ((int)f_idx == 0)
            {
                // std::cout << "eval smoothness obj" << std::endl;
                // std::cout << w_bound << " " << w_smooth_vector << " " << w_smooth << " " << w_curl << " " << w_attenuate << std::endl;
            }


            
            if ( f_idx >= appState.cur_tet_mesh->nTets() )
            {
                return T(0);
            }

            


    // A bit hacky but exit early if on a boundary element.  Should really do it same as in matlab mint and make boundary elements distict from the mesh. 
            // Eigen::VectorXi bound_face_idx = appState.bound_face_idx;
            // if (bound_face_idx(f_idx) == 1)
            // {
            //     return T(0);
            // }

            e.setNeighborData(appState, f_idx, element);



///////////////////
//// Curl Term 
///////////////////

            T w_curl_new = e.w_curl; // std::min(1e8, 1./e.w_attenuate) * e.w_curl;

            

            int num_neighbors = e.num_neighbors;
            T curl_term = T(0);
            for (int i = 0; i < num_neighbors; i++)
            {
                // if enforce boundary curl, then should multiply the weight by 2, but 
                // this doens't do much, maybe useful for improving boundary coupling?  Can try on and off 
                T boundary_multiplier = 1;
                int n_idx = e.neighbor_data.at(i).curr_idx;
                if (n_idx >= appState.cur_tet_mesh->nTets())
                    boundary_multiplier = 2; // could also be zero? 
                curl_term += boundary_multiplier * e.neighbor_data.at(i).Lk_edge_contract_diff.transpose() * weights * e.neighbor_data.at(i).Lk_edge_contract_diff;
            }

            if ( f_idx < appState.cur_tet_mesh->nTets() )
                appState.os->curls_Lks[term](f_idx) = TinyAD::to_passive(curl_term);

            // hacky, should ideally have a curl_Lk for each k.
            // appState.os->curl_L2(f_idx) = TinyAD::to_passive(curl_term);


            if (w_curl_new > 0)
            {
                ret = ret + 1./e.w_attenuate * w_curl_new * curl_term;
            }
        }
          
          return ret;

    });
}


// deps primal vars 
static void addCurlTerm_L2(ADFunc& func, AppState& appState) {

    std::cout << "add L2 curl obj" << std::endl;

    func.template add_elements<5>(TinyAD::range(appState.nelem), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
    {

     // Evaluate element using either double or TinyAD::Double
        using T = TINYAD_SCALAR_TYPE(element);
        using VAR = std::decay_t<decltype(element)>; 



        // Get variable 2D vertex positions
        Eigen::Index f_idx = element.handle;
        // Eigen::VectorX<T> s_curr = element.variables(f_idx);

        ProcElement<T,VAR> e(ElementLiftType::L2_facets); 
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
        // if (bound_face_idx(f_idx) == 1)
        // {
        //     return T(0);
        // }

        e.setNeighborData(appState, f_idx, element);



///////////////////
//// Curl Term 
///////////////////





          T w_curl_new = e.w_curl; // std::min(1e8, 1./e.w_attenuate) * e.w_curl;

          int num_neighbors = e.num_neighbors;
          T curl_term = T(0);
          for (int i = 0; i < num_neighbors; i++)
          {
                curl_term += (e.neighbor_data.at(i).L_2_facet_diff).squaredNorm();
          }

          appState.os->curl_L2(f_idx) = TinyAD::to_passive(curl_term);


          

          T ret = T(0);

          if (w_curl_new > 0)
          {
            ret = ret + 1./e.w_attenuate * w_curl_new * curl_term;
          }
          
          return ret;

    });
}

    // deps primal vars 
static void addCurlTerm_L4(ADFunc& func, AppState& appState) {

    std::cout << "add L4 curl obj" << std::endl;

    func.template add_elements<5>(TinyAD::range(appState.nelem), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
    {

     // Evaluate element using either double or TinyAD::Double
        using T = TINYAD_SCALAR_TYPE(element);
        using VAR = std::decay_t<decltype(element)>; 



        // Get variable 2D vertex positions
        Eigen::Index f_idx = element.handle;
        // Eigen::VectorX<T> s_curr = element.variables(f_idx);

        ProcElement<T,VAR> e(ElementLiftType::L4_facets); 
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
        // if (bound_face_idx(f_idx) == 1)
        // {
        //     return T(0);
            
        // }

        e.setNeighborData(appState, f_idx, element);



///////////////////
//// Curl Term 
///////////////////

          T w_curl_new = e.w_curl; // std::min(1e8, 1./e.w_attenuate) * e.w_curl;

          int num_neighbors = e.num_neighbors;
          T curl_term = T(0);
          for (int i = 0; i < num_neighbors; i++)
          {
                curl_term += (e.neighbor_data.at(i).L_4_facet_diff).squaredNorm();
          }

          appState.os->curl_L4(f_idx) = TinyAD::to_passive(curl_term);


          T ret = T(0);

          if (w_curl_new > 0)
          {
            ret = ret + 1./e.w_attenuate * w_curl_new * curl_term;
          }
          
          return ret;


    });
}



 }; // namespace OptZoo

#endif // OPTZOO_H



