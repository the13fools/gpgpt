#ifndef OPTZOO_H
#define OPTZOO_H

#include <string>


#include <TinyAD/ScalarFunction.hh>
#include "Surface.h"
#include "AppState.h"

#include "UtilsMisc.h"

/**
 * Enumeration for the different field views.
 */

 namespace OptZoo 
 {

    using SF6 = decltype(TinyAD::scalar_function<6>(TinyAD::range(1)));


// The pinned boundary condition.  
// TODO,  add in other boundary conditions like mixed neumann and free for meshing examples.  
void addPinnedBoundaryTerm(SF6& func, const AppState& appState) {

    func.add_elements<4>(TinyAD::range(appState.F.rows()), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
    {

     // Evaluate element using either double or TinyAD::Double
        using T = TINYAD_SCALAR_TYPE(element);



        // Get variable 2D vertex positions
        Eigen::Index f_idx = element.handle;
        Eigen::VectorX<T> s_curr = element.variables(f_idx);
        Eigen::Vector2<T> curr =  s_curr.head(2);
        Eigen::Vector4<T> delta = s_curr.tail(4);

        Eigen::Matrix2<T> currcurr = curr*curr.transpose();
        Eigen::Vector4<T> currcurrt = flatten(currcurr);

        Surface cur_surf = *(appState.cur_surf);

        T w_bound = appState.config->w_bound;

        Eigen::VectorXi bound_face_idx = appState.bound_face_idx;

        // Surface* cur_surf = appState;
          // metadata field
          // Eigen::Vector2<T> metadata = s_curr.segment(2, 2);


        if (bound_face_idx(f_idx) == 1)
        {
            Eigen::Vector2<T> targ = appState.frames_orig.row(f_idx);
            return w_bound*(curr-targ).squaredNorm() + w_bound*delta.squaredNorm();
        }

        if (bound_face_idx(f_idx) == -1)
        {
            T ret = w_bound*delta.squaredNorm();
            for(int i = 0; i < 3; i++)
            {
            
              int neighbor_edge_idx = cur_surf.data().faceNeighbors(f_idx, i);
              if(neighbor_edge_idx > -1)
              {
                Eigen::VectorX<T> s_n = element.variables(neighbor_edge_idx);
                Eigen::Vector2<T> n_i = s_n.head(2);
                // ret = ret + (n_i-curr).squaredNorm() * w_smooth;
                Eigen::Matrix2<T> nini = n_i*n_i.transpose();
                Eigen::Vector4<T> ninit = flatten(nini);
                ret = ret + (ninit-currcurrt).squaredNorm() * w_bound; // * w_smooth * w_attenuate;
                // ret = ret + (n_i*n_i.transpose()-currcurr).norm() * w_smooth;
              }
            }
            return ret;
        }

    });

 

}


void addSmoothnessTerm(SF6& func, const AppState& appState) {

    func.add_elements<4>(TinyAD::range(appState.F.rows()), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
    {
        using T = TINYAD_SCALAR_TYPE(element);
        // Init element state from solver DOFs + appState optimization data.  

         // Get variable 2D vertex positions
        Eigen::Index f_idx = element.handle;
        Eigen::VectorX<T> s_curr = element.variables(f_idx);
        Eigen::Vector2<T> curr =  s_curr.head(2);
        Eigen::Vector4<T> delta = s_curr.tail(4);

        Eigen::Matrix2<T> currcurr = curr*curr.transpose();
        Eigen::Vector4<T> currcurrt = flatten(currcurr);

        Surface cur_surf = *appState.cur_surf;

        double w_bound = appState.config->w_bound;
        double w_smooth_vector = appState.config->w_bound;

        double w_smooth = appState.config->w_bound;
        double w_curl = appState.config->w_curl;
        double w_attenuate = appState.config->w_attenuate;
        // double w_bound = appState.config->w_bound;




        // 

          Eigen::VectorX<T> s_a = element.variables(cur_surf.data().faceNeighbors(f_idx, 0));
          Eigen::VectorX<T> s_b = element.variables(cur_surf.data().faceNeighbors(f_idx, 1));
          Eigen::VectorX<T> s_c = element.variables(cur_surf.data().faceNeighbors(f_idx, 2));



          Eigen::Vector2<T> a = s_a.head(2);
          Eigen::Vector2<T> b = s_b.head(2);
          Eigen::Vector2<T> c = s_c.head(2);

          Eigen::Matrix2<T> aa = a*a.transpose();
          Eigen::Matrix2<T> bb = b*b.transpose();
          Eigen::Matrix2<T> cc = c*c.transpose();


          Eigen::Vector4<T> a_delta = s_a.tail(4);
          Eigen::Vector4<T> b_delta = s_b.tail(4);
          Eigen::Vector4<T> c_delta = s_c.tail(4);

          Eigen::Vector4<T> aat = flatten(aa);
          Eigen::Vector4<T> bbt = flatten(bb);
          Eigen::Vector4<T> cct = flatten(cc);

          aat = aat + a_delta;
          bbt = bbt + b_delta; 
          cct = cct + c_delta;
          currcurrt = currcurrt + delta;



                    // return (T) 0.;
         






          Eigen::Vector2<T> curr_normalized = curr.normalized();
          Eigen::Vector2<T> curr_perp; // = curr_normalized;
          curr_perp(0) = curr_normalized(1);
          curr_perp(1) = -curr_normalized(0);

          // T s_perp_term = pow(a.dot(curr_perp),2) + pow(b.dot(curr_perp),2) + pow(c.dot(curr_perp), 2);

          // T s_perp_term = ((a.dot(curr_perp) + b.dot(curr_perp) + c.dot(curr_perp)) * (curr_perp * curr_perp.transpose())).norm();
          T s_perp_term = ((a.dot(curr_perp) + b.dot(curr_perp) + c.dot(curr_perp)) * (currcurrt)).squaredNorm();

          // T primal_dirichlet_term = (a + b + c - 3*curr).squaredNorm();
          // T dirichlet_term = (aat+bbt+cct-3*currcurrt).squaredNorm();

          T primal_dirichlet_term = (a - curr).squaredNorm() + (b - curr).squaredNorm() + (c - curr).squaredNorm();
          T dirichlet_term = (aat-currcurrt).squaredNorm() + (bbt-currcurrt).squaredNorm() + (cct-currcurrt).squaredNorm();



          // T dirichlet_term = (aa + bb + cc - 3*currcurr).norm();
  // dirichlet_term += 1e-5*abs(dirichlet_term - metadata(0));
          // T delta_rescale = std::max(frames.row(f_idx).squaredNorm(), 1e-8);
          // delta_rescale = (.0001 + 1./delta_rescale);
          T delta_rescale = 1.;
          // std::cout << delta_rescale << std::endl;

        //   appState.smoothness_primal(f_idx) = TinyAD::to_passive(primal_dirichlet_term);
        //   appState.smoothness_sym(f_idx) = TinyAD::to_passive(dirichlet_term);

          // T delta_dirichlet = (a_delta+b_delta+c_delta-3*delta).squaredNorm()*delta_rescale;

          T delta_norm_term = delta_rescale * delta.squaredNorm();// + delta_dirichlet;

    

          // curl_term += 1e-5*abs(curl_term - metadata(1));


          // T atten = 1./(cur_iter*cur_iter + 1);
          // T atten = 1./(cur_iter + 1);

          // T atten = 1.;

          T delta_weight = .1; // std::min(w_curl/100., 1./w_attenuate);

          T ret = delta_norm_term * delta_weight;
          if (w_smooth_vector > 0)
            return w_smooth_vector * primal_dirichlet_term + ret;
          if (w_smooth > 0)
            ret = ret + w_attenuate * w_smooth * dirichlet_term;
     

          return ret;

    });


// (w_smooth * dirichlet_term + 
//                   w_s_perp * s_perp_term) *  + 
//                   w_curl*curl_term  + 
//                  delta_norm_term * delta_weight;

          // return (w_smooth * dirichlet_term + 
          //         w_s_perp * s_perp_term) * atten + 
          //        w_curl*curl_term  + 
          //        delta_norm_term;



}

void addCurlTerm(SF6& func, const AppState& appState) {

    //   Eigen::Vector4d ea = e_projs2.row(cur_surf.data().faceEdges(f_idx, 0));
    //       Eigen::Vector4d eb = e_projs2.row(cur_surf.data().faceEdges(f_idx, 1));
    //       Eigen::Vector4d ec = e_projs2.row(cur_surf.data().faceEdges(f_idx, 2));

    //       // Eigen::Vector4<T> ea = e_projs.at(cur_surf.data().faceEdges(f_idx, 0));
    //       // Eigen::Vector4<T> eb = e_projs.at(cur_surf.data().faceEdges(f_idx, 1));
    //       // Eigen::Vector4<T> ec = e_projs.at(cur_surf.data().faceEdges(f_idx, 2));

    //       // T curl_term = pow(ea.dot(aat + a_delta) - ea.dot(currcurrt + delta),2);
    //       // curl_term +=  pow(eb.dot(bbt + b_delta) - eb.dot(currcurrt + delta),2);
    //       // curl_term +=  pow(ec.dot(cct + c_delta) - ec.dot(currcurrt + delta),2);

    //       T curl_term = pow(ea.dot(aat ) - ea.dot(currcurrt ),2);
    //       curl_term +=  pow(eb.dot(bbt ) - eb.dot(currcurrt ),2);
    //       curl_term +=  pow(ec.dot(cct ) - ec.dot(currcurrt ),2);

    //       appState->curls_sym(f_idx) = TinyAD::to_passive(curl_term);


    //       Eigen::Vector2<T> ea_primal = e_projs_primal.at(cur_surf.data().faceEdges(f_idx, 0));
    //       Eigen::Vector2<T> eb_primal = e_projs_primal.at(cur_surf.data().faceEdges(f_idx, 1));
    //       Eigen::Vector2<T> ec_primal = e_projs_primal.at(cur_surf.data().faceEdges(f_idx, 2));

    //       T curl_term_primal = pow(ea_primal.dot(a) - ea_primal.dot(curr),2);
    //       curl_term_primal +=  pow(eb_primal.dot(b) - eb_primal.dot(curr),2);
    //       curl_term_primal +=  pow(ec_primal.dot(c) - ec_primal.dot(curr),2);

    //       appState->curls_primal(f_idx) = TinyAD::to_passive(curl_term_primal);

    //       T w_curl_new = std::min(1e8, 1./w_attenuate) * w_curl;

    // T ret = T(0);

    //  if (w_s_perp > 0)
    //         ret = ret + w_attenuate * w_s_perp * s_perp_term;
    //       if (w_curl_new > 0)
    //         ret = ret + w_curl_new * curl_term;


}

 } // namespace OptZoo

#endif // OPTZOO_H
