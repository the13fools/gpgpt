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
 * Enumeration for the different field views. Here we assume that the vector fields are stored barycentrically
 */

template<int N>
class OptZoo
{
public:

    using ADFunc = TinyAD::ScalarFunction<N, double, Eigen::Index>;

    // The pinned boundary condition.  
    // TODO,  add in other boundary conditions like mixed neumann and free for meshing examples.  
    static void addPinnedBoundaryTerm(ADFunc& func, AppState& appState) {

        std::cout << "add boundary obj: TODO replace with hard constraint" << std::endl;

        func.template add_elements<4>(TinyAD::range(appState.F.rows()), [&](auto& element)->TINYAD_SCALAR_TYPE(element)
        {
            // Evaluate element using either double or TinyAD::Double
            using T = TINYAD_SCALAR_TYPE(element);
            using VAR = std::decay_t<decltype(element)>; // TINYAD_VARIABLES_TYPE(element);

            // Get boundary annotations 
            Eigen::VectorXi& bound_face_idx = appState.bound_face_idx;
            Eigen::Index f_idx = element.handle;

            Eigen::Matrix<T, 2, 2> metric = appState.cur_surf->data().Bs[f_idx].transpose() * appState.cur_surf->data().Bs[f_idx];

            // Dirichlet (i.e. pinned) boundary condition
            if (bound_face_idx(f_idx) == 1)
            {
                ProcElement<T, VAR> e(ElementLiftType::primal);
                e.SetSelfData(appState, f_idx, element);

                Eigen::VectorX<T> curr = e.self_data.dofs_curr_elem;

                Eigen::VectorX<T> targ = appState.frames_orig.row(f_idx);

                T eucclidean_norm = 0;
                int nvecs = curr.rows() / 2;
                for (int i = 0; i < nvecs; i++) {
                    Eigen::Matrix<T, 2, 1> diff;
                    diff << curr[2 * i] - targ[2 * i], curr[2 * i + 1] - targ[2 * i + 1];
                    eucclidean_norm += diff.dot(metric * diff);
                }
                return e.w_bound * eucclidean_norm; 
            }
            return T(0);

        });
    }


    // Deps primal vars 
    static void addSmoothness_L4_Term(ADFunc& func, AppState& appState) {

        std::cout << "add L4 smoonthess obj" << std::endl;

        func.template add_elements<4>(TinyAD::range(appState.F.rows()), [&](auto& element)->TINYAD_SCALAR_TYPE(element)
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

            ProcElement<T, VAR> e(ElementLiftType::L4_krushkal);
            // e.setElementVars(appState, f_idx, s_curr);
            e.SetSelfData(appState, f_idx, element);
            e.SetNeighborData(appState, f_idx, element);


            ///////////////////
            //// Initlaize elementwise objectives 
            ///////////////////

            T dirichlet_term = e.GetL4DirichletTerm();

            appState.os->smoothness_L4(f_idx) = TinyAD::to_passive(dirichlet_term);

            T ret = T(0);//  delta_norm_term * delta_weight;
            ret = ret + e.w_attenuate * e.w_smooth * dirichlet_term;

            return ret;
        });

    }

    // Deps primal vars 
    static void addSmoothness_L2_Term(ADFunc& func, AppState& appState) {

        std::cout << "add L2 smoonthess obj" << std::endl;

        func.template add_elements<4>(TinyAD::range(appState.F.rows()), [&](auto& element)->TINYAD_SCALAR_TYPE(element)
        {

            // Evaluate element using either double or TinyAD::Double
            using T = TINYAD_SCALAR_TYPE(element);
            using VAR = std::decay_t<decltype(element)>;

            // Get variable 2D vertex positions
            Eigen::Index f_idx = element.handle;
            // Eigen::VectorX<T> s_curr = element.variables(f_idx);

            ProcElement<T, VAR> e(ElementLiftType::L2_krushkal);
            // e.setElementVars(appState, f_idx, s_curr);
            e.SetSelfData(appState, f_idx, element);



            Eigen::VectorXi& bound_face_idx = appState.bound_face_idx;


            // A bit hacky but exit early if on a boundary element.  Should really do it same as in matlab mint and make boundary elements distict from the mesh. 
            if (bound_face_idx(f_idx) == 1)
            {
                return T(0);
            }

            e.SetNeighborData(appState, f_idx, element);


            ///////////////////
            //// Initlaize elementwise objectives 
            ///////////////////
            T dirichlet_term = e.GetL2DirichletTerm();

            appState.os->smoothness_L2(f_idx) = TinyAD::to_passive(dirichlet_term);

            T delta_weight = 1.; // std::min(w_curl/100., 1./w_attenuate);

            T ret = T(0);//  delta_norm_term * delta_weight;
            if (e.w_smooth > 0)
            {
                ret = ret + e.w_attenuate * e.w_smooth * dirichlet_term;
            }
            return ret;
        });

    }

    // deps primal vars 
    static void addCurlTerm_L2(ADFunc& func, AppState& appState) {

        std::cout << "add L2 curl obj" << std::endl;

        func.template add_elements<4>(TinyAD::range(appState.F.rows()), [&](auto& element)->TINYAD_SCALAR_TYPE(element)
        {

            // Evaluate element using either double or TinyAD::Double
            using T = TINYAD_SCALAR_TYPE(element);
            using VAR = std::decay_t<decltype(element)>;



            // Get variable 2D vertex positions
            Eigen::Index f_idx = element.handle;

            ProcElement<T, VAR> e(ElementLiftType::L2_tensor);
            e.SetSelfData(appState, f_idx, element);



            Eigen::VectorXi bound_face_idx = appState.bound_face_idx;

            // A bit hacky but exit early if on a boundary element.  Should really do it same as in matlab mint and make boundary elements distict from the mesh. 
                    // Eigen::VectorXi bound_face_idx = appState.bound_face_idx;
            if (bound_face_idx(f_idx) == 1)
            {
                return T(0);
            }

            e.SetNeighborData(appState, f_idx, element);



            ///////////////////
            //// Curl Term 
            ///////////////////



            Eigen::Vector4d ea = appState.C_sym_2.row(appState.cur_surf->data().faceEdges(f_idx, 0));
            Eigen::Vector4d eb = appState.C_sym_2.row(appState.cur_surf->data().faceEdges(f_idx, 1));
            Eigen::Vector4d ec = appState.C_sym_2.row(appState.cur_surf->data().faceEdges(f_idx, 2));

            T curl_term = pow(ea.dot(e.neighbor_data.at(0).L_2_primals) - ea.dot(e.self_data.L_2_primals), 2);
            curl_term += pow(eb.dot(e.neighbor_data.at(1).L_2_primals) - eb.dot(e.self_data.L_2_primals), 2);
            curl_term += pow(ec.dot(e.neighbor_data.at(2).L_2_primals) - ec.dot(e.self_data.L_2_primals), 2);

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

        func.template add_elements<4>(TinyAD::range(appState.F.rows()), [&](auto& element)->TINYAD_SCALAR_TYPE(element)
        {

            // Evaluate element using either double or TinyAD::Double
            using T = TINYAD_SCALAR_TYPE(element);
            using VAR = std::decay_t<decltype(element)>;



            // Get variable 2D vertex positions
            Eigen::Index f_idx = element.handle;
            // Eigen::VectorX<T> s_curr = element.variables(f_idx);

            ProcElement<T, VAR> e(ElementLiftType::L4_tensor);
            // e.setElementVars(appState, f_idx, s_curr);
            e.SetSelfData(appState, f_idx, element);



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

            e.SetNeighborData(appState, f_idx, element);



            ///////////////////
            //// Curl Term 
            ///////////////////



            Eigen::VectorXd ea = appState.C_sym_4.row(appState.cur_surf->data().faceEdges(f_idx, 0));
            Eigen::VectorXd eb = appState.C_sym_4.row(appState.cur_surf->data().faceEdges(f_idx, 1));
            Eigen::VectorXd ec = appState.C_sym_4.row(appState.cur_surf->data().faceEdges(f_idx, 2));

            T curl_term = pow(ea.dot(e.neighbor_data.at(0).L_4_primals) - ea.dot(e.self_data.L_4_primals), 2);
            curl_term += pow(eb.dot(e.neighbor_data.at(1).L_4_primals) - eb.dot(e.self_data.L_4_primals), 2);
            curl_term += pow(ec.dot(e.neighbor_data.at(2).L_4_primals) - ec.dot(e.self_data.L_4_primals), 2);

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



