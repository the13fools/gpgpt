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

            Eigen::Index f_idx = element.handle;

            Eigen::VectorXi& bound_face_idx = appState.bound_face_idx;

            // A bit hacky but exit early if on a boundary element.  Should really do it same as in matlab mint and make boundary elements distict from the mesh. 
            if (bound_face_idx(f_idx) == 1)
            {
                return T(0);
            }

            auto& V = appState.cur_surf->data().V;
            auto& E = appState.cur_surf->data().E;

            Eigen::VectorX<T> cur_face_dofs = element.variables(f_idx);    // the vector fields on the current face

            T curl_term = T(0);

            for (int i = 0; i < 3; i++) {
                int eid = appState.cur_surf->data().faceEdges(f_idx, i);
                assert(E(eid, 0) != -1 || E(eid, 1) != -1); // should not be the boundary edge

                int efid = E(eid, 0) == f_idx ? 0 : 1;
                int n_idx = E(eid, 1 - efid);       // neighboring face id

                Eigen::VectorX<T> nei_face_dofs = element.variables(n_idx);

                Eigen::Vector<T, 3> edge = V.row(appState.cur_surf->data().edgeVerts(eid, 0)) - V.row(appState.cur_surf->data().edgeVerts(eid, 1));

                assert(cur_face_dofs.size() == nei_face_dofs.size() && cur_face_dofs.rows() % 2 == 0);

                int nvecs = cur_face_dofs.rows() / 2;
                T diff = T(0);
                for (int vi = 0; vi < nvecs; vi++) {
                    Eigen::Vector<T, 3> ext_v = appState.cur_surf->data().Bs[f_idx] * cur_face_dofs.template segment<2>(2 * vi);
                    Eigen::Vector<T, 3> ext_nv = appState.cur_surf->data().Bs[n_idx] * nei_face_dofs.template segment<2>(2 * vi);
                    diff += pow(ext_v.dot(edge), 2) - pow(ext_nv.dot(edge), 2);
                }

                curl_term += pow(diff, 2);
            }

            T ret = T(0);
            if (appState.config->w_curl > 0)
            {
                ret = ret + appState.config->w_curl * curl_term;
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

            Eigen::Index f_idx = element.handle;

            Eigen::VectorXi& bound_face_idx = appState.bound_face_idx;

            // A bit hacky but exit early if on a boundary element.  Should really do it same as in matlab mint and make boundary elements distict from the mesh. 
            if (bound_face_idx(f_idx) == 1)
            {
                return T(0);
            }

            auto& V = appState.cur_surf->data().V;
            auto& E = appState.cur_surf->data().E;

            Eigen::VectorX<T> cur_face_dofs = element.variables(f_idx);    // the vector fields on the current face

            T curl_term = T(0);

            for (int i = 0; i < 3; i++) {
                int eid = appState.cur_surf->data().faceEdges(f_idx, i);
                assert(E(eid, 0) != -1 || E(eid, 1) != -1); // should not be the boundary edge

                int efid = E(eid, 0) == f_idx ? 0 : 1;
                int n_idx = E(eid, 1 - efid);       // neighboring face id

                Eigen::VectorX<T> nei_face_dofs = element.variables(n_idx);

                Eigen::Vector<T, 3> edge = V.row(appState.cur_surf->data().edgeVerts(eid, 0)) - V.row(appState.cur_surf->data().edgeVerts(eid, 1));

                assert(cur_face_dofs.size() == nei_face_dofs.size() && cur_face_dofs.rows() % 2 == 0);

                int nvecs = cur_face_dofs.rows() / 2;
                T diff = T(0);
                for (int vi = 0; vi < nvecs; vi++) {
                    Eigen::Vector<T, 3> ext_v = appState.cur_surf->data().Bs[f_idx] * cur_face_dofs.template segment<2>(2 * vi);
                    Eigen::Vector<T, 3> ext_nv = appState.cur_surf->data().Bs[n_idx] * nei_face_dofs.template segment<2>(2 * vi);
                    diff += pow(ext_v.dot(edge), 4) - pow(ext_nv.dot(edge), 4);
                }

                curl_term += pow(diff, 2);
            }

            T ret = T(0);
            if (appState.config->w_curl > 0)
            {
                ret = ret + appState.config->w_curl * curl_term;
            }
            return ret;

        });
    }

}; // namespace OptZoo

#endif // OPTZOO_H



