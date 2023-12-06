#ifndef ELEMENTVARSTINYAD_H
#define ELEMENTVARSTINYAD_H

#include <string>

#include <TinyAD/ScalarFunction.hh>

// #include "UtilsMisc.h"

// This is a helper header for the TinyAD interface of gpgpt.
// We suggest that the pattern modeled here is a convenient way to structure 
// generic geometric optimization problems.  

// As you'll see this enables a significant amount of code reuse in the OptZoo downstream.

/// Can make this more efficient to cache these values per element each step
// Not going to worry about doing this for now.  
// There's a non-zero performance gain to be had here though.  Pretty low hanging thing to hack in maybe?


template <typename T_active>
class ElementData {
public:

    // Eigen::Index f_idx; 
    Eigen::Index n_idx; //  faceNeighbors idx
    Eigen::VectorX<T_active> dofs_curr_elem;

    // Eigen::Vector2<T> curr;
    // Eigen::Vector4<T> delta;
    Eigen::VectorX<T_active> curr;
    Eigen::VectorX<T_active> delta;

    Eigen::MatrixX<T_active> currcurr;
    Eigen::VectorX<T_active> currcurrt;

        // Eigen::VectorXi bound_face_idx;

};



/**
 * ElementVars class 
 */
template <typename T_active, typename ELEM>
class ProcElement {
public:

    ProcElement() : w_bound(0) {};

    // ElementVars(AppState& appState) {}

    double w_bound;
    double w_smooth_vector;

    double w_smooth;
    double w_curl;
    double w_attenuate;

    int num_neighbors = 0;

    Surface* cur_surf;

    // Replace this with ElementData
    Eigen::Index f_idx; 
    // Eigen::VectorX<T_active> s_curr;

    // // Eigen::Vector2<T> curr;
    // // Eigen::Vector4<T> delta;
    // Eigen::VectorX<T_active> curr;
    // Eigen::VectorX<T_active> delta;

    // Eigen::MatrixX<T_active> currcurr;
    // Eigen::VectorX<T_active> currcurrt;


    ElementData<T_active> self_data;

    std::vector<ElementData<T_active>> neighbor_data;

    
    // T w_bound;


    void setSelfData(AppState& appState, const Eigen::Index f_idx, ELEM& element) { 

        setElemState(appState, f_idx);

        self_data.dofs_curr_elem = element.variables(f_idx);
        setL2Vars(appState, f_idx, self_data.dofs_curr_elem, self_data);

    }

    void setNeighborData(AppState& appState, const Eigen::Index f_idx, ELEM& element) { 

        neighbor_data.clear();
        int max_neighbors = 3;

        num_neighbors = 0;

        for (int i = 0; i < max_neighbors; i++)
        {
            ElementData<T_active> neighbor_data_i;
            int n_idx = cur_surf->data().faceNeighbors(f_idx, i);
            neighbor_data_i.n_idx = n_idx;
            if (n_idx == -1) continue; // Make sure this is the correct convention.  Do more advance boundary handling later.
            
            num_neighbors += 1;
            neighbor_data_i.dofs_curr_elem = element.variables(cur_surf->data().faceNeighbors(f_idx, i));
            setL2Vars(appState, f_idx, neighbor_data_i.dofs_curr_elem, neighbor_data_i);
            neighbor_data.push_back(neighbor_data_i);
        }

    }



    void setElemState(AppState& appState, const Eigen::Index f_idx) { 
        cur_surf = appState.cur_surf;

        w_bound = appState.config->w_bound;
        w_smooth_vector = appState.config->w_smooth_vector;

        w_smooth = appState.config->w_smooth;
        w_curl = appState.config->w_curl;
        w_attenuate = appState.config->w_attenuate;

        this->f_idx = f_idx;
    }


    // void setElementVars(AppState& appState, const Eigen::Index& f_idx, const Eigen::VectorX<T_active>& s_curr) { 
    void setL2Vars(AppState& appState, const Eigen::Index f_idx, const Eigen::VectorX<T_active>& s_curr, ElementData<T_active>& data) { 
        
        int primals_size = appState.primals_layout.size;
        int deltas_size = appState.deltas_layout.size;


        data.curr.resize(primals_size);
        data.delta.resize(deltas_size);

        
        data.curr =  s_curr.segment(appState.primals_layout.start, primals_size); // head(2);
        data.delta =  s_curr.segment(appState.deltas_layout.start, deltas_size); // head(2);


        data.currcurr = data.curr*data.curr.transpose();

        data.currcurrt.resize(primals_size*primals_size);

        // flatten(currcurr, currcurrt);
        // TODO fix this later, not sure why this function isn't having it.  
        // maybe just copy over the one from UtilsMisc
        for (int i = 0; i < primals_size; i++)
        {
            for (int j = 0; j < primals_size; j++)
            {
                data.currcurrt(i*primals_size + j) = data.curr(i)*data.curr(j);
            }
        }

    }

//     void setNeighborVars()
//     {
//         ///////////////////
// //// Initialize the neighbor meta-data 
// ///////////////////

//           Eigen::VectorX<T> s_a = element.variables(cur_surf->data().faceNeighbors(f_idx, 0));
//           Eigen::VectorX<T> s_b = element.variables(cur_surf->data().faceNeighbors(f_idx, 1));
//           Eigen::VectorX<T> s_c = element.variables(cur_surf->data().faceNeighbors(f_idx, 2));



//           Eigen::Vector2<T> a = s_a.head(2);
//           Eigen::Vector2<T> b = s_b.head(2);
//           Eigen::Vector2<T> c = s_c.head(2);

//           Eigen::Matrix2<T> aa = a*a.transpose();
//           Eigen::Matrix2<T> bb = b*b.transpose();
//           Eigen::Matrix2<T> cc = c*c.transpose();


//           Eigen::Vector4<T> a_delta = s_a.tail(4);
//           Eigen::Vector4<T> b_delta = s_b.tail(4);
//           Eigen::Vector4<T> c_delta = s_c.tail(4);

//         //   Eigen::Vector4<T> aat = flatten(aa);
//         //   Eigen::Vector4<T> bbt = flatten(bb);
//         //   Eigen::Vector4<T> cct = flatten(cc);
//                 currcurrt.resize(primals_size*primals_size);

//         // flatten(currcurr, currcurrt);
//         // TODO fix this later, not sure why this function isn't having it.  
//         // maybe just copy over the one from UtilsMisc
//         for (int i = 0; i < primals_size; i++)
//         {
//             for (int j = 0; j < primals_size; j++)
//             {
//                 currcurrt(i*primals_size + j) = curr(i)*curr(j);
//             }
//         }

//           aat = aat + a_delta;
//           bbt = bbt + b_delta; 
//           cct = cct + c_delta;
//     }

    // inline void flatten(const Eigen::MatrixX<T_active>& xxt, Eigen::VectorX<T_active>& xxt_flattened)
    // {
    //     int xdim = xxt.rows();

    //     xxt_flattened.resize(xdim*xdim);

    //     for (int i = 0; i < xdim; i++)
    //     {
    //         for (int j = 0; j < xdim; j++)
    //         {
    //             xxt_flattened(i*xdim + j) = xxt(i)*xxt(j);
    //         }
    //     }


    // }



/**
 * Convenience Struct for computing the element metadata for the objective functions.
 * 
 * @param AppData, type, element number  
 * @return ElementVars struct.  
 */
// std::string fieldViewToString(Field_View view);
// std::string fieldViewToFileStub(Field_View view);


//  };

};

#endif // ELEMENTVARSTINYAD_H
