#ifndef ELEMENTVARSTINYAD_H
#define ELEMENTVARSTINYAD_H

#include <string>

// #include "UtilsMisc.h"

// This is a helper header for the TinyAD interface of gpgpt.
// We suggest that the pattern modeled here is a convenient way to structure 
// generic geometric optimization problems.  

// As you'll see this enables a significant amount of code reuse in the OptZoo downstream.

/// Can make this more efficient to cache these values per element each step
// Not going to worry about doing this for now.  
// There's a non-zero performance gain to be had here though.  Pretty low hanging thing to hack in maybe?



/**
 * ElementVars class 
 */
template <typename T_active>
class ElementVars {
public:

    // ElementVars(AppState& appState) {}

    Eigen::Index f_idx; 
    Eigen::VectorX<T_active> s_curr;

    // Eigen::Vector2<T> curr;
    // Eigen::Vector4<T> delta;
    Eigen::VectorX<T_active> curr;
    Eigen::VectorX<T_active> delta;

    Eigen::MatrixX<T_active> currcurr;
    Eigen::VectorX<T_active> currcurrt;

    Surface* cur_surf;

    // T w_bound;

    Eigen::VectorXi bound_face_idx;

    double w_bound;
    double w_smooth_vector;

    double w_smooth;
    double w_curl;
    double w_attenuate;

    inline void flatten(const Eigen::MatrixX<T_active>& xxt, Eigen::VectorX<T_active>& xxt_flattened)
    {
        int xdim = xxt.rows();

        xxt_flattened.resize(xdim*xdim);

        for (int i = 0; i < xdim; i++)
        {
            for (int j = 0; j < xdim; j++)
            {
                xxt_flattened(i*xdim + j) = xxt(i)*xxt(j);
            }
        }


    }



    void setElementVars(AppState& appState, const Eigen::Index& f_idx, const Eigen::VectorX<T_active>& s_curr) { 

        int primals_size = appState.primals_layout.size;
        int deltas_size = appState.deltas_layout.size;

        curr.resize(primals_size);
        delta.resize(deltas_size);

        curr =  s_curr.segment(appState.primals_layout.start, primals_size); // head(2);
        delta =  s_curr.segment(appState.deltas_layout.start, deltas_size); // head(2);

        // std::cout << "s_curr " << s_curr << std::endl;
        // std::cout << "curr " << curr << "curr.rows(): " <<  curr.rows() << std::endl;
        // std::cout << "delta " << delta << "delta.rows(): " << delta.rows() <<  std::endl;

        // f_idx = element.handle;
        // s_curr = element.variables(f_idx);
        // curr =  s_curr.head(2);
        // delta = s_curr.tail(4);

        currcurr = curr*curr.transpose();

        // std::cout << "currcurr " << currcurr << std::endl;
        // std::cout << "currcurr rows " << currcurr.rows() << std::endl;
        // std::cout << "currcurr cols " << currcurr.cols() << std::endl;



        currcurrt.resize(primals_size*primals_size);

        flatten(currcurr, currcurrt);

        // for (int i = 0; i < primals_size; i++)
        // {
        //     for (int j = 0; j < primals_size; j++)
        //     {
        //         currcurrt(i*primals_size + j) = curr(i)*curr(j);
        //     }
        // }
        // // currcurrt = flatten(currcurr);

        // std::cout << "currcurrt " << currcurrt << std::endl;

        cur_surf = appState.cur_surf;

        // T w_bound = appState.config->w_bound;

        bound_face_idx = appState.bound_face_idx;

        w_bound = appState.config->w_bound;
        w_smooth_vector = appState.config->w_smooth_vector;

        w_smooth = appState.config->w_smooth;
        w_curl = appState.config->w_curl;
        w_attenuate = appState.config->w_attenuate;
    }

//     void setNeighborVars()
//     {
//         ///////////////////
// //// Initialize the neighbor meta-data 
// ///////////////////

//           Eigen::VectorX<T> s_a = element.variables(e.cur_surf->data().faceNeighbors(f_idx, 0));
//           Eigen::VectorX<T> s_b = element.variables(e.cur_surf->data().faceNeighbors(f_idx, 1));
//           Eigen::VectorX<T> s_c = element.variables(e.cur_surf->data().faceNeighbors(f_idx, 2));



//           Eigen::Vector2<T> a = s_a.head(2);
//           Eigen::Vector2<T> b = s_b.head(2);
//           Eigen::Vector2<T> c = s_c.head(2);

//           Eigen::Matrix2<T> aa = a*a.transpose();
//           Eigen::Matrix2<T> bb = b*b.transpose();
//           Eigen::Matrix2<T> cc = c*c.transpose();


//           Eigen::Vector4<T> a_delta = s_a.tail(4);
//           Eigen::Vector4<T> b_delta = s_b.tail(4);
//           Eigen::Vector4<T> c_delta = s_c.tail(4);

//           Eigen::Vector4<T> aat = flatten(aa);
//           Eigen::Vector4<T> bbt = flatten(bb);
//           Eigen::Vector4<T> cct = flatten(cc);

//           aat = aat + a_delta;
//           bbt = bbt + b_delta; 
//           cct = cct + c_delta;
//     }



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
