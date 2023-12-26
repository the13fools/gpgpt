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

enum class ElementLoadType {
    L2, L4, Element_COUNT
};



template <typename T_active>
class ElementData {
public:

    Eigen::Index curr_idx; //  faceNeighbors idx
    Eigen::Index faceNeighbor_idx; //  faceNeighbors idx
    
    Eigen::VectorX<T_active> dofs_curr_elem;
    std::vector< Eigen::VectorX<T_active> > primals_rank1;
    std::vector< T_active > primal_norms;
    
    T_active frame_norm_euclidian;

    Eigen::VectorX<T_active> primals;
    Eigen::VectorX<double> moments; // Moments only get updated during local step for ADMM style algos.  
    Eigen::VectorX<T_active> deltas;

    // These are the lifted primals to moment space.  
    // i.e. if primals_rank1 = [u, v];
    // then L_2_primals = flatten(u*u^T + v*v^T);
    // and L_4_primals = flatten(u^4 + v^4 = flatten(u*u^T)*flatten(u*u^T)^T + flatten(v*v^T)*flatten(v*v^T)^T);

    // For now we don't do the extra optimization of removing the symmetry from these 
    // entries, just as a base line, will add this in momentarily.  
    Eigen::VectorX<T_active> L_2_primals;
    Eigen::VectorX<T_active> L_2x2_primals;
    Eigen::VectorX<T_active> L_4_primals;

    Eigen::VectorX<T_active> L_2_krushkal;
    Eigen::VectorX<T_active> L_4_krushkal;

    void set_primals_rank1(DOFMemoryLayout& primals_layout)
    {
        int primals_size = primals_layout.size;
        int primal_rank = primals_layout.rank;

        int rank1_size = primals_size/primal_rank; // Elsewhere should check that this is an integer
        
        frame_norm_euclidian = 0;
        primals_rank1.resize(primal_rank);
        primal_norms.resize(primal_rank);
        for (int v_i = 0; v_i < primal_rank; v_i++)
        {
            primals_rank1[v_i] = dofs_curr_elem.segment(primals_layout.start + v_i*rank1_size, rank1_size);
            primal_norms[v_i] = primals_rank1[v_i].squaredNorm();
            frame_norm_euclidian += primals_rank1[v_i].squaredNorm();
        }

        primals.resize(primals_size);
        primals = dofs_curr_elem.segment(primals_layout.start, primals_size);

    };
    

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

    // Surface* cur_surf;
    CubeCover::TetMeshConnectivity* cur_mesh;

    // Eigen::Index f_idx; // maybe get rid of this it's redundant. 
    Eigen::Index t_idx; // maybe get rid of this it's redundant. 


    ElementData<T_active> self_data;

    std::vector<ElementData<T_active>> neighbor_data;

    // static const Eigen::VectorX<int> L_2_weights << 1, 2, 1;
    // static const Eigen::VectorX<int> L_4_weights << 1, 4, 6, 4, 1;
    // static const Eigen::VectorX<int> L_6_weights << 1, 6, 15, 20, 15, 6, 1;

    // Eigen::MatrixX<T_active> currcurr;
    // Eigen::VectorX<T_active> currcurrt;


    void setSelfData(AppState& appState, const Eigen::Index t_idx, ELEM& element) { 

        setElemState(appState, t_idx);

        self_data.dofs_curr_elem = element.variables(t_idx);

        self_data.set_primals_rank1(appState.primals_layout);
        
        self_data.curr_idx = t_idx;

        




        L2_primals(appState, self_data.dofs_curr_elem, self_data);
        // L4_primals(appState, f_idx, self_data.dofs_curr_elem, self_data);

      

        self_data.frame_norm_euclidian = 1;

        if (self_data.frame_norm_euclidian < 1e-3) {
            // std::cout << "Warning: frame_norm_euclidian is small: " << self_data.frame_norm_euclidian << std::endl;
            
            // std::cout << L_2_primals.rows() << std::endl;
            // if (L_2_primals.rows() == 0)
            // {
            //     throw std::runtime_error("L_2_primals.rows() == 0");
            // }
            // std::cout << L_2_primals(0) << " " << L_2_primals(1) << " " << L_2_primals(2) << " " << L_2_primals(3) << std::endl;
            // frame_norm_euclidian = 1e-3;
        }


    }

    void setNeighborData(AppState& appState, const Eigen::Index t_idx, ELEM& element) { 

        neighbor_data.clear();
        int max_neighbors = 4;

        num_neighbors = 0;

        for (int i = 0; i < max_neighbors; i++)
        {
            ElementData<T_active> neighbor_data_i;
            // int n_idx = cur_surf->data().faceNeighbors(f_idx, i); // 2d change 
            int n_idx = cur_mesh->tetOppositeVertex(t_idx, i);

            neighbor_data_i.curr_idx = n_idx;
            neighbor_data_i.faceNeighbor_idx = i;
            if (n_idx == -1) continue; // Make sure this is the correct convention.  Do more advance boundary handling later.
            
            num_neighbors += 1;
            neighbor_data_i.dofs_curr_elem = element.variables(cur_mesh->tetOppositeVertex(t_idx, i));
            neighbor_data_i.set_primals_rank1(appState.primals_layout);

            L2_primals(appState, neighbor_data_i.dofs_curr_elem, neighbor_data_i);
            // L4_primals(appState, f_idx, neighbor_data_i.dofs_curr_elem, neighbor_data_i);
            neighbor_data.push_back(neighbor_data_i);
        }

    }



    void setElemState(AppState& appState, const Eigen::Index t_idx) { 
        // cur_surf = appState.cur_surf.get();
        cur_mesh = appState.cur_tet_mesh.get();

        w_bound = appState.config->w_bound;
        w_smooth_vector = appState.config->w_smooth_vector;

        w_smooth = appState.config->w_smooth;
        w_curl = appState.config->w_curl;
        w_attenuate = appState.config->w_attenuate;

        this->t_idx = t_idx;
    }


    // void setElementVars(AppState& appState, const Eigen::Index& f_idx, const Eigen::VectorX<T_active>& s_curr) { 
    void L2_primals(AppState& appState, const Eigen::VectorX<T_active>& s_curr, ElementData<T_active>& data) { 
        
        // Seperate out the primal dofs into rank-1 components
        // Maybe make this a seperate function. 
        int nprimals = data.primals_rank1.size();
        int primals_size = data.primals_rank1[0].size();
        data.L_2_primals.resize(primals_size*primals_size);
        data.L_2_primals.setZero();
        data.L_2_krushkal = data.L_2_primals;

        // update for n vectors per frame.  
        for (int v_i = 0; v_i < nprimals; v_i++)
        {
            Eigen::VectorX<T_active> cur = data.primals_rank1[v_i];
            T_active cur_norm = cur.norm() + 1e-10;
            Eigen::VectorX<T_active> cur_normalized = cur / cur_norm;

            // Eigen::MatrixX<T_active> curcurt = cur*cur.transpose();
            Eigen::VectorX<T_active> curcurt_flattened;
            curcurt_flattened.resize(cur.rows()*cur.rows());// = Eigen::Zeros(cur.rows()*cur.rows()); // TODO: make this compressed 
            Eigen::VectorX<T_active> curcurt_normalized_flattened = curcurt_flattened;
            // TODO fix this later
            for (int i = 0; i < primals_size; i++)
            {
                for (int j = 0; j < primals_size; j++)
                {
                    curcurt_flattened(i*primals_size + j) = cur(i)*cur(j);
                    curcurt_normalized_flattened(i*primals_size + j) = cur_normalized(i)*cur_normalized(j);
                }
            }

            data.L_2_primals = data.L_2_primals + curcurt_flattened;
            data.L_2_krushkal = data.L_2_krushkal + curcurt_normalized_flattened * cur_norm;
        
        }

        // This sets the squared values which scales in the same way as L_4 to make the energy scale invariant. 
        L2_to_L2x2(data);

    }


    // There's some cleverer thing to do here to save some computation.  
    // Look at the matlab mint repo.  Do this for now cause it's easy. 
    void L2_to_L2x2(ElementData<T_active>& data)
    {
        int L2_size = data.L_2_primals.rows();
        int L2x2_size = L2_size*L2_size;
        data.L_2x2_primals.resize(L2x2_size);
        data.L_2x2_primals.setZero();

        Eigen::VectorX<T_active> L_2 = data.L_2_primals;
            // Eigen::MatrixX<T_active> curcurt = cur*cur.transpose();
        Eigen::VectorX<T_active> L_2x2 = data.L_2x2_primals;

        for (int i = 0; i < L2_size; i++)
        {
            for (int j = 0; j < L2_size; j++)
            {
                L_2x2(i*L2_size + j) = L_2(i)*L_2(j);
            }
        }

        data.L_2x2_primals = L_2x2;

    }

        // void setElementVars(AppState& appState, const Eigen::Index& f_idx, const Eigen::VectorX<T_active>& s_curr) { 
    void L4_primals(AppState& appState, const Eigen::Index f_idx, const Eigen::VectorX<T_active>& s_curr, ElementData<T_active>& data) { 
        
        // Seperate out the primal dofs into rank-1 components
        // Maybe make this a seperate function. 
        int nprimals = data.primals_rank1.size();
        int primals_size = data.primals_rank1[0].size();
        int L4_size = primals_size*primals_size*primals_size*primals_size;
        data.L_4_primals.resize(primals_size*primals_size*primals_size*primals_size);
        data.L_4_primals.setZero();
        data.L_4_krushkal = data.L_4_primals;


        // update for n vectors per frame.  
        for (int v_i = 0; v_i < nprimals; v_i++)
        {
            Eigen::VectorX<T_active> cur = data.primals_rank1[v_i];
            
            T_active cur_norm = cur.norm() + 1e-10;
            Eigen::VectorX<T_active> cur_normalized = cur / cur_norm;


            // Eigen::MatrixX<T_active> curcurt = cur*cur.transpose();
            Eigen::VectorX<T_active> curcurt_flattened;
            curcurt_flattened.resize(L4_size);// = Eigen::Zeros(cur.rows()*cur.rows()); // TODO: make this compressed 
            Eigen::VectorX<T_active> curcurt_normalized_flattened = curcurt_flattened;
            
            // TODO fix this later
            int psize = primals_size;
            int psize2 = psize*psize;
            int psize3 = psize2*psize;
            for (int i = 0; i < primals_size; i++)
            {
                for (int j = 0; j < primals_size; j++)
                {
                    for (int k = 0; k < primals_size; k++)
                    {
                        for (int l = 0; l < primals_size; l++)
                        {
                            curcurt_flattened(i*psize3 + j*psize2 + k*psize + l) = cur(i)*cur(j)*cur(k)*cur(l);
                            curcurt_normalized_flattened(i*psize3 + j*psize2 + k*psize + l)  = cur_normalized(i)*cur_normalized(j)*cur_normalized(k)*cur_normalized(l); 
                        }
                    }
                }
            }

            data.L_4_primals = data.L_4_primals + curcurt_flattened;
            data.L_4_krushkal = data.L_4_krushkal + curcurt_normalized_flattened * cur_norm;
        
        }

        // This sets the squared values which scales in the same way as L_4 to make the energy scale invariant. 
        // L2_to_L2x2(data);

    }



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
