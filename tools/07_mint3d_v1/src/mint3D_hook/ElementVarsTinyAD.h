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


// Borrowing the notation from here: https://www.tensortoolbox.org/tensor_types.html
// Primal means initilize vector field variables without doing any lifting 
// Tensor is a multidimensional array, i.e. v*v'.  
// Krushkal is a normalized represenation which is more suitable for geometric optimization, i.e. |v| \hat{v}\hat{v}'
// TODO sym_tensor and sym_krushkal, which would both offer a non-trivial performance improvement 
enum class ElementLiftType {
    primal, 
    Lk_edge_contract,
    L2_krushkal, 
    L2_sym_krushkal,
    L2_facets, 
    L4_krushkal, 
    L4_sym_krushkal,
    L4_facets, 
    Element_COUNT
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
    // Eigen::VectorX<T_active> L_2_primals;
    // Eigen::VectorX<T_active> L_2x2_primals;
    // Eigen::VectorX<T_active> L_4_primals;
    Eigen::VectorX<T_active> L_2_facet_diff;
    Eigen::VectorX<T_active> L_4_facet_diff;
    // Eigen::VectorX<T_active> L_2x2_primals;
    // Eigen::VectorX<T_active> L_4_primals;

    Eigen::VectorX<T_active> L_2_krushkal;
    Eigen::VectorX<T_active> L_2_sym_krushkal;

    Eigen::VectorX<T_active> L_4_krushkal;
    Eigen::VectorX<T_active> L_4_sym_krushkal;


    Eigen::VectorX<T_active> Lk_edge_contract_diff;

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

    ProcElement(ElementLiftType t) : curr_lift(t), w_bound(0) {};

    // ElementVars(AppState& appState) {}

    ElementLiftType curr_lift; 

    double w_bound;
    double w_smooth_vector;

    double w_smooth;
    double w_curl;
    double w_attenuate;

    int num_neighbors = 0;

    // Surface* cur_surf;
    CubeCover::TetMeshConnectivity* cur_mesh;

    Eigen::Index t_idx; // maybe get rid of this, it's redundant. 

    int edge_contract_order = -1;


    ElementData<T_active> self_data;

    std::vector<ElementData<T_active>> neighbor_data;


    void setSelfData(AppState& appState, const Eigen::Index t_idx, ELEM& element) { 

        setElemState(appState, t_idx);

        self_data.dofs_curr_elem = element.variables(t_idx);

        self_data.set_primals_rank1(appState.primals_layout);
        
        self_data.curr_idx = t_idx;

        

        switch(curr_lift)
        {
            // Can seperate these two out to make it more granular if it's necessary for a bit of a speed boost
            case(ElementLiftType::primal):
                break;

            case(ElementLiftType::L2_krushkal):
                L2_krushkal(appState, self_data);
                break;  
            case(ElementLiftType::L2_sym_krushkal):
                L2_sym_krushkal(appState, self_data);
                break;  
            case(ElementLiftType::Lk_edge_contract):
            case(ElementLiftType::L2_facets): 
            case(ElementLiftType::L4_facets):  
            // Computation for this happens per neighbor      
                break;       

            case(ElementLiftType::L4_krushkal):
                L4_krushkal(appState, self_data);
                break;
            case(ElementLiftType::L4_sym_krushkal):
                L4_sym_krushkal(appState, self_data);
                break;
        
            default:
                std::cout << "Error: LiftType not implemented" << std::endl;
        }


    }

    void setNeighborData(AppState& appState, const Eigen::Index t_idx, ELEM& element) { 

        neighbor_data.clear();
        int max_neighbors = 4;

        num_neighbors = 0;

        int pinned_neighbors = 0;

        for (int i = 0; i < max_neighbors; i++)
        {
            ElementData<T_active> neighbor_data_i;
            // int n_idx = cur_surf->data().faceNeighbors(f_idx, i); // 2d change 

            int n_idx;

            if (t_idx >= cur_mesh->nTets())
            {
                n_idx = cur_mesh->boundaryElementTet(t_idx - cur_mesh->nTets());
                max_neighbors = 1;
            }
            else
            {
                n_idx = cur_mesh->tetOppositeVertex(t_idx, i);
                if (n_idx == -1) // continue; 
                {
                    // std::cout <<" current face " << cur_mesh->tetFace( t_idx, i ) << std::endl;
                    // std::cout <<" current face orientation " << cur_mesh->tetFaceOrientation( t_idx, i ) << std::endl;

                    n_idx = cur_mesh->faceBoundaryElement( cur_mesh->tetFace( t_idx, i ), (cur_mesh->tetFaceOrientation(t_idx, i) + 1) % 2  ); 
                    // std::cout << n_idx << " current n_idx" << std::endl;
                }

            }


            neighbor_data_i.curr_idx = n_idx;
            neighbor_data_i.faceNeighbor_idx = i;
            
            // Make sure this is the correct convention.  Do more advance boundary handling later.
            
            // if(appState.bound_face_idx(n_idx) == 1 && curr_lift == ElementLiftType::Lk_edge_contract) // && edge_contract_order != 2)
            // {
            //     // continue;
            //     pinned_neighbors++;
            //     if (pinned_neighbors > 2)
            //         continue;
            // }

            // if(appState.bound_face_idx(n_idx) == 1 && curr_lift == ElementLiftType::Lk_edge_contract && edge_contract_order != 2)
            // {
            //     continue;
            //     pinned_neighbors++;
            //     if (pinned_neighbors > 2)
            //         continue;
            // }


            num_neighbors += 1;
            neighbor_data_i.dofs_curr_elem = element.variables(n_idx);
            
            // if(appState.bound_face_idx(n_idx) == 1)
            // {
            //     neighbor_data_i.dofs_curr_elem = appState.frames_orig.row(n_idx);
            // }

            neighbor_data_i.set_primals_rank1(appState.primals_layout);

            switch(curr_lift)
            {
                case(ElementLiftType::primal):
                    break;
                // Can seperate these two out to make it more granular if it's necessary for a bit of a speed boost
                case(ElementLiftType::Lk_edge_contract):
                    // std::cout << "Lk_edge_contract" << std::endl;
                    // std::cout << edge_contract_order << std::endl;
                    Lk_edge_contract(appState, i, neighbor_data_i);

                    // if(appState.bound_face_idx(n_idx) == 1)
                    // {
                    //     neighbor_data_i.Lk_edge_contract_diff *= appState.config->w_attenuate / appState.config->w_curl * 10;
                    // }
                    break;
                
                case(ElementLiftType::L2_krushkal):
                    L2_krushkal(appState, neighbor_data_i);
                    break;
                case(ElementLiftType::L2_sym_krushkal):
                    L2_sym_krushkal(appState, neighbor_data_i);
                    break;
                case(ElementLiftType::L2_facets):
                    L2_facet_diff(appState, i, neighbor_data_i);
                    break;

                case(ElementLiftType::L4_sym_krushkal):
                    L4_sym_krushkal(appState, neighbor_data_i);
                    break;
                case(ElementLiftType::L4_krushkal):
                    L4_krushkal(appState, neighbor_data_i);
                    break;
                case(ElementLiftType::L4_facets):
                    L4_facet_diff(appState, i, neighbor_data_i);
                    
                break;

                default:
                    std::cout << "Error: LiftType not implemented" << std::endl;
            }

            // L2_primals(appState, neighbor_data_i.dofs_curr_elem, neighbor_data_i);
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
    void L2_sym_krushkal(AppState& appState, ElementData<T_active>& data) { 
        
        // Seperate out the primal dofs into rank-1 components
        // Maybe make this a seperate function. 
        int nprimals = data.primals_rank1.size();
        int primals_size = data.primals_rank1[0].size();
        int sym_size = 6;
        data.L_2_sym_krushkal.resize(sym_size);
        data.L_2_sym_krushkal.setZero();

        // update for n vectors per frame.  
        for (int v_i = 0; v_i < nprimals; v_i++)
        {
            Eigen::VectorX<T_active> cur = data.primals_rank1[v_i];
            T_active cur_norm = cur.norm() + 1e-10;
            Eigen::VectorX<T_active> cur_normalized = cur / cur_norm;

            // Eigen::MatrixX<T_active> curcurt = cur*cur.transpose();
            Eigen::VectorX<T_active> curcurt_normalized_flattened;
            curcurt_normalized_flattened.resize(sym_size);// = Eigen::Zeros(cur.rows()*cur.rows()); // TODO: make this compressed 
            // Eigen::VectorX<T_active> curcurt_normalized_flattened = curcurt_flattened;
            // TODO fix this later
            int cur_idx = 0;
            for (int i = 0; i < primals_size; i++)
            {
                for (int j = i; j < primals_size; j++)
                {
                    // curcurt_flattened(i*primals_size + j) = cur(i)*cur(j);
                    curcurt_normalized_flattened(cur_idx) = cur_normalized(i)*cur_normalized(j);
                    cur_idx++;
                }
            }

            // data.L_2_primals = data.L_2_primals + curcurt_flattened;
            data.L_2_sym_krushkal = data.L_2_sym_krushkal + curcurt_normalized_flattened * cur_norm;
        
        }

    }

    // void setElementVars(AppState& appState, const Eigen::Index& f_idx, const Eigen::VectorX<T_active>& s_curr) { 
    void L2_krushkal(AppState& appState, ElementData<T_active>& data) { 
        
        // Seperate out the primal dofs into rank-1 components
        // Maybe make this a seperate function. 
        int nprimals = data.primals_rank1.size();
        int primals_size = data.primals_rank1[0].size();
        data.L_2_krushkal.resize(primals_size*primals_size);
        data.L_2_krushkal.setZero();

        // update for n vectors per frame.  
        for (int v_i = 0; v_i < nprimals; v_i++)
        {
            Eigen::VectorX<T_active> cur = data.primals_rank1[v_i];
            T_active cur_norm = cur.norm() + 1e-10;
            Eigen::VectorX<T_active> cur_normalized = cur / cur_norm;

            // Eigen::MatrixX<T_active> curcurt = cur*cur.transpose();
            Eigen::VectorX<T_active> curcurt_normalized_flattened;
            curcurt_normalized_flattened.resize(cur.rows()*cur.rows());// = Eigen::Zeros(cur.rows()*cur.rows()); // TODO: make this compressed 
            // Eigen::VectorX<T_active> curcurt_normalized_flattened = curcurt_flattened;
            // TODO fix this later
            for (int i = 0; i < primals_size; i++)
            {
                for (int j = 0; j < primals_size; j++)
                {
                    // curcurt_flattened(i*primals_size + j) = cur(i)*cur(j);
                    curcurt_normalized_flattened(i*primals_size + j) = cur_normalized(i)*cur_normalized(j);
                }
            }

            // data.L_2_primals = data.L_2_primals + curcurt_flattened;
            data.L_2_krushkal = data.L_2_krushkal + curcurt_normalized_flattened * cur_norm;
        
        }

    }

// This is being called twice as much as it needs to be, there's a significant savings that could be
// had from doing this per dual edge instead of per tet... 
    void Lk_edge_contract(AppState& appState, const int n_idx, ElementData<T_active>& data) { 
        
        // Seperate out the primal dofs into rank-1 components
        // Maybe make this a seperate function. 
        int nprimals = data.primals_rank1.size();
        int facet_dim = data.primals_rank1[0].size() - 1;
        data.Lk_edge_contract_diff.resize(edge_contract_order+1);
        data.Lk_edge_contract_diff.setZero();

        // if (t_idx >= cur_mesh->nTets())
        // {
        //     return;
        // }
        
        Eigen::MatrixXd tet_facet_basis = appState.tet_facet_basis.at(t_idx).at(n_idx);

        Eigen::VectorX<T_active> neighbor_edge_proj_vals = Eigen::VectorX<T_active>::Zero(edge_contract_order+1);
        Eigen::VectorX<T_active> self_edge_proj_vals = Eigen::VectorX<T_active>::Zero(edge_contract_order+1);
       
        Eigen::VectorXd e1 = tet_facet_basis.row(0);
        Eigen::VectorXd e2 = tet_facet_basis.row(1);

        // std::cout << "n_idx: " << n_idx << std::endl;
        // std::cout << "e1: " << e1.transpose() << std::endl;
        // std::cout << "e2: " << e2.transpose() << std::endl;


        for (int v_i = 0; v_i < nprimals; v_i++)
        {
            Eigen::VectorX<T_active> cur_neighbor_vec = data.primals_rank1[v_i];
            Eigen::VectorX<T_active> cur_self_vec = self_data.primals_rank1[v_i];

            for (int term = 0; term < edge_contract_order+1; term++)
            {
                T_active neighbor_e1_contract = cur_neighbor_vec.dot(e1);// + 1e-12;
                T_active neighbor_e2_contract = cur_neighbor_vec.dot(e2);// + 1e-12;
                T_active self_e1_contract = cur_self_vec.dot(e1);// + 1e-12;
                T_active self_e2_contract = cur_self_vec.dot(e2);// + 1e-12;

                // for (int i = 0; i<3; i++)
                // {
                //     neighbor_e1_contract += cur_neighbor_vec(i)*e1(i);
                //     neighbor_e2_contract += cur_neighbor_vec(i)*e2(i);
                //     self_e1_contract += cur_self_vec(i)*e1(i);
                //     self_e2_contract += cur_self_vec(i)*e2(i);
                // }

                // T_active pow_n_e1 = 1;
                // T_active pow_n_e2 = 1;
                // T_active pow_s_e1 = 1;
                // T_active pow_s_e2 = 1;

                // for (int i = 0; i < term; i++)
                // {
                //     pow_n_e1 *= neighbor_e1_contract;
                //     pow_s_e1 *= self_e1_contract;
                // }
                // for (int i = 0; i < edge_contract_order-term; i++)
                // {
                //     pow_n_e2 *= neighbor_e2_contract;
                //     pow_s_e2 *= self_e2_contract;
                // }

                // neighbor_edge_proj_vals(term) += pow_n_e1 * pow_n_e2;
                // self_edge_proj_vals(term) += pow_s_e1 * pow_s_e2;

                neighbor_edge_proj_vals(term) += pow(neighbor_e1_contract, term)*
                                                 pow(neighbor_e2_contract, edge_contract_order-term);
                self_edge_proj_vals(term) += pow(self_e1_contract, term)*
                                             pow(self_e2_contract, edge_contract_order-term);
            }
        
        }

        data.Lk_edge_contract_diff = neighbor_edge_proj_vals - self_edge_proj_vals;



        for (int term = 0; term < edge_contract_order+1; term++)
        {
            data.Lk_edge_contract_diff(term) = data.Lk_edge_contract_diff(term);
        }

        
        // std::cout << "Lk_edge_contract_diff: " << data.Lk_edge_contract_diff.transpose() << std::endl;

    }

// This function rotates a given neighbor and element l2 tensors to their shared facet 
// using a precomputed rotation and saves their difference.  
    void L2_facet_diff(AppState& appState, const int n_idx, ElementData<T_active>& data) { 
        
        // Seperate out the primal dofs into rank-1 components
        // Maybe make this a seperate function. 
        int nprimals = data.primals_rank1.size();
        int facet_dim = data.primals_rank1[0].size() - 1;
        data.L_2_facet_diff.resize(facet_dim*facet_dim);
        data.L_2_facet_diff.setZero();
        
        Eigen::Matrix3d R_to_template = appState.R_facet_to_template.at(t_idx).at(n_idx);

        Eigen::VectorX<T_active> neighbor_L2_facet = Eigen::VectorX<T_active>::Zero(facet_dim*facet_dim);
        Eigen::VectorX<T_active> self_L2_facet = Eigen::VectorX<T_active>::Zero(facet_dim*facet_dim);
       
        

        for (int v_i = 0; v_i < nprimals; v_i++)
        {
            Eigen::VectorX<T_active> cur_neighbor_vec = data.primals_rank1[v_i];
            Eigen::VectorX<T_active> cur_self_vec = self_data.primals_rank1[v_i];

            Eigen::VectorX<T_active> rot_neighbor = R_to_template * cur_neighbor_vec;
            Eigen::VectorX<T_active> rot_self = R_to_template * cur_self_vec;

            // std::cout << rot_neighbor.transpose() << " | " << rot_self.transpose() << std::endl;


            // Eigen::MatrixX<T_active> curcurt = cur*cur.transpose();
            // Eigen::VectorX<T_active> curcurt_normalized_flattened;
            // curcurt_normalized_flattened.resize(cur.rows()*cur.rows());// = Eigen::Zeros(cur.rows()*cur.rows()); // TODO: make this compressed 
            // Eigen::VectorX<T_active> curcurt_normalized_flattened = curcurt_flattened;
            // TODO fix this later
            for (int i = 0; i < facet_dim; i++)
            {
                for (int j = 0; j < facet_dim; j++)
                {
                    // std::cout << "i: " << i << " j: " << j << std::endl;

                    neighbor_L2_facet(i*facet_dim + j) += rot_neighbor(i)*rot_neighbor(j);
                    self_L2_facet(i*facet_dim + j) += rot_self(i)*rot_self(j);
                }
            }
        
        }

        data.L_2_facet_diff = neighbor_L2_facet - self_L2_facet;

    }


    void L4_sym_krushkal(AppState& appState, ElementData<T_active>& data) { 
        
        // Seperate out the primal dofs into rank-1 components
        // Maybe make this a seperate function. 
        int nprimals = data.primals_rank1.size();
        int primals_size = data.primals_rank1[0].size();
        int L4_size = 15;
        data.L_4_sym_krushkal.resize(L4_size);
        data.L_4_sym_krushkal.setZero();
        // data.L_4_krushkal = data.L_4_primals;


        // update for n vectors per frame.  
        for (int v_i = 0; v_i < nprimals; v_i++)
        {
            Eigen::VectorX<T_active> cur = data.primals_rank1[v_i];
            
            T_active cur_norm = cur.norm() + 1e-10;
            Eigen::VectorX<T_active> cur_normalized = cur / cur_norm;


            // Eigen::MatrixX<T_active> curcurt = cur*cur.transpose();
            // Eigen::VectorX<T_active> curcurt_flattened;
            // curcurt_flattened.resize(L4_size);// = Eigen::Zeros(cur.rows()*cur.rows()); // TODO: make this compressed 
            Eigen::VectorX<T_active> curcurt_normalized_flattened; // = curcurt_flattened;
            curcurt_normalized_flattened.resize(L4_size);
            
            int cur_idx = 0;
            for (int i = 0; i < primals_size; i++)
            {
                for (int j = i; j < primals_size; j++)
                {
                    for (int k = j; k < primals_size; k++)
                    {
                        for (int l = k; l < primals_size; l++)
                        {
                            // std::cout << "i: " << i << " j: " << j << " k: " << k << " l: " << l << std::endl;
                            // curcurt_flattened(i*psize3 + j*psize2 + k*psize + l) = cur(i)*cur(j)*cur(k)*cur(l);
                            curcurt_normalized_flattened(cur_idx)  = cur_normalized(i)*cur_normalized(j)*cur_normalized(k)*cur_normalized(l); 
                            cur_idx++;
                        }
                    }
                }
            }

            // data.L_4_primals = data.L_4_primals + curcurt_flattened;
            data.L_4_sym_krushkal = data.L_4_sym_krushkal + curcurt_normalized_flattened * cur_norm;
        
        }

        // This sets the squared values which scales in the same way as L_4 to make the energy scale invariant. 
        // L2_to_L2x2(data);

    }
    
    
    
    
    void L4_krushkal(AppState& appState, ElementData<T_active>& data) { 
        
        // Seperate out the primal dofs into rank-1 components
        // Maybe make this a seperate function. 
        int nprimals = data.primals_rank1.size();
        int primals_size = data.primals_rank1[0].size();
        int L4_size = primals_size*primals_size*primals_size*primals_size;
        data.L_4_krushkal.resize(primals_size*primals_size*primals_size*primals_size);
        data.L_4_krushkal.setZero();
        // data.L_4_krushkal = data.L_4_primals;


        // update for n vectors per frame.  
        for (int v_i = 0; v_i < nprimals; v_i++)
        {
            Eigen::VectorX<T_active> cur = data.primals_rank1[v_i];
            
            T_active cur_norm = cur.norm() + 1e-10;
            Eigen::VectorX<T_active> cur_normalized = cur / cur_norm;


            // Eigen::MatrixX<T_active> curcurt = cur*cur.transpose();
            // Eigen::VectorX<T_active> curcurt_flattened;
            // curcurt_flattened.resize(L4_size);// = Eigen::Zeros(cur.rows()*cur.rows()); // TODO: make this compressed 
            Eigen::VectorX<T_active> curcurt_normalized_flattened; // = curcurt_flattened;
            curcurt_normalized_flattened.resize(L4_size);
            
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
                            // curcurt_flattened(i*psize3 + j*psize2 + k*psize + l) = cur(i)*cur(j)*cur(k)*cur(l);
                            curcurt_normalized_flattened(i*psize3 + j*psize2 + k*psize + l)  = cur_normalized(i)*cur_normalized(j)*cur_normalized(k)*cur_normalized(l); 
                        }
                    }
                }
            }

            // data.L_4_primals = data.L_4_primals + curcurt_flattened;
            data.L_4_krushkal = data.L_4_krushkal + curcurt_normalized_flattened * cur_norm;
        
        }

        // This sets the squared values which scales in the same way as L_4 to make the energy scale invariant. 
        // L2_to_L2x2(data);

    }


    // This function rotates a given neighbor and element l2 tensors to their shared facet 
// using a precomputed rotation and saves their difference.  
    void L4_facet_diff(AppState& appState, const int n_idx, ElementData<T_active>& data) { 
        
        // Seperate out the primal dofs into rank-1 components
        // Maybe make this a seperate function. 
        int nprimals = data.primals_rank1.size();
        int facet_dim = data.primals_rank1[0].size() - 1;
        data.L_4_facet_diff.resize(facet_dim*facet_dim*facet_dim*facet_dim);
        data.L_4_facet_diff.setZero();
        
        Eigen::Matrix3d R_to_template = appState.R_facet_to_template.at(t_idx).at(n_idx);

        Eigen::VectorX<T_active> neighbor_L4_facet = Eigen::VectorX<T_active>::Zero(facet_dim*facet_dim*facet_dim*facet_dim);
        Eigen::VectorX<T_active> self_L4_facet = Eigen::VectorX<T_active>::Zero(facet_dim*facet_dim*facet_dim*facet_dim);
       
        

        for (int v_i = 0; v_i < nprimals; v_i++)
        {
            Eigen::VectorX<T_active> cur_neighbor_vec = data.primals_rank1[v_i];
            Eigen::VectorX<T_active> cur_self_vec = self_data.primals_rank1[v_i];

            Eigen::VectorX<T_active> rot_neighbor = R_to_template * cur_neighbor_vec;
            Eigen::VectorX<T_active> rot_self = R_to_template * cur_self_vec;


            // TODO fix this later
            int fsize = facet_dim;
            int fsize2 = fsize*fsize;
            int fsize3 = fsize2*fsize;
            for (int i = 0; i < facet_dim; i++)
            {
                for (int j = 0; j < facet_dim; j++)
                {
                    for (int k = 0; k < facet_dim; k++)
                    {
                        for (int l = 0; l < facet_dim; l++)
                        {
                            // curcurt_flattened(i*psize3 + j*psize2 + k*psize + l) = cur(i)*cur(j)*cur(k)*cur(l);
                            neighbor_L4_facet(i*fsize3 + j*fsize2 + k*fsize + l) += rot_neighbor(i)*rot_neighbor(j)*rot_neighbor(k)*rot_neighbor(l);
                            self_L4_facet(i*fsize3 + j*fsize2 + k*fsize + l) += rot_self(i)*rot_self(j)*rot_self(k)*rot_self(l);
                            
                        }
                    }
                }
            }

        }

        data.L_4_facet_diff = neighbor_L4_facet - self_L4_facet;

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
