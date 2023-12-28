#ifndef ELEMENTVARSTINYAD_H
#define ELEMENTVARSTINYAD_H

#include <string>

#include <TinyAD/ScalarFunction.hh>

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
	primal, L2_tensor, L2_krushkal, L4_tensor, L4_krushkal, Element_COUNT
};

template <typename T_active>
class ElementData {
public:
	ElementData() : curr_idx(-1) {}
	ElementData(const Eigen::Index idx, const Eigen::VectorX<T_active>& dofs) : curr_idx(idx), dofs_curr_elem(dofs) {}

	Eigen::Index curr_idx;                                      // face idx
	
	Eigen::VectorX<T_active> dofs_curr_elem;                    // DOFs per element. For single vector field, its size is 2, and for cross fields, its size is 4
	std::vector< Eigen::VectorX<T_active> > primals_rank1;      // the vector of vector fields on the element, where each element is stored barycentrically  

	// These are the lifted primals to moment space.  
	// i.e. if primal_rank1, denoted as v = [u0, u1] barycentrically (v = u0 * b0 + v0 * b1, bi in R^{3x1})
	// then:
	// 1. L_2_primals(v) := \sum_i\sum_j ui * uj * bi \otimes bj (the definition of v \otimes v), where i, j takes the summation from 0 to 1.
	// Notice that bi \otimes bj is determined by the basis. We can encode this rank2 tensor as a 2x2 matrix (4x1 vector) = [u0 * u0, u0 * u1, u1 * u0, u1 * u1] 
	// Remark for any tensor T = \sum_i\sum_j Tij bi bj^T, if we treat bi \otimes bj as a 3x3 matrix bi bj^T := B_{ij}:
	//  ||T||^2 := FN(T)^2 = Tij * Tkl * Tr(B_{ij}^T B_{kl}) = Tij * Tkl * (bi^T bk) * (bj^T bl), where we use einstein summation notation and FN is the standard frobenius norm of a matrix
	// Tr(B_{ij}^T B_{kl}) = Tr(bj bi^T bk bl^T) = Tr(bi^T bk bl^T bj) = (bi^T bk) * (bl^T bj) (both expressions are scalar)
	// 
	// 1.1 L_2_krushkal(v) := |v| * L_2_primals(\hat(v))
	// 
	// 2. L_4_primals(v) := (v \otimes v) \otimes (v \otimes v) = ui * uj * uk * ul * (bi \otimes bj) \otimes (bk \otimes bl), where einstein summation notation again is applied.
	// Similar to the L_2 version, (bi \otimes bj) \otimes (bk \otimes bl) is all determined by the basis, we can encode this rank-4 tensor by a 16x1 vector L4. For each entry n = i * 2^3 + j * 2^2 + k * 2^1 + l * 2^0, L4[n] = ui * uj * uk * ul 
	 
	// Remark for any tenor T = T_{ijkl} (bi \otimes bj) \otimes (bk \otimes bl), we can treat (bi \otimes bj) \otimes (bk \otimes bl) as a 9x9 matrix.
	// More specifically:
	// (bi \otimes bj) \otimes (bk \otimes bl) = (bi bj^T) \otimes (bk bl^T) = B_{ijkl}, where the (m, n)-th 3x3 block of B_{ijkl} is (bi bj^T)[m, n] * bk bl^T.
	// In this way
	// ||T||^2 := FN(T)^2 = T_{ijkl} * T_{mnpq} * Tr(B_{ijkl}^T * B_{mnpq}) = T_{ijkl} * T_{mnpq} * Tr(B_{ij}^T B_{mn}) * Tr(B_{kl}^T B_{pq}) = T_{ijkl} * T_{mnpq} * (bi^T bm) * (bj^T bn) * (bk^T bp) * (bl^T bq) (with some calculation) 
	//
	// 2.1 L_4_krushkal(v) := |v| * L_4_primals(\hat(v))
	Eigen::VectorX<T_active> L_2_primals;
	Eigen::VectorX<T_active> L_4_primals;
	Eigen::VectorX<T_active> L_2_krushkal;
	Eigen::VectorX<T_active> L_4_krushkal;


	Eigen::Matrix<T_active, 3, 2> B;        // basis functions, where b0 = B.col(0) and b1 = B.col(1) 
	Eigen::Matrix<T_active, 2, 2> g;         // metric derived from the basis functions g = B^T B 
	Eigen::Matrix<T_active, 4, 1> Bij;      // If n = i * 2 + j, then Bij[n] = bi^T bj

	// set basis
	void SetBasis(const Eigen::Matrix<double, 3, 2>& basis) {
		B = basis;
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				int idx = (i << 1) | j;
				Bij[idx] = (B.col(i)).dot(B.col(j));
			}
		}
		g = B.transpose() * B;
	}

	// Seperate out the primal dofs into rank-1 components
	void SetPrimalsRank1(DOFMemoryLayout& primals_layout) {
		int primals_size = primals_layout.size;     // #dofs
		int primal_rank = primals_layout.rank;      // #vectors per elements

		int rank1_size = primals_size / primal_rank; // Elsewhere should check that this is an integer
		
		primals_rank1.resize(primal_rank);
		for (int v_i = 0; v_i < primal_rank; v_i++) {
			primals_rank1[v_i] = dofs_curr_elem.segment(primals_layout.start + v_i*rank1_size, rank1_size);
		}
	}

	int GetIndex(int i, int j) {
		// i, j \in {0, 1}
		return (i << 1) | j;  // 2 * i + j
	}

	int GetIndex(int i, int j, int k, int l) {
		// i, j, k, l \in {0, 1}
		return (i << 3) | (j << 2) | (k << 1) | l; // 8 * i + 4 * j + 2 * k + l;
	}

	// Set the L2 tensors
	void SetL2Primals() {
		assert(!primals_rank1.empty());

		// Maybe make this a seperate function for L2_primals and L_2_krushkal. 
		int nprimals = primals_rank1.size();
		int primals_size = primals_rank1[0].size();

		assert(primals_size == 2);      // we are handling 2d manifold surfaces

		L_2_primals.resize(primals_size * primals_size);
		L_2_primals.setZero();
		L_2_krushkal = L_2_primals;

		// if we have n vector perface, we define
		// L_2_primals(u^1, u^2, ..., u^n) = \sum_i\sum_j (u^1_i * u^1_j + ... + u^n_i * u^n_j)  +  * bi \otimes bj

		// update for n vectors per frame.  
		for (int v_i = 0; v_i < nprimals; v_i++) {
			Eigen::VectorX<T_active> cur = primals_rank1[v_i];
			Eigen::VectorX<T_active> cur_euclidean = cur[0] * B.col(0) + cur[1] * B.col(1);
			T_active cur_norm = cur_euclidean.norm() + 1e-10;
			Eigen::VectorX<T_active> cur_normalized = cur / cur_norm;

			Eigen::VectorX<T_active> curcurt_flattened;
			curcurt_flattened.setZero(cur.rows() * cur.rows());

			Eigen::VectorX<T_active> curcurt_normalized_flattened = curcurt_flattened;
			for (int i = 0; i < primals_size; i++) {
				for (int j = 0; j < primals_size; j++) {
					int idx = GetIndex(i, j);
					curcurt_flattened(idx) = cur(i) * cur(j);
					curcurt_normalized_flattened(idx) = cur_normalized(i) * cur_normalized(j);
				}
			}

			L_2_primals = L_2_primals + curcurt_flattened;
			L_2_krushkal = L_2_krushkal + curcurt_normalized_flattened * cur_norm;
		}
	}

	// Set L4 tensors
	void SetL4Primals() {
		assert(!primals_rank1.empty());

		// Maybe make this a seperate function. 
		int nprimals = primals_rank1.size();
		int primals_size = primals_rank1[0].size();

		assert(primals_size == 2);      // we are handling 2d manifold surfaces

		int L4_size = primals_size * primals_size * primals_size * primals_size;
		L_4_primals.resize(primals_size * primals_size * primals_size * primals_size);
		L_4_primals.setZero();
		L_4_krushkal = L_4_primals;


		// update for n vectors per frame.  
		for (int v_i = 0; v_i < nprimals; v_i++) {
			Eigen::VectorX<T_active> cur = primals_rank1[v_i];
			Eigen::VectorX<T_active> cur_euclidean = cur[0] * B.col(0) + cur[1] * B.col(1);
			T_active cur_norm = cur_euclidean.norm() + 1e-10;
			Eigen::VectorX<T_active> cur_normalized = cur / cur_norm;

			Eigen::VectorX<T_active> curcurt_flattened;
			curcurt_flattened.setZero(L4_size);
			Eigen::VectorX<T_active> curcurt_normalized_flattened = curcurt_flattened;

			for (int i = 0; i < primals_size; i++) {
				for (int j = 0; j < primals_size; j++) {
					for (int k = 0; k < primals_size; k++) {
						for (int l = 0; l < primals_size; l++) {
							int idx = GetIndex(i, j, k, l);
							curcurt_flattened(idx) = cur(i) * cur(j) * cur(k) * cur(l);
							curcurt_normalized_flattened(idx) = cur_normalized(i) * cur_normalized(j) * cur_normalized(k) * cur_normalized(l);
						}
					}
				}
			}

			L_4_primals = L_4_primals + curcurt_flattened;
			L_4_krushkal = L_4_krushkal + curcurt_normalized_flattened * cur_norm;
		}
	}
};



/**
 * ElementVars class 
 */
template <typename T_active, typename ELEM>
class ProcElement {
public:
	ProcElement(ElementLiftType t) : w_bound(0), curr_lift(t) {};

	// ElementVars(AppState& appState) {}
	ElementLiftType curr_lift; 
	double w_bound;
	double w_smooth_vector;

	double w_smooth;
	double w_curl;
	double w_attenuate;

	int num_neighbors = 0;

	Surface* cur_surf;

	Eigen::Index f_idx; // maybe get rid of this it's redundant. 

	ElementData<T_active> self_data;		// current element data

	std::vector<ElementData<T_active>> neighbor_data;	// neighboring element data, where we convert the barycentric vector fields into the current basis
	std::vector<Eigen::Index> shared_eids;	// the shared edge ids with neighboring faces 


	void SetSelfData(AppState& appState, const Eigen::Index f_idx, ELEM& element) { 
		SetElemState(appState, f_idx);

		self_data.curr_idx = f_idx;
		self_data.dofs_curr_elem = element.variables(f_idx);
		self_data.SetPrimalsRank1(appState.primals_layout);
		self_data.SetBasis(appState.cur_surf->data().Bs[f_idx]);
		
		switch(curr_lift) {
			case(ElementLiftType::primal):
                break;
			// Can seperate these two out to make it more granular if it's necessary for a bit of a speed boost
			case(ElementLiftType::L2_krushkal):
			case(ElementLiftType::L2_tensor):
				self_data.SetL2Primals();
			break;

			case(ElementLiftType::L4_krushkal):
			case(ElementLiftType::L4_tensor):
				self_data.SetL4Primals();
			break;
			
			default:
                std::cout << "Error: LiftType not implemented" << std::endl;

		}
	}

	void SetNeighborData(AppState& appState, const Eigen::Index f_idx, ELEM& element) { 
		neighbor_data.clear();
		shared_eids.clear();

		num_neighbors = 0;

		auto& E = cur_surf->data().E;
		auto& Ts = cur_surf->data().Ts;

		for (int i = 0; i < 3; i++) {
			int eid = cur_surf->data().faceEdges(f_idx, i);

			if (E(eid, 0) == -1 || E(eid, 1) == -1) {
				continue;
			} // Make sure this is the correct convention.  Do more advance boundary handling later.
			
			int efid = E(eid, 0) == f_idx ? 0 : 1;
			
			assert(E(eid, efid) == f_idx);
			
			int n_idx = E(eid, 1 - efid);	// neighboring face id
			Eigen::Matrix<T_active, 2, 2> T = Ts.block<2, 2>(2 * i, 2 * (1 - efid)); // the transition matrix maps vectors from barycentric coordinates of face E(i, 1 - efid) (= n_idx) to barycentric coordinates of E(i, efid) (= f_idx)

			ElementData<T_active> neighbor_data_i;
			num_neighbors += 1;
			neighbor_data_i.curr_idx = n_idx;
			neighbor_data_i.dofs_curr_elem = element.variables(n_idx);

			int nvecs = neighbor_data_i.dofs_curr_elem.rows() / 2;
			
			// convert the vectors into the same basis
			for (int n = 0; n < nvecs; n++) {
				neighbor_data_i.dofs_curr_elem.template segment<2>(2 * n) = T * neighbor_data_i.dofs_curr_elem.template segment<2>(2 * n);
			}

			neighbor_data_i.SetPrimalsRank1(appState.primals_layout);
			neighbor_data_i.SetBasis(appState.cur_surf->data().Bs[f_idx]);

			switch(curr_lift) {
				case(ElementLiftType::primal):
                    break;
				// Can seperate these two out to make it more granular if it's necessary for a bit of a speed boost
				case(ElementLiftType::L2_krushkal):
				case(ElementLiftType::L2_tensor):
					neighbor_data_i.SetL2Primals();
				break;

				case(ElementLiftType::L4_krushkal):
				case(ElementLiftType::L4_tensor):
					neighbor_data_i.SetL4Primals();
				break;

				default:
                    std::cout << "Error: LiftType not implemented" << std::endl;
      
			}
			neighbor_data.push_back(neighbor_data_i);
			shared_eids.push_back(eid);
		}

	}

	void SetElemState(AppState& appState, const Eigen::Index f_idx) { 
		cur_surf = appState.cur_surf.get();

		w_bound = appState.config->w_bound;
		w_smooth_vector = appState.config->w_smooth_vector;

		w_smooth = appState.config->w_smooth;
		w_curl = appState.config->w_curl;
		w_attenuate = appState.config->w_attenuate;

		this->f_idx = f_idx;
	}

	T_active GetL2DirichletTerm() {
		T_active dirichlet_term = T_active(0);
		T_active face_area = cur_surf->faceArea(f_idx);

		auto get_face_centroid = [](Surface* surf, int face_id) {
			Eigen::VectorX<T_active> face_centroid = surf->data().V.row(surf->data().F(face_id, 0));
			face_centroid += surf->data().V.row(surf->data().F(face_id, 1));
			face_centroid += surf->data().V.row(surf->data().F(face_id, 2));
			face_centroid /= 3;
			return face_centroid;
			};

		Eigen::VectorX<T_active> c = get_face_centroid(cur_surf, f_idx);

		for (int nid = 0; nid < neighbor_data.size(); nid++) {
			int eid = shared_eids[nid];
			int nfid = neighbor_data[nid].curr_idx;
			int vid0 = cur_surf->data().edgeVerts(eid, 0);
			int vid1 = cur_surf->data().edgeVerts(eid, 1);

			Eigen::VectorX<T_active> nc = get_face_centroid(cur_surf, nfid);
			Eigen::VectorX<T_active> edge_midpt = (cur_surf->data().V.row(vid0) + cur_surf->data().V.row(vid1)) / 2;

			T_active barycentric_len = (c - edge_midpt).norm() + (nc - edge_midpt).norm();

			T_active nei_face_area = cur_surf->faceArea(nfid);
			T_active edge_area = (face_area + nei_face_area) / 2;

			Eigen::VectorX<T_active> L_2_krushkal_fd = (self_data.L_2_krushkal - neighbor_data[nid].L_2_krushkal) / barycentric_len;

			assert(L_2_krushkal_fd.rows() == 4);
			// ||T||^2 := FN(T)^2 = Tij * Tkl * Tr(B_{ij}^T B_{kl}) = Tij * Tkl * (bi^T bk) * (bj^T bl)

			T_active norm_sq = T_active(0);
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						for (int l = 0; l < 2; l++) {
							int idx_ij = self_data.GetIndex(i, j);
							int idx_kl = self_data.GetIndex(k, l);
							int idx_ik = self_data.GetIndex(i, k);
							int idx_jl = self_data.GetIndex(j, l);

							norm_sq += L_2_krushkal_fd(idx_ij) * L_2_krushkal_fd(idx_kl) * self_data.Bij[idx_ik] * self_data.Bij[idx_jl];
						}
					}
				}
			}

			dirichlet_term += norm_sq * edge_area;
		}
		return dirichlet_term;
	}

	T_active GetL4DirichletTerm() {
		T_active dirichlet_term = T_active(0);
		T_active face_area = cur_surf->faceArea(f_idx);

		auto get_face_centroid = [](Surface* surf, int face_id) {
			Eigen::VectorX<T_active> face_centroid = surf->data().V.row(surf->data().F(face_id, 0));
			face_centroid += surf->data().V.row(surf->data().F(face_id, 1));
			face_centroid += surf->data().V.row(surf->data().F(face_id, 2));
			face_centroid /= 3;
			return face_centroid;
			};

		Eigen::VectorX<T_active> c = get_face_centroid(cur_surf, f_idx);

		for (int nid = 0; nid < neighbor_data.size(); nid++) {
			int eid = shared_eids[nid];
			int nfid = neighbor_data[nid].curr_idx;
			int vid0 = cur_surf->data().edgeVerts(eid, 0);
			int vid1 = cur_surf->data().edgeVerts(eid, 1);

			Eigen::VectorX<T_active> nc = get_face_centroid(cur_surf, nfid);
			Eigen::VectorX<T_active> edge_midpt = (cur_surf->data().V.row(vid0) + cur_surf->data().V.row(vid1)) / 2;

			T_active barycentric_len = (c - edge_midpt).norm() + (nc - edge_midpt).norm();

			T_active nei_face_area = cur_surf->faceArea(nfid);
			T_active edge_area = (face_area + nei_face_area) / 2;

			Eigen::VectorX<T_active> L_4_krushkal_fd = (self_data.L_4_krushkal - neighbor_data[nid].L_4_krushkal) / barycentric_len;

			assert(L_4_krushkal_fd.rows() == 16);

			T_active norm_sq = T_active(0);

			//  ||T||^2 := FN(T)^2 = T_{ijkl} * T_{mnpq} * (bi^T bm) * (bj^T bn) * (bk^T bp) * (bl^T bq). This can be modified!
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						for (int l = 0; l < 2; l++) {
							int idx_ijkl = self_data.GetIndex(i, j, k, l);

							for (int m = 0; m < 2; m++) {
								for (int n = 0; n < 2; n++) {
									for (int p = 0; p < 2; p++) {
										for (int q = 0; q < 2; q++) {
											int idx_mnpq = self_data.GetIndex(m, n, p, q);

											int idx_im = self_data.GetIndex(i, m);
											int idx_jn = self_data.GetIndex(j, n);
											int idx_kp = self_data.GetIndex(k, p);
											int idx_lq = self_data.GetIndex(l, q);
											
											norm_sq += L_4_krushkal_fd(idx_ijkl) * L_4_krushkal_fd(idx_mnpq) * self_data.Bij[idx_im] * self_data.Bij[idx_jn] * self_data.Bij[idx_kp] * self_data.Bij[idx_lq];
										}
									}
								}
							}
						}
					}
				}
			}

			dirichlet_term += norm_sq * edge_area;
		}

		return dirichlet_term;
	}


/**
 * Convenience Struct for computing the element metadata for the objective functions.
 * 
 * @param AppData, type, element number  
 * @return ElementVars struct.  
 */
};

#endif // ELEMENTVARSTINYAD_H
