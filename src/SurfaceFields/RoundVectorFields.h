#pragma once
#include <Eigen/Dense>
#include <vector>
#include "GlobalFieldIntegration.h"

// Round vector fields into a scalar field. 
// mesh_pts:			the vertex positions of the input mesh 
// mesh_faces:			the face connectivity of the input mesh
// vecs:			the vector fields on the mesh. Defined per face. vecs.rows() = k * nfaces, where the first k rows is for the vector fields on the first face, and so on.
//				Moreover, the vector fields is stored intrinsically. Assume for face fid, the basis is b0 and b1 \in R^3, then the extrinsic vector (for j-th vector) is vecs(fid * k + j, 0) * b0 + vecs(fid * k + j, 1) * b1
// splitted_pts:		the vertex positions of the output mesh, which contains 2 * k copies of the input mesh
// splitted_faces:		the face connvectivity of the output mesh
// round_op:			the round operation.
// theta:		        the theta value whose gradient is equal to vecs on the splitted mesh
// grad_theta;                  the gradient of theta
// splitted_vecs:		the extrinsic vector fields on the splitted mesh
// global_rescaling:	        the global rescaling ratio
// cut_pts:			the cutted pts. RM: we do NOT actually cut the splitted meshes
// cut_edges:			the connectivity of the cut points
void roundVectorFields(const Eigen::MatrixXd& mesh_pts, const Eigen::MatrixXi& mesh_faces, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& splitted_pts, Eigen::MatrixXi& splitted_faces, Eigen::VectorXd& theta, SurfaceFields::GlobalFieldIntegration* round_op, Eigen::MatrixXd* grad_theta = nullptr, double global_rescaling = 1, Eigen::MatrixXd* splitted_vecs = nullptr, std::vector<Eigen::Vector3d>* cut_pts = nullptr, std::vector<std::vector<int>>* cut_edges = nullptr);