#pragma once

#include "../Surface.h"

namespace SurfaceFields {
	/*
	* Convert poly vector fields on the surface into corresponding vectors. Any vector field v, has been turned to v and w, where the angle between v and w is pi
	*/
	std::vector<Eigen::Vector2d> getVectors(Surface& s, int fid, const std::vector<Eigen::Vector2d>& poly_vecs);

	double vectorAngle(Surface& s, int face, const Eigen::Vector2d& v);


	// Comb poly vectors. Assume the input polyvectors are stored as vector (barycentrically) fields. More specifically, 
	// assume the vec_fields.rows() = k * s.nFaces(). Then for each face, we will have k poly vectors. 
	// The return will be a 2k * nFaces() x 2 matrix, the rows i, i + k, i + 2 * k, ..., i + (nfaces - 1) * k are the 
	// corresponding combed vectors
	Eigen::MatrixXd CombPolyVectors(Surface& s, const Eigen::MatrixXd& vec_fields);
	
}