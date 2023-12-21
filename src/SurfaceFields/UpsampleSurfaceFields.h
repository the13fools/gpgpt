#pragma once
#include "../Surface.h"

namespace SurfaceFields {
  void upsampleSurface(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int up_level, Eigen::MatrixXd& up_V, Eigen::MatrixXi& up_F, std::vector<std::pair<int, Eigen::Vector3d>>* bary = nullptr);

  // The rendering technique mentioned in the paper "Stripe Patterns on Surfaces"
  void upsampleScalarFields(const Surface& surf, const Eigen::VectorXd& edge_one_form, const Eigen::VectorXd& theta, const std::vector<std::pair<int, Eigen::Vector3d>>& bary, Eigen::VectorXd& up_theta);
}