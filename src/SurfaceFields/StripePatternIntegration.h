#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "../Surface.h"


namespace SurfaceFields {
// Stripe patterns on the surface
class StripePatternsGlobalIntegration
{
public:
  void globallyIntegrateOneComponent(const Surface& surf, const Eigen::MatrixXd& v, Eigen::VectorXd& theta, Eigen::VectorXd* edge_omega = nullptr);

private:
  // Get edge one forms from the input face vectors, where the face vectors will be rescaled by "scales", and averaged by face area (if provided)
  Eigen::VectorXd GetEdgeOneForms(const Surface&surf, const Eigen::MatrixXd& v, const Eigen::VectorXd* face_area = nullptr);

  // Build the consistnce matrix form the edge one forms. Section 4.1 in "Stripe Patterns on Surfaces"
  Eigen::SparseMatrix<double> BuildConsistenceMatrix(const Surface& surf, const Eigen::VectorXd& omega, const Eigen::VectorXd& edge_metric);

  // Build vertex lumped mass matrix
  Eigen::SparseMatrix<double> BuildVertLumpedMassMatrix(const Surface& surf, const Eigen::VectorXd& face_area);
};
}