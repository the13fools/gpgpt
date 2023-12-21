#include "StripePatternIntegration.h"

#include <igl/doublearea.h>
#include <igl/cotmatrix_entries.h>

#include <SymGEigsShiftSolver.h>
#include <MatOp/SparseSymShiftSolve.h>



namespace SurfaceFields {
void StripePatternsGlobalIntegration::globallyIntegrateOneComponent(const Surface& surf, const Eigen::MatrixXd& v, Eigen::VectorXd& theta, Eigen::VectorXd* edge_omega) {
  size_t nedges = surf.data().E.rows();

  // Step 0: prepare the neccessary data
  // face mass matrix
  Eigen::VectorXd face_areas;
  igl::doublearea(surf.data().V, surf.data().F, face_areas);
  face_areas *= 0.5;

  // build edge metric matrix and inverse (cotan weights)
  Eigen::MatrixXd C;
  igl::cotmatrix_entries(surf.data().V, surf.data().F, C);
  Eigen::VectorXd edge_metric = Eigen::VectorXd::Zero(nedges);
  size_t nfaces = surf.nFaces();
  for (int i = 0; i < nfaces; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      int eidx = surf.data().faceEdges(i, j);
      edge_metric[eidx] += C(i, j);
    }
  }

  // Step 1: We need to convert the face vector fields to edge one forms
  Eigen::VectorXd edge_one_forms = GetEdgeOneForms(surf, v, &face_areas);

  // Step 2: form the problem. For each edge,
  Eigen::SparseMatrix<double> A = BuildConsistenceMatrix(surf, edge_one_forms, edge_metric);
  Eigen::SparseMatrix<double> B = BuildVertLumpedMassMatrix(surf, face_areas);

  // Step 3: solve the Eigen value problem: min 0.5 * x^T A * x, s.t. ||x||_B = 1
  Spectra::SymShiftInvert<double> op(A, B);
  Spectra::SparseSymMatProd<double> Bop(B);
  const double eps = 1e-10;
  Spectra::SymGEigsShiftSolver<Spectra::SymShiftInvert<double>, Spectra::SparseSymMatProd<double>, Spectra::GEigsMode::ShiftInvert> geigs(op, Bop, 1, 6, -2 * eps);

  geigs.init();
  int nconv = geigs.compute(Spectra::SortRule::LargestMagn, 1e6);

  Eigen::VectorXd evalues;
  Eigen::MatrixXd evecs;

  evalues = geigs.eigenvalues();
  evecs = geigs.eigenvectors();
  if (nconv != 1 || geigs.info() != Spectra::CompInfo::Successful) {
    throw "Eigensolver failed to converge!!";
  }

  std::cout << "Eigenvalue is " << evalues[0] << std::endl;

  size_t nverts = surf.data().V.rows();
  theta.setZero(nverts);
  for (int i = 0; i < nverts; i++)
  {
    std::complex<double> z = std::complex<double>(evecs(2 * i, 0), evecs(2 * i + 1, 0));
    theta[i] = std::arg(z);
  }

  if(edge_omega) {
    *edge_omega = std::move(edge_one_forms);
  }
}

// Get edge one forms from the input face vectors, where the face vectors will be rescaled by "scales", and averaged by face area
Eigen::VectorXd StripePatternsGlobalIntegration::GetEdgeOneForms(const Surface& surf, const Eigen::MatrixXd& v, const Eigen::VectorXd* face_areas) {
  size_t nedges = surf.data().E.rows();
  Eigen::VectorXd omega = Eigen::VectorXd::Zero(nedges);
  for (size_t eid = 0; eid < nedges; eid++) {
    Eigen::Vector3d e = (surf.data().V.row(surf.data().edgeVerts(eid, 1)) - surf.data().V.row(surf.data().edgeVerts(eid, 0))).transpose();
    double area_sum = 0;
    for (size_t id = 0; id < 2; id++) {
      int fid = surf.data().E(eid, id);
      if (fid == -1) {
        continue;
      }
      double area = face_areas ? (*face_areas)[fid] : 1;
      area_sum += area;
      Eigen::Vector3d Euclidean_vec = surf.data().Bs[fid] * (v.row(fid).transpose()); // The original vec
      omega[eid] += Euclidean_vec.dot(e) * area;
    }
    omega[eid] /= area_sum;
  }
  return omega;
}

// Build the consistnce matrix form the edge one forms. Section 4.1 in "Stripe Patterns on Surfaces"
Eigen::SparseMatrix<double> StripePatternsGlobalIntegration::BuildConsistenceMatrix(const Surface& surf, const Eigen::VectorXd& omega, const Eigen::VectorXd& edge_metric) {
  size_t nverts = surf.data().V.rows();
  size_t nedges = surf.data().E.rows();
  std::vector<Eigen::Triplet<double>> AT;
  for (size_t i = 0; i < nedges; i++) {
    int vid0 = surf.data().edgeVerts(i, 0);
    int vid1 = surf.data().edgeVerts(i, 1);

    std::complex<double> expw0 = std::complex<double>(std::cos(omega(i)), std::sin(omega(i)));


    AT.push_back({ 2 * vid0, 2 * vid0, edge_metric(i) });
    AT.push_back({ 2 * vid0 + 1, 2 * vid0 + 1, edge_metric(i) });

    AT.push_back({ 2 * vid1, 2 * vid1, edge_metric(i) });
    AT.push_back({ 2 * vid1 + 1, 2 * vid1 + 1, edge_metric(i) });


    AT.push_back({ 2 * vid0, 2 * vid1, -edge_metric(i) * (expw0.real()) });
    AT.push_back({ 2 * vid0 + 1, 2 * vid1, -edge_metric(i) * (-expw0.imag()) });
    AT.push_back({ 2 * vid0, 2 * vid1 + 1, -edge_metric(i) * (expw0.imag()) });
    AT.push_back({ 2 * vid0 + 1, 2 * vid1 + 1, -edge_metric(i) * (expw0.real()) });

    AT.push_back({ 2 * vid1, 2 * vid0, -edge_metric(i) * (expw0.real()) });
    AT.push_back({ 2 * vid1, 2 * vid0 + 1, -edge_metric(i) * (-expw0.imag()) });
    AT.push_back({ 2 * vid1 + 1, 2 * vid0, -edge_metric(i) * (expw0.imag()) });
    AT.push_back({ 2 * vid1 + 1, 2 * vid0 + 1, -edge_metric(i) * (expw0.real()) });
  }

  Eigen::SparseMatrix<double> A(2 * nverts, 2 * nverts);
  A.setFromTriplets(AT.begin(), AT.end());
  return A;
}

// Build vertex lumped mass matrix
Eigen::SparseMatrix<double> StripePatternsGlobalIntegration::BuildVertLumpedMassMatrix(const Surface& surf, const Eigen::VectorXd& face_area) {
  size_t nfaces = surf.data().F.rows();
  size_t nverts = surf.data().V.rows();

  std::vector<Eigen::Triplet<double>> BT;
  for (size_t fid = 0; fid < nfaces; fid++) {
    for (size_t i = 0; i < 3; i++) {
      int vid = surf.data().F(fid, i);
      BT.push_back(Eigen::Triplet<double>(2 * vid, 2 * vid, face_area[fid] / 3));
      BT.push_back(Eigen::Triplet<double>(2 * vid + 1, 2 * vid + 1, face_area[fid] / 3));
    }
  }

  Eigen::SparseMatrix<double>B(2 * nverts, 2 * nverts);
  B.setFromTriplets(BT.begin(), BT.end());

  return B;
}
}