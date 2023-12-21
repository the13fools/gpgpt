#include "UpsampleSurfaceFields.h"

#include <Eigen/Sparse>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/adjacency_list.h>

namespace SurfaceFields {

static double lArg(const long &n, const Eigen::Vector3d &bary)
{
  double larg = 0;
  double ti = bary(0), tj = bary(1), tk = bary(2);
  if(tk <= ti && tk <= tj)
    larg = M_PI / 3 * n * (1 + (tj - ti) / (1 - 3 * tk));
  else if (ti <= tj && ti <= tk)
    larg = M_PI / 3 * n * (3 + (tk - tj) / (1 - 3 * ti));
  else
    larg = M_PI / 3 * n * (5 + (ti - tk) / (1 - 3 * tj));
  return larg;
}

static void midPoint(const int n_verts, const Eigen::MatrixXi& F, Eigen::SparseMatrix<double>& S, Eigen::MatrixXi& NF, std::vector<int>* faceTracing)
{

  typedef Eigen::SparseMatrix<double> SparseMat;
  typedef Eigen::Triplet<double> Triplet_t;

  //Ref. igl::loop
  Eigen::MatrixXi FF, FFi;
  igl::triangle_triangle_adjacency(F, FF, FFi);
  std::vector<std::vector<typename Eigen::MatrixXi::Scalar>> adjacencyList;
  igl::adjacency_list(F, adjacencyList, true);
  //Compute the number and positions of the vertices to insert (on edges)
  Eigen::MatrixXi NI = Eigen::MatrixXi::Constant(FF.rows(), FF.cols(), -1);
  Eigen::MatrixXi NIdoubles = Eigen::MatrixXi::Zero(FF.rows(), FF.cols());
  Eigen::VectorXi vertIsOnBdry = Eigen::VectorXi::Zero(n_verts);
  int counter = 0;
  for (int i = 0; i < FF.rows(); ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      if (NI(i, j) == -1)
      {
        NI(i, j) = counter;
        NIdoubles(i, j) = 0;
        if (FF(i, j) != -1)
        {
          //If it is not a boundary
          NI(FF(i, j), FFi(i, j)) = counter;
          NIdoubles(i, j) = 1;
        }
        else
        {
          //Mark boundary vertices for later
          vertIsOnBdry(F(i, j)) = 1;
          vertIsOnBdry(F(i, (j + 1) % 3)) = 1;
        }
        ++counter;
      }
    }
  }

  const int& n_odd = n_verts;
  const int& n_even = counter;
  const int n_newverts = n_odd + n_even;

  //Construct vertex positions
  std::vector<Triplet_t> tripletList;
  for (int i = 0; i < n_odd; ++i)
  {
    //Old vertices
    tripletList.emplace_back(i, i, 1.);
  }
  for (int i = 0; i < FF.rows(); ++i)
  {
    //New vertices
    for (int j = 0; j < 3; ++j)
    {
      if (NIdoubles(i, j) == 0)
      {
        if (FF(i, j) == -1)
        {
          //Boundary vertex
          tripletList.emplace_back(NI(i, j) + n_odd, F(i, j), 1. / 2.);
          tripletList.emplace_back(NI(i, j) + n_odd, F(i, (j + 1) % 3), 1. / 2.);
        }
        else
        {
          tripletList.emplace_back(NI(i, j) + n_odd, F(i, j), 1. / 2.);
          tripletList.emplace_back(NI(i, j) + n_odd, F(i, (j + 1) % 3), 1. / 2.);
        }
      }
    }
  }
  S.resize(n_newverts, n_verts);
  S.setFromTriplets(tripletList.begin(), tripletList.end());

  // Build the new topology (Every face is replaced by four)
  if (faceTracing)
    faceTracing->resize(F.rows() * 4);
  NF.resize(F.rows() * 4, 3);

  for (int i = 0; i < F.rows(); ++i)
  {
    Eigen::VectorXi VI(6);
    VI << F(i, 0), F(i, 1), F(i, 2), NI(i, 0) + n_odd, NI(i, 1) + n_odd, NI(i, 2) + n_odd;

    Eigen::VectorXi f0(3), f1(3), f2(3), f3(3);
    f0 << VI(0), VI(3), VI(5);
    f1 << VI(1), VI(4), VI(3);
    f2 << VI(2), VI(5), VI(4);
    f3 << VI(3), VI(4), VI(5);

    NF.row((i * 4) + 0) = f0;
    NF.row((i * 4) + 1) = f1;
    NF.row((i * 4) + 2) = f2;
    NF.row((i * 4) + 3) = f3;

    if (faceTracing) {
      for (int j = 0; j < 4; j++) {
        (*faceTracing)[4 * i + j] = i;
      }
    }

  }
}


void upsampleSurface(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int up_level, Eigen::MatrixXd& NV, Eigen::MatrixXi& NF, std::vector<std::pair<int, Eigen::Vector3d>>* bary) {
  NV = V;
  NF = F;
  // midpoint subdivision
  std::vector<int> tmp_facemap;
  Eigen::SparseMatrix<double> tmp_mat(V.rows(), V.rows());
  tmp_mat.setIdentity();

  tmp_facemap.resize(F.rows());
  for (int i = 0; i < F.rows(); i++) {
    tmp_facemap[i] = i;
  }

  for(int i = 0; i < up_level; ++i) {
    Eigen::MatrixXi tempF = NF;
    Eigen::SparseMatrix<double> S;
    std::vector<int> faceTracing;
    midPoint(NV.rows(), tempF, S, NF, &faceTracing);
    // This .eval is super important
    NV = (S*NV).eval();

    tmp_mat = S * tmp_mat;

    std::vector<int> tmp_facemap1;
    tmp_facemap1.resize(NF.rows());
    for (int j = 0; j < NF.rows(); j++) {
      tmp_facemap1[j] = tmp_facemap[faceTracing[j]];
    }
    tmp_facemap = std::move(tmp_facemap1);
  }

  if(bary) {
    bary->resize(NV.rows());
    Eigen::VectorXi isVisited = Eigen::VectorXi::Zero(NV.rows());

    std::vector<std::vector<std::pair<int, double>>> nonzeroTracing(NV.rows());

    for (int k=0; k < tmp_mat.outerSize(); ++k) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(tmp_mat,k); it; ++it) {
        nonzeroTracing[it.row()].push_back({it.col(), it.value()});
      }
    }

    for(int i = 0; i < NF.rows(); i++) {
      for(int j = 0; j < 3; j++) {
        int vid = NF(i, j);
        if(isVisited(vid))
          continue;
        std::vector<std::pair<int, double>> perFaceBary = nonzeroTracing[vid];
        if(perFaceBary.size() >3) {
          std::cerr << "some error in the upsampling matrix." << std::endl;
          exit(1);
        }
        int preFace = tmp_facemap[i];
        Eigen::Vector3d ptBary = Eigen::Vector3d::Zero();

        for(int k = 0; k < perFaceBary.size(); k++) {
          int preVid = perFaceBary[k].first;
          for(int n = 0; n < 3; n++) {
            if(F(preFace, n) == preVid) {
              ptBary(n) = perFaceBary[k].second;
            }
          }
        }
        bary->at(vid) = std::pair<int, Eigen::Vector3d>(preFace, ptBary);
        isVisited(vid) = 1;
      }
    }
  }
}

void upsampleScalarFields(const Surface& surf, const Eigen::VectorXd& edge_one_form, const Eigen::VectorXd& theta, const std::vector<std::pair<int, Eigen::Vector3d>>& bary, Eigen::VectorXd& up_theta) {
  int upsize = bary.size();
  up_theta.setZero(upsize);

  for(int i = 0; i < bary.size(); i++)
  {
    int fid = bary[i].first;
    double omegaIJ = edge_one_form(surf.data().faceEdges(fid, 2));
    double omegaJK = edge_one_form(surf.data().faceEdges(fid, 0));
    double omegaKI = edge_one_form(surf.data().faceEdges(fid, 1));

    double cIJ = surf.data().F(fid, 0) == surf.data().edgeVerts(surf.data().faceEdges(fid, 2), 0) ? 1 : -1;
    double cJK = surf.data().F(fid, 1) == surf.data().edgeVerts(surf.data().faceEdges(fid, 0), 0) ? 1 : -1;
    double cKI = surf.data().F(fid, 2) == surf.data().edgeVerts(surf.data().faceEdges(fid, 1), 0) ? 1 : -1;

    omegaIJ *= cIJ;
    omegaJK *= cJK;
    omegaKI *= cKI;

    std::complex<double> rij( std::cos(omegaIJ), std::sin(omegaIJ) );
    std::complex<double> rjk( std::cos(omegaJK), std::sin(omegaJK) );
    std::complex<double> rki( std::cos(omegaKI), std::sin(omegaKI) );

    std::complex<double> psiI = std::complex<double>(std::cos(theta[surf.data().F(fid, 0)]), std::sin(theta[surf.data().F(fid, 0)]));
    std::complex<double> psiJ = std::complex<double>(std::cos(theta[surf.data().F(fid, 1)]), std::sin(theta[surf.data().F(fid, 1)]));
    std::complex<double> psiK = std::complex<double>(std::cos(theta[surf.data().F(fid, 2)]), std::sin(theta[surf.data().F(fid, 2)]));


    double alphaI = std::arg(psiI);
    double alphaJ = alphaI + omegaIJ - std::arg(rij*psiI/psiJ); //fmodPI((varphiI + omegaIJ) - varphiJ); // could do this in terms of angles instead of complex numbers...
    double alphaK = alphaJ + omegaJK - std::arg(rjk*psiJ/psiK); //fmodPI((varphiJ + omegaJK) - varphiK); // mostly a matter of taste---possibly a matter of performance?
    double alphaL = alphaK + omegaKI - std::arg(rki*psiK/psiI); //fmodPI((varphiK + omegaKI) - varphiI);

    // adjust triangles containing zeros
    long n = std::lround((alphaL-alphaI)/(2.*M_PI));
    alphaJ -= 2.*M_PI*n/3.;
    alphaK -= 4.*M_PI*n/3.;

    double theta = lArg(n, bary[i].second);
    up_theta(i) = theta + bary[i].second(0) * alphaI + bary[i].second(1) * alphaJ + bary[i].second(2) * alphaK;
  }
}
}