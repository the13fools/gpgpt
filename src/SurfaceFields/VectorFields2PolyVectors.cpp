#include "VectorFields2PolyVectors.h"

#include <Eigen/Geometry>
#include <deque>
#include <iostream>

namespace SurfaceFields {
double vectorAngle(Surface& s, int face, const Eigen::Vector2d& v)
{
  Eigen::Vector3d n = s.faceNormal(face);
  Eigen::Vector3d extv = s.data().Bs[face] * v;
  Eigen::Vector3d u = s.data().Bs[face].col(0);
  return 2 * std::atan2(u.cross(extv).dot(n), u.norm() * extv.norm() + u.dot(extv));
}

std::vector<Eigen::Vector2d> getVectors(Surface& s, int fid, const std::vector<Eigen::Vector2d>& poly_vecs) {
  int nvecs = poly_vecs.size();
  Eigen::Vector3d n = s.faceNormal(fid);
  Eigen::Matrix<double, 3, 2> B = s.data().Bs[fid];
  Eigen::Vector3d u = B.col(0);
  u.normalize();
  Eigen::Matrix2d BTBinv = (B.transpose() * B).inverse();

  std::vector<std::pair<double, Eigen::Vector2d>> vecs = {};
  for (int i = 0; i < nvecs; i++) {
    auto& v = poly_vecs[i];
    double theta = vectorAngle(s, fid, v); // theta should between -pi and pi
    if (theta < 0) {
      theta += 2 * M_PI;
    } // convert it into [0, 2 * pi]

    vecs.push_back({ theta, v });
    double angle = theta + M_PI; // angle between [pi, 3 * pi]

    if (angle >= 2 * M_PI) {
      angle = angle - 2 * M_PI;
    }

    Eigen::Vector2d oppv = BTBinv * B.transpose() * Eigen::AngleAxisd(angle, n).toRotationMatrix() * u;

    double angle_1 = vectorAngle(s, fid, oppv);

    oppv *= (v.norm() / oppv.norm());
    vecs.push_back({ angle, oppv });
  }

  std::sort(vecs.begin(), vecs.end(), [&](const std::pair<double, Eigen::Vector2d>& a, const std::pair<double, Eigen::Vector2d>& b) {
    return a.first < b.first;
  });

  std::vector<Eigen::Vector2d> sorted_vecs;
  for (auto& angle_vec : vecs) {
    sorted_vecs.push_back(angle_vec.second);
  }

  return sorted_vecs;
}

Eigen::MatrixXd CombPolyVectors(Surface& s, const Eigen::MatrixXd& vec_fields) {
  int nfaces = s.nFaces();
  std::vector<bool> visited(nfaces, false);
  assert((vec_fields.rows() % nfaces) == 0);
  int npoly = vec_fields.rows() / nfaces;

  Eigen::MatrixXd comb_vecs(2 * npoly * nfaces, 2);

  // BFS of faces
  for (int fid = 0; fid < nfaces; fid++) {
    if (visited[fid]) {
      continue;
    }
    struct Visit
    {
      int from;
      int to;
    };

    std::deque<Visit> q;
    q.push_back(Visit{ -1, fid });
    while (!q.empty()) {
      Visit vis = q.front();
      q.pop_front();
      if (visited[vis.to]) {
        continue;
      }
      visited[vis.to] = true;
      std::vector<Eigen::Vector2d> poly_vecs;
      for (int pid = 0; pid < npoly; pid++) {
        poly_vecs.push_back(vec_fields.row(npoly * vis.to + pid).transpose());
      }
      std::vector<Eigen::Vector2d> vecs = getVectors(s, vis.to, poly_vecs);

      if (vis.from == -1) {
        for (int vid = 0; vid < vecs.size(); vid++) {
          comb_vecs.row(2 * npoly * vis.to + vid) = vecs[vid].transpose();
        }
      }
      else {
        // since the vectors are sorted anticlockwisely, we only need to find the circular permutation
        // vector 0 on neighbor face
        Eigen::Vector2d v0 = comb_vecs.row(2 * npoly * vis.from).transpose();
        // tranport to current face
        int edge = -1;
        int fromside = -1;
        for (int j = 0; j < 3; j++) {
          int eid = s.data().faceEdges(vis.to, j);
          int opp = 0;
          if (s.data().E(eid, opp) == vis.to) {
            opp = 1;
          }

          if (s.data().E(eid, opp) == vis.from) {
            edge = eid;
            fromside = opp;
          }
        }
        assert(edge != -1);
        assert(s.data().E(edge, fromside) == vis.from);
        assert(s.data().E(edge, 1 - fromside) == vis.to);
        Eigen::Vector2d xportv0 = s.data().Ts.block<2, 2>(2 * edge, 2 * fromside) * v0;
        int bestidx = 0;
        double bestdot = -1.0;
        for (int j = 0; j < vecs.size(); j++) {
          Eigen::Vector3d vec1 = s.data().Bs[vis.to] * xportv0;
          vec1.normalize();
          Eigen::Vector3d vec2 = s.data().Bs[vis.to] * vecs[j];
          vec2.normalize();
          double curdot = (vec1).dot(vec2);
          if (curdot > bestdot) {
            bestidx = j;
            bestdot = curdot;
          }
        }

        // set vectors
        for (int j = 0; j < vecs.size(); j++) {
          comb_vecs.row(2 * npoly * vis.to + j) = vecs[(j + bestidx) % (vecs.size())].transpose();
        }
      }

      // queue neighbors
      for (int j = 0; j < 3; j++) {
        int nb = s.data().faceNeighbors(vis.to, j);
        if (nb != -1 && !visited[nb])
          q.push_back(Visit{ vis.to, nb });
      }

    }
  }

  return comb_vecs;
}

}