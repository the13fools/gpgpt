#include "CoverMesh.h"

#include <queue>
#include <set>
#include <Eigen/Dense>

#include <igl/writeOBJ.h>
#include <igl/is_vertex_manifold.h>
#include <igl/is_edge_manifold.h>
#include <igl/cotmatrix_entries.h>
#include <igl/facet_components.h>
#include <igl/remove_unreferenced.h>
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>

#include "FieldSurface.h"
#include "FieldIntegration.h"
#include "GNGlobalIntegration.h"
#include "SpectralLocalIntegration.h"
#include "CurlLocalIntegration.h"

#include "../Surface.h"


namespace SurfaceFields {
    CoverMesh::CoverMesh(const Surface& originalSurf, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::VectorXi& oldToNewVertMap, const Eigen::MatrixXd& field, int ncovers)
    {
        originalSurf_ = new Surface(originalSurf);
        fs = new FieldSurface(V, F, 1);
        int nfaces = F.rows();
        ncovers_ = ncovers;
        for (int i = 0; i < nfaces; i++)
        {
            fs->vectorFields.segment<2>(fs->vidx(i, 0)) = field.row(i).transpose();
        }

        theta.resize(fs->nVerts());
        theta.setZero();

        scales.resize(fs->nFaces());
        scales.setZero();

        renderScale_ = 1.0;

        initializeSplitMesh(oldToNewVertMap);
    }

    CoverMesh::~CoverMesh()
    {
        delete fs;
        if (data_.splitMesh)
            delete data_.splitMesh;
        delete originalSurf_;
    }

    double CoverMesh::barycentric(double val1, double val2, double target)
    {
        return (target - val1) / (val2 - val1);
    }

    bool CoverMesh::crosses(double isoval, double val1, double val2, double minval, double maxval, double& bary)
    {
        double halfperiod = 0.5 * (maxval - minval);
        if (fabs(val2 - val1) <= halfperiod)
        {
            bary = barycentric(val1, val2, isoval);
            if (bary >= 0 && bary < 1)
                return true;
            return false;
        }
        if (val1 < val2)
        {
            double wrapval1 = val1 + (maxval - minval);
            bary = barycentric(wrapval1, val2, isoval);
            if (bary >= 0 && bary < 1)
                return true;
            double wrapval2 = val2 + (minval - maxval);
            bary = barycentric(val1, wrapval2, isoval);
            if (bary >= 0 && bary < 1)
                return true;
        }
        else
        {
            double wrapval1 = val1 + (minval - maxval);
            bary = barycentric(wrapval1, val2, isoval);
            if (bary >= 0 && bary < 1)
                return true;
            double wrapval2 = val2 + (maxval - minval);
            bary = barycentric(val1, wrapval2, isoval);
            if (bary >= 0 && bary < 1)
                return true;
        }
        return false;
    }

    void CoverMesh::createVisualization(Eigen::MatrixXd& V, Eigen::MatrixXi& F,
        Eigen::MatrixXd& edgePts, Eigen::MatrixXi& edgeSegs, Eigen::MatrixXd& colors,
        Eigen::MatrixXd& cutPts1, Eigen::MatrixXd& cutPts2, Eigen::MatrixXd& cutColors,
        bool hideVectors, double vectorLength)
    {
        int splitFace = data_.splitMesh->nFaces();
        int origverts = originalSurf_->nVerts();
        int origfaces = originalSurf_->nFaces();
        V = data_.splitMesh->data().V;
        F = data_.splitMesh->data().F;

        if (hideVectors)
        {
            edgePts.resize(0, 3);
            edgeSegs.resize(0, 3);
            colors.resize(0, 3);
        }
        else
        {
            edgePts.resize(2 * splitFace, 3);
            edgePts.setZero();
            edgeSegs.resize(splitFace, 2);
            colors.resize(splitFace, 3);

            for (int c = 0; c < ncovers_; c++)
            {
                for (int i = 0; i < origfaces; i++)
                {
                    Eigen::Vector3d centroid;
                    centroid.setZero();
                    for (int j = 0; j < 3; j++)
                        centroid += renderScale_ * originalSurf_->data().V.row(originalSurf_->data().F(i, j));
                    centroid /= 3.0;
                    centroid += data_.splitOffsets[c];

                    edgePts.row(2 * c * origfaces + 2 * i) = centroid.transpose();
                    Eigen::Vector3d vec = originalSurf_->data().Bs[i] * fs->v(c * origfaces + i, 0);
                    vec *= vectorLength * renderScale_ * fs->data().averageEdgeLength / vec.norm() * sqrt(3.0) / 6.0 * 0.75;
                    edgePts.row(2 * c * origfaces + 2 * i + 1) = (centroid + vec).transpose();
                    edgeSegs(c * origfaces + i, 0) = 2 * (c * origfaces + i);
                    edgeSegs(c * origfaces + i, 1) = 2 * (c * origfaces + i) + 1;
                    colors.row(c * origfaces + i) = Eigen::Vector3d(0, 0, 0).transpose();
                }
            }
        }

        int ncutedges = data_.splitMeshCuts.size();
        int nsliceedges = slicedEdges.size();
        cutPts1.resize(ncutedges + nsliceedges, 3);
        cutPts2.resize(ncutedges + nsliceedges, 3);
        cutColors.resize(ncutedges + nsliceedges, 3);
        for (int i = 0; i < ncutedges; i++)
        {
            int edgeid = data_.splitMeshCuts[i];
            int v0 = data_.splitMesh->data().edgeVerts(edgeid, 0);
            int v1 = data_.splitMesh->data().edgeVerts(edgeid, 1);
            Eigen::Vector3d n(0, 0, 0);
            int f0 = data_.splitMesh->data().E(edgeid, 0);
            int f1 = data_.splitMesh->data().E(edgeid, 1);
            if (f0 != -1)
                n += data_.splitMesh->faceNormal(f0);
            if (f1 != -1)
                n += data_.splitMesh->faceNormal(f1);
            Eigen::Vector3d offset = 0.0001 * n / n.norm();
            cutPts1.row(i) = data_.splitMesh->data().V.row(v0) + offset.transpose();
            cutPts2.row(i) = data_.splitMesh->data().V.row(v1) + offset.transpose();
            cutColors.row(i) = Eigen::RowVector3d(0.9, .1, .9);
        }
        for (int i = 0; i < nsliceedges; i++)
        {
            int v0 = data_.splitMesh->data().edgeVerts(slicedEdges[i], 0);
            int v1 = data_.splitMesh->data().edgeVerts(slicedEdges[i], 1);
            Eigen::Vector3d n(0, 0, 0);
            int f0 = data_.splitMesh->data().E(slicedEdges[i], 0);
            int f1 = data_.splitMesh->data().E(slicedEdges[i], 1);
            if (f0 != -1)
                n += data_.splitMesh->faceNormal(f0);
            if (f1 != -1)
                n += data_.splitMesh->faceNormal(f1);
            Eigen::Vector3d offset = 0.0001 * n / n.norm();
            cutPts1.row(i + ncutedges) = data_.splitMesh->data().V.row(v0) + offset.transpose();
            cutPts2.row(i + ncutedges) = data_.splitMesh->data().V.row(v1) + offset.transpose();
            cutColors.row(i + ncutedges) = Eigen::RowVector3d(0.1, .9, .9);
        }
    }

    int CoverMesh::visMeshToCoverMesh(int vertid)
    {
        return data_.splitToCoverVerts[vertid];
    }


    void CoverMesh::integrateField(LocalFieldIntegration* lmethod, GlobalFieldIntegration* gmethod, double globalScale)
    {
        int globalverts = fs->nVerts();
        theta.resize(globalverts);
        theta.setZero();
        scales.resize(fs->nFaces());
        scales.setZero();

        // create a mesh without deleted faces
        int undeletedFaces = fs->numUndeletedFaces();
        Eigen::MatrixXi undelF(undeletedFaces, 3);
        Eigen::VectorXi undelFaceMap(undeletedFaces);
        int fid = 0;
        int globalfaces = fs->nFaces();
        for (int i = 0; i < globalfaces; i++)
        {
            if (!fs->isFaceDeleted(i))
            {
                undelFaceMap[fid] = i;
                undelF.row(fid) = fs->data().F.row(i);
                fid++;
            }
        }

        // separate cut mesh into connected components
        Eigen::VectorXi components;
        igl::facet_components(undelF, components);
        int ncomponents = 0;
        for (int i = 0; i < components.size(); i++)
            ncomponents = std::max(ncomponents, components[i]);
        ncomponents++;
        std::vector<int> componentsizes;
        for (int i = 0; i < ncomponents; i++)
            componentsizes.push_back(0);
        for (int i = 0; i < components.size(); i++)
            componentsizes[components[i]]++;
        std::cout << "Covering mesh has " << ncomponents << " connected components" << std::endl;
        // loop over the connected components
        for (int component = 0; component < ncomponents; component++)
        {
            std::cout << "Component " << component << ": " << componentsizes[component] << " faces" << std::endl;
            // faces for just this connected component
            Eigen::VectorXi compFacesToGlobal(componentsizes[component]);
            Eigen::MatrixXi compF(componentsizes[component], 3);
            Eigen::MatrixXd compField(componentsizes[component], 2);
            int idx = 0;
            for (int i = 0; i < components.size(); i++)
            {
                if (components[i] == component)
                {
                    compFacesToGlobal[idx] = undelFaceMap[i];
                    compF.row(idx) = undelF.row(i);
                    Eigen::Vector2d vec = fs->v(undelFaceMap[i], 0);
                    double vecnorm = (fs->data().Bs[undelFaceMap[i]] * vec).norm();
                    compField.row(idx) = vec.transpose() / vecnorm;
                    idx++;
                }
            }

            Eigen::MatrixXd prunedV;
            Eigen::MatrixXi prunedF;
            Eigen::VectorXi I;
            igl::remove_unreferenced(fs->data().V, compF, prunedV, prunedF, I);
            // connected component surface
            Surface surf(prunedV, prunedF);

            std::cout << "Built connected component surface" << std::endl;

            // component theta and s
            Eigen::VectorXd compS;
            Eigen::VectorXd compTheta;
            lmethod->locallyIntegrateOneComponent(surf, compField, compS);

            double maxS = 0;
            for (int i = 0; i < compS.size(); i++)
            {
                if (fabs(compS[i]) > maxS)
                {
                    maxS = fabs(compS[i]);
                }
            }

            double s_scale = 3.1415 / fs->data().averageEdgeLength / maxS;
            compS *= globalScale * s_scale;

            for (int i = 0; i < componentsizes[component]; i++)
            {
                scales[compFacesToGlobal[i]] = compS[i];
            }

            gmethod->globallyIntegrateOneComponent(surf, compField, compS, compTheta);


            // map component theta to the global theta vector
            for (int i = 0; i < globalverts; i++)
            {
                if (I[i] != -1)
                    theta[i] = compTheta[I[i]];
            }
        }
    }

    double CoverMesh::inversePowerIteration(Eigen::SparseMatrix<double>& M, Eigen::VectorXd& evec, int iters)
    {
        evec.resize(M.cols());
        evec.setRandom();
        evec /= evec.norm();
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver(M);
        for (int i = 0; i < iters; i++)
        {
            Eigen::VectorXd newvec = solver.solve(evec);
            evec = newvec / newvec.norm();
        }
        return evec.transpose() * M * evec;
    }

    void CoverMesh::initializeSplitMesh(const Eigen::VectorXi& oldToNewVertMap)
    {
        data_.splitToCoverVerts = oldToNewVertMap;
        int facespercover = fs->nFaces() / ncovers_;
        int rows = 2;
        int meshesperrow = ncovers_ / rows + (ncovers_ % rows == 0 ? 0 : 1);
        data_.splitOffsets.clear();
        for (int i = 0; i < ncovers_; i++)
        {
            int row = i / meshesperrow;
            int col = i % meshesperrow;
            double dy = (-1.1 * row + (1.1) * (rows - row - 1)) / double(rows);
            double dx = (1.1 * col + (-1.1) * (meshesperrow - col - 1)) / double(meshesperrow);
            data_.splitOffsets.push_back(Eigen::Vector3d(dx, dy, 0.0));
        }

        int origverts = originalSurf_->nVerts();
        int origfaces = originalSurf_->nFaces();
        int newverts = ncovers_ * origverts;
        int newfaces = ncovers_ * origfaces;
        Eigen::MatrixXd V(newverts, 3);
        Eigen::MatrixXi F(newfaces, 3);
        renderScale_ = 1.0 / std::max(rows, meshesperrow);
        for (int i = 0; i < ncovers_; i++)
        {
            for (int j = 0; j < origverts; j++)
            {
                V.row(i * origverts + j) = data_.splitOffsets[i].transpose() + renderScale_ * originalSurf_->data().V.row(j);
            }
            for (int j = 0; j < origfaces; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    F(i * origfaces + j, k) = i * origverts + originalSurf_->data().F(j, k);
                }
            }
        }
        data_.splitMesh = new Surface(V, F);

        data_.coverToSplitVerts.clear();
        for (int i = 0; i < oldToNewVertMap.rows(); i++)
        {
            data_.coverToSplitVerts[oldToNewVertMap[i]].push_back(i);
        }

        data_.splitMeshCuts.clear();
        for (int i = 0; i < newfaces; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                int edge = fs->data().faceEdges(i, j);
                if (edge == -1)
                    continue;
                int f0 = fs->data().E(edge, 0);
                int f1 = fs->data().E(edge, 1);
                if (f0 == -1 || f1 == -1)
                    continue;
                int face0copy = f0 / origfaces;
                int face1copy = f1 / origfaces;
                if (face0copy != face1copy)
                {
                    data_.splitMeshCuts.push_back(data_.splitMesh->data().faceEdges(i, j));
                }
            }
        }
    }

    const Surface& CoverMesh::splitMesh() const
    {
        return *data_.splitMesh;
    }

    static double periodicDiff(double a, double b)
    {
        double diff = b - a;
        if (diff < -M_PI)
            diff += 2.0 * M_PI;
        else if (diff > M_PI)
            diff -= 2.0 * M_PI;
        return diff;
    }

    void CoverMesh::gradThetaDeviation(Eigen::VectorXd& error) const
    {
        int nfaces = fs->nFaces();
        error.resize(nfaces);
        error.setZero();
        for (int i = 0; i < nfaces; i++)
        {
            if (fs->isFaceDeleted(i))
                error[i] = 0;
            else
            {
                double vals[3];
                for (int j = 0; j < 3; j++)
                {
                    vals[j] = theta[fs->data().F(i, j)];
                    Eigen::Vector2d diffs;
                    diffs[0] = periodicDiff(vals[0], vals[1]);
                    diffs[1] = periodicDiff(vals[0], vals[2]);
                    Eigen::Matrix<double, 3, 2> B = fs->data().Bs[i];
                    Eigen::Matrix2d BTB = B.transpose() * B;
                    Eigen::Matrix2d BTBinv = BTB.inverse();
                    Eigen::Vector2d grad = BTBinv * diffs;
                    Eigen::Vector3d grademb = B * grad;
                    Eigen::Vector3d vemb = B * fs->data().Js.block<2, 2>(2 * i, 0) * fs->v(i, 0);
                    Eigen::Vector3d n = fs->faceNormal(i);
                    double theta = asin(grademb.cross(vemb).dot(n) / grademb.norm() / vemb.norm());
                    error[i] = fabs(theta) / (0.5 * M_PI);
                }
            }
        }
    }

}

