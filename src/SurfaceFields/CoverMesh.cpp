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
#include "CoMISoWrapper.h"
#include "Weave.h"



using namespace std;

// Serialize Eigen matrix to a binary file
static bool serializeMatrix(const Eigen::MatrixXd& mat, const std::string& filepath) {
    try {
        std::ofstream outFile(filepath, std::ios::binary);
        if (!outFile.is_open()) {
            std::cerr << "Error: Unable to open file for writing: " << filepath << std::endl;
            return false;
        }

        // Write matrix rows and cols
        int rows = static_cast<int>(mat.rows());
        int cols = static_cast<int>(mat.cols());
        outFile.write(reinterpret_cast<char*>(&rows), sizeof(int));
        outFile.write(reinterpret_cast<char*>(&cols), sizeof(int));

        // Write matrix data
        outFile.write(reinterpret_cast<const char*>(mat.data()), rows * cols * sizeof(double));
        outFile.close();
        return true;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: Unable to serialize matrix: " << e.what() << std::endl;
        return false;
    }
}

CoverMesh::CoverMesh(const Surface &originalSurf, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXi &oldToNewVertMap, const Eigen::MatrixXd &field, int ncovers)
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
   

    initializeSplitMesh(oldToNewVertMap);
}

CoverMesh::~CoverMesh()
{
    delete fs;
    if (data_.splitMesh)
        delete data_.splitMesh;
    delete originalSurf_;
}


void CoverMesh::createVisualization(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXd& splitted_theta,
    std::vector<Eigen::Vector3d>& face_vectors,
    std::vector<Eigen::Vector3d>& cut_pts, std::vector<std::vector<int>>& cut_edges)
{
    int splitFace = data_.splitMesh->nFaces();
    int origverts = originalSurf_->nVerts();
    int origfaces = originalSurf_->nFaces();
    V = data_.splitMesh->data().V;
    F = data_.splitMesh->data().F;
    splitted_theta.setZero(V.rows());

    face_vectors.resize(ncovers_ * origfaces);

    for (int c = 0; c < ncovers_; c++)
    {
        for (int i = 0; i < origfaces; i++)
        {
            face_vectors[c * origfaces + i] = originalSurf_->data().Bs[i] * fs->v(c * origfaces + i, 0);
        }
    }

    for (int i = 0; i < V.rows(); i++) {
        int cover_vid = visMeshToCoverMesh(i);
        splitted_theta[i] = theta[cover_vid];
    }

    int ncutedges = data_.splitMeshCuts.size();
    int nsliceedges = slicedEdges.size();
    for (int i = 0; i < ncutedges; i++)
    {
        int edgeid = data_.splitMeshCuts[i];
        int v0 = data_.splitMesh->data().edgeVerts(edgeid, 0);
        int v1 = data_.splitMesh->data().edgeVerts(edgeid, 1);
       
        cut_pts.push_back(data_.splitMesh->data().V.row(v0).transpose());
        cut_pts.push_back(data_.splitMesh->data().V.row(v1).transpose());

        cut_edges.push_back({ int(cut_pts.size()) - 2, int(cut_pts.size()) - 1 });

    }
    for (int i = 0; i < nsliceedges; i++)
    {
        int v0 = data_.splitMesh->data().edgeVerts(slicedEdges[i], 0);
        int v1 = data_.splitMesh->data().edgeVerts(slicedEdges[i], 1);
       
        cut_pts.push_back(data_.splitMesh->data().V.row(v0).transpose());
        cut_pts.push_back(data_.splitMesh->data().V.row(v1).transpose());

        cut_edges.push_back({ int(cut_pts.size()) - 2, int(cut_pts.size()) - 1 });
    }
}

int CoverMesh::visMeshToCoverMesh(int vertid)
{
    return data_.splitToCoverVerts[vertid];
}

void CoverMesh::roundAntipodalCovers(int numISOLines)
{
    // create a mesh without deleted faces
    int undeletedFaces = fs->numUndeletedFaces();
    Eigen::MatrixXi undelF(undeletedFaces, 3);
    Eigen::VectorXi undelFaceMap(undeletedFaces);
    int globalfaces = fs->nFaces();
    Eigen::VectorXi coverFacesToUndelFaces(globalfaces);
    int fid=0;
    for(int i=0; i<globalfaces; i++)
    {
        if(!fs->isFaceDeleted(i))
        {
            undelFaceMap[fid] = i;
            undelF.row(fid) = fs->data().F.row(i);
            coverFacesToUndelFaces[i] = fid;
            fid++;
        }        
        else
        {
            coverFacesToUndelFaces[i] = -1;
        }
    }

    Eigen::MatrixXd prunedV;
    Eigen::MatrixXi prunedF;
    Eigen::VectorXi I;
    igl::remove_unreferenced(fs->data().V, undelF, prunedV, prunedF, I);

    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(prunedV, prunedF, L);
    
    int nfields = ncovers_ / 2;
    int newverts = prunedV.rows();
    
    // map vertices on the cover mesh to those on the antipodal cover
    std::map<int, int> correspondences;
    int prunedfaces = prunedF.rows();
    int origfaces = originalSurf_->nFaces();
    for(int i=0; i<prunedfaces; i++)
    {    
        // face on the covering mesh    
        int origface = undelFaceMap[i];
        int cover = origface / origfaces;
        // only consider faces on positively-oriented covers
        if(cover < nfields)
        {            
            int corrface = origface + nfields*origfaces;
            if(coverFacesToUndelFaces[corrface] == -1)
                continue; // shouldn't happen probably
            for(int j=0; j<3; j++)
            {
                correspondences[fs->data().F(origface,j)] = fs->data().F(corrface, j);
            }
        }
    }
    int ncorrs = correspondences.size();
    std::vector<std::pair<int, int> > corrs;
    for(auto it : correspondences)
        corrs.push_back(std::pair<int,int>(it.first, it.second));
    
    double phase = 2.0 * M_PI / double(numISOLines);
    double offset = M_PI / numISOLines;

    Eigen::VectorXd result(newverts + ncorrs*nfields);
    result.setZero();

    std::vector<Eigen::Triplet<double> > Ccoeffs;
    int row = 0;
    for (int i = 0; i < ncorrs; i++)
    {
        for (int j = 0; j < nfields; j++)
        {
            int v1 = corrs[i].first;
            int v2 = corrs[i].second;
            assert(I[v1] != -1 && I[v2] != -1);
            double thetadiff = theta[v1] + theta[v2];
            Ccoeffs.push_back(Eigen::Triplet<double>(row, I[v1], 1.0));
            Ccoeffs.push_back(Eigen::Triplet<double>(row, I[v2], 1.0));
            Ccoeffs.push_back(Eigen::Triplet<double>(row, newverts + j * ncorrs + i, phase));
            Ccoeffs.push_back(Eigen::Triplet<double>(row, newverts + nfields * ncorrs, offset + thetadiff));
            row++;
        }
    }
    Eigen::SparseMatrix<double> C(row, newverts + nfields * ncorrs + 1);
    C.setFromTriplets(Ccoeffs.begin(), Ccoeffs.end());

    std::vector<Eigen::Triplet<double> > Acoeffs;
    for (int k = 0; k < L.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(L, k); it; ++it)
        {
            Acoeffs.push_back(Eigen::Triplet<double>(it.row(), it.col(), -it.value()));
        }
    }
    for (int i = 0; i < newverts + nfields*ncorrs; i++)
        Acoeffs.push_back(Eigen::Triplet<double>(i, i, 1e-4));
    Eigen::SparseMatrix<double> A(newverts + nfields*ncorrs, newverts + nfields*ncorrs);
    A.setFromTriplets(Acoeffs.begin(), Acoeffs.end());

    Eigen::VectorXd rhs(newverts + nfields*ncorrs);
    rhs.setZero();

    Eigen::VectorXi toRound(nfields*ncorrs);
    for (int i = 0; i < nfields*ncorrs; i++)
        toRound[i] = newverts + i;

    ComisoWrapper(C, A, result, rhs, toRound, 1e-6);
    std::cout << "Residual: " << (A*result - rhs).norm() << std::endl;
    Eigen::VectorXd ctest(newverts + nfields * ncorrs + 1);
    ctest.segment(0, newverts + nfields * ncorrs) = result;
    ctest[newverts + nfields * ncorrs] = 1.0;
    std::cout << "Constraint residual: " << (C*ctest).norm() << std::endl;
    for (int i = 0; i < fs->nVerts(); i++)
    {
        if (I[i] != -1)
        {
            theta[i] += result[I[i]];
            theta[i] = std::remainder(theta[i], 2.0*M_PI);
        }
    }
}

void CoverMesh::integrateField(SurfaceFields::GlobalFieldIntegration* gmethod, double globalScale)
{
    int globalverts = fs->nVerts();
    std::cout << "global verts: " << globalverts << std::endl;
    theta.resize(globalverts);
    theta.setZero();

    //comp_pts.clear();
    //comp_faces.clear();
    //comp_thetas.clear();

    // create a mesh without deleted faces
    int undeletedFaces = fs->numUndeletedFaces();
    Eigen::MatrixXi undelF(undeletedFaces, 3);
    Eigen::VectorXi undelFaceMap(undeletedFaces);
    int fid=0;
    int globalfaces = fs->nFaces();
    for(int i=0; i<globalfaces; i++)
    {
        if(!fs->isFaceDeleted(i))
        {
            undelFaceMap[fid] = i;
            undelF.row(fid) = fs->data().F.row(i);
            fid++;
        }        
    }

    //Eigen::MatrixXd global_vecs(globalfaces, 2);
    //for (int i = 0; i < globalfaces; i++) {
    //    global_vecs.row(i) = fs->v(i, 0).transpose();
    //}
    //igl::writeOBJ("global_mesh.obj", fs->data().V, fs->data().F);
    //serializeMatrix(global_vecs, "global_vec.brfa");
    
    // separate cut mesh into connected components
    Eigen::VectorXi components;    
    igl::facet_components(undelF, components);
    int ncomponents = 0;
    for(int i=0; i<components.size(); i++)    
        ncomponents = std::max(ncomponents, components[i]);
    ncomponents++;
    std::vector<int> componentsizes;
    for(int i=0; i<ncomponents; i++)
        componentsizes.push_back(0);
    for(int i=0; i<components.size(); i++)
        componentsizes[components[i]]++;    
    std::cout << "Covering mesh has " << ncomponents << " connected components" << std::endl;
    // loop over the connected components
    for(int component = 0; component < ncomponents; component++)
    {
        std::cout << "Component " << component << ": " << componentsizes[component] << " faces" << std::endl;
        // faces for just this connected component
        Eigen::VectorXi compFacesToGlobal(componentsizes[component]);
        Eigen::MatrixXi compF(componentsizes[component], 3);
        Eigen::MatrixXd compField(componentsizes[component], 2);
        int idx=0;
        for(int i=0; i<components.size(); i++)
        {
            if(components[i] == component)
            {
                compFacesToGlobal[idx] = undelFaceMap[i];
                compF.row(idx) = undelF.row(i);
                Eigen::Vector2d vec = fs->v(undelFaceMap[i], 0); 
                compField.row(idx) = vec.transpose();
                idx++;
            }
        }

        std::cout << "comp faces: " << compF.row(0) << std::endl;
        
        Eigen::MatrixXd prunedV;
        Eigen::MatrixXi prunedF;
        Eigen::VectorXi I;
        igl::remove_unreferenced(fs->data().V, compF, prunedV, prunedF, I);
        // connected component surface
        Surface surf(prunedV, prunedF);
        
        //igl::writeOBJ("pruned_mesh_" + std::to_string(component) + ".obj", prunedV, prunedF);
        //serializeMatrix(compField, "pruned_fields_" + std::to_string(component) + ".bfra");

        std::cout << "Built connected component surface" << std::endl;
        Eigen::VectorXd compTheta;

        compField *= globalScale;
        gmethod->globallyIntegrateOneComponent(surf, compField, compTheta);
        std::cout << "comp field's first row: " << compField.row(0) << std::endl;
        std::cout << "compTheta's norm: " << compTheta.norm() << std::endl;

        //comp_pts.push_back(prunedV);
        //comp_faces.push_back(prunedF);
        //comp_thetas.push_back(compTheta);
        
        int ncount = 0;
        // map component theta to the global theta vector
        for (int i = 0; i < globalverts; i++)
        {
            if (I[i] != -1) {
                theta[i] = compTheta[I[i]];
                ncount++;
                if ((fs->data().V.row(i) - prunedV.row(I[i])).norm()) {
                    std::cout << fs->data().V.row(i) << ", " << prunedV.row(I[i]) << std::endl;
                }
            }
                           
        }    
        std::cout << compTheta.size() << ", " << ncount << std::endl;
    }
}

void CoverMesh::initializeSplitMesh(const Eigen::VectorXi &oldToNewVertMap)
{
    data_.splitToCoverVerts = oldToNewVertMap;
    int facespercover = fs->nFaces() / ncovers_;
    int rows = 2;
    int meshesperrow = ncovers_ / rows + (ncovers_ % rows == 0 ? 0 : 1);
    data_.splitOffsets.clear();

    double ratio = (originalSurf_->data().V.colwise().maxCoeff() - originalSurf_->data().V.colwise().minCoeff()).norm();

    for (int i = 0; i < ncovers_; i++)
    {
        int row = i / meshesperrow;
        int col = i%meshesperrow;
        double dy = (-ratio * row + ratio * (rows - row - 1)) / double(rows);
        double dx = (ratio * col + (-ratio) * (meshesperrow - col - 1)) / double(meshesperrow);
        data_.splitOffsets.push_back(Eigen::Vector3d(dx, dy, 0.0));
    }

    int origverts = originalSurf_->nVerts();
    int origfaces = originalSurf_->nFaces();
    int newverts = ncovers_*origverts;
    int newfaces = ncovers_*origfaces;
    Eigen::MatrixXd V(newverts, 3);
    Eigen::MatrixXi F(newfaces, 3);

    for (int i = 0; i < ncovers_; i++)
    {
        for (int j = 0; j < origverts; j++)
        {
            V.row(i*origverts + j) = data_.splitOffsets[i].transpose() +  originalSurf_->data().V.row(j);
        }
        for (int j = 0; j < origfaces; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                F(i*origfaces + j, k) = i*origverts + originalSurf_->data().F(j, k);
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
            if(edge == -1)
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

const Surface &CoverMesh::splitMesh() const
{
    return *data_.splitMesh;
}


static double periodicDiff(double a, double b)
{
    double diff = b-a;
    if(diff < -M_PI)
        diff += 2.0*M_PI;
    else if(diff > M_PI)
        diff -= 2.0*M_PI;
    return diff;
}

void CoverMesh::gradThetaDeviation(Eigen::VectorXd &error, double globalScale) const
{
    int nfaces = fs->nFaces();
    error.resize(nfaces);
    error.setZero();
    for(int i=0; i<nfaces; i++)
    {
        if(fs->isFaceDeleted(i))
            error[i] = 0;
        else
        {
            double vals[3];
            for(int j=0; j<3; j++)
            {
                vals[j] = theta[fs->data().F(i,j)];
                Eigen::Vector2d diffs;
                diffs[0] = periodicDiff(vals[0], vals[1]);
                diffs[1] = periodicDiff(vals[0], vals[2]);
                Eigen::Matrix<double, 3, 2> B = fs->data().Bs[i];
                Eigen::Matrix2d BTB = B.transpose()*B;
                Eigen::Matrix2d BTBinv = BTB.inverse();
                Eigen::Vector2d grad = BTBinv * diffs;
                Eigen::Vector3d grademb = B * grad;
                Eigen::Vector3d vemb = B*fs->data().Js.block<2,2>(2*i,0)*fs->v(i,0) * globalScale;
                error[i] = (grademb - vemb).norm();
                //Eigen::Vector3d n = fs->faceNormal(i);
                //double theta = asin(grademb.cross(vemb).dot(n) / grademb.norm() / vemb.norm());
                //error[i] = fabs( theta ) / (0.5 * M_PI);
            }
        }
    }
}
