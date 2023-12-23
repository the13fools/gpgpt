#ifndef COVERMESH_H
#define COVERMESH_H

#include <vector>
#include <Eigen/Core>
#include <Eigen/Sparse>

#include "GlobalFieldIntegration.h"

class FieldSurface;
class Weave;
class Surface;

struct CoverData
{
    Surface *splitMesh; // the covering mesh split into 2*m copies of the original mesh
    
    std::vector<Eigen::Vector3d> splitOffsets; // translation of each split mesh
    std::map<int, std::vector<int> > coverToSplitVerts; // map from vertex indices on the covering mesh to their "child" vertices on the split mesh
    Eigen::VectorXi splitToCoverVerts; // inverse of coverToSplitVerts

    std::vector<int> splitMeshCuts; // edges of the split mesh that are cuts
};

class CoverMesh
{
public:
    CoverMesh(const Surface &originalSurf, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXi &oldToNewVertMap, const Eigen::MatrixXd &field, int ncovers);
    ~CoverMesh();

    CoverMesh(const CoverMesh &) = delete;
    CoverMesh &operator=(const CoverMesh &) = delete;

    FieldSurface *fs;
    Eigen::VectorXd theta;
    
    void createVisualization(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXd& splitted_theta,
        std::vector<Eigen::Vector3d>& face_vectors,
        std::vector<Eigen::Vector3d>& cutPts, std::vector<std::vector<int>>& cut_edges);

    void integrateField(SurfaceFields::GlobalFieldIntegration* gmethod, double globalScale);
    void roundAntipodalCovers(int numISOLines);
    const Surface &splitMesh() const;
    void gradThetaDeviation(Eigen::VectorXd &error, double globalScale) const;
    
    // maps indices of vertices on the visualization mesh to corresponding "parent" vertices on the cover mesh
    int visMeshToCoverMesh(int vertid);

   
private:
    void initializeSplitMesh(const Eigen::VectorXi &oldToNewVertMap);
   
    std::vector<std::vector<Eigen::Vector3d> > isoNormal;

    CoverData data_;
    int ncovers_;
    Surface *originalSurf_;
    // edges along which the multiple cover is cut to create a topological disk
    std::vector<int> slicedEdges;
};

#endif
