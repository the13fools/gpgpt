#ifndef WEAVE_H
#define WEAVE_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include "FieldSurface.h"


// handles on the mesh
struct Handle
{
    int face; // which face the handle is on
    int field; // which of the m fields we're prescribing
    Eigen::Vector2d dir; // the vector itself in triangle barycentric coordinates
};

struct Cut
{ 
    std::vector<std::pair<int, int> > path; // (edge, orientation) list        
};

class CoverMesh;

class Weave
{
public:
    Weave(const std::string &objname, int m);
    Weave(const Eigen::MatrixXd Vtmp, const Eigen::MatrixXi Ftmp, int m);
    Weave(const Weave &w);
    ~Weave();

    Weave &operator=(const Weave &w);
   
    FieldSurface *fs;
    
    std::vector<Handle> handles; // handles on the vector fields
 
    std::vector<Cut> cuts; // list of cuts 

    bool addHandle(Handle h);  // this method will add a handle, also taking care to normalize the handle vector length

    int nHandles() const { return handles.size(); }    
    
    CoverMesh *createCover(const std::vector<std::pair<int, int> > &singularities) const;
    std::pair<std::vector<Eigen::Vector3d>, std::vector<std::vector<int>>> getPermuatedEdges();

    void combFieldsOnFieldSurface();

private:

    std::vector<long> _BFS_adj_list(std::vector<std::vector<long> > & relaxadj_list, int i) const;
    std::vector<Eigen::MatrixXd> _augmentPs() const;
};

#endif
