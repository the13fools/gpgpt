#include "FieldSurface.h"
#include "../Surface.h"
#include <set>
#include <map>
#include <igl/remove_unreferenced.h>
#include <Eigen/Dense>

FieldSurface::FieldSurface(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, int numFields) : Surface(V,F), nFields_(numFields)
{
    int nfaces = data().F.rows();
    // initialize vector fields
    vectorFields.resize(2 * nfaces * numFields);
    vectorFields.setZero();
    vectorFields.segment(0, 2 * nfaces * numFields).setRandom();

    // initialize permutation matrices
    int nedges = nEdges();
    Ps_.resize(nedges);
    for (int i = 0; i < nedges; i++)
    {
        Ps_[i].resize(numFields, numFields);
        Ps_[i].setIdentity();
    }
    
    faceDeleted_.resize(nfaces);
    for (int i = 0; i < nfaces; i++) {
        faceDeleted_[i] = false;
    }
        
}


int FieldSurface::vidx(int face, int field) const
{
    return (2 * nFields()*face + 2 * field);
}

Eigen::Vector2d FieldSurface::v(int face, int field) const
{
    return vectorFields.segment<2>(vidx(face,field));
}

FieldSurface *FieldSurface::removeDeletedFacesFromMesh(std::map<int, int> &faceMap, std::map<int, int> &vertMap) const
{
    faceMap.clear();
    vertMap.clear();
    
    std::set<int> facesToDelete;
    
    std::map<std::pair<int, int>, int> edgeMap;
    for (int e = 0; e < nEdges(); e++) 
    { 
        std::pair<int, int> p(data().edgeVerts(e,0), data().edgeVerts(e,1));
        edgeMap[p] = e; 
    }

    for(int i=0; i<nFaces(); i++)
    {
        if(isFaceDeleted(i))
            facesToDelete.insert(i);
    }

    std::vector<int> faceIds;
    for (std::set<int>::iterator it = facesToDelete.begin(); it != facesToDelete.end(); ++it)
        faceIds.push_back(*it);

    if (faceIds.empty())
    {
        // nothing to do
        FieldSurface *ret = new FieldSurface(data().V, data().F, nFields());
        ret->vectorFields = vectorFields;
        ret->Ps_ = Ps_;
        for(int i=0; i<nFaces(); i++)
            faceMap[i] = i;
        for(int i=0; i<nVerts(); i++)
            vertMap[i] = i;
        return ret;
    }
    
    int fieldIdx = 0;
    int newNFaces = nFaces() - faceIds.size();
    Eigen::VectorXd vectorFields_clean = Eigen::VectorXd::Zero( 5*nFields()*newNFaces );
    Eigen::MatrixXi F_temp = Eigen::MatrixXi::Zero(newNFaces, 3); 

    std::vector<bool> newdeleted(newNFaces);
    
    for (int i = 0; i < nFaces(); i++)   
    { 
        if (facesToDelete.count(i))
            continue;
        
        // vec field
        vectorFields_clean.segment(2*fieldIdx*nFields(), 2*nFields()) = vectorFields.segment(2*i*nFields(), 2*nFields());
        // beta
        vectorFields_clean.segment(2*fieldIdx*nFields() + 2*newNFaces*nFields(), 2*nFields()) 
            = vectorFields.segment(2*i*nFields() + 2*nFaces()*nFields(), 2*nFields() );
        // alpha
        vectorFields_clean.segment(fieldIdx*nFields() + 4*newNFaces*nFields(), nFields()) 
            = vectorFields.segment(i * nFields() + 4*nFaces()*nFields(), nFields() );
        // faces 
        F_temp.row(fieldIdx) = data().F.row(i);
        newdeleted[fieldIdx] = faceDeleted_[i];
        faceMap[i] = fieldIdx;
        fieldIdx++;
    }

    Eigen::MatrixXd V_new;
    Eigen::MatrixXi F_new;
    Eigen::VectorXi marked; 
    Eigen::VectorXi vertMapVec; 

    igl::remove_unreferenced(data().V, F_temp, V_new, F_new, marked, vertMapVec);

    FieldSurface *result = new FieldSurface(V_new, F_new, nFields());
    result->vectorFields = vectorFields_clean;
    result->faceDeleted_ = newdeleted;

    for(int i=0; i<vertMapVec.rows(); i++)
    {
        vertMap[vertMapVec(i)] = i;
    }
           
    std::vector<Eigen::MatrixXi> Ps_new;
    for( int i = 0; i < result->nEdges(); i++) 
    {
        int v0 = vertMapVec( result->data().edgeVerts(i, 0) );
        int v1 = vertMapVec( result->data().edgeVerts(i, 1) );
        if ( v0 > v1 ) 
        { 
            std::swap(v0, v1);
        }
        std::pair<int, int> p(v0, v1);
        int oldEdge = edgeMap[p];
        Ps_new.push_back(Ps_[oldEdge]);
    }
    result->Ps_ = Ps_new;

    return result;
}

const Eigen::MatrixXi FieldSurface::Ps(int edge) const
{
    return Ps_[edge];
}

void FieldSurface::deleteVertex(int vid)
{
    int nfaces = data().F.rows();
    for(int i=0; i<nfaces; i++)
    {
        for(int j=0; j<3; j++)
        {
            if(data().F(i,j) == vid)
                faceDeleted_[i] = true;
        }
    }
}


void FieldSurface::undeleteAllFaces()
{
    int nfaces = data().F.rows();
    for(int i=0; i<nfaces; i++)
    {
        faceDeleted_[i] = false;
    }
}

void FieldSurface::setFaceDeleted(int fid, bool newstatus)
{
    faceDeleted_[fid] = newstatus;
}

int FieldSurface::numUndeletedFaces() const
{
    int ret = 0;
    int nfaces = data().F.rows();
    for(int i=0; i<nfaces; i++)
        if(!faceDeleted_[i])
            ret++;
    return ret;
}
