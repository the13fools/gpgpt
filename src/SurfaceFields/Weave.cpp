#include "Weave.h"

#include <math.h>
#include <deque>
#include <queue>
#include <algorithm>
#include <set>
#include <map>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <Eigen/Eigenvalues> 
#include <Eigen/Geometry>

#include <igl/remove_unreferenced.h>

#include "../Surface.h"
#include "CoverMesh.h"
#include "VectorFields2PolyVectors.h"


Weave::Weave(Eigen::MatrixXd Vtmp, Eigen::MatrixXi Ftmp, int m)
{
    fs = new FieldSurface(Vtmp, Ftmp, m);
}

Weave::Weave(const Weave& w)
{
    fs = nullptr;
    operator=(w);
}

Weave& Weave::operator=(const Weave& w)
{
    delete fs;
    handles = w.handles;
    cuts = w.cuts;
    fs = new FieldSurface(w.fs->data().V, w.fs->data().F, w.fs->nFields());
    fs->vectorFields = w.fs->vectorFields;
    fs->Ps_ = w.fs->Ps_;

    return *this;
}

Weave::~Weave()
{
    delete fs;
}

bool Weave::addHandle(Handle h)
{
    if (h.face < 0 || h.face > fs->nFaces())
        return false;
    if (h.field < 0 || h.field > fs->nFields())
        return false;

    Eigen::Vector3d extrinsic = fs->data().Bs[h.face] * h.dir;
    double mag = extrinsic.norm();
    h.dir /= mag;
    handles.push_back(h);
    return true;
}

using namespace std;


static Eigen::Vector3d parallelTransport(const Eigen::Vector3d& v, const Eigen::Vector3d& e1, const Eigen::Vector3d& e2)
{
    Eigen::Vector3d t1 = e1 / e1.norm();
    Eigen::Vector3d t2 = e2 / e2.norm();
    Eigen::Vector3d n = t1.cross(t2);
    if (n.norm() < 1e-8)
        return v;
    n /= n.norm();
    Eigen::Vector3d p1 = n.cross(t1);
    Eigen::Vector3d p2 = n.cross(t2);
    return v.dot(n) * n + v.dot(t1) * t2 + v.dot(p1) * p2;
}

static double angle(Eigen::Vector3d v, Eigen::Vector3d w, Eigen::Vector3d n)
{
    return 2.0 * atan2(v.cross(w).dot(n), v.norm() * w.norm() + v.dot(w));
}

static void alignFrame(std::vector<Eigen::Vector3d>& from, const std::vector<Eigen::Vector3d>& to, const Eigen::Vector3d& n)
{
    int m = to.size();
    std::vector<int> cand;
    std::vector<int> bestperm;
    double bestangle = std::numeric_limits<double>::infinity();
    for (int i = 0; i < m; i++)
        cand.push_back(i);
    do
    {
        double totangle = 0;
        for (int j = 0; j < m; j++)
        {
            Eigen::Vector3d v1 = to[j];
            Eigen::Vector3d v2 = from[cand[j]];
            double theta = angle(v1, v2, n);
            totangle += theta * theta;
        }
        if (totangle < bestangle)
        {
            bestangle = totangle;
            bestperm = cand;
        }
    } while (std::next_permutation(cand.begin(), cand.end()));

    std::vector<Eigen::Vector3d> newvecs;
    for (int i = 0; i < m; i++)
        newvecs.push_back(from[bestperm[i]]);
    from = newvecs;
}

std::pair<std::vector<Eigen::Vector3d>, std::vector<std::vector<int>>> Weave::getPermuatedEdges() {
    int nCover = fs->nFields() * 2;
    int nfaces = fs->nFaces();
    int nverts = fs->nVerts();
    std::vector<Eigen::MatrixXd> perms;
    Eigen::MatrixXd perm;
    perms = _augmentPs();
    Eigen::MatrixXd eye = Eigen::MatrixXd::Identity(nCover, nCover);
    int not_identity = 0;

    std::vector<Eigen::Vector3d> pts;
    std::vector<std::vector<int>> edges;

    for (int e = 0; e < fs->nEdges(); e++)
    {
        perm = perms[e];
        int f1Id = fs->data().E(e, 0);
        int f2Id = fs->data().E(e, 1);
        if (f1Id == -1 || f2Id == -1)
            continue;
        int v1ID = fs->data().edgeVerts(e, 0);
        int v2ID = fs->data().edgeVerts(e, 1);
        int v1f1 = -1, v2f1 = -1, v1f2 = -1, v2f2 = -1;
        for (int i = 0; i < 3; i++)
        { // find the vid at face (0,1,or 2)
            if (fs->data().F(f1Id, i) == v1ID) v1f1 = i;
            if (fs->data().F(f1Id, i) == v2ID) v2f1 = i;
            if (fs->data().F(f2Id, i) == v1ID) v1f2 = i;
            if (fs->data().F(f2Id, i) == v2ID) v2f2 = i;
        }
        assert((v1f1 != -1) && (v2f1 != -1) && (v1f2 != -1) && (v2f2 != -1));
        if ((perm - eye).norm() != 0)
        {   // perm == I case
            pts.push_back(fs->data().V.row(v1ID));
            pts.push_back(fs->data().V.row(v2ID));
            edges.push_back({ (int)pts.size() - 2, (int)pts.size() - 1 });
        }
    }

    return { pts, edges };
}


CoverMesh* Weave::createCover(const std::vector<std::pair<int, int> >& singularities) const
{
    int nCover = fs->nFields() * 2;
    int nfaces = fs->nFaces();
    int nverts = fs->nVerts();
    std::vector<Eigen::MatrixXd> perms;
    Eigen::MatrixXd perm;
    perms = _augmentPs();
    Eigen::MatrixXd eye = Eigen::MatrixXd::Identity(nCover, nCover);
    // Compute points to glue
    vector<vector<long> > adj_list(nCover * nfaces * 3);

    int not_identity = 0;

    for (int e = 0; e < fs->nEdges(); e++)
    {
        perm = perms[e];
        int f1Id = fs->data().E(e, 0);
        int f2Id = fs->data().E(e, 1);
        if (f1Id == -1 || f2Id == -1)
            continue;
        int v1ID = fs->data().edgeVerts(e, 0);
        int v2ID = fs->data().edgeVerts(e, 1);
        int v1f1 = -1, v2f1 = -1, v1f2 = -1, v2f2 = -1;
        for (int i = 0; i < 3; i++)
        { // find the vid at face (0,1,or 2)
            if (fs->data().F(f1Id, i) == v1ID) v1f1 = i;
            if (fs->data().F(f1Id, i) == v2ID) v2f1 = i;
            if (fs->data().F(f2Id, i) == v1ID) v1f2 = i;
            if (fs->data().F(f2Id, i) == v2ID) v2f2 = i;
        }
        assert((v1f1 != -1) && (v2f1 != -1) && (v1f2 != -1) && (v2f2 != -1));
        if ((perm - eye).norm() == 0)
        { // perm == I case
            for (int l = 0; l < nCover; l++)
            {
                long v1f1_idx = v1f1 + f1Id * 3 + l * 3 * nfaces;
                long v2f1_idx = v2f1 + f1Id * 3 + l * 3 * nfaces;
                long v1f2_idx = v1f2 + f2Id * 3 + l * 3 * nfaces;
                long v2f2_idx = v2f2 + f2Id * 3 + l * 3 * nfaces;
                adj_list[v1f1_idx].push_back(v1f2_idx);
                adj_list[v1f2_idx].push_back(v1f1_idx);
                adj_list[v2f1_idx].push_back(v2f2_idx);
                adj_list[v2f2_idx].push_back(v2f1_idx);
            }
        }
        else
        { // perm != I case
            not_identity++;
            for (int l1 = 0; l1 < nCover; l1++)
            {
                int l2 = -1;
                for (int j = 0; j < nCover; j++)
                    if (perm(l1, j) == 1) { l2 = j; break; }
                long v1f1_idx = v1f1 + f1Id * 3 + l1 * 3 * nfaces;
                long v2f1_idx = v2f1 + f1Id * 3 + l1 * 3 * nfaces;
                long v1f2_idx = v1f2 + f2Id * 3 + l2 * 3 * nfaces;
                long v2f2_idx = v2f2 + f2Id * 3 + l2 * 3 * nfaces;
                adj_list[v1f1_idx].push_back(v1f2_idx);
                adj_list[v1f2_idx].push_back(v1f1_idx);
                adj_list[v2f1_idx].push_back(v2f2_idx);
                adj_list[v2f2_idx].push_back(v2f1_idx);
            }
        }
    }

    std::cout << "not identity number: " << not_identity << std::endl;

    // Do some glueing
    vector<vector<long> > gluePointList;
    vector<bool> toSearchFlag(nCover * nfaces * 3, 1);
    for (int i = 0; i < nCover * nfaces * 3; i++)
    {
        if (toSearchFlag[i] == 0)
            continue;
        vector<long> gluePoint = _BFS_adj_list(adj_list, i);
        gluePointList.push_back(gluePoint);
        for (int j = 0; j < gluePoint.size(); j++)
            toSearchFlag[gluePoint[j]] = 0;
    }
    int nNewPoints = gluePointList.size();
    Eigen::MatrixXd VAug = Eigen::MatrixXd::Zero(nNewPoints, 3); // |gluePointList| x 3
    Eigen::VectorXi oldId2NewId(nCover * nverts);
    oldId2NewId.setConstant(-1);
    vector<long> encodeDOldId2NewId(nCover * 3 * nfaces);
    for (int i = 0; i < nNewPoints; i++)
    { // Assign a new Vertex for each group of glue vetices
        long encodedVid = gluePointList[i][0];
        int layerId = floor(encodedVid / (nfaces * 3));
        int atFace = floor((encodedVid - layerId * nfaces * 3) / 3);
        int atVid = encodedVid - layerId * nfaces * 3 - 3 * atFace;
        int vid = fs->data().F(atFace, atVid);
        for (int j = 0; j < 3; j++)
            VAug(i, j) = fs->data().V(vid, j);
        for (int j = 0; j < gluePointList[i].size(); j++)
        { // Maintain a vid mapping
            encodedVid = gluePointList[i][j];
            layerId = floor(encodedVid / (nfaces * 3));
            atFace = floor((encodedVid - layerId * nfaces * 3) / 3);
            atVid = encodedVid - layerId * nfaces * 3 - 3 * atFace;
            oldId2NewId[vid + layerId * nverts] = i;
            encodeDOldId2NewId[gluePointList[i][j]] = i;
        }
    }

    Eigen::MatrixXi FAug = Eigen::MatrixXi::Zero(nCover * nfaces, 3);; // |gluePointList| x 3
    for (int cId = 0; cId < nCover; cId++)
    {
        for (int fId = 0; fId < nfaces; fId++)
        {
            int id0 = (fId + cId * nfaces) * 3;
            int id1 = (fId + cId * nfaces) * 3 + 1;
            int id2 = (fId + cId * nfaces) * 3 + 2;
            FAug(fId + cId * nfaces, 0) = encodeDOldId2NewId[id0];
            FAug(fId + cId * nfaces, 1) = encodeDOldId2NewId[id1];
            FAug(fId + cId * nfaces, 2) = encodeDOldId2NewId[id2];
        }
    }
    Eigen::MatrixXd flattenedField(nCover * nfaces, 2);
    for (int cId = 0; cId < nCover; cId++)
    {
        for (int fId = 0; fId < nfaces; fId++)
        {
            int field = cId % fs->nFields();
            double sign = (cId < fs->nFields() ? 1.0 : -1.0);
            Eigen::Vector2d vec = sign * fs->v(fId, field).transpose();
            Eigen::Vector3d embvec = fs->data().Bs[fId] * vec;
            //double norm = embvec.norm();
            //flattenedField.row(fId + cId * nfaces) = vec / norm;
            flattenedField.row(fId + cId * nfaces) = vec;
        }
    }


    CoverMesh* ret = new CoverMesh(*this->fs, VAug, FAug, oldId2NewId, flattenedField, nCover);
    // transfer deleted faces
    for (int i = 0; i < nfaces; i++)
    {
        for (int j = 0; j < nCover; j++)
        {
            ret->fs->setFaceDeleted(j * nfaces + i, fs->isFaceDeleted(i));
        }
    }
    // also remove faces around singularities
    std::map<int, std::set<int> > delvertcovers;
    for (auto it : singularities)
    {
        delvertcovers[it.first].insert(it.second);
    }
    for (int i = 0; i < nfaces; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int vert = fs->data().F(i, j);
            for (int k = 0; k < fs->nFields(); k++)
            {
                if (delvertcovers[vert].count(k))
                {
                    ret->fs->setFaceDeleted(k * nfaces + i, true);
                    ret->fs->setFaceDeleted((k + fs->nFields()) * nfaces + i, true);
                }
            }
        }
    }
    return ret;
}

void Weave::combFieldsOnFieldSurface() {
    int m = fs->nFields();
    FieldSurface* new_fs = new FieldSurface(fs->data().V, fs->data().F, m);

    // set vector fields    
    int nfaces = fs->nFaces();
    std::vector<bool> visited(nfaces, false);
    // bread-first search of faces
    for (int i = 0; i < nfaces; i++)
    {
        if (visited[i])
            continue;

        struct Visit
        {
            int from;
            int to;
        };

        std::deque<Visit> q;
        q.push_back(Visit{ -1,i });
        while (!q.empty())
        {
            Visit vis = q.front();
            q.pop_front();
            if (visited[vis.to])
                continue;
            visited[vis.to] = true;

            std::vector<Eigen::Vector2d> poly_vecs;
            for (int pid = 0; pid < m; pid++) {
                poly_vecs.push_back(fs->v(vis.to, pid));
            }
            std::vector<Eigen::Vector2d> vecs = SurfaceFields::getVectors(*fs, vis.to, poly_vecs);

            if (vis.from == -1)
            {
                for (int j = 0; j < m; j++) // only half of the vectors are pushed
                {
                    int vidx = fs->vidx(vis.to, j);
                    new_fs->vectorFields.segment<2>(vidx) = vecs[j];
                }
            }
            else
            {
                // find optimal cyclic permutation
                // vector 0 on neighbor face
                Eigen::Vector2d v0 = new_fs->v(vis.from, 0);
                // tranport to current face
                int edge = -1;
                int fromside = -1;
                for (int j = 0; j < 3; j++)
                {
                    int eid = new_fs->data().faceEdges(vis.to, j);
                    int opp = 0;
                    if (new_fs->data().E(eid, opp) == vis.to)
                        opp = 1;
                    if (new_fs->data().E(eid, opp) == vis.from)
                    {
                        edge = eid;
                        fromside = opp;
                    }
                }
                assert(edge != -1);
                assert(new_fs->data().E(edge, fromside) == vis.from);
                assert(new_fs->data().E(edge, 1 - fromside) == vis.to);
                Eigen::Vector2d xportv0 = new_fs->data().Ts.block<2, 2>(2 * edge, 2 * fromside) * v0;
                int bestidx = 0;
                double bestdot = -1.0;
                for (int j = 0; j < 2 * m; j++)
                {
                    Eigen::Vector3d vec1 = new_fs->data().Bs[vis.to] * xportv0;
                    vec1.normalize();
                    Eigen::Vector3d vec2 = new_fs->data().Bs[vis.to] * vecs[j];
                    vec2.normalize();
                    double curdot = (vec1).dot(vec2);
                    if (curdot > bestdot)
                    {
                        bestidx = j;
                        bestdot = curdot;
                    }
                }

                // set vectors
                for (int j = 0; j < m; j++)
                {
                    int vidx = new_fs->vidx(vis.to, j);
                    new_fs->vectorFields.segment<2>(vidx) = vecs[(j + bestidx) % (2 * m)];

                }
            }

            // queue neighbors
            for (int j = 0; j < 3; j++)
            {
                int nb = new_fs->data().faceNeighbors(vis.to, j);
                if (nb != -1 && !visited[nb])
                    q.push_back(Visit{ vis.to, nb });
            }
        }
    }

    std::swap(fs->vectorFields, new_fs->vectorFields);
    delete new_fs;  // Honestly speaking, unique_ptr should be used
}


std::vector<long> Weave::_BFS_adj_list(std::vector<std::vector<long> >& adj_list, int startPoint) const
{
    vector<long> traversed;
    queue<long> que;
    traversed.push_back(startPoint);
    que.push(startPoint);
    while (que.size() > 0)
    {
        long curPoint = que.front();
        que.pop();
        for (int j = 0; j < adj_list[curPoint].size(); j++)
        {
            long to_add = adj_list[curPoint][j];
            bool visited = false;
            for (int i = 0; i < traversed.size(); i++)
            {
                if (traversed[i] == to_add) {
                    visited = true;
                    break;
                }
            }
            if (visited)
                continue;
            traversed.push_back(to_add);
            que.push(to_add);
        }
    }
    return traversed;
}

std::vector<Eigen::MatrixXd> Weave::_augmentPs() const
{
    int nCover = fs->nFields() * 2;
    int nfaces = fs->nFaces();
    int nverts = fs->nVerts();
    int nfields = fs->nFields();
    std::vector<Eigen::MatrixXd> perms;
    Eigen::MatrixXd perm;
    for (int e = 0; e < fs->nEdges(); e++)
    {
        perm = Eigen::MatrixXd::Zero(nCover, nCover);
        for (int j = 0; j < nfields; j++)
        {
            for (int k = 0; k < nfields; k++)
            {
                if (fs->Ps(e)(j, k) == 1)
                {
                    perm(j, k) = 1;
                    perm(j + nfields, k + nfields) = 1;
                }
                if (fs->Ps(e)(j, k) == -1)
                {
                    perm(j, k + nfields) = 1;
                    perm(j + nfields, k) = 1;
                }
            }
        }
        perms.push_back(perm);
    }
    return perms;
}

