#include "RoundVectorFields.h"

#include <unordered_set>

#include "Weave.h"
#include "StripePatternIntegration.h"
#include "FieldSurface.h"
#include "Permutations.h"
#include "CoverMesh.h"

#include "../Surface.h"

void roundVectorFields(const Eigen::MatrixXd& mesh_pts, const Eigen::MatrixXi& mesh_faces, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& splitted_pts, Eigen::MatrixXi& splitted_faces, Eigen::VectorXd& theta, double global_rescaling, std::vector<Eigen::Vector3d>* splitted_vecs, std::vector<Eigen::Vector3d>* cut_pts, std::vector<std::vector<int>>* cut_edges, Eigen::VectorXd* err) {
	using namespace SurfaceFields;
    int nfields = vecs.rows() / mesh_faces.rows();
    int nfaces = mesh_faces.rows();
    // step 1: build weave
    std::unique_ptr<Weave> weave = std::make_unique<Weave>(mesh_pts, mesh_faces, nfields);
    for (int i = 0; i < nfaces; i++) {
        for (int j = 0; j < nfields; j++) {
            int vidx = weave->fs->vidx(i, j);
            weave->fs->vectorFields.segment<2>(vidx) = vecs.row(nfields * i + j).transpose();
        }
    }

    // step 2: comb the vector fields
    weave->combFieldsOnFieldSurface();
    reassignAllPermutations(*weave);


    // step 3: get the singularities
    std::vector<std::pair<int, int> > topsingularities;
    std::vector<std::pair<int, int> > geosingularities;
    findSingularVertices(*weave, topsingularities, geosingularities);

    struct pair_hash {
        inline std::size_t operator()(const std::pair<int, int>& v) const {
            return std::hash<int>()(v.first) ^ std::hash<int>()(v.second);
        }
    };

    std::unordered_set<std::pair<int, int>, pair_hash> todelete_set;
    todelete_set.insert(topsingularities.begin(), topsingularities.end());
    todelete_set.insert(geosingularities.begin(), geosingularities.end());

    std::vector<std::pair<int, int>> todelete;
    todelete.insert(todelete.end(), todelete_set.begin(), todelete_set.end());

    // step 4: cover mesh
    CoverMesh* cover_mesh = weave->createCover(todelete);

    // step 5: integration
    std::unique_ptr<SurfaceFields::StripePatternsGlobalIntegration> ptr = std::make_unique<SurfaceFields::StripePatternsGlobalIntegration>();
    cover_mesh->integrateField(ptr.get(), global_rescaling);

    // step 6: get the splitted mesh
    std::vector<Eigen::Vector3d> tmp_splitted_vecs;
    std::vector<Eigen::Vector3d> tmp_cut_pts;
    std::vector<std::vector<int>> tmp_cut_edges;
    cover_mesh->createVisualization(splitted_pts, splitted_faces, theta, tmp_splitted_vecs, tmp_cut_pts, tmp_cut_edges);

    


    if (splitted_vecs) {
        *splitted_vecs = std::move(tmp_splitted_vecs);
    }
    if (cut_edges) {
        *cut_edges = std::move(tmp_cut_edges);
    }
    if (cut_pts) {
        *cut_pts = std::move(tmp_cut_pts);
    }

    if (err) {
        Eigen::VectorXd cover_err;
        cover_mesh->gradThetaDeviation(cover_err, global_rescaling);
        *err = std::move(cover_err);
    }

    delete cover_mesh;
}