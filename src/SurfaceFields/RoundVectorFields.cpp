#include "RoundVectorFields.h"

#include <unordered_set>

#include "Weave.h"
#include "StripePatternIntegration.h"
#include "MIGlobalIntegration.h"
#include "FieldSurface.h"
#include "Permutations.h"
#include "CoverMesh.h"
#include <memory>

#include "../Surface.h"

void roundVectorFields(const Eigen::MatrixXd& mesh_pts, const Eigen::MatrixXi& mesh_faces, const Eigen::MatrixXd& vecs, Eigen::MatrixXd& splitted_pts, Eigen::MatrixXi& splitted_faces, Eigen::VectorXd& theta, SurfaceFields::GlobalFieldIntegration* round_op, Eigen::MatrixXd* grad_theta, double global_rescaling, Eigen::MatrixXd* splitted_vecs, std::vector<Eigen::Vector3d>* cut_pts, std::vector<std::vector<int>>* cut_edges) {
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
    std::cout << "comb fields done" << std::endl;


    // step 3: get the singularities
    std::vector<std::pair<int, int> > topsingularities;
    std::vector<std::pair<int, int> > geosingularities;
    findSingularVertices(*weave, topsingularities, geosingularities);

    struct pair_hash {
        inline std::size_t operator()(const std::pair<int, int>& v) const {
            return std::hash<int>()(v.first) ^ std::hash<int>()(v.second);
        }
    };

    std::cout << "Found singularity done: " << topsingularities.size() << " are topological singularities, " << geosingularities.size() << " are geometrical ones" << std::endl;

    std::unordered_set<std::pair<int, int>, pair_hash> todelete_set;
    todelete_set.insert(topsingularities.begin(), topsingularities.end());
    todelete_set.insert(geosingularities.begin(), geosingularities.end());

    std::vector<std::pair<int, int>> todelete;
    todelete.insert(todelete.end(), todelete_set.begin(), todelete_set.end());

    // step 4: cover mesh
    CoverMesh* cover_mesh = weave->createCover(todelete);

    // step 5: integration
    //std::unique_ptr<SurfaceFields::StripePatternsGlobalIntegration> ptr = std::make_unique<SurfaceFields::StripePatternsGlobalIntegration>();
    cover_mesh->integrateField(round_op, global_rescaling);

    std::cout << "Integration done" << std::endl;

    // step 6: get the splitted mesh
    Eigen::MatrixXd tmp_splitted_vecs;
    std::vector<Eigen::Vector3d> tmp_cut_pts;
    std::vector<std::vector<int>> tmp_cut_edges;
    cover_mesh->createVisualization(splitted_pts, splitted_faces, theta, tmp_splitted_vecs, tmp_cut_pts, tmp_cut_edges);

    std::cout << "create visualization done" << std::endl;

    if (cut_edges) {
        *cut_edges = std::move(tmp_cut_edges);
    }
    if (cut_pts) {
        *cut_pts = std::move(tmp_cut_pts);
    }
    if(grad_theta) {
        *grad_theta = cover_mesh->getGradTheta(global_rescaling);

        for(int i = 0; i < grad_theta->rows(); i++) {
            if(cover_mesh->fs->isFaceDeleted(i)) {
              grad_theta->row(i) = tmp_splitted_vecs.row(i);
            }
        }
        std::cout << "Compute grad theta done" << std::endl;
    }
    if (splitted_vecs) {
        *splitted_vecs = std::move(tmp_splitted_vecs);
    }

    delete cover_mesh;
}