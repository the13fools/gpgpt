#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"

#include <igl/doublearea.h>
#include <igl/cotmatrix_entries.h>
#include <igl/cotmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/file_dialog_open.h>
#include <igl/file_dialog_save.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/principal_curvature.h>
#include <igl/hsv_to_rgb.h>

#include <unordered_set>
#include <sstream>

#include "../src/SurfaceFields/VectorFields2PolyVectors.h"
#include "../src/SurfaceFields/StripePatternIntegration.h"
#include "../src/SurfaceFields/MIGlobalIntegration.h"
#include "../src/SurfaceFields/GlobalFieldIntegration.h"
#include "../src/SurfaceFields/UpsampleSurfaceFields.h"
#include "../src/SurfaceFields/RoundVectorFields.h"

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


static bool serializeMatrixNew(const Eigen::MatrixXd& mat, const std::string& filepath, int vector_per_element) {
    try {
        std::ofstream outFile(filepath, std::ios::binary);
        if (!outFile.is_open()) {
            std::cerr << "Error: Unable to open file for writing: " << filepath << std::endl;
            return false;
        }

        // Write matrix rows and cols
        int rows = static_cast<int>(mat.rows());
        int cols = static_cast<int>(mat.cols());
        int vpe = static_cast<int>(vector_per_element);


        outFile.write("FRA 2", sizeof("FRA 2"));
        outFile.write(reinterpret_cast<char*>(&rows), sizeof(int));
        outFile.write(reinterpret_cast<char*>(&cols), sizeof(int));
        outFile.write(reinterpret_cast<char*>(&vpe), sizeof(int));

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

// Deserialize Eigen matrix from a binary file
static bool deserializeMatrix(Eigen::MatrixXd& mat, const std::string& filepath) {
    try {
        std::ifstream inFile(filepath, std::ios::binary);
        if (!inFile.is_open()) {
            std::cerr << "Error: Unable to open file for reading: " << filepath << std::endl;
            return false;
        }

        // Read matrix rows and cols
        int rows, cols;
        inFile.read(reinterpret_cast<char*>(&rows), sizeof(int));
        inFile.read(reinterpret_cast<char*>(&cols), sizeof(int));

        // Read matrix data
        mat.resize(rows, cols);
        inFile.read(reinterpret_cast<char*>(mat.data()), rows * cols * sizeof(double));
        inFile.close();
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error: Unable to deserialize matrix: " << e.what() << std::endl;
        return false;
    }
}

// Deserialize Eigen matrix from a binary file
static bool deserializeMatrixNew(Eigen::MatrixXd& mat, const std::string& filepath) {
    try {
        std::ifstream inFile(filepath, std::ios::binary);
        if (!inFile.is_open()) {
            std::cerr << "Error: Unable to open file for reading: " << filepath << std::endl;
            return false;
        }

        char* filename = new char[5];
        inFile.read(reinterpret_cast<char*>(&filename), sizeof("FRA 2"));
        // Read matrix rows and cols
        int rows, cols, vpe;
        inFile.read(reinterpret_cast<char*>(&rows), sizeof(int));
        inFile.read(reinterpret_cast<char*>(&cols), sizeof(int));
        inFile.read(reinterpret_cast<char*>(&vpe), sizeof(int));

        // Read matrix data
        mat.resize(rows, cols);
        inFile.read(reinterpret_cast<char*>(mat.data()), rows * cols * sizeof(double));
        inFile.close();
        return true;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: Unable to deserialize matrix: " << e.what() << std::endl;
        return false;
    }
}

enum IntegrationType {
  kBomme = 0,
  kKnoppel = 1,
};

enum IntegrabilityErrType {
  kStandard = 0,
  kL2 = 1,
  kL4 = 2,
};


static std::vector<Eigen::VectorXd> checkLocalIntegrability(Surface& surf, const Eigen::MatrixXd& vec, IntegrabilityErrType err_type) {
    int nedges = surf.nEdges();
    std::vector<Eigen::Vector3d> pts;
    int nfaces = surf.nFaces();
    int nvecs = vec.rows() / nfaces;

    std::vector<Eigen::VectorXd> edge_errs;

    if(err_type == kStandard) {
        edge_errs.resize(nvecs, Eigen::VectorXd::Zero(nedges));
    } else {
        edge_errs.resize(1, Eigen::VectorXd::Zero(nedges));
    }

    for (int i = 0; i < nedges; i++) {
        int f0 = surf.data().E(i, 0);
        int f1 = surf.data().E(i, 1);

        pts.push_back((surf.data().V.row(surf.data().edgeVerts(i, 1)) + surf.data().V.row(surf.data().edgeVerts(i, 0))) / 2.0);

        if (f0 != -1 && f1 != -1) {
            Eigen::Vector3d edge = (surf.data().V.row(surf.data().edgeVerts(i, 1)) - surf.data().V.row(surf.data().edgeVerts(i, 0))).transpose();
            std::vector<Eigen::Vector3d> v0_list(nvecs), v1_list(nvecs);
            for(int vi= 0; vi < nvecs; vi++) {
              v0_list[vi] = surf.data().Bs[f0] * (vec.row(nvecs * f0 + vi).transpose());
              v1_list[vi] = surf.data().Bs[f1] * (vec.row(nvecs * f1 + vi).transpose());
            }

            switch(err_type) {
            case kStandard: {
              for(int vi = 0; vi < nvecs; vi++) {
                auto& v0 = v0_list[vi];
                auto& v1 = v1_list[vi];
                edge_errs[vi][i] = std::abs(v1.dot(edge) - v0.dot(edge));

                if(i == 4630) {
                  std::cout << "v0: " << v0.transpose() << ", v1" << v1.transpose() << std::endl;
                  std::cout << "edge err:  " << edge_errs[vi][i] << std::endl;
                }
              }
              break;
            }
            case kL2:{
              for(int vi = 0; vi < nvecs; vi++) {
                auto& v0 = v0_list[vi];
                auto& v1 = v1_list[vi];
                edge_errs[0][i] += std::pow(v1.dot(edge), 2.0) - std::pow(v0.dot(edge), 2.0);
              }
              edge_errs[0][i] = std::abs(edge_errs[0][i]);
              break;
            }
            case kL4:{
              for(int vi = 0; vi < nvecs; vi++) {
                auto& v0 = v0_list[vi];
                auto& v1 = v1_list[vi];
                edge_errs[0][i] += std::pow(v1.dot(edge), 4.0) - std::pow(v0.dot(edge), 4.0);
              }
              edge_errs[0][i] = std::abs(edge_errs[0][i]);
              break;
            }
            }
        }
    }
    auto pt_surf = polyscope::registerPointCloud("err pts", pts);

    switch(err_type) {
    case kStandard: {
        for(int vi = 0; vi < nvecs; vi++) {
            auto err_plot = pt_surf->addScalarQuantity("Standard curl err vec " + std::to_string(vi), edge_errs[vi]);
            err_plot->setMapRange({edge_errs[vi].minCoeff(), edge_errs[vi].maxCoeff()});
        }
        break;
    }
    case kL2:{
        auto err_plot = pt_surf->addScalarQuantity("L2 curl err", edge_errs[0]);
        err_plot->setMapRange({edge_errs[0].minCoeff(), edge_errs[0].maxCoeff()});
        break;
    }
    case kL4:{
        auto err_plot = pt_surf->addScalarQuantity("L4 curl err", edge_errs[0]);
        err_plot->setMapRange({edge_errs[0].minCoeff(), edge_errs[0].maxCoeff()});
        break;
    }
    }

    return edge_errs;
}


static Eigen::VectorXd checkLocalIntegrabilityEuclidean(Surface& surf, const std::vector<Eigen::Vector3d>& face_vec) {
    int nedges = surf.nEdges();
    Eigen::VectorXd edge_err = Eigen::VectorXd::Zero(nedges);
    for (int i = 0; i < nedges; i++) {
        int f0 = surf.data().E(i, 0);
        int f1 = surf.data().E(i, 1);

        if (f0 != -1 && f1 != -1) {
            Eigen::Vector3d edge = (surf.data().V.row(surf.data().edgeVerts(i, 1)) - surf.data().V.row(surf.data().edgeVerts(i, 0))).transpose();
            Eigen::Vector3d v0 = face_vec[f0];
            Eigen::Vector3d v1 = face_vec[f1];
            edge_err[i] = std::abs((v1 - v0).dot(edge));
        }
    }
    return edge_err;
}

static Eigen::MatrixXd getExtVecs(const Surface& surf, const Eigen::MatrixXd& vec) {
    Eigen::MatrixXd face_vecs(surf.data().F.rows(), 3);
    for(int fid = 0; fid < surf.data().F.rows(); fid++) {
        Eigen::Vector3d v_ext = surf.data().Bs[fid] * (vec.row(fid).transpose());
        face_vecs.row(fid) = v_ext.transpose();
    }
    return face_vecs;
}

static void renderSingleVectorFields(polyscope::SurfaceMesh* ps_surf_mesh, const Surface& surf, const Eigen::MatrixXd& vec, std::string name, double plot_rescaling) {
    Eigen::MatrixXd face_vecs = getExtVecs(surf, vec) * plot_rescaling;
    ps_surf_mesh->addFaceVectorQuantity(name, face_vecs, polyscope::VectorType::AMBIENT);
}

Eigen::MatrixXd getGradient(const Surface& surf, const Eigen::MatrixXd& init_vec, const Eigen::VectorXd& scalars, double scaling) {
    int nfaces = surf.data().F.rows();
    Eigen::MatrixXd grad(nfaces, 2);

    for(int i = 0; i < nfaces; i++) {
        int v0 = surf.data().F(i, 0);
        int v1 = surf.data().F(i, 1);
        int v2 = surf.data().F(i, 2);

        Eigen::Vector2d dphi = {scalars[v1] - scalars[v0], scalars[v2] - scalars[v0]};
        Eigen::Matrix2d BTB = surf.data().Bs[i].transpose() * surf.data().Bs[i];

        Eigen::Vector2d int_dphi = BTB * init_vec.row(i).transpose();
        Eigen::Vector2i jumps;
        for(int k = 0; k < 2; k++) {
            jumps[k] = std::round((int_dphi[k] - dphi[k]) / (2 * M_PI));
            dphi[k] = jumps[k] * 2 * M_PI + dphi[k];
        }

        dphi = BTB.inverse() * dphi;
        grad.row(i) = dphi.transpose() / scaling;
    }

    std::cout << "gradient computation done" << std::endl;

    return grad;
}

Eigen::MatrixXd calculateErr(const Surface& surf, const Eigen::MatrixXd& vec, const Eigen::VectorXd& scalars, double scaling) {
    Eigen::MatrixXd err = getGradient(surf, vec, scalars, scaling);
    return err - vec;
}

static void renderErrFields(polyscope::SurfaceMesh* ps_surf_mesh, const Surface& surf, const std::vector<Eigen::MatrixXd>& vecs, const std::vector<Eigen::VectorXd>& scalars, double vec_rescaling_ratio, std::string name) {
    int nscalar = scalars.size();

    int nfaces = surf.nFaces();
    int nvecs = vecs.size();

    if(nscalar != nvecs) {
        return;
    }

    for(int i = 0; i < nvecs; i++) {
        Eigen::MatrixXd err = calculateErr(surf, vecs[i], scalars[i], vec_rescaling_ratio);
        Eigen::MatrixXd ext_err = getExtVecs(surf, err);
        Eigen::VectorXd err_norm(nfaces);
        for(int fid = 0; fid < nfaces; fid++) {
            err_norm[fid] = ext_err.row(fid).norm();
        }
        ps_surf_mesh->addFaceVectorQuantity(name + " err vec " + std::to_string(i), ext_err);
        ps_surf_mesh->addFaceScalarQuantity(name + " err " + std::to_string(i), err_norm);
    }
}

void renderIntegrabilityErr(const Eigen::MatrixXd& face_vecs, Surface& cur_surf, const IntegrabilityErrType& err_type, const std::string vec_name, const std::string mesh_name) {
    int nvecs = face_vecs.rows() / cur_surf.nFaces();
    auto surf_mesh = polyscope::getSurfaceMesh(mesh_name);

    std::vector<Eigen::VectorXd> edge_errs = checkLocalIntegrability(cur_surf, face_vecs, err_type);
    std::cout << "start to render integrability" << std::endl;
    std::cout << "edge err size: " << edge_errs.size() << std::endl;

    for(int i = 0; i < edge_errs.size(); i++) {
        auto& edge_err = edge_errs[i];
        std::cout << vec_name;
        if(err_type == kStandard) {
            std::cout << ": " << i;
        }
        std::cout << ", min: " << edge_err.minCoeff() << ", max: " << edge_err.maxCoeff();
    }

    // rendering
    for (int i = 0; i < edge_errs.size(); i++) {
        auto& edge_err = edge_errs[i];
        Eigen::VectorXd vert_err;
        vert_err.setZero(cur_surf.nVerts());
        for (int eid = 0; eid < cur_surf.nEdges(); eid++) {
            for (int j = 0; j < 2; j++) {
              int vid = cur_surf.data().edgeVerts(eid, j);
              vert_err[vid] += edge_err[eid] / 2;
            }
        }

        std::string name = "";
        switch(err_type) {
        case kStandard: {
            for (int vi = 0; vi < nvecs; vi++) {
              auto err_plot = surf_mesh->addVertexScalarQuantity(vec_name + " " + std::to_string(vi) +
                                                                     " local integrability err",
                                                                 vert_err);
              err_plot->setMapRange({vert_err.minCoeff(), vert_err.maxCoeff()});
            }
            break;
        }
        case kL2: {
            auto err_plot = surf_mesh->addVertexScalarQuantity(vec_name + " L2 local integrability err",
                                                               vert_err);
            err_plot->setMapRange({vert_err.minCoeff(), vert_err.maxCoeff()});
            break;
        }
        case kL4: {
            auto err_plot = surf_mesh->addVertexScalarQuantity(vec_name + " L4 local integrability err",
                                                               vert_err);
            err_plot->setMapRange({vert_err.minCoeff(), vert_err.maxCoeff()});
            break;
        }
        }
    }
}

Eigen::MatrixXd paintPhi(const Eigen::VectorXd& phi, Eigen::VectorXd* brightness = nullptr)      // brightness should between 0 and 1
{
    int nverts = phi.size();
    Eigen::MatrixXd color(nverts, 3);
    for (int i = 0; i < nverts; i++)
    {
        double r, g, b;
        //            double h = 360.0 * phi[i] / 2.0 / M_PI + 120;
        double h = 360.0 * phi[i] / 2.0 / M_PI;
        h = 360 + ((int)h % 360); // fix for libigl bug
        double s = 1.0;
        double v = 0.5;
        if(brightness)
        {
            double r = (*brightness)(i);
            v = r * r / (r * r + 1);
        }
        //                v = (*brightness)(i);
        igl::hsv_to_rgb(h, s, v, r, g, b);
        color(i, 0) = r;
        color(i, 1) = g;
        color(i, 2) = b;
    }
    return color;
}


IntegrationType int_type = kKnoppel;
IntegrabilityErrType int_err_type = kStandard;

Eigen::MatrixXd mesh_pts;
Eigen::MatrixXi mesh_faces;
Surface surf;
Eigen::MatrixXd poly_vecs;
Eigen::MatrixXd vector_fields;

std::vector<Eigen::MatrixXd> poly_vecs_list;
std::vector<Eigen::MatrixXd> vector_field_list;

std::vector<Eigen::VectorXd> edge_one_forms;
std::vector<Eigen::VectorXd> scalar_fields;

double rescale_ratio = 1;
int up_level = 0;

double vec_plot_scaling = 1;

void loadMesh(const std::string& mesh_path, Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
  if(!igl::readOBJ(mesh_path, V, F)) {
    std::string new_path = igl::file_dialog_open();
    igl::readOBJ(new_path, V, F);
  }
}

void renderVectorFields() {
    int nvecs = vector_field_list.size();
    auto surf_mesh = polyscope::getSurfaceMesh("mesh");

    for (int i = 0; i < nvecs; i++) {
        renderSingleVectorFields(surf_mesh, surf, vector_field_list[i], "vec" + std::to_string(i), vec_plot_scaling);
    }
}

void renderPolyFields() {
    int nvecs = poly_vecs_list.size();
    auto surf_mesh = polyscope::getSurfaceMesh("mesh");

    for (int i = 0; i < nvecs; i++) {
        renderSingleVectorFields(surf_mesh, surf, poly_vecs_list[i], "poly vec " + std::to_string(i) + " part 1", vec_plot_scaling);
        renderSingleVectorFields(surf_mesh, surf, -poly_vecs_list[i], "poly vec " + std::to_string(i) + " part 2", vec_plot_scaling);
    }
}

void renderScalarFields() {
  int nscalar = scalar_fields.size();

  if(up_level) {
    std::vector<std::pair<int, Eigen::Vector3d>> bary;
    Eigen::MatrixXd NV;
    Eigen::MatrixXi NF;
   SurfaceFields::upsampleSurface(mesh_pts, mesh_faces, up_level, NV, NF, &bary);
   auto up_surf_mesh =  polyscope::registerSurfaceMesh("up mesh", NV, NF);

   for(int i = 0; i < nscalar; i++) {
      Eigen::VectorXd up_theta;
      SurfaceFields::upsampleScalarFields(surf, edge_one_forms[i], scalar_fields[i], bary, up_theta);
      Eigen::MatrixXd phi = paintPhi(up_theta);
      up_surf_mesh->addVertexColorQuantity("Scalar " + std::to_string(i), phi);
    }

  } else {
    auto surf_mesh = polyscope::getSurfaceMesh("mesh");

    for(int i = 0; i < nscalar; i++) {
      Eigen::MatrixXd phi = paintPhi(scalar_fields[i]);
      surf_mesh->addVertexColorQuantity("Scalar " + std::to_string(i), phi);
    }
  }
}


void loadVectorFields() {
  Eigen::MatrixXd A;
  std::string load_vec_name = igl::file_dialog_open();
  if (!deserializeMatrixNew(A, load_vec_name)) {
    deserializeMatrix(A, load_vec_name);
  }

  std::cout << "start to register mesh" << std::endl;
  auto surf_mesh = polyscope::registerSurfaceMesh("mesh", mesh_pts, mesh_faces);
  surf_mesh->setEnabled(true);
  std::cout << "register finished" << std::endl;

  double ave_edge_len = surf.data().averageEdgeLength;

  // we need to turn this into intrinsic surface poly fields
  if (mesh_faces.rows() != A.rows()) {
    std::cout << "vector field size and nfaces do not match" << std::endl;
  } else {
    int npoly = A.cols() / 2;
    poly_vecs_list.resize(npoly, Eigen::MatrixXd::Zero(mesh_faces.rows(), 2));
    poly_vecs.resize(npoly * mesh_faces.rows(), 2);
    double max_vec_len = 0;
    for (int i = 0; i < mesh_faces.rows(); i++) {
      for (int j = 0; j < npoly; j++) {
              Eigen::Vector3d v;
              v << A(i, npoly * j + 0), A(i, npoly * j + 1), 0;
              max_vec_len = std::max(max_vec_len, v.norm());
              Eigen::Matrix<double, 3, 2> B = surf.data().Bs[i];
              Eigen::Vector2d newvec =
                  (B.transpose() * B).inverse() * B.transpose() * v;
              poly_vecs.row(npoly * i + j) = newvec.transpose();
              poly_vecs_list[j].row(i) = newvec.transpose();
      }
    }
    std::cout << "load finished" << std::endl;
    vec_plot_scaling = ave_edge_len / max_vec_len;

    vector_fields = SurfaceFields::CombPolyVectors(surf, poly_vecs);

    int nvecs = vector_fields.rows() / mesh_faces.rows();

    for (int i = 0; i < nvecs; i++) {
      Eigen::MatrixXd face_vecs(mesh_faces.rows(), 2);
      for (int fid = 0; fid < mesh_faces.rows(); fid++) {
              face_vecs.row(fid) << vector_fields(nvecs * fid + i, 0),
                  vector_fields(nvecs * fid + i, 1);
      }
      vector_field_list.push_back(face_vecs);
    }

    renderPolyFields();
    renderVectorFields();
  }
}


void load() {
  std::string load_file_name = igl::file_dialog_open();
  std::cout << "load file in: " << load_file_name << std::endl;
  loadMesh(load_file_name, mesh_pts, mesh_faces);
  surf = Surface(mesh_pts, mesh_faces);
  std::cout << "loading finished, #V: " << mesh_pts.rows() << ", #F: " << mesh_faces.rows() << std::endl;
  loadVectorFields();
}

void myCallback() {
  if (ImGui::Button("Load")) {
    load();
  }
  ImGui::SameLine();
  if (ImGui::Button("Load vector fields")) {
    loadVectorFields();
  }

  if (ImGui::Button("Get Simple PolyVectors"))
  {
    int nfaces = mesh_faces.rows();
    poly_vecs.resize(2 * nfaces, 2);
    for(int i = 0; i < nfaces; i++) {
      poly_vecs.row(2 * i) << 1, 0;
      poly_vecs.row(2 * i + 1) << 0, 1;
    }

    vector_fields = SurfaceFields::CombPolyVectors(surf, poly_vecs);
    renderPolyFields();
    renderVectorFields();
  }

  if (ImGui::Button("Get Principal PolyVectors"))
  {
    Eigen::VectorXd PV1, PV2;
    Eigen::MatrixXd PD1, PD2;

    igl::principal_curvature(mesh_pts, mesh_faces, PD1, PD2, PV1, PV2);
    auto surf_mesh = polyscope::getSurfaceMesh("mesh");
    surf_mesh->addVertexVectorQuantity("PD1", PD1);
    surf_mesh->addVertexVectorQuantity("PD2", PD2);
    surf_mesh->addVertexScalarQuantity("PV1", PV1);
    surf_mesh->addVertexScalarQuantity("PV2", PV2);

    int nfaces = mesh_faces.rows();

    poly_vecs.resize(nfaces, 2);
    for(int i = 0; i < nfaces; i++) {
      int v0 = mesh_faces(i, 0);
      int v1 = mesh_faces(i, 1);
      int v2 = mesh_faces(i, 2);

      Eigen::Vector2d tmp_v = {PV1[v1] - PV1[v0], PV1[v2] - PV1[v0]};
      Eigen::Matrix2d BTB = surf.data().Bs[i].transpose() * surf.data().Bs[i];
      tmp_v = BTB.inverse() * tmp_v;
      poly_vecs.row(i) = tmp_v.transpose();

//      tmp_v = {PV2[v1] - PV2[v0], PV2[v2] - PV2[v0]};
//      tmp_v = BTB.inverse() * tmp_v;
//      poly_vecs.row(2 * i + 1) = tmp_v.transpose();
    }

    vector_fields = SurfaceFields::CombPolyVectors(surf, poly_vecs);
    renderPolyFields();
    renderVectorFields();
  }

  if(ImGui::InputDouble("Rescaling ratio", &rescale_ratio)) {}

  if(ImGui::InputInt("Upsampling Level", &up_level)) {
    if(up_level < 0) {
      up_level = 0;
    }
    if(edge_one_forms.size() == scalar_fields.size() && !edge_one_forms.empty()) {
      renderScalarFields();
    }
  }

  ImGui::Combo("Integrability Error Type", (int*)&int_err_type, "Standard\0L2\0L4\0");

  if (ImGui::Button("Check Integrability Before Comb")) {
    renderIntegrabilityErr(poly_vecs, surf, int_err_type, "input vec", "mesh");
  }

  if (ImGui::Button("Check Integrability")) {
    renderIntegrabilityErr(vector_fields, surf, int_err_type, "combed vec", "mesh");
  }

  ImGui::Combo("Integration Method", (int*)&int_type, "Bomme\0Knoppel\0");

  if (ImGui::Button("Integrate Vector Fields")) {
    scalar_fields.clear();
    edge_one_forms.clear();

    int nvecs = vector_fields.rows() / mesh_faces.rows();
    for(int i = 0; i < nvecs; i++) {
      Eigen::MatrixXd vec(mesh_faces.rows(), 2);
      for(int j = 0; j < mesh_faces.rows(); j++) {
        vec.row(j) << vector_fields(nvecs * j + i, 0), vector_fields(nvecs * j + i, 1);
      }
      vec *= rescale_ratio;


      std::unique_ptr<SurfaceFields::GlobalFieldIntegration> int_ptr;
      switch (int_type) {
      case kBomme:
      {
          int_ptr = std::make_unique<SurfaceFields::MIGlobalIntegration>(10, 10);
          break;
      }
      case kKnoppel:
      {
          int_ptr = std::make_unique<SurfaceFields::StripePatternsGlobalIntegration>();
          break;
      }
      default:
      {
          int_ptr = std::make_unique<SurfaceFields::StripePatternsGlobalIntegration>();
          break;
      }
      }
      Eigen::VectorXd theta, edge_omega;
      int_ptr->globallyIntegrateOneComponent(surf, vec, theta, &edge_omega);
      scalar_fields.emplace_back(theta);
      edge_one_forms.emplace_back(edge_omega);
    }
    renderScalarFields();
  }

  if(ImGui::Button("Get grad phi")) {
    int nvecs = vector_fields.rows() / surf.nFaces();
    for(int i = 0; i < scalar_fields.size(); i++) {
      Eigen::MatrixXd face_vec(surf.nFaces(), 2);
      for(int fid = 0; fid < surf.nFaces(); fid++) {
          face_vec.row(fid) = vector_fields.row(nvecs * fid + i);
      }

      Eigen::MatrixXd dphi = getGradient(surf, face_vec, scalar_fields[i], rescale_ratio);

      auto surf_mesh = polyscope::getSurfaceMesh("mesh");

      renderSingleVectorFields(surf_mesh, surf, dphi, "dphi " + std::to_string(i), vec_plot_scaling);
    }
  }

  if (ImGui::Button("Integrate Uncombed Vector Fields")) {
      scalar_fields.clear();
      edge_one_forms.clear();

      int nvecs = poly_vecs.rows() / mesh_faces.rows();
      for (int i = 0; i < nvecs; i++) {
          Eigen::MatrixXd vec(mesh_faces.rows(), 2);
          for (int j = 0; j < mesh_faces.rows(); j++) {
              vec.row(j) << poly_vecs(nvecs * j + i, 0), poly_vecs(nvecs * j + i, 1);
          }
          vec *= rescale_ratio;

          std::unique_ptr<SurfaceFields::StripePatternsGlobalIntegration> ptr = std::make_unique<SurfaceFields::StripePatternsGlobalIntegration>();
          Eigen::VectorXd theta, edge_omega;
          ptr->globallyIntegrateOneComponent(surf, vec, theta, &edge_omega);
          scalar_fields.emplace_back(theta);
          edge_one_forms.emplace_back(edge_omega);
      }
      renderScalarFields();
  }

  if(ImGui::Button("Get grad uncombed phi")) {
      int nvecs = poly_vecs.rows() / surf.nFaces();
      for(int i = 0; i < scalar_fields.size(); i++) {
          Eigen::MatrixXd face_vec(surf.nFaces(), 2);
          for(int fid = 0; fid < surf.nFaces(); fid++) {
            face_vec.row(fid) = poly_vecs.row(nvecs * fid + i);
          }

          Eigen::MatrixXd dphi = getGradient(surf, face_vec, scalar_fields[i], rescale_ratio);

          auto surf_mesh = polyscope::getSurfaceMesh("mesh");
          renderSingleVectorFields(surf_mesh, surf, dphi, "d_uncombed_phi " + std::to_string(i), vec_plot_scaling);
      }
  }

  if (ImGui::Button("Round Vector Fields")) {
      Eigen::MatrixXd render_QCover, grad_theta, face_vectors;
      Eigen::MatrixXi render_FCover;
      Eigen::VectorXd render_theta, err;

      std::vector<Eigen::Vector3d> cut_pts;
      std::vector<std::vector<int>> cut_edges;

      std::unique_ptr<SurfaceFields::GlobalFieldIntegration> int_ptr;
      switch (int_type) {
      case kBomme:
      {
          int_ptr = std::make_unique<SurfaceFields::MIGlobalIntegration>(1.0, 1e-4);
          break;
      } 
      case kKnoppel:
      {
          int_ptr = std::make_unique<SurfaceFields::StripePatternsGlobalIntegration>();
          break;
      }
      default:
      {
          int_ptr = std::make_unique<SurfaceFields::StripePatternsGlobalIntegration>();
          break;
      }
      }

      roundVectorFields(mesh_pts, mesh_faces, poly_vecs, render_QCover, render_FCover, render_theta, int_ptr.get(), &grad_theta, rescale_ratio, &face_vectors, &cut_pts, &cut_edges);
      err.setZero(face_vectors.rows());
      Eigen::MatrixXd err_vec = face_vectors - grad_theta;
      for(int i = 0; i < face_vectors.rows(); i++) {
          err[i] = err_vec.row(i).norm();
      }

      auto cover_surface = polyscope::registerSurfaceMesh("splitted cover mesh", render_QCover, render_FCover);
      cover_surface->addVertexColorQuantity("phi", paintPhi(render_theta));
      cover_surface->addFaceVectorQuantity("vecs", face_vectors * vec_plot_scaling, polyscope::VectorType::AMBIENT);
      cover_surface->addFaceVectorQuantity("grad phi", grad_theta * vec_plot_scaling, polyscope::VectorType::AMBIENT);
      cover_surface->addFaceVectorQuantity("err vec", err_vec * vec_plot_scaling, polyscope::VectorType::AMBIENT);
      cover_surface->addFaceScalarQuantity("grad phi error", err);

      std::cout << "global rescaling ratio: " << vec_plot_scaling << std::endl;

      auto cut_poly = polyscope::registerCurveNetwork("cuts", cut_pts, cut_edges);

  }

}

int main(int argc, char** argv)
{
  Eigen::MatrixXd test_V(4, 3);
  Eigen::MatrixXi test_F(2, 3);

  test_V << 0, 0, 0,
            1, 0, 0,
            1, 1, 0,
            0, 1, 0;
  test_F << 0, 1, 2,
            0, 2, 3;

  Eigen::MatrixXd test_vec(2, 3);
  test_vec << M_PI, 0, 0,
              M_PI, 0, 0;
  Eigen::VectorXd test_theta(4);

  test_theta << 0, M_PI, M_PI, 0;

  Surface test_surf(test_V, test_F);

  Eigen::MatrixXd test_vec_bary(2, 2);
  for(int i = 0; i < 2; i++) {
      Eigen::Matrix<double, 3, 2> B = test_surf.data().Bs[i];
      test_vec_bary.row(i) = ((B.transpose() * B).inverse() * B.transpose() * (test_vec.row(i)).transpose()).transpose();
  }

  std::unique_ptr<SurfaceFields::MIGlobalIntegration> ptr = std::make_unique<SurfaceFields::MIGlobalIntegration>(0 ,0);
  ptr->globallyIntegrateOneComponent(test_surf, test_vec_bary, test_theta, nullptr);




  // Initialize polyscope
  polyscope::init();

  // Set the callback function
  polyscope::state::userCallback = myCallback;
  load();

  polyscope::view::upDir = polyscope::view::UpDir::ZUp;

  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}