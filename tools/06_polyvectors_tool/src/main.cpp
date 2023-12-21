#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

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

#include <sstream>

#include "../src/SurfaceFields/VectorFields2PolyVectors.h"
#include "../src/SurfaceFields/StripePatternIntegration.h"
#include "../src/SurfaceFields/UpsampleSurfaceFields.h"


Eigen::MatrixXd mesh_pts;
Eigen::MatrixXi mesh_faces;
Surface surf;
Eigen::MatrixXd poly_vecs;
Eigen::MatrixXd vector_fields;

std::vector<Eigen::VectorXd> edge_one_forms;
std::vector<Eigen::VectorXd> scalar_fields;

double rescale_ratio = 1;
int up_level = 0;

void loadMesh(const std::string& mesh_path, Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
  if(!igl::readOBJ(mesh_path, V, F)) {
    std::string new_path = igl::file_dialog_open();
    igl::readOBJ(new_path, V, F);
  }
}

// Deserialize Eigen matrix from a binary file
bool deserializeMatrix(Eigen::MatrixXd& mat, const std::string& filepath) {
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

void renderVectorFields()
{
  int nfaces = mesh_faces.rows();
  int nvecs = vector_fields.rows() / nfaces;
  auto surf_mesh = polyscope::getSurfaceMesh("mesh");

  for(int i = 0; i < nvecs; i++) {
    std::vector<Eigen::Vector3d> face_vecs;
    for(int fid = 0; fid < nfaces; fid++) {
      Eigen::Vector2d v;
      v << vector_fields(nvecs * fid + i, 0), vector_fields(nvecs * fid + i, 1);
      Eigen::Vector3d v_ext = surf.data().Bs[fid] * v;
      face_vecs.push_back(v_ext);
    }
    surf_mesh->addFaceVectorQuantity("vec " + std::to_string(i), face_vecs);
  }
}

void renderPolyFields() {
  int nfaces = mesh_faces.rows();
  int nvecs = poly_vecs.rows() / nfaces;
  auto surf_mesh = polyscope::getSurfaceMesh("mesh");

  for(int i = 0; i < nvecs; i++) {
    std::vector<Eigen::Vector3d> face_vecs_1, face_vecs_2;
    for(int fid = 0; fid < nfaces; fid++) {
      Eigen::Vector2d v;
      v << poly_vecs(nvecs * fid + i, 0), poly_vecs(nvecs * fid + i, 1);
      Eigen::Vector3d v_ext = surf.data().Bs[fid] * v;
      face_vecs_1.push_back(v_ext);
      face_vecs_2.push_back(-v_ext);
    }
    surf_mesh->addFaceVectorQuantity("poly vec " + std::to_string(i) + " part 1", face_vecs_1);
    surf_mesh->addFaceVectorQuantity("poly vec " + std::to_string(i) + " part 2", face_vecs_2);
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



void load() {
  std::string load_file_name = igl::file_dialog_open();
  std::cout << "load file in: " << load_file_name << std::endl;
  loadMesh(load_file_name, mesh_pts, mesh_faces);
  surf = Surface(mesh_pts, mesh_faces);
  std::cout << "loading finished, #V: " << mesh_pts.rows() << ", #F: " << mesh_faces.rows() << std::endl;

  Eigen::MatrixXd A;
  std::string load_vec_name = igl::file_dialog_open();
  deserializeMatrix(A, load_vec_name);

  std::cout << "start to register mesh" << std::endl;
  auto surf_mesh = polyscope::registerSurfaceMesh("mesh", mesh_pts, mesh_faces);
  surf_mesh->setEnabled(true);
  std::cout << "register finished" << std::endl;

  // we need to turn this into intrinsic surface poly fields
  if(mesh_faces.rows() != A.rows()) {
    std::cout << "vector field size and nfaces do not match" << std::endl;
  } else {
    int npoly = A.cols() / 2;
    poly_vecs.resize(npoly * mesh_faces.rows(), 2);
    for(int i = 0; i < mesh_faces.rows(); i++) {
      for(int j = 0; j < npoly; j++) {
        Eigen::Vector3d v;
        v << A(i, npoly * j + 0), A(i, npoly * j + 1), 0;
        Eigen::Matrix<double, 3, 2> B = surf.data().Bs[i];
        Eigen::Vector2d newvec = (B.transpose()*B).inverse() * B.transpose() * v;
        poly_vecs.row(npoly * i + j) = newvec.transpose();
      }
    }
    std::cout << "load finished" << std::endl;

    vector_fields = SurfaceFields::CombPolyVectors(surf, poly_vecs);
    renderPolyFields();
    renderVectorFields();
  }
}



void myCallback() {
  if (ImGui::Button("Load")) {
    load();
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

      std::unique_ptr<SurfaceFields::StripePatternsGlobalIntegration> ptr = std::make_unique<SurfaceFields::StripePatternsGlobalIntegration>();
      Eigen::VectorXd theta, edge_omega;
      ptr->globallyIntegrateOneComponent(surf, vec, theta, &edge_omega);
      scalar_fields.emplace_back(theta);
      edge_one_forms.emplace_back(edge_omega);
    }
    renderScalarFields();
  }

}

int main(int argc, char** argv)
{
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