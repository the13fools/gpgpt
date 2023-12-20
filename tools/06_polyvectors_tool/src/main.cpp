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

#include <sstream>

#include "../src/SurfaceFields/VectorFields2PolyVectors.h"

Eigen::MatrixXd mesh_pts;
Eigen::MatrixXi mesh_faces;
Surface surf;
Eigen::MatrixXd poly_vecs;
Eigen::MatrixXd vector_fields;

void loadMesh(const std::string& mesh_path, Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
  if(!igl::readOBJ(mesh_path, V, F)) {
    std::string new_path = igl::file_dialog_open();
    igl::readOBJ(new_path, V, F);
  }
}

void load() {
  std::string load_file_name = igl::file_dialog_open();
  std::cout << "load file in: " << load_file_name << std::endl;
  loadMesh(load_file_name, mesh_pts, mesh_faces);
  surf = Surface(mesh_pts, mesh_faces);
  std::cout << "loading finished, #V: " << mesh_pts.rows() << ", #F: " << mesh_faces.rows() << std::endl;

  std::cout << "start to register mesh" << std::endl;
  auto surf_mesh = polyscope::registerSurfaceMesh("mesh", mesh_pts, mesh_faces);
  surf_mesh->setEnabled(true);
  std::cout << "register finished" << std::endl;
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


void myCallback() {
  ImGui::PushItemWidth(100);
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