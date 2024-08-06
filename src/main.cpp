#include "polyscope/polyscope.h"

#include <igl/PI.h>
#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/boundary_loop.h>
#include <igl/exact_geodesic.h>
#include <igl/gaussian_curvature.h>
#include <igl/invert_diag.h>
#include <igl/lscm.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/readOBJ.h>

#include "polyscope/messages.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"

#include <iostream>
#include <unordered_set>
#include <utility>

#include "args/args.hxx"
#include "json/json.hpp"

// The mesh, Eigen representation
Eigen::MatrixXd meshV1;
Eigen::MatrixXi meshF1;
Eigen::MatrixXd meshV2;
Eigen::MatrixXi meshF2;

// Options for algorithms
int iVertexSource = 7;

void addCurvatureScalar() {
  using namespace Eigen;
  using namespace std;

  VectorXd K;
  igl::gaussian_curvature(meshV1, meshF1, K);
  SparseMatrix<double> M, Minv;
  igl::massmatrix(meshV1, meshF1, igl::MASSMATRIX_TYPE_DEFAULT, M);
  igl::invert_diag(M, Minv);
  K = (Minv * K).eval();

  polyscope::getSurfaceMesh("input mesh")
      ->addVertexScalarQuantity("gaussian curvature", K,
                                polyscope::DataType::SYMMETRIC);
}

void computeDistanceFrom() {
  Eigen::VectorXi VS, FS, VT, FT;
  // The selected vertex is the source
  VS.resize(1);
  VS << iVertexSource;
  // All vertices are the targets
  VT.setLinSpaced(meshV1.rows(), 0, meshV1.rows() - 1);
  Eigen::VectorXd d;
  igl::exact_geodesic(meshV1, meshF1, VS, FS, VT, FT, d);

  polyscope::getSurfaceMesh("input mesh")
      ->addVertexDistanceQuantity(
          "distance from vertex " + std::to_string(iVertexSource), d);
}

// void computeParameterization() {
//   using namespace Eigen;
//   using namespace std;
// 
//   // Fix two points on the boundary
//   VectorXi bnd, b(2, 1);
//   igl::boundary_loop(meshF1, bnd);
// 
//   if (bnd.size() == 0) {
//     polyscope::warning("mesh has no boundary, cannot parameterize");
//     return;
//   }
// 
//   b(0) = bnd(0);
//   b(1) = bnd(round(bnd.size() / 2));
//   MatrixXd bc(2, 2);
//   bc << 0, 0, 1, 0;
// 
//   // LSCM parametrization
//   Eigen::MatrixXd V_uv;
//   igl::lscm(meshV1, meshF1, b, bc, V_uv);
// 
//   polyscope::getSurfaceMesh("input mesh")
//       ->addVertexParameterizationQuantity("LSCM parameterization", V_uv);
// }

void computeNormals() {
  Eigen::MatrixXd N_vertices;
  igl::per_vertex_normals(meshV1, meshF1, N_vertices);

  polyscope::getSurfaceMesh("input mesh")
      ->addVertexVectorQuantity("libIGL vertex normals", N_vertices);
}

void callback() {

  static int numPoints = 2000;
  static float param = 3.14;

  ImGui::PushItemWidth(100);

  // Curvature
  // if (ImGui::Button("add curvature")) {
  //   addCurvatureScalar();
  // }
  
  // Normals 
  // if (ImGui::Button("add normals")) {
  //   computeNormals();
  // }

  // Param
  // if (ImGui::Button("add parameterization")) {
  //  computeParameterization();
  //}

  // Geodesics
  // if (ImGui::Button("compute distance")) {
  //   computeDistanceFrom();
  // }
  ImGui::SameLine();
  ImGui::InputInt("source vertex", &iVertexSource);

  ImGui::PopItemWidth();
}

int main(int argc, char **argv) {
  // Configure the argument parser
  args::ArgumentParser parser("A simple demo of Polyscope with libIGL.\nBy "
                              "Nick Sharp (nsharp@cs.cmu.edu)",
                              "");
  // args::Positional<std::string> inFile(parser, "mesh", "input mesh");
  args::PositionalList<std::string> inFiles(parser, "meshes", "input mesh files");

  // Parse args
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError e) {
    std::cerr << e.what() << std::endl;

    std::cerr << parser;
    return 1;
  }

  // Options
  polyscope::options::autocenterStructures = true;
  polyscope::view::windowWidth = 1024;
  polyscope::view::windowHeight = 1024;

  // Initialize polyscope
  polyscope::init();

  // std::string filename = args::get(inFile);
  // std::cout << "loading: " << filename << std::endl;

  // Read the mesh
  // igl::readOBJ(filename, meshV1, meshF1);

  // Register the mesh with Polyscope
  // polyscope::registerSurfaceMesh("input mesh", meshV1, meshF1);

  // Read the first mesh
  std::string filename1 = inFiles.Get()[0];
  std::cout << "loading: " << filename1 << std::endl;
  igl::readOBJ(filename1, meshV1, meshF1);
  polyscope::registerSurfaceMesh("input mesh 1", meshV1, meshF1);

  // Read the second mesh
  std::string filename2 = inFiles.Get()[1];
  std::cout << "loading: " << filename2 << std::endl;
  igl::readOBJ(filename2, meshV2, meshF2);
  polyscope::registerSurfaceMesh("input mesh 2", meshV2, meshF2);

  // Add the callback
  polyscope::state::userCallback = callback;

  // Show the gui
  polyscope::show();

  return 0;
}
