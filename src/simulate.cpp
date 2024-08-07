#include <iostream>

#include "prettyprint.hpp"

#include "igl/edges.h"
#include "igl/readOBJ.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

void RunSimulation() {
  ;
}

void CallbackFunction() {

  ImGui::PushItemWidth(100);

  if (ImGui::Button("Run simulation")) {
    RunSimulation();
  }

  ImGui::SameLine();
  ImGui::PopItemWidth();
}

int main(int argc, char **argv) {

  Eigen::MatrixXd meshV; // V x 3
  Eigen::MatrixXi meshF; // F x 3

  std::string baseDir = "../data/";

  std::string meshfilename = "plane.obj";

  // Concatenate the base directory path and the filename
  std::string meshPath = baseDir + meshfilename;
  // Read the first mesh
  igl::readOBJ(meshPath, meshV, meshF);

  // Initialize polyscope with some options
  polyscope::view::setUpDir(polyscope::UpDir::ZUp);
  polyscope::view::setFrontDir(polyscope::FrontDir::XFront);
  polyscope::init();

  // Register the mesh with Polyscope
  polyscope::registerSurfaceMesh("Cloth", meshV, meshF);

  // Specify the callback
  polyscope::state::userCallback = CallbackFunction;

  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
