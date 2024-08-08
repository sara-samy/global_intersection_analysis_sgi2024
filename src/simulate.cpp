#include <Eigen/Dense>
#include <iostream>

#include "igl/doublearea.h"
#include "igl/edge_lengths.h"
#include "igl/edges.h"
#include "igl/readOBJ.h"
#include "igl/triangle_triangle_adjacency.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::VectorXi;

class Cloth {
public:
  int numParticles;
  int numTris;

  MatrixXd pos;
  MatrixXd prevPos;
  MatrixXd restPos;
  MatrixXd vel;

  // will be used later to define the constraints.
  MatrixXi edgeIds;  // #E x 2
  MatrixXi triPairs; // #F x 3

  VectorXd invMass;
  VectorXd bendingLengths;
  VectorXd stretchingLengths;

  float stretchingCompliance;
  float bendingCompliance;

  // will store gradient of constraint functions later.
  VectorXd grads;

  void initPhysics(const MatrixXi &F);

  Cloth(const MatrixXd &V, const MatrixXi &F, float bendingCompliance);
};
Cloth::Cloth(const MatrixXd &V, const MatrixXi &F, float bendingCompliance) {
  numParticles = V.rows();
  numTris = F.rows();
  pos = V;
  prevPos = pos;
  restPos = pos;
  vel = MatrixXd::Zero(numParticles, 3);

  std::cout << "# Triangles in cloth mesh: " << numTris << std::endl;
  std::cout << "# Particles in cloth mesh: " << numParticles << std::endl;

  // Get edges
  igl::edges(F, edgeIds);

  // Get neighbouring triangle pairs
  // The (i, j) entry contains id of triangle adjacent to triangle i w.r.t edge
  // j. If the entry is -1 then there are no adjacent triangles w.r.t that edge.
  igl::triangle_triangle_adjacency(F, triPairs);

  initPhysics(F);
}
void Cloth::initPhysics(const MatrixXi &F) {
  // Compute the edge lengths
  igl::edge_lengths(this->pos, F, stretchingLengths);
  std::cout << "stretchingLengths has size " << stretchingLengths.rows()
            << " x " << stretchingLengths.cols() << std::endl;

  // Compute inverse masses of particles
  invMass = VectorXd::Zero(numParticles);
  // Vector to store the double areas
  VectorXd triAreas;
  igl::doublearea(this->pos, F, triAreas);
  for (int i = 0; i < F.rows(); ++i) {
    double area = triAreas(i) / 2.0;
    // Distribute the area equally to the three vertices of the triangle
    double vertexArea = area / 3.0;
    for (int j = 0; j < 3; ++j) {
      int vertex_index = F(i, j);
      invMass(vertex_index) += 1.0 / vertexArea;
    }
  }
  std::cout << "invMass has size " << invMass.rows() << " x " << invMass.cols()
            << std::endl;
}

void RunSimulation() { ; }

void CallbackFunction() {

  ImGui::PushItemWidth(100);

  if (ImGui::Button("Run simulation")) {
    RunSimulation();
  }

  ImGui::SameLine();
  ImGui::PopItemWidth();
}

int main(int argc, char **argv) {

  MatrixXd meshV; // #V x 3
  MatrixXi meshF; // #F x 3

  std::string baseDir = "../data/";
  std::string meshfilename = "plane.obj";
  std::string meshPath = baseDir + meshfilename;

  // Read the mesh
  igl::readOBJ(meshPath, meshV, meshF);

  Cloth cloth(meshV, meshF, 1.0f);

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
