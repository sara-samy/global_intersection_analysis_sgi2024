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

  MatrixXd initialVertices;
  MatrixXd pos;
  MatrixXd prevPos;
  MatrixXd restPos;
  MatrixXd vel;

  // will be used later to define the constraints.
  MatrixXi edgeIds;  // #E x 2
  MatrixXi triPairs; // #F x 3

  VectorXd invMass;
  VectorXd bendingLengths;
  MatrixXd stretchingLengths;

  float stretchingCompliance;
  float bendingCompliance;

  float thickness;

  bool handleCollisions;

  // will store gradient of constraint functions later.
  VectorXd grads;
  

  void initPhysics(const MatrixXi &F);
  void Simulate(double frameDt, int numSubSteps, Eigen::Vector3d gravity);
  void solveConstraints(double dt);
  void Cloth::solveGroundCollisions();
  void Cloth::solveCollisions(double dt);


  Cloth(const MatrixXd &V, const MatrixXi &F, float bendingCompliance);
};
Cloth::Cloth(const MatrixXd &V, const MatrixXi &F, float bendingCompliance) {


	numParticles = V.rows();
  numTris = F.rows();
  pos = V;
  prevPos = pos;
  restPos = pos;
  vel = MatrixXd::Zero(numParticles, 3);
  initialVertices = V;

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

	thickness = 0.01f;
	handleCollisions = true;

  // Compute the edge lengths
	
  igl::edge_lengths(pos, F, stretchingLengths);
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

void Cloth::Simulate(double frameDt, int numSubSteps, Eigen::Vector3d gravity)
{

	double dt = frameDt / numSubSteps;
	double maxVelocity = 0.2 * thickness / dt;

	
	if (handleCollisions) {
		//hash.create(pos);
		double maxTravelDist = maxVelocity * frameDt;
		//hash.queryAll(pos, maxTravelDist);
	}
	
	for (int step = 0; step < numSubSteps; step++) 
	{
		// integrate 
		for (int i = 0; i < numParticles; i++) 
		{
			if (invMass(i) > 0.0) {

				Eigen::Vector3d velocity = vel.row(i);
				velocity += gravity;
				double vNorm = velocity.norm();
				double maxV = 0.2 * thickness / dt;
				if (vNorm > maxV) {
					velocity=velocity * maxV / vNorm;
				}
				prevPos.row(i) = pos.row(i);
				pos.row(i) += velocity*dt;
			}
		}

		// solve

		solveGroundCollisions();

		solveConstraints(dt);
		if (handleCollisions)
			solveCollisions(dt);

		// update velocities

		for (int i = 0; i < numParticles; i++) {
			if (invMass(i) > 0.0)
			{
				vel.row(i) = (pos.row(i) - prevPos.row(i)) / dt;
			}
		}
	}
}

void Cloth::solveConstraints(double dt) {
	// Iterate over each constraint
	for (int i = 0; i < edgeIds.rows(); i++) {
		int id0 = edgeIds(i, 0);
		int id1 = edgeIds(i, 1);

		double w0 = invMass(id0);
		double w1 = invMass(id1);
		double w = w0 + w1;
		if (w == 0.0)
			continue;

		// Compute the vector difference between the two points
		Eigen::Vector3d vec = pos.row(id0) - pos.row(id1);
		double len = vec.norm();
		if (len == 0.0)
			continue;

		// Normalize the vector
		vec /= len;

		double restLen = (initialVertices.row(id0) - initialVertices.row(id1)).norm();
		double C = len - restLen;
		double alpha = stretchingCompliance / (dt * dt);
		double s = -C / (w + alpha);

		// Update positions
		pos.row(id0) += s * w0 * vec;
		pos.row(id1) -= s * w1 * vec;
	}
}
void Cloth::solveGroundCollisions() {
	for (int i = 0; i < numParticles; i++) {
		if (invMass(i) == 0.0)
			continue;

		// Check the y-coordinate of the particle
		double y = pos(i, 1);
		if (y < 0.5 * thickness) {
			// Apply damping
			double damping = 1.0;
			Eigen::Vector3d displacement = pos.row(i) - prevPos.row(i);
			pos.row(i) += -damping * displacement;

			// Ensure the particle stays above the ground plane
			pos(i, 1) = 0.5 * thickness;
		}
	}
}

void Cloth::solveCollisions(double dt) {
	double thickness2 = thickness * thickness;

	for (int i = 0; i < numParticles; i++) {
		if (invMass(i) == 0.0)
			continue;

		int id0 = i;
		int first = 0;//hash.firstAdjId[i];
		int last = 0;// hash.firstAdjId[i + 1];

		for (int j = first; j < last; j++) {
			int id1 = 0;// hash.adjIds[j];
			if (invMass(id1) == 0.0)
				continue;

			// Vector difference between particles
			Eigen::Vector3d vec = pos.row(id1) - pos.row(id0);
			double dist2 = vec.squaredNorm();

			if (dist2 > thickness2 || dist2 == 0.0)
				continue;

			double restDist2 = (restPos.row(id0) - restPos.row(id1)).squaredNorm();

			double minDist = thickness;
			if (dist2 > restDist2)
				continue;

			if (restDist2 < thickness2)
				minDist = sqrt(restDist2);

			// Position correction
			double dist = sqrt(dist2);
			vec *= (minDist - dist) / dist;
			pos.row(id0) -= 0.5 * vec;
			pos.row(id1) += 0.5 * vec;

			// Velocity corrections
			Eigen::Vector3d vel0 = pos.row(id0) - prevPos.row(id0);
			Eigen::Vector3d vel1 = pos.row(id1) - prevPos.row(id1);

			// Average velocity
			Eigen::Vector3d avgVel = 0.5 * (vel0 + vel1);

			// Correct velocities
			vel0 -= avgVel;
			vel1 -= avgVel;

			// Apply corrections 
			double friction = 0.0;
			pos.row(id0) += friction * vel0;
			pos.row(id1) += friction * vel1;
		}
	}
}



void RunSimulation() { ; }
bool run = false;
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
  //polyscope::state::userCallback = CallbackFunction;

  Eigen::Vector3d gravity(0,-9.8,0);
  double dt = 0.01;
  int subSteps = 5;
  auto polyscope_callback = [&]() mutable
  {

	  ImGui::Begin("Simulator");


	  if (ImGui::Button(run?"Stop simulation":"Run simulation")) {
		  run = !run;
	  }
	  ImGui::InputDouble("dt", &dt);
	  ImGui::SliderInt("Substeps", &subSteps, 1, 20);
	  ImGui::SliderFloat("Bending compliance", &cloth.bendingCompliance, 0, 1);
	  ImGui::SliderFloat("Stretching compliance", &cloth.stretchingCompliance, 0, 1);

	  if (run)
	  {
	  	//cloth.Simulate(dt, subSteps, gravity);
	  	//polyscope::getSurfaceMesh("Cloth")->updateVertexPositions(cloth.pos);

	  }
	  ImGui::End();
	

  };
  polyscope::state::userCallback = polyscope_callback;

  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
