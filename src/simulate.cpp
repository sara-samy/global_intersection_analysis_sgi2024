#include <Eigen/Dense>
#include <iostream>

#include "igl/doublearea.h"
#include "igl/edge_lengths.h"
#include "igl/edges.h"
#include "igl/readOBJ.h"
#include "igl/triangle_triangle_adjacency.h"
#include "polyscope/point_cloud.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::Vector3d;
using Eigen::VectorXi;

#include <Eigen/Core>
#include <vector>
#include <cmath>
#include <algorithm>

class Hash {
public:
	Hash(double spacing, int maxNumObjects)
		: spacing(spacing), tableSize(5 * maxNumObjects),
		cellStart(tableSize + 1, 0),
		cellEntries(maxNumObjects, 0),
		queryIds(maxNumObjects, 0),
		querySize(0),
		maxNumObjects(maxNumObjects),
		firstAdjId(maxNumObjects + 1, 0),
		adjIds(10 * maxNumObjects, 0) {}

	int hashCoords(int xi, int yi, int zi) const {
		int h = (xi * 92837111) ^ (yi * 689287499) ^ (zi * 283923481);
		return std::abs(h) % tableSize;
	}

	int intCoord(double coord) const {
		return static_cast<int>(std::floor(coord / spacing));
	}

	int hashPos(const Eigen::MatrixXd& pos, int nr) const {
		return hashCoords(
			intCoord(pos(nr, 0)),
			intCoord(pos(nr, 1)),
			intCoord(pos(nr, 2))
		);
	}

	void create(const Eigen::MatrixXd& pos) {
		int numObjects = std::min(static_cast<int>(pos.rows()), static_cast<int>(cellEntries.size()));

		// determine cell sizes
		std::fill(cellStart.begin(), cellStart.end(), 0);
		std::fill(cellEntries.begin(), cellEntries.end(), 0);

		for (int i = 0; i < numObjects; i++) {
			int h = hashPos(pos, i);
			cellStart[h]++;
		}

		// determine cell starts
		int start = 0;
		for (int i = 0; i < tableSize; i++) {
			start += cellStart[i];
			cellStart[i] = start;
		}
		cellStart[tableSize] = start;  // guard

		// fill in objects ids
		for (int i = 0; i < numObjects; i++) {
			int h = hashPos(pos, i);
			cellStart[h]--;
			cellEntries[cellStart[h]] = i;
		}
	}

	void query(const Eigen::MatrixXd& pos, int nr, double maxDist) {
		int x0 = intCoord(pos(nr, 0) - maxDist);
		int y0 = intCoord(pos(nr, 1) - maxDist);
		int z0 = intCoord(pos(nr, 2) - maxDist);

		int x1 = intCoord(pos(nr, 0) + maxDist);
		int y1 = intCoord(pos(nr, 1) + maxDist);
		int z1 = intCoord(pos(nr, 2) + maxDist);

		querySize = 0;

		for (int xi = x0; xi <= x1; xi++) {
			for (int yi = y0; yi <= y1; yi++) {
				for (int zi = z0; zi <= z1; zi++) {
					int h = hashCoords(xi, yi, zi);
					int start = cellStart[h];
					int end = cellStart[h + 1];

					for (int i = start; i < end; i++) {
						queryIds[querySize] = cellEntries[i];
						querySize++;
					}
				}
			}
		}
	}

	void queryAll(const Eigen::MatrixXd& pos, double maxDist) {
		int num = 0;
		double maxDist2 = maxDist * maxDist;

		for (int i = 0; i < maxNumObjects; i++) {
			int id0 = i;
			firstAdjId[id0] = num;
			query(pos, id0, maxDist);

			for (int j = 0; j < querySize; j++) {
				int id1 = queryIds[j];
				if (id1 >= id0)
					continue;
				double dist2 = (pos.row(id0) - pos.row(id1)).squaredNorm();
				if (dist2 > maxDist2)
					continue;

				if (num >= adjIds.size()) {
					adjIds.resize(2 * num);  // dynamic array resizing
				}
				adjIds[num++] = id1;
			}
		}

		firstAdjId[maxNumObjects] = num;
	}

	const std::vector<int>& getFirstAdjId() const {
		return firstAdjId;
	}

	const std::vector<int>& getAdjIds() const {
		return adjIds;
	}

private:
	double spacing;
	int tableSize;
	std::vector<int> cellStart;
	std::vector<int> cellEntries;
	std::vector<int> queryIds;
	int querySize;
	int maxNumObjects;
	std::vector<int> firstAdjId;
	std::vector<int> adjIds;
};


class RigidObject {
public:
    int numParticles;
    MatrixXd pos;

    VectorXd invMass;

    RigidObject(const MatrixXd &V) {
        numParticles = V.rows();
        pos = V;

        // Set inverse masses of particles to 0
        invMass = VectorXd::Zero(numParticles);
    }
};

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

  int topLeftCorner = -1;
  int topRightCorner = -1;

  // will store gradient of constraint functions later.
  VectorXd grads;

  Hash hash;
  float spacing=0.01;

  void initPhysics(const MatrixXi &F);
  void Simulate(double frameDt, int numSubSteps, Eigen::Vector3d gravity, float groundHeight);
  void solveConstraints(double dt);
  void solveGroundCollisions(float groundHeight);
  void solveCollisions(double dt);


  Cloth(const MatrixXd &V, const MatrixXi &F, float bendingCompliance);
};

Cloth::Cloth(const MatrixXd &V, const MatrixXi &F, float bendingCompliance): hash(spacing,V.rows())
{

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

  invMass = VectorXd::Zero(numParticles);
  

  initPhysics(F);
}
void Cloth::initPhysics(const MatrixXi &F) {

	thickness = 0.01f;
	stretchingCompliance = 0.0;
	bendingCompliance = 1.0;
	handleCollisions = true;

  // Compute the edge lengths

  igl::edge_lengths(pos, F, stretchingLengths);
  std::cout << "stretchingLengths has size " << stretchingLengths.rows()
            << " x " << stretchingLengths.cols() << std::endl;

  // Compute inverse masses of particles
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

  // Finding the top corners
  double maxY = initialVertices.col(1).maxCoeff();

  for (int i = 0; i < initialVertices.rows(); ++i) {

	  if (initialVertices(i, 1) == maxY) {
		  if (topLeftCorner == -1 || initialVertices(i, 0) < initialVertices(topLeftCorner, 0)) {
			  topLeftCorner = i;
			  invMass(i) = 0;

		  }
		  if (topRightCorner == -1 || initialVertices(i, 0) > initialVertices(topRightCorner, 0)) {
			  topRightCorner = i;
			  invMass(i) = 0;

		  }
	  }
  }


}

void Cloth::Simulate(double frameDt, int numSubSteps, Eigen::Vector3d gravity, float groundHeight=0)
{

	double dt = frameDt / numSubSteps;
	double maxVelocity = 0.2 * thickness / dt;

	
	if (handleCollisions) {
		hash.create(pos);
		double maxTravelDist = maxVelocity * frameDt;
		hash.queryAll(pos, maxTravelDist);
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

		solveGroundCollisions(groundHeight);

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
void Cloth::solveGroundCollisions(float groundHeight=0) {
	for (int i = 0; i < numParticles; i++) 
	{
		if (invMass(i) == 0.0)
			continue;

		// Check the y-coordinate of the particle
		double y = pos(i, 1);
		//float ground = -23;
		float ground = groundHeight;
		if (y < ground) {
			// Apply damping
			double damping = 1.0;
			Eigen::Vector3d displacement = pos.row(i) - prevPos.row(i);
			pos.row(i) += -damping * displacement;

			// Ensure the particle stays above the ground plane
			pos(i, 1) = ground + thickness;
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



bool run = false;

int main() {


  // Read the mesh of cloth
  MatrixXd meshV; // #V x 3
  MatrixXi meshF; // #F x 3

  std::string baseDir = "../data/";
  std::string meshfilename = "cloth.obj";
  std::string meshPath = baseDir + meshfilename;
  igl::readOBJ(meshPath, meshV, meshF);

  // Rotate the plane mesh
  //for (int i = 0; i < meshV.rows(); ++i) {
  //    double temp = meshV(i, 0);
  //    meshV(i, 0) = meshV(i, 1);
  //    meshV(i, 1) = temp;
  //}

  Cloth cloth(meshV, meshF, 1.0f);


  // Read the mesh for the rigid plane
  //MatrixXd rigidMeshV;
  //MatrixXi rigidMeshF;
  //std::string rigidMeshfilename = "rigid_plane.obj";
  //std::string rigidMeshPath = baseDir + rigidMeshfilename;
  //igl::readOBJ(rigidMeshPath, rigidMeshV, rigidMeshF);

  // Rotate the rigid mesh
  //for (int i = 0; i < rigidMeshV.rows(); ++i) {
  //    double temp = rigidMeshV(i, 2);
  //    rigidMeshV(i, 2) = rigidMeshV(i, 1);
  //    rigidMeshV(i, 1) = temp;
  //}
  // Translate the rigid plane below cloth plane.
  //double planeLevel = 0.0;
  //double minY = std::numeric_limits<double>::min();
  //for (int i = 0; i < meshV.rows(); ++i) {
  //    double y = meshV(i, 1);

  //    if (y < minY) {
  //        minY = y;
  //        planeLevel = meshV(i, 1);
  //    }
  //}
  //for (int i = 0; i < rigidMeshV.rows(); ++i) {
  //  rigidMeshV(i, 1) = planeLevel;
  //}

  // Initialize polyscope with some options
  polyscope::view::setUpDir(polyscope::UpDir::ZUp);
  polyscope::view::setFrontDir(polyscope::FrontDir::XFront);
  polyscope::init();

  // Register the mesh with Polyscope
  polyscope::registerSurfaceMesh("Cloth", meshV, meshF);
  //polyscope::registerSurfaceMesh("Fake ground", rigidMeshV, rigidMeshF);

  


  Vector3d gravity(0,-9.8,0); 
  double dt = 0.01;
  int subSteps = 5;
  float planeHeight = 0;
  float collisionPlaneOffset = 0;
  auto polyscope_callback = [&]() mutable
  {
      ImGui::PushItemWidth(100);

	  ImGui::Begin("Simulator");
	  if (ImGui::Button(run ? "Stop simulation" : "Run simulation")) {
		  run = !run;
	  }
  	
	  ImGui::InputDouble("dt", &dt);
	  ImGui::SliderInt("Substeps", &subSteps, 1, 20);
	  ImGui::SliderFloat("Bending compliance", &cloth.bendingCompliance, 0, 1);
	  ImGui::SliderFloat("Stretching compliance", &cloth.stretchingCompliance, 0, 1);
	  ImGui::SliderFloat("Cloth Thickness", &cloth.thickness, 0, 20, "%.2f");
	  ImGui::Checkbox("Handle Collisions", &cloth.handleCollisions);
	  ImGui::SliderFloat("Ground Height", &planeHeight, -20, 20);
	  ImGui::SliderFloat("Collision Plane Offset", &collisionPlaneOffset, -25, 25);
	//ImGui::Checkbox("G")
	  polyscope::options::groundPlaneHeightFactor = polyscope::absoluteValue(planeHeight);
	  if (ImGui::Button("Reset Simulation"))
	  {
		  planeHeight = 0;
		  cloth.pos = cloth.initialVertices;
		  polyscope::getSurfaceMesh("Cloth")->updateVertexPositions(cloth.pos);

	  }
	if(ImGui::Button("Fix point"))
	{
		
	}

	  if (run)
	  {
	  	cloth.Simulate(dt, subSteps, gravity,-planeHeight-collisionPlaneOffset);
		polyscope::getSurfaceMesh("Cloth")->updateVertexPositions(cloth.pos);
	  }

	  ImGui::End();
  };
  polyscope::view::setUpDir(polyscope::UpDir::YUp);
  polyscope::state::userCallback = polyscope_callback;


  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
