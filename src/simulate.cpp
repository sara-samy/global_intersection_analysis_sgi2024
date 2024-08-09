#include <Eigen/Dense>
#include <iostream>

#include "igl/doublearea.h"
#include "igl/edge_lengths.h"
#include "igl/edges.h"
#include "igl/readOBJ.h"
#include "igl/triangle_triangle_adjacency.h"
#include "igl/predicates/predicates.h"
#include "ipc/broad_phase/bvh.hpp"
#include "ipc/distance/edge_edge.hpp"
#include "ipc/distance/point_edge.hpp"
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


///////////////////

bool custom_is_edge_intersecting_triangle(const Eigen::Vector3d& e0,
    const Eigen::Vector3d& e1,
    const Eigen::Vector3d& t0,
    const Eigen::Vector3d& t1,
    const Eigen::Vector3d& t2, double& u,
    double& v, double& t

) {
    igl::predicates::exactinit();
    const auto ori1 = igl::predicates::orient3d(t0, t1, t2, e0);
    const auto ori2 = igl::predicates::orient3d(t0, t1, t2, e1);

    if (ori1 != igl::predicates::Orientation::COPLANAR &&
        ori2 != igl::predicates::Orientation::COPLANAR && ori1 == ori2) {
        // edge is completly on one side of the plane that triangle is in
        return false;
    }

    Eigen::Matrix3d M;
    M.col(0) = t1 - t0;
    M.col(1) = t2 - t0;
    M.col(2) = e0 - e1;
    Eigen::Vector3d uvt = M.fullPivLu().solve(e0 - t0);

    if (uvt[0] >= 0.0 && uvt[1] >= 0.0 && uvt[0] + uvt[1] <= 1.0 &&
        uvt[2] >= 0.0 && uvt[2] <= 1.0) {
        u = uvt[0];
        v = uvt[1];
        t = uvt[2];
        return true;
    }
    else
        return false;

    return uvt[0] >= 0.0 && uvt[1] >= 0.0 && uvt[0] + uvt[1] <= 1.0 &&
        uvt[2] >= 0.0 && uvt[2] <= 1.0;
}


void floodVertex(int vertexID, Eigen::MatrixXd meshV, Eigen::MatrixXi meshE, const Eigen::VectorXi& contourIntersectingEdgeIds,
    Eigen::VectorXi& color, int color1, int color2)
{

    std::vector<int> otherVertices;
    std::vector<int> otherColor;
    int otherVertexNumber = 0;

    //color current vertex
    //color(vertexID) = startingColor;

    for (int i = 0; i < meshE.rows(); i++)
    {
        bool connecting = false;
        int otherVertex;
        //for all vertices connected to v
        if (meshE(i, 0) == vertexID)
        {
            connecting = true;
            otherVertex = meshE(i, 1);

        }
        else if (meshE(i, 1) == vertexID)
        {
            connecting = true;
            otherVertex = meshE(i, 0);
        }

        if (connecting && color(otherVertex, 0) == 0) // if they are connected and the next vertex has no assigned color, add them to list
        {
            otherVertices.push_back(otherVertex);

            bool intersecting = false;
            for (int j = 0; j < contourIntersectingEdgeIds.rows(); j++)
            {
                if (contourIntersectingEdgeIds(j) == i)
                {
                    intersecting = true;
                    break;
                }
            }
            if (intersecting)
            {
                //swap the starting color if it intersects a contour during flood fill
                if (color(vertexID) == color2)otherColor.push_back(color1);
                else if (color(vertexID) == color1)otherColor.push_back(color2);
            }
            else
                otherColor.push_back(color(vertexID));
            otherVertexNumber++;

        }

    }

    if (otherVertexNumber > 0) // color all vertices in list if list size greater than 0
    {
        //color vertices
        for (int i = 0; i < otherVertexNumber; i++)
        {
            color(otherVertices[i]) = otherColor[i];
        }

        //flood fill vertices
        for (int i = 0; i < otherVertexNumber; i++)
        {
            floodVertex(otherVertices[i], meshV, meshE, contourIntersectingEdgeIds, color, color1, color2);
        }
    }
}


void colorVertices(Eigen::MatrixXd meshV, Eigen::MatrixXi meshE, Eigen::MatrixXi meshF, std::string polyscopeName, Eigen::MatrixXi contourEdges, Eigen::MatrixXd intersectingVertices, Eigen::VectorXi& color, int color1, int color2)
{
    std::vector<ipc::EdgeEdgeCandidate> edge_edge_candidates;


    Eigen::MatrixXd allVertices(meshV.rows() + intersectingVertices.rows(), 3);
    allVertices << meshV, intersectingVertices;

    Eigen::MatrixXi allEdges(meshE.rows() + contourEdges.rows(), 2);

    Eigen::MatrixXi indexshiftedContourEdges(contourEdges.rows(), 2);
    indexshiftedContourEdges << contourEdges;
    indexshiftedContourEdges = indexshiftedContourEdges.rowwise() + Eigen::Vector2i(meshV.rows(), meshV.rows()).transpose();

    allVertices << meshV, intersectingVertices;
    allEdges << meshE, indexshiftedContourEdges;


    ipc::BVH bvh;
    bvh.build(allVertices, allEdges, meshF, 1e-4);

    bvh.can_vertices_collide = [meshV](size_t vi, size_t vj) {
        return (vi <= meshV.rows() && vj > meshV.rows()) || (vi > meshV.rows() && vj <= meshV.rows());
    };



    bvh.detect_edge_edge_candidates(edge_edge_candidates);


    color = color.setZero();
    int firstFloodVertex = 0;

    bool chosen = false;

    for (auto candidate : edge_edge_candidates)
    {


        Eigen::Vector3d a = allVertices.row(allEdges(candidate.edge0_id, 0));
        Eigen::Vector3d b = allVertices.row(allEdges(candidate.edge0_id, 1));
        Eigen::Vector3d c = allVertices.row(allEdges(candidate.edge1_id, 0));
        Eigen::Vector3d d = allVertices.row(allEdges(candidate.edge1_id, 1));

        // Check for degenerate edges
        if ((a == b) || (c == d)) {
            // One of the edges is degenerate, skip this candidate
            edge_edge_candidates.erase(std::remove(edge_edge_candidates.begin(), edge_edge_candidates.end(), candidate), edge_edge_candidates.end());
            continue;
        }
        //check if an edge has a vertex near the countour edge, remove the edge if this is the case
        if (ipc::point_edge_distance(a, c, d) < 1e-4)
            if (allEdges(candidate.edge0_id, 0) < meshV.rows())
            {
                color(allEdges(candidate.edge0_id, 0)) = color2;
                if (!chosen)
                {
                    chosen = true;
                    firstFloodVertex = allEdges(candidate.edge0_id, 0);

                }
            }
        if (ipc::point_edge_distance(b, c, d) < 1e-4)
            if (allEdges(candidate.edge0_id, 1) < meshV.rows())
            {
                color(allEdges(candidate.edge0_id, 1)) = color2;
                if (!chosen)
                {
                    chosen = true;
                    firstFloodVertex = allEdges(candidate.edge0_id, 1);

                }
            }
        if (ipc::point_edge_distance(c, a, b) < 1e-4)
            if (allEdges(candidate.edge1_id, 0) < meshV.rows())
            {
                color(allEdges(candidate.edge1_id, 0)) = color2;
                if (!chosen)
                {
                    chosen = true;
                    firstFloodVertex = allEdges(candidate.edge1_id, 0);

                }
            }
        if (ipc::point_edge_distance(d, a, b) < 1e-4)
            if (allEdges(candidate.edge1_id, 1) < meshV.rows())
            {

                color(allEdges(candidate.edge1_id, 1)) = color2;
                if (!chosen)
                {
                    chosen = true;
                    firstFloodVertex = allEdges(candidate.edge1_id, 1);

                }
            }

    }

    Eigen::MatrixXd contourIntersectingEdges(0, 2);
    Eigen::VectorXi contourIntersectingEdgesIDs(0);
    //take each edge, compare to each contour edge
    for (auto candidate : edge_edge_candidates)
    {

        Eigen::Vector3d a = allVertices.row(allEdges(candidate.edge0_id, 0));
        Eigen::Vector3d b = allVertices.row(allEdges(candidate.edge0_id, 1));
        Eigen::Vector3d c = allVertices.row(allEdges(candidate.edge1_id, 0));
        Eigen::Vector3d d = allVertices.row(allEdges(candidate.edge1_id, 1));
        if (ipc::edge_edge_distance(a, b, c, d) < 1e-5)
        {
            if (candidate.edge0_id < meshE.rows())//if first edge belongs to first mesh 
            {
                int k = contourIntersectingEdges.rows();
                contourIntersectingEdges.conservativeResize(k + 1, 2);
                contourIntersectingEdges(k, 0) = allEdges(candidate.edge0_id, 0);
                contourIntersectingEdges(k, 1) = allEdges(candidate.edge0_id, 1);

                contourIntersectingEdgesIDs.conservativeResize(k + 1);
                contourIntersectingEdgesIDs(k, 0) = candidate.edge0_id;

            }
            else if (candidate.edge1_id < meshE.rows())//if second edge belongs to first mesh 
            {
                int k = contourIntersectingEdges.rows();
                contourIntersectingEdges.conservativeResize(k + 1, 2);
                contourIntersectingEdges(k, 0) = allEdges(candidate.edge1_id, 0);
                contourIntersectingEdges(k, 1) = allEdges(candidate.edge1_id, 1);

                contourIntersectingEdgesIDs.conservativeResize(k + 1);
                contourIntersectingEdgesIDs(k, 0) = candidate.edge1_id;
            }
        }
    }


    if (firstFloodVertex == 0) color(0) = color1;
    floodVertex(firstFloodVertex, meshV, meshE, contourIntersectingEdgesIDs, color, color1, color2);
    polyscope::getSurfaceMesh(polyscopeName)->addVertexScalarQuantity("coloring", color);



}


void createContour(Eigen::MatrixXd meshV, Eigen::MatrixXi meshF, Eigen::MatrixXd meshVc, Eigen::MatrixXi meshFc, Eigen::VectorXi& color1, Eigen::VectorXi& color2)
{
    Eigen::MatrixXd intersectingVertices;

    Eigen::VectorXi intersection(meshF.rows(), 1);


    // find intersecting triangles
    for (int i = 0; i < meshF.rows(); i++) // faces of surface mesh
    {
        for (int j = 0; j < meshFc.rows(); j++) // faces of cloth
        {
            // vertex indices
            int vIndex1 = meshFc(j, 0);
            int vIndex2 = meshFc(j, 1);
            int vIndex3 = meshFc(j, 2);

            // cloth vertices
            Eigen::Vector3d v1 = meshVc.row(vIndex1);
            Eigen::Vector3d v2 = meshVc.row(vIndex2);
            Eigen::Vector3d v3 = meshVc.row(vIndex3);

            // surface vertices
            Eigen::Vector3d sv1 = meshV.row(meshF(i, 0));
            Eigen::Vector3d sv2 = meshV.row(meshF(i, 1));
            Eigen::Vector3d sv3 = meshV.row(meshF(i, 2));

            // go through each edge of cloth and check if intersecting with a face of
            // surface mesh

            double u = 0, v = 0, t = 0;

            bool i1; // = ipc::is_edge_intersecting_triangle(v1, v2, sv1, sv2, sv3);
            bool i2; // = ipc::is_edge_intersecting_triangle(v2, v3, sv1, sv2, sv3);
            bool i3; // = ipc::is_edge_intersecting_triangle(v3, v1, sv1, sv2, sv3);

            i1 = custom_is_edge_intersecting_triangle(v1, v2, sv1, sv2, sv3, u, v, t);
            i2 = custom_is_edge_intersecting_triangle(v2, v3, sv1, sv2, sv3, u, v, t);
            i3 = custom_is_edge_intersecting_triangle(v3, v1, sv1, sv2, sv3, u, v, t);

            // does edge orientation matter?

            if (i1 || i2 || i3) {
                // std::cout << "Intersection";
                intersection(i, 0) = 1;

                // find points of intersection
                Eigen::Vector3d iVector = u * (sv2 - sv1) + v * (sv3 - sv1) + sv1;

                intersectingVertices.conservativeResize(intersectingVertices.rows() + 1,
                    3);
                intersectingVertices.row(intersectingVertices.rows() - 1) = iVector;
            }
        }
    }

    // create path
    // using intersecting vertices- find closest point- mark it using a boolean
    Eigen::VectorXi inPath(intersectingVertices.rows(), 1);
    inPath = inPath.setZero();

    if (intersectingVertices.rows() > 0)
    {
        Eigen::Vector3d currentVertex = intersectingVertices.row(0);
        Eigen::MatrixXi contourEdges(intersectingVertices.rows(), 2);

        inPath(0, 0) = 1;
        int current = 0;
        int next = 0;
        double shortest = 9999;
        int k = 0;
        contourEdges(k, 0) = current;
        while (inPath.minCoeff() == 0) {
            // get all vertices in triangle

            for (int i = 0; i < intersectingVertices.rows(); i++) {
                if (i == current || inPath(i, 0) == 1)
                    continue;

                Eigen::Vector3d comparingVertex = intersectingVertices.row(i);
                double currentDistance = glm::distance(
                    glm::vec3(currentVertex.x(), currentVertex.y(), currentVertex.z()),
                    glm::vec3(comparingVertex.x(), comparingVertex.y(),
                        comparingVertex.z()));

                if (currentDistance < shortest) {
                    shortest = currentDistance;
                    next = i;
                }
            }
            inPath(next, 0) = 1;
            contourEdges(k, 0) = current;
            contourEdges(k, 1) = next;
            current = next;
            currentVertex = intersectingVertices.row(current);
            shortest = 9999;
            k++;
        }
        contourEdges(k, 0) = contourEdges(k - 1, 1);
        contourEdges(k, 1) = contourEdges(0, 0);


        Eigen::MatrixXi meshE;
        Eigen::MatrixXi meshEc;


        igl::edges(meshF, meshE);
        igl::edges(meshFc, meshEc);



        int colorOutsideContour = 1;
        int colorInsideContour = 2;
        int colorInsideContour2 = 3;

        colorVertices(meshV, meshE, meshF, "Cloth", contourEdges, intersectingVertices, color1, colorOutsideContour, colorInsideContour);
        colorVertices(meshVc, meshEc, meshFc, "Ball", contourEdges, intersectingVertices, color2, colorOutsideContour, colorInsideContour2);



        //if (contourEdges.rows() > 0)
          //  polyscope::registerCurveNetwork("contour", intersectingVertices, contourEdges);

    }

    //if(intersection.rows()>0)
      // polyscope::getSurfaceMesh("input1")->addFaceScalarQuantity("intersection",intersection);
    //if (intersectingVertices.rows() > 0)
      //  polyscope::registerPointCloud("points", intersectingVertices);

}

Eigen::Vector3d getAttractiveForceDirection(Eigen::MatrixXd meshV, Eigen::MatrixXi meshF, Eigen::MatrixXd meshVc, Eigen::MatrixXi meshFc, Eigen::VectorXi color1, Eigen::VectorXi color2)
{


    createContour(meshV, meshF, meshVc, meshFc, color1, color2);
    //using both positions, find attractive force direction
    Eigen::Vector3d avg1 = { 0,0,0 };
    int colorCount = 0;
    for (int i = 0; i < meshV.rows(); i++)
    {
        if (color1(i) == 2)//we use value 2 for color for mesh 1 for now
        {
            avg1 += meshV.row(i);
            colorCount++;
        }
    }
    if (colorCount == 0)avg1.setZero();
    else avg1 /= colorCount;


    Eigen::Vector3d avg2 = { 0,0,0 };
    colorCount = 0;
    for (int i = 0; i < meshVc.rows(); i++)
    {
        if (color2(i) == 3)//we use value 2 for color for mesh 1 for now
        {
            avg2 += meshVc.row(i);
            colorCount++;
        }
    }
    if (colorCount == 0)avg2.setZero();
    else avg2 /= colorCount;



    Eigen::MatrixXd vn(2, 3);
    vn.row(0) = avg1;
    vn.row(1) = avg2;



    polyscope::registerPointCloud("averages", vn);

    return avg2 - avg1;
}



///////////////////
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

    //for now we copy this from the contouring function
  VectorXi coloredVertices;
  Vector3d attractiveForceDir;//direction of GIA intersection resolving force

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
	double maxVelocity = thickness / dt;

	
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

		//solveGroundCollisions(groundHeight);

		solveConstraints(dt);

        //add in cloth collision resolution
        for (int i = 0; i < numParticles; i++)
        {
            if (invMass(i) > 0.0) 
            {
                Eigen::Vector3d velocity = vel.row(i);
                pos.row(i) += attractiveForceDir * dt * 0.5;
            }
        }


		//if (handleCollisions)
			//solveCollisions(dt);

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

void Cloth::solveCollisions(double dt){
	double thickness2 = thickness * thickness;

	for (int i = 0; i < numParticles; i++) {
		if (invMass(i) == 0.0)
			continue;

		int id0 = i;
		int first = hash.firstAdjId[i];
		int last = hash.firstAdjId[i + 1];

		for (int j = first; j < last; j++) {
			int id1 = hash.adjIds[j];
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
			double friction = 0.5;
			pos.row(id0) += friction * vel0;
			pos.row(id1) += friction * vel1;
		}
	}
}



bool run = false;

int main() {


  // Read the mesh of cloth
  MatrixXd meshV; // #V x 3
  MatrixXd meshVb; // #V x 3
  MatrixXi meshF; // #F x 3
  MatrixXi meshFb; // #F x 3

  std::string baseDir = "../data/";
  std::string meshfilename1 = "plane.obj";
  std::string meshfilename2 = "sphere3.obj";


  std::string meshPath1 = baseDir + meshfilename1;
  igl::readOBJ(meshPath1, meshV, meshF);
  std::string meshPath2 = baseDir + meshfilename2;
  igl::readOBJ(meshPath2, meshVb, meshFb);

  // Rotate the plane mesh
  //for (int i = 0; i < meshV.rows(); ++i) {
  //    double temp = meshV(i, 0);
  //    meshV(i, 0) = meshV(i, 1);
  //    meshV(i, 1) = temp;
  //}

  for (int i = 0; i < meshV.rows(); i++)
  {
	  meshV.row(i) = meshV.row(i) * 0.05;
  }
  meshV.rowwise() +=
	  Eigen::Vector3d(0.0, 0.5, 0.7)
	  .transpose();
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
  polyscope::view::setUpDir(polyscope::UpDir::YUp);
  polyscope::view::setFrontDir(polyscope::FrontDir::XFront);
  polyscope::init();

  // Register the mesh with Polyscope
  polyscope::registerSurfaceMesh("Cloth", meshV, meshF);
  polyscope::registerSurfaceMesh("Ball", meshVb, meshFb);
  //polyscope::registerSurfaceMesh("Fake ground", rigidMeshV, rigidMeshF);



  //////////////////////////////////////////////GIA
  
  Eigen::VectorXi color1(meshV.rows());
  Eigen::VectorXi color2(meshVb.rows());

  Eigen::Vector3d forceDir = getAttractiveForceDirection(meshV, meshF, meshVb, meshFb, color1, color2).normalized();
  cloth.coloredVertices = color1;
  cloth.attractiveForceDir = forceDir;



  /////////////////////////////////////////////SIMULATION
  
  
  Vector3d gravity(0,-9.8,0); 
  double dt = 0.01;
  int subSteps = 5;
  float planeHeight = 0;
  float collisionPlaneOffset = 23;
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
	  if (run)
	  {

		  //ADD STEP TO DETECT COLLISIONS HERE
		  
		  //get attractive forces from both given meshes


          forceDir = getAttractiveForceDirection(meshV, meshF, meshVb, meshFb, color1, color2).normalized();
          cloth.coloredVertices = color1;
          cloth.attractiveForceDir = forceDir;

	  	cloth.Simulate(dt, subSteps, gravity, -planeHeight - collisionPlaneOffset);

		polyscope::getSurfaceMesh("Cloth")->updateVertexPositions(cloth.pos);
	  }

	  ImGui::End();
  };
  polyscope::state::userCallback = polyscope_callback;
  ////////////////////////////////////////////

  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
