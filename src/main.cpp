#include <iostream>

//#include "args/args.hxx"

#include "igl/edges.h"
#include "igl/predicates/predicates.h"
#include "igl/readOBJ.h"
#include "ipc/collision_mesh.hpp"
#include "ipc/collisions/collisions.hpp"
#include "ipc/potentials/barrier_potential.hpp"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"
#include <polyscope/polyscope.h>

#include "ipc/distance/edge_edge.hpp"
#include "ipc/utils/intersection.hpp"

#include "ipc/utils/save_obj.hpp"
#include "polyscope/curve_network.h"

bool custom_is_edge_intersecting_triangle(const Eigen::Vector3d &e0,
                                          const Eigen::Vector3d &e1,
                                          const Eigen::Vector3d &t0,
                                          const Eigen::Vector3d &t1,
                                          const Eigen::Vector3d &t2, double &u,
                                          double &v, double &t

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
  } else
    return false;

  return uvt[0] >= 0.0 && uvt[1] >= 0.0 && uvt[0] + uvt[1] <= 1.0 &&
         uvt[2] >= 0.0 && uvt[2] <= 1.0;
}

void add_edge(std::set<std::pair<int, int>>& edge_set, int v1, int v2)
{

    if (v1 > v2) std::swap(v1, v2);  
	edge_set.insert(std:: make_pair(v1, v2));
}

// Function to find all unique edges from a faces matrix
Eigen::MatrixXi find_unique_edges(const Eigen::MatrixXi& meshF)
{
    std::set<std::pair<int, int>> edge_set;

    for (int i = 0; i < meshF.rows(); ++i) {
        int v0 = meshF(i, 0);
        int v1 = meshF(i, 1);
        int v2 = meshF(i, 2);
    	add_edge(edge_set, v0, v1);
        add_edge(edge_set, v1, v2);
        add_edge(edge_set, v2, v0);
    }
	Eigen::MatrixXi meshE(edge_set.size(), 2);
    int row = 0;
    for (const auto edge : edge_set) {
        meshE(row, 0) = edge.first;
        meshE(row, 1) = edge.second;
        ++row;
	}

    //std::cout << meshE;

    return meshE;
}

int main(int argc, char **argv) {

    // Read the mesh
    Eigen::MatrixXd meshV;  // V x 3
    Eigen::MatrixXd meshVc; // Vx3
    Eigen::MatrixXi meshF;  // F x 3
    Eigen::MatrixXi meshFc;

    std::string baseDir = "../data/";

    std::string meshfilename1 = "sphere3.obj";
    std::string meshfilename2 = "cloth.obj";

    // Concatenate the base directory path and the filename
    std::string mesh1Path = baseDir + meshfilename1;
	std::string mesh2Path = baseDir + meshfilename2;
    // Read the first mesh
    igl::readOBJ(mesh1Path, meshV, meshF);

    // Read the second mesh
    igl::readOBJ(mesh2Path, meshVc, meshFc);

#if __APPLE__
    
  // Configure the argument parser
  args::ArgumentParser parser("Global intersection analysis demo", "");
  // args::Positional<std::string> inFile(parser, "mesh", "input mesh");
  args::PositionalList<std::string> inFiles(parser, "meshes",
                                            "input mesh files");

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

  std::string filename1 = inFiles.Get()[0];
  std::cout << "loading: " << filename1 << std::endl;
  igl::readOBJ(filename1, meshV, meshF);


  std::string filename2 = inFiles.Get()[1];
  std::cout << "loading: " << filename2 << std::endl;
  igl::readOBJ(filename2, meshVc, meshFc);
#endif


  

  polyscope::init();

  meshVc.rowwise() +=
      Eigen::Vector3d(0, 0, 1.3)
          .transpose(); //+ Eigen::Vector3d(-1, -1, 0).transpose();

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

  // Eigen::Vector3d

  // take a point, make it current point, mark it
  //
  // loop
  // for the triangle the point is in
  // find closest any unmarked vertex in triangle, add that to list, mark it
  //  if no closest in triangle, look at joined triangles, and see if those have
  //  points closest which are unmarked, choose that, mark it, add tp list, make
  //  it current point repeat all points are marked
  //

  // while (inPath.maxCoeff()!=2)
  //{

  // }

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

  //get mesh edges
	Eigen::MatrixXi meshE = find_unique_edges(meshF);

  // get the flood fill
  //precompute step
    // find intersecting edges with a contour
    // for now do for a single contour
    //later- multiple contours

    //take each edge, compare to each contour edge
	


    //pick an arbitrary point
    //color it
    //pick all points next to current point
    //if a next point has an edge that intersects with our contour, change argument color to a different one
    //color them with current color value
    //repeat algo














  // Register the mesh with Polyscope
  polyscope::registerSurfaceMesh("input1", meshV, meshF);

  polyscope::registerSurfaceMesh("input2", meshVc, meshFc);

  // polyscope::getSurfaceMesh("input2")->addVertexVectorQuantity("gradient",a);
  // polyscope::getSurfaceMesh("input2")->translate(glm::vec3(3, 0, 0));

  polyscope::getSurfaceMesh("input1")->addFaceScalarQuantity("intersection",
                                                             intersection);

  polyscope::registerPointCloud("points", intersectingVertices);

  polyscope::registerCurveNetwork("contour", meshV, meshE);

  // Show the GUI
  polyscope::show();
  std::cout << "Hello world";
}
