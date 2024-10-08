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

#include "ipc/broad_phase/bvh.hpp"
#include "ipc/distance/edge_edge.hpp"
#include "ipc/distance/point_edge.hpp"
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


void floodVertex(int vertexID, Eigen::MatrixXd meshV, Eigen::MatrixXi meshE, const Eigen::VectorXi& contourIntersectingEdgeIds, 
    Eigen::VectorXi& color, int color1, int color2)
{

    std::vector<int> otherVertices;
    std::vector<int> otherColor;
    int otherVertexNumber = 0;

    //color current vertex
    //color(vertexID) = startingColor;

    for(int i=0;i<meshE.rows(); i++)
    {
        bool connecting = false;
        int otherVertex;
        //for all vertices connected to v
        if(meshE(i,0)==vertexID)
	    {
            connecting = true;
        	otherVertex = meshE(i, 1);

	    }
        else if(meshE(i, 1) == vertexID)
        {
            connecting = true;
            otherVertex = meshE(i, 0);
        }

        if(connecting && color(otherVertex,0)==0) // if they are connected and the next vertex has no assigned color, add them to list
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
    
    if(otherVertexNumber>0) // color all vertices in list if list size greater than 0
    {
	    //color vertices
        for(int i=0;i<otherVertexNumber;i++)
        {
            color(otherVertices[i]) = otherColor[i];
        }

        //flood fill vertices
        for(int i=0;i<otherVertexNumber;i++)
        {
            floodVertex(otherVertices[i], meshV, meshE, contourIntersectingEdgeIds, color,color1,color2);
        }
    }
}


void colorVertices(Eigen::MatrixXd meshV, Eigen::MatrixXi meshE, Eigen::MatrixXi meshF,std::string polyscopeName, Eigen::MatrixXi contourEdges, Eigen::MatrixXd intersectingVertices, Eigen::VectorXi &color, int color1, int color2)
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
    floodVertex(firstFloodVertex, meshV, meshE, contourIntersectingEdgesIDs, color,color1,color2);
    //polyscope::getSurfaceMesh(polyscopeName)->addVertexScalarQuantity("coloring", color);



}


void createContour(Eigen::MatrixXd meshV, Eigen::MatrixXi meshF, Eigen::MatrixXd meshVc, Eigen::MatrixXi meshFc, Eigen::VectorXi &color1, Eigen::VectorXi& color2)
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

        colorVertices(meshV, meshE, meshF, "input1", contourEdges, intersectingVertices, color1, colorOutsideContour, colorInsideContour);
        colorVertices(meshVc, meshEc, meshFc, "input2", contourEdges, intersectingVertices, color2, colorOutsideContour, colorInsideContour2);



        //if (contourEdges.rows() > 0)
          //  polyscope::registerCurveNetwork("contour", intersectingVertices, contourEdges);

    }

    //if(intersection.rows()>0)
       //polyscope::getSurfaceMesh("input1")->addFaceScalarQuantity("intersection",intersection);
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





int main(int argc, char** argv) {

    // Read the mesh
    Eigen::MatrixXd meshV;  // V x 3
    Eigen::MatrixXd meshVc; // Vx3
    Eigen::MatrixXi meshF;  // F x 3
    Eigen::MatrixXi meshFc;

    std::string baseDir = "../data/";

    std::string meshfilename1 = "high_res_plane.obj";
    std::string meshfilename2 = "sphere3.obj";

    // Concatenate the base directory path and the filename
    std::string mesh1Path = baseDir + meshfilename1;
    std::string mesh2Path = baseDir + meshfilename2;
    // Read the first mesh
    igl::readOBJ(mesh1Path, meshV, meshF);

    // Read the second mesh
    igl::readOBJ(mesh2Path, meshVc, meshFc);


    polyscope::init();

    
    
    for (int i = 0; i<meshV.rows(); i++)
    {
        meshV.row(i) = meshV.row(i) * 0.05;
    }

    meshV.rowwise() +=
        Eigen::Vector3d(0.0, 0.5, 5.7)
        .transpose(); //+ Eigen::Vector3d(-1, -1, 0).transpose();


    // Register the mesh with Polyscope
    polyscope::registerSurfaceMesh("input1", meshV, meshF);

    polyscope::registerSurfaceMesh("input2", meshVc, meshFc);

    

    //use this color data for , something
    Eigen::VectorXi color1(meshV.rows());
    Eigen::VectorXi color2(meshVc.rows());

    Eigen::Vector3d forceDir= getAttractiveForceDirection(meshV, meshF, meshVc, meshFc, color1, color2);
    //1 blue
    // 2 black
    // 3 white

    std::cout << forceDir;


  
  // Show the GUI
  //polyscope::state::userCallback = CallbackFunction;

  polyscope::show();
  std::cout << "Hello world";
}
