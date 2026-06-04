//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2025.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation version 3.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// contact: delft3d.support@deltares.nl
// Stichting Deltares
// P.O. Box 177
// 2600 MH Delft, The Netherlands
//
// All indications and logos of, and references to, "Delft3D" and "Deltares"
// are registered trademarks of Stichting Deltares, and remain the property of
// Stichting Deltares. All rights reserved.
//
//------------------------------------------------------------------------------

#include <algorithm>
#include <chrono>
#include <execution>
#include <gtest/gtest.h>
#include <random>

#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/MeshRefinement.hpp>
#include <MeshKernel/Parameters.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/Utilities/Utilities.hpp>
#include <TestUtils/Definitions.hpp>
#include <TestUtils/MakeMeshes.hpp>
#include <TestUtils/PolygonReader.hpp>

TEST(Mesh, OneQuadTestConstructor)
{
    // 1 Setup
    std::vector<meshkernel::Point> nodes;
    nodes.push_back({0.0, 0.0});
    nodes.push_back({0.0, 10.0});
    nodes.push_back({10.0, 0.0});
    nodes.push_back({10.0, 10.0});
    std::vector<meshkernel::Edge> edges;
    edges.push_back({0, 2});
    edges.push_back({1, 3});
    edges.push_back({0, 1});
    edges.push_back({2, 3});

    // 2 Execution
    const auto mesh = meshkernel::Mesh2D(edges, nodes, meshkernel::Projection::cartesian);

    // 3 Validation
    // expect nodesEdges to be sorted ccw
    ASSERT_EQ(0, mesh.m_nodesEdges[0][0]);
    ASSERT_EQ(2, mesh.m_nodesEdges[0][1]);

    ASSERT_EQ(1, mesh.m_nodesEdges[1][0]);
    ASSERT_EQ(2, mesh.m_nodesEdges[1][1]);

    ASSERT_EQ(0, mesh.m_nodesEdges[2][0]);
    ASSERT_EQ(3, mesh.m_nodesEdges[2][1]);

    ASSERT_EQ(1, mesh.m_nodesEdges[3][0]);
    ASSERT_EQ(3, mesh.m_nodesEdges[3][1]);

    // each node has two edges int this case
    ASSERT_EQ(2, mesh.GetNumNodesEdges(0));
    ASSERT_EQ(2, mesh.GetNumNodesEdges(1));
    ASSERT_EQ(2, mesh.GetNumNodesEdges(2));
    ASSERT_EQ(2, mesh.GetNumNodesEdges(3));

    // the nodes composing the face, in ccw order
    ASSERT_EQ(0, mesh.m_facesNodes[0][0]);
    ASSERT_EQ(2, mesh.m_facesNodes[0][1]);
    ASSERT_EQ(3, mesh.m_facesNodes[0][2]);
    ASSERT_EQ(1, mesh.m_facesNodes[0][3]);

    // the edges composing the face, in ccw order
    ASSERT_EQ(0, mesh.m_facesEdges[0][0]);
    ASSERT_EQ(3, mesh.m_facesEdges[0][1]);
    ASSERT_EQ(1, mesh.m_facesEdges[0][2]);
    ASSERT_EQ(2, mesh.m_facesEdges[0][3]);

    // // the found circumcenter for the face
    // ASSERT_DOUBLE_EQ(5.0, mesh.m_facesCircumcenters[0].x);
    // ASSERT_DOUBLE_EQ(5.0, mesh.m_facesCircumcenters[0].y);

    // each edge has only one face in this case
    ASSERT_EQ(1, mesh.GetNumEdgesFaces(0));
    ASSERT_EQ(1, mesh.GetNumEdgesFaces(1));
    ASSERT_EQ(1, mesh.GetNumEdgesFaces(2));
    ASSERT_EQ(1, mesh.GetNumEdgesFaces(3));

    // each edge is a boundary edge, so the second entry of edgesFaces is an invalid index (meshkernel::constants::missing::sizetValue)
    ASSERT_EQ(meshkernel::constants::missing::uintValue, mesh.m_edgesFaces[0][1]);
    ASSERT_EQ(meshkernel::constants::missing::uintValue, mesh.m_edgesFaces[1][1]);
    ASSERT_EQ(meshkernel::constants::missing::uintValue, mesh.m_edgesFaces[2][1]);
    ASSERT_EQ(meshkernel::constants::missing::uintValue, mesh.m_edgesFaces[3][1]);
}

TEST(Mesh2D, TriangulateSamplesWithSkinnyTriangle)
{
    // Prepare
    std::vector<meshkernel::Point> nodes;

    nodes.push_back({302.002502, 472.130371});
    nodes.push_back({144.501526, 253.128174});
    nodes.push_back({368.752930, 112.876755});
    nodes.push_back({707.755005, 358.879242});
    nodes.push_back({301.252502, 471.380371});
    nodes.push_back({302.002502, 472.130371});

    meshkernel::Polygons polygons(nodes, meshkernel::Projection::cartesian);

    // Execute
    const auto generatedPoints = polygons.ComputePointsInPolygons();

    meshkernel::Mesh2D mesh(generatedPoints[0], polygons, meshkernel::Projection::cartesian);

    // Assert
    ASSERT_EQ(5, mesh.GetNumNodes());
    ASSERT_EQ(6, mesh.GetNumEdges());

    ASSERT_EQ(3, mesh.GetEdge(0).first);
    ASSERT_EQ(0, mesh.GetEdge(0).second);

    ASSERT_EQ(0, mesh.GetEdge(1).first);
    ASSERT_EQ(1, mesh.GetEdge(1).second);

    ASSERT_EQ(1, mesh.GetEdge(2).first);
    ASSERT_EQ(3, mesh.GetEdge(2).second);

    ASSERT_EQ(4, mesh.GetEdge(3).first);
    ASSERT_EQ(1, mesh.GetEdge(3).second);

    ASSERT_EQ(1, mesh.GetEdge(4).first);
    ASSERT_EQ(2, mesh.GetEdge(4).second);

    ASSERT_EQ(2, mesh.GetEdge(5).first);
    ASSERT_EQ(4, mesh.GetEdge(5).second);
}

// extern "C"
// {

//     extern void mshoce_(
//         const int* jnew,   // Input: Reset indicator (Fortran Logical pointer)
//         double* coor,      // Output: Flattened array of node coordinates
//         int* kmeshc,       // Output: Connectivity/topology grid matrix
//         const int* inpelm, // Input: Element type identifier
//         const int* nbound, // Input: Number of boundary elements
//         double* bcord,     // Input: Coordinates of boundary control nodes

//         int* kbndpt,            // Input: Type flags for boundary nodes
//         int* boundary,          // Input: Edge-to-node connectivity map
//         const int* numcurvboun, // Input: Total count of curved boundary segments
//         int* npoint,            // Output: Count of generated points
//         int* nelem,             // Output: Count of generated elements

//         int* holeinfo,     // Input: Structural layout parameters for holes
//         const int* nholes, // Input: Total count of internal holes
//         const int* ncoar,  // Input: Quantity of sizing descriptors passed
//         double* coar,      // Input: Target element sizing arrays
//         int* userpoints,   // Input: Fixed target internal points

//         int* isurnr,             // Output: Surface structural adjacency mapping register
//         const int* numextcurves, // Input: Auxiliary alignment curve flags
//         int* numnodextcurvs,     // Input: Node mappings for alignment paths

//         int* curvenumbers, // Input: Curve curvature flags
//         double* rinput,    // Input: Supplementary sizing/weighting matrices
//         const int* nuspnt, // Input: Count of forced target control nodes
//         const int* ndim    // Input: Domain spatial dimension identifier
//     );
// }

// double minimumEdgeDelta(const std::vector<meshkernel::Point>& polygonNodes)
// {
//     double delta = 1.0e20;

//     for (size_t i = 0; i + 1 < polygonNodes.size(); ++i)
//     {
//         double dx = polygonNodes[i + 1].x - polygonNodes[i].x;
//         double dy = polygonNodes[i + 1].y - polygonNodes[i].y;

//         delta = std::min(delta, std::sqrt(dx * dx + dy * dy));
//     }

//     return delta;
// }

// std::vector<std::reference_wrapper<const meshkernel::Polygon>> generatePolygonReferences(const meshkernel::Polygons& polygon)
// {
//     const auto& enclosure = polygon.Enclosure(0);

//     std::vector<std::reference_wrapper<const meshkernel::Polygon>> boundaryLoops;
//     boundaryLoops.reserve(1 + enclosure.NumberOfInner());
//     boundaryLoops.emplace_back(enclosure.Outer());

//     for (meshkernel::UInt i = 0; i < enclosure.NumberOfInner(); ++i)
//     {
//         boundaryLoops.emplace_back(enclosure.Inner(i));
//     }

//     return boundaryLoops;
// }

// std::vector<meshkernel::Point> pointsFromFlatArray(const std::vector<double>& coordinates, const int numberOfPoints)
// {
//     std::vector<meshkernel::Point> meshNodes(numberOfPoints);

//     for (int i = 0; i < numberOfPoints; ++i)
//     {
//         meshNodes[i].x = coordinates[2 * i];
//         meshNodes[i].y = coordinates[2 * i + 1];
//     }

//     return meshNodes;
// }

// std::tuple<std::vector<meshkernel::Edge>, std::vector<std::vector<meshkernel::UInt>>, std::vector<std::uint8_t>> gatherEdgesAndFaces(const std::vector<int>& kmeshc, const int numberOfElements)
// {

//     std::vector<meshkernel::Edge> edges(3 * numberOfElements);
//     std::vector<std::vector<meshkernel::UInt>> faceNodes(numberOfElements);
//     std::vector<std::uint8_t> numFaceNodes(numberOfElements, 0);

//     // Gather all edges in the mesh
//     // The first stage will find shared edges twice, once for each triangle. The duplicate will be removed later.
//     for (int i = 0; i < numberOfElements; ++i)
//     {
//         int index = 3 * i;

//         meshkernel::UInt n1 = static_cast<meshkernel::UInt>(kmeshc[index] - 1);
//         meshkernel::UInt n2 = static_cast<meshkernel::UInt>(kmeshc[index + 1] - 1);
//         meshkernel::UInt n3 = static_cast<meshkernel::UInt>(kmeshc[index + 2] - 1);

//         // Which is better, reserve and push back or allocate and assign?
//         edges[index] = n1 < n2 ? meshkernel::Edge(n1, n2) : meshkernel::Edge(n2, n1);
//         edges[index + 1] = n2 < n3 ? meshkernel::Edge(n2, n3) : meshkernel::Edge(n3, n2);
//         edges[index + 2] = n3 < n1 ? meshkernel::Edge(n3, n1) : meshkernel::Edge(n1, n3);

//         // Set the nodes of the element.
//         faceNodes[i].resize(3);
//         faceNodes[i][0] = n1;
//         faceNodes[i][1] = n2;
//         faceNodes[i][2] = n3;

//         // Indicate that there are 3 nodes for this element
//         numFaceNodes[i] = 3;
//     }

//     // Remove the duplicated edges.
//     // std::pair has a predefined less-than (Edge is a std;:pair<UInt>)
//     std::sort(std::execution::par, edges.begin(), edges.end());
//     auto [first, last] = std::ranges::unique(edges);
//     edges.erase(first, last);

//     return {edges, faceNodes, numFaceNodes};
// }

// meshkernel::Mesh2D generateMesh(const meshkernel::Polygons& poly)
// {
//     std::vector<std::reference_wrapper<const meshkernel::Polygon>> boundaryLoops(generatePolygonReferences(poly));

//     const int numberOfBoundaryNodes = std::accumulate(boundaryLoops.begin(), boundaryLoops.end(), 0, [](int sum, const auto& poly)
//                                                       { return sum + static_cast<int>(poly.get().Size()); });

//     const int dimension = 2;
//     const int elementIdentifier = 3; // 3-node linear triangles
//     const int numberOfPolygons = static_cast<int>(boundaryLoops.size());
//     const int numberOfHoles = numberOfPolygons - 1;

//     // 1. Interleave coordinates into bcord: [x1, y1, x2, y2...]
//     // Flatten the outer polygon followed by each inner polygon, without separators.
//     std::vector<double> boundaryCoordinates(2 * numberOfBoundaryNodes);

//     // 2. Map edgeNodeConnectivity: Fortran uses 1-based indexing
//     std::vector<int> edgeNodeConnectivity(numberOfBoundaryNodes, 0);

//     // 3. Map boundary segments: Fortran boundary(2, numberOfPolygons)
//     // One curve entry per closed loop, flattened in column-major order.

//     std::vector<int> boundaryConnectivity(2 * numberOfPolygons);

//     // 4. Compute elementSizing array for local polyline point matching
//     const int numberElementSizing = 0;         // numberOfBoundaryNodes;
//     std::vector<double> elementSizing(1, 0.0); // 3 * numberOfBoundaryNodes);

//     int pointOffset = 0;

//     for (int loopIndex = 0; loopIndex < numberOfPolygons; ++loopIndex)
//     {
//         const auto& loop = boundaryLoops[loopIndex].get().Nodes();

//         boundaryConnectivity[2 * loopIndex] = loopIndex + 1;
//         boundaryConnectivity[2 * loopIndex + 1] = pointOffset + 1;

//         for (int i = 0; i < static_cast<int>(loop.size()); ++i)
//         {
//             boundaryCoordinates[2 * (pointOffset + i)] = loop[i].x;
//             boundaryCoordinates[2 * (pointOffset + i) + 1] = loop[i].y;
//             edgeNodeConnectivity[pointOffset + i] = pointOffset + i + 1;
//         }

//         pointOffset += static_cast<int>(loop.size());
//     }

//     // 5. Dynamic Memory Allocation for Output Buffers
//     // Get estimate of area covered by polygon, subtracting area covered by holes
//     const double estimatedArea = poly.Enclosure(0).ComputeSurfaceArea();
//     // Compute the minimum spacing between points of the polygon, both outer and all inner polygons
//     const auto [minimumDelta, _] = poly.Enclosure(0).SegmentLengthExtrema();

//     const int estimatedNumberOfElements = static_cast<int>((estimatedArea / (0.433 * minimumDelta * minimumDelta)) * 3.5);
//     const int maximumNumberOfNodes = numberOfBoundaryNodes + 3 * estimatedNumberOfElements;
//     const int maximumNumberOfElements = 2 * maximumNumberOfNodes;

//     std::vector<double> triangulationNodes(dimension * maximumNumberOfNodes, 0.0);
//     std::vector<int> triangulationElementNodes(elementIdentifier * maximumNumberOfElements, 0);

//     // 6. Dummies and Placeholders
//     std::vector<int> holeinfo(2 * (numberOfHoles + 2), 0);
//     int forcedControlPoints = 0;
//     int surfaceSequenceNumber = 1;
//     int auxiliaryAlignment = 0;

//     std::vector<int> numnodextcurvs(1, 0);
//     std::vector<int> curvenumbers(1, 0);
//     std::vector<double> rinput(1, 0.0);
//     std::vector<int> userpoints(1, 0);

//     // 7. Make the Call
//     int newMesh = 1;
//     int numberOfPoints = maximumNumberOfNodes;
//     int numberOfElements = maximumNumberOfElements;

//     mshoce_(&newMesh, triangulationNodes.data(), triangulationElementNodes.data(), &elementIdentifier, &numberOfBoundaryNodes, boundaryCoordinates.data(),
//             edgeNodeConnectivity.data(), boundaryConnectivity.data(), &numberOfPolygons, &numberOfPoints, &numberOfElements,
//             holeinfo.data(), &numberOfHoles, &numberElementSizing, elementSizing.data(), userpoints.data(),
//             &surfaceSequenceNumber, &auxiliaryAlignment, numnodextcurvs.data(), curvenumbers.data(),
//             rinput.data(), &forcedControlPoints, &dimension);

//     // Recover array of Points
//     std::vector<meshkernel::Point> meshNodes(pointsFromFlatArray(triangulationNodes, numberOfPoints));

//     // Recover arrays of edges, face-node connectivity and number of nodes per face.
//     auto [edges, faceNodes, numFaceNodes] = gatherEdgesAndFaces(triangulationElementNodes, numberOfElements);

//     return meshkernel::Mesh2D(edges, meshNodes, faceNodes, numFaceNodes, poly.GetProjection());
// }

#if 0

TEST(Mesh, TriangulateSamples)
{
    // Prepare

    // std::vector<meshkernel::Point> nodes{{0.0, 0.0}, {100.0, 0.0}, {100.0, 100.0}, {0.0, 100.0}, {0.0, 0.0} };
    // // std::vector<meshkernel::Point> nodes{{0.0, 0.0}, {10.0, 0.0}, {10.0, 10.0}, {0.0, 10.0}, {0.0, 0.0} };

    // // nodes.push_back({498.503152894023, 1645.82297461613});
    // // nodes.push_back({-5.90937355559299, 814.854361678898});
    // // nodes.push_back({851.30035347439, 150.079471329115});
    // // nodes.push_back({1411.11078745316, 1182.22995897746});
    // // nodes.push_back({501.418832237663, 1642.90729527249});
    // // nodes.push_back({498.503152894023, 1645.82297461613});

    // meshkernel::Polygons polygons(nodes, meshkernel::Projection::cartesian);

    // nodes = polygons.RefineFirstPolygon (0, 4, 25);
    // // nodes = polygons.RefineFirstPolygon (0, 4, 2.5);
    // // nodes = polygons.RefineFirstPolygon (0, 1, 1.0);
    // meshkernel::Polygons polygons1(nodes, meshkernel::Projection::cartesian);

    // for (size_t i = 0; i < nodes.size (); ++i)
    // {
    //     std::cout << " nodes " << i << "  " << nodes [i].x << ", " << nodes[i].y << std::endl;
    // }

    // nodes = polygons1.RefineFirstPolygon (4, 8, 5.0);
    // // nodes = polygons1.RefineFirstPolygon (4, 8, 0.5);
    // // nodes = polygons1.RefineFirstPolygon (10, 11, 2.0);
    // meshkernel::Polygons polygons2(nodes, meshkernel::Projection::cartesian);

    // // Execute
    // const auto generatedPoints = polygons2.ComputePointsInPolygons();

    // for (size_t i = 0; i < generatedPoints[0].size (); ++i)
    // {
    //     std::cout << " pont " << i << "  " << generatedPoints[0] [i].x << ", " << generatedPoints[0][i].y << std::endl;
    // }

    // auto polys = ReadPolygons ("/home/wcs1/MeshKernel/MeshKernel01/build_deb/northbank_001b.pol", meshkernel::Projection::cartesian);

    // std::vector<meshkernel::Point> nodes{{0.0, 0.0}, {0.5, 0.0}, {1.0, 0.0}, {1.5, 0.0}, {2.0, 0}, {2.5, 0.0}, {3.0, 0.0}, {3.5, 0.0}, {4.0, 0.0}, {4.5, 0.0}, {5.0, 0}, {5.5, 0.0}, {6.0, 0.0}, {6.5, 0.0}, {7.0, 0.0}, {7.5, 0.0}, {8.0, 0.0}, {8.5, 0.0}, {9.0, 0.0}, {9.5, 0.0}, {10.0, 0.0}, {10.0, 2.5}, {10.0, 5.0}, {10.0, 7.5}, {10.0, 10.0}, {7.5, 10.0}, {5.0, 10.0}, {2.5, 10.0}, {0.0, 10.0}, {0.0, 7.5}, {0.0, 5.0}, {0.0, 2.5}, {0.0, 0.0}};

    std::vector<meshkernel::Point> nodes{{0.0, 0.0}, {0.5, 0.0}, {1.0, 0.0}, {1.5, 0.0}, {2.0, 0}, {2.5, 0.0}, {3.0, 0.0}, {3.5, 0.0}, {4.0, 0.0}, {4.5, 0.0}, {5.0, 0}, {5.5, 0.0}, {6.0, 0.0}, {6.5, 0.0}, {7.0, 0.0}, {7.5, 0.0}, {8.0, 0.0}, {8.5, 0.0}, {9.0, 0.0}, {9.5, 0.0}, {10.0, 0.0}, {10.0, 2.5}, {10.0, 5.0}, {10.0, 7.5}, {10.0, 10.0}, {7.5, 10.0}, {5.0, 10.0}, {2.5, 10.0}, {0.0, 10.0}, {0.0, 7.5}, {0.0, 5.0}, {0.0, 2.5}, {0.0, 0.0}, {-998.0, -998.0}, {2.0, 2.0}, {2.5, 2.0}, {3.0, 2.0}, {3.5, 2.0}, {4.0, 2.0}, {4.5, 2.0}, {5.0, 2.0}, {5.0, 3.5}, {5.0, 5.0}, {2.0, 5.0}, {2.0, 2.0}, {-998.0, -998.0}, {6.0, 6.0}, {8.0, 6.0}, {8.0, 8.0}, {6.0, 8.0}, {6.0, 6.0}};

    // std::vector<meshkernel::Point> nodes{{0.0, 0.0}, {1.25, 0.0}, {2.5, 0}, {5.0, 0.0}, {6.25, 0.0}, {7.5, 0}, {8.75, 0.0}, {10.0, 0.0}, {10.0, 2.5}, {10.0, 5.0}, {10.0, 7.5}, {10.0, 10.0}, {7.5, 10.0}, {5.0, 10.0}, {2.5, 10.0}, {0.0, 10.0}, {0.0, 7.5}, {0.0, 5.0}, {0.0, 2.5}, {0.0, 0.0}};

    // std::vector<meshkernel::Point> nodes{{0.0, 0.0}, {2.5, 0}, {5.0, 0.0}, {7.5, 0}, {10.0, 0.0}, {10.0, 2.5}, {10.0, 5.0}, {10.0, 7.5}, {10.0, 10.0}, {7.5, 10.0}, {5.0, 10.0}, {2.5, 10.0}, {0.0, 10.0}, {0.0, 7.5}, {0.0, 5.0}, {0.0, 2.5}, {0.0, 0.0}, {-998.0, -998.0}, {2.0, 2.0}, {5.0, 2.0}, {5.0, 5.0}, {2.0, 5.0}, {2.0, 2.0}, {-998.0, -998.0}, {6.0, 6.0}, {8.0, 6.0}, {8.0, 8.0}, {6.0, 8.0}, {6.0, 6.0}};

    // // std::vector<meshkernel::Point> nodes{{0.0, 0.0}, {5.0, 0.0}, {10.0, 0.0},
    // //                                      {10.0, 5.0}, {10.0, 10.0}, {5.0, 10.0},
    // //                                      {0.0, 10.0}, {0.0, 5.0}, {0.0, 0.0}};

    meshkernel::Polygons polygons(nodes, meshkernel::Projection::cartesian);
    auto polys = &polygons;
    // auto mesh2 = generateMesh (polygons);
    auto mesh2 = generateMesh(polygons);

    meshkernel::SaveVtk(mesh2.Nodes(), mesh2.m_facesNodes, "trianglemesh.vtu");

    // auto polys = ReadPolygons ("/home/wcs1/MeshKernel/MeshKernel01/build_deb/northbank_001b.pol", meshkernel::Projection::cartesian);
    // generateMesh (*polys);
    return;

    const auto generatedPoints = polys->ComputePointsInPolygons();

    auto pnts = polys->GatherAllEnclosureNodes();

    for (size_t i = 0; i < pnts.size(); ++i)
    {
        for (size_t j = i + 1; j < pnts.size(); ++j)
        {

            if (pnts[i] == pnts[j])
            {
                std::cout << " equal points " << i << "  " << j << "  " << pnts[i].x << ", " << pnts[i].y << std::endl;
            }
        }
    }

    meshkernel::Mesh2D mesh(generatedPoints[0], *polys, meshkernel::Projection::cartesian);

    // std::cout << std::endl;
    // std::cout << std::endl;

    // for (size_t i = 0; i < mesh.Nodes ().size (); ++i)
    // {
    //     std::cout << " pont " << i << "  " << mesh.Nodes ()[i].x << ", " << mesh.Nodes ()[i].y << std::endl;
    // }

    meshkernel::SaveVtk(mesh.Nodes(), mesh.m_facesNodes, "trianglemesh.vtu");
}

#endif

TEST(Mesh, TwoTrianglesDuplicatedEdges)
{
    // 1 Setup
    std::vector<meshkernel::Point> nodes;
    nodes.push_back({0.0, 0.0});
    nodes.push_back({5.0, -5.0});
    nodes.push_back({10.0, 0.0});
    nodes.push_back({5.0, 5.0});
    std::vector<meshkernel::Edge> edges;
    edges.push_back({0, 3});
    edges.push_back({0, 2});
    edges.push_back({2, 3});
    edges.push_back({0, 1});
    edges.push_back({2, 1});

    // 2 Execution
    const auto mesh = meshkernel::Mesh2D(edges, nodes, meshkernel::Projection::cartesian);

    // 3 Validation
    ASSERT_EQ(2, mesh.GetNumFaces());
}

TEST(Mesh, MeshBoundaryToPolygon)
{
    // 1 Setup
    std::vector<meshkernel::Point> nodes;
    nodes.push_back({0.0, 0.0});
    nodes.push_back({5.0, -5.0});
    nodes.push_back({10.0, 0.0});
    nodes.push_back({5.0, 5.0});
    std::vector<meshkernel::Edge> edges;
    edges.push_back({0, 3});
    edges.push_back({0, 2});
    edges.push_back({2, 3});
    edges.push_back({0, 1});
    edges.push_back({2, 1});

    auto mesh = meshkernel::Mesh2D(edges, nodes, meshkernel::Projection::cartesian);

    std::vector<meshkernel::Point> polygonNodes;
    const auto meshBoundaryPolygon = mesh.ComputeBoundaryPolygons(polygonNodes);

    const double tolerance = 1e-5;
    ASSERT_NEAR(0.0, meshBoundaryPolygon[0].x, tolerance);
    ASSERT_NEAR(5.0, meshBoundaryPolygon[1].x, tolerance);
    ASSERT_NEAR(10.0, meshBoundaryPolygon[2].x, tolerance);
    ASSERT_NEAR(5.0, meshBoundaryPolygon[3].x, tolerance);
    ASSERT_NEAR(0.0, meshBoundaryPolygon[4].x, tolerance);

    ASSERT_NEAR(0.0, meshBoundaryPolygon[0].y, tolerance);
    ASSERT_NEAR(5.0, meshBoundaryPolygon[1].y, tolerance);
    ASSERT_NEAR(0.0, meshBoundaryPolygon[2].y, tolerance);
    ASSERT_NEAR(-5.0, meshBoundaryPolygon[3].y, tolerance);
    ASSERT_NEAR(0.0, meshBoundaryPolygon[4].y, tolerance);
}

TEST(Mesh, HangingEdge)
{
    // 1 Setup
    std::vector<meshkernel::Point> nodes;
    nodes.push_back({0.0, 0.0});
    nodes.push_back({5.0, 0.0});
    nodes.push_back({3.0, 2.0});
    nodes.push_back({3.0, 4.0});

    std::vector<meshkernel::Edge> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 3});
    edges.push_back({3, 0});
    edges.push_back({2, 1});

    auto mesh = meshkernel::Mesh2D(edges, nodes, meshkernel::Projection::cartesian);

    ASSERT_EQ(1, mesh.GetNumFaces());
}

TEST(Mesh, NodeMerging)
{
    // 1. Setup
    const int n = 10; // x
    const int m = 10; // y

    std::vector<std::vector<int>> indicesValues(n, std::vector<int>(m));
    std::vector<meshkernel::Point> nodes(n * m);
    meshkernel::UInt nodeIndex = 0;
    for (auto j = 0; j < m; ++j)
    {
        for (auto i = 0; i < n; ++i)
        {
            indicesValues[i][j] = i + j * n;
            nodes[nodeIndex] = {(double)i, (double)j};
            nodeIndex++;
        }
    }

    std::vector<meshkernel::Edge> edges((n - 1) * m + (m - 1) * n);
    std::size_t edgeIndex = 0;
    for (auto j = 0; j < m; ++j)
    {
        for (auto i = 0; i < n - 1; ++i)
        {
            edges[edgeIndex] = {indicesValues[i][j], indicesValues[i + 1][j]};
            edgeIndex++;
        }
    }

    for (auto j = 0; j < m - 1; ++j)
    {
        for (auto i = 0; i < n; ++i)
        {
            edges[edgeIndex] = {indicesValues[i][j + 1], indicesValues[i][j]};
            edgeIndex++;
        }
    }

    auto mesh = std::make_unique<meshkernel::Mesh2D>(edges, nodes, meshkernel::Projection::cartesian);

    // Add overlapping nodes
    double generatingDistance = std::sqrt(std::pow(0.001 * 0.9, 2) / 2.0);
    std::uniform_real_distribution<double> x_distribution(0.0, generatingDistance);
    std::uniform_real_distribution<double> y_distribution(0.0, generatingDistance);
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());

    nodes.resize(mesh->GetNumNodes() * 2);
    edges.resize(mesh->GetNumEdges() + mesh->GetNumNodes() * 2);
    meshkernel::UInt originalNodeIndex = 0;
    for (meshkernel::UInt j = 0; j < m; ++j)
    {
        for (meshkernel::UInt i = 0; i < n; ++i)
        {
            nodes[nodeIndex] = {i + x_distribution(generator), j + y_distribution(generator)};

            edges[edgeIndex] = {nodeIndex, originalNodeIndex};
            edgeIndex++;

            nodeIndex++;
            originalNodeIndex++;
        }
    }

    nodes.resize(nodeIndex);
    edges.resize(edgeIndex);

    // re set with augmented nodes
    mesh = std::make_unique<meshkernel::Mesh2D>(edges, nodes, meshkernel::Projection::cartesian);

    // 2. Act
    meshkernel::Polygons polygon;
    [[maybe_unused]] auto action = mesh->MergeNodesInPolygon(polygon, 0.001);

    // 3. Assert
    ASSERT_EQ(mesh->GetNumValidNodes(), n * m);
    ASSERT_EQ(mesh->GetNumValidEdges(), (n - 1) * m + (m - 1) * n);
}

TEST(Mesh, MillionQuads)
{
    const int n = 4; // x
    const int m = 4; // y

    std::vector<std::vector<meshkernel::UInt>> indicesValues(n, std::vector<meshkernel::UInt>(m));
    std::vector<meshkernel::Point> nodes(n * m);
    std::size_t nodeIndex = 0;
    for (meshkernel::UInt j = 0; j < m; ++j)
    {
        for (meshkernel::UInt i = 0; i < n; ++i)
        {
            indicesValues[i][j] = i + j * n;
            nodes[nodeIndex] = {static_cast<double>(i), static_cast<double>(j)};
            nodeIndex++;
        }
    }

    std::vector<meshkernel::Edge> edges((n - 1) * m + (m - 1) * n);
    std::size_t edgeIndex = 0;
    for (meshkernel::UInt j = 0; j < m; ++j)
    {
        for (meshkernel::UInt i = 0; i < n - 1; ++i)
        {
            edges[edgeIndex] = {indicesValues[i][j], indicesValues[i + 1][j]};
            edgeIndex++;
        }
    }

    for (meshkernel::UInt j = 0; j < m - 1; ++j)
    {
        for (meshkernel::UInt i = 0; i < n; ++i)
        {
            edges[edgeIndex] = {indicesValues[i][j + 1], indicesValues[i][j]};
            edgeIndex++;
        }
    }

    // now build node-edge mapping
    auto start(std::chrono::steady_clock::now());
    const auto mesh = meshkernel::Mesh2D(edges, nodes, meshkernel::Projection::cartesian);
    auto end(std::chrono::steady_clock::now());

    double elapsedTime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    // std::cout << "Elapsed time " << elapsedTime << " s " << std::endl;
    // std::cout << "Number of found cells " << mesh.GetNumFaces() << std::endl;

    EXPECT_LE(elapsedTime, 5.0);
}

TEST(Mesh, InsertNodeInMeshWithExistingNodesRtreeTriggersRTreeReBuild)
{
    // Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);
    mesh->BuildTree(meshkernel::Location::Nodes);

    // insert nodes modifies the number of nodes, m_nodesRTreeRequiresUpdate is set to true
    meshkernel::Point newPoint{10.0, 10.0};

    auto [newNodeIndex, insertAction] = mesh->InsertNode(newPoint);

    [[maybe_unused]] auto connectAction = mesh->ConnectNodes(0, newNodeIndex);

    // when m_nodesRTreeRequiresUpdate = true m_nodesRTree is not empty the mesh.m_nodesRTree is re-build
    mesh->Administrate();

    // builds edges RTree
    mesh->BuildTree(meshkernel::Location::Edges);
    const auto& rtreeEdges = mesh->GetRTree(meshkernel::Location::Edges);

    // even if m_edgesRTreeRequiresUpdate = true, m_edgesRTree is initially empty, so it is assumed that is not needed for searches
    ASSERT_EQ(5, rtreeEdges.Size());
}

TEST(Mesh, DeleteNodeInMeshWithExistingNodesRtreeTriggersRTreeReBuild)
{
    // Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    meshkernel::Point newPoint{10.0, 10.0};
    mesh->BuildTree(meshkernel::Location::Nodes);
    auto& rtree = mesh->GetRTree(meshkernel::Location::Nodes);

    [[maybe_unused]] auto [nodeId, indertAction] = mesh->InsertNode(newPoint);

    // delete nodes modifies the number of nodes, m_nodesRTreeRequiresUpdate is set to true
    [[maybe_unused]] auto deleteAction = mesh->DeleteNode(0);

    // when m_nodesRTreeRequiresUpdate
    mesh->DeleteInvalidNodesAndEdges();
    mesh->Administrate();

    // building a tree based on nodes
    rtree.BuildTree(mesh->Nodes());

    // After deleting a node, the nodes RTree is reduced
    ASSERT_EQ(3, rtree.Size());
}

TEST(Mesh, ConnectNodesInMeshWithExistingEdgesRtreeTriggersRTreeReBuild)
{
    // 1 Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    meshkernel::Point newPoint{10.0, 10.0};

    auto [newNodeIndex, insertAction] = mesh->InsertNode(newPoint);

    // connect nodes modifies the number of edges, m_nodesRTreeRequiresUpdate is set to true
    [[maybe_unused]] auto connectAction = mesh->ConnectNodes(0, newNodeIndex);

    // re-do mesh adminstration
    mesh->Administrate();

    // re-build tree
    mesh->BuildTree(meshkernel::Location::Edges);
    const auto& rtree = mesh->GetRTree(meshkernel::Location::Edges);

    // even if m_nodesRTreeRequiresUpdate = true, m_nodesRTree is initially empty, so it is assumed that is not needed for searches
    ASSERT_EQ(5, rtree.Size());
}

TEST(Mesh, DeleteEdgeInMeshWithExistingEdgesRtreeTriggersRTreeReBuild)
{
    // 1 Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);
    mesh->BuildTree(meshkernel::Location::Edges);
    const auto& rtree = mesh->GetRTree(meshkernel::Location::Edges);

    // DeleteEdge modifies the number of edges, m_edgesRTreeRequiresUpdate is set to true
    [[maybe_unused]] auto action = mesh->DeleteEdge(0);

    // re-do mesh administration
    mesh->Administrate();

    // re-build tree
    mesh->BuildTree(meshkernel::Location::Edges);

    // deleting an edge produces an edges rtree of size 3
    ASSERT_EQ(3, rtree.Size());
}

TEST(Mesh, InsertUnconnectedNodeInMeshIsSetToInvalid)
{
    // Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);
    mesh->BuildTree(meshkernel::Location::Nodes);
    const auto& rtreesNodes = mesh->GetRTree(meshkernel::Location::Nodes);

    // insert nodes modifies the number of nodes, m_nodesRTreeRequiresUpdate is set to true
    meshkernel::Point newPoint{10.0, 10.0};

    [[maybe_unused]] auto [nodeId, action] = mesh->InsertNode(newPoint);

    // when m_nodesRTreeRequiresUpdate = true m_nodesRTree is not empty the mesh.m_nodesRTree is re-build
    mesh->Administrate();

    // building a tree based on nodes
    mesh->BuildTree(meshkernel::Location::Nodes);

    // building the edges rtree
    mesh->BuildTree(meshkernel::Location::Edges);
    const auto& rtreeEdges = mesh->GetRTree(meshkernel::Location::Edges);

    // building a tree based on edges
    EXPECT_EQ(5, mesh->GetNumNodes());
    EXPECT_EQ(4, mesh->GetNumValidNodes());
    EXPECT_EQ(4, rtreesNodes.Size());
    EXPECT_EQ(4, rtreeEdges.Size());
    // Administrate should set the unconnected node to be invalid.
    EXPECT_FALSE(mesh->Node(4).IsValid());
}

TEST(Mesh, EdgeConnectedToInvalidNodeInMeshIsSetToInvalid)
{
    // Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);
    mesh->BuildTree(meshkernel::Location::Nodes);
    const auto& nodesRtree = mesh->GetRTree(meshkernel::Location::Nodes);
    const auto& edgesRtree = mesh->GetRTree(meshkernel::Location::Edges);

    meshkernel::Point newPoint{meshkernel::constants::missing::doubleValue,
                               meshkernel::constants::missing::doubleValue};
    auto [nodeIndex, nodeAction] = mesh->InsertNode(newPoint);
    auto [edgeIndex, edgeAction] = mesh->ConnectNodes(0, nodeIndex);

    EXPECT_EQ(mesh->GetEdge(edgeIndex).first, 0);
    EXPECT_EQ(mesh->GetEdge(edgeIndex).second, nodeIndex);

    // when m_nodesRTreeRequiresUpdate = true m_nodesRTree is not empty the mesh.m_nodesRTree is re-build
    mesh->Administrate();

    // building a tree based on nodes
    mesh->BuildTree(meshkernel::Location::Nodes);

    // building a tree based on edges
    mesh->BuildTree(meshkernel::Location::Edges);

    EXPECT_EQ(5, mesh->GetNumNodes());
    EXPECT_EQ(4, mesh->GetNumValidNodes());

    EXPECT_EQ(5, mesh->GetNumEdges());
    EXPECT_EQ(4, mesh->GetNumValidEdges());

    EXPECT_EQ(4, nodesRtree.Size());
    EXPECT_EQ(4, edgesRtree.Size());

    // Administrate should set the unconnected node to be invalid.
    EXPECT_FALSE(mesh->Node(4).IsValid());

    // Administrate should set the edge connecting an invalid node to be invalid.
    EXPECT_EQ(mesh->GetEdge(4).first, meshkernel::constants::missing::uintValue);
    EXPECT_EQ(mesh->GetEdge(4).second, meshkernel::constants::missing::uintValue);
}

TEST(Mesh, GetNodeIndexShouldTriggerNodesRTreeBuild)
{
    // 1 Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    // By default, no nodesRTree is build
    const auto& edgesRTree = mesh->GetRTree(meshkernel::Location::Edges);
    const auto& nodesRTree = mesh->GetRTree(meshkernel::Location::Nodes);

    ASSERT_EQ(0, nodesRTree.Size());
    ASSERT_EQ(0, edgesRTree.Size());

    // FindNodeCloseToAPoint builds m_nodesRTree for searching the nodes
    const auto index = mesh->FindNodeCloseToAPoint({1.5, 1.5}, 10.0);
    ASSERT_EQ(3, index);

    // m_nodesRTree is build
    ASSERT_EQ(4, nodesRTree.Size());

    // m_edgesRTree is not build when searching for nodes
    ASSERT_EQ(0, edgesRTree.Size());
}

TEST(Mesh, FindEdgeCloseToAPointShouldTriggerEdgesRTreeBuild)
{
    // 1 Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    // FindEdgeCloseToAPoint builds m_edgesRTree for searching the edges

    mesh->BuildTree(meshkernel::Location::Edges);
    const auto& edgesRTree = mesh->GetRTree(meshkernel::Location::Edges);
    const auto& nodesRTree = mesh->GetRTree(meshkernel::Location::Nodes);

    const auto index = mesh->FindLocationIndex({1.5, 1.5}, meshkernel::Location::Edges);
    ASSERT_EQ(1, index);

    // m_nodesRTree is not build when searching for edges
    ASSERT_EQ(0, nodesRTree.Size());

    // m_edgesRTree is build
    ASSERT_EQ(4, edgesRTree.Size());
}

TEST(Mesh, GetObtuseTriangles)
{
    // Setup a mesh with two triangles, one obtuse
    std::vector<meshkernel::Point> nodes{
        {0.0, 0.0},
        {3.0, 0.0},
        {-1.0, 2.0},
        {1.5, -2.0}};

    std::vector<meshkernel::Edge> edges{
        {0, 1},
        {1, 2},
        {2, 0},
        {0, 3},
        {3, 1}};

    auto mesh = std::make_unique<meshkernel::Mesh2D>(edges, nodes, meshkernel::Projection::cartesian);

    // execute, only one obtuse triangle should be found
    const auto obtuseTrianglesCount = mesh->GetObtuseTrianglesCenters().size();

    // assert a small flow edge is found
    ASSERT_EQ(1, obtuseTrianglesCount);
}

TEST(Mesh, GetSmallFlowEdgeCenters)
{
    // Setup a mesh with two triangles
    std::vector<meshkernel::Point> nodes{
        {0.0, 0.0},
        {1.0, 0.0},
        {1.0, 0.3},
        {1.0, -0.3}};

    std::vector<meshkernel::Edge> edges{
        {0, 3},
        {3, 1},
        {1, 0},
        {1, 2},
        {2, 0},
    };

    meshkernel::Mesh2D mesh = meshkernel::Mesh2D(edges, nodes, meshkernel::Projection::cartesian);

    // execute, by setting the smallFlowEdgesThreshold high, a small flow edge will be found
    const auto numSmallFlowEdgeFirstQuery = mesh.GetEdgesCrossingSmallFlowEdges(100).size();

    // execute, by setting the smallFlowEdgesThreshold low, no small flow edge will be found
    const auto numSmallFlowEdgeSecondQuery = mesh.GetEdgesCrossingSmallFlowEdges(0.0).size();

    // assert a small flow edge is found
    ASSERT_EQ(1, numSmallFlowEdgeFirstQuery);
    ASSERT_EQ(0, numSmallFlowEdgeSecondQuery);
}

TEST(Mesh, DeleteSmallFlowEdge)
{
    // Setup a mesh with eight triangles
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/RemoveSmallFlowEdgesTests/remove_small_flow_edges_net.nc");

    ASSERT_EQ(8, mesh->GetNumFaces());

    // After merging the number of faces is reduced
    auto undoAction = mesh->DeleteSmallFlowEdges(1.0);

    ASSERT_EQ(3, mesh->GetNumFaces());

    // Test the undo action has been computed correctly
    undoAction->Restore();
    // Recompute faces
    mesh->Administrate();
    ASSERT_EQ(8, mesh->GetNumFaces());
}

TEST(Mesh, DeleteSmallTrianglesAtBoundaries)
{
    // Setup a mesh with two triangles
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/RemoveSmallFlowEdgesTests/remove_small_flow_edges_quad_net.nc");

    ASSERT_EQ(2, mesh->GetNumFaces());

    // After merging
    auto undoAction = mesh->DeleteSmallTrianglesAtBoundaries(0.6);

    ASSERT_EQ(1, mesh->GetNumFaces());

    const double tolerance = 1e-8;
    ASSERT_NEAR(364.17013549804688, mesh->Node(0).x, tolerance);
    ASSERT_NEAR(meshkernel::constants::missing::doubleValue, mesh->Node(1).x, tolerance);
    ASSERT_NEAR(295.21142578125000, mesh->Node(2).x, tolerance);
    ASSERT_NEAR(421.46209716796875, mesh->Node(3).x, tolerance);
    ASSERT_NEAR(359.79510498046875, mesh->Node(4).x, tolerance);

    ASSERT_NEAR(374.00662231445313, mesh->Node(0).y, tolerance);
    ASSERT_NEAR(meshkernel::constants::missing::doubleValue, mesh->Node(1).y, tolerance);
    ASSERT_NEAR(300.48181152343750, mesh->Node(2).y, tolerance);
    ASSERT_NEAR(295.33038330078125, mesh->Node(3).y, tolerance);
    ASSERT_NEAR(398.59295654296875, mesh->Node(4).y, tolerance);

    // Test the undo action has been computed correctly
    undoAction->Restore();
    // Recompute faces
    mesh->Administrate();
    ASSERT_EQ(2, mesh->GetNumFaces());
}

TEST(Mesh, DeleteHangingEdge)
{
    // 1 Setup
    std::vector<meshkernel::Point> nodes;
    nodes.push_back({0.0, 0.0});
    nodes.push_back({5.0, 0.0});
    nodes.push_back({3.0, 4.0});

    std::vector<meshkernel::Edge> edges;
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    edges.push_back({2, 0});

    // Execute
    auto mesh = meshkernel::Mesh2D(edges, nodes, meshkernel::Projection::cartesian);

    // Add new node and connect with existing node, creating hanging edge
    nodes.push_back({3.0, 2.0});
    edges.push_back({3, 1});

    [[maybe_unused]] auto undoInsertNode = mesh.InsertNode(nodes[3]);
    [[maybe_unused]] auto undoConnectNodes = mesh.ConnectNodes(3, 1);

    // Assert
    ASSERT_EQ(1, mesh.GetNumFaces());
    ASSERT_EQ(4, mesh.GetNumEdges());

    // Execute
    auto hangingEdges = mesh.GetHangingEdges();

    // Assert
    ASSERT_EQ(1, hangingEdges.size());

    // Execute
    auto undoAction = mesh.DeleteHangingEdges();
    hangingEdges = mesh.GetHangingEdges();

    // Assert
    ASSERT_EQ(0, hangingEdges.size());

    // Test the undo action has been computed correctly
    undoAction->Restore();
    // Recompute faces
    mesh.Administrate();
    // Assert
    ASSERT_EQ(1, mesh.GetNumFaces());
    ASSERT_EQ(4, mesh.GetNumEdges());
}

class MeshDeletion : public ::testing::TestWithParam<std::tuple<meshkernel::Mesh2D::DeleteMeshOptions, bool, int>>
{
public:
    [[nodiscard]] static std::vector<std::tuple<meshkernel::Mesh2D::DeleteMeshOptions, bool, int>> GetData()
    {
        return {
            {meshkernel::Mesh2D::DeleteMeshOptions::InsideNotIntersected, false, 16},
            {meshkernel::Mesh2D::DeleteMeshOptions::InsideAndIntersected, false, 14},
            {meshkernel::Mesh2D::DeleteMeshOptions::InsideNotIntersected, true, 0},
            {meshkernel::Mesh2D::DeleteMeshOptions::InsideAndIntersected, true, 6}

        };
    }
};

TEST_P(MeshDeletion, expected_results)
{
    // Get the test parameters
    auto const [deleteOption, invertSelection, numNodes] = GetParam();

    // Setup
    auto mesh = MakeRectangularMeshForTesting(4, 4, 1.0, meshkernel::Projection::cartesian);

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    std::vector<meshkernel::Point> polygonNodes{
        {-0.5, -1.0},
        {0.8, -1.0},
        {0.8, 1.8},
        {-0.5, 1.8},
        {-0.5, -1.0}};

    const meshkernel::Polygons polygon(polygonNodes, meshkernel::Projection::cartesian);

    // Execute
    auto undoAction = mesh->DeleteMesh(polygon, deleteOption, invertSelection);

    // Assert
    ASSERT_EQ(numNodes, mesh->GetNumValidNodes());

    // Test the undo action has been computed correctly
    undoAction->Restore();
    // Recompute faces
    mesh->Administrate();

    ASSERT_EQ(originalNodes.size(), mesh->Nodes().size());
    ASSERT_EQ(originalEdges.size(), mesh->Edges().size());

    for (meshkernel::UInt i = 0; i < mesh->Nodes().size(); ++i)
    {
        EXPECT_EQ(originalNodes[i].x, mesh->Node(i).x);
        EXPECT_EQ(originalNodes[i].y, mesh->Node(i).y);
    }

    for (meshkernel::UInt i = 0; i < mesh->Edges().size(); ++i)
    {
        EXPECT_EQ(originalEdges[i].first, mesh->GetEdge(i).first) << "first node of edge " << i;
        EXPECT_EQ(originalEdges[i].second, mesh->GetEdge(i).second) << "second node of edge " << i;
    }
}

INSTANTIATE_TEST_SUITE_P(Mesh, MeshDeletion, ::testing::ValuesIn(MeshDeletion::GetData()));

class MeshDeletionWithInnerPolygons : public ::testing::TestWithParam<std::tuple<meshkernel::Mesh2D::DeleteMeshOptions, bool, std::vector<meshkernel::Point>, int>>
{

    static inline std::vector<meshkernel::Point> firstPolygon_{
        {-0.5, -0.5},
        {7.5, -0.5},
        {7.5, 7.5},
        {-0.5, 7.5},
        {-0.5, -0.5},
        {meshkernel::constants::missing::innerOuterSeparator, meshkernel::constants::missing::innerOuterSeparator},
        {1.5, 1.5},
        {4.5, 1.5},
        {4.5, 4.5},
        {1.5, 4.5},
        {1.5, 1.5},
    };

    static inline std::vector<meshkernel::Point> secondPolygon_{
        {-0.5, -0.5},
        {7.5, -0.5},
        {7.5, 7.5},
        {-0.5, 7.5},
        {-0.5, -0.5},
        {meshkernel::constants::missing::innerOuterSeparator, meshkernel::constants::missing::innerOuterSeparator},
        {1.5, 1.5},
        {4.5, 1.5},
        {4.5, 4.5},
        {2.7, 4.5},
        {2.7, 3.3},
        {1.5, 3.3},
        {1.5, 1.5}};

public:
    [[nodiscard]] static std::vector<std::tuple<meshkernel::Mesh2D::DeleteMeshOptions, bool, std::vector<meshkernel::Point>, int>> GetData()
    {
        return {
            {meshkernel::Mesh2D::DeleteMeshOptions::InsideAndIntersected, false, firstPolygon_, 9},
            {meshkernel::Mesh2D::DeleteMeshOptions::InsideAndIntersected, true, firstPolygon_, 48},
            {meshkernel::Mesh2D::DeleteMeshOptions::InsideAndIntersected, false, secondPolygon_, 8},
            {meshkernel::Mesh2D::DeleteMeshOptions::InsideAndIntersected, true, secondPolygon_, 49}};
    }
};

TEST_P(MeshDeletionWithInnerPolygons, expected_results)
{
    // Get the test parameters
    auto const& [deleteOption, invertSelection, polygonNodes, numNodes] = GetParam();

    // Setup
    auto mesh = MakeRectangularMeshForTesting(7, 7, 1.0, meshkernel::Projection::cartesian);

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    const meshkernel::Polygons polygon(polygonNodes, meshkernel::Projection::cartesian);

    // Execute
    auto undoAction = mesh->DeleteMesh(polygon, deleteOption, invertSelection);

    // Assert
    const auto nodes = mesh->Nodes();
    ASSERT_EQ(numNodes, mesh->GetNumValidNodes());

    // Test the undo action has been computed correctly
    undoAction->Restore();
    // Recompute faces
    mesh->Administrate();

    ASSERT_EQ(originalNodes.size(), mesh->Nodes().size());
    ASSERT_EQ(originalEdges.size(), mesh->Edges().size());

    for (meshkernel::UInt i = 0; i < mesh->Nodes().size(); ++i)
    {
        EXPECT_EQ(originalNodes[i].x, mesh->Node(i).x);
        EXPECT_EQ(originalNodes[i].y, mesh->Node(i).y);
    }

    for (meshkernel::UInt i = 0; i < mesh->Edges().size(); ++i)
    {
        EXPECT_EQ(originalEdges[i].first, mesh->GetEdge(i).first);
        EXPECT_EQ(originalEdges[i].second, mesh->GetEdge(i).second);
    }
}

INSTANTIATE_TEST_SUITE_P(Mesh, MeshDeletionWithInnerPolygons, ::testing::ValuesIn(MeshDeletionWithInnerPolygons::GetData()));

TEST(Mesh, ManageMeshHoles)
{
    // Setup
    auto mesh = MakeRectangularMeshForTesting(4, 4, 1.0, meshkernel::Projection::cartesian);
    mesh->BuildTree(meshkernel::Location::Nodes);

    std::vector<meshkernel::Point> polygonPoints{{0.99, 0.99}, {2.01, 0.99}, {2.01, 2.01}, {0.99, 2.01}, {0.99, 0.99}};
    meshkernel::Polygons polygon(polygonPoints, mesh->m_projection);

    // Create holes in mesh
    auto undoDeleteFaces = mesh->DeleteMeshFacesInPolygon(polygon);

    meshkernel::Polygons refinementPolygon;

    meshkernel::MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 1;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;

    meshkernel::MeshRefinement meshRefinement(*mesh, refinementPolygon, meshRefinementParameters);
    auto refinementUndoAction = meshRefinement.Compute();

    ASSERT_EQ(mesh->GetNumNodes(), 48);
    ASSERT_EQ(mesh->GetNumEdges(), 80);
    ASSERT_EQ(mesh->GetNumFaces(), 32);

    std::vector<double> expectedNodeX{0.0, 0.0, 0.0, 0.0, 1.0, 1.0,
                                      1.0, 1.0, 2.0, 2.0, 2.0, 2.0,
                                      3.0, 3.0, 3.0, 3.0, 0.5, 0.5,
                                      0.5, 0.5, 1.5, 1.5, 1.5, 1.5,
                                      2.5, 2.5, 2.5, 2.5, 0.0, 0.0,
                                      0.0, 1.0, 1.0, 1.0, 2.0, 2.0,
                                      2.0, 3.0, 3.0, 3.0, 0.5, 0.5,
                                      0.5, 1.5, 1.5, 2.5, 2.5, 2.5};

    std::vector<double> expectedNodeY{0.0, 1.0, 2.0, 3.0, 0.0, 1.0,
                                      2.0, 3.0, 0.0, 1.0, 2.0, 3.0,
                                      0.0, 1.0, 2.0, 3.0, 0.0, 1.0,
                                      2.0, 3.0, 0.0, 1.0, 2.0, 3.0,
                                      0.0, 1.0, 2.0, 3.0, 0.5, 1.5,
                                      2.5, 0.5, 1.5, 2.5, 0.5, 1.5,
                                      2.5, 0.5, 1.5, 2.5, 0.5, 1.5,
                                      2.5, 0.5, 2.5, 0.5, 1.5, 2.5};

    std::vector<meshkernel::UInt> expectedEdgesFirst{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                                                     11, 1, 2, 3, 5, 6, 7, 9, 10, 11,
                                                     13, 14, 15, 16, 31, 17, 28, 17,
                                                     32, 18, 29, 18, 33, 19, 30, 20,
                                                     34, 21, 31, 22, 36, 23, 33, 24,
                                                     37, 25, 34, 25, 38, 26, 35, 26,
                                                     39, 27, 36, 16, 17, 18, 19, 20,
                                                     21, 22, 23, 24, 25, 26, 27, 28,
                                                     29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39};

    std::vector<meshkernel::UInt> expectedEdgesSecond{16, 17, 18, 19, 20, 21, 22, 23,
                                                      24, 25, 26, 27, 28, 29, 30, 31,
                                                      32, 33, 34, 35, 36, 37, 38, 39,
                                                      40, 40, 40, 40, 41, 41, 41, 41,
                                                      42, 42, 42, 42, 43, 43, 43, 43,
                                                      44, 44, 44, 44, 45, 45, 45, 45,
                                                      46, 46, 46, 46, 47, 47, 47, 47,
                                                      4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
                                                      14, 15, 0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14};

    auto nodes = mesh->Nodes();
    auto edges = mesh->Edges();

    for (size_t i = 0; i < nodes.size(); ++i)
    {
        EXPECT_EQ(expectedNodeX[i], nodes[i].x);
        EXPECT_EQ(expectedNodeY[i], nodes[i].y);
    }

    for (size_t i = 0; i < edges.size(); ++i)
    {
        EXPECT_EQ(expectedEdgesFirst[i], edges[i].first);
        EXPECT_EQ(expectedEdgesSecond[i], edges[i].second);
    }

    refinementUndoAction->Restore();
    undoDeleteFaces->Restore();
    mesh->Administrate();

    ASSERT_EQ(mesh->GetNumValidNodes(), 16);
    ASSERT_EQ(mesh->GetNumValidEdges(), 24);
    ASSERT_EQ(mesh->GetNumFaces(), 9);
}
