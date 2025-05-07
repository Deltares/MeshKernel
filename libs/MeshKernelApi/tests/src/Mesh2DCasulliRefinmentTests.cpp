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
#include <gtest/gtest.h>
#include <random>

#include "CartesianApiTestFixture.hpp"
#include "MeshKernel/Parameters.hpp"
#include "MeshKernelApi/BoundingBox.hpp"
#include "MeshKernelApi/Mesh2D.hpp"
#include "MeshKernelApi/MeshKernel.hpp"
#include "TestUtils/MakeMeshes.hpp"

// namespace aliases
namespace mk = meshkernel;
namespace mkapi = meshkernelapi;

TEST(CasulliRefinement, RefineGridUsingApi_CasulliRefinement_ShouldRefine)
{
    // Prepare
    int meshKernelId = -1;
    const int projectionType = 1;

    auto errorCode = meshkernelapi::mkernel_allocate_state(projectionType, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::MakeGridParameters makeGridParameters;

    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = 10.0;
    makeGridParameters.block_size_y = 10.0;
    makeGridParameters.num_columns = 11;
    makeGridParameters.num_rows = 11;

    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_convert_to_mesh2d(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Refine using Casulli algorithm
    errorCode = meshkernelapi::mkernel_mesh2d_casulli_refinement(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    // Check number of nodes and edges is correct.
    EXPECT_EQ(576, mesh2d.num_valid_nodes);
    EXPECT_EQ(1104, mesh2d.num_valid_edges);
}

TEST(CasulliRefinement, ConnectingTwoMeshesAfterCasulliRefinementDoesNotCrash)
{
    // This test is really just to check that the sequence of steps involved does not cause a crash.

    // Prepare
    int meshKernelId1;
    int meshKernelId2;
    constexpr int isSpherical = 0;
    int errorCode = meshkernelapi::mkernel_allocate_state(isSpherical, meshKernelId1);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_allocate_state(isSpherical, meshKernelId2);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::MakeGridParameters gridParameters;
    gridParameters.num_columns = 3;
    gridParameters.num_rows = 3;
    gridParameters.block_size_x = 10.0;
    gridParameters.block_size_y = 10.0;
    gridParameters.origin_x = 0.0;
    gridParameters.origin_y = 0.0;
    gridParameters.angle = 0.0;

    errorCode = meshkernelapi::mkernel_mesh2d_make_rectangular_mesh(meshKernelId1, gridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    gridParameters.origin_x = 30.5;
    gridParameters.origin_y = 0.0;

    errorCode = meshkernelapi::mkernel_mesh2d_make_rectangular_mesh(meshKernelId2, gridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    //--------------------------------

    std::vector<double> polygonPointsX({-5.0, 15.0, 15.0, -5.0, -5.0});
    std::vector<double> polygonPointsY({15.0, 15.0, -5.0, -5.0, 15.0});
    meshkernelapi::GeometryList polygon;
    polygon.num_coordinates = 5;
    polygon.coordinates_x = polygonPointsX.data();
    polygon.coordinates_y = polygonPointsY.data();

    errorCode = meshkernelapi::mkernel_mesh2d_casulli_refinement_on_polygon(meshKernelId1, polygon);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    //--------------------------------

    meshkernelapi::Mesh2D mesh2d1{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId2, mesh2d1);

    std::vector<double> node_x(mesh2d1.num_nodes);
    std::vector<double> node_y(mesh2d1.num_nodes);
    std::vector<int> edge_nodes(mesh2d1.num_edges * 2);

    mesh2d1.node_x = node_x.data();
    mesh2d1.node_y = node_y.data();
    mesh2d1.edge_nodes = edge_nodes.data();

    errorCode = mkernel_mesh2d_get_node_edge_data(meshKernelId2, mesh2d1);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Connect the two meshes
    // to note: the mesh associated with meshKernelId1 has gaps in the node and edge data arrays (with invalid data)
    // after the Casulli refinement. The mesh data, mesh2d1, does not have any such gaps.
    // The mesh with gaps in the data must appear first in the mesh connection.
    errorCode = meshkernelapi::mkernel_mesh2d_connect_meshes(meshKernelId1, mesh2d1, 0.1, true);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    //--------------------------------

    meshkernelapi::Mesh2D mesh2d2{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId1, mesh2d2);

    std::vector<double> node_x2(mesh2d2.num_nodes);
    std::vector<double> node_y2(mesh2d2.num_nodes);
    std::vector<int> edge_nodes2(mesh2d2.num_edges * 2);

    mesh2d2.node_x = node_x2.data();
    mesh2d2.node_y = node_y2.data();
    mesh2d2.edge_nodes = edge_nodes2.data();

    errorCode = mkernel_mesh2d_get_node_edge_data(meshKernelId1, mesh2d2);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    EXPECT_EQ(mesh2d2.num_nodes, 44);
    EXPECT_EQ(mesh2d2.num_edges, 72);
}

TEST(CasulliRefinement, CasulliRefinementErrorCases)
{
    // Prepare
    int meshKernelId;
    const int isGeographic = 0;
    auto errorCode = meshkernelapi::mkernel_allocate_state(isGeographic, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_casulli_refinement(meshKernelId + 1);
    ASSERT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_casulli_refinement_on_polygon(meshKernelId + 1,
                                                                            meshkernelapi::GeometryList{});
    ASSERT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    meshkernelapi::mkernel_deallocate_state(meshKernelId);
}

TEST(CasulliRefinement, CasulliRefinementWholeMesh)
{
    // Prepare
    int meshKernelId;
    const int isGeographic = 0;
    meshkernelapi::mkernel_allocate_state(isGeographic, meshKernelId);

    // Set-up new mesh
    auto [num_nodes, num_edges, node_x, node_y, edge_nodes] = MakeRectangularMeshForApiTesting(10, 10, 1.0);

    std::vector<double> originalNodesX(node_x);
    std::vector<double> originalNodesY(node_y);
    std::vector<int> originalEdges(edge_nodes);

    meshkernelapi::Mesh2D mesh2d{};
    mesh2d.num_edges = static_cast<int>(num_edges);
    mesh2d.num_nodes = static_cast<int>(num_nodes);
    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_nodes = edge_nodes.data();

    auto errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_casulli_refinement(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D refinedMesh2d{};

    // Just do a rudimentary check that the number of noeds and edges is greater in the refined mesh.

    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, refinedMesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    EXPECT_GT(refinedMesh2d.num_nodes, mesh2d.num_nodes);
    EXPECT_GT(refinedMesh2d.num_edges, mesh2d.num_edges);
}

TEST(CasulliRefinement, CasulliRefinementMeshRegion)
{
    // Prepare
    int meshKernelId;
    const int isGeographic = 0;
    meshkernelapi::mkernel_allocate_state(isGeographic, meshKernelId);

    // Set-up new mesh
    auto [num_nodes, num_edges, node_x, node_y, edge_nodes] = MakeRectangularMeshForApiTesting(10, 10, 1.0);

    std::vector<double> originalNodesX(node_x);
    std::vector<double> originalNodesY(node_y);
    std::vector<int> originalEdges(edge_nodes);

    meshkernelapi::Mesh2D mesh2d{};
    mesh2d.num_edges = static_cast<int>(num_edges);
    mesh2d.num_nodes = static_cast<int>(num_nodes);
    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_nodes = edge_nodes.data();

    std::vector<double> polygonPointsX({2.5, 7.5, 5.5, 2.5});
    std::vector<double> polygonPointsY({2.5, 4.5, 8.5, 2.5});
    meshkernelapi::GeometryList polygon;
    polygon.num_coordinates = 4;
    polygon.coordinates_x = polygonPointsX.data();
    polygon.coordinates_y = polygonPointsY.data();

    auto errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_casulli_refinement_on_polygon(meshKernelId, polygon);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D refinedMesh2d{};

    // Just do a rudimentary check that the number of noeds and edges is greater in the refined mesh.

    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, refinedMesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    EXPECT_GT(refinedMesh2d.num_nodes, mesh2d.num_nodes);
    EXPECT_GT(refinedMesh2d.num_edges, mesh2d.num_edges);
}

TEST(CasulliRefinement, CasulliDeRefinementErrorCases)
{
    // Prepare
    int meshKernelId;
    const int isGeographic = 0;
    auto errorCode = meshkernelapi::mkernel_allocate_state(isGeographic, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_casulli_derefinement(meshKernelId + 1);
    ASSERT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_casulli_derefinement_on_polygon(meshKernelId + 1,
                                                                              meshkernelapi::GeometryList{});
    ASSERT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    meshkernelapi::mkernel_deallocate_state(meshKernelId);
}

TEST(CasulliRefinement, CasulliDeRefinementWholeMesh)
{
    // Prepare
    int meshKernelId;
    const int isGeographic = 0;
    meshkernelapi::mkernel_allocate_state(isGeographic, meshKernelId);

    // Set-up new mesh
    auto [num_nodes, num_edges, node_x, node_y, edge_nodes] = MakeRectangularMeshForApiTesting(20, 20, 1.0);

    std::vector<double> originalNodesX(node_x);
    std::vector<double> originalNodesY(node_y);
    std::vector<int> originalEdges(edge_nodes);

    std::vector<double> edgeCentresX(num_edges);
    std::vector<double> edgeCentresY(num_edges);
    std::vector<int> edgeFaces(2 * num_edges);
    std::vector<double> faceCentresX(20 * 20);
    std::vector<double> faceCentresY(20 * 20);
    std::vector<int> faceNodes(20 * 20 * 4);
    std::vector<int> faceEdges(20 * 20 * 4);
    std::vector<int> nodesPerFace(20 * 20);

    meshkernelapi::Mesh2D mesh2d{};
    mesh2d.num_edges = static_cast<int>(num_edges);
    mesh2d.num_nodes = static_cast<int>(num_nodes);
    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_nodes = edge_nodes.data();

    mesh2d.edge_x = edgeCentresX.data();
    mesh2d.edge_y = edgeCentresY.data();
    mesh2d.edge_faces = edgeFaces.data();
    mesh2d.face_x = faceCentresX.data();
    mesh2d.face_y = faceCentresY.data();
    mesh2d.face_nodes = faceNodes.data();
    mesh2d.face_edges = faceEdges.data();
    mesh2d.nodes_per_face = nodesPerFace.data();

    auto errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_casulli_derefinement(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D derefinedMesh2d{};

    // Just do a rudimentary check that there are fewer nodes and edges in the de-refined mesh.

    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, derefinedMesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    EXPECT_LT(derefinedMesh2d.num_nodes, mesh2d.num_nodes);
    EXPECT_LT(derefinedMesh2d.num_edges, mesh2d.num_edges);

    // Now check undo

    bool didUndo = false;
    int undoId = meshkernel::constants::missing::intValue;
    errorCode = meshkernelapi::mkernel_undo_state(didUndo, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndo);
    ASSERT_EQ(meshKernelId, undoId);

    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = mkernel_mesh2d_get_data(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    ASSERT_EQ(static_cast<int>(originalNodesX.size()), mesh2d.num_nodes);
    ASSERT_EQ(static_cast<int>(originalEdges.size()), 2 * mesh2d.num_edges);

    constexpr double tolerance = 1.0e-12;

    for (size_t i = 0; i < originalNodesX.size(); ++i)
    {
        EXPECT_NEAR(originalNodesX[i], mesh2d.node_x[i], tolerance);
        EXPECT_NEAR(originalNodesY[i], mesh2d.node_y[i], tolerance);
    }

    for (size_t i = 0; i < static_cast<size_t>(2 * mesh2d.num_edges); ++i)
    {
        EXPECT_EQ(originalEdges[i], mesh2d.edge_nodes[i]);
    }
}

TEST(CasulliRefinement, CasullDeRefinementMeshRegion)
{
    // Prepare
    int meshKernelId;
    const int isGeographic = 0;
    meshkernelapi::mkernel_allocate_state(isGeographic, meshKernelId);

    // Set-up new mesh
    auto [num_nodes, num_edges, node_x, node_y, edge_nodes] = MakeRectangularMeshForApiTesting(10, 10, 1.0);

    std::vector<double> originalNodesX(node_x);
    std::vector<double> originalNodesY(node_y);
    std::vector<int> originalEdges(edge_nodes);

    meshkernelapi::Mesh2D mesh2d{};
    mesh2d.num_edges = static_cast<int>(num_edges);
    mesh2d.num_nodes = static_cast<int>(num_nodes);
    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_nodes = edge_nodes.data();

    std::vector<double> polygonPointsX({2.5, 7.5, 5.5, 2.5});
    std::vector<double> polygonPointsY({2.5, 4.5, 8.5, 2.5});
    meshkernelapi::GeometryList polygon;
    polygon.num_coordinates = 4;
    polygon.coordinates_x = polygonPointsX.data();
    polygon.coordinates_y = polygonPointsY.data();

    auto errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_casulli_derefinement_on_polygon(meshKernelId, polygon);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D refinedMesh2d{};

    // Just do a rudimentary check that the number of noeds and edges is lesser in the refined mesh.

    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, refinedMesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    EXPECT_LT(refinedMesh2d.num_nodes, mesh2d.num_nodes);
    EXPECT_LT(refinedMesh2d.num_edges, mesh2d.num_edges);
}

TEST(CasulliRefinement, CasulliDeRefinementElementsWholeMesh)
{
    // Prepare
    int meshKernelId;
    const int isGeographic = 0;
    meshkernelapi::mkernel_allocate_state(isGeographic, meshKernelId);

    const int numElementsX = 10;
    const int numElementsY = 10;
    const int numElements = numElementsX * numElementsY;

    // Set-up new mesh
    auto [num_nodes, num_edges, node_x, node_y, edge_nodes] = MakeRectangularMeshForApiTesting(numElementsX, numElementsY, 1.0);

    std::vector<double> originalNodesX(node_x);
    std::vector<double> originalNodesY(node_y);
    std::vector<int> originalEdges(edge_nodes);

    std::vector<double> edgeCentresX(num_edges);
    std::vector<double> edgeCentresY(num_edges);
    std::vector<int> edgeFaces(2 * num_edges);
    std::vector<double> faceCentresX(numElements);
    std::vector<double> faceCentresY(numElements);
    std::vector<int> faceNodes(numElements * 4);
    std::vector<int> faceEdges(numElements * 4);
    std::vector<int> nodesPerFace(numElements);

    std::vector<double> removedElementCentresX(num_nodes);
    std::vector<double> removedElementCentresY(num_nodes);
    meshkernelapi::GeometryList elementsToRemove;

    elementsToRemove.coordinates_x = removedElementCentresX.data();
    elementsToRemove.coordinates_y = removedElementCentresY.data();

    meshkernelapi::Mesh2D mesh2d{};
    mesh2d.num_edges = static_cast<int>(num_edges);
    mesh2d.num_nodes = static_cast<int>(num_nodes);
    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_nodes = edge_nodes.data();

    mesh2d.edge_x = edgeCentresX.data();
    mesh2d.edge_y = edgeCentresY.data();
    mesh2d.edge_faces = edgeFaces.data();
    mesh2d.face_x = faceCentresX.data();
    mesh2d.face_y = faceCentresY.data();
    mesh2d.face_nodes = faceNodes.data();
    mesh2d.face_edges = faceEdges.data();
    mesh2d.nodes_per_face = nodesPerFace.data();

    auto errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_casulli_derefinement_elements(meshKernelId, elementsToRemove);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    const std::vector<meshkernel::Point> removedElementCentres{{1.5, 0.5}, {1.5, 2.5}, {1.5, 4.5}, {1.5, 6.5}, {1.5, 8.5}, {3.5, 0.5}, {3.5, 2.5}, {3.5, 4.5}, {3.5, 6.5}, {3.5, 8.5}, {5.5, 0.5}, {5.5, 2.5}, {5.5, 4.5}, {5.5, 6.5}, {5.5, 8.5}, {7.5, 0.5}, {7.5, 2.5}, {7.5, 4.5}, {7.5, 6.5}, {7.5, 8.5}, {9.5, 0.5}, {9.5, 2.5}, {9.5, 4.5}, {9.5, 6.5}, {9.5, 8.5}};

    ASSERT_EQ(elementsToRemove.num_coordinates, static_cast<int>(removedElementCentres.size()));
    constexpr double tolerance = 1.0e-12;

    for (int i = 0; i < elementsToRemove.num_coordinates; ++i)
    {
        EXPECT_NEAR(removedElementCentres[i].x, elementsToRemove.coordinates_x[i], tolerance);
        EXPECT_NEAR(removedElementCentres[i].y, elementsToRemove.coordinates_y[i], tolerance);
    }
}

TEST(CasulliRefinement, CasulliDeRefinementElementsMeshRegion)
{
    // Prepare
    int meshKernelId;
    const int isGeographic = 0;
    meshkernelapi::mkernel_allocate_state(isGeographic, meshKernelId);

    const int numElementsX = 10;
    const int numElementsY = 10;
    const int numElements = numElementsX * numElementsY;

    // Set-up new mesh
    auto [num_nodes, num_edges, node_x, node_y, edge_nodes] = MakeRectangularMeshForApiTesting(numElementsX, numElementsY, 1.0);

    std::vector<double> originalNodesX(node_x);
    std::vector<double> originalNodesY(node_y);
    std::vector<int> originalEdges(edge_nodes);

    std::vector<double> edgeCentresX(num_edges);
    std::vector<double> edgeCentresY(num_edges);
    std::vector<int> edgeFaces(2 * num_edges);
    std::vector<double> faceCentresX(numElements);
    std::vector<double> faceCentresY(numElements);
    std::vector<int> faceNodes(numElements * 4);
    std::vector<int> faceEdges(numElements * 4);
    std::vector<int> nodesPerFace(numElements);

    std::vector<double> polygonPointsX({2.5, 7.5, 5.5, 2.5});
    std::vector<double> polygonPointsY({2.5, 4.5, 8.5, 2.5});
    meshkernelapi::GeometryList polygon;
    polygon.num_coordinates = 4;
    polygon.coordinates_x = polygonPointsX.data();
    polygon.coordinates_y = polygonPointsY.data();

    std::vector<double> removedElementCentresX(num_nodes);
    std::vector<double> removedElementCentresY(num_nodes);
    meshkernelapi::GeometryList elementsToRemove;

    elementsToRemove.coordinates_x = removedElementCentresX.data();
    elementsToRemove.coordinates_y = removedElementCentresY.data();

    meshkernelapi::Mesh2D mesh2d{};
    mesh2d.num_edges = static_cast<int>(num_edges);
    mesh2d.num_nodes = static_cast<int>(num_nodes);
    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_nodes = edge_nodes.data();

    mesh2d.edge_x = edgeCentresX.data();
    mesh2d.edge_y = edgeCentresY.data();
    mesh2d.edge_faces = edgeFaces.data();
    mesh2d.face_x = faceCentresX.data();
    mesh2d.face_y = faceCentresY.data();
    mesh2d.face_nodes = faceNodes.data();
    mesh2d.face_edges = faceEdges.data();
    mesh2d.nodes_per_face = nodesPerFace.data();

    auto errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_casulli_derefinement_elements_on_polygon(meshKernelId, polygon, elementsToRemove);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    const std::vector<meshkernel::Point> removedElementCentres{{2.5, 2.5}, {4.5, 4.5}, {4.5, 6.5}, {6.5, 4.5}, {6.5, 6.5}};

    ASSERT_EQ(elementsToRemove.num_coordinates, static_cast<int>(removedElementCentres.size()));
    constexpr double tolerance = 1.0e-12;

    for (int i = 0; i < elementsToRemove.num_coordinates; ++i)
    {
        EXPECT_NEAR(removedElementCentres[i].x, elementsToRemove.coordinates_x[i], tolerance);
        EXPECT_NEAR(removedElementCentres[i].y, elementsToRemove.coordinates_y[i], tolerance);
    }
}

TEST(CasulliRefinement, CurvilinearFullMeshRefinement)
{
    int meshKernelId = -1;
    int errorCode = meshkernel::ExitCode::Success;

    const int mRefinement = 2;
    const int nRefinement = 3;

    const double tolerance = 1.0e-10;
    const double deltaX = 2.0;
    const double deltaY = 3.0;
    const double origin = 0.0;

    const double refinedDeltaX = deltaX / static_cast<double>(mRefinement);
    const double refinedDeltaY = deltaY / static_cast<double>(nRefinement);

    const int isGeographic = 0;
    errorCode = meshkernelapi::mkernel_allocate_state(isGeographic, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::MakeGridParameters makeMeshParameters;
    makeMeshParameters.num_columns = 3;
    makeMeshParameters.num_rows = 3;
    makeMeshParameters.block_size_x = deltaX;
    makeMeshParameters.block_size_y = deltaY;
    makeMeshParameters.origin_x = origin;
    makeMeshParameters.origin_y = origin;
    makeMeshParameters.angle = 0.0;

    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeMeshParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    //--------------------------------
    // Compute refinement

    // Should succeed
    errorCode = meshkernelapi::mkernel_curvilinear_full_refine(meshKernelId, mRefinement, nRefinement);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGridOut{};
    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGridOut);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    ASSERT_EQ(7, curvilinearGridOut.num_m);
    ASSERT_EQ(10, curvilinearGridOut.num_n);

    std::vector<double> xNodes(curvilinearGridOut.num_m * curvilinearGridOut.num_n);
    std::vector<double> yNodes(curvilinearGridOut.num_m * curvilinearGridOut.num_n);

    curvilinearGridOut.node_x = xNodes.data();
    curvilinearGridOut.node_y = yNodes.data();

    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, curvilinearGridOut);

    size_t count = 0;

    double expectedY = origin;

    for (int n = 0; n < curvilinearGridOut.num_n; ++n)
    {
        double expectedX = origin;

        for (int m = 0; m < curvilinearGridOut.num_m; ++m)
        {
            EXPECT_NEAR(expectedX, xNodes[count], tolerance);
            EXPECT_NEAR(expectedY, yNodes[count], tolerance);
            ++count;
            expectedX += refinedDeltaX;
        }

        expectedY += refinedDeltaY;
    }

    errorCode = meshkernelapi::mkernel_deallocate_state(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}

TEST(CasulliRefinement, CurvilinearFullMeshRefinementFailureTests)
{
    int meshKernelId = -1;
    int errorCode = meshkernel::ExitCode::Success;

    const double deltaX = 2.0;
    const double deltaY = 3.0;
    const double origin = 0.0;

    // Should fail, meshkernelId not defined
    errorCode = meshkernelapi::mkernel_curvilinear_full_refine(meshKernelId, 1, 2);
    ASSERT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    const int isGeographic = 0;
    errorCode = meshkernelapi::mkernel_allocate_state(isGeographic, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Should fail, curvilinear grid not valid (i.e. not grid.IsValid)
    errorCode = meshkernelapi::mkernel_curvilinear_full_refine(meshKernelId, 2, 3);
    ASSERT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    meshkernel::MakeGridParameters makeMeshParameters;
    makeMeshParameters.num_columns = 3;
    makeMeshParameters.num_rows = 3;
    makeMeshParameters.block_size_x = deltaX;
    makeMeshParameters.block_size_y = deltaY;
    makeMeshParameters.origin_x = origin;
    makeMeshParameters.origin_y = origin;
    makeMeshParameters.angle = 0.0;

    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeMeshParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    //--------------------------------

    // Should fail, invalid refinement factor
    errorCode = meshkernelapi::mkernel_curvilinear_full_refine(meshKernelId, -1, 0);
    ASSERT_EQ(meshkernel::ExitCode::ConstraintErrorCode, errorCode);
    errorCode = meshkernelapi::mkernel_curvilinear_full_refine(meshKernelId, 0, -2);
    ASSERT_EQ(meshkernel::ExitCode::ConstraintErrorCode, errorCode);
    errorCode = meshkernelapi::mkernel_curvilinear_full_refine(meshKernelId, -3, -4);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_deallocate_state(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}

TEST(CasulliRefinement, CurvilinearRefinementBasedOnDepths)
{
    // Prepare
    int meshKernelId = -1;
    const int projectionType = 0;

    auto errorCode = meshkernelapi::mkernel_allocate_state(projectionType, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::MakeGridParameters makeGridParameters;

    const int numberOfNodes = 11;
    const double delta = 30.0;

    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = delta;
    makeGridParameters.block_size_y = delta;
    makeGridParameters.num_columns = numberOfNodes;
    makeGridParameters.num_rows = numberOfNodes;

    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_convert_to_mesh2d(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    //--------------------------------
    // Prepare depth interpolator

    const int overSampleSize = 3;
    const int numberOfSamplePoints = overSampleSize * (numberOfNodes + 1) * overSampleSize * (numberOfNodes + 1);

    std::vector<double> sampleXNodes(static_cast<size_t>(numberOfSamplePoints));
    std::vector<double> sampleYNodes(static_cast<size_t>(numberOfSamplePoints));
    std::vector<double> sampleDataValues(static_cast<size_t>(numberOfSamplePoints));

    meshkernel::InterpolationParameters interpolationParameters{.interpolation_type = 0};
    int propertyId = -1;
    meshkernelapi::GeometryList sampleData;

    sampleData.num_coordinates = numberOfSamplePoints;
    sampleData.coordinates_x = sampleXNodes.data();
    sampleData.coordinates_y = sampleYNodes.data();
    sampleData.values = sampleDataValues.data();

    size_t sampleCount = 0;

    double sampleDelta = delta / static_cast<double>(overSampleSize);
    double yCoord = -0.5 * sampleDelta;

    for (size_t i = 0; i < static_cast<size_t>(overSampleSize * (numberOfNodes + 1)); ++i)
    {
        double xCoord = -0.5 * sampleDelta;

        for (size_t j = 0; j < static_cast<size_t>(overSampleSize * (numberOfNodes + 1)); ++j)
        {
            sampleXNodes[sampleCount] = xCoord;
            sampleYNodes[sampleCount] = yCoord;
            // Should create a zero line double the middle of the sample data set
            sampleDataValues[sampleCount] = 0.5 * (xCoord - 0.5 * static_cast<double>(numberOfNodes - 1) * delta) * 3.0;
            ++sampleCount;
            xCoord += sampleDelta;
        }

        yCoord += sampleDelta;
    }

    errorCode = meshkernelapi::mkernel_mesh2d_set_property(meshKernelId, interpolationParameters, sampleData, propertyId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    //--------------------------------

    meshkernelapi::GeometryList polygon;
    meshkernelapi::Mesh2D mesh2d{};
    meshkernel::MeshRefinementParameters refinementParameters;
    refinementParameters.min_edge_size = 0.5 * delta; // 12.5;
    refinementParameters.max_courant_time = 2.0;
    const double minimumRefinementDepth = -2.0 * delta;

    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Check number of nodes and edges is correct before refinement
    EXPECT_EQ(144, mesh2d.num_valid_nodes);
    EXPECT_EQ(264, mesh2d.num_valid_edges);

    // Refine using Casulli algorithm
    errorCode = meshkernelapi::mkernel_mesh2d_casulli_refinement_based_on_depths(meshKernelId, polygon, propertyId, refinementParameters, minimumRefinementDepth);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Check number of nodes and edges is correct.
    EXPECT_EQ(252, mesh2d.num_valid_nodes);
    EXPECT_EQ(489, mesh2d.num_valid_edges);

    errorCode = meshkernelapi::mkernel_expunge_state(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}
