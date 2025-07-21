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

#include "MeshKernel/Parameters.hpp"
#include "MeshKernelApi/BoundingBox.hpp"
#include "MeshKernelApi/Mesh2D.hpp"
#include "MeshKernelApi/MeshKernel.hpp"
#include "TestUtils/MakeMeshes.hpp"

#include "TestUtils/Definitions.hpp"

#include "TestMeshGeneration.hpp"

// namespace aliases
namespace mk = meshkernel;
namespace mkapi = meshkernelapi;

TEST(Mesh2DTests, Mesh2dApiNodeEdgeDataTest)
{
    const int clgSize = 2;

    int meshKernelId = meshkernel::constants::missing::intValue;
    int errorCode;

    // Clear the meshkernel state and undo stack before starting the test.
    errorCode = mkapi::mkernel_clear_state();
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    errorCode = mkapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    meshkernel::MakeGridParameters makeGridParameters;

    // num_columns and num_rows indicate number of elements in each direction, so value = clgSize - 1
    makeGridParameters.num_columns = clgSize - 1;
    makeGridParameters.num_rows = clgSize - 1;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = 1.0;
    makeGridParameters.block_size_y = 1.0;

    // Generate curvilinear grid.
    errorCode = mkapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    // Convert the curvilinear grid to an unstructured grid.
    errorCode = mkapi::mkernel_curvilinear_convert_to_mesh2d(meshKernelId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    int node1 = meshkernel::constants::missing::intValue;
    int node2 = meshkernel::constants::missing::intValue;
    int node3 = meshkernel::constants::missing::intValue;
    int node4 = meshkernel::constants::missing::intValue;

    mkapi::BoundingBox boundingBox{-0.1, -0.1, 1.1, 1.1};

    errorCode = mkapi::mkernel_mesh2d_get_location_index(meshKernelId, 0.0, 0.0,
                                                         1 /*Location::Node*/,
                                                         boundingBox,
                                                         node1);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    errorCode = mkapi::mkernel_mesh2d_get_location_index(meshKernelId, 1.0, 0.0,
                                                         1 /*Location::Node*/,
                                                         boundingBox,
                                                         node2);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    errorCode = mkapi::mkernel_mesh2d_get_location_index(meshKernelId, 1.0, 1.0,
                                                         1 /*Location::Node*/,
                                                         boundingBox,
                                                         node3);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    errorCode = mkapi::mkernel_mesh2d_get_location_index(meshKernelId, 0.0, 1.0,
                                                         1 /*Location::Node*/,
                                                         boundingBox,
                                                         node4);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    int newNodeId;
    errorCode = mkapi::mkernel_mesh2d_insert_node(meshKernelId, 0.5, 0.5, newNodeId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    int edge1 = meshkernel::constants::missing::intValue;
    int edge2 = meshkernel::constants::missing::intValue;
    int edge3 = meshkernel::constants::missing::intValue;
    int edge4 = meshkernel::constants::missing::intValue;

    errorCode = mkapi::mkernel_mesh2d_insert_edge(meshKernelId, node1, newNodeId, edge1);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    errorCode = mkapi::mkernel_mesh2d_insert_edge(meshKernelId, node2, newNodeId, edge2);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    errorCode = mkapi::mkernel_mesh2d_insert_edge(meshKernelId, node3, newNodeId, edge3);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    errorCode = mkapi::mkernel_mesh2d_insert_edge(meshKernelId, node4, newNodeId, edge4);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    // All of the newly added edges should be removed.
    // But, there should still be space for them in the mesh node and edge lists.
    errorCode = mkapi::mkernel_mesh2d_delete_node(meshKernelId, newNodeId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    mkapi::Mesh2D mesh2d{};
    errorCode = mkapi::mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    std::vector<double> nodeX(mesh2d.num_nodes);
    std::vector<double> nodeY(mesh2d.num_nodes);
    std::vector<int> edges(2 * mesh2d.num_edges);

    mesh2d.node_x = nodeX.data();
    mesh2d.node_y = nodeY.data();
    mesh2d.edge_nodes = edges.data();

    errorCode = mkapi::mkernel_mesh2d_get_node_edge_data(meshKernelId, mesh2d);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    std::vector<double> expectedX{0.0, 1.0, 0.0, 1.0, meshkernel::constants::missing::doubleValue};
    std::vector<double> expectedY{0.0, 0.0, 1.0, 1.0, meshkernel::constants::missing::doubleValue};
    std::vector<int> expectedEdges{0, 2, 1, 3, 0, 1, 2, 3, -1, -1, -1, -1, -1, -1, -1, -1};

    ASSERT_EQ(mesh2d.num_nodes, clgSize * clgSize + 1);
    ASSERT_EQ(mesh2d.num_edges, 2 * (clgSize - 1) * clgSize + 4);

    for (size_t i = 0; i < nodeX.size(); ++i)
    {
        EXPECT_EQ(nodeX[i], expectedX[i]);
        EXPECT_EQ(nodeY[i], expectedY[i]);
    }

    for (size_t i = 0; i < edges.size(); ++i)
    {
        EXPECT_EQ(edges[i], expectedEdges[i]);
    }
}

TEST(Mesh2DTests, Mesh2DGetPropertyTest)
{
    std::vector<double> nodesX{57.0, 49.1, 58.9, 66.7, 48.8, 65.9, 67.0, 49.1};
    std::vector<double> nodesY{23.6, 14.0, 6.9, 16.2, 23.4, 24.0, 7.2, 6.7};

    std::vector edges{
        0, 1,
        1, 2,
        2, 3,
        0, 3,
        1, 4,
        0, 4,
        0, 5,
        3, 5,
        3, 6,
        2, 6,
        2, 7,
        1, 7};

    meshkernelapi::Mesh2D mesh2d;
    mesh2d.edge_nodes = edges.data();
    mesh2d.node_x = nodesX.data();
    mesh2d.node_y = nodesY.data();
    mesh2d.num_nodes = static_cast<int>(nodesX.size());
    mesh2d.num_edges = static_cast<int>(edges.size() * 0.5);

    int meshKernelId = -1;
    auto errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    int orthogonalityId = -1;
    errorCode = meshkernelapi::mkernel_mesh2d_get_orthogonality_property_type(orthogonalityId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    int geometryListDimension = -1;
    errorCode = meshkernelapi::mkernel_mesh2d_get_property_dimension(meshKernelId, orthogonalityId, geometryListDimension);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    int locationId = static_cast<int>(meshkernel::Location::Edges);
    meshkernelapi::GeometryList propertyvalues{};
    propertyvalues.num_coordinates = geometryListDimension;
    propertyvalues.geometry_separator = meshkernel::constants::missing::doubleValue;
    std::vector<double> values(geometryListDimension);
    propertyvalues.values = values.data();
    errorCode = mkernel_mesh2d_get_property(meshKernelId, orthogonalityId, locationId, propertyvalues);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    EXPECT_EQ(propertyvalues.num_coordinates, 12);
    const double tolerance = 1e-4;

    EXPECT_NEAR(values[0], 0.053115002097392186, tolerance);
    EXPECT_NEAR(values[1], 0.053447414050301602, tolerance);
    EXPECT_NEAR(values[2], 0.053830968680574229, tolerance);
    EXPECT_NEAR(values[3], 0.059369643517766427, tolerance);
}

TEST(Mesh2DTests, Mesh2DGetCircumcenterPropertyTest)
{
    std::vector<double> nodesX{57.0, 49.1, 58.9, 66.7, 48.8, 65.9, 67.0, 49.1};
    std::vector<double> nodesY{23.6, 14.0, 6.9, 16.2, 23.4, 24.0, 7.2, 6.7};

    std::vector edges{
        0, 1,
        1, 2,
        2, 3,
        0, 3,
        1, 4,
        0, 4,
        0, 5,
        3, 5,
        3, 6,
        2, 6,
        2, 7,
        1, 7};

    meshkernelapi::Mesh2D mesh2d;
    mesh2d.edge_nodes = edges.data();
    mesh2d.node_x = nodesX.data();
    mesh2d.node_y = nodesY.data();
    mesh2d.num_nodes = static_cast<int>(nodesX.size());
    mesh2d.num_edges = static_cast<int>(edges.size() * 0.5);

    int meshKernelId = -1;
    auto errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    int circumcenterId = -1;
    errorCode = meshkernelapi::mkernel_mesh2d_get_face_circumcenter_property_type(circumcenterId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    int geometryListDimension = -1;
    errorCode = meshkernelapi::mkernel_mesh2d_get_property_dimension(meshKernelId, circumcenterId, geometryListDimension);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    int locationId = static_cast<int>(meshkernel::Location::Faces);
    meshkernelapi::GeometryList propertyvalues{};
    propertyvalues.num_coordinates = geometryListDimension;
    propertyvalues.geometry_separator = meshkernel::constants::missing::doubleValue;
    std::vector<double> xCoords(geometryListDimension);
    std::vector<double> yCoords(geometryListDimension);
    propertyvalues.coordinates_x = xCoords.data();
    propertyvalues.coordinates_y = yCoords.data();
    errorCode = mkernel_mesh2d_get_property(meshKernelId, circumcenterId, locationId, propertyvalues);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    EXPECT_EQ(propertyvalues.num_coordinates, 5);
    const double tolerance = 1e-4;

    EXPECT_NEAR(xCoords[0], 61.8801441733, tolerance);
    EXPECT_NEAR(xCoords[1], 53.0139097744, tolerance);
    EXPECT_NEAR(xCoords[2], 53.9275510204, tolerance);
    EXPECT_NEAR(xCoords[3], 62.7984959122, tolerance);
    EXPECT_NEAR(xCoords[4], 57.9082712446, tolerance);

    EXPECT_NEAR(yCoords[0], 19.8770034142, tolerance);
    EXPECT_NEAR(yCoords[1], 18.8296992481, tolerance);
    EXPECT_NEAR(yCoords[2], 10.35, tolerance);
    EXPECT_NEAR(yCoords[3], 11.5482066645, tolerance);
    EXPECT_NEAR(yCoords[4], 15.2402740785, tolerance);
}

TEST(Mesh2DTests, Mesh2DGetEdgeLengthPropertyTest)
{
    std::vector<double> nodesX{57.0, 49.1, 58.9, 66.7, 48.8, 65.9, 67.0, 49.1};
    std::vector<double> nodesY{23.6, 14.0, 6.9, 16.2, 23.4, 24.0, 7.2, 6.7};

    std::vector edges{
        0, 1,
        1, 2,
        2, 3,
        0, 3,
        1, 4,
        0, 4,
        0, 5,
        3, 5,
        3, 6,
        2, 6,
        2, 7,
        1, 7};

    meshkernelapi::Mesh2D mesh2d;
    mesh2d.edge_nodes = edges.data();
    mesh2d.node_x = nodesX.data();
    mesh2d.node_y = nodesY.data();
    mesh2d.num_nodes = static_cast<int>(nodesX.size());
    mesh2d.num_edges = static_cast<int>(edges.size() * 0.5);

    int meshKernelId = -1;
    auto errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    int circumcenterId = -1;
    errorCode = meshkernelapi::mkernel_mesh2d_get_edge_length_property_type(circumcenterId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    int geometryListDimension = -1;
    errorCode = meshkernelapi::mkernel_mesh2d_get_property_dimension(meshKernelId, circumcenterId, geometryListDimension);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    int locationId = static_cast<int>(meshkernel::Location::Edges);
    meshkernelapi::GeometryList propertyvalues{};
    propertyvalues.num_coordinates = geometryListDimension;
    propertyvalues.geometry_separator = meshkernel::constants::missing::doubleValue;
    std::vector<double> edgeLengths(geometryListDimension);
    propertyvalues.values = edgeLengths.data();
    errorCode = mkernel_mesh2d_get_property(meshKernelId, circumcenterId, locationId, propertyvalues);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    EXPECT_EQ(propertyvalues.num_coordinates, 12);
    const double tolerance = 1e-4;

    EXPECT_NEAR(edgeLengths[0], 12.4326183887, tolerance);
    EXPECT_NEAR(edgeLengths[1], 12.1016527797, tolerance);
    EXPECT_NEAR(edgeLengths[2], 12.1379569945, tolerance);
    EXPECT_NEAR(edgeLengths[3], 12.2004098292, tolerance);
    EXPECT_NEAR(edgeLengths[4], 9.40478601564, tolerance);
    EXPECT_NEAR(edgeLengths[5], 8.20243866176, tolerance);
    EXPECT_NEAR(edgeLengths[6], 8.90898422942, tolerance);
    EXPECT_NEAR(edgeLengths[7], 7.84091831357, tolerance);
    EXPECT_NEAR(edgeLengths[8], 9.00499861188, tolerance);
    EXPECT_NEAR(edgeLengths[9], 8.10555365166, tolerance);
    EXPECT_NEAR(edgeLengths[10], 9.80204060387, tolerance);
    EXPECT_NEAR(edgeLengths[11], 7.3, tolerance);
}

TEST(Mesh2DTests, GetPolygonsOfDeletedFaces_WithPolygon_ShouldGetPolygonOfDeletedFaces)
{
    // Prepare: set a mesh with two faces sharing an high orthogonality edge. 2 polygon faces should be return
    std::vector<double> nodesX{57.0, 49.1, 58.9, 66.7, 48.8, 65.9, 67.0, 49.1};
    std::vector<double> nodesY{23.6, 14.0, 6.9, 16.2, 23.4, 24.0, 7.2, 6.7};

    std::vector edges{
        0, 1,
        1, 2,
        2, 3,
        0, 3,
        1, 4,
        0, 4,
        0, 5,
        3, 5,
        3, 6,
        2, 6,
        2, 7,
        1, 7};

    meshkernelapi::Mesh2D mesh2d;
    mesh2d.edge_nodes = edges.data();
    mesh2d.node_x = nodesX.data();
    mesh2d.node_y = nodesY.data();
    mesh2d.num_nodes = static_cast<int>(nodesX.size());
    mesh2d.num_edges = static_cast<int>(edges.size() * 0.5);

    int meshKernelId = -1;
    auto errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    int propertyType = -1;
    errorCode = meshkernelapi::mkernel_mesh2d_get_orthogonality_property_type(propertyType);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    int geometryListDimension = -1;
    errorCode = meshkernelapi::mkernel_mesh2d_get_filtered_face_polygons_dimension(meshKernelId,
                                                                                   propertyType,
                                                                                   0.04,
                                                                                   1.0,
                                                                                   geometryListDimension);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(5, geometryListDimension);

    meshkernelapi::GeometryList facePolygons{};
    facePolygons.num_coordinates = geometryListDimension;
    facePolygons.geometry_separator = meshkernel::constants::missing::doubleValue;
    std::vector<double> xfacePolygons(geometryListDimension);
    std::vector<double> yfacePolygons(geometryListDimension);
    facePolygons.coordinates_x = xfacePolygons.data();
    facePolygons.coordinates_y = yfacePolygons.data();
    errorCode = mkernel_mesh2d_get_filtered_face_polygons(meshKernelId,
                                                          propertyType,
                                                          0.04,
                                                          1.0,
                                                          facePolygons);

    // Assert
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    std::vector<double> expectedFacePolygonsX{57.0, 49.1, 58.9, 66.7, 57.0};
    std::vector<double> expectedFacePolygonsY{23.6, 14.0, 6.9, 16.2, 23.6};
    for (size_t i = 0u; i < xfacePolygons.size(); ++i)
    {
        ASSERT_NEAR(expectedFacePolygonsX[i], xfacePolygons[i], 1e-6);
        ASSERT_NEAR(expectedFacePolygonsY[i], yfacePolygons[i], 1e-6);
    }
}

TEST(Mesh2DTests, GetPolygonsOfDeletedFaces_WithPolygon_FailureTests)
{
    meshkernel::MakeGridParameters makeGridParameters;

    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = 1.0;
    makeGridParameters.block_size_y = 1.0;
    makeGridParameters.num_columns = 20;
    makeGridParameters.num_rows = 20;

    int meshKernelId = 0;
    const int projectionType = 0;
    auto errorCode = meshkernelapi::mkernel_allocate_state(projectionType, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = mkapi::mkernel_curvilinear_convert_to_mesh2d(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    //

    int propertyType = -1;
    errorCode = meshkernelapi::mkernel_mesh2d_get_orthogonality_property_type(propertyType);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    double minValue = -1.0;
    double maxValue = 1.0;

    meshkernelapi::GeometryList facePolygons{};

    // face data not yet cached
    errorCode = mkernel_mesh2d_get_filtered_face_polygons(meshKernelId,
                                                          propertyType + 1,
                                                          minValue,
                                                          maxValue,
                                                          facePolygons);

    // Assert
    ASSERT_EQ(meshkernel::ExitCode::ConstraintErrorCode, errorCode);

    // Execute
    int geometryListDimension = -1;
    errorCode = meshkernelapi::mkernel_mesh2d_get_filtered_face_polygons_dimension(meshKernelId,
                                                                                   propertyType,
                                                                                   minValue,
                                                                                   maxValue,
                                                                                   geometryListDimension);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    facePolygons.num_coordinates = geometryListDimension;
    facePolygons.geometry_separator = meshkernel::constants::missing::doubleValue;
    std::vector<double> xfacePolygons(geometryListDimension);
    std::vector<double> yfacePolygons(geometryListDimension);
    facePolygons.coordinates_x = xfacePolygons.data();
    facePolygons.coordinates_y = yfacePolygons.data();

    // Check different property type
    errorCode = mkernel_mesh2d_get_filtered_face_polygons(meshKernelId,
                                                          propertyType + 1,
                                                          minValue,
                                                          maxValue,
                                                          facePolygons);

    // Assert
    ASSERT_EQ(meshkernel::ExitCode::ConstraintErrorCode, errorCode);

    // Check different minimum value
    errorCode = mkernel_mesh2d_get_filtered_face_polygons(meshKernelId,
                                                          propertyType,
                                                          minValue + 0.5,
                                                          maxValue,
                                                          facePolygons);

    // Assert
    ASSERT_EQ(meshkernel::ExitCode::ConstraintErrorCode, errorCode);
    // Check different minimum value
    errorCode = mkernel_mesh2d_get_filtered_face_polygons(meshKernelId,
                                                          propertyType,
                                                          minValue,
                                                          maxValue + 0.5,
                                                          facePolygons);

    // Assert
    ASSERT_EQ(meshkernel::ExitCode::ConstraintErrorCode, errorCode);
}

TEST(Mesh2D, Mesh2DInitializeOrthogonalization_WithHexagon_ShouldOrthogonalize)
{
    // Prepare
    int meshKernelId;
    const int isGeographic = 0;
    meshkernelapi::mkernel_allocate_state(isGeographic, meshKernelId);

    auto [num_nodes, num_edges, node_x, node_y, node_type, edge_nodes, edge_type] =
        ReadLegacyMeshFile(TEST_FOLDER + "/data/MeshWithHexagon.nc");
    meshkernelapi::Mesh2D mesh2d;
    mesh2d.num_edges = static_cast<int>(num_edges);
    mesh2d.num_nodes = static_cast<int>(num_nodes);
    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_nodes = edge_nodes.data();

    auto errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::OrthogonalizationParameters orthogonalizationParameters{};
    orthogonalizationParameters.outer_iterations = 1;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;

    meshkernelapi::GeometryList geometryList{};
    meshkernelapi::GeometryList landBoundaries{};

    // Execute
    errorCode = mkernel_mesh2d_initialize_orthogonalization(meshKernelId,
                                                            1,
                                                            orthogonalizationParameters,
                                                            landBoundaries,
                                                            geometryList);

    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_prepare_outer_iteration_orthogonalization(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_compute_inner_ortogonalization_iteration(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_finalize_inner_ortogonalization_iteration(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_delete_orthogonalization(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}

TEST(Mesh2D, IntersectionsFromPolyline_ShouldIntersectMesh)
{
    // Setup
    int meshKernelId;
    auto errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);

    // Create a curvilinear grid in the back-end and convert to an unstructured grid
    meshkernel::MakeGridParameters makeMeshParameters;
    makeMeshParameters.num_columns = 3;
    makeMeshParameters.num_rows = 3;
    makeMeshParameters.block_size_x = 1.0;
    makeMeshParameters.block_size_y = 1.0;
    makeMeshParameters.origin_x = 0.0;
    makeMeshParameters.origin_y = 0.0;
    makeMeshParameters.angle = 0.0;

    // Creates an unstructured grid from mesh parameters
    errorCode = meshkernelapi::mkernel_mesh2d_make_rectangular_mesh(meshKernelId, makeMeshParameters);

    // Get the mesh dimensions
    meshkernelapi::Mesh2D mesh2dDimensions{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2dDimensions);

    // Set the polyLine
    std::vector xCoordinates{0.6, 0.6, 2.5, 2.5, 0.6};
    std::vector yCoordinates{2.5, 0.5, 0.5, 2.5, 2.5};

    meshkernelapi::GeometryList boundaryPolygon{};
    boundaryPolygon.geometry_separator = meshkernel::constants::missing::doubleValue;
    boundaryPolygon.coordinates_x = xCoordinates.data();
    boundaryPolygon.coordinates_y = yCoordinates.data();
    boundaryPolygon.values = nullptr;

    boundaryPolygon.num_coordinates = static_cast<int>(xCoordinates.size());

    std::vector<int> polylineSegmentIndexes(mesh2dDimensions.num_edges * 2, meshkernel::constants::missing::intValue);

    std::vector<int> edgeNodes(mesh2dDimensions.num_edges * 2, meshkernel::constants::missing::intValue);
    std::vector<int> edgeIndex(mesh2dDimensions.num_edges, meshkernel::constants::missing::intValue);
    std::vector<double> edgeDistances(mesh2dDimensions.num_edges * 2, meshkernel::constants::missing::doubleValue);
    std::vector<double> segmentDistances(mesh2dDimensions.num_edges * 2, meshkernel::constants::missing::doubleValue);
    std::vector<int> segmentIndexes(mesh2dDimensions.num_edges, meshkernel::constants::missing::intValue);

    std::vector<int> faceEdgeIndex(mesh2dDimensions.num_edges, meshkernel::constants::missing::intValue);
    std::vector<int> faceNumEdges(mesh2dDimensions.num_edges, meshkernel::constants::missing::intValue);
    std::vector<int> faceIndexes(mesh2dDimensions.num_edges, meshkernel::constants::missing::intValue);

    errorCode = mkernel_mesh2d_intersections_from_polygon(meshKernelId,
                                                          boundaryPolygon,
                                                          edgeNodes.data(),
                                                          edgeIndex.data(),
                                                          edgeDistances.data(),
                                                          segmentDistances.data(),
                                                          segmentIndexes.data(),
                                                          faceIndexes.data(),
                                                          faceNumEdges.data(),
                                                          faceEdgeIndex.data());

    /// Assert
    const double tolerance = 1e-6;

    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    ASSERT_EQ(segmentIndexes[0], 0);
    ASSERT_EQ(segmentIndexes[1], 0);
    ASSERT_EQ(segmentIndexes[2], 1);

    ASSERT_NEAR(segmentDistances[0], 0.25, tolerance);
    ASSERT_NEAR(segmentDistances[1], 0.75, tolerance);
    ASSERT_NEAR(segmentDistances[2], 0.21052631578947370, tolerance);

    ASSERT_NEAR(edgeDistances[0], 0.6, tolerance);
    ASSERT_NEAR(edgeDistances[1], 0.6, tolerance);
    ASSERT_NEAR(edgeDistances[2], 0.50000000000000000, tolerance);

    ASSERT_EQ(faceIndexes[0], 3);
    ASSERT_EQ(faceIndexes[1], 3);
    ASSERT_EQ(faceIndexes[2], 0);
    ASSERT_EQ(faceIndexes[3], 0);
    ASSERT_EQ(faceIndexes[4], 1);

    ASSERT_EQ(faceNumEdges[0], 2);
    ASSERT_EQ(faceNumEdges[1], 2);
    ASSERT_EQ(faceNumEdges[2], 2);
    ASSERT_EQ(faceNumEdges[3], 2);

    ASSERT_EQ(faceEdgeIndex[0], 18);
    ASSERT_EQ(faceEdgeIndex[1], 15);
    ASSERT_EQ(faceEdgeIndex[2], 1);
    ASSERT_EQ(faceEdgeIndex[3], 15);
    ASSERT_EQ(faceEdgeIndex[4], 1);
}
TEST(Mesh2D, CurvilinearMakeRectangularOnExtension_OnSpericalCoordinates_ShouldGenerateCurvilinearMesh)
{
    // Setup
    auto makeGridParameters = meshkernel::MakeGridParameters();
    makeGridParameters.origin_x = -1.0;
    makeGridParameters.origin_y = 49.1;
    makeGridParameters.block_size_x = 0.01;
    makeGridParameters.block_size_y = 0.01;
    makeGridParameters.upper_right_x = -0.2;
    makeGridParameters.upper_right_y = 49.6;

    // Execute
    int meshKernelId;
    const int projectionType = 1;
    auto errorCode = meshkernelapi::mkernel_allocate_state(projectionType, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_make_rectangular_mesh_on_extension(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d;
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    EXPECT_EQ(mesh2d.num_nodes, 6318);
    EXPECT_EQ(mesh2d.num_edges, 12477);
    EXPECT_EQ(mesh2d.num_faces, 6160);
}

TEST(Mesh2D, RemoveSingleIsland)
{
    // Load mesh with 2 disconnected regions, first a 10x10 and the second is a smaller 2x2 mesh
    auto [num_nodes, num_edges, node_x, node_y, node_type, edge_nodes, edge_type] =
        ReadLegacyMeshFile(TEST_FOLDER + "/data/RemoveDomainIslands/single_disconnected_region.nc");

    meshkernelapi::Mesh2D mesh;
    mesh.num_edges = static_cast<int>(num_edges);
    mesh.num_nodes = static_cast<int>(num_nodes);
    mesh.node_x = node_x.data();
    mesh.node_y = node_y.data();
    mesh.edge_nodes = edge_nodes.data();

    const int projectionType = 1;
    int meshKernelId;
    auto errorCode = meshkernelapi::mkernel_allocate_state(projectionType, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_set(meshKernelId, mesh);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_remove_disconnected_regions(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}

TEST(Mesh2D, RemoveMultipleIslands)
{
    // Load mesh with 4 disconnected regions, the main domain is a 10x10, there are 3 other much small island regions,
    // each with a different shape and number of elements.
    auto [num_nodes, num_edges, node_x, node_y, node_type, edge_nodes, edge_type] =
        ReadLegacyMeshFile(TEST_FOLDER + "/data/RemoveDomainIslands/single_disconnected_region.nc");

    meshkernelapi::Mesh2D mesh;
    mesh.num_edges = static_cast<int>(num_edges);
    mesh.num_nodes = static_cast<int>(num_nodes);
    mesh.node_x = node_x.data();
    mesh.node_y = node_y.data();
    mesh.edge_nodes = edge_nodes.data();

    const int projectionType = 1;
    int meshKernelId;
    auto errorCode = meshkernelapi::mkernel_allocate_state(projectionType, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_set(meshKernelId, mesh);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_remove_disconnected_regions(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}

TEST(Mesh2D, Mesh2DSetAndAdd)
{
    using namespace meshkernelapi;

    int errorCode = mkernel_clear_undo_state();
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::UInt const num_nodes_x = 20;
    meshkernel::UInt const num_nodes_y = 15;
    double const delta = 1.0;
    double const offset_x = 0.5;

    // create first mesh
    auto [num_nodes_1, num_edges_1, node_x_1, node_y_1, edge_nodes_1] =
        MakeRectangularMeshForApiTesting(num_nodes_x,
                                         num_nodes_y,
                                         delta,
                                         meshkernel::Point(0.0, 0.0));
    Mesh2D mesh2d_1{};
    mesh2d_1.num_nodes = static_cast<int>(num_nodes_1);
    mesh2d_1.num_edges = static_cast<int>(num_edges_1);
    mesh2d_1.node_x = node_x_1.data();
    mesh2d_1.node_y = node_y_1.data();
    mesh2d_1.edge_nodes = edge_nodes_1.data();

    // create second mesh
    auto [num_nodes_2, num_edges_2, node_x_2, node_y_2, edge_nodes_2] =
        MakeRectangularMeshForApiTesting(num_nodes_x,
                                         num_nodes_y,
                                         delta,
                                         meshkernel::Point(0.0 + offset_x, 0.0));
    Mesh2D mesh2d_2{};
    mesh2d_2.num_nodes = static_cast<int>(num_nodes_2);
    mesh2d_2.num_edges = static_cast<int>(num_edges_2);
    mesh2d_2.node_x = node_x_2.data();
    mesh2d_2.node_y = node_y_2.data();
    mesh2d_2.edge_nodes = edge_nodes_2.data();

    // allocate state
    int mk_id = 0;
    errorCode = mkernel_allocate_state(0, mk_id);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // first initialise using the first mesh, mesh2d_1
    errorCode = mkernel_mesh2d_set(mk_id, mesh2d_1);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // then add the second mesh, mesh2d_2
    errorCode = mkernel_mesh2d_add(mk_id, mesh2d_2);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // get the dimensions of the resulting mesh
    Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(mk_id, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // allocate memory for the arrays of mesh2d
    std::vector<int> edge_faces(mesh2d.num_edges * 2);
    std::vector<int> edge_nodes(mesh2d.num_edges * 2);
    std::vector<int> face_nodes(mesh2d.num_face_nodes);
    std::vector<int> face_edges(mesh2d.num_face_nodes);
    std::vector<int> nodes_per_face(mesh2d.num_faces);
    std::vector<double> node_x(mesh2d.num_nodes);
    std::vector<double> node_y(mesh2d.num_nodes);
    std::vector<double> edge_x(mesh2d.num_edges);
    std::vector<double> edge_y(mesh2d.num_edges);
    std::vector<double> face_x(mesh2d.num_faces);
    std::vector<double> face_y(mesh2d.num_faces);

    mesh2d.edge_faces = edge_faces.data();
    mesh2d.edge_nodes = edge_nodes.data();
    mesh2d.face_nodes = face_nodes.data();
    mesh2d.face_edges = face_edges.data();
    mesh2d.nodes_per_face = nodes_per_face.data();
    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_x = edge_x.data();
    mesh2d.edge_y = edge_y.data();
    mesh2d.face_x = face_x.data();
    mesh2d.face_y = face_y.data();

    // upon allocation, get the data of the resulting mesh
    errorCode = mkernel_mesh2d_get_data(mk_id, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    ASSERT_EQ(mesh2d.num_nodes, num_nodes_1 + num_nodes_2);
    ASSERT_EQ(mesh2d.num_edges, num_edges_1 + num_edges_2);

    for (meshkernel::UInt i = 0; i < num_nodes_1; ++i)
    {
        EXPECT_EQ(mesh2d.node_x[i], mesh2d_1.node_x[i]);
        EXPECT_EQ(mesh2d.node_y[i], mesh2d_1.node_y[i]);
    }

    for (meshkernel::UInt i = num_nodes_1; i < num_nodes_1 + num_nodes_2; ++i)
    {
        EXPECT_EQ(mesh2d.node_x[i], mesh2d_2.node_x[i - num_nodes_1]);
        EXPECT_EQ(mesh2d.node_y[i], mesh2d_2.node_y[i - num_nodes_1]);
    }

    errorCode = mkernel_deallocate_state(mk_id);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}

TEST(Mesh2D, Mesh2DAddEdge)
{
    using namespace meshkernelapi;

    int errorCode = mkernel_clear_state();
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::UInt const num_nodes_x = 4;
    meshkernel::UInt const num_nodes_y = 4;
    double const delta = 1.0;

    // create first mesh
    auto [num_nodes, num_edges, node_x, node_y, edge_nodes] =
        MakeRectangularMeshForApiTesting(num_nodes_x,
                                         num_nodes_y,
                                         delta,
                                         meshkernel::Point(0.0, 0.0));
    Mesh2D mesh2d{};
    mesh2d.num_nodes = static_cast<int>(num_nodes);
    mesh2d.num_edges = static_cast<int>(num_edges);
    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_nodes = edge_nodes.data();

    // allocate state
    int mk_id = 0;
    errorCode = mkernel_allocate_state(0, mk_id);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // first initialise using the first mesh, mesh2d
    errorCode = mkernel_mesh2d_set(mk_id, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    int newEdgeId = -1;
    errorCode = mkernel_mesh2d_insert_edge(mk_id, 0, 4, newEdgeId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_TRUE(newEdgeId > 0);

    // Should be only a single item on the undo action stack
    bool undoInsertEdge = false;
    int undoId = meshkernel::constants::missing::intValue;
    errorCode = mkernel_undo_state(undoInsertEdge, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_TRUE(undoInsertEdge);
    ASSERT_EQ(mk_id, undoId);

    // Undo creation of mesh2d
    undoInsertEdge = false;
    errorCode = mkernel_undo_state(undoInsertEdge, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_TRUE(undoInsertEdge);
    ASSERT_EQ(undoId, mk_id);

    // Should be no items on the undo action stack
    undoInsertEdge = false;
    errorCode = mkernel_undo_state(undoInsertEdge, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_FALSE(undoInsertEdge);
    ASSERT_EQ(undoId, meshkernel::constants::missing::intValue);
}

TEST(Mesh2D, SimpleMultiMeshUndoTest)
{
    using namespace meshkernelapi;

    int committedCount = 0;
    int restoredCount = 0;

    int errorCode = mkernel_clear_undo_state();
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::UInt const num_nodes_x = 4;
    meshkernel::UInt const num_nodes_y = 4;
    double const delta = 1.0;

    // create first mesh
    auto [num_nodes, num_edges, node_x, node_y, edge_nodes] =
        MakeRectangularMeshForApiTesting(num_nodes_x,
                                         num_nodes_y,
                                         delta,
                                         meshkernel::Point(0.0, 0.0));
    Mesh2D mesh2d{};
    mesh2d.num_nodes = static_cast<int>(num_nodes);
    mesh2d.num_edges = static_cast<int>(num_edges);
    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_nodes = edge_nodes.data();

    // allocate state
    int mkid1 = 0;
    errorCode = mkernel_allocate_state(0, mkid1);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // first initialise using the first mesh, mesh2d
    errorCode = mkernel_mesh2d_set(mkid1, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // create second mesh
    std::tie(num_nodes, num_edges, node_x, node_y, edge_nodes) =
        MakeRectangularMeshForApiTesting(num_nodes_x,
                                         num_nodes_y,
                                         delta,
                                         meshkernel::Point(10.0, 10.0));
    mesh2d.num_nodes = static_cast<int>(num_nodes);
    mesh2d.num_edges = static_cast<int>(num_edges);
    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_nodes = edge_nodes.data();

    int mkid2 = 0;
    errorCode = mkernel_allocate_state(0, mkid2);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // initialise using the second mesh, mesh2d
    errorCode = mkernel_mesh2d_set(mkid2, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Start test

    // Check validity of mkid's
    bool isValid = false;
    errorCode = mkernel_is_valid_state(mkid2, isValid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(isValid);

    isValid = false;
    errorCode = mkernel_is_valid_state(mkid1, isValid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(isValid);

    int nodeId = 0;
    errorCode = mkernel_mesh2d_insert_node(mkid1, 0.25, 0.25, nodeId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = mkernel_mesh2d_insert_node(mkid1, 0.5, 0.25, nodeId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = mkernel_mesh2d_insert_node(mkid1, 0.75, 0.25, nodeId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = mkernel_mesh2d_insert_node(mkid1, 0.25, 0.5, nodeId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Save id to delete node later
    int nodeId2;
    errorCode = mkernel_mesh2d_insert_node(mkid1, 0.5, 0.5, nodeId2);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = mkernel_mesh2d_insert_node(mkid2, 0.75, 0.5, nodeId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = mkernel_mesh2d_insert_node(mkid2, 0.25, 0.75, nodeId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = mkernel_mesh2d_insert_node(mkid2, 0.5, 0.75, nodeId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = mkernel_mesh2d_insert_node(mkid2, 0.75, 0.75, nodeId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Delete node
    errorCode = mkernel_mesh2d_delete_node(mkid1, nodeId2);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = mkernel_undo_state_count(committedCount, restoredCount);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    EXPECT_EQ(committedCount, 12);
    EXPECT_EQ(restoredCount, 0);

    bool didUndo = false;
    int undoId = meshkernel::constants::missing::intValue;
    // Undo deletion of node from mkid1
    errorCode = mkernel_undo_state(didUndo, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndo);
    ASSERT_EQ(mkid1, undoId);

    didUndo = false;
    // Undo node insertion from mkid2
    errorCode = mkernel_undo_state(didUndo, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndo);
    ASSERT_EQ(mkid2, undoId);

    errorCode = mkernel_undo_state_count(committedCount, restoredCount);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    EXPECT_EQ(committedCount, 10);
    EXPECT_EQ(restoredCount, 2);

    errorCode = mkernel_expunge_state(mkid1);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = mkernel_undo_state_count_for_id(mkid2, committedCount, restoredCount);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    EXPECT_EQ(committedCount, 4);
    EXPECT_EQ(restoredCount, 1);

    isValid = false;
    errorCode = mkernel_is_valid_state(mkid1, isValid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_FALSE(isValid);

    errorCode = mkernel_expunge_state(mkid2);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    isValid = false;
    errorCode = mkernel_is_valid_state(mkid2, isValid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_FALSE(isValid);

    errorCode = mkernel_undo_state_count_for_id(mkid2, committedCount, restoredCount);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    EXPECT_EQ(committedCount, 0);
    EXPECT_EQ(restoredCount, 0);
}

TEST(Mesh2D, Mesh2DInsertNode)
{
    using namespace meshkernelapi;

    // Clear the undo stack before starting the test.
    int errorCode = mkernel_clear_undo_state();
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::UInt const num_nodes_x = 4;
    meshkernel::UInt const num_nodes_y = 4;
    double const delta = 1.0;

    // create first mesh
    auto [num_nodes, num_edges, node_x, node_y, edge_nodes] =
        MakeRectangularMeshForApiTesting(num_nodes_x,
                                         num_nodes_y,
                                         delta,
                                         meshkernel::Point(0.0, 0.0));
    Mesh2D mesh2d{};
    mesh2d.num_nodes = static_cast<int>(num_nodes);
    mesh2d.num_edges = static_cast<int>(num_edges);
    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_nodes = edge_nodes.data();

    // allocate state
    int mk_id = 0;
    errorCode = mkernel_allocate_state(0, mk_id);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // first initialise using the first mesh, mesh2d
    errorCode = mkernel_mesh2d_set(mk_id, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    //--------------------------------

    int newNodeId = -1;
    errorCode = mkernel_mesh2d_insert_node(mk_id, 0.5, -1.0, newNodeId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_TRUE(newNodeId > 0);

    // Should be only a single item on the undo action stack
    bool undoInsertNode = false;
    int undoId = meshkernel::constants::missing::intValue;
    errorCode = mkernel_undo_state(undoInsertNode, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_TRUE(undoInsertNode);
    ASSERT_EQ(mk_id, undoId);

    // Undo creation of the mesh2d
    undoInsertNode = false;
    errorCode = mkernel_undo_state(undoInsertNode, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_TRUE(undoInsertNode);
    ASSERT_EQ(undoId, mk_id);

    // Should be zero items on the undo action stack.
    undoInsertNode = false;
    errorCode = mkernel_undo_state(undoInsertNode, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_FALSE(undoInsertNode);
    ASSERT_EQ(undoId, meshkernel::constants::missing::intValue);
}

TEST(Mesh2D, InsertEdgeFromCoordinates_OnEmptyMesh_ShouldInsertNewEdge)
{
    // Prepare
    int meshKernelId = 0;
    int errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    int firstNodeIndex;
    int secondNodeIndex;
    int edgeIndex;
    errorCode = meshkernelapi::mkernel_mesh2d_insert_edge_from_coordinates(meshKernelId,
                                                                           0,
                                                                           0,
                                                                           1,
                                                                           0,
                                                                           firstNodeIndex,
                                                                           secondNodeIndex,
                                                                           edgeIndex);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    ASSERT_EQ(0, edgeIndex);
    ASSERT_EQ(0, firstNodeIndex);
    ASSERT_EQ(1, secondNodeIndex);

    ASSERT_EQ(2, mesh2d.num_nodes);
    ASSERT_EQ(1, mesh2d.num_edges);
}
TEST(Mesh2D, ConvertToCurvilinear_ShouldConvertMeshToCurvilinear)
{
    // create first mesh
    auto [num_nodes, num_edges, node_x, node_y, edge_nodes] = MakeRectangularMeshForApiTesting(10, 10, 1.0, meshkernel::Point(0.0, 0.0));

    meshkernelapi::Mesh2D mesh2d{};
    mesh2d.num_nodes = static_cast<int>(num_nodes);
    mesh2d.num_edges = static_cast<int>(num_edges);
    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_nodes = edge_nodes.data();

    // allocate state
    int meshKernelId = 0;
    int errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // first initialise using the first mesh, mesh2d
    errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    errorCode = meshkernelapi::mkernel_mesh2d_convert_to_curvilinear(meshKernelId, 5.0, 5.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2dOut{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2dOut);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearOut{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearOut);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    ASSERT_EQ(0, mesh2dOut.num_valid_nodes);
    ASSERT_EQ(0, mesh2dOut.num_valid_edges);

    ASSERT_EQ(11, curvilinearOut.num_m);
    ASSERT_EQ(11, curvilinearOut.num_n);
}

TEST(Mesh2D, ConvertToCurvilinear_OnStructuredAndUnstructuredMesh_ShouldConvertMeshToCurvilinear)
{
    // Prepare
    int meshKernelId;
    const int isGeographic = 0;
    meshkernelapi::mkernel_allocate_state(isGeographic, meshKernelId);

    auto [num_nodes, num_edges, node_x, node_y, node_type, edge_nodes, edge_type] =
        ReadLegacyMeshFile(TEST_FOLDER + "/data/SmallCurviAndTriangularMesh.nc");
    meshkernelapi::Mesh2D mesh2d;
    mesh2d.num_edges = static_cast<int>(num_edges);
    mesh2d.num_nodes = static_cast<int>(num_nodes);
    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_nodes = edge_nodes.data();

    auto errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2dOut{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2dOut);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(152, mesh2dOut.num_valid_nodes);
    ASSERT_EQ(308, mesh2dOut.num_valid_edges);

    // Execute: convert mesh to curvilinear
    errorCode = meshkernelapi::mkernel_mesh2d_convert_to_curvilinear(meshKernelId, 5.0, 5.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2dOut);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearOut{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearOut);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert: mesh and curvilinear grid co-exist in the same instance
    ASSERT_EQ(38, mesh2dOut.num_valid_nodes);
    ASSERT_EQ(94, mesh2dOut.num_valid_edges);
    ASSERT_EQ(11, curvilinearOut.num_m);
    ASSERT_EQ(11, curvilinearOut.num_n);

    // Execute: curvilinear to mesh
    errorCode = meshkernelapi::mkernel_curvilinear_convert_to_mesh2d(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2dOut);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(159, mesh2dOut.num_valid_nodes);
    ASSERT_EQ(314, mesh2dOut.num_valid_edges);

    // Execute: Merge nodes ro restore mesh back to the original state
    meshkernelapi::GeometryList geometry_list;
    geometry_list.num_coordinates = 0;
    errorCode = mkernel_mesh2d_merge_nodes_with_merging_distance(meshKernelId, geometry_list, 0.001);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2dOut);
    ASSERT_EQ(152, mesh2dOut.num_valid_nodes);
    ASSERT_EQ(308, mesh2dOut.num_valid_edges);
}

TEST(Mesh2d, GetFacePolygons_OnAValidMesh_ShouldGetFacePolygons)
{
    int meshKernelId;
    int isGeographic = 0;
    meshkernelapi::mkernel_allocate_state(isGeographic, meshKernelId);

    meshkernelapi::Mesh2D mesh2d;
    std::vector<double> node_x{
        0,
        1,
        2,
        3,
        0,
        1,
        2,
        3,
        0,
        1,
        3,
        0,
        1};

    std::vector<double> node_y{
        0,
        0,
        0,
        0,
        1,
        1,
        1,
        1,
        2,
        2,
        2,
        3,
        3};

    std::vector edge_nodes{
        0, 4,
        1, 5,
        2, 6,
        3, 7,
        4, 8,
        5, 9,
        7, 10,
        8, 11,
        10, 12,
        0, 1,
        1, 2,
        2, 3,
        5, 6,
        6, 7,
        8, 9,
        9, 10,
        11, 12};

    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_nodes = edge_nodes.data();
    mesh2d.num_edges = static_cast<int>(edge_nodes.size() * 0.5);
    mesh2d.num_nodes = static_cast<int>(node_x.size());

    auto errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    int geometryListDimension = -1;
    errorCode = meshkernelapi::mkernel_mesh2d_get_face_polygons_dimension(meshKernelId, 5, geometryListDimension);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(13, geometryListDimension);

    meshkernelapi::GeometryList geometryList;
    auto facesXCoordinates = std::vector<double>(geometryListDimension);
    auto facesYCoordinates = std::vector<double>(geometryListDimension);
    geometryList.coordinates_x = facesXCoordinates.data();
    geometryList.coordinates_y = facesYCoordinates.data();
    geometryList.num_coordinates = geometryListDimension;

    errorCode = mkernel_mesh2d_get_face_polygons(meshKernelId, 5, geometryList);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    const auto expectedFacesXCoordinates = std::vector{1.0000,
                                                       2.0000,
                                                       3.0000,
                                                       3.0000,
                                                       1.0000,
                                                       1.0000,
                                                       -999.0,
                                                       0.0000,
                                                       1.0000,
                                                       3.0000,
                                                       1.0000,
                                                       0.0000,
                                                       0.0000};

    const auto expectedFacesYCoordinates = std::vector{1.0000,
                                                       1.0000,
                                                       1.0000,
                                                       2.0000,
                                                       2.0000,
                                                       1.0000,
                                                       -999.0,
                                                       2.0000,
                                                       2.0000,
                                                       2.0000,
                                                       3.0000,
                                                       3.0000,
                                                       2.0000};

    ASSERT_THAT(expectedFacesXCoordinates, ::testing::ContainerEq(facesXCoordinates));
    ASSERT_THAT(expectedFacesYCoordinates, ::testing::ContainerEq(facesYCoordinates));
}

TEST(Mesh2D, UndoConnectMeshes)
{
    const meshkernel::UInt num_nodes_x = 3;
    const meshkernel::UInt num_nodes_y = 3;
    const double delta = 1.0;

    // create first mesh
    auto [num_nodes_1, num_edges_1, node_x_1, node_y_1, edge_nodes_1] =
        MakeRectangularMeshForApiTesting(num_nodes_x,
                                         num_nodes_y,
                                         delta,
                                         meshkernel::Point(0.0, 0.0));
    meshkernelapi::Mesh2D mesh2d_1{};
    mesh2d_1.num_nodes = static_cast<int>(num_nodes_1);
    mesh2d_1.num_edges = static_cast<int>(num_edges_1);
    mesh2d_1.node_x = node_x_1.data();
    mesh2d_1.node_y = node_y_1.data();
    mesh2d_1.edge_nodes = edge_nodes_1.data();

    // create second mesh
    auto [num_nodes_2, num_edges_2, node_x_2, node_y_2, edge_nodes_2] =
        MakeRectangularMeshForApiTesting(2 * num_nodes_x,
                                         2 * num_nodes_y,
                                         0.5 * delta,
                                         meshkernel::Point(3.0, 0.0));
    meshkernelapi::Mesh2D mesh2d_2{};
    mesh2d_2.num_nodes = static_cast<int>(num_nodes_2);
    mesh2d_2.num_edges = static_cast<int>(num_edges_2);
    mesh2d_2.node_x = node_x_2.data();
    mesh2d_2.node_y = node_y_2.data();
    mesh2d_2.edge_nodes = edge_nodes_2.data();

    // allocate state
    int mk_id = 0;
    int errorCode = meshkernelapi::mkernel_allocate_state(0, mk_id);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // first initialise using the first mesh, mesh2d_1
    errorCode = meshkernelapi::mkernel_mesh2d_set(mk_id, mesh2d_1);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // The polygon defining the region to be connected.
    meshkernelapi::GeometryList selectedRegion;

    errorCode = meshkernelapi::mkernel_mesh2d_connect_meshes(mk_id, mesh2d_2, selectedRegion, 0.1, true);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    bool didUndo = false;
    int undoMkId = meshkernel::constants::missing::intValue;

    errorCode = meshkernelapi::mkernel_undo_state(didUndo, undoMkId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(undoMkId, mk_id);

    errorCode = meshkernelapi::mkernel_redo_state(didUndo, undoMkId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(undoMkId, mk_id);
}

TEST(Mesh2D, SimpleConnectMeshes)
{
    const meshkernel::UInt num_nodes_x = 3;
    const meshkernel::UInt num_nodes_y = 3;
    const double delta = 1.0;

    // create first mesh
    auto [num_nodes_1, num_edges_1, node_x_1, node_y_1, edge_nodes_1] =
        MakeRectangularMeshForApiTesting(num_nodes_x,
                                         num_nodes_y,
                                         delta,
                                         meshkernel::Point(0.0, 0.0));
    meshkernelapi::Mesh2D mesh2d_1{};
    mesh2d_1.num_nodes = static_cast<int>(num_nodes_1);
    mesh2d_1.num_edges = static_cast<int>(num_edges_1);
    mesh2d_1.node_x = node_x_1.data();
    mesh2d_1.node_y = node_y_1.data();
    mesh2d_1.edge_nodes = edge_nodes_1.data();

    // create second mesh
    auto [num_nodes_2, num_edges_2, node_x_2, node_y_2, edge_nodes_2] =
        MakeRectangularMeshForApiTesting(num_nodes_x,
                                         2 * num_nodes_y,
                                         0.5 * delta,
                                         meshkernel::Point(3.0, 0.0));
    meshkernelapi::Mesh2D mesh2d_2{};
    mesh2d_2.num_nodes = static_cast<int>(num_nodes_2);
    mesh2d_2.num_edges = static_cast<int>(num_edges_2);
    mesh2d_2.node_x = node_x_2.data();
    mesh2d_2.node_y = node_y_2.data();
    mesh2d_2.edge_nodes = edge_nodes_2.data();

    // allocate state
    int mk_id = 0;
    int errorCode = meshkernelapi::mkernel_allocate_state(0, mk_id);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // first initialise using the first mesh, mesh2d_1
    errorCode = meshkernelapi::mkernel_mesh2d_set(mk_id, mesh2d_1);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    // The polygon defining the region to be connected.
    meshkernelapi::GeometryList selectedRegion;

    errorCode = meshkernelapi::mkernel_mesh2d_connect_meshes(mk_id, mesh2d_2, selectedRegion, 0.1, false);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d_3{};

    errorCode = meshkernelapi::mkernel_mesh2d_get_dimensions(mk_id, mesh2d_3);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    ASSERT_EQ(mesh2d_3.num_nodes, mesh2d_1.num_nodes + mesh2d_2.num_nodes);
    ASSERT_EQ(mesh2d_3.num_edges, mesh2d_1.num_edges + mesh2d_2.num_edges);

    std::vector<double> node_x_3(mesh2d_3.num_nodes);
    std::vector<double> node_y_3(mesh2d_3.num_nodes);

    std::vector<int> edge_nodes_3(2 * mesh2d_3.num_edges);

    mesh2d_3.node_x = node_x_3.data();
    mesh2d_3.node_y = node_y_3.data();
    mesh2d_3.edge_nodes = edge_nodes_3.data();

    errorCode = meshkernelapi::mkernel_mesh2d_get_node_edge_data(mk_id, mesh2d_3);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    for (size_t i = 0; i < node_x_1.size(); ++i)
    {
        EXPECT_EQ(node_x_1[i], node_x_3[i]);
        EXPECT_EQ(node_y_1[i], node_y_3[i]);
    }

    for (size_t i = 0; i < node_x_2.size(); ++i)
    {
        EXPECT_EQ(node_x_2[i], node_x_3[i + node_x_1.size()]);
        EXPECT_EQ(node_y_2[i], node_y_3[i + node_x_1.size()]);
    }

    for (size_t i = 0; i < edge_nodes_1.size(); ++i)
    {
        EXPECT_EQ(edge_nodes_1[i], edge_nodes_3[i]);
    }

    for (size_t i = 0; i < edge_nodes_2.size(); ++i)
    {
        EXPECT_EQ(edge_nodes_2[i] + mesh2d_1.num_nodes, edge_nodes_3[i + edge_nodes_1.size()]);
    }
}

TEST(Mesh2D, InsertEdgeThroughApi)
{
    int meshkernelId = 0;
    int errorCode = meshkernelapi::mkernel_allocate_state(0, meshkernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = GenerateUnstructuredMesh(meshkernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    int newEdgeIndex;
    errorCode = meshkernelapi::mkernel_mesh2d_insert_edge(meshkernelId, 0, 4, newEdgeIndex);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(17, newEdgeIndex);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshkernelId, mesh2d);

    // Assert
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(12, mesh2d.num_nodes);
    ASSERT_EQ(18, mesh2d.num_edges);
}

TEST(Mesh2D, MergeTwoNodesThroughApi)
{
    int meshkernelId = 0;
    int errorCode = meshkernelapi::mkernel_allocate_state(0, meshkernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = GenerateUnstructuredMesh(meshkernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    errorCode = meshkernelapi::mkernel_mesh2d_merge_two_nodes(meshkernelId, 0, 4);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshkernelId, mesh2d);

    // Assert
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(11, mesh2d.num_valid_nodes);
    ASSERT_EQ(15, mesh2d.num_valid_edges);
}

TEST(Mesh2D, MergeNodesThroughApi)
{
    int meshkernelId = 0;
    int errorCode = meshkernelapi::mkernel_allocate_state(0, meshkernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = GenerateUnstructuredMesh(meshkernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::GeometryList geometry_list{};
    std::vector xCoordinates{-0.5, 2.5, 2.5, -0.5, -0.5};
    std::vector yCoordinates{-0.5, -0.5, 2.5, 2.5, -0.5};
    std::vector zCoordinates{0.0, 0.0, 0.0, 0.5, 0.5};

    geometry_list.geometry_separator = meshkernel::constants::missing::doubleValue;
    geometry_list.coordinates_x = xCoordinates.data();
    geometry_list.coordinates_y = yCoordinates.data();
    geometry_list.values = zCoordinates.data();
    geometry_list.num_coordinates = static_cast<int>(xCoordinates.size());

    // Execute
    errorCode = mkernel_mesh2d_merge_nodes(meshkernelId, geometry_list);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshkernelId, mesh2d);

    // Assert (nothing is done, just check that the api communication works)
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(12, mesh2d.num_valid_nodes);
    ASSERT_EQ(17, mesh2d.num_valid_edges);
}

TEST(Mesh2D, MergeNodesWithMergingDistanceThroughApi)
{
    int meshkernelId = 0;
    int errorCode = meshkernelapi::mkernel_allocate_state(0, meshkernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    errorCode = GenerateUnstructuredMesh(meshkernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::GeometryList geometry_list{};

    // Execute
    errorCode = mkernel_mesh2d_merge_nodes_with_merging_distance(meshkernelId, geometry_list, 0.001);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshkernelId, mesh2d);

    // Assert (nothing is done, just check that the api communication works)
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(12, mesh2d.num_valid_nodes);
    ASSERT_EQ(17, mesh2d.num_valid_edges);
}

TEST(Mesh2D, InsertNodeAndEdge_OnMesh2D_ShouldInsertNodeAndEdge)
{
    int meshkernelId = 0;
    int errorCode = meshkernelapi::mkernel_allocate_state(0, meshkernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    errorCode = GenerateUnstructuredMesh(meshkernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    // Isolated nodes are removed by the administration done in mkernel_mesh2d_get_dimensions.
    // The newly inserted node should be connected to another one to form an edge.
    // In this manner, the edge will not be removed during the administration
    int newNodeIndex;
    errorCode = meshkernelapi::mkernel_mesh2d_insert_node(meshkernelId, -0.5, -0.5, newNodeIndex);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    int newEdgeIndex;
    errorCode = meshkernelapi::mkernel_mesh2d_insert_edge(meshkernelId, newNodeIndex, 0, newEdgeIndex);

    // Assert
    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshkernelId, mesh2d);
    ASSERT_EQ(mesh2d.num_nodes, 13);
    ASSERT_EQ(mesh2d.num_edges, 18);
}

TEST(Mesh2D, MoveNode_OnMesh2D_ShouldMoveNode)
{
    int meshkernelId = 0;
    int errorCode = meshkernelapi::mkernel_allocate_state(0, meshkernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    errorCode = GenerateUnstructuredMesh(meshkernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    errorCode = meshkernelapi::mkernel_mesh2d_move_node(meshkernelId, -0.5, -0.5, 0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshkernelId, mesh2d);
    std::vector<int> edge_nodes(mesh2d.num_edges * 2);
    std::vector<int> face_nodes(mesh2d.num_face_nodes);
    std::vector<int> nodes_per_face(mesh2d.num_faces);

    std::vector<double> node_x(mesh2d.num_nodes);
    std::vector<double> node_y(mesh2d.num_nodes);

    std::vector<double> edge_x(mesh2d.num_edges);
    std::vector<double> edge_y(mesh2d.num_edges);

    std::vector<double> face_x(mesh2d.num_faces);
    std::vector<double> face_y(mesh2d.num_faces);

    std::vector<int> edge_faces(mesh2d.num_edges * 2);
    std::vector<int> face_edges(mesh2d.num_face_nodes * 2);

    mesh2d.edge_faces = edge_faces.data();
    mesh2d.edge_nodes = edge_nodes.data();
    mesh2d.face_nodes = face_nodes.data();
    mesh2d.face_edges = face_edges.data();
    mesh2d.nodes_per_face = nodes_per_face.data();
    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_x = edge_x.data();
    mesh2d.edge_y = edge_y.data();
    mesh2d.face_x = face_x.data();
    mesh2d.face_y = face_y.data();
    errorCode = mkernel_mesh2d_get_data(meshkernelId, mesh2d);

    // Assert
    ASSERT_EQ(mesh2d.node_x[0], -0.5);
    ASSERT_EQ(mesh2d.node_y[0], -0.5);
}

TEST(Mesh2D, MoveNode_OnMesh2DWithInvalidIndex_ShouldReturnAnErrorCode)
{
    int meshkernelId = 0;
    int errorCode = meshkernelapi::mkernel_allocate_state(0, meshkernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    errorCode = GenerateUnstructuredMesh(meshkernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    errorCode = meshkernelapi::mkernel_mesh2d_move_node(meshkernelId, -0.5, -0.5, -1);

    // Assert call was not successful
    ASSERT_NE(meshkernel::ExitCode::Success, errorCode);
}

TEST(Mesh2D, GetEdge_OnMesh2D_ShouldGetAnEdgeIndex)
{
    int meshkernelId = 0;
    int errorCode = meshkernelapi::mkernel_allocate_state(0, meshkernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    errorCode = GenerateUnstructuredMesh(meshkernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    int edgeIndex;
    errorCode = meshkernelapi::mkernel_mesh2d_get_edge(meshkernelId, 0.5, -0.5, 0.0, 0.0, 3.0, 3.0, edgeIndex);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    ASSERT_EQ(edgeIndex, 0);
}

TEST(Mesh2D, GetNode_OnMesh2D_ShouldGetANodeIndex)
{
    int meshkernelId = 0;
    int errorCode = meshkernelapi::mkernel_allocate_state(0, meshkernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    errorCode = GenerateUnstructuredMesh(meshkernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    int nodeIndex;
    errorCode = meshkernelapi::mkernel_mesh2d_get_node_index(meshkernelId, 3.0, 3.0, 10.0, 0.0, 0.0, 3.0, 3.0, nodeIndex);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    ASSERT_EQ(nodeIndex, 11);
}
