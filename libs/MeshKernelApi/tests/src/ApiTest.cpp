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

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "MeshKernel/MeshTransformation.hpp"
#include "MeshKernel/Parameters.hpp"

#include "MeshKernelApi/BoundingBox.hpp"
#include "MeshKernelApi/GeometryList.hpp"
#include "MeshKernelApi/Mesh1D.hpp"
#include "MeshKernelApi/Mesh2D.hpp"
#include "MeshKernelApi/MeshKernel.hpp"
#include "Version/Version.hpp"

#include "TestUtils/Definitions.hpp"
#include "TestUtils/MakeCurvilinearGrids.hpp"
#include "TestUtils/MakeMeshes.hpp"
#include "TestUtils/SampleFileReader.hpp"

#include <memory>
#include <numeric>

#include "CartesianApiTestFixture.hpp"

TEST(MeshState, AllocateState)
{
    int errorCode;
    int meshKernelId = 0;
    // cartesian
    errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    errorCode = meshkernelapi::mkernel_deallocate_state(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    // spherical
    errorCode = meshkernelapi::mkernel_allocate_state(1, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    errorCode = meshkernelapi::mkernel_deallocate_state(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    // spherical accurate
    errorCode = meshkernelapi::mkernel_allocate_state(2, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    errorCode = meshkernelapi::mkernel_deallocate_state(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    // invalid
    errorCode = meshkernelapi::mkernel_allocate_state(3, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::RangeErrorCode, errorCode);
    // invalid
    errorCode = meshkernelapi::mkernel_allocate_state(-1, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::RangeErrorCode, errorCode);
}

TEST(MeshState, MKernelGetProjection_ShouldGetProjection)
{
    // Setup
    int meshKernelId = 0;
    int setProjectionType = 0;
    auto errorCode = meshkernelapi::mkernel_allocate_state(setProjectionType, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    int getProjectionType = 0;
    errorCode = meshkernelapi::mkernel_get_projection(meshKernelId, getProjectionType);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    ASSERT_EQ(setProjectionType, getProjectionType);
}

TEST_F(CartesianApiTestFixture, Mesh2DDeleteNode_ShouldDeleteNode)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    auto errorCode = meshkernelapi::mkernel_mesh2d_delete_node(meshKernelId, 0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Get the dimensions
    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert dimensions
    ASSERT_EQ(11, mesh2d.num_valid_nodes);
    ASSERT_EQ(15, mesh2d.num_valid_edges);

    // Allocate memory and get data
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
    errorCode = mkernel_mesh2d_get_data(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    /*  1---4---7---10
        |   |   |   |
        0---3---6---9
            |   |   |
            2---5---8
    */

    /*  2---5---8---11
        |   |   |   |
        1---4---7---10
            |   |   |
            3---6---9

        Node 0 is invalid
        Edges 0 and 9 are invalid

    */

    // Assert data
    const double tolerance = 1e-6;
    // Nodes
    EXPECT_NEAR(0.0, mesh2d.node_x[1], tolerance);
    EXPECT_NEAR(1.0, mesh2d.node_y[1], tolerance);
    // Edges
    EXPECT_EQ(1, mesh2d.edge_nodes[2]);
    EXPECT_EQ(4, mesh2d.edge_nodes[3]);
    // TODO should the edge_x/y have the same number of values as there are edges
    EXPECT_EQ(meshkernel::constants::missing::doubleValue, mesh2d.edge_x[0]);
    EXPECT_EQ(meshkernel::constants::missing::doubleValue, mesh2d.edge_y[0]);

    EXPECT_NEAR(0.5, mesh2d.edge_x[1], tolerance);
    EXPECT_NEAR(1.0, mesh2d.edge_y[1], tolerance);
    // First face
    EXPECT_EQ(0, mesh2d.edge_faces[2]);
    EXPECT_EQ(-1, mesh2d.edge_faces[3]);
    EXPECT_EQ(0, mesh2d.edge_faces[4]);
    EXPECT_EQ(-1, mesh2d.edge_faces[5]);

    EXPECT_EQ(1, mesh2d.face_edges[0]);
    EXPECT_EQ(12, mesh2d.face_edges[1]);
    EXPECT_EQ(2, mesh2d.face_edges[2]);
    EXPECT_EQ(10, mesh2d.face_edges[3]);

    EXPECT_EQ(4, mesh2d.nodes_per_face[0]);
    EXPECT_NEAR(0.5, mesh2d.face_x[0], tolerance);
    EXPECT_NEAR(1.5, mesh2d.face_y[0], tolerance);
    // Second Face
    EXPECT_EQ(3, mesh2d.face_nodes[4]);
    EXPECT_EQ(6, mesh2d.face_nodes[5]);
    EXPECT_EQ(7, mesh2d.face_nodes[6]);
    EXPECT_EQ(4, mesh2d.face_nodes[7]);
    EXPECT_EQ(4, mesh2d.nodes_per_face[1]);

    EXPECT_NEAR(1.5, mesh2d.face_x[1], tolerance);
    EXPECT_NEAR(0.5, mesh2d.face_y[1], tolerance);
}

TEST_F(CartesianApiTestFixture, FlipEdges_ShouldFlipEdges)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    const int isTriangulationRequired = 1;
    const int projectToLandBoundaryOption = 1;
    meshkernelapi::GeometryList selectingPolygon{};
    meshkernelapi::GeometryList landBoundaries{};
    auto errorCode = mkernel_mesh2d_flip_edges(meshKernelId,
                                               isTriangulationRequired,
                                               projectToLandBoundaryOption,
                                               selectingPolygon,
                                               landBoundaries);

    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);

    // Assert
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(12, mesh2d.num_nodes);
    ASSERT_EQ(23, mesh2d.num_edges);
}

TEST_F(CartesianApiTestFixture, FlipEdges_WithALandBoundary_ShouldFlipEdges)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    const int isTriangulationRequired = 1;
    const int projectToLandBoundaryOption = 1;
    meshkernelapi::GeometryList selectingPolygon{};

    std::vector xCoordinates{-0.5, -0.5, 4.0, meshkernel::constants::missing::doubleValue};
    std::vector yCoordinates{3.0, -0.5, -0.5, meshkernel::constants::missing::doubleValue};
    std::vector zCoordinates{0.0, 0.0, 0.0, meshkernel::constants::missing::doubleValue};

    meshkernelapi::GeometryList landBoundaries{};
    landBoundaries.geometry_separator = meshkernel::constants::missing::doubleValue;
    landBoundaries.coordinates_x = xCoordinates.data();
    landBoundaries.coordinates_y = yCoordinates.data();
    landBoundaries.values = zCoordinates.data();
    landBoundaries.num_coordinates = static_cast<int>(xCoordinates.size());

    auto errorCode = mkernel_mesh2d_flip_edges(meshKernelId,
                                               isTriangulationRequired,
                                               projectToLandBoundaryOption,
                                               selectingPolygon,
                                               landBoundaries);

    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);

    // Assert
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(12, mesh2d.num_nodes);
    ASSERT_EQ(23, mesh2d.num_edges);
}

TEST_F(CartesianApiTestFixture, GenerateTriangularGridThroughApi)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometry_separator = meshkernel::constants::missing::doubleValue;
    std::vector xCoordinates{
        415.319672,
        390.271973,
        382.330048,
        392.715668,
        418.374268,
        453.807556,
        495.960968,
        532.005188,
        565.605774,
        590.653442,
        598.595398,
        593.708008,
        564.994812,
        514.899475,
        461.138611,
        422.039764,
        415.319672};

    std::vector yCoordinates{
        490.293762,
        464.024139,
        438.365448,
        411.484894,
        386.437103,
        366.276703,
        363.222107,
        370.553162,
        386.437103,
        412.095825,
        445.085571,
        481.129944,
        497.624817,
        504.955872,
        501.290344,
        493.348358,
        490.293762};

    std::vector zCoordinates{
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0};

    geometryListIn.coordinates_x = xCoordinates.data();
    geometryListIn.coordinates_y = yCoordinates.data();
    geometryListIn.values = zCoordinates.data();
    geometryListIn.num_coordinates = static_cast<int>(xCoordinates.size());

    // Execute
    auto errorCode = mkernel_mesh2d_make_triangular_mesh_from_polygon(meshKernelId, geometryListIn);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    // Get the new state

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    ASSERT_EQ(23, mesh2d.num_nodes);
    ASSERT_EQ(50, mesh2d.num_edges);
}

TEST_F(CartesianApiTestFixture, GenerateTriangularGridFromSamplesThroughApi)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;

    geometryListIn.geometry_separator = meshkernel::constants::missing::doubleValue;

    std::vector xCoordinates{
        0.0,
        10.0,
        10.0,
        0.0,
        0.0};

    std::vector yCoordinates{
        0.0,
        0.0,
        10.0,
        10.0,
        0.0};

    std::vector zCoordinates{
        0.0,
        0.0,
        0.0,
        0.0,
        0.0};

    geometryListIn.coordinates_x = xCoordinates.data();
    geometryListIn.coordinates_y = yCoordinates.data();
    geometryListIn.values = zCoordinates.data();

    geometryListIn.num_coordinates = static_cast<int>(xCoordinates.size());

    // Execute
    auto errorCode = mkernel_mesh2d_make_triangular_mesh_from_samples(meshKernelId, geometryListIn);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Get the new state

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    ASSERT_EQ(4, mesh2d.num_nodes);
    ASSERT_EQ(5, mesh2d.num_edges);
}

TEST_F(CartesianApiTestFixture, GetMeshBoundariesThroughApi)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();
    int numberOfpolygonNodes;
    meshkernelapi::GeometryList selectingPolygon;
    auto errorCode = mkernel_mesh2d_count_mesh_boundaries_as_polygons(meshKernelId, selectingPolygon, numberOfpolygonNodes);
    ASSERT_EQ(11, numberOfpolygonNodes);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::GeometryList geometryListOut;
    geometryListOut.geometry_separator = meshkernel::constants::missing::doubleValue;
    geometryListOut.num_coordinates = numberOfpolygonNodes;

    std::vector<double> xCoordinates(numberOfpolygonNodes);
    std::vector<double> yCoordinates(numberOfpolygonNodes);
    std::vector<double> zCoordinates(numberOfpolygonNodes);

    geometryListOut.coordinates_x = xCoordinates.data();
    geometryListOut.coordinates_y = yCoordinates.data();
    geometryListOut.values = zCoordinates.data();

    // Execute
    errorCode = mkernel_mesh2d_get_mesh_boundaries_as_polygons(meshKernelId, selectingPolygon, geometryListOut);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    const double tolerance = 1e-6;
    ASSERT_NEAR(0.0, geometryListOut.coordinates_x[0], tolerance);
    ASSERT_NEAR(0.0, geometryListOut.coordinates_y[0], tolerance);
}

TEST_F(CartesianApiTestFixture, GetMeshBoundariesAsPolygons_WithSelectedPolygon_ShouldGetSelectedMeshBoundary)
{
    // Prepare
    MakeMesh(10, 10, 1);
    auto const meshKernelId = GetMeshKernelId();

    int nodeIndex = 0;
    auto errorCode = meshkernelapi::mkernel_mesh2d_get_node_index(meshKernelId,
                                                                  4.0,
                                                                  4.0,
                                                                  0.1,
                                                                  0.0,
                                                                  0.0,
                                                                  10.0,
                                                                  10.0,
                                                                  nodeIndex);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    errorCode = meshkernelapi::mkernel_mesh2d_delete_node(meshKernelId, nodeIndex);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    int numberOfBoundaryPolygonNodes;
    meshkernelapi::GeometryList selectingPolygon;
    auto coordinates_x = std::vector<double>{2, 6, 6, 2, 2};
    auto coordinates_y = std::vector<double>{2, 2, 6, 6, 2};
    selectingPolygon.coordinates_x = coordinates_x.data();
    selectingPolygon.coordinates_y = coordinates_y.data();
    selectingPolygon.num_coordinates = static_cast<int>(coordinates_x.size());

    errorCode = mkernel_mesh2d_count_mesh_boundaries_as_polygons(meshKernelId, selectingPolygon, numberOfBoundaryPolygonNodes);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::GeometryList BoundaryPolygon;
    BoundaryPolygon.geometry_separator = meshkernel::constants::missing::doubleValue;
    BoundaryPolygon.num_coordinates = numberOfBoundaryPolygonNodes;

    std::vector<double> xCoordinatesBoundaryPolygon(numberOfBoundaryPolygonNodes);
    std::vector<double> yCoordinatesBoundaryPolygon(numberOfBoundaryPolygonNodes);
    std::vector<double> zCoordinatesBoundaryPolygon(numberOfBoundaryPolygonNodes);

    BoundaryPolygon.coordinates_x = xCoordinatesBoundaryPolygon.data();
    BoundaryPolygon.coordinates_y = yCoordinatesBoundaryPolygon.data();
    BoundaryPolygon.values = zCoordinatesBoundaryPolygon.data();

    // Execute
    errorCode = mkernel_mesh2d_get_mesh_boundaries_as_polygons(meshKernelId, selectingPolygon, BoundaryPolygon);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    const double tolerance = 1e-6;
    ASSERT_NEAR(3.0, BoundaryPolygon.coordinates_x[0], tolerance);
    ASSERT_NEAR(4.0, BoundaryPolygon.coordinates_x[1], tolerance);
    ASSERT_NEAR(5.0, BoundaryPolygon.coordinates_x[2], tolerance);
    ASSERT_NEAR(5.0, BoundaryPolygon.coordinates_x[3], tolerance);
    ASSERT_NEAR(5.0, BoundaryPolygon.coordinates_x[4], tolerance);
    ASSERT_NEAR(4.0, BoundaryPolygon.coordinates_x[5], tolerance);
    ASSERT_NEAR(3.0, BoundaryPolygon.coordinates_x[6], tolerance);
    ASSERT_NEAR(3.0, BoundaryPolygon.coordinates_x[7], tolerance);
    ASSERT_NEAR(3.0, BoundaryPolygon.coordinates_x[8], tolerance);

    ASSERT_NEAR(3.0, BoundaryPolygon.coordinates_y[0], tolerance);
    ASSERT_NEAR(3.0, BoundaryPolygon.coordinates_y[1], tolerance);
    ASSERT_NEAR(3.0, BoundaryPolygon.coordinates_y[2], tolerance);
    ASSERT_NEAR(4.0, BoundaryPolygon.coordinates_y[3], tolerance);
    ASSERT_NEAR(5.0, BoundaryPolygon.coordinates_y[4], tolerance);
    ASSERT_NEAR(5.0, BoundaryPolygon.coordinates_y[5], tolerance);
    ASSERT_NEAR(5.0, BoundaryPolygon.coordinates_y[6], tolerance);
    ASSERT_NEAR(4.0, BoundaryPolygon.coordinates_y[7], tolerance);
    ASSERT_NEAR(3.0, BoundaryPolygon.coordinates_y[8], tolerance);
}

TEST_F(CartesianApiTestFixture, OffsetAPolygonThroughApi)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometry_separator = meshkernel::constants::missing::doubleValue;

    std::vector xCoordinatesIn{0.0, 1.0, 1.0, 0.0, 0.0};
    std::vector yCoordinatesIn{0.0, 0.0, 1.0, 1.0, 0.0};
    std::vector valuesIn{0.0, 0.0, 0.0, 0.0, 0.0};

    geometryListIn.coordinates_x = xCoordinatesIn.data();
    geometryListIn.coordinates_y = yCoordinatesIn.data();
    geometryListIn.values = valuesIn.data();
    geometryListIn.num_coordinates = static_cast<int>(xCoordinatesIn.size());

    // Execute
    int numberOfpolygonNodes;
    auto errorCode = mkernel_polygon_count_offset(meshKernelId, geometryListIn, false, 0.5, numberOfpolygonNodes);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(5, numberOfpolygonNodes);

    meshkernelapi::GeometryList geometryListOut;

    geometryListOut.geometry_separator = meshkernel::constants::missing::doubleValue;

    std::vector<double> xCoordinatesOut(numberOfpolygonNodes);
    std::vector<double> yCoordinatesOut(numberOfpolygonNodes);
    std::vector<double> valuesOut(numberOfpolygonNodes);
    geometryListOut.coordinates_x = xCoordinatesOut.data();
    geometryListOut.coordinates_y = yCoordinatesOut.data();
    geometryListOut.values = valuesOut.data();
    geometryListOut.num_coordinates = static_cast<int>(xCoordinatesOut.size());

    errorCode = mkernel_polygon_get_offset(meshKernelId, geometryListIn, false, 10.0, geometryListOut);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    const double tolerance = 1e-6;

    ASSERT_NEAR(-10.0, geometryListOut.coordinates_x[0], tolerance);
    ASSERT_NEAR(11.0, geometryListOut.coordinates_x[1], tolerance);
    ASSERT_NEAR(11.0, geometryListOut.coordinates_x[2], tolerance);
    ASSERT_NEAR(-10.0, geometryListOut.coordinates_x[3], tolerance);
    ASSERT_NEAR(-10.0, geometryListOut.coordinates_x[4], tolerance);

    ASSERT_NEAR(-10.0, geometryListOut.coordinates_y[0], tolerance);
    ASSERT_NEAR(-10.0, geometryListOut.coordinates_y[1], tolerance);
    ASSERT_NEAR(11.0, geometryListOut.coordinates_y[2], tolerance);
    ASSERT_NEAR(11.0, geometryListOut.coordinates_y[3], tolerance);
    ASSERT_NEAR(-10.0, geometryListOut.coordinates_y[4], tolerance);
}

TEST_F(CartesianApiTestFixture, ComputeSingleContactsThroughApi_ShouldGenerateContacts)
{
    // Prepare
    MakeMesh(3, 3, 10);
    auto const meshKernelId = GetMeshKernelId();

    // Init 1d mesh
    meshkernelapi::Mesh1D mesh1d;
    std::vector node_x{
        1.73493900000000,
        2.35659313023165,
        5.38347452702839,
        14.2980910429074,
        22.9324017677239,
        25.3723169493137,
        25.8072280000000};
    std::vector node_y{
        -7.6626510000000,
        1.67281447902331,
        10.3513746546384,
        12.4797224193970,
        15.3007317677239,
        24.1623588554512,
        33.5111870000000};
    mesh1d.node_x = node_x.data();
    mesh1d.node_y = node_y.data();
    mesh1d.num_nodes = static_cast<int>(node_x.size());

    std::vector edge_nodes{0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6};
    mesh1d.edge_nodes = edge_nodes.data();
    mesh1d.num_edges = 6;

    auto errorCode = mkernel_mesh1d_set(meshKernelId, mesh1d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Init 1d mask
    std::vector onedNodeMask{1, 1, 1, 1, 1, 1, 1};

    // Init polygon
    meshkernelapi::GeometryList polygon;
    polygon.geometry_separator = meshkernel::constants::missing::doubleValue;

    std::vector<double> xCoordinates{-30, 40, 40, -40, -30};
    std::vector<double> yCoordinates{-20, -20, 50, 50, -20};
    std::vector<double> zCoordinates{0, 0, 0, 0, 0};
    polygon.coordinates_x = xCoordinates.data();
    polygon.coordinates_y = yCoordinates.data();
    polygon.values = zCoordinates.data();

    polygon.num_coordinates = static_cast<int>(xCoordinates.size());

    // Execute
    errorCode = mkernel_contacts_compute_single(meshKernelId, onedNodeMask.data(), polygon, 0.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Get the new state
    meshkernelapi::Contacts contacts{};
    errorCode = mkernel_contacts_get_dimensions(meshKernelId, contacts);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<int> mesh1d_indices(contacts.num_contacts);
    std::vector<int> mesh2d_indices(contacts.num_contacts);
    contacts.mesh1d_indices = mesh1d_indices.data();
    contacts.mesh2d_indices = mesh2d_indices.data();

    errorCode = mkernel_contacts_get_data(meshKernelId, contacts);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    ASSERT_EQ(5, contacts.num_contacts);

    ASSERT_EQ(1, contacts.mesh1d_indices[0]);
    ASSERT_EQ(2, contacts.mesh1d_indices[1]);
    ASSERT_EQ(3, contacts.mesh1d_indices[2]);
    ASSERT_EQ(4, contacts.mesh1d_indices[3]);
    ASSERT_EQ(5, contacts.mesh1d_indices[4]);

    ASSERT_EQ(0, contacts.mesh2d_indices[0]);
    ASSERT_EQ(1, contacts.mesh2d_indices[1]);
    ASSERT_EQ(4, contacts.mesh2d_indices[2]);
    ASSERT_EQ(7, contacts.mesh2d_indices[3]);
    ASSERT_EQ(8, contacts.mesh2d_indices[4]);
}

TEST_F(CartesianApiTestFixture, ComputeMultipleContactsThroughApi)
{
    // Prepare
    MakeMesh(3, 3, 10);
    auto const meshKernelId = GetMeshKernelId();

    // Init 1d mesh
    meshkernelapi::Mesh1D mesh1d;
    std::vector node_x{
        1.73493900000000,
        2.35659313023165,
        5.38347452702839,
        14.2980910429074,
        22.9324017677239,
        25.3723169493137,
        25.8072280000000};
    std::vector node_y{
        -7.6626510000000,
        1.67281447902331,
        10.3513746546384,
        12.4797224193970,
        15.3007317677239,
        24.1623588554512,
        33.5111870000000};
    mesh1d.node_x = node_x.data();
    mesh1d.node_y = node_y.data();
    mesh1d.num_nodes = 7;

    std::vector edge_nodes{0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6};
    mesh1d.edge_nodes = edge_nodes.data();
    mesh1d.num_edges = 6;

    auto errorCode = mkernel_mesh1d_set(meshKernelId, mesh1d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Init 1d mask
    std::vector onedNodeMask{1, 1, 1, 1, 1, 1, 1};

    // Execute
    // Why do we need to specify the namespace "meshkernelapi"?
    errorCode = meshkernelapi::mkernel_contacts_compute_multiple(meshKernelId, onedNodeMask.data());
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Get the new state
    meshkernelapi::Contacts contacts{};
    errorCode = mkernel_contacts_get_dimensions(meshKernelId, contacts);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<int> mesh1d_indices(contacts.num_contacts);
    std::vector<int> mesh2d_indices(contacts.num_contacts);
    contacts.mesh1d_indices = mesh1d_indices.data();
    contacts.mesh2d_indices = mesh2d_indices.data();

    errorCode = mkernel_contacts_get_data(meshKernelId, contacts);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    ASSERT_EQ(5, contacts.num_contacts);

    ASSERT_EQ(1, contacts.mesh1d_indices[0]);
    ASSERT_EQ(2, contacts.mesh1d_indices[1]);
    ASSERT_EQ(3, contacts.mesh1d_indices[2]);
    ASSERT_EQ(4, contacts.mesh1d_indices[3]);
    ASSERT_EQ(5, contacts.mesh1d_indices[4]);

    ASSERT_EQ(0, contacts.mesh2d_indices[0]);
    ASSERT_EQ(1, contacts.mesh2d_indices[1]);
    ASSERT_EQ(4, contacts.mesh2d_indices[2]);
    ASSERT_EQ(7, contacts.mesh2d_indices[3]);
    ASSERT_EQ(8, contacts.mesh2d_indices[4]);
}

TEST_F(CartesianApiTestFixture, ComputeContactsWithPolygonsThroughApi)
{
    // Prepare
    MakeMesh(3, 3, 10);
    auto const meshKernelId = GetMeshKernelId();

    // Init 1d mesh
    meshkernelapi::Mesh1D mesh1d;
    std::vector node_x{
        1.73493900000000,
        2.35659313023165,
        5.38347452702839,
        14.2980910429074,
        22.9324017677239,
        25.3723169493137,
        25.8072280000000};
    std::vector node_y{
        -7.6626510000000,
        1.67281447902331,
        10.3513746546384,
        12.4797224193970,
        15.3007317677239,
        24.1623588554512,
        33.5111870000000};
    mesh1d.node_x = node_x.data();
    mesh1d.node_y = node_y.data();
    mesh1d.num_nodes = 7;

    std::vector edge_nodes{
        0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6};
    mesh1d.edge_nodes = edge_nodes.data();
    mesh1d.num_edges = 6;

    auto errorCode = mkernel_mesh1d_set(meshKernelId, mesh1d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Init 1d mask
    std::vector onedNodeMask{1, 1, 1, 1, 1, 1, 1};

    // Init polygon
    meshkernelapi::GeometryList polygon;
    polygon.geometry_separator = meshkernel::constants::missing::doubleValue;

    std::vector<double> xCoordinates{25, 50, 50, 25, 25};
    std::vector<double> yCoordinates{25, 25, 50, 50, 25};
    std::vector<double> zCoordinates{0, 0, 0, 0, 0};
    polygon.coordinates_x = xCoordinates.data();
    polygon.coordinates_y = yCoordinates.data();
    polygon.values = zCoordinates.data();
    polygon.num_coordinates = static_cast<int>(xCoordinates.size());

    // Execute
    errorCode = mkernel_contacts_compute_with_polygons(meshKernelId, onedNodeMask.data(), polygon);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Get the new state
    meshkernelapi::Contacts contacts{};
    errorCode = mkernel_contacts_get_dimensions(meshKernelId, contacts);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<int> mesh1d_indices(contacts.num_contacts);
    std::vector<int> mesh2d_indices(contacts.num_contacts);
    contacts.mesh1d_indices = mesh1d_indices.data();
    contacts.mesh2d_indices = mesh2d_indices.data();

    errorCode = mkernel_contacts_get_data(meshKernelId, contacts);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    ASSERT_EQ(1, contacts.num_contacts);
    ASSERT_EQ(5, contacts.mesh1d_indices[0]);
    ASSERT_EQ(8, contacts.mesh2d_indices[0]);
}

TEST_F(CartesianApiTestFixture, ComputeContactsWithPointsThroughApi)
{
    // Prepare
    MakeMesh(3, 3, 10);
    auto const meshKernelId = GetMeshKernelId();

    // Init 1d mesh
    meshkernelapi::Mesh1D mesh1d;
    std::vector node_x{
        1.73493900000000,
        2.35659313023165,
        5.38347452702839,
        14.2980910429074,
        22.9324017677239,
        25.3723169493137,
        25.8072280000000};
    std::vector node_y{
        -7.6626510000000,
        1.67281447902331,
        10.3513746546384,
        12.4797224193970,
        15.3007317677239,
        24.1623588554512,
        33.5111870000000};
    mesh1d.node_x = node_x.data();
    mesh1d.node_y = node_y.data();
    mesh1d.num_nodes = static_cast<int>(node_x.size());

    std::vector edge_nodes{
        0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6};
    mesh1d.edge_nodes = edge_nodes.data();
    mesh1d.num_edges = 6;

    auto errorCode = mkernel_mesh1d_set(meshKernelId, mesh1d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Init 1d mask
    std::vector onedNodeMask{1, 1, 1, 1, 1, 1, 1};

    // Init polygon
    meshkernelapi::GeometryList points;
    points.geometry_separator = meshkernel::constants::missing::doubleValue;

    std::vector<double> xCoordinates{5, 15, 25, 35};
    std::vector<double> yCoordinates{5, 15, 25, 35};
    std::vector<double> zCoordinates{0, 0, 0, 0};
    points.coordinates_x = xCoordinates.data();
    points.coordinates_y = yCoordinates.data();
    points.values = zCoordinates.data();
    points.num_coordinates = static_cast<int>(xCoordinates.size());

    // Execute
    errorCode = mkernel_contacts_compute_with_points(meshKernelId, onedNodeMask.data(), points);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Get the new state
    meshkernelapi::Contacts contacts{};
    errorCode = mkernel_contacts_get_dimensions(meshKernelId, contacts);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<int> mesh1d_indices(contacts.num_contacts);
    std::vector<int> mesh2d_indices(contacts.num_contacts);
    contacts.mesh1d_indices = mesh1d_indices.data();
    contacts.mesh2d_indices = mesh2d_indices.data();

    errorCode = mkernel_contacts_get_data(meshKernelId, contacts);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    ASSERT_EQ(3, contacts.num_contacts);

    ASSERT_EQ(1, contacts.mesh1d_indices[0]);
    ASSERT_EQ(3, contacts.mesh1d_indices[1]);
    ASSERT_EQ(5, contacts.mesh1d_indices[2]);

    ASSERT_EQ(0, contacts.mesh2d_indices[0]);
    ASSERT_EQ(4, contacts.mesh2d_indices[1]);
    ASSERT_EQ(8, contacts.mesh2d_indices[2]);
}

TEST_F(CartesianApiTestFixture, ComputeBoundaryContactsThroughApi)
{
    // Prepare
    MakeMesh(3, 3, 10);
    auto const meshKernelId = GetMeshKernelId();

    // Init 1d mesh
    meshkernelapi::Mesh1D mesh1d;
    std::vector node_x{
        -16.1886410000000,
        -16.1464995876014,
        -16.1043581752028,
        -16.0622167628042,
        -15.7539488236928,
        -6.86476658679268,
        2.02441565010741,
        10.9135970000000,
    };
    std::vector node_y{
        0.89018900000000,
        9.78201442138723,
        18.6738398427745,
        27.5656652641617,
        36.1966603330179,
        36.4175095626911,
        36.6383587923643,
        36.8592080000000};
    mesh1d.node_x = node_x.data();
    mesh1d.node_y = node_y.data();
    mesh1d.num_nodes = 8;

    std::vector edge_nodes{0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7};
    mesh1d.edge_nodes = edge_nodes.data();
    mesh1d.num_edges = 7;

    auto errorCode = mkernel_mesh1d_set(meshKernelId, mesh1d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Init 1d mask
    std::vector onedNodeMask{1, 1, 1, 1, 1, 1, 1, 1};

    // Init polygon
    meshkernelapi::GeometryList polygon;
    polygon.geometry_separator = meshkernel::constants::missing::doubleValue;

    std::vector<double> xCoordinates{-30, 40, 40, -40, -30};
    std::vector<double> yCoordinates{-20, -20, 50, 50, -20};
    std::vector<double> zCoordinates{0, 0, 0, 0, 0};
    polygon.coordinates_x = xCoordinates.data();
    polygon.coordinates_y = yCoordinates.data();
    polygon.values = zCoordinates.data();
    polygon.num_coordinates = static_cast<int>(xCoordinates.size());

    // Execute
    errorCode = mkernel_contacts_compute_boundary(meshKernelId, onedNodeMask.data(), polygon, 200.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Get the new state
    meshkernelapi::Contacts contacts{};
    errorCode = mkernel_contacts_get_dimensions(meshKernelId, contacts);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector mesh1d_indices(contacts.num_contacts, 0);
    std::vector mesh2d_indices(contacts.num_contacts, 0);
    contacts.mesh1d_indices = mesh1d_indices.data();
    contacts.mesh2d_indices = mesh2d_indices.data();

    errorCode = mkernel_contacts_get_data(meshKernelId, contacts);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    ASSERT_EQ(8, contacts.num_contacts);

    ASSERT_EQ(0, contacts.mesh1d_indices[0]);
    ASSERT_EQ(2, contacts.mesh1d_indices[1]);
    ASSERT_EQ(6, contacts.mesh1d_indices[2]);
    ASSERT_EQ(0, contacts.mesh1d_indices[3]);
    ASSERT_EQ(7, contacts.mesh1d_indices[4]);
    ASSERT_EQ(7, contacts.mesh1d_indices[5]);
    ASSERT_EQ(7, contacts.mesh1d_indices[6]);
    ASSERT_EQ(7, contacts.mesh1d_indices[7]);

    ASSERT_EQ(0, contacts.mesh2d_indices[0]);
    ASSERT_EQ(1, contacts.mesh2d_indices[1]);
    ASSERT_EQ(2, contacts.mesh2d_indices[2]);
    ASSERT_EQ(3, contacts.mesh2d_indices[3]);
    ASSERT_EQ(5, contacts.mesh2d_indices[4]);
    ASSERT_EQ(6, contacts.mesh2d_indices[5]);
    ASSERT_EQ(7, contacts.mesh2d_indices[6]);
    ASSERT_EQ(8, contacts.mesh2d_indices[7]);
}

TEST(ApiStatelessTests, GetSplinesThroughApi)
{
    // Prepare
    meshkernelapi::GeometryList geometryListIn;

    std::vector<double> splineCoordinatesX{10.0, 20.0, 30.0};
    std::vector<double> splineCoordinatesY{-5.0, 5.0, -5.0};
    std::vector<double> values{0.0, 0.0, 0.0};
    geometryListIn.coordinates_x = splineCoordinatesX.data();
    geometryListIn.coordinates_y = splineCoordinatesY.data();
    geometryListIn.values = values.data();
    geometryListIn.num_coordinates = static_cast<int>(splineCoordinatesX.size());
    geometryListIn.geometry_separator = meshkernel::constants::missing::doubleValue;

    meshkernelapi::GeometryList geometryListOut;
    int const numberOfPointsBetweenNodes = 3;
    size_t totalNumPoints = (numberOfPointsBetweenNodes + 1) * 2 + 1;
    std::vector<double> CoordinatesOutX(totalNumPoints);
    std::vector<double> CoordinatesOutY(totalNumPoints);
    std::vector<double> valuesOut(totalNumPoints);
    geometryListOut.coordinates_x = CoordinatesOutX.data();
    geometryListOut.coordinates_y = CoordinatesOutY.data();
    geometryListOut.values = valuesOut.data();
    geometryListOut.num_coordinates = static_cast<int>(totalNumPoints);
    geometryListOut.geometry_separator = meshkernel::constants::missing::doubleValue;

    // Execute
    auto errorCode = mkernel_get_splines(geometryListIn, geometryListOut, numberOfPointsBetweenNodes);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert. The last value is a geometry separator  to handle the case of multiple splines
    ASSERT_EQ(totalNumPoints, geometryListOut.num_coordinates);

    std::vector computedCoordinatesX(CoordinatesOutX.data(), CoordinatesOutX.data() + totalNumPoints);
    std::vector ValidCoordinatesX{10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0};
    ASSERT_THAT(computedCoordinatesX, ::testing::ContainerEq(ValidCoordinatesX));

    std::vector computedCoordinatesY(CoordinatesOutY.data(), CoordinatesOutY.data() + totalNumPoints);
    std::vector ValidCoordinatesY{-5.000000, -1.328125, 1.8750000, 4.1406250, 5.0000000, 4.1406250, 1.8750000, -1.328125, -5.000000};
    ASSERT_THAT(computedCoordinatesY, ::testing::ContainerEq(ValidCoordinatesY));
}

TEST(ApiStatelessTests, TestGettingVersionThroughApi)
{
    auto versionFromApi = std::make_unique<char[]>(64);
    meshkernelapi::mkernel_get_version(versionFromApi.get());
    ASSERT_EQ(strcmp(versionFromApi.get(), versionString), 0);
}

TEST_F(CartesianApiTestFixture, GetClosestMeshCoordinateThroughApi)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    double xCoordinatesOut, yCoordinatesOut;
    auto errorCode = meshkernelapi::mkernel_mesh2d_get_closest_node(meshKernelId, -5.0, 5.0, 10.0, 0, 0, 0, 0, xCoordinatesOut, yCoordinatesOut);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    ASSERT_EQ(0.0, xCoordinatesOut);
}

TEST_F(CartesianApiTestFixture, CurvilinearComputeTransfiniteFromTriangle_OnValidSelection_ShouldCreateCurvilinearGrid)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometry_separator = meshkernel::constants::missing::doubleValue;
    std::vector xCoordinatesIn{
        444.504791,
        427.731781,
        405.640503,
        381.094666,
        451.050354,
        528.778931,
        593.416260,
        558.643005,
        526.733398,
        444.095703,
        444.504791};
    std::vector yCoordinatesIn{
        437.155945,
        382.745758,
        317.699005,
        262.470612,
        262.879700,
        263.288788,
        266.561584,
        324.653687,
        377.836578,
        436.746857,
        437.155945};
    std::vector valuesIn{
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0};
    geometryListIn.coordinates_x = xCoordinatesIn.data();
    geometryListIn.coordinates_y = yCoordinatesIn.data();
    geometryListIn.values = valuesIn.data();
    geometryListIn.num_coordinates = static_cast<int>(xCoordinatesIn.size());

    // Execute
    auto errorCode = mkernel_curvilinear_compute_transfinite_from_triangle(meshKernelId,
                                                                           geometryListIn,
                                                                           0,
                                                                           3,
                                                                           6);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = meshkernelapi::mkernel_curvilinear_convert_to_mesh2d(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);

    // Assert
    ASSERT_EQ(28, mesh2d.num_nodes);
    ASSERT_EQ(40, mesh2d.num_edges);
}

TEST_F(CartesianApiTestFixture, CurvilinearComputeTransfiniteFromTriangle_OnInvalidSelection_ShouldEmitErrorMessage)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();
    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometry_separator = meshkernel::constants::missing::doubleValue;
    std::vector xCoordinatesIn{

        624.810178166822,
        523.558231396829,
        422.306284626837,
        321.054337856845,
        219.802391086853,
        118.550444316861,
        48.2615856660667,
        237.80030632643,
        427.339026986793,
        616.877747647156,
        806.416468307519,
        995.955188967883,
        1185.49390962825,
        1375.03263028861,
        1564.57135094897,
        1430.31975483724,
        1296.0681587255,
        1161.81656261376,
        1027.56496650203,
        893.313370390293,
        759.061774278557,
        624.810178166822};

    std::vector yCoordinatesIn{
        1678.91491569813,
        1518.24565884755,
        1357.57640199696,
        1196.90714514637,
        1036.23788829579,
        875.568631445201,
        733.274504823614,
        745.175773330195,
        757.077041836776,
        768.978310343356,
        780.879578849938,
        792.780847356519,
        804.6821158631,
        816.58338436968,
        828.484652876261,
        949.974690422243,
        1071.46472796822,
        1192.95476551421,
        1314.44480306019,
        1435.93484060617,
        1557.42487815215,
        1678.91491569813};

    std::vector valuesIn{
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0};

    geometryListIn.coordinates_x = xCoordinatesIn.data();
    geometryListIn.coordinates_y = yCoordinatesIn.data();
    geometryListIn.values = valuesIn.data();
    geometryListIn.num_coordinates = static_cast<int>(xCoordinatesIn.size());

    // Execute
    auto errorCode = mkernel_curvilinear_compute_transfinite_from_triangle(meshKernelId,
                                                                           geometryListIn,
                                                                           0,
                                                                           6,
                                                                           9);

    // Assert
    ASSERT_EQ(meshkernel::ExitCode::StdLibExceptionCode, errorCode);

    auto exceptionMessage = std::make_unique<char[]>(512);
    errorCode = meshkernelapi::mkernel_get_error(exceptionMessage.get());
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::string exceptionMessageString(exceptionMessage.get());
    ASSERT_FALSE(exceptionMessageString.empty());
}

TEST_F(CartesianApiTestFixture, Delete_WithEmptyPolygon_ShouldDeleteMesh2D)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryList{};

    // Execute
    auto errorCode = mkernel_mesh2d_delete(meshKernelId, geometryList, 0, false);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);

    // Assert
    ASSERT_EQ(0, mesh2d.num_valid_nodes);
    ASSERT_EQ(0, mesh2d.num_valid_edges);
}

TEST_F(CartesianApiTestFixture, DeleteFaces_WithPolygon_ShouldDeleteMesh2D)
{
    // Prepare
    MakeMesh(4, 4, 2);
    auto const meshKernelId = GetMeshKernelId();

    // By using an empty list, all nodes will be selected
    meshkernelapi::GeometryList geometryList{};

    geometryList.num_coordinates = 5;
    geometryList.geometry_separator = meshkernel::constants::missing::doubleValue;
    std::vector<double> xCoordinatesOut{2, 6, 6, 2, 2};
    std::vector<double> yCoordinatesOut{2, 2, 6, 6, 2};
    geometryList.coordinates_x = xCoordinatesOut.data();
    geometryList.coordinates_y = yCoordinatesOut.data();

    // Execute
    int deletionOption = 2;
    auto errorCode = mkernel_mesh2d_delete(meshKernelId, geometryList, deletionOption, false);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);

    // Assert
    ASSERT_EQ(12, mesh2d.num_faces);
}

TEST_F(CartesianApiTestFixture, CountHangingEdgesMesh2D_WithZeroHangingEdges_ShouldCountZeroEdges)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    int numHangingEdges;
    auto const errorCode = meshkernelapi::mkernel_mesh2d_count_hanging_edges(meshKernelId, numHangingEdges);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    ASSERT_EQ(0, numHangingEdges);
}

TEST_F(CartesianApiTestFixture, GetHangingEdgesMesh2D_WithOneHangingEdges_ShouldGetOneHangingEdges)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // delete an edge at the lower left corner to create a new hanging edge
    auto errorCode = meshkernelapi::mkernel_mesh2d_delete_edge(meshKernelId, 0.5, 0.0, 0.0, 0.0, 1.0, 1.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    int numHangingEdges;
    errorCode = meshkernelapi::mkernel_mesh2d_count_hanging_edges(meshKernelId, numHangingEdges);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_GT(numHangingEdges, 0);

    // Execute
    std::vector<int> hangingEdges(numHangingEdges);
    errorCode = meshkernelapi::mkernel_mesh2d_get_hanging_edges(meshKernelId, hangingEdges.data());
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    ASSERT_EQ(hangingEdges[0], 9);
}

TEST_F(CartesianApiTestFixture, DeleteHangingEdgesMesh2D_WithOneHangingEdges_ShouldDeleteOneHangingEdges)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // delete an edge at the lower left corner to create a new hanging edge
    auto errorCode = meshkernelapi::mkernel_mesh2d_delete_edge(meshKernelId, 0.5, 0.0, 0.0, 0.0, 1.0, 1.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Before deletion
    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(mesh2d.num_valid_edges, 16);

    // Execute
    errorCode = meshkernelapi::mkernel_mesh2d_delete_hanging_edges(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);

    // Assert
    ASSERT_EQ(mesh2d.num_valid_edges, 15);
}

TEST_F(CartesianApiTestFixture, DeleteEdgeByIndexThenDeleteHangingEdges)
{
    MakeMesh();
    int const meshKernelId = GetMeshKernelId();
    meshkernelapi::Mesh2D mesh2d{};
    int errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    int const initial_num_valid_edges = mesh2d.num_valid_edges;

    // delete horizontal edge of lower left corner
    errorCode = meshkernelapi::mkernel_mesh2d_delete_edge_by_index(meshKernelId, 0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(initial_num_valid_edges - 1, mesh2d.num_valid_edges);

    // Execute
    errorCode = meshkernelapi::mkernel_mesh2d_delete_hanging_edges(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(initial_num_valid_edges - 2, mesh2d.num_valid_edges);
}

TEST_F(CartesianApiTestFixture, GetNodesInPolygonMesh2D_OnMesh2D_ShouldGetAllNodes)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // By using an empty list, all nodes will be selected
    const meshkernelapi::GeometryList geometryListIn{};

    meshkernelapi::Mesh2D mesh2d{};
    auto errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    int numberOfNodes = -1;

    errorCode = mkernel_mesh2d_count_nodes_in_polygons(meshKernelId,
                                                       geometryListIn,
                                                       1,
                                                       numberOfNodes);

    ASSERT_EQ(numberOfNodes, mesh2d.num_nodes);

    // Execute
    std::vector<int> selectedNodes(mesh2d.num_nodes);

    errorCode = mkernel_mesh2d_get_nodes_in_polygons(meshKernelId,
                                                     geometryListIn,
                                                     1,
                                                     selectedNodes.data());

    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert (all nodes indices will be selected)
    std::vector actualResult(selectedNodes.data(), selectedNodes.data() + mesh2d.num_nodes);
    std::vector<int> expectedResult(mesh2d.num_nodes);
    std::iota(expectedResult.begin(), expectedResult.end(), 0);
    ASSERT_THAT(actualResult, ::testing::ContainerEq(expectedResult));
}

TEST_F(CartesianApiTestFixture, CountNodesInPolygonMesh2D_OnMesh2D_ShouldCountAllNodes)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // By using an empty list, all nodes will be selected
    const meshkernelapi::GeometryList geometryListIn{};

    // Execute
    int numNodes;
    const auto errorCode = mkernel_mesh2d_count_nodes_in_polygons(meshKernelId,
                                                                  geometryListIn,
                                                                  1,
                                                                  numNodes);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert all nodes have been selected
    ASSERT_EQ(12, numNodes);
}

TEST_F(CartesianApiTestFixture, CountSmallFlowEdges_OnMesh2D_ShouldCountSmallFlowEdges)
{
    // Prepare a mesh with two triangles
    meshkernelapi::Mesh2D mesh2d;

    std::vector node_x{0.0, 1.0, 1.0, 1.0};
    std::vector node_y{0.0, 0.0, 0.3, -0.3};
    std::vector edge_nodes{0, 3, 3, 1, 1, 0, 1, 2, 2, 0};

    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_nodes = edge_nodes.data();
    mesh2d.num_edges = static_cast<int>(edge_nodes.size() * 0.5);
    mesh2d.num_nodes = static_cast<int>(node_x.size());

    // Get the meshkernel id
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    auto errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    int numSmallFlowEdges;
    double const smallFlowEdgesThreshold = 100.0;
    errorCode = meshkernelapi::mkernel_mesh2d_count_small_flow_edge_centers(meshKernelId, smallFlowEdgesThreshold, numSmallFlowEdges);
    ASSERT_EQ(1, numSmallFlowEdges);
}

TEST_F(CartesianApiTestFixture, GetSmallFlowEdges_OnMesh2D_ShouldGetSmallFlowEdges)
{
    // Prepare a mesh with two triangles
    meshkernelapi::Mesh2D mesh2d;

    std::vector<double> node_x{0.0, 1.0, 1.0, 1.0};
    std::vector<double> node_y{0.0, 0.0, 0.3, -0.3};
    std::vector<int> edge_nodes{0, 3, 3, 1, 1, 0, 1, 2, 2, 0};

    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_nodes = edge_nodes.data();
    mesh2d.num_edges = static_cast<int>(edge_nodes.size() * 0.5);
    mesh2d.num_nodes = static_cast<int>(node_x.size());

    // Get the meshkernel id
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    auto errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    int numSmallFlowEdges;
    double const smallFlowEdgesThreshold = 100.0;
    errorCode = meshkernelapi::mkernel_mesh2d_count_small_flow_edge_centers(meshKernelId, smallFlowEdgesThreshold, numSmallFlowEdges);

    // Assert
    std::vector<double> coordinates_x(numSmallFlowEdges);
    std::vector<double> coordinates_y(numSmallFlowEdges);
    meshkernelapi::GeometryList result{};
    result.coordinates_x = coordinates_x.data();
    result.coordinates_y = coordinates_y.data();
    result.num_coordinates = numSmallFlowEdges;

    errorCode = mkernel_mesh2d_get_small_flow_edge_centers(meshKernelId, smallFlowEdgesThreshold, result);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    const double tolerance = 1e-6;
    ASSERT_NEAR(result.coordinates_x[0], 0.5, tolerance);
    ASSERT_NEAR(result.coordinates_y[0], 0.0, tolerance);
}

TEST_F(CartesianApiTestFixture, CountObtuseTriangles_OnMesh2DWithOneObtuseTriangle_ShouldCountObtuseTriangles)
{
    // Prepare a mesh with one obtuse triangle
    meshkernelapi::Mesh2D mesh2d;
    std::vector<double> coordinatesX{0.0, 3.0, -1.0, 1.5};
    std::vector<double> coordinatesY{0.0, 0.0, 2.0, -2.0};
    std::vector<int> edge_nodes{0, 1, 1, 2, 2, 0, 0, 3, 3, 1};
    mesh2d.node_x = coordinatesX.data();
    mesh2d.node_y = coordinatesY.data();
    mesh2d.edge_nodes = edge_nodes.data();
    mesh2d.num_edges = static_cast<int>(edge_nodes.size() * 0.5);
    mesh2d.num_nodes = static_cast<int>(coordinatesX.size());
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    auto errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    int numObtuseTriangles;
    errorCode = meshkernelapi::mkernel_mesh2d_count_obtuse_triangles(meshKernelId, numObtuseTriangles);

    // Assert
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(1, numObtuseTriangles);

    // Need to clear the obtuse triangle cache for the next tests
    meshkernelapi::GeometryList geometryList{};

    std::vector<double> coordinatesObtuseTrianglesX(numObtuseTriangles);
    std::vector<double> coordinatesObtuseTrianglesY(numObtuseTriangles);
    geometryList.coordinates_x = coordinatesObtuseTrianglesX.data();
    geometryList.coordinates_y = coordinatesObtuseTrianglesY.data();
    geometryList.num_coordinates = numObtuseTriangles;
    errorCode = mkernel_mesh2d_get_obtuse_triangles_mass_centers(meshKernelId, geometryList);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}

TEST_F(CartesianApiTestFixture, Mesh2DCountObtuseTriangles_OnMesh2DWithOneObtuseTriangle_ShouldGetObtuseTriangle)
{
    // Prepare a mesh with one obtuse triangle
    meshkernelapi::Mesh2D mesh2d;
    std::vector<double> coordinatesX{0.0, 3.0, -1.0, 1.5};
    std::vector<double> coordinatesY{0.0, 0.0, 2.0, -2.0};
    std::vector<int> edge_nodes{0, 1, 1, 2, 2, 0, 0, 3, 3, 1};
    mesh2d.node_x = coordinatesX.data();
    mesh2d.node_y = coordinatesY.data();
    mesh2d.edge_nodes = edge_nodes.data();
    mesh2d.num_edges = static_cast<int>(edge_nodes.size() * 0.5);
    mesh2d.num_nodes = static_cast<int>(coordinatesX.size());
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    auto errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    int numObtuseTriangles;
    errorCode = meshkernelapi::mkernel_mesh2d_count_obtuse_triangles(meshKernelId, numObtuseTriangles);
    meshkernelapi::GeometryList geometryList{};

    std::vector<double> coordinatesObtuseTrianglesX(numObtuseTriangles);
    std::vector<double> coordinatesObtuseTrianglesY(numObtuseTriangles);
    geometryList.coordinates_x = coordinatesObtuseTrianglesX.data();
    geometryList.coordinates_y = coordinatesObtuseTrianglesY.data();
    geometryList.num_coordinates = numObtuseTriangles;
    errorCode = mkernel_mesh2d_get_obtuse_triangles_mass_centers(meshKernelId, geometryList);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    ASSERT_EQ(1, numObtuseTriangles);
    const double tolerance = 1e-6;
    std::vector computedCoordinatesX(coordinatesObtuseTrianglesX.data(), coordinatesObtuseTrianglesX.data() + numObtuseTriangles);
    ASSERT_NEAR(computedCoordinatesX[0], 0.66666666666666652, tolerance);
    std::vector computedCoordinatesY(coordinatesObtuseTrianglesY.data(), coordinatesObtuseTrianglesY.data() + numObtuseTriangles);
    ASSERT_NEAR(computedCoordinatesY[0], 0.66666666666666652, tolerance);
}

TEST_F(CartesianApiTestFixture, Mesh2DDeleteSmallFlowEdgesAndSmallTriangles_OnMesh2DWithOneObtuseTriangle_ShouldDeleteOneEdge)
{
    // Prepare a mesh with one obtuse triangle
    meshkernelapi::Mesh2D mesh2d;
    std::vector coordinatesX{0.0, 3.0, -1.0, 1.5};
    std::vector coordinatesY{0.0, 0.0, 2.0, -2.0};
    std::vector edge_nodes{0, 1, 1, 2, 2, 0, 0, 3, 3, 1};

    mesh2d.node_x = coordinatesX.data();
    mesh2d.node_y = coordinatesY.data();
    mesh2d.edge_nodes = edge_nodes.data();
    mesh2d.num_nodes = static_cast<int>(coordinatesX.size());
    mesh2d.num_edges = static_cast<int>(edge_nodes.size() * 0.5);

    auto const meshKernelId = GetMeshKernelId();

    // Execute, with large length threshold
    auto errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    errorCode = meshkernelapi::mkernel_mesh2d_delete_small_flow_edges_and_small_triangles(meshKernelId, 1.0, 0.01);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    meshkernelapi::Mesh2D newMesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, newMesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // One edge is removed
    ASSERT_EQ(4, newMesh2d.num_valid_edges);
}

TEST_F(CartesianApiTestFixture, Mesh2dAveragingInterpolation_OnMesh2D_ShouldInterpolateValues)
{
    // Setup
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList samples{};
    std::vector<double> firstGridNodeCoordinateX{1.0, 2.0, 3.0, 1.0};
    std::vector<double> firstGridNodeCoordinateY{1.0, 3.0, 2.0, 4.0};
    std::vector<double> values{3.0, 10, 4.0, 5.0};

    samples.coordinates_x = firstGridNodeCoordinateX.data();
    samples.coordinates_y = firstGridNodeCoordinateY.data();
    samples.values = values.data();
    samples.num_coordinates = static_cast<int>(firstGridNodeCoordinateX.size());

    int const locationType = 1;             // Nodes
    int const averagingMethodType = 1;      // Simple averaging
    double const relativeSearchSize = 1.01; // The relative search size

    meshkernelapi::Mesh2D mesh2d;
    auto errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::GeometryList results{};
    std::vector<double> resultsCoordinateX(mesh2d.num_nodes);
    std::vector<double> resultsCoordinateY(mesh2d.num_nodes);
    std::vector<double> resultsValues(mesh2d.num_nodes);

    results.coordinates_x = resultsCoordinateX.data();
    results.coordinates_y = resultsCoordinateY.data();
    results.values = resultsValues.data();
    results.num_coordinates = static_cast<int>(resultsCoordinateX.size());

    // Execute
    errorCode = mkernel_mesh2d_averaging_interpolation(meshKernelId,
                                                       samples,
                                                       locationType,
                                                       averagingMethodType,
                                                       relativeSearchSize,
                                                       0,
                                                       results);

    // Assert the value has been interpolated
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    const double tolerance = 1e-6;
    std::vector computedResultsValues(resultsValues.data(), resultsValues.data() + mesh2d.num_nodes);
    ASSERT_NEAR(computedResultsValues[4], 3.0, tolerance);
}

TEST_F(CartesianApiTestFixture, Mesh2dTriangulationInterpolation_ShouldInterpolateValues)
{
    // Setup
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList samples{};
    std::vector<double> firstGridNodeCoordinateX{1.0, 2.0, 3.0, 1.0};
    std::vector<double> firstGridNodeCoordinateY{1.0, 3.0, 2.0, 4.0};
    std::vector<double> values{3.0, 10, 4.0, 5.0};

    samples.coordinates_x = firstGridNodeCoordinateX.data();
    samples.coordinates_y = firstGridNodeCoordinateY.data();
    samples.values = values.data();
    samples.num_coordinates = static_cast<int>(firstGridNodeCoordinateX.size());

    int const locationType = 1; // Nodes

    meshkernelapi::Mesh2D mesh2d;
    auto errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::GeometryList results{};
    std::vector<double> resultsCoordinateX(mesh2d.num_nodes);
    std::vector<double> resultsCoordinateY(mesh2d.num_nodes);
    std::vector<double> resultsValues(mesh2d.num_nodes);

    results.coordinates_x = resultsCoordinateX.data();
    results.coordinates_y = resultsCoordinateY.data();
    results.values = resultsValues.data();
    results.num_coordinates = static_cast<int>(resultsCoordinateX.size());

    // Execute
    errorCode = mkernel_mesh2d_triangulation_interpolation(meshKernelId,
                                                           samples,
                                                           locationType,
                                                           results);

    // Assert the value has been interpolated
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    const double tolerance = 1e-6;
    std::vector computedResultsValues(resultsValues.data(), resultsValues.data() + mesh2d.num_nodes);
    ASSERT_NEAR(computedResultsValues[8], 5.6666666666666670, tolerance);
}

TEST_F(CartesianApiTestFixture, Network1DComputeFixedChainages_ShouldGenerateMesh1D)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();

    double separator = meshkernelapi::mkernel_get_separator();
    std::vector<double> polyLineXCoordinate{0.0, 10.0, 20.0, separator, 10.0, 10.0, 10.0};
    std::vector<double> polyLineYCoordinate{0.0, 0.0, 0.0, separator, -10.0, 0.0, 10.0};

    meshkernelapi::GeometryList polylines{};
    polylines.coordinates_x = polyLineXCoordinate.data();
    polylines.coordinates_y = polyLineYCoordinate.data();
    polylines.num_coordinates = static_cast<int>(polyLineXCoordinate.size());
    polylines.geometry_separator = separator;

    // Set the network
    auto errorCode = mkernel_network1d_set(meshKernelId, polylines);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute, compute fixed chainages
    int const fixedChainagesSize = 3;
    double const minFaceSize = 0.01;
    double const fixedChainagesOffset = 10.0;
    std::vector<double> fixedChainages{5.0, separator, 5.0};
    errorCode = meshkernelapi::mkernel_network1d_compute_fixed_chainages(meshKernelId, fixedChainages.data(), fixedChainagesSize, minFaceSize, fixedChainagesOffset);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    // Convert network 1d to mesh1d
    errorCode = meshkernelapi::mkernel_network1d_to_mesh1d(meshKernelId, minFaceSize);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Asserts
    meshkernelapi::Mesh1D mesh1dResults;
    errorCode = mkernel_mesh1d_get_dimensions(meshKernelId, mesh1dResults);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(6, mesh1dResults.num_nodes);
    ASSERT_EQ(4, mesh1dResults.num_edges);

    auto edge_nodes = std::vector<int>(mesh1dResults.num_edges * 2);
    auto node_x = std::vector<double>(mesh1dResults.num_nodes);
    auto node_y = std::vector<double>(mesh1dResults.num_nodes);
    mesh1dResults.edge_nodes = edge_nodes.data();
    mesh1dResults.node_x = node_x.data();
    mesh1dResults.node_y = node_y.data();

    errorCode = mkernel_mesh1d_get_data(meshKernelId, mesh1dResults);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}

TEST_F(CartesianApiTestFixture, Network1DToMesh1d_FromPolylines_ShouldGenerateMesh1D)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();

    double separator = meshkernelapi::mkernel_get_separator();
    std::vector<double> polyLineXCoordinate{0.0, 10.0, 20.0, separator, 10.0, 10.0, 10.0};
    std::vector<double> polyLineYCoordinate{0.0, 0.0, 0.0, separator, -10.0, 0.0, 10.0};
    meshkernelapi::GeometryList polylines{};
    polylines.coordinates_x = polyLineXCoordinate.data();
    polylines.coordinates_y = polyLineYCoordinate.data();
    polylines.num_coordinates = static_cast<int>(polyLineXCoordinate.size());
    polylines.geometry_separator = separator;

    // Set the network
    auto errorCode = mkernel_network1d_set(meshKernelId, polylines);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute, compute offsetted chainages
    double const offset = 1.0;
    errorCode = meshkernelapi::mkernel_network1d_compute_offsetted_chainages(meshKernelId, offset);

    // Convert network 1d to mesh1d
    errorCode = meshkernelapi::mkernel_network1d_to_mesh1d(meshKernelId, 0.01);

    // Asserts
    meshkernelapi::Mesh1D mesh1dResults;
    errorCode = mkernel_mesh1d_get_dimensions(meshKernelId, mesh1dResults);
    ASSERT_EQ(41, mesh1dResults.num_valid_nodes);
    ASSERT_EQ(40, mesh1dResults.num_valid_edges);
}

TEST_F(CartesianApiTestFixture, ContactsComputeSingle_OnMesh2D_ShouldComputeContacts)
{
    auto [nodes_x, nodes_y, edges, face_nodes, num_face_nodes] = MakeMeshWithFaceNodesForApiTesting();
    const auto meshKernelId = GetMeshKernelId();

    meshkernelapi::Mesh2D mesh2d;
    mesh2d.node_x = &nodes_x[0];
    mesh2d.node_y = &nodes_y[0];

    mesh2d.edge_nodes = &edges[0];
    mesh2d.face_nodes = &face_nodes[0];
    mesh2d.nodes_per_face = &num_face_nodes[0];

    mesh2d.num_nodes = static_cast<int>(nodes_x.size());
    mesh2d.num_edges = static_cast<int>(edges.size() / 2);
    mesh2d.num_faces = static_cast<int>(num_face_nodes.size());

    auto errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Init 1d mesh
    meshkernelapi::Mesh1D mesh1d;
    std::vector<double> node_x{
        1.73493900000000,
        2.35659313023165,
        5.38347452702839,
        14.2980910429074,
        22.9324017677239,
        25.3723169493137,
        25.8072280000000};
    std::vector<double> node_y{
        -7.6626510000000,
        1.67281447902331,
        10.3513746546384,
        12.4797224193970,
        15.3007317677239,
        24.1623588554512,
        33.5111870000000};
    mesh1d.node_x = node_x.data();
    mesh1d.node_y = node_y.data();
    mesh1d.num_nodes = 7;

    std::vector<int> edge_nodes{0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6};
    mesh1d.edge_nodes = edge_nodes.data();
    mesh1d.num_edges = 6;

    errorCode = mkernel_mesh1d_set(meshKernelId, mesh1d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Init 1d mask
    std::vector<int> onedNodeMask{1, 1, 1, 1, 1, 1, 1};

    // Init polygon
    meshkernelapi::GeometryList polygon;
    polygon.geometry_separator = meshkernel::constants::missing::doubleValue;

    std::vector<double> xCoordinates{-30, 40, 40, -40, -30};
    std::vector<double> yCoordinates{-20, -20, 50, 50, -20};
    std::vector<double> zCoordinates{0, 0, 0, 0, 0};
    polygon.coordinates_x = xCoordinates.data();
    polygon.coordinates_y = yCoordinates.data();
    polygon.values = zCoordinates.data();
    polygon.num_coordinates = static_cast<int>(xCoordinates.size());

    // Execute
    errorCode = mkernel_contacts_compute_single(meshKernelId, onedNodeMask.data(), polygon, 0.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Get the new state
    meshkernelapi::Contacts contacts{};
    errorCode = mkernel_contacts_get_dimensions(meshKernelId, contacts);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<int> mesh1d_indices(contacts.num_contacts);
    std::vector<int> mesh2d_indices(contacts.num_contacts);
    contacts.mesh1d_indices = mesh1d_indices.data();
    contacts.mesh2d_indices = mesh2d_indices.data();

    errorCode = mkernel_contacts_get_data(meshKernelId, contacts);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    ASSERT_EQ(5, contacts.num_contacts);

    ASSERT_EQ(1, contacts.mesh1d_indices[0]);
    ASSERT_EQ(2, contacts.mesh1d_indices[1]);
    ASSERT_EQ(3, contacts.mesh1d_indices[2]);
    ASSERT_EQ(4, contacts.mesh1d_indices[3]);
    ASSERT_EQ(5, contacts.mesh1d_indices[4]);

    ASSERT_EQ(0, contacts.mesh2d_indices[0]);
    ASSERT_EQ(3, contacts.mesh2d_indices[1]);
    ASSERT_EQ(4, contacts.mesh2d_indices[2]);
    ASSERT_EQ(5, contacts.mesh2d_indices[3]);
    ASSERT_EQ(8, contacts.mesh2d_indices[4]);
}

TEST_F(CartesianApiTestFixture, GenerateTriangularGridThroughApi_OnClockWisePolygon_ShouldComputeValidTriangles)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometry_separator = meshkernel::constants::missing::doubleValue;
    std::vector xCoordinates{0.0, 0.0, 1.0, 1.0, 0.0};
    std::vector yCoordinates{0.0, 1.0, 1.0, 0.0, 0.0};
    std::vector zCoordinates(yCoordinates.size(), 0.0);

    geometryListIn.coordinates_x = xCoordinates.data();
    geometryListIn.coordinates_y = yCoordinates.data();
    geometryListIn.values = zCoordinates.data();
    geometryListIn.num_coordinates = static_cast<int>(xCoordinates.size());

    // Execute
    auto errorCode = mkernel_mesh2d_make_triangular_mesh_from_polygon(meshKernelId, geometryListIn);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    // Get the new state

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    ASSERT_EQ(5, mesh2d.num_nodes);
    ASSERT_EQ(8, mesh2d.num_edges);
}

TEST_F(CartesianApiTestFixture, TranslateMesh)
{
    // Prepare

    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::Mesh2D mesh2d{};
    [[maybe_unused]] int errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);

    // Get copy of original mesh
    // Will be used later to compare with translated mesh
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
    errorCode = mkernel_mesh2d_get_data(meshKernelId, mesh2d);

    double translationX = 10.0;
    double translationY = 5.0;

    meshkernelapi::mkernel_mesh2d_translate(meshKernelId, translationX, translationY);

    meshkernelapi::Mesh2D mesh2dTranslated{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2dTranslated);

    // Data for translated mesh
    std::vector<int> edge_faces_t(mesh2dTranslated.num_edges * 2);
    std::vector<int> edge_nodes_t(mesh2dTranslated.num_edges * 2);
    std::vector<int> face_nodes_t(mesh2dTranslated.num_face_nodes);
    std::vector<int> face_edges_t(mesh2dTranslated.num_face_nodes);
    std::vector<int> nodes_per_face_t(mesh2dTranslated.num_faces);
    std::vector<double> node_x_t(mesh2dTranslated.num_nodes);
    std::vector<double> node_y_t(mesh2dTranslated.num_nodes);
    std::vector<double> edge_x_t(mesh2dTranslated.num_edges);
    std::vector<double> edge_y_t(mesh2dTranslated.num_edges);
    std::vector<double> face_x_t(mesh2dTranslated.num_faces);
    std::vector<double> face_y_t(mesh2dTranslated.num_faces);

    mesh2dTranslated.edge_faces = edge_faces_t.data();
    mesh2dTranslated.edge_nodes = edge_nodes_t.data();
    mesh2dTranslated.face_nodes = face_nodes_t.data();
    mesh2dTranslated.face_edges = face_edges_t.data();
    mesh2dTranslated.nodes_per_face = nodes_per_face_t.data();
    mesh2dTranslated.node_x = node_x_t.data();
    mesh2dTranslated.node_y = node_y_t.data();
    mesh2dTranslated.edge_x = edge_x_t.data();
    mesh2dTranslated.edge_y = edge_y_t.data();
    mesh2dTranslated.face_x = face_x_t.data();
    mesh2dTranslated.face_y = face_y_t.data();

    errorCode = mkernel_mesh2d_get_data(meshKernelId, mesh2dTranslated);

    const double tolerance = 1e-12;

    for (int i = 0; i < mesh2d.num_nodes; ++i)
    {
        EXPECT_NEAR(mesh2d.node_x[i] + translationX, mesh2dTranslated.node_x[i], tolerance);
        EXPECT_NEAR(mesh2d.node_y[i] + translationY, mesh2dTranslated.node_y[i], tolerance);
    }
}

TEST_F(CartesianApiTestFixture, RotateMesh)
{
    // Prepare

    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::Mesh2D mesh2d{};
    [[maybe_unused]] int errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);

    // Get copy of original mesh
    // Will be used later to compare with translated mesh
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
    errorCode = mkernel_mesh2d_get_data(meshKernelId, mesh2d);

    double centreX = 10.0;
    double centreY = 5.0;
    double angle = 45.0;

    meshkernelapi::mkernel_mesh2d_rotate(meshKernelId, centreX, centreY, angle);

    meshkernelapi::Mesh2D mesh2dTranslated{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2dTranslated);

    // Data for translated mesh
    std::vector<int> edge_faces_t(mesh2dTranslated.num_edges * 2);
    std::vector<int> edge_nodes_t(mesh2dTranslated.num_edges * 2);
    std::vector<int> face_nodes_t(mesh2dTranslated.num_face_nodes);
    std::vector<int> face_edges_t(mesh2dTranslated.num_face_nodes);
    std::vector<int> nodes_per_face_t(mesh2dTranslated.num_faces);
    std::vector<double> node_x_t(mesh2dTranslated.num_nodes);
    std::vector<double> node_y_t(mesh2dTranslated.num_nodes);
    std::vector<double> edge_x_t(mesh2dTranslated.num_edges);
    std::vector<double> edge_y_t(mesh2dTranslated.num_edges);
    std::vector<double> face_x_t(mesh2dTranslated.num_faces);
    std::vector<double> face_y_t(mesh2dTranslated.num_faces);

    mesh2dTranslated.edge_faces = edge_faces_t.data();
    mesh2dTranslated.edge_nodes = edge_nodes_t.data();
    mesh2dTranslated.face_nodes = face_nodes_t.data();
    mesh2dTranslated.face_edges = face_edges_t.data();
    mesh2dTranslated.nodes_per_face = nodes_per_face_t.data();
    mesh2dTranslated.node_x = node_x_t.data();
    mesh2dTranslated.node_y = node_y_t.data();
    mesh2dTranslated.edge_x = edge_x_t.data();
    mesh2dTranslated.edge_y = edge_y_t.data();
    mesh2dTranslated.face_x = face_x_t.data();
    mesh2dTranslated.face_y = face_y_t.data();

    errorCode = mkernel_mesh2d_get_data(meshKernelId, mesh2dTranslated);

    const double tolerance = 1e-12;

    meshkernel::RigidBodyTransformation transformation;

    // First translate
    meshkernel::Translation translation(meshkernel::Vector(-centreX, -centreY));
    transformation.compose(translation);

    // Second rotate
    transformation.compose(meshkernel::Rotation(angle));

    // Third translate back
    translation.reset(meshkernel::Vector(centreX, centreY));
    transformation.compose(translation);

    for (int i = 0; i < mesh2d.num_nodes; ++i)
    {
        // Apply expected transformed to the original mesh node
        meshkernel::Point expected = transformation(meshkernel::Point(mesh2d.node_x[i], mesh2d.node_y[i]));

        // The compare the expected node with the transformed mesh node
        EXPECT_NEAR(expected.x, mesh2dTranslated.node_x[i], tolerance);
        EXPECT_NEAR(expected.y, mesh2dTranslated.node_y[i], tolerance);
    }
}

TEST_F(CartesianApiTestFixture, ConvertProjection)
{
    // Prepare

    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::Mesh2D mesh2d{};
    int errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Get copy of original mesh
    // Will be used later to compare with translated mesh
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
    errorCode = mkernel_mesh2d_get_data(meshKernelId, mesh2d);

    {
        std::string const zoneString("+proj=utm +lat_1=0.5 +lat_2=2 +n=0.5 +zone=31");

        // The mesh above has Cartesian projection, convert it to spherical
        int target_projection = 1;
        errorCode = meshkernelapi::mkernel_mesh2d_convert_projection(meshKernelId, target_projection, zoneString.c_str());
        ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

        int projection;
        errorCode = meshkernelapi::mkernel_get_projection(meshKernelId, projection);
        ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
        ASSERT_EQ(target_projection, projection);

        // Round trip back to Cartesian
        target_projection = 0;
        errorCode = meshkernelapi::mkernel_mesh2d_convert_projection(meshKernelId, target_projection, zoneString.c_str());
        ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
        errorCode = meshkernelapi::mkernel_get_projection(meshKernelId, projection);
        ASSERT_EQ(target_projection, projection);
    }

    meshkernelapi::Mesh2D mesh2dConvertedProjection{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2dConvertedProjection);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Data for translated mesh
    std::vector<int> edge_faces_t(mesh2dConvertedProjection.num_edges * 2);
    std::vector<int> edge_nodes_t(mesh2dConvertedProjection.num_edges * 2);
    std::vector<int> face_nodes_t(mesh2dConvertedProjection.num_face_nodes);
    std::vector<int> face_edges_t(mesh2dConvertedProjection.num_face_nodes);
    std::vector<int> nodes_per_face_t(mesh2dConvertedProjection.num_faces);
    std::vector<double> node_x_t(mesh2dConvertedProjection.num_nodes);
    std::vector<double> node_y_t(mesh2dConvertedProjection.num_nodes);
    std::vector<double> edge_x_t(mesh2dConvertedProjection.num_edges);
    std::vector<double> edge_y_t(mesh2dConvertedProjection.num_edges);
    std::vector<double> face_x_t(mesh2dConvertedProjection.num_faces);
    std::vector<double> face_y_t(mesh2dConvertedProjection.num_faces);

    mesh2dConvertedProjection.edge_faces = edge_faces_t.data();
    mesh2dConvertedProjection.edge_nodes = edge_nodes_t.data();
    mesh2dConvertedProjection.face_nodes = face_nodes_t.data();
    mesh2dConvertedProjection.face_edges = face_edges_t.data();
    mesh2dConvertedProjection.nodes_per_face = nodes_per_face_t.data();
    mesh2dConvertedProjection.node_x = node_x_t.data();
    mesh2dConvertedProjection.node_y = node_y_t.data();
    mesh2dConvertedProjection.edge_x = edge_x_t.data();
    mesh2dConvertedProjection.edge_y = edge_y_t.data();
    mesh2dConvertedProjection.face_x = face_x_t.data();
    mesh2dConvertedProjection.face_y = face_y_t.data();

    errorCode = mkernel_mesh2d_get_data(meshKernelId, mesh2dConvertedProjection);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    constexpr double tolerance = 1e-6;

    for (int i = 0; i < mesh2d.num_nodes; ++i)
    {
        EXPECT_NEAR(mesh2d.node_x[i], mesh2dConvertedProjection.node_x[i], tolerance);
        EXPECT_NEAR(mesh2d.node_y[i], mesh2dConvertedProjection.node_y[i], tolerance);
    }
}

TEST_F(CartesianApiTestFixture, InsertEdgeFromCoordinates_OnNonEmptyMesh_ShouldInsertNewEdge)
{
    // Prepare
    MakeMesh(4, 4, 2);
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    int firstNodeIndex;
    int secondNodeIndex;
    int edgeIndex;
    auto errorCode = meshkernelapi::mkernel_mesh2d_insert_edge_from_coordinates(meshKernelId,
                                                                                0,
                                                                                0,
                                                                                -1.0,
                                                                                -1.0,
                                                                                firstNodeIndex,
                                                                                secondNodeIndex,
                                                                                edgeIndex);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    ASSERT_EQ(40, edgeIndex);
    ASSERT_EQ(0, firstNodeIndex);
    ASSERT_EQ(25, secondNodeIndex);

    ASSERT_EQ(26, mesh2d.num_nodes);
    ASSERT_EQ(41, mesh2d.num_edges);
}

TEST_F(CartesianApiTestFixture, InsertEdgeFromCoordinates_OnMeshWithAnHole_ShouldInsertNewEdge)
{
    // Prepare
    MakeMesh(4, 4, 2);
    auto const meshKernelId = GetMeshKernelId();

    // Delete internal node
    int nodeIndex = 0;
    auto errorCode = meshkernelapi::mkernel_mesh2d_get_node_index(meshKernelId,
                                                                  4.0,
                                                                  4.0,
                                                                  0.1,
                                                                  0.0,
                                                                  0.0,
                                                                  10.0,
                                                                  10.0,
                                                                  nodeIndex);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    errorCode = meshkernelapi::mkernel_mesh2d_delete_node(meshKernelId, nodeIndex);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    int firstNodeIndex;
    int secondNodeIndex;
    int edgeIndex;
    errorCode = meshkernelapi::mkernel_mesh2d_insert_edge_from_coordinates(meshKernelId,
                                                                           4,
                                                                           4,
                                                                           3,
                                                                           3,
                                                                           firstNodeIndex,
                                                                           secondNodeIndex,
                                                                           edgeIndex);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    ASSERT_EQ(40, edgeIndex);
    ASSERT_EQ(25, firstNodeIndex);
    ASSERT_EQ(26, secondNodeIndex);

    ASSERT_EQ(27, mesh2d.num_nodes);
    ASSERT_EQ(41, mesh2d.num_edges);
}

class MeshLocationIndexTests : public ::testing::TestWithParam<std::tuple<meshkernel::Point, meshkernel::Location, int>>
{
public:
    [[nodiscard]] static std::vector<std::tuple<meshkernel::Point, meshkernel::Location, int>> GetData()
    {
        return {
            std::make_tuple<meshkernel::Point, meshkernel::Location, int>(meshkernel::Point{-1.0, -1.0}, meshkernel::Location::Nodes, 0),
            std::make_tuple<meshkernel::Point, meshkernel::Location, int>(meshkernel::Point{18.0, 18.0}, meshkernel::Location::Nodes, 10),
            std::make_tuple<meshkernel::Point, meshkernel::Location, int>(meshkernel::Point{5.0, -1.0}, meshkernel::Location::Edges, 12),
            std::make_tuple<meshkernel::Point, meshkernel::Location, int>(meshkernel::Point{11.0, 13.0}, meshkernel::Location::Edges, 5),
            std::make_tuple<meshkernel::Point, meshkernel::Location, int>(meshkernel::Point{0.5, 0.5}, meshkernel::Location::Faces, 0),
            std::make_tuple<meshkernel::Point, meshkernel::Location, int>(meshkernel::Point{18.0, 18.0}, meshkernel::Location::Faces, 4),
            std::make_tuple<meshkernel::Point, meshkernel::Location, int>(meshkernel::Point{7.0, 14.0}, meshkernel::Location::Faces, 3)};
    }
};

TEST_P(MeshLocationIndexTests, GetLocationIndex_OnAMesh_ShouldGetTheLocationIndex)
{
    // Prepare
    auto const& [point, location, expectedIndex] = GetParam();

    meshkernel::MakeGridParameters makeGridParameters;

    makeGridParameters.origin_x = 0.0;

    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = 10.0;
    makeGridParameters.block_size_y = 10.0;
    makeGridParameters.num_columns = 3;
    makeGridParameters.num_rows = 3;

    int meshKernelId = 0;
    const int projectionType = 0;
    auto errorCode = meshkernelapi::mkernel_allocate_state(projectionType, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_convert_to_mesh2d(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    int locationIndex = -1;
    const meshkernelapi::BoundingBox boundingBox;
    auto const locationInt = static_cast<int>(location);
    errorCode = mkernel_mesh2d_get_location_index(meshKernelId, point.x, point.y, locationInt, boundingBox, locationIndex);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    ASSERT_EQ(locationIndex, expectedIndex);
}
INSTANTIATE_TEST_SUITE_P(LocationIndexParametrizedTests, MeshLocationIndexTests, ::testing::ValuesIn(MeshLocationIndexTests::GetData()));
