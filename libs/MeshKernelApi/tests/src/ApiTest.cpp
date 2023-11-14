#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <MeshKernel/MeshTransformation.hpp>
#include <MeshKernel/Parameters.hpp>

#include <MeshKernelApi/CurvilinearGrid.hpp>
#include <MeshKernelApi/GeometryList.hpp>
#include <MeshKernelApi/Mesh1D.hpp>
#include <MeshKernelApi/Mesh2D.hpp>
#include <MeshKernelApi/MeshKernel.hpp>
#include <Version/Version.hpp>

#include <TestUtils/Definitions.hpp>
#include <TestUtils/MakeCurvilinearGrids.hpp>
#include <TestUtils/MakeMeshes.hpp>
#include <TestUtils/SampleFileReader.hpp>

#include <memory>
#include <numeric>

#include "CartesianApiTestFixture.hpp"

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
    ASSERT_EQ(11, mesh2d.num_nodes);
    ASSERT_EQ(15, mesh2d.num_edges);

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
    // Assert data
    const double tolerance = 1e-6;
    // Nodes
    ASSERT_NEAR(0.0, mesh2d.node_x[0], tolerance);
    ASSERT_NEAR(1.0, mesh2d.node_y[0], tolerance);
    // Edges
    ASSERT_EQ(0, mesh2d.edge_nodes[0]);
    ASSERT_EQ(3, mesh2d.edge_nodes[1]);
    ASSERT_NEAR(0.5, mesh2d.edge_x[0], tolerance);
    ASSERT_NEAR(1.0, mesh2d.edge_y[0], tolerance);
    // First face
    ASSERT_EQ(0, mesh2d.edge_faces[0]);
    ASSERT_EQ(-1, mesh2d.edge_faces[1]);
    ASSERT_EQ(0, mesh2d.edge_faces[2]);
    ASSERT_EQ(-1, mesh2d.edge_faces[3]);

    ASSERT_EQ(0, mesh2d.face_edges[0]);
    ASSERT_EQ(10, mesh2d.face_edges[1]);
    ASSERT_EQ(1, mesh2d.face_edges[2]);
    ASSERT_EQ(8, mesh2d.face_edges[3]);

    ASSERT_EQ(4, mesh2d.nodes_per_face[0]);
    ASSERT_NEAR(0.5, mesh2d.face_x[0], tolerance);
    ASSERT_NEAR(1.5, mesh2d.face_y[0], tolerance);
    // Second Face
    ASSERT_EQ(2, mesh2d.face_nodes[4]);
    ASSERT_EQ(5, mesh2d.face_nodes[5]);
    ASSERT_EQ(6, mesh2d.face_nodes[6]);
    ASSERT_EQ(3, mesh2d.face_nodes[7]);
    ASSERT_EQ(4, mesh2d.nodes_per_face[1]);

    ASSERT_NEAR(1.5, mesh2d.face_x[1], tolerance);
    ASSERT_NEAR(0.5, mesh2d.face_y[1], tolerance);
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

TEST_F(CartesianApiTestFixture, InsertEdgeThroughApi)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    int newEdgeIndex;
    auto errorCode = meshkernelapi::mkernel_mesh2d_insert_edge(meshKernelId, 0, 4, newEdgeIndex);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(17, newEdgeIndex);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);

    // Assert
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(12, mesh2d.num_nodes);
    ASSERT_EQ(18, mesh2d.num_edges);
}

TEST_F(CartesianApiTestFixture, MergeTwoNodesThroughApi)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    auto errorCode = meshkernelapi::mkernel_mesh2d_merge_two_nodes(meshKernelId, 0, 4);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);

    // Assert
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(11, mesh2d.num_nodes);
    ASSERT_EQ(15, mesh2d.num_edges);
}

TEST_F(CartesianApiTestFixture, MergeNodesThroughApi)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();
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
    auto errorCode = mkernel_mesh2d_merge_nodes(meshKernelId, geometry_list);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);

    // Assert (nothing is done, just check that the api communication works)
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(12, mesh2d.num_nodes);
    ASSERT_EQ(17, mesh2d.num_edges);
}

TEST_F(CartesianApiTestFixture, MergeNodesWithMergingDistanceThroughApi)
{
    // Prepare
    MakeMesh(2, 3, 1.0);
    auto const meshKernelId = GetMeshKernelId();
    meshkernelapi::GeometryList geometry_list{};

    // Execute
    auto errorCode = mkernel_mesh2d_merge_nodes_with_merging_distance(meshKernelId, geometry_list, 0.001);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);

    // Assert (nothing is done, just check that the api communication works)
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(12, mesh2d.num_nodes);
    ASSERT_EQ(17, mesh2d.num_edges);
}

TEST_F(CartesianApiTestFixture, OrthogonalizationThroughApi)
{
    // Set a new mesh in mesh
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Prepare
    meshkernel::OrthogonalizationParameters orthogonalizationParameters{};
    orthogonalizationParameters.outer_iterations = 1;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;

    meshkernelapi::GeometryList geometryList{};
    meshkernelapi::GeometryList landBoundaries{};

    // Execute
    auto errorCode = mkernel_mesh2d_initialize_orthogonalization(meshKernelId,
                                                                 1,
                                                                 orthogonalizationParameters,
                                                                 landBoundaries,
                                                                 geometryList);

    // Assert (nothing is done, just check that the api communication works)
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_prepare_outer_iteration_orthogonalization(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_compute_inner_ortogonalization_iteration(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_finalize_inner_ortogonalization_iteration(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_delete_orthogonalization(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);

    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(12, mesh2d.num_nodes);
    ASSERT_EQ(17, mesh2d.num_edges);
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
    ASSERT_EQ(42, mesh2d.num_nodes);
    ASSERT_EQ(107, mesh2d.num_edges);
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
    auto errorCode = meshkernelapi::mkernel_mesh2d_count_mesh_boundaries_as_polygons(meshKernelId, numberOfpolygonNodes);
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
    errorCode = mkernel_mesh2d_get_mesh_boundaries_as_polygons(meshKernelId, geometryListOut);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    const double tolerance = 1e-6;
    ASSERT_NEAR(0.0, geometryListOut.coordinates_x[0], tolerance);
    ASSERT_NEAR(0.0, geometryListOut.coordinates_y[0], tolerance);
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

TEST(ApiStatelessTests, Orthogonalize_OnInvaliMesh_ShouldThrowAMeshGeometryError)
{
    // Prepare
    int meshKernelId;
    int isGeographic = 0;
    meshkernelapi::mkernel_allocate_state(isGeographic, meshKernelId);

    auto [num_nodes, num_edges, node_x, node_y, node_type, edge_nodes, edge_type] =
        ReadLegacyMeshFile(TEST_FOLDER + "/data/InvalidMeshes/invalid_orthogonalization_net.nc");
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

    // Assert
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert there is a geometry error
    errorCode = meshkernelapi::mkernel_mesh2d_prepare_outer_iteration_orthogonalization(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::MeshGeometryErrorCode, errorCode);

    // Delete orthogonalization instance
    errorCode = meshkernelapi::mkernel_mesh2d_delete_orthogonalization(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Get the message
    auto exceptionMessage = std::make_unique<char[]>(512);
    errorCode = meshkernelapi::mkernel_get_error(exceptionMessage.get());
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Get the index of the invalid location
    int invalidIndex;
    int type;
    errorCode = meshkernelapi::mkernel_get_geometry_error(invalidIndex, type);
    ASSERT_EQ(static_cast<int>(meshkernel::Mesh::Location::Nodes), type);
    ASSERT_EQ(478, invalidIndex);
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

TEST_F(CartesianApiTestFixture, MakeCurvilinearGridFromTriangleThroughApi)
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

TEST_F(CartesianApiTestFixture, DeleteMesh2D_WithEmptyPolygon_ShouldDeleteMesh2D)
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
    ASSERT_EQ(0, mesh2d.num_nodes);
    ASSERT_EQ(0, mesh2d.num_edges);
}

TEST_F(CartesianApiTestFixture, GetDimensionsMesh1D_WithMesh1D_ShouldGetDimensionsMesh1D)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    std::vector nodesX{-16.1886410000000,
                       -16.1464995876014,
                       -16.1043581752028,
                       -16.0622167628042,
                       -15.7539488236928,
                       -6.86476658679268,
                       2.02441565010741,
                       10.9135970000000};

    std::vector nodesY{0.89018900000000,
                       9.78201442138723,
                       18.6738398427745,
                       27.5656652641617,
                       36.1966603330179,
                       36.4175095626911,
                       36.6383587923643,
                       36.8592080000000};

    std::vector edges{0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7};

    meshkernelapi::Mesh1D mesh1d;
    mesh1d.edge_nodes = edges.data();
    mesh1d.node_x = nodesX.data();
    mesh1d.node_y = nodesY.data();
    mesh1d.num_nodes = 8;
    mesh1d.num_edges = 7;

    auto errorCode = mkernel_mesh1d_set(GetMeshKernelId(), mesh1d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    meshkernelapi::Mesh1D mesh1dResults;
    errorCode = mkernel_mesh1d_get_dimensions(meshKernelId, mesh1dResults);

    // Assert
    ASSERT_EQ(8, mesh1dResults.num_nodes);
    ASSERT_EQ(7, mesh1dResults.num_edges);
}

TEST_F(CartesianApiTestFixture, GetDataMesh1D_WithMesh1D_ShouldGetDataMesh1D)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    std::vector nodesX{-16.1886410000000,
                       -16.1464995876014,
                       -16.1043581752028,
                       -16.0622167628042,
                       -15.7539488236928,
                       -6.86476658679268,
                       2.02441565010741,
                       10.9135970000000};

    std::vector nodesY{0.89018900000000,
                       9.78201442138723,
                       18.6738398427745,
                       27.5656652641617,
                       36.1966603330179,
                       36.4175095626911,
                       36.6383587923643,
                       36.8592080000000};

    std::vector edges{0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7};

    meshkernelapi::Mesh1D mesh1d;
    mesh1d.edge_nodes = edges.data();
    mesh1d.node_x = nodesX.data();
    mesh1d.node_y = nodesY.data();
    mesh1d.num_nodes = 8;
    mesh1d.num_edges = 7;

    auto errorCode = mkernel_mesh1d_set(GetMeshKernelId(), mesh1d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    meshkernelapi::Mesh1D mesh1dResults;
    errorCode = mkernel_mesh1d_get_dimensions(meshKernelId, mesh1dResults);
    std::vector<double> meshNodesX(mesh1dResults.num_nodes);
    std::vector<double> meshNodesY(mesh1dResults.num_nodes);
    std::vector<int> meshEdges(mesh1dResults.num_edges * 2);

    mesh1dResults.node_x = meshNodesX.data();
    mesh1dResults.node_y = meshNodesY.data();
    mesh1dResults.edge_nodes = meshEdges.data();
    errorCode = mkernel_mesh1d_get_data(meshKernelId, mesh1dResults);

    // Assert
    std::vector validMeshNodesX(nodesX.data(), nodesX.data() + mesh1d.num_nodes);
    std::vector computedMeshNodesX(meshNodesX.data(), meshNodesX.data() + mesh1dResults.num_nodes);
    ASSERT_THAT(computedMeshNodesX, ::testing::ContainerEq(validMeshNodesX));

    std::vector validMeshNodesY(nodesY.data(), nodesY.data() + mesh1d.num_nodes);
    std::vector computedMeshNodesY(meshNodesY.data(), meshNodesY.data() + mesh1dResults.num_nodes);
    ASSERT_THAT(computedMeshNodesY, ::testing::ContainerEq(validMeshNodesY));

    std::vector<double> validEdges(edges.data(), edges.data() + mesh1d.num_edges);
    std::vector<double> computedEdges(meshEdges.data(), meshEdges.data() + mesh1dResults.num_edges);
    ASSERT_THAT(computedEdges, ::testing::ContainerEq(validEdges));
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
    ASSERT_EQ(hangingEdges[0], 8);
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
    ASSERT_EQ(mesh2d.num_edges, 16);

    // Execute
    errorCode = meshkernelapi::mkernel_mesh2d_delete_hanging_edges(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);

    // Assert
    ASSERT_EQ(mesh2d.num_edges, 15);
}

TEST_F(CartesianApiTestFixture, ComputeOrthogonalizationMesh2D_WithOrthogonalMesh2D_ShouldOrthogonalize)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    meshkernel::OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor_at_boundary = 1.0;
    orthogonalizationParameters.areal_to_angle_smoothing_factor = 1.0;
    // By using an empty polygon the entire mesh will be orthogonalized
    meshkernelapi::GeometryList polygons{};
    // No land boundaries accounted
    meshkernelapi::GeometryList landBoundaries{};

    auto errorCode = mkernel_mesh2d_compute_orthogonalization(meshKernelId, 1, orthogonalizationParameters, polygons, landBoundaries);

    // Assert
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}

TEST_F(CartesianApiTestFixture, GetOrthogonalityMesh2D_OnMesh2D_ShouldGetOrthogonality)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::Mesh2D mesh2d{};
    auto errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::GeometryList edgeOrthogonality;
    std::vector<double> values(mesh2d.num_edges);
    edgeOrthogonality.values = values.data();
    edgeOrthogonality.num_coordinates = mesh2d.num_edges;

    // Execute
    errorCode = mkernel_mesh2d_get_orthogonality(meshKernelId, edgeOrthogonality);

    // Assert
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}

TEST_F(CartesianApiTestFixture, GetSmoothnessMesh2D_OnMesh2D_ShouldGetSmoothness)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::Mesh2D mesh2d{};
    auto errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::GeometryList edgeSmoothness;
    std::vector<double> values(mesh2d.num_edges);
    edgeSmoothness.values = values.data();
    edgeSmoothness.num_coordinates = mesh2d.num_edges;

    // Execute
    errorCode = mkernel_mesh2d_get_smoothness(meshKernelId, edgeSmoothness);

    // Assert
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
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

TEST_F(CartesianApiTestFixture, InsertNodeAndEdge_OnMesh2D_ShouldInsertNodeAndEdge)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    // Isolated nodes are removed by the administration done in mkernel_mesh2d_get_dimensions.
    // The newly inserted node should be connected to another one to form an edge.
    // In this manner, the edge will not be removed during the administration
    int newNodeIndex;
    auto errorCode = meshkernelapi::mkernel_mesh2d_insert_node(meshKernelId, -0.5, -0.5, newNodeIndex);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    int newEdgeIndex;
    errorCode = meshkernelapi::mkernel_mesh2d_insert_edge(meshKernelId, newNodeIndex, 0, newEdgeIndex);

    // Assert
    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(mesh2d.num_nodes, 13);
    ASSERT_EQ(mesh2d.num_edges, 18);
}

TEST_F(CartesianApiTestFixture, MoveNode_OnMesh2D_ShouldMoveNode)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    auto errorCode = meshkernelapi::mkernel_mesh2d_move_node(meshKernelId, -0.5, -0.5, 0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
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
    errorCode = mkernel_mesh2d_get_data(meshKernelId, mesh2d);

    // Assert
    ASSERT_EQ(mesh2d.node_x[0], -0.5);
    ASSERT_EQ(mesh2d.node_y[0], -0.5);
}

TEST_F(CartesianApiTestFixture, MoveNode_OnMesh2DWithInvalidIndex_ShouldReturnAnErrorCode)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    const auto errorCode = meshkernelapi::mkernel_mesh2d_move_node(meshKernelId, -0.5, -0.5, -1);

    // Assert call was not successful
    ASSERT_NE(meshkernel::ExitCode::Success, errorCode);
}

TEST_F(CartesianApiTestFixture, GetEdge_OnMesh2D_ShouldGetAnEdgeIndex)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    int edgeIndex;
    const auto errorCode = meshkernelapi::mkernel_mesh2d_get_edge(meshKernelId, 0.5, -0.5, 0.0, 0.0, 3.0, 3.0, edgeIndex);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    ASSERT_EQ(edgeIndex, 0);
}

TEST_F(CartesianApiTestFixture, GetNode_OnMesh2D_ShouldGetANodeIndex)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    int nodeIndex;
    const auto errorCode = meshkernelapi::mkernel_mesh2d_get_node_index(meshKernelId, 3.0, 3.0, 10.0, 0.0, 0.0, 3.0, 3.0, nodeIndex);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    ASSERT_EQ(nodeIndex, 11);
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
    ASSERT_EQ(4, newMesh2d.num_edges);
}

TEST_F(CartesianApiTestFixture, CurvilinearSetFrozenLinesOrthogonalize_ShouldSetFrozenLines)
{
    // Setup
    MakeRectangularCurvilinearGrid();
    auto const meshKernelId = GetMeshKernelId();
    meshkernel::OrthogonalizationParameters const orthogonalizationParameters{};

    auto errorCode = meshkernelapi::mkernel_curvilinear_initialize_orthogonalize(meshKernelId, orthogonalizationParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    errorCode = meshkernelapi::mkernel_curvilinear_set_frozen_lines_orthogonalize(meshKernelId, 20.0, 0.0, 20.0, 10.0);

    // Asset
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}

TEST_F(CartesianApiTestFixture, CurvilinearFinalizeOrthogonalize_ShouldFinalize)
{
    // Setup
    MakeRectangularCurvilinearGrid();
    auto const meshKernelId = GetMeshKernelId();
    meshkernel::OrthogonalizationParameters const orthogonalizationParameters{};

    auto errorCode = meshkernelapi::mkernel_curvilinear_initialize_orthogonalize(meshKernelId, orthogonalizationParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    errorCode = meshkernelapi::mkernel_curvilinear_finalize_orthogonalize(meshKernelId);

    // Assert
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
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
    // Convert network 1d to mesh1d
    errorCode = meshkernelapi::mkernel_network1d_to_mesh1d(meshKernelId, minFaceSize);

    // Asserts
    meshkernelapi::Mesh1D mesh1dResults;
    errorCode = mkernel_mesh1d_get_dimensions(meshKernelId, mesh1dResults);
    ASSERT_EQ(6, mesh1dResults.num_nodes);
    ASSERT_EQ(4, mesh1dResults.num_edges);
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
    ASSERT_EQ(41, mesh1dResults.num_nodes);
    ASSERT_EQ(40, mesh1dResults.num_edges);
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
    ASSERT_EQ(mesh2d.num_nodes, 8343);
    ASSERT_EQ(mesh2d.num_edges, 16502);
    ASSERT_EQ(mesh2d.num_faces, 8160);
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

TEST(MeshState, MKernelSnapSplineToLandBoundary_ShouldSnap)
{
    const double tolerance = 1e-6;

    // Setup
    int meshKernelId = 0;
    int setProjectionType = 0; // Cartesian
    auto errorCode = meshkernelapi::mkernel_allocate_state(setProjectionType, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // The land boundary to which the spline is to be snapped.
    std::vector<double> landBoundaryPointsX{257.002197, 518.753845, 938.006470};
    std::vector<double> landBoundaryPointsY{442.130066, 301.128662, 416.629822};

    // The original spline points.
    std::vector<double> splinePointsX{281.0023, 367.2529, 461.7534, 517.2538, 614.0045, 720.5051, 827.7558, 923.7563};
    std::vector<double> splinePointsY{447.3801, 401.6296, 354.3792, 318.3788, 338.629, 377.6294, 417.3798, 424.1299};

    // The expected spline values after snapping to land boundary.
    std::vector<double> expectedSplinePointsX{273.5868719643935, 359.5998304717778, 451.5303458337523, 517.7962262926076,
                                              616.7325138813335, 725.7358644094627, 836.2627853156330, 923.5001778441060};

    std::vector<double> expectedSplinePointsY{434.2730022174478, 386.1712239047134, 338.3551703843473, 306.3259738916997,
                                              327.9627689164845, 358.0902879743862, 388.6415116416172, 412.5818685325169};

    meshkernelapi::GeometryList landBoundaryGeometry{};
    landBoundaryGeometry.geometry_separator = meshkernel::constants::missing::doubleValue;
    landBoundaryGeometry.coordinates_x = landBoundaryPointsX.data();
    landBoundaryGeometry.coordinates_y = landBoundaryPointsY.data();
    landBoundaryGeometry.num_coordinates = static_cast<int>(landBoundaryPointsX.size());

    meshkernelapi::GeometryList splineGeometry{};
    splineGeometry.geometry_separator = meshkernel::constants::missing::doubleValue;
    splineGeometry.coordinates_x = splinePointsX.data();
    splineGeometry.coordinates_y = splinePointsY.data();
    splineGeometry.num_coordinates = static_cast<int>(splinePointsX.size());

    errorCode = mkernel_splines_snap_to_landboundary(meshKernelId, landBoundaryGeometry, splineGeometry, 0, static_cast<int>(splinePointsX.size() - 1));
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    for (size_t i = 0; i < splinePointsX.size(); ++i)
    {
        EXPECT_NEAR(splineGeometry.coordinates_x[i], expectedSplinePointsX[i], tolerance);
    }

    for (size_t i = 0; i < splinePointsX.size(); ++i)
    {
        EXPECT_NEAR(splineGeometry.coordinates_y[i], expectedSplinePointsY[i], tolerance);
    }
}

TEST(MeshState, MKernelSnapSplineToLandBoundary_ShouldThrowException)
{

    // Setup
    int meshKernelId = 0;
    int setProjectionType = 0; // Cartesian
    auto errorCode = meshkernelapi::mkernel_allocate_state(setProjectionType, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // The land boundary to which the spline is to be snapped.
    std::vector<double> landBoundaryPointsX{257.002197, 518.753845, 938.006470};
    std::vector<double> landBoundaryPointsY{442.130066, 301.128662, 416.629822};

    // The original spline points.
    std::vector<double> splinePointsX{281.0023, 367.2529, 461.7534, 517.2538, 614.0045, 720.5051, 827.7558, 923.7563};
    std::vector<double> splinePointsY{447.3801, 401.6296, 354.3792, 318.3788, 338.629, 377.6294, 417.3798, 424.1299};

    meshkernelapi::GeometryList landBoundaryGeometry{};
    landBoundaryGeometry.geometry_separator = meshkernel::constants::missing::doubleValue;

    meshkernelapi::GeometryList splineGeometry{};
    splineGeometry.geometry_separator = meshkernel::constants::missing::doubleValue;

    //--------------------------------
    // Start index is less than 0
    errorCode = mkernel_splines_snap_to_landboundary(meshKernelId,
                                                     landBoundaryGeometry,
                                                     splineGeometry,
                                                     -2,
                                                     1);
    EXPECT_EQ(meshkernel::ExitCode::ConstraintErrorCode, errorCode);

    //--------------------------------
    // Start index is greater than end index
    errorCode = mkernel_splines_snap_to_landboundary(meshKernelId,
                                                     landBoundaryGeometry,
                                                     splineGeometry,
                                                     2,
                                                     1);
    EXPECT_EQ(meshkernel::ExitCode::ConstraintErrorCode, errorCode);

    //--------------------------------
    // The land boundary is not set
    errorCode = mkernel_splines_snap_to_landboundary(meshKernelId,
                                                     landBoundaryGeometry,
                                                     splineGeometry,
                                                     0,
                                                     static_cast<int>(splinePointsX.size() - 1));
    EXPECT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    // First define the number of land boundary points
    landBoundaryGeometry.num_coordinates = static_cast<int>(landBoundaryPointsX.size());

    // The land boundary points are null
    errorCode = mkernel_splines_snap_to_landboundary(meshKernelId,
                                                     landBoundaryGeometry,
                                                     splineGeometry, 0,
                                                     static_cast<int>(splinePointsX.size() - 1));
    EXPECT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    //--------------------------------
    // Now define the land boundary
    landBoundaryGeometry.coordinates_x = landBoundaryPointsX.data();
    landBoundaryGeometry.coordinates_y = landBoundaryPointsY.data();

    // The number of spline points is 0
    errorCode = mkernel_splines_snap_to_landboundary(meshKernelId,
                                                     landBoundaryGeometry,
                                                     splineGeometry,
                                                     0,
                                                     static_cast<int>(splinePointsX.size() - 1));
    EXPECT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    // define the number of spline points
    splineGeometry.num_coordinates = static_cast<int>(splinePointsX.size());

    // The spline values are null
    errorCode = mkernel_splines_snap_to_landboundary(meshKernelId,
                                                     landBoundaryGeometry,
                                                     splineGeometry,
                                                     0,
                                                     static_cast<int>(splinePointsX.size() - 1));
    EXPECT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    splineGeometry.coordinates_x = splinePointsX.data();
    splineGeometry.coordinates_y = splinePointsY.data();

    // Start spline index is greater than the number of spline points
    errorCode = mkernel_splines_snap_to_landboundary(meshKernelId,
                                                     landBoundaryGeometry,
                                                     splineGeometry,
                                                     static_cast<int>(splinePointsX.size()) + 1,
                                                     static_cast<int>(splinePointsX.size()) + 2);
    EXPECT_EQ(meshkernel::ExitCode::ConstraintErrorCode, errorCode);

    // End spline index is greater than the number of spline points
    errorCode = mkernel_splines_snap_to_landboundary(meshKernelId,
                                                     landBoundaryGeometry,
                                                     splineGeometry,
                                                     0,
                                                     static_cast<int>(splinePointsX.size()));
    EXPECT_EQ(meshkernel::ExitCode::ConstraintErrorCode, errorCode);
}

TEST(MeshState, PolygonSnapToLandboundary_ShouldSnapPolygonToLandBoundary)
{
    const double tolerance = 1e-6;

    // Setup
    int meshKernelId = 0;
    int setProjectionType = 0; // Cartesian
    auto errorCode = meshkernelapi::mkernel_allocate_state(setProjectionType, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> landBoundaryPointsX{139.251465, 527.753906, 580.254211, 194.001801};
    std::vector<double> landBoundaryPointsY{497.630615, 499.880676, 265.878296, 212.627762};

    std::vector<double> polygonPointsX{170.001648, 263.002228, 344.002747,
                                       458.753448, 515.753845, 524.753906,
                                       510.503754, 557.754089, 545.004028,
                                       446.003387, 340.252716, 242.752106,
                                       170.001648};
    std::vector<double> polygonPointsY{472.880371, 472.880371, 475.130432,
                                       482.630493, 487.130554, 434.630005,
                                       367.129333, 297.378601, 270.378357,
                                       259.128235, 244.128067, 226.877884,
                                       472.880371};

    // The expected polygon values after snapping to land boundary.
    std::vector<double> expectedSnappedPointX = {169.8572772242283, 262.8547378163090, 343.8655709877979,
                                                 458.6558591358565, 515.6804060372598, 541.5480568270806,
                                                 555.2836667233159, 572.4472626165707, 546.2703464583593,
                                                 447.5942143903486, 341.7865993173012, 243.7707524316129,
                                                 169.8572772242283};

    std::vector<double> expectedSnappedPointY = {497.8078724305628, 498.3464789799546, 498.8156634613377,
                                                 499.4804859264834, 499.8107507986815, 438.3979070214996,
                                                 377.1760644727631, 300.6751319852315, 261.1931241088368,
                                                 247.5891786326750, 233.0020541046851, 219.4891385810638,
                                                 497.8078724305628};

    meshkernelapi::GeometryList landBoundaryGeometry{};
    landBoundaryGeometry.geometry_separator = meshkernel::constants::missing::doubleValue;
    landBoundaryGeometry.coordinates_x = landBoundaryPointsX.data();
    landBoundaryGeometry.coordinates_y = landBoundaryPointsY.data();
    landBoundaryGeometry.num_coordinates = static_cast<int>(landBoundaryPointsX.size());

    meshkernelapi::GeometryList polygonGeometry{};
    polygonGeometry.geometry_separator = meshkernel::constants::missing::doubleValue;
    polygonGeometry.coordinates_x = polygonPointsX.data();
    polygonGeometry.coordinates_y = polygonPointsY.data();
    polygonGeometry.num_coordinates = static_cast<int>(polygonPointsX.size());

    errorCode = meshkernelapi::mkernel_polygon_snap_to_landboundary(meshKernelId,
                                                                    landBoundaryGeometry,
                                                                    polygonGeometry,
                                                                    0,
                                                                    static_cast<int>(polygonPointsX.size()) - 1);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    for (size_t i = 0; i < polygonPointsX.size(); ++i)
    {
        EXPECT_NEAR(polygonGeometry.coordinates_x[i], expectedSnappedPointX[i], tolerance);
    }

    for (size_t i = 0; i < polygonPointsX.size(); ++i)
    {
        EXPECT_NEAR(polygonGeometry.coordinates_y[i], expectedSnappedPointY[i], tolerance);
    }
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

TEST(Mesh1D, Mesh1DSetAndAdd)
{
    using namespace meshkernelapi;

    // Init 1d mesh
    meshkernelapi::Mesh1D mesh1d_1;
    meshkernelapi::Mesh1D mesh1d_2;
    {
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

        std::vector edge_nodes{0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6};

        mesh1d_1.node_x = node_x.data();
        mesh1d_1.node_y = node_y.data();
        mesh1d_1.num_nodes = static_cast<int>(node_x.size());
        mesh1d_1.edge_nodes = edge_nodes.data();
        mesh1d_1.num_edges = static_cast<int>(edge_nodes.size()) / 2;

        double const offset = node_x.back() + 1.0;
        std::transform(node_x.begin(),
                       node_x.end(),
                       node_x.begin(),
                       [offset](double const val)
                       { return val + offset; });
        mesh1d_2.node_x = node_x.data();
        mesh1d_2.node_y = node_y.data();
        mesh1d_2.num_nodes = static_cast<int>(node_x.size());
        mesh1d_2.edge_nodes = edge_nodes.data();
        mesh1d_2.num_edges = static_cast<int>(edge_nodes.size()) / 2;
    }

    // allocate state
    int mk_id = 0;
    int errorCode = mkernel_allocate_state(0, mk_id);

    // first initialise using the first mesh, mesh1d_1
    errorCode = mkernel_mesh1d_set(mk_id, mesh1d_1);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // then add the second mesh, mesh1d_2
    errorCode = mkernel_mesh1d_add(mk_id, mesh1d_2);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // get the dimensions and data of the resulting mesh
    Mesh1D mesh1d;
    errorCode = mkernel_mesh1d_get_dimensions(mk_id, mesh1d);
    std::vector<double> meshNodesX(mesh1d.num_nodes);
    std::vector<double> meshNodesY(mesh1d.num_nodes);
    std::vector<int> meshEdges(mesh1d.num_edges * 2);

    mesh1d.node_x = meshNodesX.data();
    mesh1d.node_y = meshNodesY.data();
    mesh1d.edge_nodes = meshEdges.data();
    errorCode = mkernel_mesh1d_get_data(mk_id, mesh1d);

    EXPECT_EQ(mesh1d.num_nodes, mesh1d_1.num_nodes + mesh1d_2.num_nodes);
    EXPECT_EQ(mesh1d.num_edges, mesh1d_1.num_edges + mesh1d_1.num_edges);

    for (int i = 0; i < mesh1d_1.num_nodes; ++i)
    {
        EXPECT_EQ(mesh1d.node_x[i], mesh1d_1.node_x[i]);
        EXPECT_EQ(mesh1d.node_y[i], mesh1d_1.node_y[i]);
    }

    for (int i = mesh1d_1.num_nodes; i < mesh1d_1.num_nodes + mesh1d_2.num_nodes; ++i)
    {
        EXPECT_EQ(mesh1d.node_x[i], mesh1d_2.node_x[i - mesh1d_1.num_nodes]);
        EXPECT_EQ(mesh1d.node_y[i], mesh1d_2.node_y[i - mesh1d_1.num_nodes]);
    }

    errorCode = mkernel_deallocate_state(mk_id);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}

TEST(Mesh2D, Mesh2DSetAndAdd)
{
    using namespace meshkernelapi;

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
    int errorCode = mkernel_allocate_state(0, mk_id);

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

    EXPECT_EQ(mesh2d.num_nodes, num_nodes_1 + num_nodes_2);
    EXPECT_EQ(mesh2d.num_edges, num_edges_1 + num_edges_2);

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