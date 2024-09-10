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

    int geometryListDimension = -1;
    errorCode = meshkernelapi::mkernel_mesh2d_get_property_dimension(meshKernelId, 0, geometryListDimension);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    meshkernelapi::GeometryList propertyvalues{};
    propertyvalues.num_coordinates = geometryListDimension;
    propertyvalues.geometry_separator = meshkernel::constants::missing::doubleValue;
    std::vector<double> values(geometryListDimension);
    propertyvalues.values = values.data();
    errorCode = mkernel_mesh2d_get_property(meshKernelId, 0, propertyvalues);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    EXPECT_EQ(propertyvalues.num_coordinates, 12);
    const double tolerance = 1e-4;

    EXPECT_NEAR(values[0], 0.055751274056612614, tolerance);
    EXPECT_NEAR(values[1], 0.056220640190527582, tolerance);
    EXPECT_NEAR(values[2], 0.051193798544321531, tolerance);
    EXPECT_NEAR(values[3], 0.056591641726992326, tolerance);
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
