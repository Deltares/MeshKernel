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

// Create mesh and return the meshkernel id
int MakeUnstructuredMesh(int meshKernelId, meshkernel::UInt numRows = 2, meshkernel::UInt numColumns = 3, double delta = 1.0)
{

    // Set-up new mesh
    auto [num_nodes, num_edges, node_x, node_y, edge_nodes] = MakeRectangularMeshForApiTesting(numRows, numColumns, delta);
    meshkernelapi::Mesh2D mesh2d{};
    mesh2d.num_edges = static_cast<int>(num_edges);
    mesh2d.num_nodes = static_cast<int>(num_nodes);
    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_nodes = edge_nodes.data();
    int errorCode = meshkernelapi::mkernel_mesh2d_set(meshKernelId, mesh2d);

    return errorCode;
}

TEST(ApiCacheTest, GetHangingEdgesMesh2D_WithOneHangingEdges_GetOneHangingEdgesFailures)
{
    int meshKernelId = 0;

    int errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Prepare
    errorCode = MakeUnstructuredMesh(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // delete an edge at the lower left corner to create a new hanging edge
    errorCode = meshkernelapi::mkernel_mesh2d_delete_edge(meshKernelId, 0.5, 0.0, 0.0, 0.0, 1.0, 1.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    int numHangingEdges;
    errorCode = meshkernelapi::mkernel_mesh2d_count_hanging_edges(meshKernelId, numHangingEdges);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Already cached
    errorCode = meshkernelapi::mkernel_mesh2d_count_hanging_edges(meshKernelId, numHangingEdges);
    ASSERT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    // Re-cache
    errorCode = meshkernelapi::mkernel_mesh2d_count_hanging_edges(meshKernelId, numHangingEdges);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<int> hangingEdges(numHangingEdges);
    errorCode = meshkernelapi::mkernel_mesh2d_get_hanging_edges(meshKernelId, hangingEdges.data());
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_get_hanging_edges(meshKernelId, hangingEdges.data());
    ASSERT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    errorCode = meshkernelapi::mkernel_expunge_state(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}

TEST(ApiCacheTest, GetNodesInPolygonMesh2D_OnMesh2D_NodeInPolygonFailures)
{
    int meshKernelId = 0;

    int errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = MakeUnstructuredMesh(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // By using an empty list, all nodes will be selected
    const meshkernelapi::GeometryList geometryListIn{};

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    int numberOfNodes = -1;
    // Execute
    std::vector<int> selectedNodes(mesh2d.num_nodes);

    // No daata has been cached.
    errorCode = mkernel_mesh2d_get_nodes_in_polygons(meshKernelId, geometryListIn, 0, selectedNodes.data());
    ASSERT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    errorCode = mkernel_mesh2d_count_nodes_in_polygons(meshKernelId, geometryListIn, 1, numberOfNodes);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(numberOfNodes, mesh2d.num_nodes);

    // Values have been cached already
    errorCode = mkernel_mesh2d_count_nodes_in_polygons(meshKernelId, geometryListIn, 1, numberOfNodes);
    ASSERT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    // Re-cache data
    errorCode = mkernel_mesh2d_count_nodes_in_polygons(meshKernelId, geometryListIn, 1, numberOfNodes);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Incorrect parameters, should be the same as the call to count (which does the cahcing)
    errorCode = mkernel_mesh2d_get_nodes_in_polygons(meshKernelId, geometryListIn, 0, selectedNodes.data());
    ASSERT_EQ(meshkernel::ExitCode::ConstraintErrorCode, errorCode);

    // Re-cache data
    errorCode = mkernel_mesh2d_count_nodes_in_polygons(meshKernelId, geometryListIn, 1, numberOfNodes);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Incorrect parameters, should be the same as the call to count (which does the cahcing)
    int* nullArray = nullptr;
    errorCode = mkernel_mesh2d_get_nodes_in_polygons(meshKernelId, geometryListIn, 1, nullArray);
    ASSERT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    errorCode = meshkernelapi::mkernel_expunge_state(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}

TEST(ApiCacheTest, GetSmallFlowEdges_OnMesh2D_GetSmallFlowEdgesFailures)
{
    // Prepare a mesh with two triangles
    meshkernelapi::Mesh2D mesh2d;
    int meshKernelId = 0;

    int errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = MakeUnstructuredMesh(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> node_x{0.0, 1.0, 1.0, 1.0};
    std::vector<double> node_y{0.0, 0.0, 0.3, -0.3};
    std::vector<int> edge_nodes{0, 3, 3, 1, 1, 0, 1, 2, 2, 0};

    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_nodes = edge_nodes.data();
    mesh2d.num_edges = static_cast<int>(edge_nodes.size() * 0.5);
    mesh2d.num_nodes = static_cast<int>(node_x.size());

    const double smallFlowEdgesThreshold = 100.0;
    meshkernelapi::GeometryList result{};

    errorCode = mkernel_mesh2d_get_small_flow_edge_centers(meshKernelId, smallFlowEdgesThreshold, result);
    ASSERT_EQ(meshkernel::ExitCode::ConstraintErrorCode, errorCode);

    // Execute
    errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    int numSmallFlowEdges;
    errorCode = meshkernelapi::mkernel_mesh2d_count_small_flow_edge_centers(meshKernelId, smallFlowEdgesThreshold, numSmallFlowEdges);

    // Assert
    std::vector<double> coordinates_x(numSmallFlowEdges);
    std::vector<double> coordinates_y(numSmallFlowEdges);
    result.coordinates_x = coordinates_x.data();
    result.coordinates_y = coordinates_y.data();
    result.num_coordinates = numSmallFlowEdges;

    errorCode = mkernel_mesh2d_get_small_flow_edge_centers(meshKernelId, smallFlowEdgesThreshold, result);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Cache has been deleted
    errorCode = mkernel_mesh2d_get_small_flow_edge_centers(meshKernelId, smallFlowEdgesThreshold, result);
    ASSERT_EQ(meshkernel::ExitCode::ConstraintErrorCode, errorCode);

    // Re-cache values
    errorCode = meshkernelapi::mkernel_mesh2d_count_small_flow_edge_centers(meshKernelId, smallFlowEdgesThreshold, numSmallFlowEdges);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Get values with different set of options
    errorCode = mkernel_mesh2d_get_small_flow_edge_centers(meshKernelId, 2.0 * smallFlowEdgesThreshold, result);
    ASSERT_EQ(meshkernel::ExitCode::ConstraintErrorCode, errorCode);

    // Re-cache values
    errorCode = meshkernelapi::mkernel_mesh2d_count_small_flow_edge_centers(meshKernelId, smallFlowEdgesThreshold, numSmallFlowEdges);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    // Attempt at getting the size again will cause an error
    errorCode = meshkernelapi::mkernel_mesh2d_count_small_flow_edge_centers(meshKernelId, smallFlowEdgesThreshold, numSmallFlowEdges);
    ASSERT_EQ(meshkernel::ExitCode::ConstraintErrorCode, errorCode);

    errorCode = meshkernelapi::mkernel_expunge_state(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}

TEST(ApiCacheTest, CountObtuseTriangles_OnMesh2DWithOneObtuseTriangle_ObtuseTrianglesFailures)
{
    // Prepare a mesh with one obtuse triangle
    int meshKernelId = 0;

    int errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d;
    std::vector<double> coordinatesX{0.0, 3.0, -1.0, 1.5};
    std::vector<double> coordinatesY{0.0, 0.0, 2.0, -2.0};
    std::vector<int> edge_nodes{0, 1, 1, 2, 2, 0, 0, 3, 3, 1};
    mesh2d.node_x = coordinatesX.data();
    mesh2d.node_y = coordinatesY.data();
    mesh2d.edge_nodes = edge_nodes.data();
    mesh2d.num_edges = static_cast<int>(edge_nodes.size() * 0.5);
    mesh2d.num_nodes = static_cast<int>(coordinatesX.size());

    // Execute
    errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Need to clear the obtuse triangle cache for the next tests
    meshkernelapi::GeometryList geometryList{};

    // Data has not yet been cached
    errorCode = mkernel_mesh2d_get_obtuse_triangles_mass_centers(meshKernelId, geometryList);
    ASSERT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    int numObtuseTriangles;
    errorCode = meshkernelapi::mkernel_mesh2d_count_obtuse_triangles(meshKernelId, numObtuseTriangles);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Already cached, cached data will be deleted
    errorCode = meshkernelapi::mkernel_mesh2d_count_obtuse_triangles(meshKernelId, numObtuseTriangles);
    ASSERT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    // Re-cache data.
    errorCode = meshkernelapi::mkernel_mesh2d_count_obtuse_triangles(meshKernelId, numObtuseTriangles);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> coordinatesObtuseTrianglesX(numObtuseTriangles);
    std::vector<double> coordinatesObtuseTrianglesY(numObtuseTriangles);
    geometryList.coordinates_x = coordinatesObtuseTrianglesX.data();
    geometryList.coordinates_y = coordinatesObtuseTrianglesY.data();
    geometryList.num_coordinates = numObtuseTriangles;
    errorCode = mkernel_mesh2d_get_obtuse_triangles_mass_centers(meshKernelId, geometryList);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Cache has been deleted in the last call
    errorCode = mkernel_mesh2d_get_obtuse_triangles_mass_centers(meshKernelId, geometryList);
    ASSERT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    errorCode = meshkernelapi::mkernel_expunge_state(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}
