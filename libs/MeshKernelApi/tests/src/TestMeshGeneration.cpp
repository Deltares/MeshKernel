#include "TestMeshGeneration.hpp"

#include "MeshKernelApi/MeshKernel.hpp"

int GenerateCurvilinearMesh(const int meshKernelId,
                            const int nodesX, const int nodesY,
                            const double deltaX, const double deltaY,
                            const double originX, const double originY)
{
    meshkernel::MakeGridParameters makeGridParameters;

    // num_columns and num_rows indicate number of elements in each direction, so value = nodes - 1
    makeGridParameters.num_columns = nodesX - 1;
    makeGridParameters.num_rows = nodesY - 1;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = originX;
    makeGridParameters.origin_y = originY;
    makeGridParameters.block_size_x = deltaX;
    makeGridParameters.block_size_y = deltaY;

    // Generate curvilinear grid.
    int errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    return errorCode;
}

int GenerateUnstructuredMesh(int meshKernelId, meshkernel::UInt numRows, meshkernel::UInt numColumns, double delta)
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
