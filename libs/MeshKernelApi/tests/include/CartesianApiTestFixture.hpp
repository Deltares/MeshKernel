#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <MeshKernel/Parameters.hpp>

#include <MeshKernelApi/Mesh2D.hpp>
#include <MeshKernelApi/MeshKernel.hpp>

#include <TestUtils/MakeMeshes.hpp>

class CartesianApiTestFixture : public testing::Test
{
public:
    /// Constructor for allocating state
    CartesianApiTestFixture()
    {
        int isGeographic = 0;
        const auto errorCode = meshkernelapi::mkernel_allocate_state(isGeographic, m_meshKernelId);
        if (errorCode != 0)
        {
            throw std::runtime_error("Could not allocate state");
        }
    }

    /// Destructor for deallocating state
    ~CartesianApiTestFixture()
    {
        meshkernelapi::mkernel_deallocate_state(m_meshKernelId);
    }

    /// @brief Make a mesh
    /// @param[in]  numRows            Number of rows
    /// @param[in]  numColumns            Number of columns
    /// @param[in]  delta        Distance between neighboring nodes
    void MakeMesh(meshkernel::UInt numRows = 2, meshkernel::UInt numColumns = 3, double delta = 1.0)
    {
        // Set-up new mesh
        auto [num_nodes, num_edges, node_x, node_y, edge_nodes] = MakeRectangularMeshForApiTesting(numRows, numColumns, delta);
        meshkernelapi::Mesh2D mesh2d{};
        mesh2d.num_edges = static_cast<int>(num_edges);
        mesh2d.num_nodes = static_cast<int>(num_nodes);
        mesh2d.node_x = node_x.data();
        mesh2d.node_y = node_y.data();
        mesh2d.edge_nodes = edge_nodes.data();
        const auto errorCode = mkernel_mesh2d_set(m_meshKernelId, mesh2d);
        if (errorCode != 0)
        {
            throw std::runtime_error("Could not set mesh2d");
        }
    }

    void MakeRectangularCurvilinearGrid(meshkernel::UInt numberOfColumns = 4,
                                        meshkernel::UInt numberOfRows = 4,
                                        double blockSizeX = 10.0,
                                        double blockSizeY = 10.0,
                                        double originX = 0.0,
                                        double originY = 0.0) const
    {
        meshkernel::MakeGridParameters makeGridParameters{};
        makeGridParameters.num_columns = static_cast<int>(numberOfColumns);
        makeGridParameters.num_rows = static_cast<int>(numberOfRows);
        makeGridParameters.angle = 0.0;
        makeGridParameters.origin_x = originX;
        makeGridParameters.origin_y = originY;
        makeGridParameters.block_size_x = blockSizeX;
        makeGridParameters.block_size_y = blockSizeY;

        auto const errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(m_meshKernelId, makeGridParameters);
        if (errorCode != 0)
        {
            throw std::runtime_error("Could not create rectangular curvilinear grid");
        }
    }

    [[nodiscard]] int GetMeshKernelId() const
    {
        return m_meshKernelId;
    }

private:
    int m_meshKernelId{};
};

static auto GebcoMakeGridParameters()
{

    double lonMin = -1;
    double lonMax = -0.2;
    double latMin = 49.1;
    double latMax = 49.6;
    double lonRes = 0.1;
    double latRes = 0.1;
    int numX = static_cast<int>(std::ceil((lonMax - lonMin) / lonRes));
    int numY = static_cast<int>(std::ceil((latMax - latMin) / latRes));

    meshkernel::MakeGridParameters makeGridParameters{};

    makeGridParameters.num_columns = numX;
    makeGridParameters.num_rows = numY;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = lonMin;
    makeGridParameters.origin_y = latMin;
    makeGridParameters.block_size_x = 0.1;
    makeGridParameters.block_size_y = 0.1;

    return makeGridParameters;
}
