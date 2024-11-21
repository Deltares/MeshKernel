#include <algorithm>
#include <gtest/gtest.h>
#include <random>
#include <vector>

#include "CartesianApiTestFixture.hpp"
#include "MeshKernel/Parameters.hpp"
#include "MeshKernelApi/BoundingBox.hpp"
#include "MeshKernelApi/Mesh2D.hpp"
#include "MeshKernelApi/MeshKernel.hpp"
#include "TestUtils/MakeMeshes.hpp"

// namespace aliases
namespace mk = meshkernel;
namespace mkapi = meshkernelapi;

TEST(MeshPropertyTests, BathymetryTest)
{

    const int numXCoords = 4;
    const int numYCoords = 4;
    const int numberOfCoordinates = numXCoords * numYCoords;

    const int numXBathyCoords = 6;
    const int numYBathyCoords = 6;
    const int numberOfBathyCoordinates = numXBathyCoords * numYBathyCoords;

    const double origin = 0.0;
    const double delta = 1.0;

    int meshKernelId = meshkernel::constants::missing::intValue;
    int errorCode;

    // Clear the meshkernel state and undo stack before starting the test.
    errorCode = mkapi::mkernel_clear_state();
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    errorCode = mkapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    meshkernel::MakeGridParameters makeGridParameters;

    // num_columns and num_rows indicate number of elements in each direction
    makeGridParameters.num_columns = numXCoords - 1;
    makeGridParameters.num_rows = numYCoords - 1;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = origin;
    makeGridParameters.origin_y = origin;
    makeGridParameters.block_size_x = delta;
    makeGridParameters.block_size_y = delta;

    // Generate curvilinear grid.
    errorCode = mkapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    // Convert the curvilinear grid to an unstructured grid.
    errorCode = mkapi::mkernel_curvilinear_convert_to_mesh2d(meshKernelId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    //--------------------------------

    mkapi::Mesh2D meshData{};
    mkapi::mkernel_mesh2d_get_dimensions(meshKernelId, meshData);

    std::vector<double> xs(meshData.num_nodes);
    std::vector<double> ys(meshData.num_nodes);
    std::vector<int> edges(meshData.num_edges * 2);

    std::cout << "meshData.num_nodes " << meshData.num_nodes << std::endl;

    meshData.node_x = xs.data();
    meshData.node_y = ys.data();
    meshData.edge_nodes = edges.data();

    mkapi::mkernel_mesh2d_get_node_edge_data(meshKernelId, meshData);

    for (int i = 0; i < meshData.num_nodes; ++i)
    {
        std::cout << "node " << i << " = " << xs[i] << ", " << ys[i] << std::endl;
    }

    // return;

    //--------------------------------

    std::vector<double> bathymetryData(numberOfBathyCoordinates, 1.0);
    std::vector<double> bathymetryXNodes(numberOfBathyCoordinates);
    std::vector<double> bathymetryYNodes(numberOfBathyCoordinates);

    double bathyY = origin - 0.5 * delta;
    size_t count = 0;

    for (size_t i = 0; i < numXBathyCoords; ++i)
    {
        double bathyX = origin - 0.5 * delta;

        for (size_t j = 0; j < numYBathyCoords; ++j)
        {
            bathymetryXNodes[count] = bathyX;
            bathymetryYNodes[count] = bathyY;

            bathyX += delta;
            ++count;
        }

        bathyY += delta;
    }

    mkapi::GeometryList sampleData{};

    sampleData.num_coordinates = numberOfCoordinates;
    sampleData.values = bathymetryData.data();
    sampleData.coordinates_x = bathymetryXNodes.data();
    sampleData.coordinates_y = bathymetryYNodes.data();

    // mkapi::GeometryList sampleData{.num_coordinates = numberOfCoordinates,
    //                                .coordinates_x = bathymetryXNodes.data(),
    //                                .coordinates_y = bathymetryYNodes.data(),
    //                                .values = bathymetryData.data()};

    errorCode = mkapi::mkernel_mesh2d_set_bathymetry_data(meshKernelId, sampleData);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    int sampleDataSize = -1;
    errorCode = mkapi::mkernel_mesh2d_get_bathymetry_dimension(meshKernelId, sampleDataSize);
    ASSERT_EQ(sampleDataSize, numberOfCoordinates);

    std::vector<double> retrievedData(numberOfCoordinates, -1.0);
    sampleData.values = retrievedData.data();

    sampleData.coordinates_x = nullptr;
    sampleData.coordinates_y = nullptr;

    errorCode = mkapi::mkernel_mesh2d_get_bathymetry_data(meshKernelId, sampleData);

    // const double tolerance = 1.0e-13;

    for (size_t i = 0; i < retrievedData.size(); ++i)
    {
        std::cout << " retrievedData[i] " << i << " = " << retrievedData[i] << std::endl;

        // EXPECT_NEAR(retrievedData[i], 1.0, tolerance);
    }
}
