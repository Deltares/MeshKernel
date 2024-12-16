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
    const int projectionType = 0;

    const int numXCoords = 4;
    const int numYCoords = 4;
    const int numberOfEdges = (numXCoords - 1) * numYCoords + numXCoords * (numYCoords - 1);

    const int numXBathyCoords = 5;
    const int numYBathyCoords = 5;
    const int numberOfBathyCoordinates = numXBathyCoords * numYBathyCoords;

    const double origin = 0.0;
    const double delta = 1.0;

    int bathymetryPropertyId = -1;

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
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    std::vector<double> xs(meshData.num_nodes);
    std::vector<double> ys(meshData.num_nodes);
    std::vector<int> edges(meshData.num_edges * 2);

    meshData.node_x = xs.data();
    meshData.node_y = ys.data();
    meshData.edge_nodes = edges.data();

    mkapi::mkernel_mesh2d_get_node_edge_data(meshKernelId, meshData);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

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
            bathymetryData[count] = bathyX;

            bathyX += delta;
            ++count;
        }

        bathyY += delta;
    }

    mkapi::GeometryList sampleData{};

    sampleData.num_coordinates = numberOfBathyCoordinates;
    sampleData.values = bathymetryData.data();
    sampleData.coordinates_x = bathymetryXNodes.data();
    sampleData.coordinates_y = bathymetryYNodes.data();

    errorCode = mkapi::mkernel_mesh2d_set_property(projectionType, 0 /*use interpolation based on triangulation*/, sampleData, bathymetryPropertyId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    const double tolerance = 1.0e-13;

    std::vector<double> expectedInterpolatedData{0.0, 1.0, 2.0, 3.0,
                                                 0.0, 1.0, 2.0, 3.0,
                                                 0.0, 1.0, 2.0, 3.0,
                                                 0.5, 1.5, 2.5, 0.5,
                                                 1.5, 2.5, 0.5, 1.5,
                                                 2.5, 0.5, 1.5, 2.5};

    int sampleDataSize = -1;
    errorCode = mkapi::mkernel_mesh2d_get_property_dimension(meshKernelId, bathymetryPropertyId, sampleDataSize);

    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    ASSERT_EQ(sampleDataSize, numberOfEdges);

    mkapi::GeometryList propertyData{};
    std::vector<double> retrievedPropertyData(numberOfEdges, -1.0);
    propertyData.num_coordinates = numberOfEdges;
    propertyData.values = retrievedPropertyData.data();

    errorCode = mkapi::mkernel_mesh2d_get_property(meshKernelId, bathymetryPropertyId, propertyData);

    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    for (size_t i = 0; i < retrievedPropertyData.size(); ++i)
    {
        EXPECT_NEAR(retrievedPropertyData[i], expectedInterpolatedData[i], tolerance);
    }

    //--------------------------------
    // Remove kernel and property id from state.

    errorCode = mkapi::mkernel_expunge_state(meshKernelId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    errorCode = mkapi::mkernel_deallocate_property(bathymetryPropertyId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
}

TEST(MeshPropertyTests, PropertyFailureTest)
{
    const int projectionType = 0;

    const int numXCoords = 4;
    const int numYCoords = 4;
    const int numberOfCoordinates = numXCoords * numYCoords;

    const int numXBathyCoords = 5;
    const int numYBathyCoords = 5;
    const int numberOfBathyCoordinates = numXBathyCoords * numYBathyCoords;

    const double origin = 0.0;
    const double delta = 1.0;

    int bathymetryPropertyId = -1;

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
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    std::vector<double> xs(meshData.num_nodes);
    std::vector<double> ys(meshData.num_nodes);
    std::vector<int> edges(meshData.num_edges * 2);

    meshData.node_x = xs.data();
    meshData.node_y = ys.data();
    meshData.edge_nodes = edges.data();

    mkapi::mkernel_mesh2d_get_node_edge_data(meshKernelId, meshData);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

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
            bathymetryData[count] = bathyX;

            bathyX += delta;
            ++count;
        }

        bathyY += delta;
    }

    mkapi::GeometryList sampleData{};

    sampleData.num_coordinates = numberOfBathyCoordinates;
    sampleData.values = bathymetryData.data();
    sampleData.coordinates_x = bathymetryXNodes.data();
    sampleData.coordinates_y = bathymetryYNodes.data();

    bool hasBathymetryData = false;

    errorCode = mkapi::mkernel_mesh2d_is_valid_property(meshKernelId, bathymetryPropertyId, hasBathymetryData);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_FALSE(hasBathymetryData);

    errorCode = mkapi::mkernel_mesh2d_set_property(projectionType, 0 /*use interpolation based on triangulation*/, sampleData, bathymetryPropertyId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    errorCode = mkapi::mkernel_mesh2d_is_valid_property(meshKernelId, bathymetryPropertyId, hasBathymetryData);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_TRUE(hasBathymetryData);

    //--------------------------------
    // Start test

    int sampleDataSize = -1;

    // Expected test failure due to invalid meshkernel id
    errorCode = mkapi::mkernel_mesh2d_get_property_dimension(meshKernelId + 100, bathymetryPropertyId, sampleDataSize);
    ASSERT_EQ(mk::ExitCode::MeshKernelErrorCode, errorCode);

    mkapi::GeometryList propertyData{};
    std::vector<double> retrievedPropertyData(numberOfCoordinates, -1.0);
    propertyData.num_coordinates = numberOfCoordinates;

    // Expected test failure due to values set to null
    propertyData.values = nullptr;
    errorCode = mkapi::mkernel_mesh2d_get_property(meshKernelId, bathymetryPropertyId, propertyData);
    ASSERT_EQ(mk::ExitCode::ConstraintErrorCode, errorCode);

    propertyData.values = retrievedPropertyData.data();

    // Expected test failure due to invalid property id
    errorCode = mkapi::mkernel_mesh2d_get_property(meshKernelId, bathymetryPropertyId + 100, propertyData);
    ASSERT_EQ(mk::ExitCode::MeshKernelErrorCode, errorCode);

    // Expected test failure due to incorrect size
    propertyData.num_coordinates = numberOfCoordinates - 1;
    errorCode = mkapi::mkernel_mesh2d_get_property(meshKernelId, bathymetryPropertyId, propertyData);
    ASSERT_EQ(mk::ExitCode::ConstraintErrorCode, errorCode);

    //--------------------------------
    // Remove kernel and property id from state.

    errorCode = mkapi::mkernel_deallocate_property(bathymetryPropertyId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);

    errorCode = mkapi::mkernel_mesh2d_is_valid_property(meshKernelId, bathymetryPropertyId, hasBathymetryData);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
    EXPECT_FALSE(hasBathymetryData);

    errorCode = mkapi::mkernel_expunge_state(meshKernelId);
    ASSERT_EQ(mk::ExitCode::Success, errorCode);
}
