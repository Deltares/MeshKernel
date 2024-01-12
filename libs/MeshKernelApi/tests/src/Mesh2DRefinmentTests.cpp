#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <MeshKernel/Parameters.hpp>

#include <MeshKernelApi/GeometryList.hpp>
#include <MeshKernelApi/Mesh2D.hpp>
#include <MeshKernelApi/MeshKernel.hpp>

#include <TestUtils/Definitions.hpp>
#include <TestUtils/MakeCurvilinearGrids.hpp>
#include <TestUtils/MakeMeshes.hpp>
#include <TestUtils/SampleFileReader.hpp>

#include "CartesianApiTestFixture.hpp"

TEST_F(CartesianApiTestFixture, RefineAPolygonThroughApi)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometry_separator = meshkernel::constants::missing::doubleValue;

    std::vector xCoordinatesIn{76.251099, 498.503723, 505.253784, 76.251099};
    std::vector yCoordinatesIn{92.626556, 91.126541, 490.130554, 92.626556};
    std::vector valuesIn{0.0, 0.0, 0.0, 0.0};

    geometryListIn.coordinates_x = xCoordinatesIn.data();
    geometryListIn.coordinates_y = yCoordinatesIn.data();
    geometryListIn.values = valuesIn.data();
    geometryListIn.num_coordinates = static_cast<int>(xCoordinatesIn.size());

    // Execute
    int numberOfpolygonNodes;
    auto errorCode = mkernel_polygon_count_refine(meshKernelId, geometryListIn, 0, 2, 40, numberOfpolygonNodes);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(23, numberOfpolygonNodes);

    meshkernelapi::GeometryList geometryListOut;
    geometryListOut.num_coordinates = numberOfpolygonNodes;
    geometryListOut.geometry_separator = meshkernel::constants::missing::doubleValue;
    std::vector<double> xCoordinatesOut(numberOfpolygonNodes);
    std::vector<double> yCoordinatesOut(numberOfpolygonNodes);
    std::vector<double> valuesOut(numberOfpolygonNodes);
    geometryListOut.coordinates_x = xCoordinatesOut.data();
    geometryListOut.coordinates_y = yCoordinatesOut.data();
    geometryListOut.values = valuesOut.data();
    errorCode = mkernel_polygon_refine(meshKernelId, geometryListIn, false, 0, 2, geometryListOut);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    const double tolerance = 1e-6;
    ASSERT_NEAR(76.251099, geometryListOut.coordinates_x[0], tolerance);
    ASSERT_NEAR(92.626556, geometryListOut.coordinates_y[0], tolerance);
}

TEST_F(CartesianApiTestFixture, RefineBasedOnSamples_OnAUniformMesh_shouldRefineMesh)
{
    // Prepare
    MakeMesh(10, 10, 25);
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometry_separator = meshkernel::constants::missing::doubleValue;
    std::vector xCoordinatesIn{
        50.0,
        150.0,
        250.0,
        50.0,
        150.0,
        250.0,
        50.0,
        150.0,
        250.0};

    std::vector yCoordinatesIn{
        50.0,
        50.0,
        50.0,
        150.0,
        150.0,
        150.0,
        250.0,
        250.0,
        250.0};

    std::vector valuesIn{
        2.0,
        2.0,
        2.0,
        3.0,
        3.0,
        3.0,
        4.0,
        4.0,
        4.0};

    geometryListIn.coordinates_x = xCoordinatesIn.data();
    geometryListIn.coordinates_y = yCoordinatesIn.data();
    geometryListIn.values = valuesIn.data();
    geometryListIn.num_coordinates = static_cast<int>(valuesIn.size());

    meshkernel::MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 2;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.min_edge_size = 0.5;
    meshRefinementParameters.refinement_type = 1;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.account_for_samples_outside = 0;
    meshRefinementParameters.max_courant_time = 0.1;

    // Get the old state
    meshkernelapi::Mesh2D mesh2d{};
    auto errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(121, mesh2d.num_nodes);
    ASSERT_EQ(220, mesh2d.num_edges);

    // Execute
    errorCode = mkernel_mesh2d_refine_based_on_samples(meshKernelId, geometryListIn, 1.0, 1, meshRefinementParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Get the new state
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    ASSERT_EQ(1626, mesh2d.num_nodes);
    ASSERT_EQ(3225, mesh2d.num_edges);
}
TEST_F(CartesianApiTestFixture, RefineAGridBasedOnPolygonThroughApi)
{
    // Prepare
    MakeMesh(10, 10, 25);
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList polygon;
    polygon.geometry_separator = meshkernel::constants::missing::doubleValue;
    std::vector xCoordinatesIn{
        50.0,
        150.0,
        250.0,
        50.0,
        50.0};

    std::vector yCoordinatesIn{
        50.0,
        50.0,
        150.0,
        150.0,
        50.0};

    polygon.coordinates_x = xCoordinatesIn.data();
    polygon.coordinates_y = yCoordinatesIn.data();
    polygon.num_coordinates = static_cast<int>(xCoordinatesIn.size());

    meshkernel::MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 2;
    meshRefinementParameters.refine_intersected = 0;

    // Get the old state
    meshkernelapi::Mesh2D mesh2d{};
    auto errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(121, mesh2d.num_nodes);
    ASSERT_EQ(220, mesh2d.num_edges);

    // Execute
    errorCode = mkernel_mesh2d_refine_based_on_polygon(meshKernelId, polygon, meshRefinementParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Get the new state
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    ASSERT_EQ(279, mesh2d.num_nodes);
    ASSERT_EQ(595, mesh2d.num_edges);
}

TEST(MeshRefinement, Mesh2DRefineBasedOnGriddedSamples_WithGriddedSamples_ShouldRefineMesh)
{
    // Prepare
    int meshKernelId;
    constexpr int isGeographic = 0;
    meshkernelapi::mkernel_allocate_state(isGeographic, meshKernelId);

    auto [num_nodes, num_edges, node_x, node_y, node_type, edge_nodes, edge_type] =
        ReadLegacyMeshFile(TEST_FOLDER + "/data/MeshRefinementTests/gebco.nc");
    meshkernelapi::Mesh2D mesh2d;
    mesh2d.num_edges = static_cast<int>(num_edges);
    mesh2d.num_nodes = static_cast<int>(num_nodes);
    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_nodes = edge_nodes.data();

    auto errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    auto [numX, numY, xllCenter, yllCenter, cellSize, nodatavalue, values] = ReadAscFile<double>(TEST_FOLDER + "/data/MeshRefinementTests/gebco.asc");
    int interpolationType;
    errorCode = meshkernelapi::mkernel_get_interpolation_type_double(interpolationType);
    meshkernelapi::GriddedSamples griddedSamples;

    griddedSamples.num_x = numX;
    griddedSamples.num_y = numY;
    griddedSamples.x_origin = xllCenter;
    griddedSamples.y_origin = yllCenter;
    griddedSamples.cell_size = cellSize;
    griddedSamples.value_type = interpolationType;
    griddedSamples.values = values.data();
    griddedSamples.x_coordinates = nullptr;
    griddedSamples.y_coordinates = nullptr;

    meshkernel::MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 5;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.min_edge_size = 0.01;
    meshRefinementParameters.refinement_type = 1;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.account_for_samples_outside = 0;
    meshRefinementParameters.smoothing_iterations = 0;

    errorCode = mkernel_mesh2d_refine_based_on_gridded_samples(meshKernelId, griddedSamples, meshRefinementParameters, true);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2dResults;
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2dResults);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    ASSERT_EQ(417, mesh2dResults.num_nodes);
    ASSERT_EQ(988, mesh2dResults.num_edges);
    ASSERT_EQ(572, mesh2dResults.num_faces);
    ASSERT_EQ(1936, mesh2dResults.num_face_nodes);
}

TEST_F(CartesianApiTestFixture, Mesh2DRefineBasedOnGriddedSamples_WithNotUniformlySpacedSamples_ShouldRefineMesh)
{
    // Prepare
    meshkernel::UInt nRows{5};
    meshkernel::UInt nCols{4};
    MakeMesh(nRows, nCols, 100.0);
    auto const meshKernelId = GetMeshKernelId();

    int interpolationType;
    auto errorCode = meshkernelapi::mkernel_get_interpolation_type_double(interpolationType);

    meshkernelapi::GriddedSamples griddedSamples;
    griddedSamples.num_y = 6;
    griddedSamples.num_x = 7;
    std::vector<double> x_coordinates(griddedSamples.num_x);
    std::vector<double> y_coordinates(griddedSamples.num_y);

    double coordinate = -50.0;
    const double dx = 100.0;
    for (size_t i = 0; i < x_coordinates.size(); ++i)
    {
        x_coordinates[i] = coordinate + i * dx;
    }
    coordinate = -50.0;
    const double dy = 100.0;
    for (size_t i = 0; i < y_coordinates.size(); ++i)
    {
        y_coordinates[i] = coordinate + i * dy;
    }

    std::vector<double> values(griddedSamples.num_y * griddedSamples.num_x);
    for (size_t i = 0; i < values.size(); ++i)
    {
        values[i] = -0.05;
    }

    griddedSamples.x_coordinates = x_coordinates.data();
    griddedSamples.y_coordinates = y_coordinates.data();
    griddedSamples.values = values.data();
    griddedSamples.value_type = interpolationType;

    meshkernel::MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;
    meshRefinementParameters.min_edge_size = 2.0;
    meshRefinementParameters.refinement_type = 1;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.account_for_samples_outside = 0;
    meshRefinementParameters.max_num_refinement_iterations = 5;
    meshRefinementParameters.smoothing_iterations = 0;
    meshRefinementParameters.max_courant_time = 120.0;
    meshRefinementParameters.directional_refinement = 0;

    errorCode = mkernel_mesh2d_refine_based_on_gridded_samples(meshKernelId, griddedSamples, meshRefinementParameters, true);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2dResults;
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2dResults);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    ASSERT_EQ(86, mesh2dResults.num_nodes);
    ASSERT_EQ(161, mesh2dResults.num_edges);
    ASSERT_EQ(76, mesh2dResults.num_faces);
}

TEST(MeshRefinement, RefineBasedOnGriddedSamples_WithUniformSamplesAndSphericalCoordinates_ShouldRefineMesh2d)
{
    // Prepare
    int meshKernelId;
    constexpr int isSpherical = 1;
    meshkernelapi::mkernel_allocate_state(isSpherical, meshKernelId);

    const auto makeGridParameters = GebcoMakeGridParameters();

    auto errorCode = meshkernelapi::mkernel_mesh2d_make_rectangular_mesh(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    int interpolationType;
    errorCode = meshkernelapi::mkernel_get_interpolation_type_double(interpolationType);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    auto [ncols, nrows, xllcenter, yllcenter, cellsize, nodata_value, values] = ReadAscFile<double>(TEST_FOLDER + "/data/MeshRefinementTests/gebco.asc");
    meshkernelapi::GriddedSamples griddedSamples;
    griddedSamples.num_x = ncols;
    griddedSamples.num_y = nrows;
    griddedSamples.x_origin = xllcenter;
    griddedSamples.y_origin = yllcenter;
    griddedSamples.cell_size = cellsize;
    griddedSamples.value_type = interpolationType;
    griddedSamples.values = values.data();
    griddedSamples.x_coordinates = nullptr;
    griddedSamples.y_coordinates = nullptr;

    meshkernel::MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 1;
    meshRefinementParameters.min_edge_size = 0.01;
    meshRefinementParameters.refinement_type = 1;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.account_for_samples_outside = 0;
    meshRefinementParameters.max_num_refinement_iterations = 5;
    meshRefinementParameters.smoothing_iterations = 5;
    meshRefinementParameters.max_courant_time = 120;
    meshRefinementParameters.directional_refinement = 0;

    errorCode = mkernel_mesh2d_refine_based_on_gridded_samples(meshKernelId, griddedSamples, meshRefinementParameters, true);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2dResults;
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2dResults);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    ASSERT_EQ(5223, mesh2dResults.num_nodes);
    ASSERT_EQ(10745, mesh2dResults.num_edges);
    ASSERT_EQ(5523, mesh2dResults.num_faces);
    ASSERT_EQ(21212, mesh2dResults.num_face_nodes);
}

TEST(MeshRefinement, RefineBasedOnGriddedSamples_WithUniformSamplesAndSphericalCoordinatesAndLargeMinEdgeSize_ShouldNotRefineMesh2d)
{
    // Prepare
    int meshKernelId;
    constexpr int isSpherical = 1;
    meshkernelapi::mkernel_allocate_state(isSpherical, meshKernelId);

    const auto makeGridParameters = GebcoMakeGridParameters();

    auto errorCode = meshkernelapi::mkernel_mesh2d_make_rectangular_mesh(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    int interpolationType;
    errorCode = meshkernelapi::mkernel_get_interpolation_type_double(interpolationType);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    auto [ncols, nrows, xllcenter, yllcenter, cellsize, nodata_value, values] = ReadAscFile<double>(TEST_FOLDER + "/data/MeshRefinementTests/gebco.asc");
    meshkernelapi::GriddedSamples griddedSamples;
    griddedSamples.num_x = ncols;
    griddedSamples.num_y = nrows;
    griddedSamples.x_origin = xllcenter;
    griddedSamples.y_origin = yllcenter;
    griddedSamples.value_type = interpolationType;
    griddedSamples.cell_size = cellsize;
    griddedSamples.values = values.data();
    griddedSamples.x_coordinates = nullptr;
    griddedSamples.y_coordinates = nullptr;

    meshkernel::MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 1;
    meshRefinementParameters.refinement_type = 1;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.account_for_samples_outside = 0;
    meshRefinementParameters.max_num_refinement_iterations = 5;
    meshRefinementParameters.smoothing_iterations = 5;
    meshRefinementParameters.max_courant_time = 120;
    meshRefinementParameters.directional_refinement = 0;

    // Set a large value of min edge size
    meshRefinementParameters.min_edge_size = 1e9;

    errorCode = mkernel_mesh2d_refine_based_on_gridded_samples(meshKernelId, griddedSamples, meshRefinementParameters, true);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2dResults;
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2dResults);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert no refinement took place
    ASSERT_EQ(54, (makeGridParameters.num_columns + 1) * (makeGridParameters.num_rows + 1));
}

TEST(MeshRefinement, RefineAGridBasedOnPolygonThroughApi_OnSpericalCoordinateWithLargeMinEdgeSize_ShouldNotRefine)
{
    // Prepare
    int isGeographic = 1;
    int meshKernelId = -1;
    auto errorCode = meshkernelapi::mkernel_allocate_state(isGeographic, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::MakeGridParameters makeGridParameters;

    makeGridParameters.origin_x = -6.0;
    makeGridParameters.origin_y = 48.5;
    makeGridParameters.upper_right_x = 2;
    makeGridParameters.upper_right_y = 51.2;
    makeGridParameters.block_size_x = 0.5;
    makeGridParameters.block_size_y = 0.5;

    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid_on_extension(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_convert_to_mesh2d(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Get the old state
    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(221, mesh2d.num_nodes);
    ASSERT_EQ(412, mesh2d.num_edges);

    std::vector xCoordinatesIn{-5.0, -4.0, 0.0, -5.0};
    std::vector yCoordinatesIn{49.0, 51.0, 49.5, 49.0};
    std::vector valuesIn{1.0, 1.0, 1.0, 1.0};

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometry_separator = meshkernel::constants::missing::doubleValue;
    geometryListIn.coordinates_x = xCoordinatesIn.data();
    geometryListIn.coordinates_y = yCoordinatesIn.data();
    geometryListIn.values = valuesIn.data();
    geometryListIn.num_coordinates = static_cast<int>(xCoordinatesIn.size());

    meshkernel::MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 10;
    meshRefinementParameters.min_edge_size = 200000.0;

    // Execute
    errorCode = mkernel_mesh2d_refine_based_on_polygon(meshKernelId, geometryListIn, meshRefinementParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Get the new state
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    ASSERT_EQ(221, mesh2d.num_nodes);
    ASSERT_EQ(412, mesh2d.num_edges);
}

TEST(MeshRefinement, RefineAGridBasedOnPolygonThroughApi_OnSpericalCoordinateWithSmallMinEdgeSize_ShouldRefine)
{
    // Prepare
    int isGeographic = 1;
    int meshKernelId = -1;
    auto errorCode = meshkernelapi::mkernel_allocate_state(isGeographic, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::MakeGridParameters makeGridParameters;

    makeGridParameters.origin_x = -6.0;
    makeGridParameters.origin_y = 48.5;
    makeGridParameters.upper_right_x = 2;
    makeGridParameters.upper_right_y = 51.2;
    makeGridParameters.block_size_x = 0.5;
    makeGridParameters.block_size_y = 0.5;

    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid_on_extension(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_convert_to_mesh2d(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Get the new state
    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(221, mesh2d.num_nodes);
    ASSERT_EQ(412, mesh2d.num_edges);

    std::vector xCoordinatesIn{-5.0, -4.0, 0.0, -5.0};
    std::vector yCoordinatesIn{49.0, 51.0, 49.5, 49.0};
    std::vector valuesIn{1.0, 1.0, 1.0, 1.0};
    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometry_separator = meshkernel::constants::missing::doubleValue;
    geometryListIn.coordinates_x = xCoordinatesIn.data();
    geometryListIn.coordinates_y = yCoordinatesIn.data();
    geometryListIn.values = valuesIn.data();
    geometryListIn.num_coordinates = static_cast<int>(xCoordinatesIn.size());

    meshkernel::MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 10;

    // The minimum edge size
    meshRefinementParameters.min_edge_size = 2000.0;

    // Execute
    errorCode = mkernel_mesh2d_refine_based_on_polygon(meshKernelId, geometryListIn, meshRefinementParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(1570, mesh2d.num_nodes);
    ASSERT_EQ(3361, mesh2d.num_edges);
}

class MeshRefinementValueTypes : public ::testing::TestWithParam<meshkernel::InterpolationValuesTypes>
{
public:
    [[nodiscard]] static std::vector<meshkernel::InterpolationValuesTypes> GetData()
    {
        return {meshkernel::InterpolationValuesTypes::shortType,
                meshkernel::InterpolationValuesTypes::intType,
                meshkernel::InterpolationValuesTypes::floatType,
                meshkernel::InterpolationValuesTypes::doubleType};
    }
};

TEST_P(MeshRefinementValueTypes, parameters)
{
    // Get the test parameters
    auto const interpolationValueType = GetParam();

    // Prepare
    int meshKernelId;
    constexpr int isSpherical = 1;
    meshkernelapi::mkernel_allocate_state(isSpherical, meshKernelId);

    const auto makeGridParameters = GebcoMakeGridParameters();

    auto errorCode = meshkernelapi::mkernel_mesh2d_make_rectangular_mesh(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::GriddedSamples griddedSamples;

    std::vector<short> valuesShort;
    std::vector<int> valuesInt;
    std::vector<float> valuesFloat;
    std::vector<double> valuesDouble;

    int interpolationType;
    if (interpolationValueType == meshkernel::InterpolationValuesTypes::shortType)
    {

        errorCode = meshkernelapi::mkernel_get_interpolation_type_short(interpolationType);
        ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
        auto [ncols, nrows, xllcenter, yllcenter, cellsize, nodata_value, values] = ReadAscFile<short>(TEST_FOLDER + "/data/MeshRefinementTests/gebcoIntegers.asc");
        valuesShort = values;
        griddedSamples.num_x = ncols;
        griddedSamples.num_y = nrows;
        griddedSamples.x_origin = xllcenter;
        griddedSamples.y_origin = yllcenter;
        griddedSamples.cell_size = cellsize;
        griddedSamples.value_type = interpolationType;
        griddedSamples.values = valuesShort.data();
        griddedSamples.x_coordinates = nullptr;
        griddedSamples.y_coordinates = nullptr;
    }
    else if (interpolationValueType == meshkernel::InterpolationValuesTypes::intType)
    {
        errorCode = meshkernelapi::mkernel_get_interpolation_type_int(interpolationType);
        ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
        auto [ncols, nrows, xllcenter, yllcenter, cellsize, nodata_value, values] = ReadAscFile<int>(TEST_FOLDER + "/data/MeshRefinementTests/gebcoIntegers.asc");
        valuesInt = values;
        griddedSamples.num_x = ncols;
        griddedSamples.num_y = nrows;
        griddedSamples.x_origin = xllcenter;
        griddedSamples.y_origin = yllcenter;
        griddedSamples.cell_size = cellsize;
        griddedSamples.value_type = interpolationType;
        griddedSamples.values = valuesInt.data();
        griddedSamples.x_coordinates = nullptr;
        griddedSamples.y_coordinates = nullptr;
    }
    else if (interpolationValueType == meshkernel::InterpolationValuesTypes::floatType)
    {
        errorCode = meshkernelapi::mkernel_get_interpolation_type_float(interpolationType);
        ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
        auto [ncols, nrows, xllcenter, yllcenter, cellsize, nodata_value, values] = ReadAscFile<float>(TEST_FOLDER + "/data/MeshRefinementTests/gebcoIntegers.asc");
        valuesFloat = values;
        griddedSamples.num_x = ncols;
        griddedSamples.num_y = nrows;
        griddedSamples.x_origin = xllcenter;
        griddedSamples.y_origin = yllcenter;
        griddedSamples.cell_size = cellsize;
        griddedSamples.value_type = interpolationType;
        griddedSamples.values = valuesFloat.data();
        griddedSamples.x_coordinates = nullptr;
        griddedSamples.y_coordinates = nullptr;
    }
    else if (interpolationValueType == meshkernel::InterpolationValuesTypes::doubleType)
    {
        errorCode = meshkernelapi::mkernel_get_interpolation_type_double(interpolationType);
        ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
        auto [ncols, nrows, xllcenter, yllcenter, cellsize, nodata_value, values] = ReadAscFile<double>(TEST_FOLDER + "/data/MeshRefinementTests/gebcoIntegers.asc");
        valuesDouble = values;
        griddedSamples.num_x = ncols;
        griddedSamples.num_y = nrows;
        griddedSamples.x_origin = xllcenter;
        griddedSamples.y_origin = yllcenter;
        griddedSamples.cell_size = cellsize;
        griddedSamples.value_type = interpolationType;
        griddedSamples.values = valuesDouble.data();
        griddedSamples.x_coordinates = nullptr;
        griddedSamples.y_coordinates = nullptr;
    }

    meshkernel::MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 1;
    meshRefinementParameters.min_edge_size = 0.01;
    meshRefinementParameters.refinement_type = 1;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.account_for_samples_outside = 0;
    meshRefinementParameters.max_num_refinement_iterations = 5;
    meshRefinementParameters.smoothing_iterations = 5;
    meshRefinementParameters.max_courant_time = 120;
    meshRefinementParameters.directional_refinement = 0;

    errorCode = mkernel_mesh2d_refine_based_on_gridded_samples(meshKernelId, griddedSamples, meshRefinementParameters, true);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2dResults;
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2dResults);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::mkernel_deallocate_state(meshKernelId);

    ASSERT_EQ(5223, mesh2dResults.num_nodes);
    ASSERT_EQ(10745, mesh2dResults.num_edges);
    ASSERT_EQ(5523, mesh2dResults.num_faces);
    ASSERT_EQ(21212, mesh2dResults.num_face_nodes);
}

INSTANTIATE_TEST_SUITE_P(MeshRefinement, MeshRefinementValueTypes, ::testing::ValuesIn(MeshRefinementValueTypes::GetData()));
