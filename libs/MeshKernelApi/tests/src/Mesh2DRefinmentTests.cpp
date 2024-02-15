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
#include "MeshKernel/AveragingInterpolation.hpp"
#include "MeshKernel/MeshRefinement.hpp"
#include "MeshKernel/SamplesHessianCalculator.hpp"
#include "SampleFileWriter.hpp"
#include "SampleGenerator.hpp"

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

TEST_F(CartesianApiTestFixture, RefineBasedOnSamplesWaveCourant_OnAUniformMesh_shouldRefineMesh)
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

TEST(MeshRefinement, Mesh2DRefineBasedOnGriddedSamplesWaveCourant_WithGriddedSamples_ShouldRefineMesh)
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

    auto [numX, numY, xllCenter, yllCenter, cellSize, nodatavalue, values] = ReadAscFile<float>(TEST_FOLDER + "/data/MeshRefinementTests/gebco.asc");
    int interpolationType;
    errorCode = meshkernelapi::mkernel_get_interpolation_type_float(interpolationType);
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

TEST_F(CartesianApiTestFixture, Mesh2DRefineBasedOnGriddedSamplesWaveCourant_WithNotUniformlySpacedSamples_ShouldRefineMesh)
{
    // Prepare
    meshkernel::UInt nRows{5};
    meshkernel::UInt nCols{4};
    MakeMesh(nRows, nCols, 100.0);
    auto const meshKernelId = GetMeshKernelId();

    int interpolationType;
    auto errorCode = meshkernelapi::mkernel_get_interpolation_type_float(interpolationType);

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

    std::vector<float> values(griddedSamples.num_y * griddedSamples.num_x);
    for (size_t i = 0; i < values.size(); ++i)
    {
        values[i] = -0.05f;
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

TEST(MeshRefinement, RefineBasedOnGriddedSamplesWaveCourant_WithUniformSamplesAndSphericalCoordinates_ShouldRefineMesh2d)
{
    // Prepare
    int meshKernelId;
    constexpr int isSpherical = 1;
    meshkernelapi::mkernel_allocate_state(isSpherical, meshKernelId);

    const auto makeGridParameters = GebcoMakeGridParameters();

    auto errorCode = meshkernelapi::mkernel_mesh2d_make_rectangular_mesh(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    int interpolationType;
    errorCode = meshkernelapi::mkernel_get_interpolation_type_float(interpolationType);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    auto [ncols, nrows, xllcenter, yllcenter, cellsize, nodata_value, values] = ReadAscFile<float>(TEST_FOLDER + "/data/MeshRefinementTests/gebco.asc");
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

TEST(MeshRefinement, RefineBasedOnGriddedSamplesWaveCourant_WithUniformSamplesAndSphericalCoordinatesAndLargeMinEdgeSize_ShouldNotRefineMesh2d)
{
    // Prepare
    int meshKernelId;
    constexpr int isSpherical = 1;
    meshkernelapi::mkernel_allocate_state(isSpherical, meshKernelId);

    const auto makeGridParameters = GebcoMakeGridParameters();

    auto errorCode = meshkernelapi::mkernel_mesh2d_make_rectangular_mesh(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    int interpolationType;
    errorCode = meshkernelapi::mkernel_get_interpolation_type_float(interpolationType);
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

class MeshRefinementSampleValueTypes : public ::testing::TestWithParam<meshkernel::InterpolationValues>
{
public:
    [[nodiscard]] static std::vector<meshkernel::InterpolationValues> GetData()
    {
        return {meshkernel::InterpolationValues::shortType,
                meshkernel::InterpolationValues::floatType};
    }
};

TEST_P(MeshRefinementSampleValueTypes, parameters)
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
    std::vector<float> valuesFloat;

    int interpolationType;
    if (interpolationValueType == meshkernel::InterpolationValues::shortType)
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
    else if (interpolationValueType == meshkernel::InterpolationValues::floatType)
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

INSTANTIATE_TEST_SUITE_P(MeshRefinement, MeshRefinementSampleValueTypes, ::testing::ValuesIn(MeshRefinementSampleValueTypes::GetData()));

TEST_F(CartesianApiTestFixture, RefineAMeshBasedOnRidgeRefinement_OnAUniformMesh_shouldRefineMesh)
{
    // Prepare
    meshkernel::UInt nRows{21};
    meshkernel::UInt nCols{41};
    const double meshDelta = 10.0;
    MakeMesh(nRows, nCols, 10.0);

    // Generate gridded samples, from a gaussian distribution
    const int numSamplesXCoordinates = (nCols - 1) * 2 + 1;
    const int numSamplesYCoordinates = (nRows * 1) * 2 + 1;
    const double deltaX = 5.0;
    const double deltaY = 5.0;
    const auto sampleData = generateSampleData(FunctionTestCase::GaussianBump, numSamplesXCoordinates, numSamplesYCoordinates, deltaX, deltaY);

    // Create an instance of gridded samples,
    meshkernelapi::GriddedSamples griddedSamples;
    griddedSamples.num_x = numSamplesXCoordinates;
    griddedSamples.num_y = numSamplesYCoordinates;
    griddedSamples.x_origin = 0.0;
    griddedSamples.y_origin = 0.0;
    griddedSamples.cell_size = meshDelta * 0.1;

    // Flatten the values before passing them to Mesh Kernel
    meshkernel::UInt index = 0;
    std::vector<float> values(numSamplesXCoordinates * numSamplesYCoordinates, 0.0);
    for (int j = 0; j < numSamplesYCoordinates; ++j)
    {
        for (int i = 0; i < numSamplesXCoordinates; ++i)
        {
            const auto griddedIndex = griddedSamples.num_y * i + j;
            values[index] = static_cast<float>(sampleData[griddedIndex].value);
            index++;
        }
    }
    int interpolationType;
    auto errorCode = meshkernelapi::mkernel_get_interpolation_type_float(interpolationType);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    griddedSamples.values = values.data();
    griddedSamples.value_type = interpolationType;

    // Create meshrefinement parameters
    meshkernel::MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 3;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;
    meshRefinementParameters.min_edge_size = 2.0;
    meshRefinementParameters.account_for_samples_outside = 0;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.smoothing_iterations = 0;
    meshRefinementParameters.refinement_type = 3;

    // Execute
    const auto meshKernelId = GetMeshKernelId();
    mkernel_mesh2d_refine_ridges_based_on_gridded_samples(meshKernelId, griddedSamples, 1.01, 1, 0, meshRefinementParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert on the mesh values
    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(956, mesh2d.num_nodes);
    ASSERT_EQ(1872, mesh2d.num_edges);
}

TEST(MeshRefinement, RefineAGridBasedOnSamplesThroughApi_OnSpericalCoordinateWithRealSamples_ShouldRefine)
{
    // Prepare
    int isGeographic = 1;
    int meshKernelId = -1;
    auto errorCode = meshkernelapi::mkernel_allocate_state(isGeographic, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::MakeGridParameters makeGridParameters;
    makeGridParameters.angle = 0;
    makeGridParameters.origin_x = -68.55;
    makeGridParameters.origin_y = 11.8;
    makeGridParameters.upper_right_x = -67.9;
    makeGridParameters.upper_right_y = 12.6;
    makeGridParameters.block_size_x = 0.05;
    makeGridParameters.block_size_y = 0.05;

    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid_on_extension(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_convert_to_mesh2d(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> lon{
        -68.54791667,
        -68.46458333,
        -68.38125,
        -68.29791667,
        -68.21458333,
        -68.13125,
        -68.04791667,
        -67.96458333,
    };

    std::vector<double> lat{
        11.80208333,
        11.88541667,
        11.96875,
        12.05208333,
        12.13541667,
        12.21875,
        12.30208333,
        12.38541667,
        12.46875,
        12.55208333,
    };

    std::vector<float> values{
        -1700.0,
        -1769.0,
        -1688.0,
        -1641.0,
        -1526.0,
        -1291.0,
        -1121.0,
        -1537.0,
        -1561.0,
        -1674.0,
        -1354.0,
        -757.0,
        -837.0,
        -838.0,
        -1080.0,
        -1466.0,
        -1630.0,
        -1390.0,
        -710.0,
        -562.0,
        -479.0,
        -753.0,
        -1246.0,
        -1703.0,
        -1553.0,
        -1446.0,
        -1147.0,
        -248.0,
        -175.0,
        -712.0,
        -1621.0,
        -1920.0,
        -1503.0,
        -1380.0,
        -1080.0,
        -305.0,
        18.0,
        -543.0,
        -1563.0,
        -2241.0,
        -1477.0,
        -1571.0,
        -3.0,
        100.0,
        11.0,
        -891.0,
        -1521.0,
        -2446.0,
        -1892.0,
        -1808.0,
        16.0,
        -3102.0,
        -2015.0,
        -1302.0,
        -1484.0,
        -2581.0,
        -2516.0,
        -2091.0,
        -1957.0,
        -2647.0,
        -1422.0,
        -1486.0,
        -2340.0,
        -2702.0,
        -2689.0,
        -2353.0,
        -2614.0,
        -3612.0,
        -3058.0,
        -3017.0,
        -3181.0,
        -2848.0,
        -3110.0,
        -3025.0,
        -3861.0,
        -3927.0,
        -3818.0,
        -4162.0,
        -4386.0,
        -4504.0,
    };

    meshkernelapi::GriddedSamples griddedSamples;

    griddedSamples.x_coordinates = lon.data();
    griddedSamples.y_coordinates = lat.data();
    griddedSamples.values = values.data();

    meshkernel::MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.min_edge_size = 300;
    meshRefinementParameters.refinement_type = static_cast<int>(meshkernel::MeshRefinement::RefinementType::WaveCourant);
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.smoothing_iterations = 2;
    meshRefinementParameters.max_courant_time = 120;

    errorCode = meshkernelapi::mkernel_mesh2d_refine_based_on_gridded_samples(meshKernelId, griddedSamples, meshRefinementParameters, true);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    // Get the new state

    // Assert
    meshkernelapi::Mesh2D mesh2d{};
    errorCode = meshkernelapi::mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(252, mesh2d.num_nodes);
    ASSERT_EQ(472, mesh2d.num_edges);
}
