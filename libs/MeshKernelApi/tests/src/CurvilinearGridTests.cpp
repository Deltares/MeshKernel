#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <MeshKernel/Parameters.hpp>
#include <MeshKernelApi/Mesh2D.hpp>
#include <MeshKernelApi/MeshKernel.hpp>

#include "CartesianApiTestFixture.hpp"
#include "MeshKernelApi/BoundingBox.hpp"

#include <TestUtils/MakeCurvilinearGrids.hpp>
#include <TestUtils/MakeMeshes.hpp>

using namespace meshkernel;

TEST(CurvilinearGrid, MakeRectangular_OnSphericalCoordinates_ShouldMakeCurvilinearGrid)
{
    MakeGridParameters makeGridParameters;

    makeGridParameters.origin_x = -1.0;
    makeGridParameters.origin_y = 49.1;
    makeGridParameters.block_size_x = 0.1;
    makeGridParameters.block_size_y = 0.1;
    makeGridParameters.num_columns = 8;
    makeGridParameters.num_rows = 11;

    int meshKernelId = 0;
    const int projectionType = 1;
    auto errorCode = meshkernelapi::mkernel_allocate_state(projectionType, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGridResults;
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGridResults);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    ASSERT_EQ(9, curvilinearGridResults.num_m);
    ASSERT_EQ(12, curvilinearGridResults.num_n);
}

TEST(CurvilinearGrid, MakeRectangular_OnSphericalCoordinatesWithpolygon_ShouldMakeCurvilinearGrid)
{
    // Setup
    const double lonMin = -6.0;
    const double lonMax = 2.0;

    const double latMin = 48.5;
    const double latMax = 51.2;

    const double lonRes = 0.2;
    const double latRes = 0.2;

    MakeGridParameters makeGridParameters;
    makeGridParameters.origin_x = -6;
    makeGridParameters.origin_y = 15.5;
    makeGridParameters.num_rows = static_cast<int>(std::ceil((latMax - latMin) / latRes));
    makeGridParameters.num_columns = static_cast<int>(std::ceil((lonMax - lonMin) / lonRes));
    makeGridParameters.block_size_x = lonRes;
    makeGridParameters.block_size_y = latRes;

    int meshKernelId = 0;
    int projectionType = 1;
    auto errorCode = meshkernelapi::mkernel_allocate_state(projectionType, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    meshkernelapi::GeometryList geometryList{};
    auto coordinates_x = std::vector{-6.0, -4.0, 0.0, -6.0};
    auto coordinates_y = std::vector{48.0, 51.0, 49.5, 48.0};
    geometryList.coordinates_x = coordinates_x.data();
    geometryList.coordinates_y = coordinates_y.data();
    geometryList.num_coordinates = static_cast<int>(coordinates_x.size());

    errorCode = mkernel_curvilinear_compute_rectangular_grid_from_polygon(meshKernelId, makeGridParameters, geometryList);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGridResults;
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGridResults);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    ASSERT_EQ(61, curvilinearGridResults.num_m);
    ASSERT_EQ(125, curvilinearGridResults.num_n);
}

TEST_F(CartesianApiTestFixture, CurvilinearComputeTransfiniteFromPolygon_ShouldComputeAValidCurvilinearGrid)
{
    /*

    Input polygon:
    6---5---4
    |       |
    7       3
    |       |
    0---1---2

    Generated curvilinearGrid:

    6---7---8
    |   |   |
    3---4---5
    |   |   |
    0---1---2

    */

    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometry_separator = constants::missing::doubleValue;
    std::vector<double> xCoordinatesIn{0, 5, 10, 10, 10, 5, 0, 0, 0};
    std::vector<double> yCoordinatesIn{0, 0, 0, 5, 10, 10, 10, 5, 0};
    std::vector<double> valuesIn{0, 0, 0, 0, 0, 0, 0, 0, 0};

    geometryListIn.coordinates_x = xCoordinatesIn.data();
    geometryListIn.coordinates_y = yCoordinatesIn.data();
    geometryListIn.values = valuesIn.data();
    geometryListIn.num_coordinates = static_cast<int>(xCoordinatesIn.size());

    // Execute
    auto errorCode = mkernel_curvilinear_compute_transfinite_from_polygon(meshKernelId,
                                                                          geometryListIn,
                                                                          0,
                                                                          2,
                                                                          4,
                                                                          false);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Get the new state
    meshkernelapi::CurvilinearGrid curvilinear_grid{};

    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinear_grid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    ASSERT_EQ(3, curvilinear_grid.num_m);
    ASSERT_EQ(3, curvilinear_grid.num_n);
}

TEST_F(CartesianApiTestFixture, MakeCurvilinearGridThroughApi)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();

    MakeGridParameters makeGridParameters;

    makeGridParameters.num_columns = 3;
    makeGridParameters.num_rows = 2;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = 1.0;
    makeGridParameters.block_size_y = 1.0;

    // Execute
    auto errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);

    // Assert
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(4, curvilinearGrid.num_m);
    ASSERT_EQ(3, curvilinearGrid.num_n);

    // Allocate memory and get data
    std::vector<double> node_x(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> node_y(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = node_x.data();
    curvilinearGrid.node_y = node_y.data();
    errorCode = mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    /*  8---9--10---11
        |   |   |   |
        4---5---6---7
        |   |   |   |
        0---1---2---3
    */
    // Assert data
    const double tolerance = 1e-6;
    // Nodes
    ASSERT_NEAR(0.0, curvilinearGrid.node_x[0], tolerance);
    ASSERT_NEAR(0.0, curvilinearGrid.node_y[0], tolerance);
    ASSERT_NEAR(1.0, curvilinearGrid.node_x[1], tolerance);
    ASSERT_NEAR(0.0, curvilinearGrid.node_y[1], tolerance);
}

TEST_F(CartesianApiTestFixture, GenerateTransfiniteCurvilinearGridThroughApi)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometry_separator = constants::missing::doubleValue;
    std::vector xCoordinates{1.340015E+02, 3.642529E+02, 6.927549E+02, constants::missing::doubleValue,
                             2.585022E+02, 4.550035E+02, 8.337558E+02, constants::missing::doubleValue,
                             1.002513E+02, 4.610035E+02, constants::missing::doubleValue,
                             6.522547E+02, 7.197551E+02};

    std::vector yCoordinates{2.546282E+02, 4.586302E+02, 5.441311E+02, constants::missing::doubleValue,
                             6.862631E+01, 2.726284E+02, 3.753794E+02, constants::missing::doubleValue,
                             4.068797E+02, 7.912642E+01, constants::missing::doubleValue,
                             6.026317E+02, 2.681283E+02};

    std::vector zCoordinates{0.0, 0.0, 0.0, constants::missing::doubleValue,
                             0.0, 0.0, 0.0, constants::missing::doubleValue,
                             0.0, 0.0, constants::missing::doubleValue,
                             0.0, 0.0};

    geometryListIn.coordinates_x = xCoordinates.data();
    geometryListIn.coordinates_y = yCoordinates.data();
    geometryListIn.values = zCoordinates.data();
    geometryListIn.num_coordinates = static_cast<int>(xCoordinates.size());

    CurvilinearParameters curvilinearParameters;

    curvilinearParameters.m_refinement = 10;
    curvilinearParameters.n_refinement = 10;
    curvilinearParameters.smoothing_iterations = 10;
    curvilinearParameters.smoothing_parameter = 0.5;
    curvilinearParameters.attraction_parameter = 0.0;

    // Execute
    auto errorCode = mkernel_curvilinear_compute_transfinite_from_splines(meshKernelId,
                                                                          geometryListIn,
                                                                          curvilinearParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);

    // Assert
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(11, curvilinearGrid.num_m);
    ASSERT_EQ(11, curvilinearGrid.num_n);
}

TEST_F(CartesianApiTestFixture, GenerateOrthogonalCurvilinearGridThroughApi)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;

    geometryListIn.geometry_separator = constants::missing::doubleValue;
    std::vector xCoordinates{1.175014E+02, 3.755030E+02, 7.730054E+02, constants::missing::doubleValue,
                             4.100089E+01, 3.410027E+02};

    std::vector yCoordinates{2.437587E+01, 3.266289E+02, 4.563802E+02, constants::missing::doubleValue,
                             2.388780E+02, 2.137584E+01};

    std::vector zCoordinates{0.0, 0.0, 0.0, constants::missing::doubleValue,
                             0.0, 0.0};

    geometryListIn.coordinates_x = xCoordinates.data();
    geometryListIn.coordinates_y = yCoordinates.data();
    geometryListIn.values = zCoordinates.data();
    geometryListIn.num_coordinates = static_cast<int>(xCoordinates.size());

    CurvilinearParameters curvilinearParameters;
    curvilinearParameters.m_refinement = 40;
    curvilinearParameters.n_refinement = 10;
    SplinesToCurvilinearParameters splinesToCurvilinearParameters;
    splinesToCurvilinearParameters.aspect_ratio = 0.1;
    splinesToCurvilinearParameters.aspect_ratio_grow_factor = 1.1;
    splinesToCurvilinearParameters.average_width = 500.0;
    splinesToCurvilinearParameters.curvature_adapted_grid_spacing = 1;
    splinesToCurvilinearParameters.maximum_num_faces_in_uniform_part = 5;
    splinesToCurvilinearParameters.nodes_on_top_of_each_other_tolerance = 0.0001;
    splinesToCurvilinearParameters.min_cosine_crossing_angles = 0.95;
    splinesToCurvilinearParameters.check_front_collisions = 0;
    splinesToCurvilinearParameters.remove_skinny_triangles = 1;
    splinesToCurvilinearParameters.grow_grid_outside = 0;

    // Execute
    auto errorCode = mkernel_curvilinear_initialize_orthogonal_grid_from_splines(meshKernelId,
                                                                                 geometryListIn,
                                                                                 curvilinearParameters,
                                                                                 splinesToCurvilinearParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Grow grid, from the second layer
    for (auto layer = 1; layer <= curvilinearParameters.n_refinement; ++layer)
    {
        errorCode = meshkernelapi::mkernel_curvilinear_iterate_orthogonal_grid_from_splines(meshKernelId, layer);
        ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    }

    // Puts the computed curvilinear mesh into the mesh state (unstructured mesh)
    errorCode = meshkernelapi::mkernel_curvilinear_refresh_orthogonal_grid_from_splines(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Delete the mesh curvilinearGridFromSplinesInstances vector entry
    errorCode = meshkernelapi::mkernel_curvilinear_delete_orthogonal_grid_from_splines(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Get the new state

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);

    // Assert
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(3, curvilinearGrid.num_n);
    ASSERT_EQ(7, curvilinearGrid.num_m);
}

TEST_F(CartesianApiTestFixture, RefineCompute_OnCurvilinearGrid_ShouldRefine)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();

    MakeRectangularCurvilinearGrid(3, 3);

    auto errorCode = meshkernelapi::mkernel_curvilinear_refine(meshKernelId, 10.0, 20.0, 20.0, 20.0, 10);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    ASSERT_EQ(13, curvilinearGrid.num_m);
    ASSERT_EQ(4, curvilinearGrid.num_n);
}

TEST_F(CartesianApiTestFixture, DerefineCompute_OnCurvilinearGrid_ShouldDeRefine)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();

    MakeRectangularCurvilinearGrid();

    // Execute
    auto errorCode = meshkernelapi::mkernel_curvilinear_derefine(meshKernelId, 10.0, 20.0, 30.0, 20.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    ASSERT_EQ(4, curvilinearGrid.num_m);
    ASSERT_EQ(5, curvilinearGrid.num_n);
}

TEST_F(CartesianApiTestFixture, Orthogonalize_CurvilinearGrid_ShouldOrthogonalize)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();

    MakeRectangularCurvilinearGrid(3, 3);
    // Move a node to make the grid non orthogonal
    auto errorCode = meshkernelapi::mkernel_curvilinear_move_node(meshKernelId, 10.0, 20.0, 18.0, 12.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    OrthogonalizationParameters orthogonalizationParameters{};
    orthogonalizationParameters.outer_iterations = 1;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;

    // Execute
    errorCode = meshkernelapi::mkernel_curvilinear_initialize_orthogonalize(meshKernelId, orthogonalizationParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    errorCode = meshkernelapi::mkernel_curvilinear_set_block_orthogonalize(meshKernelId, 0.0, 0.0, 30.0, 30.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    errorCode = meshkernelapi::mkernel_curvilinear_orthogonalize(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert (nothing changed)
    ASSERT_EQ(4, curvilinearGrid.num_m);
    ASSERT_EQ(4, curvilinearGrid.num_n);

    std::vector<double> xNodesCurvilinearGrid(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> yNodesCurvilinearGrid(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = xNodesCurvilinearGrid.data();
    curvilinearGrid.node_y = yNodesCurvilinearGrid.data();

    errorCode = mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    double const tolerance = 1e-6;
    ASSERT_NEAR(11.841396536135521, curvilinearGrid.node_x[9], tolerance);
    ASSERT_NEAR(18.158586078094562, curvilinearGrid.node_y[9], tolerance);
}

TEST_F(CartesianApiTestFixture, Smoothing_CurvilinearGrid_ShouldSmooth)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();

    MakeRectangularCurvilinearGrid();

    // Execute
    auto errorCode = meshkernelapi::mkernel_curvilinear_smoothing(meshKernelId, 10, 10.0, 20.0, 30.0, 20.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert (nothing changed)
    ASSERT_EQ(5, curvilinearGrid.num_m);
    ASSERT_EQ(5, curvilinearGrid.num_n);
}

TEST_F(CartesianApiTestFixture, ComputedDirectionalSmooth_CurvilinearGrid_ShouldSmooth)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();

    MakeRectangularCurvilinearGrid();

    // Execute
    auto errorCode = meshkernelapi::mkernel_curvilinear_smoothing_directional(meshKernelId,
                                                                              10,
                                                                              10.0,
                                                                              0.0,
                                                                              10.0,
                                                                              30.0,
                                                                              10.0,
                                                                              0.0,
                                                                              30.0,
                                                                              0.0);

    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert (nothing changed)
    ASSERT_EQ(5, curvilinearGrid.num_m);
    ASSERT_EQ(5, curvilinearGrid.num_n);
}

TEST_F(CartesianApiTestFixture, ComputedLineShift_CurvilinearGrid_ShouldShift)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();

    MakeRectangularCurvilinearGrid();

    meshkernelapi::mkernel_curvilinear_initialize_line_shift(meshKernelId);

    auto errorCode = meshkernelapi::mkernel_curvilinear_set_line_line_shift(meshKernelId, 0.0, 0.0, 0.0, 30.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_set_block_line_shift(meshKernelId, 0.0, 0.0, 30.0, 30.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    /// Move a gridline point, in this case the origin to -10.0, 0.0
    errorCode = meshkernelapi::mkernel_curvilinear_move_node_line_shift(meshKernelId, 0.0, 0.0, -10.0, 0.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    errorCode = meshkernelapi::mkernel_curvilinear_line_shift(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    errorCode = meshkernelapi::mkernel_curvilinear_finalize_line_shift(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert, the nodes along the grid line have changed
    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> xNodesCurvilinearGrid(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> yNodesCurvilinearGrid(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = xNodesCurvilinearGrid.data();
    curvilinearGrid.node_y = yNodesCurvilinearGrid.data();

    errorCode = mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(-10.0, curvilinearGrid.node_x[0]);
    ASSERT_EQ(2.5, curvilinearGrid.node_x[1]);
    ASSERT_EQ(17.5, curvilinearGrid.node_x[2]);
    ASSERT_EQ(30.0, curvilinearGrid.node_x[3]);
}

TEST_F(CartesianApiTestFixture, CurvilinearComputeOrthogonalGridFromSplines_ShouldMakeCurvilinearGrid)
{
    // Setup
    meshkernelapi::GeometryList splines{};
    double geometrySeparator = meshkernelapi::mkernel_get_separator();

    std::vector<double> coordinatesX{
        7.7979524E+04,
        7.7979524E+04,
        7.8302860E+04,
        7.9732343E+04,
        8.0889543E+04,
        8.0668314E+04,
        7.9579184E+04,
        geometrySeparator,
        7.6618112E+04,
        7.6754253E+04,
        7.7179694E+04,
        7.8404966E+04,
        7.9681290E+04,
        8.0140766E+04,
        7.9477078E+04,
        7.8779354E+04,
        geometrySeparator,
        7.7281800E+04,
        7.7366889E+04,
        7.7928471E+04,
        7.9153742E+04,
        8.0242872E+04,
        8.0481119E+04,
        7.9970590E+04,
        7.9579184E+04,
        7.9170760E+04,
        geometrySeparator,
        7.613792E+04,
        7.831719E+04,
        geometrySeparator,
        7.857202E+04,
        8.003072E+04};

    std::vector<double> coordinatesY{
        3.7127829E+05,
        3.7025723E+05,
        3.6898090E+05,
        3.6809598E+05,
        3.6698984E+05,
        3.6578158E+05,
        3.6419894E+05,
        geometrySeparator,
        3.7136337E+05,
        3.7005301E+05,
        3.6874265E+05,
        3.6780668E+05,
        3.6721107E+05,
        3.6636018E+05,
        3.6544123E+05,
        3.6452228E+05,
        geometrySeparator,
        3.7144846E+05,
        3.6984880E+05,
        3.6874265E+05,
        3.6792581E+05,
        3.6722808E+05,
        3.6641124E+05,
        3.6542421E+05,
        3.6484561E+05,
        3.6431806E+05,
        geometrySeparator,
        3.712157E+05,
        3.710751E+05,
        geometrySeparator,
        3.649151E+05,
        3.641506E+05};

    splines.coordinates_x = coordinatesX.data();
    splines.coordinates_y = coordinatesY.data();
    splines.num_coordinates = static_cast<int>(coordinatesX.size());

    SplinesToCurvilinearParameters splinesToCurvilinearParameters;
    CurvilinearParameters curvilinearParameters{};

    curvilinearParameters.m_refinement = 20;
    curvilinearParameters.n_refinement = 40;
    splinesToCurvilinearParameters.aspect_ratio = 0.5;
    splinesToCurvilinearParameters.aspect_ratio_grow_factor = 1.0;
    splinesToCurvilinearParameters.average_width = 500.0;
    splinesToCurvilinearParameters.curvature_adapted_grid_spacing = true;
    splinesToCurvilinearParameters.grow_grid_outside = 0;
    splinesToCurvilinearParameters.maximum_num_faces_in_uniform_part = 8;

    splinesToCurvilinearParameters.nodes_on_top_of_each_other_tolerance = 1e-4;
    splinesToCurvilinearParameters.min_cosine_crossing_angles = 0.95;
    splinesToCurvilinearParameters.check_front_collisions = false;
    splinesToCurvilinearParameters.remove_skinny_triangles = true;

    // Execute, with large length threshold
    auto const meshKernelId = GetMeshKernelId();
    auto errorCode = mkernel_curvilinear_compute_orthogonal_grid_from_splines(meshKernelId,
                                                                              splines,
                                                                              curvilinearParameters,
                                                                              splinesToCurvilinearParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert one curvilinear grid is produced
    ASSERT_GT(curvilinearGrid.num_m, 0);
    ASSERT_GT(curvilinearGrid.num_n, 0);
}

TEST_F(CartesianApiTestFixture, CurvilinearGridFromSplines_ShouldMakeCurvilinearGrid)
{
    // Setup
    meshkernelapi::GeometryList splines{};
    double geometrySeparator = meshkernelapi::mkernel_get_separator();

    std::vector<double> coordinatesX{
        -1.0, 11.0,
        geometrySeparator,
        10.0, 10.0,
        geometrySeparator,
        11.0, -1.0,
        geometrySeparator,
        0.0, 0.0};

    std::vector<double> coordinatesY{
        0.0, 0.0,
        geometrySeparator,
        11.0, -1.0,
        geometrySeparator,
        10.0, 10.0,
        geometrySeparator,
        11.0, -1.0};

    splines.coordinates_x = coordinatesX.data();
    splines.coordinates_y = coordinatesY.data();
    splines.num_coordinates = static_cast<int>(coordinatesX.size());

    CurvilinearParameters curvilinearParameters{};

    curvilinearParameters.m_refinement = 5;
    curvilinearParameters.n_refinement = 10;

    // Execute, with large length threshold
    auto const meshKernelId = GetMeshKernelId();
    auto errorCode = mkernel_curvilinear_compute_grid_from_splines(meshKernelId,
                                                                   splines,
                                                                   curvilinearParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Expected size of the generated grid
    EXPECT_EQ(curvilinearGrid.num_n, 11);
    EXPECT_EQ(curvilinearGrid.num_m, 6);
}

TEST_F(CartesianApiTestFixture, CurvilinearGridFromSplines_ShouldThrowAlgorithmError)
{
    // Setup
    meshkernelapi::GeometryList splines{};
    double geometrySeparator = meshkernelapi::mkernel_get_separator();

    std::vector<double> coordinatesX{
        -1.0, 11.0,
        geometrySeparator,
        10.0, 10.0,
        geometrySeparator,
        11.0, -1.0,
        geometrySeparator,
        0.0, 0.0,
        geometrySeparator,
        15.0, 15.0};

    std::vector<double> coordinatesY{
        0.0, 0.0,
        geometrySeparator,
        11.0, -1.0,
        geometrySeparator,
        10.0, 10.0,
        geometrySeparator,
        11.0, -1.0,
        geometrySeparator,
        15.0, -1.0};

    splines.coordinates_x = coordinatesX.data();
    splines.coordinates_y = coordinatesY.data();
    splines.num_coordinates = static_cast<int>(coordinatesX.size());

    CurvilinearParameters curvilinearParameters{};

    curvilinearParameters.m_refinement = 5;
    curvilinearParameters.n_refinement = 5;

    // Execute, with large length threshold
    auto const meshKernelId = GetMeshKernelId();
    auto errorCode = mkernel_curvilinear_compute_grid_from_splines(meshKernelId,
                                                                   splines,
                                                                   curvilinearParameters);
    ASSERT_EQ(meshkernel::ExitCode::AlgorithmErrorCode, errorCode);
}

TEST_F(CartesianApiTestFixture, CurvilinearInsertFace_ShouldInsertAFace)
{
    // Setup
    MakeRectangularCurvilinearGrid();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    auto errorCode = meshkernelapi::mkernel_curvilinear_insert_face(meshKernelId, -5.0, 5.0);

    // Assert
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> node_x(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> node_y(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = node_x.data();
    curvilinearGrid.node_y = node_y.data();
    errorCode = mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert two extra nodes have been inserted (before it was 5 by 5 = 25 nodes, not it is 25 + 2 = 27)
    auto const numValidNodes = CurvilinearGridCountValidNodes(curvilinearGrid);
    ASSERT_EQ(numValidNodes, 27);
}

TEST_F(CartesianApiTestFixture, CurvilinearLineMirror_ShouldInsertANewGridLine)
{
    // Setup
    MakeRectangularCurvilinearGrid();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    auto errorCode = meshkernelapi::mkernel_curvilinear_line_mirror(meshKernelId,
                                                                    1.2,
                                                                    0.0,
                                                                    0.0,
                                                                    0.0,
                                                                    50.0);

    // Assert
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert that 5 nodes have been inserted on the bottom boundary, so the total count is 30
    ASSERT_EQ(curvilinearGrid.num_m * curvilinearGrid.num_n, 30);
}

TEST_F(CartesianApiTestFixture, CurvilinearLineAttractionRepulsion_ShouldAttractGridlines)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();

    MakeRectangularCurvilinearGrid(5, 5);

    // Execute
    auto errorCode = meshkernelapi::mkernel_curvilinear_line_attraction_repulsion(meshKernelId,
                                                                                  0.5,
                                                                                  30.0, 0.0, 30.0, 50.0,
                                                                                  10.0, 20.0, 50.0, 20.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> node_x(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> node_y(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = node_x.data();
    curvilinearGrid.node_y = node_y.data();
    errorCode = mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert data
    const double tolerance = 1e-6;
    // Nodes
    ASSERT_NEAR(17.5, curvilinearGrid.node_x[2], tolerance);
    ASSERT_NEAR(0.0, curvilinearGrid.node_y[2], tolerance);
    ASSERT_NEAR(42.5, curvilinearGrid.node_x[4], tolerance);
    ASSERT_NEAR(0.0, curvilinearGrid.node_y[4], tolerance);
}

TEST_F(CartesianApiTestFixture, CurvilinearDeleteNode_ShouldDeleteNode)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();
    MakeRectangularCurvilinearGrid(5, 5);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    auto errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    std::vector<double> node_x(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> node_y(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = node_x.data();
    curvilinearGrid.node_y = node_y.data();
    errorCode = mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    auto const numValidNodesBefore = CurvilinearGridCountValidNodes(curvilinearGrid);

    // Execute
    errorCode = meshkernelapi::mkernel_curvilinear_delete_node(meshKernelId, 10.0, 0.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Asserts
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Two nodes are removed, one by delete node and one by the administration. The node is at the corner and the entire face will be removed
    auto const numValidNodesAfter = CurvilinearGridCountValidNodes(curvilinearGrid);
    ASSERT_EQ(numValidNodesBefore - 2, numValidNodesAfter);
}

TEST(CurvilinearGrid, MakeRectangularOnDefinedExtension_OnSphericalCoordinates_ShouldMakeCurvilinearGrid)
{
    // Setup
    MakeGridParameters makeGridParameters;
    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.upper_right_x = 10.0;
    makeGridParameters.upper_right_y = 10.0;
    makeGridParameters.block_size_x = 1.0;
    makeGridParameters.block_size_y = 2.0;

    int meshKernelId = 0;
    int projectionType = 0;
    auto errorCode = meshkernelapi::mkernel_allocate_state(projectionType, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid_on_extension(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGridResults;
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGridResults);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    ASSERT_EQ(11, curvilinearGridResults.num_m);
    ASSERT_EQ(6, curvilinearGridResults.num_n);
}

TEST_F(CartesianApiTestFixture, CurvilinearDeleteExterior_OnCurvilinearGrid_ShouldDeleteExterior)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();
    MakeRectangularCurvilinearGrid(10, 10, 1.0, 1.0, 0.0, 0.0);

    // Execute
    auto errorCode = meshkernelapi::mkernel_curvilinear_delete_exterior(meshKernelId, 5.0, 5.0, 15.0, 15.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> node_x(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> node_y(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = node_x.data();
    curvilinearGrid.node_y = node_y.data();
    errorCode = mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    double const tolerance = 1e-6;
    ASSERT_NEAR(meshkernel::constants::missing::doubleValue, curvilinearGrid.node_x[0], tolerance);
    ASSERT_NEAR(meshkernel::constants::missing::doubleValue, curvilinearGrid.node_y[0], tolerance);
    ASSERT_NEAR(5.0, curvilinearGrid.node_x[60], tolerance);
    ASSERT_NEAR(5.0, curvilinearGrid.node_y[60], tolerance);
}

TEST_F(CartesianApiTestFixture, CurvilinearDeleteInterior_OnCurvilinearGrid_ShouldDeleteExterior)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();
    MakeRectangularCurvilinearGrid(10, 10, 1.0, 1.0, 0.0, 0.0);

    // Execute
    auto errorCode = meshkernelapi::mkernel_curvilinear_delete_interior(meshKernelId, 5.0, 5.0, 15.0, 15.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(11, curvilinearGrid.num_m);
    ASSERT_EQ(11, curvilinearGrid.num_n);

    std::vector<double> node_x(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> node_y(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = node_x.data();
    curvilinearGrid.node_y = node_y.data();
    errorCode = mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    double const tolerance = 1e-6;
    ASSERT_NEAR(0.0, curvilinearGrid.node_x[0], tolerance);
    ASSERT_NEAR(0.0, curvilinearGrid.node_y[0], tolerance);
    ASSERT_NEAR(meshkernel::constants::missing::doubleValue, curvilinearGrid.node_x[72], tolerance);
    ASSERT_NEAR(meshkernel::constants::missing::doubleValue, curvilinearGrid.node_y[72], tolerance);
}

TEST(CurvilinearGrid, MakeRectangular_ComputeSmoothnessTest)
{
    MakeGridParameters makeGridParameters;

    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = 10.0;
    makeGridParameters.block_size_y = 10.0;
    makeGridParameters.num_columns = 4;
    makeGridParameters.num_rows = 5;

    int meshKernelId = 0;
    const int projectionType = 0;
    auto errorCode = meshkernelapi::mkernel_allocate_state(projectionType, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGridResults;
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGridResults);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> smoothness(curvilinearGridResults.num_n * curvilinearGridResults.num_m);
    std::vector<double> expectedX{-999.0, 1.0, 1.0, 1.0, -999.0, -999.0, 1.0, 1.0, 1.0, -999.0, -999.0, 1.0, 1.0, 1.0, -999.0, -999.0, 1.0, 1.0, 1.0, -999.0, -999.0, 1.0, 1.0, 1.0, -999.0, -999.0, 1.0, 1.0, 1.0, -999.0};
    std::vector<double> expectedY{-999.0, -999.0, -999.0, -999.0, -999.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -999.0, -999.0, -999.0, -999.0, -999.0};

    // Test x direction
    errorCode = meshkernelapi::mkernel_curvilinear_compute_smoothness(meshKernelId, 0, smoothness.data());
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    constexpr double tolerance = 1.0e-13;

    for (size_t i = 0; i < smoothness.size(); ++i)
    {
        EXPECT_NEAR(expectedX[i], smoothness[i], tolerance);
    }

    // Now test y direction
    errorCode = meshkernelapi::mkernel_curvilinear_compute_smoothness(meshKernelId, 1, smoothness.data());
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    for (size_t i = 0; i < smoothness.size(); ++i)
    {
        EXPECT_NEAR(expectedY[i], smoothness[i], tolerance);
    }
}

TEST(CurvilinearGrid, MakeRectangular_ComputeCurvatureTest)
{
    MakeGridParameters makeGridParameters;

    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = 10.0;
    makeGridParameters.block_size_y = 10.0;
    makeGridParameters.num_columns = 4;
    makeGridParameters.num_rows = 5;

    int meshKernelId = 0;
    const int projectionType = 0;
    auto errorCode = meshkernelapi::mkernel_allocate_state(projectionType, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGridResults;
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGridResults);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> curvature(curvilinearGridResults.num_n * curvilinearGridResults.num_m);

    std::vector<double> expectedX{-999.0, 0.001000001000001, 0.001000001000001, 0.001000001000001, -999.0,
                                  -999.0, 0.001000001000001, 0.001000001000001, 0.001000001000001, -999.0,
                                  -999.0, 0.001000001000001, 0.001000001000001, 0.001000001000001, -999.0,
                                  -999.0, 0.001000001000001, 0.001000001000001, 0.001000001000001, -999.0,
                                  -999.0, 0.001000001000001, 0.001000001000001, 0.001000001000001, -999.0,
                                  -999.0, 0.001000001000001, 0.001000001000001, 0.001000001000001, -999.0};

    std::vector<double> expectedY{-999.0, -999.0, -999.0, -999.0, -999.0,
                                  0.001000001000001, 0.001000001000001, 0.001000001000001, 0.001000001000001,
                                  0.001000001000001, 0.001000001000001, 0.001000001000001, 0.001000001000001,
                                  0.001000001000001, 0.001000001000001, 0.001000001000001, 0.001000001000001,
                                  0.001000001000001, 0.001000001000001, 0.001000001000001, 0.001000001000001,
                                  0.001000001000001, 0.001000001000001, 0.001000001000001, 0.001000001000001,
                                  -999.0, -999.0, -999.0, -999.0, -999.0};

    // Test x direction
    errorCode = meshkernelapi::mkernel_curvilinear_compute_curvature(meshKernelId, 0, curvature.data());
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    constexpr double tolerance = 1.0e-13;

    for (size_t i = 0; i < curvature.size(); ++i)
    {
        EXPECT_NEAR(expectedX[i], curvature[i], tolerance);
    }

    // Now test y direction
    errorCode = meshkernelapi::mkernel_curvilinear_compute_curvature(meshKernelId, 1, curvature.data());
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    for (size_t i = 0; i < curvature.size(); ++i)
    {
        EXPECT_NEAR(expectedY[i], curvature[i], tolerance);
    }
}

TEST(CurvilinearGrid, MakeRectangular_ConvertToMesh2D)
{
    // Prepare
    MakeGridParameters makeGridParameters;

    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = 10.0;
    makeGridParameters.block_size_y = 10.0;
    makeGridParameters.num_columns = 5;
    makeGridParameters.num_rows = 5;

    int meshKernelId = 0;
    const int projectionType = 0;
    auto errorCode = meshkernelapi::mkernel_allocate_state(projectionType, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGridResults;
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGridResults);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    errorCode = meshkernelapi::mkernel_curvilinear_convert_to_mesh2d(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Data for translated mesh
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

    // assert
    const double tolerance = 1e-8;
    ASSERT_NEAR(0.0, node_x[0], tolerance);
    ASSERT_NEAR(10.0, node_x[1], tolerance);
    ASSERT_NEAR(20.0, node_x[2], tolerance);
    ASSERT_NEAR(30.0, node_x[3], tolerance);
    ASSERT_NEAR(40.0, node_x[4], tolerance);

    ASSERT_NEAR(0.0, node_y[0], tolerance);
    ASSERT_NEAR(0.0, node_y[1], tolerance);
    ASSERT_NEAR(10.0, node_y[6], tolerance);
    ASSERT_NEAR(10.0, node_y[7], tolerance);
    ASSERT_NEAR(20.0, node_y[12], tolerance);
    ASSERT_NEAR(20.0, node_y[13], tolerance);

    ASSERT_EQ(mesh2d.num_nodes, 36);
    ASSERT_EQ(mesh2d.num_edges, 60);
}

class CurvilinearLocationIndexTests : public testing::TestWithParam<std::tuple<Point, Location, int>>
{
public:
    [[nodiscard]] static std::vector<std::tuple<Point, Location, int>> GetData()
    {
        return {
            std::make_tuple<Point, Location, int>(Point{-1.0, -1.0}, Location::Nodes, 0),
            std::make_tuple<Point, Location, int>(Point{18.0, 18.0}, Location::Nodes, 10),
            std::make_tuple<Point, Location, int>(Point{21.0, 21.0}, Location::Nodes, 10),
            std::make_tuple<Point, Location, int>(Point{20.0, 20.0}, Location::Nodes, 10),
            std::make_tuple<Point, Location, int>(Point{19.0, 19.0}, Location::Nodes, 10),
            std::make_tuple<Point, Location, int>(Point{5.0, -1.0}, Location::Edges, 12),
            std::make_tuple<Point, Location, int>(Point{11.0, 13.0}, Location::Edges, 5),
            std::make_tuple<Point, Location, int>(Point{0.5, 0.5}, Location::Faces, 0),
            std::make_tuple<Point, Location, int>(Point{18.0, 18.0}, Location::Faces, 4),
            std::make_tuple<Point, Location, int>(Point{7.0, 14.0}, Location::Faces, 3)};
    }
};

TEST_P(CurvilinearLocationIndexTests, GetLocationIndex_OnACurvilinearGrid_ShouldGetTheLocationIndex)
{
    // Prepare
    auto const& [point, location, expectedIndex] = GetParam();

    MakeGridParameters makeGridParameters;

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

    // Execute
    int locationIndex = -1;
    const meshkernelapi::BoundingBox boundingBox;
    auto const locationInt = static_cast<int>(location);
    errorCode = mkernel_curvilinear_get_location_index(meshKernelId, point.x, point.y, locationInt, boundingBox, locationIndex);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    ASSERT_EQ(locationIndex, expectedIndex);
}
INSTANTIATE_TEST_SUITE_P(LocationIndexParametrizedTests, CurvilinearLocationIndexTests, ::testing::ValuesIn(CurvilinearLocationIndexTests::GetData()));

class CurvilinearCornersTests : public testing::TestWithParam<std::tuple<Point, Point, int, int>>
{
public:
    [[nodiscard]] static std::vector<std::tuple<Point, Point, int, int>> GetData()
    {
        return {
            std::make_tuple<Point, Point, int, int>(Point{-0.5, 0.5}, Point{-0.3, 1.5}, 4, 3),
            std::make_tuple<Point, Point, int, int>(Point{-0.5, 0.5}, Point{-0.5, 1.1}, 4, 3),
            std::make_tuple<Point, Point, int, int>(Point{2.5, 0.5}, Point{2.5, 1.1}, 4, 3),
            std::make_tuple<Point, Point, int, int>(Point{2.5, 0.5}, Point{2.1, 1.5}, 4, 3),
            std::make_tuple<Point, Point, int, int>(Point{0.5, -0.5}, Point{1.1, -0.5}, 3, 4),
            std::make_tuple<Point, Point, int, int>(Point{0.5, -0.5}, Point{1.5, -0.1}, 3, 4),
            std::make_tuple<Point, Point, int, int>(Point{1.5, -0.5}, Point{0.9, -0.5}, 3, 4),
            std::make_tuple<Point, Point, int, int>(Point{1.5, -0.5}, Point{0.5, -0.1}, 3, 4)};
    }
};

TEST_P(CurvilinearCornersTests, InsertAFace_OnACurvilinearGrid_ShouldInsertFaces)
{
    // Prepare
    auto const& [firstPoint, secondPoint, expectedM, expectedN] = GetParam();

    MakeGridParameters makeGridParameters;

    makeGridParameters.origin_x = 0.0;

    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = 1.0;
    makeGridParameters.block_size_y = 1.0;
    makeGridParameters.num_columns = 2;
    makeGridParameters.num_rows = 2;

    int meshKernelId = 0;
    const int projectionType = 0;
    auto errorCode = meshkernelapi::mkernel_allocate_state(projectionType, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    errorCode = meshkernelapi::mkernel_curvilinear_insert_face(meshKernelId, firstPoint.x, firstPoint.y);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_insert_face(meshKernelId, secondPoint.x, secondPoint.y);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(expectedM, curvilinearGrid.num_m);
    ASSERT_EQ(expectedN, curvilinearGrid.num_n);
}
INSTANTIATE_TEST_SUITE_P(CurvilinearCornersParametrizedTests, CurvilinearCornersTests, ::testing::ValuesIn(CurvilinearCornersTests::GetData()));

TEST(CurvilinearGrid, SnapToLandBoundary)
{
    constexpr double tolerance = 1.0e-10;

    // Prepare
    int meshKernelId;
    auto errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    MakeGridParameters makeGridParameters;

    makeGridParameters.num_columns = 10;
    makeGridParameters.num_rows = 10;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = 10.0;
    makeGridParameters.block_size_y = 10.0;

    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> node_x(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> node_y(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = node_x.data();
    curvilinearGrid.node_y = node_y.data();
    errorCode = mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Make copy of node values.
    std::vector<double> originalNodeX(node_x);
    std::vector<double> originalNodeY(node_y);

    std::vector<double> expectedPointsX{0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0,
                                        0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0,
                                        0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0,
                                        0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0,
                                        0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0,
                                        0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0,
                                        0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0,
                                        0.01774961128485, 10.01533144852, 19.98050480592, 29.97597166392, 39.97143852191, 50.0447171313, 60.0262335099, 69.96945303784, 79.9410737541, 90.08118201565, 100.0,
                                        0.1556264747728, 10.13442431202, 19.8290684635, 29.78932241507, 39.74957636663, 50.39207447388, 60.23001228601, 69.73216787907, 79.48334170402, 90.71179423073, 100.0,
                                        0.3316761521277, 10.28648942044, 19.6357052077, 29.55099714999, 39.46628909227, 50.83560173828, 60.49020958726, 69.42918756326, 78.89888121011, 91.51699877479, 100.0,
                                        0.4105750388196, 10.35463931961, 19.54904702213, 29.44418867202, 39.3393303219, 51.03437408427, 60.60682029451, 69.29340310761, 78.63694785711, 91.87786136221, 100.0};

    std::vector<double> expectedPointsY{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                        10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,
                                        20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0,
                                        30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0,
                                        40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0,
                                        50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
                                        60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0,
                                        70.23666072723, 70.20441866004, 70.18938128575, 70.23341738295, 70.27745348015, 70.21158750721, 70.12412878028, 70.11524478944, 70.22231155967, 70.20664439596, 70.0,
                                        82.07501303012, 81.79231843049, 81.66047252617, 82.0465757743, 82.43267902243, 81.85517402742, 81.08834634086, 81.01045256871, 81.94920124036, 81.81183341756, 80.0,
                                        94.42233455747, 93.81984672774, 93.53885249278, 94.36172816257, 95.18460383235, 93.95380659904, 92.31951875183, 92.1535090375, 94.15420042167, 93.86143769632, 90.0,
                                        105.4743163503, 104.728509145, 104.3806722017, 105.3992929494, 106.4179136971, 104.8943352951, 102.8712842195, 102.6657842326, 105.1423986574, 104.7799937436, 100.0};

    // Set up land boundary

    // std::vector<double> landXCoordinates({11.0000, 11.5878, 11.9511, 11.9511, 11.5878, 11.0000, 10.4122, 10.0489, 10.0489, 10.4122, 11.0000});
    // std::vector<double> landYCoordinates({0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0});

    std::vector<double> landXCoordinates({-10.614098, -1.765560, 17.628498, 42.355930, 65.871231, 85.992554, 99.568405});
    std::vector<double> landYCoordinates({100.910210, 105.637527, 104.182968, 106.728447, 101.758698, 107.092087, 101.758698});

    // Not part of the unit test but added to ensure arrays are the same size
    ASSERT_EQ(landXCoordinates.size(), landYCoordinates.size());

    meshkernelapi::GeometryList land;
    land.geometry_separator = constants::missing::doubleValue;
    land.inner_outer_separator = constants::missing::innerOuterSeparator;
    land.num_coordinates = static_cast<int>(landXCoordinates.size());
    land.coordinates_x = landXCoordinates.data();
    land.coordinates_y = landYCoordinates.data();

    double sectionControlPoint1x = 0.0;
    double sectionControlPoint1y = 100.0;
    double sectionControlPoint2x = 90.0;
    double sectionControlPoint2y = 100.0;

    errorCode = mkernel_curvilinear_snap_to_landboundary(meshKernelId, land,
                                                         sectionControlPoint1x,
                                                         sectionControlPoint1y,
                                                         sectionControlPoint2x,
                                                         sectionControlPoint2y);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    for (int i = 0; i < curvilinearGrid.num_m * curvilinearGrid.num_n; ++i)
    {
        EXPECT_NEAR(curvilinearGrid.node_x[i], expectedPointsX[i], tolerance);
        EXPECT_NEAR(curvilinearGrid.node_y[i], expectedPointsY[i], tolerance);
    }

    // Now undo snapping.

    bool didUndoOfDeleteNode = false;
    int undoId = meshkernel::constants::missing::intValue;

    errorCode = meshkernelapi::mkernel_undo_state(didUndoOfDeleteNode, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndoOfDeleteNode);
    ASSERT_EQ(meshKernelId, undoId);

    errorCode = mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    for (int i = 0; i < curvilinearGrid.num_m * curvilinearGrid.num_n; ++i)
    {
        EXPECT_NEAR(curvilinearGrid.node_x[i], originalNodeX[i], tolerance);
        EXPECT_NEAR(curvilinearGrid.node_y[i], originalNodeY[i], tolerance);
    }
}

TEST(CurvilinearGrid, SnapToSpline)
{
    constexpr double tolerance = 2.0e-5;

    // Prepare
    int meshKernelId;
    auto errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    MakeGridParameters makeGridParameters;

    makeGridParameters.num_columns = 10;
    makeGridParameters.num_rows = 10;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = 1.0;
    makeGridParameters.block_size_y = 1.0;

    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> node_x(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> node_y(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = node_x.data();
    curvilinearGrid.node_y = node_y.data();
    errorCode = mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Make copy of node values.
    std::vector<double> originalNodeX(node_x);
    std::vector<double> originalNodeY(node_y);

    std::vector<double> expectedPointsX{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
                                        0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.014873267871414519, 8.372359722435612994, 9.916182851714186341, 11.17180940417338775,
                                        0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.021612916987076680, 8.541090219036453490, 10.33134050232082579, 11.70280126708063406,
                                        0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.021612921838043420, 8.541090340482851317, 10.33134080113699937, 11.70280164927028110,
                                        0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.014873257214348179, 8.372359455630595804, 9.916182195246349806, 11.17180856454275784,
                                        0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.009231913341891662, 8.231125581765601495, 9.568679376010450710, 10.72734808288726427,
                                        0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.004257781943209693, 8.106595706893173769, 9.262276376411159973, 10.33545478808742679,
                                        0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.000598733449153954, 8.014989592258231710, 9.036881559829700095, 10.04717197944511220,
                                        0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.000598734163229864, 8.014989610135470954, 9.036881603816251385, 10.04717203570446493,
                                        0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
                                        0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};

    // The expected (y-) points, generated by interactor.
    std::vector<double> expectedPointsY{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9907980331007956165, 0.7696241424499061790, 0.4331653038203038264, 0.2750112858449960118,
                                        2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.990363764459876617, 1.758752009170715835, 1.406414660633680480, 1.240796876357717071,
                                        3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.009636224659985260, 3.241247718439566050, 3.593584669157450229, 3.759202266436699524,
                                        4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.009201984124368678, 4.230376288790720807, 4.566835757237710958, 4.724990071261565561,
                                        5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.005581416323523491, 5.139733557614412973, 5.343811324323562317, 5.439739013178277460,
                                        6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.001994298641023740, 6.049928266931340559, 6.122847395199174869, 6.157123365388534886,
                                        7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.000114079035108183, 7.002856025872458368, 7.007027188416842911, 7.008987862473465391,
                                        8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 7.999885924712966556, 7.997144067962484115, 7.992973042461954059, 7.991012432823422884,
                                        9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0,
                                        10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0};

    // Set up the spline

    std::vector<double> splineXCoordinates({11.0000, 11.5878, 11.9511, 11.9511, 11.5878, 11.0000, 10.4122, 10.0489, 10.0489, 10.4122, 11.0000});
    std::vector<double> splineYCoordinates({0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0});

    // Not part of the unit test but added to ensure arrays are the same size
    ASSERT_EQ(splineXCoordinates.size(), splineYCoordinates.size());

    meshkernelapi::GeometryList spline;
    spline.geometry_separator = constants::missing::doubleValue;
    spline.inner_outer_separator = constants::missing::innerOuterSeparator;
    spline.num_coordinates = static_cast<int>(splineXCoordinates.size());
    spline.coordinates_x = splineXCoordinates.data();
    spline.coordinates_y = splineYCoordinates.data();

    double sectionControlPoint1x = 10.0;
    double sectionControlPoint1y = 8.0;
    double sectionControlPoint2x = 10.0;
    double sectionControlPoint2y = 1.0;

    errorCode = mkernel_curvilinear_snap_to_spline(meshKernelId, spline,
                                                   sectionControlPoint1x,
                                                   sectionControlPoint1y,
                                                   sectionControlPoint2x,
                                                   sectionControlPoint2y);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    for (int i = 0; i < curvilinearGrid.num_m * curvilinearGrid.num_n; ++i)
    {
        EXPECT_NEAR(curvilinearGrid.node_x[i], expectedPointsX[i], tolerance);
        EXPECT_NEAR(curvilinearGrid.node_y[i], expectedPointsY[i], tolerance);
    }

    // Now undo snapping.

    bool didUndoOfDeleteNode = false;
    int undoId = meshkernel::constants::missing::intValue;

    errorCode = meshkernelapi::mkernel_undo_state(didUndoOfDeleteNode, undoId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    EXPECT_TRUE(didUndoOfDeleteNode);
    ASSERT_EQ(meshKernelId, undoId);

    errorCode = mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    for (int i = 0; i < curvilinearGrid.num_m * curvilinearGrid.num_n; ++i)
    {
        EXPECT_NEAR(curvilinearGrid.node_x[i], originalNodeX[i], tolerance);
        EXPECT_NEAR(curvilinearGrid.node_y[i], originalNodeY[i], tolerance);
    }
}

class CurvilineartBoundariesAsPolygonsTests : public testing::TestWithParam<std::tuple<int, int, int, int, std::vector<double>, std::vector<double>>>
{
public:
    [[nodiscard]] static std::vector<std::tuple<int, int, int, int, std::vector<double>, std::vector<double>>> GetData()
    {
        return {
            std::make_tuple<int, int, int, int, std::vector<double>, std::vector<double>>(0, 0, 9, 9, std::vector<double>{0, 10, 20, 30}, std::vector<double>{0, 0, 0, 0}),
            std::make_tuple<int, int, int, int, std::vector<double>, std::vector<double>>(1, 1, 8, 8, std::vector<double>{10, 20, 30, 40}, std::vector<double>{10, 10, 10, 10})};
    }
};

TEST_P(CurvilineartBoundariesAsPolygonsTests, GetLocationIndex_OnACurvilinearGrid_ShouldGetTheLocationIndex)
{
    // Prepare
    auto const& [lowerLeftN, lowerLeftM, upperRightN, upperRightM, ValidCoordinatesX, ValidCoordinatesY] = GetParam();

    // Prepare
    int meshKernelId;
    auto errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    MakeGridParameters makeGridParameters;

    makeGridParameters.num_columns = 10;
    makeGridParameters.num_rows = 10;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = 10.0;
    makeGridParameters.block_size_y = 10.0;

    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    int numberOfPolygonNodes;
    errorCode = meshkernelapi::mkernel_curvilinear_count_boundaries_as_polygons(meshKernelId, lowerLeftN, lowerLeftM, upperRightN, upperRightM, numberOfPolygonNodes);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::GeometryList boundaryPolygon;
    std::vector<double> coordinates_x(numberOfPolygonNodes);
    std::vector<double> coordinates_y(numberOfPolygonNodes);

    boundaryPolygon.coordinates_x = coordinates_x.data();
    boundaryPolygon.coordinates_y = coordinates_y.data();
    boundaryPolygon.num_coordinates = numberOfPolygonNodes;
    boundaryPolygon.geometry_separator = constants::missing::doubleValue;

    errorCode = mkernel_curvilinear_get_boundaries_as_polygons(meshKernelId, lowerLeftN, lowerLeftM, upperRightN, upperRightM, boundaryPolygon);

    std::vector firstFourCoordinatesX(coordinates_x.begin(), coordinates_x.begin() + 4);
    ASSERT_THAT(firstFourCoordinatesX, ::testing::ContainerEq(ValidCoordinatesX));

    std::vector firstFourCoordinatesY(coordinates_y.begin(), coordinates_y.begin() + 4);
    ASSERT_THAT(firstFourCoordinatesY, ::testing::ContainerEq(ValidCoordinatesY));

    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}

INSTANTIATE_TEST_SUITE_P(CurvilineartBoundariesAsPolygonsTests, CurvilineartBoundariesAsPolygonsTests, ::testing::ValuesIn(CurvilineartBoundariesAsPolygonsTests::GetData()));

TEST(CurvilinearGrid, MakeCircularGrid_CartesianCoordinate_ShouldMakeCurvilinearGrid)
{

    meshkernel::MakeGridParameters parameters = {.num_columns = 14,
                                                 .num_rows = 10,
                                                 .angle = 32.0,
                                                 .origin_x = 17.0,
                                                 .origin_y = -23.0,
                                                 .radius_curvature = 10.0,
                                                 .uniform_columns_fraction = 1.0,
                                                 .uniform_rows_fraction = 0.0};

    int meshKernelId = 0;
    int projectionType = 0;
    auto errorCode = meshkernelapi::mkernel_allocate_state(projectionType, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Execute
    errorCode = meshkernelapi::mkernel_curvilinear_compute_circular_grid(meshKernelId, parameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGridResults;
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGridResults);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert
    ASSERT_EQ(parameters.num_columns + 1, curvilinearGridResults.num_m);
    ASSERT_EQ(parameters.num_rows + 1, curvilinearGridResults.num_n);

    std::vector<double> radiusValues{20.0, 33.0766048601183, 50.1763643268853,
                                     72.5370441018832, 101.777221484012, 140.013446050598,
                                     190.013446050598, 255.39647035119, 340.895267685025,
                                     452.698666560014, 598.899553470658};

    std::vector<double> thetaValues{2.12930168743308, 1.68050273692025, 1.23170378640743, 0.782904835894599,
                                    0.334105885381772, -0.114693065131056, -0.563492015643884, -1.01229096615671,
                                    -1.46108991666954, -1.90988886718237, -2.35868781769519, -2.80748676820802,
                                    -3.25628571872085, -3.70508466923368, -4.1538836197465};

    std::vector<double> xValues(curvilinearGridResults.num_m * curvilinearGridResults.num_n);
    std::vector<double> yValues(curvilinearGridResults.num_m * curvilinearGridResults.num_n);

    meshkernelapi::CurvilinearGrid gridData{};
    gridData.node_x = xValues.data();
    gridData.node_y = yValues.data();
    gridData.num_m = parameters.num_columns + 1;
    gridData.num_n = parameters.num_rows + 1;

    errorCode = mkernel_curvilinear_get_data(meshKernelId, gridData);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    const double tolerance = 1.0e-10;

    size_t count = 0;

    for (meshkernel::UInt i = 0; i < static_cast<meshkernel::UInt>(parameters.num_rows + 1); ++i)
    {
        for (meshkernel::UInt j = 0; j < static_cast<meshkernel::UInt>(parameters.num_columns + 1); ++j)
        {
            double expectedX = parameters.origin_x + radiusValues[i] * std::cos(thetaValues[j]);
            double expectedY = parameters.origin_y + radiusValues[i] * std::sin(thetaValues[j]);
            EXPECT_NEAR(expectedX, xValues[count], tolerance);
            EXPECT_NEAR(expectedY, yValues[count], tolerance);
            ++count;
        }
    }

    errorCode = meshkernelapi::mkernel_deallocate_state(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}

TEST(CurvilinearGrid, MakeCircularGrid_CartesianCoordinate_ShouldFail)
{
    int meshKernelId = 0;
    int projectionType = 0;
    auto errorCode = meshkernelapi::mkernel_allocate_state(projectionType, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::MakeGridParameters parameters = {.num_columns = 14,
                                                 .num_rows = 10,
                                                 .angle = 32.0,
                                                 .origin_x = 17.0,
                                                 .origin_y = -23.0,
                                                 .radius_curvature = 10.0,
                                                 .uniform_columns_fraction = -1.0,
                                                 .uniform_rows_fraction = 0.0};

    errorCode = meshkernelapi::mkernel_curvilinear_compute_circular_grid(meshKernelId, parameters);
    ASSERT_EQ(meshkernel::ExitCode::RangeErrorCode, errorCode);
    parameters.uniform_columns_fraction = 1.0;

    parameters.uniform_rows_fraction = -1.0;
    errorCode = meshkernelapi::mkernel_curvilinear_compute_circular_grid(meshKernelId, parameters);
    ASSERT_EQ(meshkernel::ExitCode::RangeErrorCode, errorCode);
    parameters.uniform_rows_fraction = 1.0;

    parameters.maximum_uniform_size_columns = -1.0;
    errorCode = meshkernelapi::mkernel_curvilinear_compute_circular_grid(meshKernelId, parameters);
    ASSERT_EQ(meshkernel::ExitCode::RangeErrorCode, errorCode);
    parameters.maximum_uniform_size_columns = 1.0;

    parameters.maximum_uniform_size_rows = -1.0;
    errorCode = meshkernelapi::mkernel_curvilinear_compute_circular_grid(meshKernelId, parameters);
    ASSERT_EQ(meshkernel::ExitCode::RangeErrorCode, errorCode);
    parameters.maximum_uniform_size_rows = 1.0;

    errorCode = meshkernelapi::mkernel_curvilinear_compute_circular_grid(meshKernelId + 100, parameters);
    ASSERT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_compute_circular_grid(meshKernelId, parameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Now try to generate another grid with the same mesh kernel id
    errorCode = meshkernelapi::mkernel_curvilinear_compute_circular_grid(meshKernelId, parameters);
    ASSERT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    errorCode = meshkernelapi::mkernel_deallocate_state(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}
