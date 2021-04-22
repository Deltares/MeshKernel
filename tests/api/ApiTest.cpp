#include <exception>
#include <memory>
#include <numeric>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <MeshKernelApi/CurvilinearGrid.hpp>
#include <MeshKernelApi/GeometryList.hpp>
#include <MeshKernelApi/MakeMeshParameters.hpp>
#include <MeshKernelApi/Mesh1D.hpp>
#include <MeshKernelApi/Mesh2D.hpp>
#include <MeshKernelApi/MeshKernel.hpp>

#include <TestUtils/Definitions.hpp>
#include <TestUtils/MakeMeshes.hpp>

class ApiTests : public testing::Test
{
public:
    /// Constructor for allocating state
    ApiTests()
    {
        int isGeographic = 0;
        const auto errorCode = meshkernelapi::mkernel_allocate_state(isGeographic, m_meshKernelId);
        if (errorCode != 0)
        {
            throw std::runtime_error("Could not allocate state");
        }
    }

    /// Destructor for deallocating state
    ~ApiTests()
    {
        meshkernelapi::mkernel_deallocate_state(m_meshKernelId);
    }

    /// @brief Makes a mesh
    /// @param[in]  n            Number of rows
    /// @param[in]  m            Number of columns
    /// @param[in]  delta        Distance between neighboring nodes
    void MakeMesh(int n = 4, int m = 3, double delta = 1.0)
    {
        // Set-up new mesh
        const auto mesh2d = MakeRectangularMeshForApiTesting(n, m, delta);
        const auto errorCode = mkernel_set_mesh2d(m_meshKernelId, mesh2d);
        if (errorCode != 0)
        {
            throw std::runtime_error("Could not set mesh2d");
        }
        // Clean up client-allocated memory
        DeleteRectangularMeshForApiTesting(mesh2d);
    }

    void MakeUniformCurvilinearGrid(int numberOfColumns = 4, int numberOfRows = 4, double blockSize = 10.0)
    {

        meshkernelapi::MakeMeshParameters makeMeshParameters{};
        meshkernelapi::GeometryList geometryList{};

        makeMeshParameters.grid_type = 0;
        makeMeshParameters.num_columns = numberOfColumns;
        makeMeshParameters.num_rows = numberOfRows;
        makeMeshParameters.angle = 0.0;
        makeMeshParameters.block_size = 0.0;
        makeMeshParameters.origin_x = 0.0;
        makeMeshParameters.origin_y = 0.0;
        makeMeshParameters.block_size_x = blockSize;
        makeMeshParameters.block_size_y = blockSize;

        auto errorCode = mkernel_make_uniform_curvilinear(m_meshKernelId, makeMeshParameters, geometryList);
        if (errorCode != 0)
        {
            throw std::runtime_error("Could not create uniform curvilinear grid");
        }
    }

    [[nodiscard]] int GetMeshKernelId() const
    {
        return m_meshKernelId;
    }

private:
    int m_meshKernelId;
};

TEST_F(ApiTests, DeleteNodeThroughApi)
{
    // Prepare
    MakeMesh();
    auto meshKernelId = GetMeshKernelId();

    // Execute
    auto errorCode = meshkernelapi::mkernel_delete_node_mesh2d(meshKernelId, 0);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the dimensions
    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_get_dimensions_mesh2d(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert dimensions
    ASSERT_EQ(11, mesh2d.num_nodes);
    ASSERT_EQ(15, mesh2d.num_edges);

    // Allocate memory and get data
    std::vector<int> edge_nodes(mesh2d.num_edges * 2);
    std::vector<int> face_nodes(mesh2d.num_face_nodes);
    std::vector<int> nodes_per_face(mesh2d.num_faces);
    std::vector<double> node_x(mesh2d.num_nodes);
    std::vector<double> node_y(mesh2d.num_nodes);
    std::vector<double> edge_x(mesh2d.num_edges);
    std::vector<double> edge_y(mesh2d.num_edges);
    std::vector<double> face_x(mesh2d.num_faces);
    std::vector<double> face_y(mesh2d.num_faces);

    mesh2d.edge_nodes = &edge_nodes[0];
    mesh2d.face_nodes = &face_nodes[0];
    mesh2d.nodes_per_face = &nodes_per_face[0];
    mesh2d.node_x = &node_x[0];
    mesh2d.node_y = &node_y[0];
    mesh2d.edge_x = &edge_x[0];
    mesh2d.edge_y = &edge_y[0];
    mesh2d.face_x = &face_x[0];
    mesh2d.face_y = &face_y[0];
    errorCode = mkernel_get_data_mesh2d(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

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
    ASSERT_EQ(0, mesh2d.face_nodes[0]);
    ASSERT_EQ(3, mesh2d.face_nodes[1]);
    ASSERT_EQ(4, mesh2d.face_nodes[2]);
    ASSERT_EQ(1, mesh2d.face_nodes[3]);
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

TEST_F(ApiTests, FlipEdgesThroughApi)
{
    // Prepare
    MakeMesh();
    auto meshKernelId = GetMeshKernelId();

    // Execute
    const int isTriangulationRequired = 1;
    const int projectToLandBoundaryOption = 1;
    auto errorCode = meshkernelapi::mkernel_flip_edges_mesh2d(meshKernelId, isTriangulationRequired, projectToLandBoundaryOption);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_get_dimensions_mesh2d(meshKernelId, mesh2d);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(12, mesh2d.num_nodes);
    ASSERT_EQ(23, mesh2d.num_edges);
}

TEST_F(ApiTests, InsertEdgeThroughApi)
{
    // Prepare
    MakeMesh();
    auto meshKernelId = GetMeshKernelId();

    // Execute
    int newEdgeIndex;
    auto errorCode = meshkernelapi::mkernel_insert_edge_mesh2d(meshKernelId, 0, 4, newEdgeIndex);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(17, newEdgeIndex);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_get_dimensions_mesh2d(meshKernelId, mesh2d);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(12, mesh2d.num_nodes);
    ASSERT_EQ(18, mesh2d.num_edges);
}

TEST_F(ApiTests, MergeTwoNodesThroughApi)
{
    // Prepare
    MakeMesh();
    auto meshKernelId = GetMeshKernelId();

    // Execute
    auto errorCode = meshkernelapi::mkernel_merge_two_nodes_mesh2d(meshKernelId, 0, 4);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_get_dimensions_mesh2d(meshKernelId, mesh2d);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(11, mesh2d.num_nodes);
    ASSERT_EQ(15, mesh2d.num_edges);
}

TEST_F(ApiTests, MergeNodesThroughApi)
{
    // Prepare
    MakeMesh();
    auto meshKernelId = GetMeshKernelId();
    meshkernelapi::GeometryList geometry_list{};

    // Execute
    auto errorCode = mkernel_merge_nodes_mesh2d(meshKernelId, geometry_list);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_get_dimensions_mesh2d(meshKernelId, mesh2d);

    // Assert (nothing is done, just check that the api communication works)
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(12, mesh2d.num_nodes);
    ASSERT_EQ(17, mesh2d.num_edges);
}

TEST_F(ApiTests, OrthogonalizationThroughApi)
{
    // Set a new mesh in mesh
    MakeMesh();
    auto meshKernelId = GetMeshKernelId();

    // Prepare
    meshkernelapi::OrthogonalizationParameters orthogonalizationParameters{};
    orthogonalizationParameters.outer_iterations = 1;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;

    meshkernelapi::GeometryList geometryList{};
    meshkernelapi::GeometryList landBoundaries{};

    // Execute
    auto errorCode = mkernel_initialize_orthogonalization_mesh2d(meshKernelId,
                                                                 1,
                                                                 orthogonalizationParameters,
                                                                 landBoundaries,
                                                                 geometryList);

    // Assert (nothing is done, just check that the api communication works)
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    errorCode = meshkernelapi::mkernel_prepare_outer_iteration_orthogonalization_mesh2d(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    errorCode = meshkernelapi::mkernel_compute_inner_ortogonalization_iteration_mesh2d(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    errorCode = meshkernelapi::mkernel_finalize_inner_ortogonalization_iteration_mesh2d(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    errorCode = meshkernelapi::mkernel_delete_orthogonalization_mesh2d(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_get_dimensions_mesh2d(meshKernelId, mesh2d);

    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(12, mesh2d.num_nodes);
    ASSERT_EQ(17, mesh2d.num_edges);
}

TEST_F(ApiTests, GenerateTriangularGridThroughApi)
{
    // Prepare
    auto meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometry_separator = meshkernel::doubleMissingValue;
    std::unique_ptr<double> xCoordinates(new double[17]{
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
        415.319672});

    std::unique_ptr<double> yCoordinates(new double[17]{
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
        490.293762});

    std::unique_ptr<double> zCoordinates(new double[17]{
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
        0.0});

    geometryListIn.coordinates_x = xCoordinates.get();
    geometryListIn.coordinates_y = yCoordinates.get();
    geometryListIn.values = zCoordinates.get();
    geometryListIn.num_coordinates = 17;

    // Execute
    auto errorCode = mkernel_make_mesh_from_polygon_mesh2d(meshKernelId, geometryListIn);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    // Get the new state

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_get_dimensions_mesh2d(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(38, mesh2d.num_nodes);
    ASSERT_EQ(95, mesh2d.num_edges);
}

TEST_F(ApiTests, GenerateTriangularGridFromSamplesThroughApi)
{
    // Prepare
    auto meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;

    geometryListIn.geometry_separator = meshkernel::doubleMissingValue;

    std::unique_ptr<double> xCoordinates(new double[5]{
        0.0,
        10.0,
        10.0,
        0.0,
        0.0});

    std::unique_ptr<double> yCoordinates(new double[5]{
        0.0,
        0.0,
        10.0,
        10.0,
        0.0});

    std::unique_ptr<double> zCoordinates(new double[5]{
        0.0,
        0.0,
        0.0,
        0.0,
        0.0});

    geometryListIn.coordinates_x = xCoordinates.get();
    geometryListIn.coordinates_y = yCoordinates.get();
    geometryListIn.values = zCoordinates.get();

    geometryListIn.num_coordinates = 5;

    // Execute
    auto errorCode = mkernel_make_mesh_from_samples_mesh2d(meshKernelId, geometryListIn);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_get_dimensions_mesh2d(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(4, mesh2d.num_nodes);
    ASSERT_EQ(5, mesh2d.num_edges);
}

TEST_F(ApiTests, GetMeshBoundariesThroughApi)
{
    // Prepare
    MakeMesh();
    auto meshKernelId = GetMeshKernelId();
    int numberOfpolygonNodes;
    auto errorCode = meshkernelapi::mkernel_count_mesh_boundaries_to_polygon_mesh2d(meshKernelId, numberOfpolygonNodes);
    ASSERT_EQ(11, numberOfpolygonNodes);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::GeometryList geometryListOut;
    geometryListOut.geometry_separator = meshkernel::doubleMissingValue;
    geometryListOut.num_coordinates = numberOfpolygonNodes;

    std::unique_ptr<double> xCoordinates(new double[numberOfpolygonNodes]);
    std::unique_ptr<double> yCoordinates(new double[numberOfpolygonNodes]);
    std::unique_ptr<double> zCoordinates(new double[numberOfpolygonNodes]);

    geometryListOut.coordinates_x = xCoordinates.get();
    geometryListOut.coordinates_y = yCoordinates.get();
    geometryListOut.values = zCoordinates.get();

    // Execute
    errorCode = mkernel_get_mesh_boundaries_to_polygon_mesh2d(meshKernelId, geometryListOut);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    const double tolerance = 1e-6;
    ASSERT_NEAR(0.0, geometryListOut.coordinates_x[0], tolerance);
    ASSERT_NEAR(0.0, geometryListOut.coordinates_y[0], tolerance);
}

TEST_F(ApiTests, OffsetAPolygonThroughApi)
{
    // Prepare
    MakeMesh();
    auto meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometry_separator = meshkernel::doubleMissingValue;
    geometryListIn.num_coordinates = 4;

    std::unique_ptr<double> xCoordinatesIn(new double[4]{
        0.0,
        1.0,
        1.0,
        0.0});

    std::unique_ptr<double> yCoordinatesIn(new double[4]{
        0.0,
        0.0,
        1.0,
        1.0});

    std::unique_ptr<double> zCoordinatesIn(new double[4]{
        0.0,
        0.0,
        0.0,
        0.0});

    geometryListIn.coordinates_x = xCoordinatesIn.get();
    geometryListIn.coordinates_y = yCoordinatesIn.get();
    geometryListIn.values = zCoordinatesIn.get();

    // Execute
    int numberOfpolygonNodes;
    auto errorCode = mkernel_count_offset_polygon(meshKernelId, geometryListIn, false, 0.5, numberOfpolygonNodes);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(4, numberOfpolygonNodes);

    meshkernelapi::GeometryList geometryListOut;

    geometryListOut.num_coordinates = numberOfpolygonNodes;
    geometryListOut.geometry_separator = meshkernel::doubleMissingValue;

    std::unique_ptr<double> xCoordinatesOut(new double[numberOfpolygonNodes]);
    std::unique_ptr<double> yCoordinatesOut(new double[numberOfpolygonNodes]);
    std::unique_ptr<double> zCoordinatesOut(new double[numberOfpolygonNodes]);
    geometryListOut.coordinates_x = xCoordinatesOut.get();
    geometryListOut.coordinates_y = yCoordinatesOut.get();
    geometryListOut.values = zCoordinatesOut.get();
    errorCode = mkernel_get_offset_polygon(meshKernelId, geometryListIn, false, 10.0, geometryListOut);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    const double tolerance = 1e-6;
    ASSERT_NEAR(0.0, geometryListOut.coordinates_x[0], tolerance);
    ASSERT_NEAR(-10.0, geometryListOut.coordinates_y[0], tolerance);
}

TEST_F(ApiTests, RefineAPolygonThroughApi)
{
    // Prepare
    MakeMesh();
    auto meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometry_separator = meshkernel::doubleMissingValue;
    geometryListIn.num_coordinates = 3;
    std::unique_ptr<double> xCoordinatesIn(new double[3]{
        76.251099,
        498.503723,
        505.253784});

    std::unique_ptr<double> yCoordinatesIn(new double[3]{
        92.626556,
        91.126541,
        490.130554});

    std::unique_ptr<double> zCoordinatesIn(new double[3]{
        0.0,
        0.0,
        0.0});

    geometryListIn.coordinates_x = xCoordinatesIn.get();
    geometryListIn.coordinates_y = yCoordinatesIn.get();
    geometryListIn.values = zCoordinatesIn.get();

    // Execute
    int numberOfpolygonNodes;
    auto errorCode = mkernel_count_refine_polygon(meshKernelId, geometryListIn, 0, 2, 40, numberOfpolygonNodes);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(22, numberOfpolygonNodes);

    meshkernelapi::GeometryList geometryListOut;
    geometryListOut.num_coordinates = numberOfpolygonNodes;
    geometryListOut.geometry_separator = meshkernel::doubleMissingValue;
    std::unique_ptr<double> xCoordinatesOut(new double[numberOfpolygonNodes]);
    std::unique_ptr<double> yCoordinatesOut(new double[numberOfpolygonNodes]);
    std::unique_ptr<double> zCoordinatesOut(new double[numberOfpolygonNodes]);
    geometryListOut.coordinates_x = xCoordinatesOut.get();
    geometryListOut.coordinates_y = yCoordinatesOut.get();
    geometryListOut.values = zCoordinatesOut.get();
    errorCode = mkernel_refine_polygon(meshKernelId, geometryListIn, false, 0, 2, geometryListOut);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    const double tolerance = 1e-6;
    ASSERT_NEAR(76.251099, geometryListOut.coordinates_x[0], tolerance);
    ASSERT_NEAR(92.626556, geometryListOut.coordinates_y[0], tolerance);
}

TEST_F(ApiTests, RefineAGridBasedOnSamplesThroughApi)
{
    // Prepare
    MakeMesh();
    auto meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometry_separator = meshkernel::doubleMissingValue;
    std::unique_ptr<double> xCoordinatesIn(new double[9]{
        50.0,
        150.0,
        250.0,
        50.0,
        150.0,
        250.0,
        50.0,
        150.0,
        250.0});

    std::unique_ptr<double> yCoordinatesIn(new double[9]{
        50.0,
        50.0,
        50.0,
        150.0,
        150.0,
        150.0,
        250.0,
        250.0,
        250.0});

    std::unique_ptr<double> zCoordinatesIn(new double[9]{
        2.0,
        2.0,
        2.0,
        3.0,
        3.0,
        3.0,
        4.0,
        4.0,
        4.0});

    geometryListIn.coordinates_x = xCoordinatesIn.get();
    geometryListIn.coordinates_y = yCoordinatesIn.get();
    geometryListIn.values = zCoordinatesIn.get();

    geometryListIn.num_coordinates = 9;

    meshkernelapi::InterpolationParameters interpolationParameters;
    interpolationParameters.max_num_refinement_iterations = 2;
    interpolationParameters.averaging_method = 1;
    interpolationParameters.minimum_num_points = 1;
    interpolationParameters.relative_search_radius = 1.01;
    interpolationParameters.interpolate_to = 3;
    interpolationParameters.refine_intersected = 0;

    meshkernelapi::SampleRefineParameters samplesRefineParameters;
    samplesRefineParameters.min_face_size = 0.5;
    samplesRefineParameters.refinement_type = 3;
    samplesRefineParameters.connect_hanging_nodes = 1;
    samplesRefineParameters.maximum_time_step = 0.0;
    samplesRefineParameters.account_for_samples_outside = 0;

    // Execute
    auto errorCode = mkernel_refine_based_on_samples_mesh2d(meshKernelId, geometryListIn, interpolationParameters, samplesRefineParameters);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_get_dimensions_mesh2d(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(12, mesh2d.num_nodes);
    ASSERT_EQ(17, mesh2d.num_edges);
}

TEST_F(ApiTests, RefineAGridBasedOnPolygonThroughApi)
{
    // Prepare
    MakeMesh();
    auto meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometry_separator = meshkernel::doubleMissingValue;
    std::unique_ptr<double> xCoordinatesIn(new double[9]{
        50.0,
        150.0,
        250.0,
        50.0,
        150.0,
        250.0,
        50.0,
        150.0,
        250.0});

    std::unique_ptr<double> yCoordinatesIn(new double[9]{
        50.0,
        50.0,
        50.0,
        150.0,
        150.0,
        150.0,
        250.0,
        250.0,
        250.0});

    std::unique_ptr<double> zCoordinatesIn(new double[9]{
        2.0,
        2.0,
        2.0,
        3.0,
        3.0,
        3.0,
        4.0,
        4.0,
        4.0});

    geometryListIn.coordinates_x = xCoordinatesIn.get();
    geometryListIn.coordinates_y = yCoordinatesIn.get();
    geometryListIn.values = zCoordinatesIn.get();

    geometryListIn.num_coordinates = 9;

    meshkernelapi::InterpolationParameters interpolationParameters;
    interpolationParameters.max_num_refinement_iterations = 2;
    interpolationParameters.averaging_method = 1;
    interpolationParameters.minimum_num_points = 1;
    interpolationParameters.relative_search_radius = 1.01;
    interpolationParameters.interpolate_to = 3;
    interpolationParameters.refine_intersected = 0;

    // Execute
    auto errorCode = mkernel_refine_based_on_polygon_mesh2d(meshKernelId, geometryListIn, interpolationParameters);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_get_dimensions_mesh2d(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(12, mesh2d.num_nodes);
    ASSERT_EQ(17, mesh2d.num_edges);
}

TEST_F(ApiTests, ComputeSingleContactsThroughApi)
{
    // Prepare
    MakeMesh(4, 4, 10);
    auto meshKernelId = GetMeshKernelId();

    // Init 1d mesh
    meshkernelapi::Mesh1D mesh1d;
    std::unique_ptr<double> node_x(new double[7]{
        1.73493900000000,
        2.35659313023165,
        5.38347452702839,
        14.2980910429074,
        22.9324017677239,
        25.3723169493137,
        25.8072280000000});
    std::unique_ptr<double> node_y(new double[7]{
        -7.6626510000000,
        1.67281447902331,
        10.3513746546384,
        12.4797224193970,
        15.3007317677239,
        24.1623588554512,
        33.5111870000000});
    mesh1d.node_x = node_x.get();
    mesh1d.node_y = node_y.get();
    mesh1d.num_nodes = 7;

    std::unique_ptr<int> edge_nodes(new int[12]{
        0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6});
    mesh1d.edge_nodes = edge_nodes.get();
    mesh1d.num_edges = 6;

    auto errorCode = mkernel_set_mesh1d(meshKernelId, mesh1d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Init 1d mask
    std::unique_ptr<int> onedNodeMask(new int[7]{1, 1, 1, 1, 1, 1, 1});

    // Init polygon
    meshkernelapi::GeometryList polygon;
    polygon.geometry_separator = meshkernel::doubleMissingValue;

    std::unique_ptr<double> xCoordinates(new double[5]{-30, 40, 40, -40, -30});
    std::unique_ptr<double> yCoordinates(new double[5]{-20, -20, 50, 50, -20});
    std::unique_ptr<double> zCoordinates(new double[5]{0, 0, 0, 0, 0});
    polygon.coordinates_x = xCoordinates.get();
    polygon.coordinates_y = yCoordinates.get();
    polygon.values = zCoordinates.get();

    polygon.num_coordinates = 5;

    // Execute
    errorCode = mkernel_compute_single_contacts(meshKernelId, onedNodeMask.get(), polygon);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state
    meshkernelapi::Contacts contacts{};
    errorCode = mkernel_get_dimensions_contacts(meshKernelId, contacts);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    std::unique_ptr<int> mesh1d_indices(new int[contacts.num_contacts]);
    std::unique_ptr<int> mesh2d_indices(new int[contacts.num_contacts]);
    contacts.mesh1d_indices = mesh1d_indices.get();
    contacts.mesh2d_indices = mesh2d_indices.get();

    errorCode = mkernel_get_data_contacts(meshKernelId, contacts);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

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

TEST_F(ApiTests, ComputeMultipleContactsThroughApi)
{
    // Prepare
    MakeMesh(4, 4, 10);
    auto meshKernelId = GetMeshKernelId();

    // Init 1d mesh
    meshkernelapi::Mesh1D mesh1d;
    std::unique_ptr<double> node_x(new double[7]{
        1.73493900000000,
        2.35659313023165,
        5.38347452702839,
        14.2980910429074,
        22.9324017677239,
        25.3723169493137,
        25.8072280000000});
    std::unique_ptr<double> node_y(new double[7]{
        -7.6626510000000,
        1.67281447902331,
        10.3513746546384,
        12.4797224193970,
        15.3007317677239,
        24.1623588554512,
        33.5111870000000});
    mesh1d.node_x = node_x.get();
    mesh1d.node_y = node_y.get();
    mesh1d.num_nodes = 7;

    std::unique_ptr<int> edge_nodes(new int[12]{
        0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6});
    mesh1d.edge_nodes = edge_nodes.get();
    mesh1d.num_edges = 6;

    auto errorCode = mkernel_set_mesh1d(meshKernelId, mesh1d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Init 1d mask
    std::unique_ptr<int> onedNodeMask(new int[7]{
        1, 1, 1, 1, 1, 1, 1});

    // Execute
    // Why do we need to specify the namespace "meshkernelapi"?
    errorCode = meshkernelapi::mkernel_compute_multiple_contacts(meshKernelId, onedNodeMask.get());
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state
    meshkernelapi::Contacts contacts{};
    errorCode = mkernel_get_dimensions_contacts(meshKernelId, contacts);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    std::unique_ptr<int> mesh1d_indices(new int[contacts.num_contacts]);
    std::unique_ptr<int> mesh2d_indices(new int[contacts.num_contacts]);
    contacts.mesh1d_indices = mesh1d_indices.get();
    contacts.mesh2d_indices = mesh2d_indices.get();

    errorCode = mkernel_get_data_contacts(meshKernelId, contacts);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

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

TEST_F(ApiTests, ComputeContactsWithPolygonsThroughApi)
{
    // Prepare
    MakeMesh(4, 4, 10);
    auto meshKernelId = GetMeshKernelId();

    // Init 1d mesh
    meshkernelapi::Mesh1D mesh1d;
    std::unique_ptr<double> node_x(new double[7]{
        1.73493900000000,
        2.35659313023165,
        5.38347452702839,
        14.2980910429074,
        22.9324017677239,
        25.3723169493137,
        25.8072280000000});
    std::unique_ptr<double> node_y(new double[7]{
        -7.6626510000000,
        1.67281447902331,
        10.3513746546384,
        12.4797224193970,
        15.3007317677239,
        24.1623588554512,
        33.5111870000000});
    mesh1d.node_x = node_x.get();
    mesh1d.node_y = node_y.get();
    mesh1d.num_nodes = 7;

    std::unique_ptr<int> edge_nodes(new int[12]{
        0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6});
    mesh1d.edge_nodes = edge_nodes.get();
    mesh1d.num_edges = 6;

    auto errorCode = mkernel_set_mesh1d(meshKernelId, mesh1d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Init 1d mask
    std::unique_ptr<int> onedNodeMask(new int[7]{
        1, 1, 1, 1, 1, 1, 1});

    // Init polygon
    meshkernelapi::GeometryList polygon;
    polygon.geometry_separator = meshkernel::doubleMissingValue;

    std::unique_ptr<double> xCoordinates(new double[5]{
        25,
        50,
        50,
        25,
        25});
    std::unique_ptr<double> yCoordinates(new double[5]{
        25,
        25,
        50,
        50,
        25});
    std::unique_ptr<double> zCoordinates(new double[5]{
        0, 0, 0, 0, 0});
    polygon.coordinates_x = xCoordinates.get();
    polygon.coordinates_y = yCoordinates.get();
    polygon.values = zCoordinates.get();

    polygon.num_coordinates = 5;

    // Execute
    errorCode = mkernel_compute_with_polygons_contacts(meshKernelId, onedNodeMask.get(), polygon);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state
    meshkernelapi::Contacts contacts{};
    errorCode = mkernel_get_dimensions_contacts(meshKernelId, contacts);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    std::unique_ptr<int> mesh1d_indices(new int[contacts.num_contacts]);
    std::unique_ptr<int> mesh2d_indices(new int[contacts.num_contacts]);
    contacts.mesh1d_indices = mesh1d_indices.get();
    contacts.mesh2d_indices = mesh2d_indices.get();

    errorCode = mkernel_get_data_contacts(meshKernelId, contacts);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(1, contacts.num_contacts);
    ASSERT_EQ(5, contacts.mesh1d_indices[0]);
    ASSERT_EQ(8, contacts.mesh2d_indices[0]);
}

TEST_F(ApiTests, ComputeContactsWithPointsThroughApi)
{
    // Prepare
    MakeMesh(4, 4, 10);
    auto meshKernelId = GetMeshKernelId();

    // Init 1d mesh
    meshkernelapi::Mesh1D mesh1d;
    std::unique_ptr<double> node_x(new double[7]{
        1.73493900000000,
        2.35659313023165,
        5.38347452702839,
        14.2980910429074,
        22.9324017677239,
        25.3723169493137,
        25.8072280000000});
    std::unique_ptr<double> node_y(new double[7]{
        -7.6626510000000,
        1.67281447902331,
        10.3513746546384,
        12.4797224193970,
        15.3007317677239,
        24.1623588554512,
        33.5111870000000});
    mesh1d.node_x = node_x.get();
    mesh1d.node_y = node_y.get();
    mesh1d.num_nodes = 7;

    std::unique_ptr<int> edge_nodes(new int[12]{
        0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6});
    mesh1d.edge_nodes = edge_nodes.get();
    mesh1d.num_edges = 6;

    auto errorCode = mkernel_set_mesh1d(meshKernelId, mesh1d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Init 1d mask
    std::unique_ptr<int> onedNodeMask(new int[7]{
        1, 1, 1, 1, 1, 1, 1});

    // Init polygon
    meshkernelapi::GeometryList points;
    points.geometry_separator = meshkernel::doubleMissingValue;

    std::unique_ptr<double> xCoordinates(new double[4]{
        5,
        15,
        25,
        35});
    std::unique_ptr<double> yCoordinates(new double[4]{
        5,
        15,
        25,
        35});
    std::unique_ptr<double> zCoordinates(new double[4]{
        0, 0, 0, 0});
    points.coordinates_x = xCoordinates.get();
    points.coordinates_y = yCoordinates.get();
    points.values = zCoordinates.get();

    points.num_coordinates = 4;

    // Execute
    errorCode = mkernel_compute_with_points_contacts(meshKernelId, onedNodeMask.get(), points);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state
    meshkernelapi::Contacts contacts{};
    errorCode = mkernel_get_dimensions_contacts(meshKernelId, contacts);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    std::unique_ptr<int> mesh1d_indices(new int[contacts.num_contacts]);
    std::unique_ptr<int> mesh2d_indices(new int[contacts.num_contacts]);
    contacts.mesh1d_indices = mesh1d_indices.get();
    contacts.mesh2d_indices = mesh2d_indices.get();

    errorCode = mkernel_get_data_contacts(meshKernelId, contacts);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(3, contacts.num_contacts);

    ASSERT_EQ(1, contacts.mesh1d_indices[0]);
    ASSERT_EQ(3, contacts.mesh1d_indices[1]);
    ASSERT_EQ(5, contacts.mesh1d_indices[2]);

    ASSERT_EQ(0, contacts.mesh2d_indices[0]);
    ASSERT_EQ(4, contacts.mesh2d_indices[1]);
    ASSERT_EQ(8, contacts.mesh2d_indices[2]);
}

TEST_F(ApiTests, ComputeBoundaryContactsThroughApi)
{
    // Prepare
    MakeMesh(4, 4, 10);
    auto meshKernelId = GetMeshKernelId();

    // Init 1d mesh
    meshkernelapi::Mesh1D mesh1d;
    std::unique_ptr<double> node_x(new double[8]{
        -16.1886410000000,
        -16.1464995876014,
        -16.1043581752028,
        -16.0622167628042,
        -15.7539488236928,
        -6.86476658679268,
        2.02441565010741,
        10.9135970000000,
    });
    std::unique_ptr<double> node_y(new double[8]{
        0.89018900000000,
        9.78201442138723,
        18.6738398427745,
        27.5656652641617,
        36.1966603330179,
        36.4175095626911,
        36.6383587923643,
        36.8592080000000});
    mesh1d.node_x = node_x.get();
    mesh1d.node_y = node_y.get();
    mesh1d.num_nodes = 8;

    std::unique_ptr<int> edge_nodes(new int[14]{
        0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7});
    mesh1d.edge_nodes = edge_nodes.get();
    mesh1d.num_edges = 7;

    auto errorCode = mkernel_set_mesh1d(meshKernelId, mesh1d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Init 1d mask
    std::unique_ptr<int> onedNodeMask(new int[8]{
        1, 1, 1, 1, 1, 1, 1, 1});

    // Init polygon
    meshkernelapi::GeometryList polygon;
    polygon.geometry_separator = meshkernel::doubleMissingValue;

    std::unique_ptr<double> xCoordinates(new double[5]{
        -30,
        40,
        40,
        -40,
        -30});
    std::unique_ptr<double> yCoordinates(new double[5]{
        -20,
        -20,
        50,
        50,
        -20});
    std::unique_ptr<double> zCoordinates(new double[5]{
        0, 0, 0, 0, 0});
    polygon.coordinates_x = xCoordinates.get();
    polygon.coordinates_y = yCoordinates.get();
    polygon.values = zCoordinates.get();

    polygon.num_coordinates = 5;

    // Execute
    errorCode = mkernel_compute_boundary_contacts(meshKernelId, onedNodeMask.get(), polygon, 200.0);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state
    meshkernelapi::Contacts contacts{};
    errorCode = mkernel_get_dimensions_contacts(meshKernelId, contacts);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    std::unique_ptr<int> mesh1d_indices(new int[contacts.num_contacts]);
    std::unique_ptr<int> mesh2d_indices(new int[contacts.num_contacts]);
    contacts.mesh1d_indices = mesh1d_indices.get();
    contacts.mesh2d_indices = mesh2d_indices.get();

    errorCode = mkernel_get_data_contacts(meshKernelId, contacts);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

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
    std::unique_ptr<double> xCoordinatesIn(new double[3]{10.0, 20.0, 30.0});
    std::unique_ptr<double> yCoordinatesIn(new double[3]{-5.0, 5.0, -5.0});
    std::unique_ptr<double> zCoordinatesIn(new double[3]{0.0, 0.0, 0.0});
    geometryListIn.coordinates_x = xCoordinatesIn.get();
    geometryListIn.coordinates_y = yCoordinatesIn.get();
    geometryListIn.values = zCoordinatesIn.get();
    geometryListIn.num_coordinates = 3;
    geometryListIn.geometry_separator = meshkernel::doubleMissingValue;

    meshkernelapi::GeometryList geometryListOut;
    int numberOfPointsBetweenNodes = 20;
    std::unique_ptr<double> xCoordinatesOut(new double[(numberOfPointsBetweenNodes + 1) * 2 + 1]);
    std::unique_ptr<double> yCoordinatesOut(new double[(numberOfPointsBetweenNodes + 1) * 2 + 1]);
    std::unique_ptr<double> zCoordinatesOut(new double[(numberOfPointsBetweenNodes + 1) * 2 + 1]);
    geometryListOut.coordinates_x = xCoordinatesOut.get();
    geometryListOut.coordinates_y = yCoordinatesOut.get();
    geometryListOut.values = zCoordinatesOut.get();

    // Execute
    auto errorCode = mkernel_get_splines(geometryListIn, geometryListOut, numberOfPointsBetweenNodes);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ((numberOfPointsBetweenNodes + 1) * 2, geometryListOut.num_coordinates);
}

TEST(ApiStatelessTests, OrthogonalizingAnInvaliMeshShouldThrowAMeshGeometryError)
{
    // Prepare
    int meshKernelId;
    int isGeographic = 0;
    meshkernelapi::mkernel_allocate_state(isGeographic, meshKernelId);

    auto mesh2d = ReadLegacyMeshFromFileForApiTesting(TEST_FOLDER + "/data/InvalidMeshes/invalid_orthogonalization_net.nc");
    auto errorCode = mkernel_set_mesh2d(meshKernelId, mesh2d);
    DeleteRectangularMeshForApiTesting(mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::OrthogonalizationParameters orthogonalizationParameters{};
    orthogonalizationParameters.outer_iterations = 1;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;

    meshkernelapi::GeometryList geometryList{};
    meshkernelapi::GeometryList landBoundaries{};

    // Execute
    errorCode = mkernel_initialize_orthogonalization_mesh2d(meshKernelId,
                                                            1,
                                                            orthogonalizationParameters,
                                                            landBoundaries,
                                                            geometryList);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert there is a geometry error
    errorCode = meshkernelapi::mkernel_prepare_outer_iteration_orthogonalization_mesh2d(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::InvalidGeometry, errorCode);

    //Delete orthogonalization instance
    errorCode = meshkernelapi::mkernel_delete_orthogonalization_mesh2d(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the message
    const char* exceptionMessage;
    errorCode = meshkernelapi::mkernel_get_error(exceptionMessage);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the index of the invalid location
    int invalidIndex;
    int type;
    errorCode = meshkernelapi::mkernel_get_geometry_error(invalidIndex, type);
    ASSERT_EQ(static_cast<int>(meshkernel::MeshLocations::Nodes), type);
    ASSERT_EQ(478, invalidIndex);
}

TEST_F(ApiTests, MakeCurvilinearGridFromPolygonThroughApi)
{
    // Prepare
    MakeMesh();
    auto meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometry_separator = meshkernel::doubleMissingValue;
    std::unique_ptr<double> xCoordinatesIn(new double[9]{
        273.502319,
        274.252319,
        275.002350,
        458.003479,
        719.005127,
        741.505249,
        710.755066,
        507.503784,
        305.002533});

    std::unique_ptr<double> yCoordinatesIn(new double[9]{
        478.880432,
        325.128906,
        172.127350,
        157.127213,
        157.127213,
        328.128937,
        490.880554,
        494.630615,
        493.130615});

    std::unique_ptr<double> zCoordinatesIn(new double[9]{
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0});

    geometryListIn.coordinates_x = xCoordinatesIn.get();
    geometryListIn.coordinates_y = yCoordinatesIn.get();
    geometryListIn.values = zCoordinatesIn.get();
    geometryListIn.num_coordinates = 9;

    // Execute
    auto errorCode = mkernel_compute_transfinite_from_polygon_curvilinear(meshKernelId,
                                                                          geometryListIn,
                                                                          0,
                                                                          2,
                                                                          4,
                                                                          true);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = meshkernelapi::mkernel_convert_curvilinear_to_mesh2d(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    errorCode = mkernel_get_dimensions_mesh2d(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(21, mesh2d.num_nodes);
    ASSERT_EQ(29, mesh2d.num_edges);
}

TEST_F(ApiTests, GetClosestMeshCoordinateThroughApi)
{
    // Prepare
    MakeMesh();
    auto meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometry_separator = meshkernel::doubleMissingValue;
    std::unique_ptr<double> xCoordinatesIn(new double[1]{-5.0});
    std::unique_ptr<double> yCoordinatesIn(new double[1]{5.0});
    std::unique_ptr<double> zCoordinatesIn(new double[1]{0.0});
    geometryListIn.coordinates_x = xCoordinatesIn.get();
    geometryListIn.coordinates_y = yCoordinatesIn.get();
    geometryListIn.values = zCoordinatesIn.get();
    geometryListIn.num_coordinates = 1;

    meshkernelapi::GeometryList geometryListOut;
    geometryListOut.geometry_separator = meshkernel::doubleMissingValue;
    std::unique_ptr<double> xCoordinatesOut(new double[1]{meshkernel::doubleMissingValue});
    std::unique_ptr<double> yCoordinatesOut(new double[1]{meshkernel::doubleMissingValue});
    std::unique_ptr<double> zCoordinatesOut(new double[1]{meshkernel::doubleMissingValue});
    geometryListOut.coordinates_x = xCoordinatesOut.get();
    geometryListOut.coordinates_y = yCoordinatesOut.get();
    geometryListOut.values = zCoordinatesOut.get();
    geometryListOut.num_coordinates = 1;

    // Execute
    auto errorCode = mkernel_get_closest_node_mesh2d(meshKernelId, geometryListIn, 10.0, geometryListOut);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(0.0, geometryListOut.coordinates_x[0]);
}

TEST_F(ApiTests, MakeCurvilinearGridFromTriangleThroughApi)
{
    // Prepare
    MakeMesh();
    auto meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometry_separator = meshkernel::doubleMissingValue;
    std::unique_ptr<double> xCoordinatesIn(new double[10]{
        444.504791,
        427.731781,
        405.640503,
        381.094666,
        451.050354,
        528.778931,
        593.416260,
        558.643005,
        526.733398,
        444.095703});
    std::unique_ptr<double> yCoordinatesIn(new double[10]{
        437.155945,
        382.745758,
        317.699005,
        262.470612,
        262.879700,
        263.288788,
        266.561584,
        324.653687,
        377.836578,
        436.746857});
    std::unique_ptr<double> zCoordinatesIn(new double[10]{
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0});
    geometryListIn.coordinates_x = xCoordinatesIn.get();
    geometryListIn.coordinates_y = yCoordinatesIn.get();
    geometryListIn.values = zCoordinatesIn.get();
    geometryListIn.num_coordinates = 10;

    // Execute
    auto errorCode = mkernel_compute_transfinite_from_triangle_curvilinear(meshKernelId,
                                                                           geometryListIn,
                                                                           0,
                                                                           3,
                                                                           6);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = meshkernelapi::mkernel_convert_curvilinear_to_mesh2d(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    errorCode = mkernel_get_dimensions_mesh2d(meshKernelId, mesh2d);

    // Assert
    ASSERT_EQ(28, mesh2d.num_nodes);
    ASSERT_EQ(40, mesh2d.num_edges);
}

TEST_F(ApiTests, MakeCurvilinearGridThroughApi)
{
    // Prepare
    auto meshKernelId = GetMeshKernelId();

    meshkernelapi::MakeMeshParameters makeMeshParameters{};
    meshkernelapi::GeometryList geometryList{};

    makeMeshParameters.grid_type = 0;
    makeMeshParameters.num_columns = 3;
    makeMeshParameters.num_rows = 2;
    makeMeshParameters.angle = 0.0;
    makeMeshParameters.block_size = 0.0;
    makeMeshParameters.origin_x = 0.0;
    makeMeshParameters.origin_y = 0.0;
    makeMeshParameters.block_size_x = 1.0;
    makeMeshParameters.block_size_y = 1.0;

    // Execute
    auto errorCode = mkernel_make_uniform_curvilinear(meshKernelId,
                                                      makeMeshParameters,
                                                      geometryList);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_get_dimensions_curvilinear(meshKernelId, curvilinearGrid);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(12, curvilinearGrid.num_nodes);
    ASSERT_EQ(17, curvilinearGrid.num_edges);

    // Allocate memory and get data
    std::unique_ptr<int> edge_nodes(new int[curvilinearGrid.num_edges * 2]);
    std::unique_ptr<double> node_x(new double[curvilinearGrid.num_nodes]);
    std::unique_ptr<double> node_y(new double[curvilinearGrid.num_nodes]);
    std::unique_ptr<double> edge_x(new double[curvilinearGrid.num_edges]);
    std::unique_ptr<double> edge_y(new double[curvilinearGrid.num_edges]);
    curvilinearGrid.edge_nodes = edge_nodes.get();
    curvilinearGrid.node_x = node_x.get();
    curvilinearGrid.node_y = node_y.get();
    curvilinearGrid.edge_x = edge_x.get();
    curvilinearGrid.edge_y = edge_y.get();
    errorCode = mkernel_get_data_curvilinear(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

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
    // Edges
    ASSERT_EQ(0, curvilinearGrid.edge_nodes[0]);
    ASSERT_EQ(4, curvilinearGrid.edge_nodes[1]);
    ASSERT_NEAR(0.0, curvilinearGrid.edge_x[0], tolerance);
    ASSERT_NEAR(0.5, curvilinearGrid.edge_y[0], tolerance);
}

TEST_F(ApiTests, GenerateTransfiniteCurvilinearGridThroughApi)
{
    // Prepare
    auto meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometry_separator = meshkernel::doubleMissingValue;
    std::unique_ptr<double> xCoordinates(new double[13]{1.340015E+02, 3.642529E+02, 6.927549E+02, meshkernel::doubleMissingValue,
                                                        2.585022E+02, 4.550035E+02, 8.337558E+02, meshkernel::doubleMissingValue,
                                                        1.002513E+02, 4.610035E+02, meshkernel::doubleMissingValue,
                                                        6.522547E+02, 7.197551E+02});

    std::unique_ptr<double> yCoordinates(new double[13]{2.546282E+02, 4.586302E+02, 5.441311E+02, meshkernel::doubleMissingValue,
                                                        6.862631E+01, 2.726284E+02, 3.753794E+02, meshkernel::doubleMissingValue,
                                                        4.068797E+02, 7.912642E+01, meshkernel::doubleMissingValue,
                                                        6.026317E+02, 2.681283E+02});

    std::unique_ptr<double> zCoordinates(new double[13]{0.0, 0.0, 0.0, meshkernel::doubleMissingValue,
                                                        0.0, 0.0, 0.0, meshkernel::doubleMissingValue,
                                                        0.0, 0.0, meshkernel::doubleMissingValue,
                                                        0.0, 0.0});

    geometryListIn.coordinates_x = xCoordinates.get();
    geometryListIn.coordinates_y = yCoordinates.get();
    geometryListIn.values = zCoordinates.get();

    geometryListIn.num_coordinates = 13;
    meshkernelapi::CurvilinearParameters curvilinearParameters;
    curvilinearParameters.m_refinement = 10;
    curvilinearParameters.n_refinement = 10;
    curvilinearParameters.smoothing_iterations = 10;
    curvilinearParameters.smoothing_parameter = 0.5;
    curvilinearParameters.attraction_parameter = 0.0;

    // Execute
    auto errorCode = mkernel_compute_transfinite_from_splines_curvilinear(meshKernelId,
                                                                          geometryListIn,
                                                                          curvilinearParameters);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_get_dimensions_curvilinear(meshKernelId, curvilinearGrid);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(121, curvilinearGrid.num_nodes);
    ASSERT_EQ(220, curvilinearGrid.num_edges);
}

TEST_F(ApiTests, GenerateOrthogonalCurvilinearGridThroughApi)
{
    // Prepare
    auto meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;

    geometryListIn.geometry_separator = meshkernel::doubleMissingValue;
    std::unique_ptr<double> xCoordinates(new double[6]{1.175014E+02, 3.755030E+02, 7.730054E+02, meshkernel::doubleMissingValue,
                                                       4.100089E+01, 3.410027E+02});

    std::unique_ptr<double> yCoordinates(new double[6]{2.437587E+01, 3.266289E+02, 4.563802E+02, meshkernel::doubleMissingValue,
                                                       2.388780E+02, 2.137584E+01});

    std::unique_ptr<double> zCoordinates(new double[6]{0.0, 0.0, 0.0, meshkernel::doubleMissingValue,
                                                       0.0, 0.0});

    geometryListIn.coordinates_x = xCoordinates.get();
    geometryListIn.coordinates_y = yCoordinates.get();
    geometryListIn.values = zCoordinates.get();
    geometryListIn.num_coordinates = 6;

    meshkernelapi::CurvilinearParameters curvilinearParameters;
    curvilinearParameters.m_refinement = 40;
    curvilinearParameters.n_refinement = 10;
    meshkernelapi::SplinesToCurvilinearParameters splinesToCurvilinearParameters;
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
    auto errorCode = mkernel_initialize_orthogonal_grid_from_splines_curvilinear(meshKernelId,
                                                                                 geometryListIn,
                                                                                 curvilinearParameters,
                                                                                 splinesToCurvilinearParameters);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Grow grid, from the second layer
    for (auto layer = 1; layer <= curvilinearParameters.n_refinement; ++layer)
    {
        errorCode = meshkernelapi::mkernel_iterate_orthogonal_grid_from_splines_curvilinear(meshKernelId, layer);
        ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    }

    // Puts the computed curvilinear mesh into the mesh state (unstructured mesh)
    errorCode = meshkernelapi::mkernel_refresh_orthogonal_grid_from_splines_curvilinear(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Delete the mesh curvilinearGridFromSplinesInstances vector entry
    errorCode = meshkernelapi::mkernel_delete_orthogonal_grid_from_splines_curvilinear(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_get_dimensions_curvilinear(meshKernelId, curvilinearGrid);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(21, curvilinearGrid.num_nodes);
    ASSERT_EQ(32, curvilinearGrid.num_edges);
}

TEST_F(ApiTests, RefineCompute_OnCurvilinearGrid_ShouldRefine)
{
    // Prepare
    auto meshKernelId = GetMeshKernelId();

    MakeUniformCurvilinearGrid(3, 3, 10);

    meshkernelapi::GeometryList firstPoint{};
    std::unique_ptr<double> xCoordinatesFirstPoint(new double[1]{10.0});
    std::unique_ptr<double> yCoordinatesFirstPoint(new double[1]{20.0});
    firstPoint.coordinates_x = xCoordinatesFirstPoint.get();
    firstPoint.coordinates_y = yCoordinatesFirstPoint.get();
    firstPoint.num_coordinates = 1;

    meshkernelapi::GeometryList secondPoint{};
    std::unique_ptr<double> xCoordinateSecondPoint(new double[1]{20.0});
    std::unique_ptr<double> yCoordinatesSecondPoint(new double[1]{20.0});
    secondPoint.coordinates_x = xCoordinateSecondPoint.get();
    secondPoint.coordinates_y = yCoordinatesSecondPoint.get();
    secondPoint.num_coordinates = 1;

    auto errorCode = mkernel_refine_curvilinear(meshKernelId, firstPoint, secondPoint, 10);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_get_dimensions_curvilinear(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    ASSERT_EQ(52, curvilinearGrid.num_nodes);
    ASSERT_EQ(87, curvilinearGrid.num_edges);
}

TEST_F(ApiTests, DerefineCompute_OnCurvilinearGrid_ShouldDeRefine)
{
    // Prepare
    auto meshKernelId = GetMeshKernelId();

    MakeUniformCurvilinearGrid();

    meshkernelapi::GeometryList firstPoint{};
    std::unique_ptr<double> xCoordinatesFirstPoint(new double[1]{10.0});
    std::unique_ptr<double> yCoordinatesFirstPoint(new double[1]{20.0});
    firstPoint.coordinates_x = xCoordinatesFirstPoint.get();
    firstPoint.coordinates_y = yCoordinatesFirstPoint.get();
    firstPoint.num_coordinates = 1;

    meshkernelapi::GeometryList secondPoint{};
    std::unique_ptr<double> xCoordinateSecondPoint(new double[1]{30.0});
    std::unique_ptr<double> yCoordinatesSecondPoint(new double[1]{20.0});
    secondPoint.coordinates_x = xCoordinateSecondPoint.get();
    secondPoint.coordinates_y = yCoordinatesSecondPoint.get();
    secondPoint.num_coordinates = 1;

    // Execute
    auto errorCode = mkernel_derefine_curvilinear(meshKernelId, firstPoint, secondPoint);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_get_dimensions_curvilinear(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    ASSERT_EQ(20, curvilinearGrid.num_nodes);
    ASSERT_EQ(31, curvilinearGrid.num_edges);
}

TEST_F(ApiTests, Orthogonalize_CurvilinearGrid_ShouldOrthogonalize)
{
    // Prepare
    auto meshKernelId = GetMeshKernelId();

    MakeUniformCurvilinearGrid();

    meshkernelapi::GeometryList firstPoint{};
    std::unique_ptr<double> xCoordinatesFirstPoint(new double[1]{10.0});
    std::unique_ptr<double> yCoordinatesFirstPoint(new double[1]{20.0});
    firstPoint.coordinates_x = xCoordinatesFirstPoint.get();
    firstPoint.coordinates_y = yCoordinatesFirstPoint.get();
    firstPoint.num_coordinates = 1;

    meshkernelapi::GeometryList secondPoint{};
    std::unique_ptr<double> xCoordinateSecondPoint(new double[1]{30.0});
    std::unique_ptr<double> yCoordinatesSecondPoint(new double[1]{20.0});
    secondPoint.coordinates_x = xCoordinateSecondPoint.get();
    secondPoint.coordinates_y = yCoordinatesSecondPoint.get();
    secondPoint.num_coordinates = 1;

    meshkernelapi::OrthogonalizationParameters orthogonalizationParameters{};
    orthogonalizationParameters.outer_iterations = 1;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;

    // Execute
    auto errorCode = mkernel_initialize_orthogonalize_curvilinear(meshKernelId, orthogonalizationParameters);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    errorCode = mkernel_set_block_orthogonalize_curvilinear(meshKernelId, firstPoint, secondPoint);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    errorCode = meshkernelapi::mkernel_orthogonalize_curvilinear(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_get_dimensions_curvilinear(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert (nothing changed)
    ASSERT_EQ(25, curvilinearGrid.num_nodes);
    ASSERT_EQ(40, curvilinearGrid.num_edges);
}

TEST_F(ApiTests, Smoothing_CurvilinearGrid_ShouldSmooth)
{
    // Prepare
    auto meshKernelId = GetMeshKernelId();

    MakeUniformCurvilinearGrid();

    meshkernelapi::GeometryList firstPoint{};
    std::unique_ptr<double> xCoordinatesFirstPoint(new double[1]{10.0});
    std::unique_ptr<double> yCoordinatesFirstPoint(new double[1]{20.0});
    firstPoint.coordinates_x = xCoordinatesFirstPoint.get();
    firstPoint.coordinates_y = yCoordinatesFirstPoint.get();
    firstPoint.num_coordinates = 1;

    meshkernelapi::GeometryList secondPoint{};
    std::unique_ptr<double> xCoordinateSecondPoint(new double[1]{30.0});
    std::unique_ptr<double> yCoordinatesSecondPoint(new double[1]{20.0});
    secondPoint.coordinates_x = xCoordinateSecondPoint.get();
    secondPoint.coordinates_y = yCoordinatesSecondPoint.get();
    secondPoint.num_coordinates = 1;

    // Execute
    auto errorCode = mkernel_smoothing_curvilinear(meshKernelId, 10, firstPoint, secondPoint);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_get_dimensions_curvilinear(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert (nothing changed)
    ASSERT_EQ(25, curvilinearGrid.num_nodes);
    ASSERT_EQ(40, curvilinearGrid.num_edges);
}

TEST_F(ApiTests, ComputedDirectionalSmooth_CurvilinearGrid_ShouldSmooth)
{
    // Prepare
    auto meshKernelId = GetMeshKernelId();

    MakeUniformCurvilinearGrid();

    meshkernelapi::GeometryList firstSegmentVertex{};
    std::unique_ptr<double> xFirstSegmentVertexPtr(new double[1]{10.0});
    std::unique_ptr<double> yFirstSegmentVertexPtr(new double[1]{0.0});
    firstSegmentVertex.coordinates_x = xFirstSegmentVertexPtr.get();
    firstSegmentVertex.coordinates_y = yFirstSegmentVertexPtr.get();
    firstSegmentVertex.num_coordinates = 1;

    meshkernelapi::GeometryList secondPointOnTheLine{};
    std::unique_ptr<double> xSecondPointOnTheLinePtr(new double[1]{10.0});
    std::unique_ptr<double> ySecondPointOnTheLinePtr(new double[1]{30.0});
    secondPointOnTheLine.coordinates_x = xSecondPointOnTheLinePtr.get();
    secondPointOnTheLine.coordinates_y = ySecondPointOnTheLinePtr.get();
    secondPointOnTheLine.num_coordinates = 1;

    meshkernelapi::GeometryList lowerLeftCornerSmoothingArea{};
    std::unique_ptr<double> xLowerLeftCornerSmoothingAreaPtr(new double[1]{10.0});
    std::unique_ptr<double> yLowerLeftCornerSmoothingAreaPtr(new double[1]{0.0});
    lowerLeftCornerSmoothingArea.coordinates_x = xLowerLeftCornerSmoothingAreaPtr.get();
    lowerLeftCornerSmoothingArea.coordinates_y = yLowerLeftCornerSmoothingAreaPtr.get();
    lowerLeftCornerSmoothingArea.num_coordinates = 1;

    meshkernelapi::GeometryList upperRightCornerSmootingArea{};
    std::unique_ptr<double> xUpperRightCornerSmootingAreaPtr(new double[1]{30.0});
    std::unique_ptr<double> yUpperRightCornerSmootingAreaPtr(new double[1]{0.0});
    upperRightCornerSmootingArea.coordinates_x = xUpperRightCornerSmootingAreaPtr.get();
    upperRightCornerSmootingArea.coordinates_y = yUpperRightCornerSmootingAreaPtr.get();
    upperRightCornerSmootingArea.num_coordinates = 1;

    // Execute
    auto errorCode = mkernel_smoothing_directional_curvilinear(meshKernelId,
                                                               10,
                                                               firstSegmentVertex,
                                                               secondPointOnTheLine,
                                                               lowerLeftCornerSmoothingArea,
                                                               upperRightCornerSmootingArea);

    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_get_dimensions_curvilinear(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert (nothing changed)
    ASSERT_EQ(25, curvilinearGrid.num_nodes);
    ASSERT_EQ(40, curvilinearGrid.num_edges);
}

TEST_F(ApiTests, ComputedLineShift_CurvilinearGrid_ShouldShift)
{
    // Prepare
    auto meshKernelId = GetMeshKernelId();

    MakeUniformCurvilinearGrid();

    meshkernelapi::mkernel_initialize_line_shift_curvilinear(meshKernelId);

    /// Sets the line to shift
    meshkernelapi::GeometryList firstGridLineNode{};
    std::unique_ptr<double> xFirstGridLineNodePtr(new double[1]{0.0});
    std::unique_ptr<double> yFirstGridLineNodePtr(new double[1]{0.0});
    firstGridLineNode.coordinates_x = xFirstGridLineNodePtr.get();
    firstGridLineNode.coordinates_y = yFirstGridLineNodePtr.get();
    firstGridLineNode.num_coordinates = 1;

    meshkernelapi::GeometryList secondGridLineNode{};
    std::unique_ptr<double> xSecondGridLineNodePtr(new double[1]{0.0});
    std::unique_ptr<double> ySecondGridLineNodePtr(new double[1]{30.0});
    secondGridLineNode.coordinates_x = xSecondGridLineNodePtr.get();
    secondGridLineNode.coordinates_y = ySecondGridLineNodePtr.get();
    secondGridLineNode.num_coordinates = 1;

    auto errorCode = mkernel_set_line_line_shift_curvilinear(meshKernelId, firstGridLineNode, secondGridLineNode);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    /// Sets the block where the shifting will be distributed
    meshkernelapi::GeometryList lowerLeftCorner{};
    std::unique_ptr<double> xLowerLeftCornerPtr(new double[1]{0.0});
    std::unique_ptr<double> yLowerLeftCornerPtr(new double[1]{0.0});
    lowerLeftCorner.coordinates_x = xLowerLeftCornerPtr.get();
    lowerLeftCorner.coordinates_y = yLowerLeftCornerPtr.get();
    lowerLeftCorner.num_coordinates = 1;

    meshkernelapi::GeometryList upperRightCorner{};
    std::unique_ptr<double> xUpperRightCorner(new double[1]{30.0});
    std::unique_ptr<double> yUpperRightCorner(new double[1]{30.0});
    upperRightCorner.coordinates_x = xUpperRightCorner.get();
    upperRightCorner.coordinates_y = yUpperRightCorner.get();
    upperRightCorner.num_coordinates = 1;

    errorCode = mkernel_set_block_line_shift_curvilinear(meshKernelId, lowerLeftCorner, upperRightCorner);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    /// Move a gridline point, in this case the origin to -10.0, 0.0
    meshkernelapi::GeometryList fromPoint{};
    std::unique_ptr<double> xfromPoint(new double[1]{0.0});
    std::unique_ptr<double> yfromPoint(new double[1]{0.0});
    fromPoint.coordinates_x = xfromPoint.get();
    fromPoint.coordinates_y = yfromPoint.get();
    fromPoint.num_coordinates = 1;

    meshkernelapi::GeometryList toPoint{};
    std::unique_ptr<double> xToPoint(new double[1]{-10.0});
    std::unique_ptr<double> yToPoint(new double[1]{0.0});
    toPoint.coordinates_x = xToPoint.get();
    toPoint.coordinates_y = yToPoint.get();
    toPoint.num_coordinates = 1;

    errorCode = mkernel_move_node_line_shift_curvilinear(meshKernelId, fromPoint, toPoint);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Execute
    errorCode = meshkernelapi::mkernel_line_shift_curvilinear(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    errorCode = meshkernelapi::mkernel_finalize_line_shift_curvilinear(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert, the nodes along the grid line have changed
    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_get_dimensions_curvilinear(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    std::unique_ptr<double> xNodesCurvilinearGrid(new double[curvilinearGrid.num_nodes]);
    std::unique_ptr<double> yNodesCurvilinearGrid(new double[curvilinearGrid.num_nodes]);
    std::unique_ptr<int> edge_nodes(new int[curvilinearGrid.num_edges * 2]);
    std::unique_ptr<double> edge_x(new double[curvilinearGrid.num_edges]);
    std::unique_ptr<double> edge_y(new double[curvilinearGrid.num_edges]);
    curvilinearGrid.node_x = xNodesCurvilinearGrid.get();
    curvilinearGrid.node_y = yNodesCurvilinearGrid.get();
    curvilinearGrid.edge_nodes = edge_nodes.get();
    curvilinearGrid.edge_x = edge_x.get();
    curvilinearGrid.edge_y = edge_y.get();

    errorCode = mkernel_get_data_curvilinear(meshKernelId, curvilinearGrid);
    ASSERT_EQ(-10.0, curvilinearGrid.node_x[0]);
    ASSERT_EQ(2.5, curvilinearGrid.node_x[1]);
    ASSERT_EQ(17.5, curvilinearGrid.node_x[2]);
    ASSERT_EQ(30.0, curvilinearGrid.node_x[3]);
}

TEST_F(ApiTests, DeleteMesh2D_WithEmptyPolygon_ShouldDeleteMesh2D)
{
    // Prepare
    MakeMesh();
    auto meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryList{};

    // Execute
    auto errorCode = mkernel_delete_mesh2d(meshKernelId, geometryList, 0, false);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_get_dimensions_mesh2d(meshKernelId, mesh2d);

    // Assert
    ASSERT_EQ(0, mesh2d.num_nodes);
    ASSERT_EQ(0, mesh2d.num_edges);
}

TEST_F(ApiTests, GetDimensionsMesh1D_WithMesh1D_ShouldGetDimensionsMesh1D)
{
    // Prepare
    MakeMesh();
    auto meshKernelId = GetMeshKernelId();

    std::vector<double> nodes_x{-16.1886410000000,
                                -16.1464995876014,
                                -16.1043581752028,
                                -16.0622167628042,
                                -15.7539488236928,
                                -6.86476658679268,
                                2.02441565010741,
                                10.9135970000000};

    std::vector<double> nodes_y{0.89018900000000,
                                9.78201442138723,
                                18.6738398427745,
                                27.5656652641617,
                                36.1966603330179,
                                36.4175095626911,
                                36.6383587923643,
                                36.8592080000000};

    std::vector<int> edges{0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7};

    meshkernelapi::Mesh1D mesh1d;
    mesh1d.edge_nodes = &edges[0];
    mesh1d.node_x = &nodes_x[0];
    mesh1d.node_y = &nodes_y[0];
    mesh1d.num_nodes = static_cast<int>(nodes_x.size());
    mesh1d.num_edges = static_cast<int>(edges.size() * 0.5);

    auto errorCode = mkernel_set_mesh1d(GetMeshKernelId(), mesh1d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Execute
    meshkernelapi::Mesh1D mesh1dResults;
    errorCode = mkernel_get_dimensions_mesh1d(meshKernelId, mesh1dResults);

    // Assert
    ASSERT_EQ(8, mesh1dResults.num_nodes);
    ASSERT_EQ(7, mesh1dResults.num_edges);
}

TEST_F(ApiTests, GetDataMesh1D_WithMesh1D_ShouldGetDataMesh1D)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    std::vector<double> nodes_x{-16.1886410000000,
                                -16.1464995876014,
                                -16.1043581752028,
                                -16.0622167628042,
                                -15.7539488236928,
                                -6.86476658679268,
                                2.02441565010741,
                                10.9135970000000};

    std::vector<double> nodes_y{0.89018900000000,
                                9.78201442138723,
                                18.6738398427745,
                                27.5656652641617,
                                36.1966603330179,
                                36.4175095626911,
                                36.6383587923643,
                                36.8592080000000};

    std::vector<int> edges{0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7};

    meshkernelapi::Mesh1D mesh1d;
    mesh1d.edge_nodes = &edges[0];
    mesh1d.node_x = &nodes_x[0];
    mesh1d.node_y = &nodes_y[0];
    mesh1d.num_nodes = static_cast<int>(nodes_x.size());
    mesh1d.num_edges = static_cast<int>(edges.size() * 0.5);

    auto errorCode = mkernel_set_mesh1d(GetMeshKernelId(), mesh1d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Execute
    meshkernelapi::Mesh1D mesh1dResults;
    errorCode = mkernel_get_dimensions_mesh1d(meshKernelId, mesh1dResults);
    std::vector<double> nodes_x_results(mesh1dResults.num_nodes);
    std::vector<double> nodes_y_results(mesh1dResults.num_nodes);
    std::vector<int> edges_results(mesh1dResults.num_edges * 2);
    mesh1dResults.node_x = &nodes_x_results[0];
    mesh1dResults.node_y = &nodes_y_results[0];
    mesh1dResults.edge_nodes = &edges_results[0];
    errorCode = mkernel_get_data_mesh1d(meshKernelId, mesh1dResults);

    // Assert
    ASSERT_THAT(nodes_x_results, ::testing::ContainerEq(nodes_x));
    ASSERT_THAT(nodes_y_results, ::testing::ContainerEq(nodes_y));
    ASSERT_THAT(edges_results, ::testing::ContainerEq(edges));
}

TEST_F(ApiTests, CountHangingEdgesMesh2D_WithZeroHangingEdges_ShouldCountZeroEdges)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    int numHangingEdges;
    auto const errorCode = meshkernelapi::mkernel_count_hanging_edges_mesh2d(meshKernelId, numHangingEdges);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    ASSERT_EQ(0, numHangingEdges);
}

TEST_F(ApiTests, GetHangingEdgesMesh2D_WithOneHangingEdges_ShouldGetOneHangingEdges)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // delete an edge at the lower left corner to create a new hanging edge
    meshkernelapi::GeometryList geometryList{};
    std::vector<double> coordinates_x(1, 0.5);
    std::vector<double> coordinates_y(1, 0.0);
    geometryList.coordinates_x = &coordinates_x[0];
    geometryList.coordinates_y = &coordinates_y[0];
    geometryList.num_coordinates = 1;
    auto errorCode = mkernel_delete_edge_mesh2d(meshKernelId, geometryList);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    int numHangingEdges;
    errorCode = meshkernelapi::mkernel_count_hanging_edges_mesh2d(meshKernelId, numHangingEdges);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_GT(numHangingEdges, 0);

    // Execute
    std::vector<int> hangingEdges(numHangingEdges);
    errorCode = meshkernelapi::mkernel_get_hanging_edges_mesh2d(meshKernelId, &hangingEdges[0]);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(hangingEdges[0], 8);
}

TEST_F(ApiTests, DeleteHangingEdgesMesh2D_WithOneHangingEdges_ShouldDeleteOneHangingEdges)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // delete an edge at the lower left corner to create a new hanging edge
    meshkernelapi::GeometryList geometryList{};
    std::vector<double> coordinates_x(1, 0.5);
    std::vector<double> coordinates_y(1, 0.0);
    geometryList.coordinates_x = &coordinates_x[0];
    geometryList.coordinates_y = &coordinates_y[0];
    geometryList.num_coordinates = 1;
    auto errorCode = mkernel_delete_edge_mesh2d(meshKernelId, geometryList);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Before deletion
    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_get_dimensions_mesh2d(meshKernelId, mesh2d);
    ASSERT_EQ(mesh2d.num_edges, 16);

    // Execute
    errorCode = meshkernelapi::mkernel_delete_hanging_edges_mesh2d(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    errorCode = mkernel_get_dimensions_mesh2d(meshKernelId, mesh2d);

    // Assert
    ASSERT_EQ(mesh2d.num_edges, 15);
}

TEST_F(ApiTests, ComputeOrthogonalizationMesh2D_WithOrthogonalMesh2D_ShouldOrthogonalize)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    meshkernelapi::OrthogonalizationParameters orthogonalizationParameters;
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

    auto errorCode = mkernel_compute_orthogonalization_mesh2d(meshKernelId, 1, orthogonalizationParameters, polygons, landBoundaries);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
}

TEST_F(ApiTests, GetOrthogonalityMesh2D_OnMesh2D_ShouldGetOrthogonality)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::Mesh2D mesh2d{};
    auto errorCode = mkernel_get_dimensions_mesh2d(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::GeometryList edgeOrthogonality;
    edgeOrthogonality.values = &std::vector<double>(mesh2d.num_edges)[0];
    edgeOrthogonality.num_coordinates = mesh2d.num_edges;

    // Execute
    errorCode = mkernel_get_orthogonality_mesh2d(meshKernelId, edgeOrthogonality);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
}

TEST_F(ApiTests, GetSmoothnessMesh2D_OnMesh2D_ShouldGetSmoothness)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::Mesh2D mesh2d{};
    auto errorCode = mkernel_get_dimensions_mesh2d(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::GeometryList edgeSmoothness;
    edgeSmoothness.values = &std::vector<double>(mesh2d.num_edges)[0];
    edgeSmoothness.num_coordinates = mesh2d.num_edges;

    // Execute
    errorCode = mkernel_get_smoothness_mesh2d(meshKernelId, edgeSmoothness);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
}

TEST_F(ApiTests, GetNodesInPolygonMesh2D_OnMesh2D_ShouldGetAllNodes)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // By using an empty list, all nodes will be selected
    const meshkernelapi::GeometryList geometryListIn{};

    meshkernelapi::Mesh2D mesh2d{};
    auto errorCode = mkernel_get_dimensions_mesh2d(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Execute
    std::vector<int> selectedNodes(mesh2d.num_nodes, -1);
    errorCode = mkernel_nodes_in_polygons_mesh2d(meshKernelId,
                                                 geometryListIn,
                                                 1,
                                                 &selectedNodes[0]);

    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert (all nodes indices will be selected)
    std::vector<int> expectedResult(mesh2d.num_nodes);
    std::iota(expectedResult.begin(), expectedResult.end(), 0);
    ASSERT_THAT(selectedNodes, ::testing::ContainerEq(expectedResult));
}

TEST_F(ApiTests, CountNodesInPolygonMesh2D_OnMesh2D_ShouldCountAllNodes)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // By using an empty list, all nodes will be selected
    const meshkernelapi::GeometryList geometryListIn{};

    // Execute
    int numNodes;
    const auto errorCode = mkernel_count_nodes_in_polygons_mesh2d(meshKernelId,
                                                                  geometryListIn,
                                                                  1,
                                                                  numNodes);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert all nodes have been selected
    ASSERT_EQ(12, numNodes);
}

TEST_F(ApiTests, InsertNodeAndEdge_OnMesh2D_ShouldInsertNodeAndEdge)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    meshkernelapi::GeometryList geometryList{};
    std::vector<double> coordinates_x(1, -0.5);
    std::vector<double> coordinates_y(1, -0.5);
    geometryList.coordinates_x = &coordinates_x[0];
    geometryList.coordinates_y = &coordinates_y[0];
    geometryList.num_coordinates = 1;

    // Isolated nodes are removed by the administration done in mkernel_get_dimensions_mesh2d.
    // The newly inserted node should be connected to another one to form an edge.
    // In this manner, the edge will not be removed during the administration
    int newNodeIndex;
    auto errorCode = mkernel_insert_node_mesh2d(meshKernelId, geometryList, newNodeIndex);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    int newEdgeIndex;
    errorCode = meshkernelapi::mkernel_insert_edge_mesh2d(meshKernelId, newNodeIndex, 0, newEdgeIndex);

    // Assert
    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_get_dimensions_mesh2d(meshKernelId, mesh2d);
    ASSERT_EQ(mesh2d.num_nodes, 13);
    ASSERT_EQ(mesh2d.num_edges, 18);
}

TEST_F(ApiTests, MoveNode_OnMesh2D_ShouldMoveNode)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    meshkernelapi::GeometryList geometryList{};
    std::vector<double> coordinates_x(1, -0.5);
    std::vector<double> coordinates_y(1, -0.5);
    geometryList.coordinates_x = &coordinates_x[0];
    geometryList.coordinates_y = &coordinates_y[0];
    geometryList.num_coordinates = 1;

    auto errorCode = mkernel_move_node_mesh2d(meshKernelId, geometryList, 0);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_get_dimensions_mesh2d(meshKernelId, mesh2d);
    std::vector<int> edge_nodes(mesh2d.num_edges * 2);
    std::vector<int> face_nodes(mesh2d.num_face_nodes);
    std::vector<int> nodes_per_face(mesh2d.num_faces);
    std::vector<double> node_x(mesh2d.num_nodes);
    std::vector<double> node_y(mesh2d.num_nodes);
    std::vector<double> edge_x(mesh2d.num_edges);
    std::vector<double> edge_y(mesh2d.num_edges);
    std::vector<double> face_x(mesh2d.num_faces);
    std::vector<double> face_y(mesh2d.num_faces);

    mesh2d.edge_nodes = &edge_nodes[0];
    mesh2d.face_nodes = &face_nodes[0];
    mesh2d.nodes_per_face = &nodes_per_face[0];
    mesh2d.node_x = &node_x[0];
    mesh2d.node_y = &node_y[0];
    mesh2d.edge_x = &edge_x[0];
    mesh2d.edge_y = &edge_y[0];
    mesh2d.face_x = &face_x[0];
    mesh2d.face_y = &face_y[0];
    errorCode = mkernel_get_data_mesh2d(meshKernelId, mesh2d);

    // Assert
    ASSERT_EQ(mesh2d.node_x[0], -0.5);
    ASSERT_EQ(mesh2d.node_y[0], -0.5);
}

TEST_F(ApiTests, GetEdge_OnMesh2D_ShouldGetAnEdgeIndex)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    meshkernelapi::GeometryList geometryList{};
    std::vector<double> coordinates_x(1, -0.5);
    std::vector<double> coordinates_y(1, -0.5);
    geometryList.coordinates_x = &coordinates_x[0];
    geometryList.coordinates_y = &coordinates_y[0];
    geometryList.num_coordinates = 1;

    int edgeIndex;
    const auto errorCode = mkernel_get_edge_mesh2d(meshKernelId, geometryList, edgeIndex);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(edgeIndex, 9);
}

TEST_F(ApiTests, GetNode_OnMesh2D_ShouldGetANodeIndex)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    meshkernelapi::GeometryList geometryList{};
    std::vector<double> coordinates_x(1, 3.0);
    std::vector<double> coordinates_y(1, 3.0);
    geometryList.coordinates_x = &coordinates_x[0];
    geometryList.coordinates_y = &coordinates_y[0];
    geometryList.num_coordinates = 1;

    int nodeIndex;
    const auto errorCode = mkernel_get_node_index_mesh2d(meshKernelId, geometryList, 10.0, nodeIndex);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(nodeIndex, 11);
}

TEST_F(ApiTests, CountSmallFlowEdges_OnMesh2D_ShouldCountSmallFlowEdges)
{
    // Prepare a mesh with two triangles
    meshkernelapi::Mesh2D mesh2d;
    std::vector<double> node_x{0.0, 1.0, 1.0, 1.0};
    std::vector<double> node_y{0.0, 0.0, 0.3, -0.3};
    std::vector<int> edge_nodes{0, 3, 3, 1, 1, 0, 1, 2, 2, 0};
    mesh2d.node_x = &node_x[0];
    mesh2d.node_y = &node_y[0];
    mesh2d.edge_nodes = &edge_nodes[0];
    mesh2d.num_edges = static_cast<int>(edge_nodes.size() * 0.5);
    mesh2d.num_nodes = node_x.size();

    const auto meshKernelId = GetMeshKernelId();

    // Execute
    auto errorCode = mkernel_set_mesh2d(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    int numSmallFlowEdges;
    errorCode = meshkernelapi::mkernel_count_small_flow_edge_centers_mesh2d(meshKernelId, 100, numSmallFlowEdges);
    ASSERT_EQ(1, numSmallFlowEdges);
}