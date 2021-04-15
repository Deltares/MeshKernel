#include <exception>
#include <memory>

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
        auto mesh2d = MakeRectangularMeshForApiTesting(n, m, delta);
        auto errorCode = mkernel_set_mesh2d(m_meshKernelId, mesh2d);
        if (errorCode != 0)
        {
            throw std::runtime_error("Could not set mesh2d");
        }
        // Clean up client-allocated memory
        DeleteRectangularMeshForApiTesting(mesh2d);
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
    errorCode = meshkernelapi::mkernel_get_dimensions_mesh2d(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert dimensions
    ASSERT_EQ(11, mesh2d.num_nodes);
    ASSERT_EQ(15, mesh2d.num_edges);

    // Allocate memory and get data
    std::unique_ptr<int> edge_nodes(new int[mesh2d.num_edges * 2]);
    std::unique_ptr<int> face_nodes(new int[mesh2d.num_face_nodes]);
    std::unique_ptr<int> nodes_per_face(new int[mesh2d.num_faces]);
    std::unique_ptr<double> node_x(new double[mesh2d.num_nodes]);
    std::unique_ptr<double> node_y(new double[mesh2d.num_nodes]);
    std::unique_ptr<double> edge_x(new double[mesh2d.num_edges]);
    std::unique_ptr<double> edge_y(new double[mesh2d.num_edges]);
    std::unique_ptr<double> face_x(new double[mesh2d.num_faces]);
    std::unique_ptr<double> face_y(new double[mesh2d.num_faces]);
    mesh2d.edge_nodes = edge_nodes.get();
    mesh2d.face_nodes = face_nodes.get();
    mesh2d.nodes_per_face = nodes_per_face.get();
    mesh2d.node_x = node_x.get();
    mesh2d.node_y = node_y.get();
    mesh2d.edge_x = edge_x.get();
    mesh2d.edge_y = edge_y.get();
    mesh2d.face_x = face_x.get();
    mesh2d.face_y = face_y.get();
    errorCode = meshkernelapi::mkernel_get_data_mesh2d(meshKernelId, mesh2d);
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
    errorCode = meshkernelapi::mkernel_get_dimensions_mesh2d(meshKernelId, mesh2d);

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
    errorCode = meshkernelapi::mkernel_get_dimensions_mesh2d(meshKernelId, mesh2d);

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
    errorCode = meshkernelapi::mkernel_get_dimensions_mesh2d(meshKernelId, mesh2d);

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
    errorCode = meshkernelapi::mkernel_get_dimensions_mesh2d(meshKernelId, mesh2d);

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
    orthogonalizationParameters.OuterIterations = 1;
    orthogonalizationParameters.BoundaryIterations = 25;
    orthogonalizationParameters.InnerIterations = 25;
    orthogonalizationParameters.OrthogonalizationToSmoothingFactor = 0.975;

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
    geometryListIn.geometrySeparator = meshkernel::doubleMissingValue;
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

    geometryListIn.xCoordinates = xCoordinates.get();
    geometryListIn.yCoordinates = yCoordinates.get();
    geometryListIn.zCoordinates = zCoordinates.get();
    geometryListIn.numberOfCoordinates = 17;

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

    geometryListIn.geometrySeparator = meshkernel::doubleMissingValue;

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

    geometryListIn.xCoordinates = xCoordinates.get();
    geometryListIn.yCoordinates = yCoordinates.get();
    geometryListIn.zCoordinates = zCoordinates.get();

    geometryListIn.numberOfCoordinates = 5;

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
    geometryListOut.geometrySeparator = meshkernel::doubleMissingValue;
    geometryListOut.numberOfCoordinates = numberOfpolygonNodes;

    std::unique_ptr<double> xCoordinates(new double[numberOfpolygonNodes]);
    std::unique_ptr<double> yCoordinates(new double[numberOfpolygonNodes]);
    std::unique_ptr<double> zCoordinates(new double[numberOfpolygonNodes]);

    geometryListOut.xCoordinates = xCoordinates.get();
    geometryListOut.yCoordinates = yCoordinates.get();
    geometryListOut.zCoordinates = zCoordinates.get();

    // Execute
    errorCode = mkernel_get_mesh_boundaries_to_polygon_mesh2d(meshKernelId, geometryListOut);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    const double tolerance = 1e-6;
    ASSERT_NEAR(0.0, geometryListOut.xCoordinates[0], tolerance);
    ASSERT_NEAR(0.0, geometryListOut.yCoordinates[0], tolerance);
}

TEST_F(ApiTests, OffsetAPolygonThroughApi)
{
    // Prepare
    MakeMesh();
    auto meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometrySeparator = meshkernel::doubleMissingValue;
    geometryListIn.numberOfCoordinates = 4;

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

    geometryListIn.xCoordinates = xCoordinatesIn.get();
    geometryListIn.yCoordinates = yCoordinatesIn.get();
    geometryListIn.zCoordinates = zCoordinatesIn.get();

    // Execute
    int numberOfpolygonNodes;
    auto errorCode = mkernel_count_offset_polygon(meshKernelId, geometryListIn, false, 0.5, numberOfpolygonNodes);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(4, numberOfpolygonNodes);

    meshkernelapi::GeometryList geometryListOut;

    geometryListOut.numberOfCoordinates = numberOfpolygonNodes;
    geometryListOut.geometrySeparator = meshkernel::doubleMissingValue;

    std::unique_ptr<double> xCoordinatesOut(new double[numberOfpolygonNodes]);
    std::unique_ptr<double> yCoordinatesOut(new double[numberOfpolygonNodes]);
    std::unique_ptr<double> zCoordinatesOut(new double[numberOfpolygonNodes]);
    geometryListOut.xCoordinates = xCoordinatesOut.get();
    geometryListOut.yCoordinates = yCoordinatesOut.get();
    geometryListOut.zCoordinates = zCoordinatesOut.get();
    errorCode = mkernel_get_offset_polygon(meshKernelId, geometryListIn, false, 10.0, geometryListOut);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    const double tolerance = 1e-6;
    ASSERT_NEAR(0.0, geometryListOut.xCoordinates[0], tolerance);
    ASSERT_NEAR(-10.0, geometryListOut.yCoordinates[0], tolerance);
}

TEST_F(ApiTests, RefineAPolygonThroughApi)
{
    // Prepare
    MakeMesh();
    auto meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometrySeparator = meshkernel::doubleMissingValue;
    geometryListIn.numberOfCoordinates = 3;
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

    geometryListIn.xCoordinates = xCoordinatesIn.get();
    geometryListIn.yCoordinates = yCoordinatesIn.get();
    geometryListIn.zCoordinates = zCoordinatesIn.get();

    // Execute
    int numberOfpolygonNodes;
    auto errorCode = mkernel_count_refine_polygon(meshKernelId, geometryListIn, 0, 2, 40, numberOfpolygonNodes);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(22, numberOfpolygonNodes);

    meshkernelapi::GeometryList geometryListOut;
    geometryListOut.numberOfCoordinates = numberOfpolygonNodes;
    geometryListOut.geometrySeparator = meshkernel::doubleMissingValue;
    std::unique_ptr<double> xCoordinatesOut(new double[numberOfpolygonNodes]);
    std::unique_ptr<double> yCoordinatesOut(new double[numberOfpolygonNodes]);
    std::unique_ptr<double> zCoordinatesOut(new double[numberOfpolygonNodes]);
    geometryListOut.xCoordinates = xCoordinatesOut.get();
    geometryListOut.yCoordinates = yCoordinatesOut.get();
    geometryListOut.zCoordinates = zCoordinatesOut.get();
    errorCode = mkernel_refine_polygon(meshKernelId, geometryListIn, false, 0, 2, geometryListOut);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    const double tolerance = 1e-6;
    ASSERT_NEAR(76.251099, geometryListOut.xCoordinates[0], tolerance);
    ASSERT_NEAR(92.626556, geometryListOut.yCoordinates[0], tolerance);
}

TEST_F(ApiTests, RefineAGridBasedOnSamplesThroughApi)
{
    // Prepare
    MakeMesh();
    auto meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometrySeparator = meshkernel::doubleMissingValue;
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

    geometryListIn.xCoordinates = xCoordinatesIn.get();
    geometryListIn.yCoordinates = yCoordinatesIn.get();
    geometryListIn.zCoordinates = zCoordinatesIn.get();

    geometryListIn.numberOfCoordinates = 9;

    meshkernelapi::InterpolationParameters interpolationParameters;
    interpolationParameters.InterpolationType = 1;
    interpolationParameters.DisplayInterpolationProcess = 0;
    interpolationParameters.MaxNumberOfRefinementIterations = 2;
    interpolationParameters.AveragingMethod = 1;
    interpolationParameters.MinimumNumberOfPoints = 1;
    interpolationParameters.RelativeSearchRadius = 1.01;
    interpolationParameters.InterpolateTo = 3;
    interpolationParameters.RefineIntersected = 0;

    meshkernelapi::SampleRefineParameters samplesRefineParameters;
    samplesRefineParameters.SampleVectorDimension = 1;
    samplesRefineParameters.MinimumCellSize = 0.5;
    samplesRefineParameters.DirectionalRefinement = 0;
    samplesRefineParameters.RefinementType = 3;
    samplesRefineParameters.ConnectHangingNodes = 1;
    samplesRefineParameters.MaximumTimeStepInCourantGrid = 0.0;
    samplesRefineParameters.AccountForSamplesOutside = 0;

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
    geometryListIn.geometrySeparator = meshkernel::doubleMissingValue;
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

    geometryListIn.xCoordinates = xCoordinatesIn.get();
    geometryListIn.yCoordinates = yCoordinatesIn.get();
    geometryListIn.zCoordinates = zCoordinatesIn.get();

    geometryListIn.numberOfCoordinates = 9;

    meshkernelapi::InterpolationParameters interpolationParameters;
    interpolationParameters.InterpolationType = 1;
    interpolationParameters.DisplayInterpolationProcess = 0;
    interpolationParameters.MaxNumberOfRefinementIterations = 2;
    interpolationParameters.AveragingMethod = 1;
    interpolationParameters.MinimumNumberOfPoints = 1;
    interpolationParameters.RelativeSearchRadius = 1.01;
    interpolationParameters.InterpolateTo = 3;
    interpolationParameters.RefineIntersected = 0;

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
    polygon.geometrySeparator = meshkernel::doubleMissingValue;

    std::unique_ptr<double> xCoordinates(new double[5]{-30, 40, 40, -40, -30});
    std::unique_ptr<double> yCoordinates(new double[5]{-20, -20, 50, 50, -20});
    std::unique_ptr<double> zCoordinates(new double[5]{0, 0, 0, 0, 0});
    polygon.xCoordinates = xCoordinates.get();
    polygon.yCoordinates = yCoordinates.get();
    polygon.zCoordinates = zCoordinates.get();

    polygon.numberOfCoordinates = 5;

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
    polygon.geometrySeparator = meshkernel::doubleMissingValue;

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
    polygon.xCoordinates = xCoordinates.get();
    polygon.yCoordinates = yCoordinates.get();
    polygon.zCoordinates = zCoordinates.get();

    polygon.numberOfCoordinates = 5;

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
    points.geometrySeparator = meshkernel::doubleMissingValue;

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
    points.xCoordinates = xCoordinates.get();
    points.yCoordinates = yCoordinates.get();
    points.zCoordinates = zCoordinates.get();

    points.numberOfCoordinates = 4;

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
    polygon.geometrySeparator = meshkernel::doubleMissingValue;

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
    polygon.xCoordinates = xCoordinates.get();
    polygon.yCoordinates = yCoordinates.get();
    polygon.zCoordinates = zCoordinates.get();

    polygon.numberOfCoordinates = 5;

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
    geometryListIn.xCoordinates = xCoordinatesIn.get();
    geometryListIn.yCoordinates = yCoordinatesIn.get();
    geometryListIn.zCoordinates = zCoordinatesIn.get();
    geometryListIn.numberOfCoordinates = 3;
    geometryListIn.geometrySeparator = meshkernel::doubleMissingValue;

    meshkernelapi::GeometryList geometryListOut;
    int numberOfPointsBetweenNodes = 20;
    std::unique_ptr<double> xCoordinatesOut(new double[(numberOfPointsBetweenNodes + 1) * 2 + 1]);
    std::unique_ptr<double> yCoordinatesOut(new double[(numberOfPointsBetweenNodes + 1) * 2 + 1]);
    std::unique_ptr<double> zCoordinatesOut(new double[(numberOfPointsBetweenNodes + 1) * 2 + 1]);
    geometryListOut.xCoordinates = xCoordinatesOut.get();
    geometryListOut.yCoordinates = yCoordinatesOut.get();
    geometryListOut.zCoordinates = zCoordinatesOut.get();

    // Execute
    auto errorCode = mkernel_get_splines(geometryListIn, geometryListOut, numberOfPointsBetweenNodes);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ((numberOfPointsBetweenNodes + 1) * 2, geometryListOut.numberOfCoordinates);
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
    orthogonalizationParameters.OuterIterations = 1;
    orthogonalizationParameters.BoundaryIterations = 25;
    orthogonalizationParameters.InnerIterations = 25;
    orthogonalizationParameters.OrthogonalizationToSmoothingFactor = 0.975;

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
    geometryListIn.geometrySeparator = meshkernel::doubleMissingValue;
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

    geometryListIn.xCoordinates = xCoordinatesIn.get();
    geometryListIn.yCoordinates = yCoordinatesIn.get();
    geometryListIn.zCoordinates = zCoordinatesIn.get();
    geometryListIn.numberOfCoordinates = 9;

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
    geometryListIn.geometrySeparator = meshkernel::doubleMissingValue;
    std::unique_ptr<double> xCoordinatesIn(new double[1]{-5.0});
    std::unique_ptr<double> yCoordinatesIn(new double[1]{5.0});
    std::unique_ptr<double> zCoordinatesIn(new double[1]{0.0});
    geometryListIn.xCoordinates = xCoordinatesIn.get();
    geometryListIn.yCoordinates = yCoordinatesIn.get();
    geometryListIn.zCoordinates = zCoordinatesIn.get();
    geometryListIn.numberOfCoordinates = 1;

    meshkernelapi::GeometryList geometryListOut;
    geometryListOut.geometrySeparator = meshkernel::doubleMissingValue;
    std::unique_ptr<double> xCoordinatesOut(new double[1]{meshkernel::doubleMissingValue});
    std::unique_ptr<double> yCoordinatesOut(new double[1]{meshkernel::doubleMissingValue});
    std::unique_ptr<double> zCoordinatesOut(new double[1]{meshkernel::doubleMissingValue});
    geometryListOut.xCoordinates = xCoordinatesOut.get();
    geometryListOut.yCoordinates = yCoordinatesOut.get();
    geometryListOut.zCoordinates = zCoordinatesOut.get();
    geometryListOut.numberOfCoordinates = 1;

    // Execute
    auto errorCode = mkernel_get_closest_node_mesh2d(meshKernelId, geometryListIn, 10.0, geometryListOut);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(0.0, geometryListOut.xCoordinates[0]);
}

TEST_F(ApiTests, MakeCurvilinearGridFromTriangleThroughApi)
{
    // Prepare
    MakeMesh();
    auto meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometrySeparator = meshkernel::doubleMissingValue;
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
    geometryListIn.xCoordinates = xCoordinatesIn.get();
    geometryListIn.yCoordinates = yCoordinatesIn.get();
    geometryListIn.zCoordinates = zCoordinatesIn.get();
    geometryListIn.numberOfCoordinates = 10;

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

    makeMeshParameters.GridType = 0;
    makeMeshParameters.NumberOfColumns = 3;
    makeMeshParameters.NumberOfRows = 2;
    makeMeshParameters.GridAngle = 0.0;
    makeMeshParameters.GridBlockSize = 0.0;
    makeMeshParameters.OriginXCoordinate = 0.0;
    makeMeshParameters.OriginYCoordinate = 0.0;
    makeMeshParameters.OriginZCoordinate = 0.0;
    makeMeshParameters.XGridBlockSize = 1.0;
    makeMeshParameters.YGridBlockSize = 1.0;

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
    errorCode = meshkernelapi::mkernel_get_data_curvilinear(meshKernelId, curvilinearGrid);
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
    geometryListIn.geometrySeparator = meshkernel::doubleMissingValue;
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

    geometryListIn.xCoordinates = xCoordinates.get();
    geometryListIn.yCoordinates = yCoordinates.get();
    geometryListIn.zCoordinates = zCoordinates.get();

    geometryListIn.numberOfCoordinates = 13;
    meshkernelapi::CurvilinearParameters curvilinearParameters;
    curvilinearParameters.MRefinement = 10;
    curvilinearParameters.NRefinement = 10;
    curvilinearParameters.SmoothingIterations = 10;
    curvilinearParameters.SmoothingParameter = 0.5;
    curvilinearParameters.AttractionParameter = 0.0;

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

    geometryListIn.geometrySeparator = meshkernel::doubleMissingValue;
    std::unique_ptr<double> xCoordinates(new double[6]{1.175014E+02, 3.755030E+02, 7.730054E+02, meshkernel::doubleMissingValue,
                                                       4.100089E+01, 3.410027E+02});

    std::unique_ptr<double> yCoordinates(new double[6]{2.437587E+01, 3.266289E+02, 4.563802E+02, meshkernel::doubleMissingValue,
                                                       2.388780E+02, 2.137584E+01});

    std::unique_ptr<double> zCoordinates(new double[6]{0.0, 0.0, 0.0, meshkernel::doubleMissingValue,
                                                       0.0, 0.0});

    geometryListIn.xCoordinates = xCoordinates.get();
    geometryListIn.yCoordinates = yCoordinates.get();
    geometryListIn.zCoordinates = zCoordinates.get();
    geometryListIn.numberOfCoordinates = 6;

    meshkernelapi::CurvilinearParameters curvilinearParameters;
    curvilinearParameters.MRefinement = 40;
    curvilinearParameters.NRefinement = 10;
    meshkernelapi::SplinesToCurvilinearParameters splinesToCurvilinearParameters;
    splinesToCurvilinearParameters.AspectRatio = 0.1;
    splinesToCurvilinearParameters.AspectRatioGrowFactor = 1.1;
    splinesToCurvilinearParameters.AverageWidth = 500.0;
    splinesToCurvilinearParameters.CurvatureAdaptedGridSpacing = 1;
    splinesToCurvilinearParameters.GrowGridOutside = 1;
    splinesToCurvilinearParameters.MaximumNumberOfGridCellsInTheUniformPart = 5;
    splinesToCurvilinearParameters.GridsOnTopOfEachOtherTolerance = 0.0001;
    splinesToCurvilinearParameters.MinimumCosineOfCrossingAngles = 0.95;
    splinesToCurvilinearParameters.CheckFrontCollisions = 0;
    splinesToCurvilinearParameters.UniformGridSize = 0.0;
    splinesToCurvilinearParameters.DeleteSkinnyTriangles = 1;
    splinesToCurvilinearParameters.GrowGridOutside = 0;

    // Execute
    auto errorCode = mkernel_initialize_orthogonal_grid_from_splines_curvilinear(meshKernelId,
                                                                                 geometryListIn,
                                                                                 curvilinearParameters,
                                                                                 splinesToCurvilinearParameters);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Grow grid, from the second layer
    for (auto layer = 1; layer <= curvilinearParameters.NRefinement; ++layer)
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

    meshkernelapi::MakeMeshParameters makeMeshParameters{};
    meshkernelapi::GeometryList geometryList{};

    makeMeshParameters.GridType = 0;
    makeMeshParameters.NumberOfColumns = 3;
    makeMeshParameters.NumberOfRows = 3;
    makeMeshParameters.GridAngle = 0.0;
    makeMeshParameters.GridBlockSize = 0.0;
    makeMeshParameters.OriginXCoordinate = 0.0;
    makeMeshParameters.OriginYCoordinate = 0.0;
    makeMeshParameters.OriginZCoordinate = 0.0;
    makeMeshParameters.XGridBlockSize = 10.0;
    makeMeshParameters.YGridBlockSize = 10.0;

    // Execute
    auto errorCode = mkernel_make_uniform_curvilinear(meshKernelId,
                                                      makeMeshParameters,
                                                      geometryList);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::GeometryList firstPoint{};
    std::unique_ptr<double> xCoordinatesFirstPoint(new double[1]{10.0});
    std::unique_ptr<double> yCoordinatesFirstPoint(new double[1]{20.0});
    firstPoint.xCoordinates = xCoordinatesFirstPoint.get();
    firstPoint.yCoordinates = yCoordinatesFirstPoint.get();
    firstPoint.numberOfCoordinates = 1;

    meshkernelapi::GeometryList secondPoint{};
    std::unique_ptr<double> xCoordinateSecondPoint(new double[1]{20.0});
    std::unique_ptr<double> yCoordinatesSecondPoint(new double[1]{20.0});
    secondPoint.xCoordinates = xCoordinateSecondPoint.get();
    secondPoint.yCoordinates = yCoordinatesSecondPoint.get();
    secondPoint.numberOfCoordinates = 1;

    errorCode = mkernel_refine_curvilinear(meshKernelId, firstPoint, secondPoint, 10);
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

    meshkernelapi::MakeMeshParameters makeMeshParameters{};
    meshkernelapi::GeometryList geometryList{};

    makeMeshParameters.GridType = 0;
    makeMeshParameters.NumberOfColumns = 4;
    makeMeshParameters.NumberOfRows = 4;
    makeMeshParameters.GridAngle = 0.0;
    makeMeshParameters.GridBlockSize = 0.0;
    makeMeshParameters.OriginXCoordinate = 0.0;
    makeMeshParameters.OriginYCoordinate = 0.0;
    makeMeshParameters.OriginZCoordinate = 0.0;
    makeMeshParameters.XGridBlockSize = 10.0;
    makeMeshParameters.YGridBlockSize = 10.0;

    auto errorCode = mkernel_make_uniform_curvilinear(meshKernelId,
                                                      makeMeshParameters,
                                                      geometryList);

    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::GeometryList firstPoint{};
    std::unique_ptr<double> xCoordinatesFirstPoint(new double[1]{10.0});
    std::unique_ptr<double> yCoordinatesFirstPoint(new double[1]{20.0});
    firstPoint.xCoordinates = xCoordinatesFirstPoint.get();
    firstPoint.yCoordinates = yCoordinatesFirstPoint.get();
    firstPoint.numberOfCoordinates = 1;

    meshkernelapi::GeometryList secondPoint{};
    std::unique_ptr<double> xCoordinateSecondPoint(new double[1]{30.0});
    std::unique_ptr<double> yCoordinatesSecondPoint(new double[1]{20.0});
    secondPoint.xCoordinates = xCoordinateSecondPoint.get();
    secondPoint.yCoordinates = yCoordinatesSecondPoint.get();
    secondPoint.numberOfCoordinates = 1;

    // Execute
    errorCode = mkernel_derefine_curvilinear(meshKernelId, firstPoint, secondPoint);
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

    meshkernelapi::MakeMeshParameters makeMeshParameters{};
    meshkernelapi::GeometryList geometryList{};

    makeMeshParameters.GridType = 0;
    makeMeshParameters.NumberOfColumns = 4;
    makeMeshParameters.NumberOfRows = 4;
    makeMeshParameters.GridAngle = 0.0;
    makeMeshParameters.GridBlockSize = 0.0;
    makeMeshParameters.OriginXCoordinate = 0.0;
    makeMeshParameters.OriginYCoordinate = 0.0;
    makeMeshParameters.OriginZCoordinate = 0.0;
    makeMeshParameters.XGridBlockSize = 10.0;
    makeMeshParameters.YGridBlockSize = 10.0;

    auto errorCode = mkernel_make_uniform_curvilinear(meshKernelId,
                                                      makeMeshParameters,
                                                      geometryList);

    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::GeometryList firstPoint{};
    std::unique_ptr<double> xCoordinatesFirstPoint(new double[1]{10.0});
    std::unique_ptr<double> yCoordinatesFirstPoint(new double[1]{20.0});
    firstPoint.xCoordinates = xCoordinatesFirstPoint.get();
    firstPoint.yCoordinates = yCoordinatesFirstPoint.get();
    firstPoint.numberOfCoordinates = 1;

    meshkernelapi::GeometryList secondPoint{};
    std::unique_ptr<double> xCoordinateSecondPoint(new double[1]{30.0});
    std::unique_ptr<double> yCoordinatesSecondPoint(new double[1]{20.0});
    secondPoint.xCoordinates = xCoordinateSecondPoint.get();
    secondPoint.yCoordinates = yCoordinatesSecondPoint.get();
    secondPoint.numberOfCoordinates = 1;

    meshkernelapi::OrthogonalizationParameters orthogonalizationParameters{};
    orthogonalizationParameters.OuterIterations = 1;
    orthogonalizationParameters.BoundaryIterations = 25;
    orthogonalizationParameters.InnerIterations = 25;
    orthogonalizationParameters.OrthogonalizationToSmoothingFactor = 0.975;

    // Execute
    errorCode = mkernel_initialize_orthogonalize_curvilinear(meshKernelId, orthogonalizationParameters);
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

    meshkernelapi::MakeMeshParameters makeMeshParameters{};
    meshkernelapi::GeometryList geometryList{};

    makeMeshParameters.GridType = 0;
    makeMeshParameters.NumberOfColumns = 4;
    makeMeshParameters.NumberOfRows = 4;
    makeMeshParameters.GridAngle = 0.0;
    makeMeshParameters.GridBlockSize = 0.0;
    makeMeshParameters.OriginXCoordinate = 0.0;
    makeMeshParameters.OriginYCoordinate = 0.0;
    makeMeshParameters.OriginZCoordinate = 0.0;
    makeMeshParameters.XGridBlockSize = 10.0;
    makeMeshParameters.YGridBlockSize = 10.0;

    auto errorCode = mkernel_make_uniform_curvilinear(meshKernelId,
                                                      makeMeshParameters,
                                                      geometryList);

    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::GeometryList firstPoint{};
    std::unique_ptr<double> xCoordinatesFirstPoint(new double[1]{10.0});
    std::unique_ptr<double> yCoordinatesFirstPoint(new double[1]{20.0});
    firstPoint.xCoordinates = xCoordinatesFirstPoint.get();
    firstPoint.yCoordinates = yCoordinatesFirstPoint.get();
    firstPoint.numberOfCoordinates = 1;

    meshkernelapi::GeometryList secondPoint{};
    std::unique_ptr<double> xCoordinateSecondPoint(new double[1]{30.0});
    std::unique_ptr<double> yCoordinatesSecondPoint(new double[1]{20.0});
    secondPoint.xCoordinates = xCoordinateSecondPoint.get();
    secondPoint.yCoordinates = yCoordinatesSecondPoint.get();
    secondPoint.numberOfCoordinates = 1;

    // Execute
    errorCode = mkernel_smoothing_curvilinear(meshKernelId, 10, firstPoint, secondPoint);
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

    meshkernelapi::MakeMeshParameters makeMeshParameters{};
    meshkernelapi::GeometryList geometryList{};

    makeMeshParameters.GridType = 0;
    makeMeshParameters.NumberOfColumns = 4;
    makeMeshParameters.NumberOfRows = 4;
    makeMeshParameters.GridAngle = 0.0;
    makeMeshParameters.GridBlockSize = 0.0;
    makeMeshParameters.OriginXCoordinate = 0.0;
    makeMeshParameters.OriginYCoordinate = 0.0;
    makeMeshParameters.OriginZCoordinate = 0.0;
    makeMeshParameters.XGridBlockSize = 10.0;
    makeMeshParameters.YGridBlockSize = 10.0;

    auto errorCode = mkernel_make_uniform_curvilinear(meshKernelId,
                                                      makeMeshParameters,
                                                      geometryList);

    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::GeometryList firstSegmentVertex{};
    std::unique_ptr<double> xCoordinatesFirstSegmentVertex(new double[1]{10.0});
    std::unique_ptr<double> yCoordinatesFirstSegmentVertex(new double[1]{0.0});
    firstSegmentVertex.xCoordinates = xCoordinatesFirstSegmentVertex.get();
    firstSegmentVertex.yCoordinates = yCoordinatesFirstSegmentVertex.get();
    firstSegmentVertex.numberOfCoordinates = 1;

    meshkernelapi::GeometryList secondPointOnTheLine{};
    std::unique_ptr<double> xCoordinateSecondSegmentVertex(new double[1]{10.0});
    std::unique_ptr<double> yCoordinatesSecondSegmentVertex(new double[1]{30.0});
    secondPointOnTheLine.xCoordinates = xCoordinateSecondSegmentVertex.get();
    secondPointOnTheLine.yCoordinates = yCoordinatesSecondSegmentVertex.get();
    secondPointOnTheLine.numberOfCoordinates = 1;

    meshkernelapi::GeometryList lowerLeftCornerSmoothingArea{};
    std::unique_ptr<double> xCoordinatesLowerLeftCornerSmoothingArea(new double[1]{10.0});
    std::unique_ptr<double> yCoordinatesLowerLeftCornerSmoothingArea(new double[1]{0.0});
    lowerLeftCornerSmoothingArea.xCoordinates = xCoordinatesLowerLeftCornerSmoothingArea.get();
    lowerLeftCornerSmoothingArea.yCoordinates = yCoordinatesLowerLeftCornerSmoothingArea.get();
    lowerLeftCornerSmoothingArea.numberOfCoordinates = 1;

    meshkernelapi::GeometryList upperRightCornerSmootingArea{};
    std::unique_ptr<double> xCoordinateSecondUpperRightCornerSmoothingArea(new double[1]{30.0});
    std::unique_ptr<double> yCoordinatesSecondUpperRightCornerSmoothingArea(new double[1]{0.0});
    upperRightCornerSmootingArea.xCoordinates = xCoordinateSecondUpperRightCornerSmoothingArea.get();
    upperRightCornerSmootingArea.yCoordinates = yCoordinatesSecondUpperRightCornerSmoothingArea.get();
    upperRightCornerSmootingArea.numberOfCoordinates = 1;

    // Execute
    errorCode = mkernel_smoothing_directional_curvilinear(meshKernelId,
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
