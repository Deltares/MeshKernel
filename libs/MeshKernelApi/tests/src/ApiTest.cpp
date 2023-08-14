#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <MeshKernel/Parameters.hpp>

#include <MeshKernelApi/CurvilinearGrid.hpp>
#include <MeshKernelApi/GeometryList.hpp>
#include <MeshKernelApi/Mesh1D.hpp>
#include <MeshKernelApi/Mesh2D.hpp>
#include <MeshKernelApi/MeshKernel.hpp>
#include <Version/Version.hpp>

#include <TestUtils/Definitions.hpp>
#include <TestUtils/MakeCurvilinearGrids.hpp>
#include <TestUtils/MakeMeshes.hpp>
#include <TestUtils/SampleFileReader.hpp>

#include <numeric>

class CartesianApiTests : public testing::Test
{
public:
    /// Constructor for allocating state
    CartesianApiTests()
    {
        int isGeographic = 0;
        const auto errorCode = meshkernelapi::mkernel_allocate_state(isGeographic, m_meshKernelId);
        if (errorCode != 0)
        {
            throw std::runtime_error("Could not allocate state");
        }
    }

    /// Destructor for deallocating state
    ~CartesianApiTests()
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

    void MakeUniformCurvilinearGrid(meshkernel::UInt numberOfColumns = 4,
                                    meshkernel::UInt numberOfRows = 4,
                                    double blockSizeX = 10.0,
                                    double blockSizeY = 10.0,
                                    double originX = 0.0,
                                    double originY = 0.0) const
    {
        meshkernel::MakeGridParameters makeGridParameters{};
        meshkernelapi::GeometryList geometryList{};

        makeGridParameters.num_columns = static_cast<int>(numberOfColumns);
        makeGridParameters.num_rows = static_cast<int>(numberOfRows);
        makeGridParameters.angle = 0.0;
        makeGridParameters.origin_x = originX;
        makeGridParameters.origin_y = originY;
        makeGridParameters.block_size_x = blockSizeX;
        makeGridParameters.block_size_y = blockSizeY;

        auto const errorCode = mkernel_curvilinear_make_uniform(m_meshKernelId, makeGridParameters, geometryList);
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

TEST_F(CartesianApiTests, Mesh2DDeleteNode_ShouldDeleteNode)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    auto errorCode = meshkernelapi::mkernel_mesh2d_delete_node(meshKernelId, 0);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the dimensions
    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
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

    mesh2d.edge_nodes = edge_nodes.data();
    mesh2d.face_nodes = face_nodes.data();
    mesh2d.nodes_per_face = nodes_per_face.data();
    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_x = edge_x.data();
    mesh2d.edge_y = edge_y.data();
    mesh2d.face_x = face_x.data();
    mesh2d.face_y = face_y.data();
    errorCode = mkernel_mesh2d_get_data(meshKernelId, mesh2d);
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

TEST_F(CartesianApiTests, FlipEdges_ShouldFlipEdges)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    const int isTriangulationRequired = 1;
    const int projectToLandBoundaryOption = 1;
    meshkernelapi::GeometryList selectingPolygon{};
    meshkernelapi::GeometryList landBoundaries{};
    auto errorCode = mkernel_mesh2d_flip_edges(meshKernelId,
                                               isTriangulationRequired,
                                               projectToLandBoundaryOption,
                                               selectingPolygon,
                                               landBoundaries);

    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(12, mesh2d.num_nodes);
    ASSERT_EQ(23, mesh2d.num_edges);
}

TEST_F(CartesianApiTests, FlipEdges_WithALandBoundary_ShouldFlipEdges)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    const int isTriangulationRequired = 1;
    const int projectToLandBoundaryOption = 1;
    meshkernelapi::GeometryList selectingPolygon{};

    std::vector xCoordinates{-0.5, -0.5, 4.0, meshkernel::constants::missing::doubleValue};
    std::vector yCoordinates{3.0, -0.5, -0.5, meshkernel::constants::missing::doubleValue};
    std::vector zCoordinates{0.0, 0.0, 0.0, meshkernel::constants::missing::doubleValue};

    meshkernelapi::GeometryList landBoundaries{};
    landBoundaries.geometry_separator = meshkernel::constants::missing::doubleValue;
    landBoundaries.coordinates_x = xCoordinates.data();
    landBoundaries.coordinates_y = yCoordinates.data();
    landBoundaries.values = zCoordinates.data();
    landBoundaries.num_coordinates = static_cast<int>(xCoordinates.size());

    auto errorCode = mkernel_mesh2d_flip_edges(meshKernelId,
                                               isTriangulationRequired,
                                               projectToLandBoundaryOption,
                                               selectingPolygon,
                                               landBoundaries);

    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(12, mesh2d.num_nodes);
    ASSERT_EQ(23, mesh2d.num_edges);
}

TEST_F(CartesianApiTests, InsertEdgeThroughApi)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    int newEdgeIndex;
    auto errorCode = meshkernelapi::mkernel_mesh2d_insert_edge(meshKernelId, 0, 4, newEdgeIndex);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(17, newEdgeIndex);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(12, mesh2d.num_nodes);
    ASSERT_EQ(18, mesh2d.num_edges);
}

TEST_F(CartesianApiTests, MergeTwoNodesThroughApi)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    auto errorCode = meshkernelapi::mkernel_mesh2d_merge_two_nodes(meshKernelId, 0, 4);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(11, mesh2d.num_nodes);
    ASSERT_EQ(15, mesh2d.num_edges);
}

TEST_F(CartesianApiTests, MergeNodesThroughApi)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();
    meshkernelapi::GeometryList geometry_list{};

    // Execute
    auto errorCode = mkernel_mesh2d_merge_nodes(meshKernelId, geometry_list, 0.001);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);

    // Assert (nothing is done, just check that the api communication works)
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(12, mesh2d.num_nodes);
    ASSERT_EQ(17, mesh2d.num_edges);
}

TEST_F(CartesianApiTests, OrthogonalizationThroughApi)
{
    // Set a new mesh in mesh
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Prepare
    meshkernel::OrthogonalizationParameters orthogonalizationParameters{};
    orthogonalizationParameters.outer_iterations = 1;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;

    meshkernelapi::GeometryList geometryList{};
    meshkernelapi::GeometryList landBoundaries{};

    // Execute
    auto errorCode = mkernel_mesh2d_initialize_orthogonalization(meshKernelId,
                                                                 1,
                                                                 orthogonalizationParameters,
                                                                 landBoundaries,
                                                                 geometryList);

    // Assert (nothing is done, just check that the api communication works)
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_prepare_outer_iteration_orthogonalization(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_compute_inner_ortogonalization_iteration(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_finalize_inner_ortogonalization_iteration(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_delete_orthogonalization(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);

    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(12, mesh2d.num_nodes);
    ASSERT_EQ(17, mesh2d.num_edges);
}

TEST_F(CartesianApiTests, GenerateTriangularGridThroughApi)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometry_separator = meshkernel::constants::missing::doubleValue;
    std::vector xCoordinates{
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
        415.319672};

    std::vector yCoordinates{
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
        490.293762};

    std::vector zCoordinates{
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
        0.0};

    geometryListIn.coordinates_x = xCoordinates.data();
    geometryListIn.coordinates_y = yCoordinates.data();
    geometryListIn.values = zCoordinates.data();
    geometryListIn.num_coordinates = static_cast<int>(xCoordinates.size());

    // Execute
    auto errorCode = mkernel_mesh2d_make_mesh_from_polygon(meshKernelId, geometryListIn);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    // Get the new state

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(42, mesh2d.num_nodes);
    ASSERT_EQ(107, mesh2d.num_edges);
}

TEST_F(CartesianApiTests, GenerateTriangularGridFromSamplesThroughApi)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;

    geometryListIn.geometry_separator = meshkernel::constants::missing::doubleValue;

    std::vector xCoordinates{
        0.0,
        10.0,
        10.0,
        0.0,
        0.0};

    std::vector yCoordinates{
        0.0,
        0.0,
        10.0,
        10.0,
        0.0};

    std::vector zCoordinates{
        0.0,
        0.0,
        0.0,
        0.0,
        0.0};

    geometryListIn.coordinates_x = xCoordinates.data();
    geometryListIn.coordinates_y = yCoordinates.data();
    geometryListIn.values = zCoordinates.data();

    geometryListIn.num_coordinates = static_cast<int>(xCoordinates.size());

    // Execute
    auto errorCode = mkernel_mesh2d_make_mesh_from_samples(meshKernelId, geometryListIn);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(4, mesh2d.num_nodes);
    ASSERT_EQ(5, mesh2d.num_edges);
}

TEST_F(CartesianApiTests, GetMeshBoundariesThroughApi)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();
    int numberOfpolygonNodes;
    auto errorCode = meshkernelapi::mkernel_mesh2d_count_mesh_boundaries_as_polygons(meshKernelId, numberOfpolygonNodes);
    ASSERT_EQ(11, numberOfpolygonNodes);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::GeometryList geometryListOut;
    geometryListOut.geometry_separator = meshkernel::constants::missing::doubleValue;
    geometryListOut.num_coordinates = numberOfpolygonNodes;

    std::vector<double> xCoordinates(numberOfpolygonNodes);
    std::vector<double> yCoordinates(numberOfpolygonNodes);
    std::vector<double> zCoordinates(numberOfpolygonNodes);

    geometryListOut.coordinates_x = xCoordinates.data();
    geometryListOut.coordinates_y = yCoordinates.data();
    geometryListOut.values = zCoordinates.data();

    // Execute
    errorCode = mkernel_mesh2d_get_mesh_boundaries_as_polygons(meshKernelId, geometryListOut);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    const double tolerance = 1e-6;
    ASSERT_NEAR(0.0, geometryListOut.coordinates_x[0], tolerance);
    ASSERT_NEAR(0.0, geometryListOut.coordinates_y[0], tolerance);
}

TEST_F(CartesianApiTests, OffsetAPolygonThroughApi)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometry_separator = meshkernel::constants::missing::doubleValue;

    std::vector xCoordinatesIn{0.0, 1.0, 1.0, 0.0};
    std::vector yCoordinatesIn{0.0, 0.0, 1.0, 1.0};
    std::vector valuesIn{0.0, 0.0, 0.0, 0.0};

    geometryListIn.coordinates_x = xCoordinatesIn.data();
    geometryListIn.coordinates_y = yCoordinatesIn.data();
    geometryListIn.values = valuesIn.data();
    geometryListIn.num_coordinates = static_cast<int>(xCoordinatesIn.size());

    // Execute
    int numberOfpolygonNodes;
    auto errorCode = mkernel_polygon_count_offset(meshKernelId, geometryListIn, false, 0.5, numberOfpolygonNodes);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(4, numberOfpolygonNodes);

    meshkernelapi::GeometryList geometryListOut;

    geometryListOut.geometry_separator = meshkernel::constants::missing::doubleValue;

    std::vector<double> xCoordinatesOut(numberOfpolygonNodes);
    std::vector<double> yCoordinatesOut(numberOfpolygonNodes);
    std::vector<double> valuesOut(numberOfpolygonNodes);
    geometryListOut.coordinates_x = xCoordinatesOut.data();
    geometryListOut.coordinates_y = yCoordinatesOut.data();
    geometryListOut.values = valuesOut.data();
    geometryListOut.num_coordinates = static_cast<int>(xCoordinatesOut.size());

    errorCode = mkernel_polygon_get_offset(meshKernelId, geometryListIn, false, 10.0, geometryListOut);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    const double tolerance = 1e-6;
    ASSERT_NEAR(0.0, geometryListOut.coordinates_x[0], tolerance);
    ASSERT_NEAR(-10.0, geometryListOut.coordinates_y[0], tolerance);
}

TEST_F(CartesianApiTests, RefineAPolygonThroughApi)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometry_separator = meshkernel::constants::missing::doubleValue;

    std::vector xCoordinatesIn{76.251099, 498.503723, 505.253784};
    std::vector yCoordinatesIn{92.626556, 91.126541, 490.130554};
    std::vector valuesIn{0.0, 0.0, 0.0};

    geometryListIn.coordinates_x = xCoordinatesIn.data();
    geometryListIn.coordinates_y = yCoordinatesIn.data();
    geometryListIn.values = valuesIn.data();
    geometryListIn.num_coordinates = static_cast<int>(xCoordinatesIn.size());

    // Execute
    int numberOfpolygonNodes;
    auto errorCode = mkernel_polygon_count_refine(meshKernelId, geometryListIn, 0, 2, 40, numberOfpolygonNodes);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(22, numberOfpolygonNodes);

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
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    const double tolerance = 1e-6;
    ASSERT_NEAR(76.251099, geometryListOut.coordinates_x[0], tolerance);
    ASSERT_NEAR(92.626556, geometryListOut.coordinates_y[0], tolerance);
}

TEST_F(CartesianApiTests, RefineBasedOnSamples_OnAUniformMesh_shouldRefineMesh)
{
    // Prepare
    MakeMesh();
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
    meshRefinementParameters.refinement_type = 3;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.account_for_samples_outside = 0;

    // Execute
    auto errorCode = mkernel_mesh2d_refine_based_on_samples(meshKernelId, geometryListIn, 1.0, 1, meshRefinementParameters);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(12, mesh2d.num_nodes);
    ASSERT_EQ(17, mesh2d.num_edges);
}

TEST_F(CartesianApiTests, RefineBasedOnGebcoSamples_OnAUniformMesh_shouldRefineMesh)
{
    // Prepare
    MakeMesh();
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
    meshRefinementParameters.refinement_type = 3;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.account_for_samples_outside = 0;

    // Execute
    auto errorCode = mkernel_mesh2d_refine_based_on_samples(meshKernelId, geometryListIn, 1.0, 1, meshRefinementParameters);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(12, mesh2d.num_nodes);
    ASSERT_EQ(17, mesh2d.num_edges);
}

TEST_F(CartesianApiTests, RefineAGridBasedOnPolygonThroughApi)
{
    // Prepare
    MakeMesh();
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
    geometryListIn.num_coordinates = static_cast<int>(xCoordinatesIn.size());

    meshkernel::MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 2;
    meshRefinementParameters.refine_intersected = 0;

    // Execute
    auto errorCode = mkernel_mesh2d_refine_based_on_polygon(meshKernelId, geometryListIn, meshRefinementParameters);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(12, mesh2d.num_nodes);
    ASSERT_EQ(17, mesh2d.num_edges);
}

TEST_F(CartesianApiTests, ComputeSingleContactsThroughApi_ShouldGenerateContacts)
{
    // Prepare
    MakeMesh(3, 3, 10);
    auto const meshKernelId = GetMeshKernelId();

    // Init 1d mesh
    meshkernelapi::Mesh1D mesh1d;
    std::vector node_x{
        1.73493900000000,
        2.35659313023165,
        5.38347452702839,
        14.2980910429074,
        22.9324017677239,
        25.3723169493137,
        25.8072280000000};
    std::vector node_y{
        -7.6626510000000,
        1.67281447902331,
        10.3513746546384,
        12.4797224193970,
        15.3007317677239,
        24.1623588554512,
        33.5111870000000};
    mesh1d.node_x = node_x.data();
    mesh1d.node_y = node_y.data();
    mesh1d.num_nodes = static_cast<int>(node_x.size());

    std::vector edge_nodes{0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6};
    mesh1d.edge_nodes = edge_nodes.data();
    mesh1d.num_edges = 6;

    auto errorCode = mkernel_mesh1d_set(meshKernelId, mesh1d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Init 1d mask
    std::vector onedNodeMask{1, 1, 1, 1, 1, 1, 1};

    // Init polygon
    meshkernelapi::GeometryList polygon;
    polygon.geometry_separator = meshkernel::constants::missing::doubleValue;

    std::vector<double> xCoordinates{-30, 40, 40, -40, -30};
    std::vector<double> yCoordinates{-20, -20, 50, 50, -20};
    std::vector<double> zCoordinates{0, 0, 0, 0, 0};
    polygon.coordinates_x = xCoordinates.data();
    polygon.coordinates_y = yCoordinates.data();
    polygon.values = zCoordinates.data();

    polygon.num_coordinates = static_cast<int>(xCoordinates.size());

    // Execute
    errorCode = mkernel_contacts_compute_single(meshKernelId, onedNodeMask.data(), polygon, 0.0);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state
    meshkernelapi::Contacts contacts{};
    errorCode = mkernel_contacts_get_dimensions(meshKernelId, contacts);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    std::vector<int> mesh1d_indices(contacts.num_contacts);
    std::vector<int> mesh2d_indices(contacts.num_contacts);
    contacts.mesh1d_indices = mesh1d_indices.data();
    contacts.mesh2d_indices = mesh2d_indices.data();

    errorCode = mkernel_contacts_get_data(meshKernelId, contacts);
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

TEST_F(CartesianApiTests, ComputeMultipleContactsThroughApi)
{
    // Prepare
    MakeMesh(3, 3, 10);
    auto const meshKernelId = GetMeshKernelId();

    // Init 1d mesh
    meshkernelapi::Mesh1D mesh1d;
    std::vector node_x{
        1.73493900000000,
        2.35659313023165,
        5.38347452702839,
        14.2980910429074,
        22.9324017677239,
        25.3723169493137,
        25.8072280000000};
    std::vector node_y{
        -7.6626510000000,
        1.67281447902331,
        10.3513746546384,
        12.4797224193970,
        15.3007317677239,
        24.1623588554512,
        33.5111870000000};
    mesh1d.node_x = node_x.data();
    mesh1d.node_y = node_y.data();
    mesh1d.num_nodes = 7;

    std::vector edge_nodes{0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6};
    mesh1d.edge_nodes = edge_nodes.data();
    mesh1d.num_edges = 6;

    auto errorCode = mkernel_mesh1d_set(meshKernelId, mesh1d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Init 1d mask
    std::vector onedNodeMask{1, 1, 1, 1, 1, 1, 1};

    // Execute
    // Why do we need to specify the namespace "meshkernelapi"?
    errorCode = meshkernelapi::mkernel_contacts_compute_multiple(meshKernelId, onedNodeMask.data());
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state
    meshkernelapi::Contacts contacts{};
    errorCode = mkernel_contacts_get_dimensions(meshKernelId, contacts);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    std::vector<int> mesh1d_indices(contacts.num_contacts);
    std::vector<int> mesh2d_indices(contacts.num_contacts);
    contacts.mesh1d_indices = mesh1d_indices.data();
    contacts.mesh2d_indices = mesh2d_indices.data();

    errorCode = mkernel_contacts_get_data(meshKernelId, contacts);
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

TEST_F(CartesianApiTests, ComputeContactsWithPolygonsThroughApi)
{
    // Prepare
    MakeMesh(3, 3, 10);
    auto const meshKernelId = GetMeshKernelId();

    // Init 1d mesh
    meshkernelapi::Mesh1D mesh1d;
    std::vector node_x{
        1.73493900000000,
        2.35659313023165,
        5.38347452702839,
        14.2980910429074,
        22.9324017677239,
        25.3723169493137,
        25.8072280000000};
    std::vector node_y{
        -7.6626510000000,
        1.67281447902331,
        10.3513746546384,
        12.4797224193970,
        15.3007317677239,
        24.1623588554512,
        33.5111870000000};
    mesh1d.node_x = node_x.data();
    mesh1d.node_y = node_y.data();
    mesh1d.num_nodes = 7;

    std::vector edge_nodes{
        0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6};
    mesh1d.edge_nodes = edge_nodes.data();
    mesh1d.num_edges = 6;

    auto errorCode = mkernel_mesh1d_set(meshKernelId, mesh1d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Init 1d mask
    std::vector onedNodeMask{1, 1, 1, 1, 1, 1, 1};

    // Init polygon
    meshkernelapi::GeometryList polygon;
    polygon.geometry_separator = meshkernel::constants::missing::doubleValue;

    std::vector<double> xCoordinates{25, 50, 50, 25, 25};
    std::vector<double> yCoordinates{25, 25, 50, 50, 25};
    std::vector<double> zCoordinates{0, 0, 0, 0, 0};
    polygon.coordinates_x = xCoordinates.data();
    polygon.coordinates_y = yCoordinates.data();
    polygon.values = zCoordinates.data();
    polygon.num_coordinates = static_cast<int>(xCoordinates.size());

    // Execute
    errorCode = mkernel_contacts_compute_with_polygons(meshKernelId, onedNodeMask.data(), polygon);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state
    meshkernelapi::Contacts contacts{};
    errorCode = mkernel_contacts_get_dimensions(meshKernelId, contacts);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    std::vector<int> mesh1d_indices(contacts.num_contacts);
    std::vector<int> mesh2d_indices(contacts.num_contacts);
    contacts.mesh1d_indices = mesh1d_indices.data();
    contacts.mesh2d_indices = mesh2d_indices.data();

    errorCode = mkernel_contacts_get_data(meshKernelId, contacts);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(1, contacts.num_contacts);
    ASSERT_EQ(5, contacts.mesh1d_indices[0]);
    ASSERT_EQ(8, contacts.mesh2d_indices[0]);
}

TEST_F(CartesianApiTests, ComputeContactsWithPointsThroughApi)
{
    // Prepare
    MakeMesh(3, 3, 10);
    auto const meshKernelId = GetMeshKernelId();

    // Init 1d mesh
    meshkernelapi::Mesh1D mesh1d;
    std::vector node_x{
        1.73493900000000,
        2.35659313023165,
        5.38347452702839,
        14.2980910429074,
        22.9324017677239,
        25.3723169493137,
        25.8072280000000};
    std::vector node_y{
        -7.6626510000000,
        1.67281447902331,
        10.3513746546384,
        12.4797224193970,
        15.3007317677239,
        24.1623588554512,
        33.5111870000000};
    mesh1d.node_x = node_x.data();
    mesh1d.node_y = node_y.data();
    mesh1d.num_nodes = static_cast<int>(node_x.size());

    std::vector edge_nodes{
        0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6};
    mesh1d.edge_nodes = edge_nodes.data();
    mesh1d.num_edges = 6;

    auto errorCode = mkernel_mesh1d_set(meshKernelId, mesh1d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Init 1d mask
    std::vector onedNodeMask{1, 1, 1, 1, 1, 1, 1};

    // Init polygon
    meshkernelapi::GeometryList points;
    points.geometry_separator = meshkernel::constants::missing::doubleValue;

    std::vector<double> xCoordinates{5, 15, 25, 35};
    std::vector<double> yCoordinates{5, 15, 25, 35};
    std::vector<double> zCoordinates{0, 0, 0, 0};
    points.coordinates_x = xCoordinates.data();
    points.coordinates_y = yCoordinates.data();
    points.values = zCoordinates.data();
    points.num_coordinates = static_cast<int>(xCoordinates.size());

    // Execute
    errorCode = mkernel_contacts_compute_with_points(meshKernelId, onedNodeMask.data(), points);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state
    meshkernelapi::Contacts contacts{};
    errorCode = mkernel_contacts_get_dimensions(meshKernelId, contacts);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    std::vector<int> mesh1d_indices(contacts.num_contacts);
    std::vector<int> mesh2d_indices(contacts.num_contacts);
    contacts.mesh1d_indices = mesh1d_indices.data();
    contacts.mesh2d_indices = mesh2d_indices.data();

    errorCode = mkernel_contacts_get_data(meshKernelId, contacts);
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

TEST_F(CartesianApiTests, ComputeBoundaryContactsThroughApi)
{
    // Prepare
    MakeMesh(3, 3, 10);
    auto const meshKernelId = GetMeshKernelId();

    // Init 1d mesh
    meshkernelapi::Mesh1D mesh1d;
    std::vector node_x{
        -16.1886410000000,
        -16.1464995876014,
        -16.1043581752028,
        -16.0622167628042,
        -15.7539488236928,
        -6.86476658679268,
        2.02441565010741,
        10.9135970000000,
    };
    std::vector node_y{
        0.89018900000000,
        9.78201442138723,
        18.6738398427745,
        27.5656652641617,
        36.1966603330179,
        36.4175095626911,
        36.6383587923643,
        36.8592080000000};
    mesh1d.node_x = node_x.data();
    mesh1d.node_y = node_y.data();
    mesh1d.num_nodes = 8;

    std::vector edge_nodes{0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7};
    mesh1d.edge_nodes = edge_nodes.data();
    mesh1d.num_edges = 7;

    auto errorCode = mkernel_mesh1d_set(meshKernelId, mesh1d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Init 1d mask
    std::vector onedNodeMask{1, 1, 1, 1, 1, 1, 1, 1};

    // Init polygon
    meshkernelapi::GeometryList polygon;
    polygon.geometry_separator = meshkernel::constants::missing::doubleValue;

    std::vector<double> xCoordinates{-30, 40, 40, -40, -30};
    std::vector<double> yCoordinates{-20, -20, 50, 50, -20};
    std::vector<double> zCoordinates{0, 0, 0, 0, 0};
    polygon.coordinates_x = xCoordinates.data();
    polygon.coordinates_y = yCoordinates.data();
    polygon.values = zCoordinates.data();
    polygon.num_coordinates = static_cast<int>(xCoordinates.size());

    // Execute
    errorCode = mkernel_contacts_compute_boundary(meshKernelId, onedNodeMask.data(), polygon, 200.0);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state
    meshkernelapi::Contacts contacts{};
    errorCode = mkernel_contacts_get_dimensions(meshKernelId, contacts);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    std::vector mesh1d_indices(contacts.num_contacts, 0);
    std::vector mesh2d_indices(contacts.num_contacts, 0);
    contacts.mesh1d_indices = mesh1d_indices.data();
    contacts.mesh2d_indices = mesh2d_indices.data();

    errorCode = mkernel_contacts_get_data(meshKernelId, contacts);
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

    std::vector<double> splineCoordinatesX{10.0, 20.0, 30.0};
    std::vector<double> splineCoordinatesY{-5.0, 5.0, -5.0};
    std::vector<double> values{0.0, 0.0, 0.0};
    geometryListIn.coordinates_x = splineCoordinatesX.data();
    geometryListIn.coordinates_y = splineCoordinatesY.data();
    geometryListIn.values = values.data();
    geometryListIn.num_coordinates = static_cast<int>(splineCoordinatesX.size());
    geometryListIn.geometry_separator = meshkernel::constants::missing::doubleValue;

    meshkernelapi::GeometryList geometryListOut;
    int const numberOfPointsBetweenNodes = 3;
    size_t totalNumPoints = (numberOfPointsBetweenNodes + 1) * 2 + 1;
    std::vector<double> CoordinatesOutX(totalNumPoints);
    std::vector<double> CoordinatesOutY(totalNumPoints);
    std::vector<double> valuesOut(totalNumPoints);
    geometryListOut.coordinates_x = CoordinatesOutX.data();
    geometryListOut.coordinates_y = CoordinatesOutY.data();
    geometryListOut.values = valuesOut.data();
    geometryListOut.num_coordinates = static_cast<int>(totalNumPoints);
    geometryListOut.geometry_separator = meshkernel::constants::missing::doubleValue;

    // Execute
    auto errorCode = mkernel_get_splines(geometryListIn, geometryListOut, numberOfPointsBetweenNodes);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert. The last value is a geometry separator  to handle the case of multiple splines
    ASSERT_EQ(totalNumPoints, geometryListOut.num_coordinates);

    std::vector computedCoordinatesX(CoordinatesOutX.data(), CoordinatesOutX.data() + totalNumPoints);
    std::vector ValidCoordinatesX{10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0};
    ASSERT_THAT(computedCoordinatesX, ::testing::ContainerEq(ValidCoordinatesX));

    std::vector computedCoordinatesY(CoordinatesOutY.data(), CoordinatesOutY.data() + totalNumPoints);
    std::vector ValidCoordinatesY{-5.000000, -1.328125, 1.8750000, 4.1406250, 5.0000000, 4.1406250, 1.8750000, -1.328125, -5.000000};
    ASSERT_THAT(computedCoordinatesY, ::testing::ContainerEq(ValidCoordinatesY));
}

TEST(ApiStatelessTests, Orthogonalize_OnInvaliMesh_ShouldThrowAMeshGeometryError)
{
    // Prepare
    int meshKernelId;
    int isGeographic = 0;
    meshkernelapi::mkernel_allocate_state(isGeographic, meshKernelId);

    auto [num_nodes, num_edges, node_x, node_y, node_type, edge_nodes, edge_type] = ReadLegacyMeshFile(TEST_FOLDER + "/data/InvalidMeshes/invalid_orthogonalization_net.nc");
    meshkernelapi::Mesh2D mesh2d;
    mesh2d.num_edges = static_cast<int>(num_edges);
    mesh2d.num_nodes = static_cast<int>(num_nodes);
    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_nodes = edge_nodes.data();

    auto errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernel::OrthogonalizationParameters orthogonalizationParameters{};
    orthogonalizationParameters.outer_iterations = 1;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;

    meshkernelapi::GeometryList geometryList{};
    meshkernelapi::GeometryList landBoundaries{};

    // Execute
    errorCode = mkernel_mesh2d_initialize_orthogonalization(meshKernelId,
                                                            1,
                                                            orthogonalizationParameters,
                                                            landBoundaries,
                                                            geometryList);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert there is a geometry error
    errorCode = meshkernelapi::mkernel_mesh2d_prepare_outer_iteration_orthogonalization(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::MeshGeometryError, errorCode);

    // Delete orthogonalization instance
    errorCode = meshkernelapi::mkernel_mesh2d_delete_orthogonalization(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the message
    auto exceptionMessage = std::make_unique<char>(512);
    errorCode = meshkernelapi::mkernel_get_error(exceptionMessage.get());
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the index of the invalid location
    int invalidIndex;
    int type;
    errorCode = meshkernelapi::mkernel_get_geometry_error(invalidIndex, type);
    ASSERT_EQ(static_cast<int>(meshkernel::Mesh::Location::Nodes), type);
    ASSERT_EQ(478, invalidIndex);
}

TEST(ApiStatelessTests, TestGettingVersionThroughApi)
{
    auto versionFromApi = std::make_unique<char>(64);
    meshkernelapi::mkernel_get_version(versionFromApi.get());
    ASSERT_EQ(strcmp(versionFromApi.get(), versionString), 0);
}

TEST_F(CartesianApiTests, CurvilinearComputeTransfiniteFromPolygon_ShouldComputeAValidCurvilinearGrid)
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
    geometryListIn.geometry_separator = meshkernel::constants::missing::doubleValue;
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
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state
    meshkernelapi::CurvilinearGrid curvilinear_grid{};

    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinear_grid);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(3, curvilinear_grid.num_m);
    ASSERT_EQ(3, curvilinear_grid.num_n);
}

TEST_F(CartesianApiTests, GetClosestMeshCoordinateThroughApi)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    double xCoordinatesOut, yCoordinatesOut;
    auto errorCode = meshkernelapi::mkernel_mesh2d_get_closest_node(meshKernelId, -5.0, 5.0, 10.0, 0, 0, 0, 0, xCoordinatesOut, yCoordinatesOut);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(0.0, xCoordinatesOut);
}

TEST_F(CartesianApiTests, MakeCurvilinearGridFromTriangleThroughApi)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometry_separator = meshkernel::constants::missing::doubleValue;
    std::vector xCoordinatesIn{
        444.504791,
        427.731781,
        405.640503,
        381.094666,
        451.050354,
        528.778931,
        593.416260,
        558.643005,
        526.733398,
        444.095703};
    std::vector yCoordinatesIn{
        437.155945,
        382.745758,
        317.699005,
        262.470612,
        262.879700,
        263.288788,
        266.561584,
        324.653687,
        377.836578,
        436.746857};
    std::vector valuesIn{
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0};
    geometryListIn.coordinates_x = xCoordinatesIn.data();
    geometryListIn.coordinates_y = yCoordinatesIn.data();
    geometryListIn.values = valuesIn.data();
    geometryListIn.num_coordinates = static_cast<int>(xCoordinatesIn.size());

    // Execute
    auto errorCode = mkernel_curvilinear_compute_transfinite_from_triangle(meshKernelId,
                                                                           geometryListIn,
                                                                           0,
                                                                           3,
                                                                           6);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = meshkernelapi::mkernel_curvilinear_convert_to_mesh2d(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);

    // Assert
    ASSERT_EQ(28, mesh2d.num_nodes);
    ASSERT_EQ(40, mesh2d.num_edges);
}

TEST_F(CartesianApiTests, MakeCurvilinearGridThroughApi)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();

    meshkernel::MakeGridParameters makeGridParameters;
    meshkernelapi::GeometryList geometryList{};

    makeGridParameters.num_columns = 3;
    makeGridParameters.num_rows = 2;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.block_size_x = 1.0;
    makeGridParameters.block_size_y = 1.0;

    // Execute
    auto errorCode = mkernel_curvilinear_make_uniform(meshKernelId,
                                                      makeGridParameters,
                                                      geometryList);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(3, curvilinearGrid.num_m);
    ASSERT_EQ(4, curvilinearGrid.num_n);

    // Allocate memory and get data
    std::vector<double> node_x(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> node_y(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = node_x.data();
    curvilinearGrid.node_y = node_y.data();
    errorCode = mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
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
}

TEST_F(CartesianApiTests, GenerateTransfiniteCurvilinearGridThroughApi)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometry_separator = meshkernel::constants::missing::doubleValue;
    std::vector xCoordinates{1.340015E+02, 3.642529E+02, 6.927549E+02, meshkernel::constants::missing::doubleValue,
                             2.585022E+02, 4.550035E+02, 8.337558E+02, meshkernel::constants::missing::doubleValue,
                             1.002513E+02, 4.610035E+02, meshkernel::constants::missing::doubleValue,
                             6.522547E+02, 7.197551E+02};

    std::vector yCoordinates{2.546282E+02, 4.586302E+02, 5.441311E+02, meshkernel::constants::missing::doubleValue,
                             6.862631E+01, 2.726284E+02, 3.753794E+02, meshkernel::constants::missing::doubleValue,
                             4.068797E+02, 7.912642E+01, meshkernel::constants::missing::doubleValue,
                             6.026317E+02, 2.681283E+02};

    std::vector zCoordinates{0.0, 0.0, 0.0, meshkernel::constants::missing::doubleValue,
                             0.0, 0.0, 0.0, meshkernel::constants::missing::doubleValue,
                             0.0, 0.0, meshkernel::constants::missing::doubleValue,
                             0.0, 0.0};

    geometryListIn.coordinates_x = xCoordinates.data();
    geometryListIn.coordinates_y = yCoordinates.data();
    geometryListIn.values = zCoordinates.data();
    geometryListIn.num_coordinates = static_cast<int>(xCoordinates.size());

    meshkernel::CurvilinearParameters curvilinearParameters;

    curvilinearParameters.m_refinement = 10;
    curvilinearParameters.n_refinement = 10;
    curvilinearParameters.smoothing_iterations = 10;
    curvilinearParameters.smoothing_parameter = 0.5;
    curvilinearParameters.attraction_parameter = 0.0;

    // Execute
    auto errorCode = mkernel_curvilinear_compute_transfinite_from_splines(meshKernelId,
                                                                          geometryListIn,
                                                                          curvilinearParameters);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(11, curvilinearGrid.num_m);
    ASSERT_EQ(11, curvilinearGrid.num_n);
}

TEST_F(CartesianApiTests, GenerateOrthogonalCurvilinearGridThroughApi)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryListIn;

    geometryListIn.geometry_separator = meshkernel::constants::missing::doubleValue;
    std::vector xCoordinates{1.175014E+02, 3.755030E+02, 7.730054E+02, meshkernel::constants::missing::doubleValue,
                             4.100089E+01, 3.410027E+02};

    std::vector yCoordinates{2.437587E+01, 3.266289E+02, 4.563802E+02, meshkernel::constants::missing::doubleValue,
                             2.388780E+02, 2.137584E+01};

    std::vector zCoordinates{0.0, 0.0, 0.0, meshkernel::constants::missing::doubleValue,
                             0.0, 0.0};

    geometryListIn.coordinates_x = xCoordinates.data();
    geometryListIn.coordinates_y = yCoordinates.data();
    geometryListIn.values = zCoordinates.data();
    geometryListIn.num_coordinates = static_cast<int>(xCoordinates.size());

    meshkernel::CurvilinearParameters curvilinearParameters;
    curvilinearParameters.m_refinement = 40;
    curvilinearParameters.n_refinement = 10;
    meshkernel::SplinesToCurvilinearParameters splinesToCurvilinearParameters;
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
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Grow grid, from the second layer
    for (auto layer = 1; layer <= curvilinearParameters.n_refinement; ++layer)
    {
        errorCode = meshkernelapi::mkernel_curvilinear_iterate_orthogonal_grid_from_splines(meshKernelId, layer);
        ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    }

    // Puts the computed curvilinear mesh into the mesh state (unstructured mesh)
    errorCode = meshkernelapi::mkernel_curvilinear_refresh_orthogonal_grid_from_splines(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Delete the mesh curvilinearGridFromSplinesInstances vector entry
    errorCode = meshkernelapi::mkernel_curvilinear_delete_orthogonal_grid_from_splines(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(3, curvilinearGrid.num_m);
    ASSERT_EQ(7, curvilinearGrid.num_n);
}

TEST_F(CartesianApiTests, RefineCompute_OnCurvilinearGrid_ShouldRefine)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();

    MakeUniformCurvilinearGrid(3, 3);

    auto errorCode = meshkernelapi::mkernel_curvilinear_refine(meshKernelId, 10.0, 20.0, 20.0, 20.0, 10);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    ASSERT_EQ(4, curvilinearGrid.num_m);
    ASSERT_EQ(13, curvilinearGrid.num_n);
}

TEST_F(CartesianApiTests, DerefineCompute_OnCurvilinearGrid_ShouldDeRefine)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();

    MakeUniformCurvilinearGrid();

    // Execute
    auto errorCode = meshkernelapi::mkernel_curvilinear_derefine(meshKernelId, 10.0, 20.0, 30.0, 20.0);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    ASSERT_EQ(5, curvilinearGrid.num_m);
    ASSERT_EQ(4, curvilinearGrid.num_n);
}

TEST_F(CartesianApiTests, Orthogonalize_CurvilinearGrid_ShouldOrthogonalize)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();

    MakeUniformCurvilinearGrid(3, 3);
    // Move a node to make the grid non orthogonal
    auto errorCode = meshkernelapi::mkernel_curvilinear_move_node(meshKernelId, 10.0, 20.0, 18.0, 12.0);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernel::OrthogonalizationParameters orthogonalizationParameters{};
    orthogonalizationParameters.outer_iterations = 1;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;

    // Execute
    errorCode = meshkernelapi::mkernel_curvilinear_initialize_orthogonalize(meshKernelId, orthogonalizationParameters);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    errorCode = meshkernelapi::mkernel_curvilinear_set_block_orthogonalize(meshKernelId, 0.0, 0.0, 30.0, 30.0);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    errorCode = meshkernelapi::mkernel_curvilinear_orthogonalize(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

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

TEST_F(CartesianApiTests, Smoothing_CurvilinearGrid_ShouldSmooth)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();

    MakeUniformCurvilinearGrid();

    // Execute
    auto errorCode = meshkernelapi::mkernel_curvilinear_smoothing(meshKernelId, 10, 10.0, 20.0, 30.0, 20.0);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert (nothing changed)
    ASSERT_EQ(5, curvilinearGrid.num_m);
    ASSERT_EQ(5, curvilinearGrid.num_n);
}

TEST_F(CartesianApiTests, ComputedDirectionalSmooth_CurvilinearGrid_ShouldSmooth)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();

    MakeUniformCurvilinearGrid();

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

    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert (nothing changed)
    ASSERT_EQ(5, curvilinearGrid.num_m);
    ASSERT_EQ(5, curvilinearGrid.num_n);
}

TEST_F(CartesianApiTests, ComputedLineShift_CurvilinearGrid_ShouldShift)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();

    MakeUniformCurvilinearGrid();

    meshkernelapi::mkernel_curvilinear_initialize_line_shift(meshKernelId);

    auto errorCode = meshkernelapi::mkernel_curvilinear_set_line_line_shift(meshKernelId, 0.0, 0.0, 0.0, 30.0);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_set_block_line_shift(meshKernelId, 0.0, 0.0, 30.0, 30.0);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    /// Move a gridline point, in this case the origin to -10.0, 0.0
    errorCode = meshkernelapi::mkernel_curvilinear_move_node_line_shift(meshKernelId, 0.0, 0.0, -10.0, 0.0);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Execute
    errorCode = meshkernelapi::mkernel_curvilinear_line_shift(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    errorCode = meshkernelapi::mkernel_curvilinear_finalize_line_shift(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert, the nodes along the grid line have changed
    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

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

TEST_F(CartesianApiTests, DeleteMesh2D_WithEmptyPolygon_ShouldDeleteMesh2D)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList geometryList{};

    // Execute
    auto errorCode = mkernel_mesh2d_delete(meshKernelId, geometryList, 0, false);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);

    // Assert
    ASSERT_EQ(0, mesh2d.num_nodes);
    ASSERT_EQ(0, mesh2d.num_edges);
}

TEST_F(CartesianApiTests, GetDimensionsMesh1D_WithMesh1D_ShouldGetDimensionsMesh1D)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    std::vector nodesX{-16.1886410000000,
                       -16.1464995876014,
                       -16.1043581752028,
                       -16.0622167628042,
                       -15.7539488236928,
                       -6.86476658679268,
                       2.02441565010741,
                       10.9135970000000};

    std::vector nodesY{0.89018900000000,
                       9.78201442138723,
                       18.6738398427745,
                       27.5656652641617,
                       36.1966603330179,
                       36.4175095626911,
                       36.6383587923643,
                       36.8592080000000};

    std::vector edges{0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7};

    meshkernelapi::Mesh1D mesh1d;
    mesh1d.edge_nodes = edges.data();
    mesh1d.node_x = nodesX.data();
    mesh1d.node_y = nodesY.data();
    mesh1d.num_nodes = 8;
    mesh1d.num_edges = 7;

    auto errorCode = mkernel_mesh1d_set(GetMeshKernelId(), mesh1d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Execute
    meshkernelapi::Mesh1D mesh1dResults;
    errorCode = mkernel_mesh1d_get_dimensions(meshKernelId, mesh1dResults);

    // Assert
    ASSERT_EQ(8, mesh1dResults.num_nodes);
    ASSERT_EQ(7, mesh1dResults.num_edges);
}

TEST_F(CartesianApiTests, GetDataMesh1D_WithMesh1D_ShouldGetDataMesh1D)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    std::vector nodesX{-16.1886410000000,
                       -16.1464995876014,
                       -16.1043581752028,
                       -16.0622167628042,
                       -15.7539488236928,
                       -6.86476658679268,
                       2.02441565010741,
                       10.9135970000000};

    std::vector nodesY{0.89018900000000,
                       9.78201442138723,
                       18.6738398427745,
                       27.5656652641617,
                       36.1966603330179,
                       36.4175095626911,
                       36.6383587923643,
                       36.8592080000000};

    std::vector edges{0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7};

    meshkernelapi::Mesh1D mesh1d;
    mesh1d.edge_nodes = edges.data();
    mesh1d.node_x = nodesX.data();
    mesh1d.node_y = nodesY.data();
    mesh1d.num_nodes = 8;
    mesh1d.num_edges = 7;

    auto errorCode = mkernel_mesh1d_set(GetMeshKernelId(), mesh1d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Execute
    meshkernelapi::Mesh1D mesh1dResults;
    errorCode = mkernel_mesh1d_get_dimensions(meshKernelId, mesh1dResults);
    std::vector<double> meshNodesX(mesh1dResults.num_nodes);
    std::vector<double> meshNodesY(mesh1dResults.num_nodes);
    std::vector<int> meshEdges(mesh1dResults.num_edges * 2);

    mesh1dResults.node_x = meshNodesX.data();
    mesh1dResults.node_y = meshNodesY.data();
    mesh1dResults.edge_nodes = meshEdges.data();
    errorCode = mkernel_mesh1d_get_data(meshKernelId, mesh1dResults);

    // Assert
    std::vector validMeshNodesX(nodesX.data(), nodesX.data() + mesh1d.num_nodes);
    std::vector computedMeshNodesX(meshNodesX.data(), meshNodesX.data() + mesh1dResults.num_nodes);
    ASSERT_THAT(computedMeshNodesX, ::testing::ContainerEq(validMeshNodesX));

    std::vector validMeshNodesY(nodesY.data(), nodesY.data() + mesh1d.num_nodes);
    std::vector computedMeshNodesY(meshNodesY.data(), meshNodesY.data() + mesh1dResults.num_nodes);
    ASSERT_THAT(computedMeshNodesY, ::testing::ContainerEq(validMeshNodesY));

    std::vector<double> validEdges(edges.data(), edges.data() + mesh1d.num_edges);
    std::vector<double> computedEdges(meshEdges.data(), meshEdges.data() + mesh1dResults.num_edges);
    ASSERT_THAT(computedEdges, ::testing::ContainerEq(validEdges));
}

TEST_F(CartesianApiTests, CountHangingEdgesMesh2D_WithZeroHangingEdges_ShouldCountZeroEdges)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    int numHangingEdges;
    auto const errorCode = meshkernelapi::mkernel_mesh2d_count_hanging_edges(meshKernelId, numHangingEdges);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    ASSERT_EQ(0, numHangingEdges);
}

TEST_F(CartesianApiTests, GetHangingEdgesMesh2D_WithOneHangingEdges_ShouldGetOneHangingEdges)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // delete an edge at the lower left corner to create a new hanging edge
    auto errorCode = meshkernelapi::mkernel_mesh2d_delete_edge(meshKernelId, 0.5, 0.0, 0.0, 0.0, 1.0, 1.0);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    int numHangingEdges;
    errorCode = meshkernelapi::mkernel_mesh2d_count_hanging_edges(meshKernelId, numHangingEdges);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_GT(numHangingEdges, 0);

    // Execute
    std::vector<int> hangingEdges(numHangingEdges);
    errorCode = meshkernelapi::mkernel_mesh2d_get_hanging_edges(meshKernelId, hangingEdges.data());
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(hangingEdges[0], 8);
}

TEST_F(CartesianApiTests, DeleteHangingEdgesMesh2D_WithOneHangingEdges_ShouldDeleteOneHangingEdges)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // delete an edge at the lower left corner to create a new hanging edge
    auto errorCode = meshkernelapi::mkernel_mesh2d_delete_edge(meshKernelId, 0.5, 0.0, 0.0, 0.0, 1.0, 1.0);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Before deletion
    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(mesh2d.num_edges, 16);

    // Execute
    errorCode = meshkernelapi::mkernel_mesh2d_delete_hanging_edges(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);

    // Assert
    ASSERT_EQ(mesh2d.num_edges, 15);
}

TEST_F(CartesianApiTests, ComputeOrthogonalizationMesh2D_WithOrthogonalMesh2D_ShouldOrthogonalize)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    meshkernel::OrthogonalizationParameters orthogonalizationParameters;
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

    auto errorCode = mkernel_mesh2d_compute_orthogonalization(meshKernelId, 1, orthogonalizationParameters, polygons, landBoundaries);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
}

TEST_F(CartesianApiTests, GetOrthogonalityMesh2D_OnMesh2D_ShouldGetOrthogonality)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::Mesh2D mesh2d{};
    auto errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::GeometryList edgeOrthogonality;
    std::vector<double> values(mesh2d.num_edges);
    edgeOrthogonality.values = values.data();
    edgeOrthogonality.num_coordinates = mesh2d.num_edges;

    // Execute
    errorCode = mkernel_mesh2d_get_orthogonality(meshKernelId, edgeOrthogonality);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
}

TEST_F(CartesianApiTests, GetSmoothnessMesh2D_OnMesh2D_ShouldGetSmoothness)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::Mesh2D mesh2d{};
    auto errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::GeometryList edgeSmoothness;
    std::vector<double> values(mesh2d.num_edges);
    edgeSmoothness.values = values.data();
    edgeSmoothness.num_coordinates = mesh2d.num_edges;

    // Execute
    errorCode = mkernel_mesh2d_get_smoothness(meshKernelId, edgeSmoothness);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
}

TEST_F(CartesianApiTests, GetNodesInPolygonMesh2D_OnMesh2D_ShouldGetAllNodes)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // By using an empty list, all nodes will be selected
    const meshkernelapi::GeometryList geometryListIn{};

    meshkernelapi::Mesh2D mesh2d{};
    auto errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Execute
    std::vector<int> selectedNodes(mesh2d.num_nodes);
    errorCode = mkernel_mesh2d_get_nodes_in_polygons(meshKernelId,
                                                     geometryListIn,
                                                     1,
                                                     selectedNodes.data());

    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert (all nodes indices will be selected)
    std::vector actualResult(selectedNodes.data(), selectedNodes.data() + mesh2d.num_nodes);
    std::vector<int> expectedResult(mesh2d.num_nodes);
    std::iota(expectedResult.begin(), expectedResult.end(), 0);
    ASSERT_THAT(actualResult, ::testing::ContainerEq(expectedResult));
}

TEST_F(CartesianApiTests, CountNodesInPolygonMesh2D_OnMesh2D_ShouldCountAllNodes)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // By using an empty list, all nodes will be selected
    const meshkernelapi::GeometryList geometryListIn{};

    // Execute
    int numNodes;
    const auto errorCode = mkernel_mesh2d_count_nodes_in_polygons(meshKernelId,
                                                                  geometryListIn,
                                                                  1,
                                                                  numNodes);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert all nodes have been selected
    ASSERT_EQ(12, numNodes);
}

TEST_F(CartesianApiTests, InsertNodeAndEdge_OnMesh2D_ShouldInsertNodeAndEdge)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    // Isolated nodes are removed by the administration done in mkernel_mesh2d_get_dimensions.
    // The newly inserted node should be connected to another one to form an edge.
    // In this manner, the edge will not be removed during the administration
    int newNodeIndex;
    auto errorCode = meshkernelapi::mkernel_mesh2d_insert_node(meshKernelId, -0.5, -0.5, newNodeIndex);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    int newEdgeIndex;
    errorCode = meshkernelapi::mkernel_mesh2d_insert_edge(meshKernelId, newNodeIndex, 0, newEdgeIndex);

    // Assert
    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(mesh2d.num_nodes, 13);
    ASSERT_EQ(mesh2d.num_edges, 18);
}

TEST_F(CartesianApiTests, MoveNode_OnMesh2D_ShouldMoveNode)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    auto errorCode = meshkernelapi::mkernel_mesh2d_move_node(meshKernelId, -0.5, -0.5, 0);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    std::vector<int> edge_nodes(mesh2d.num_edges * 2);
    std::vector<int> face_nodes(mesh2d.num_face_nodes);
    std::vector<int> nodes_per_face(mesh2d.num_faces);

    std::vector<double> node_x(mesh2d.num_nodes);
    std::vector<double> node_y(mesh2d.num_nodes);

    std::vector<double> edge_x(mesh2d.num_edges);
    std::vector<double> edge_y(mesh2d.num_edges);

    std::vector<double> face_x(mesh2d.num_faces);
    std::vector<double> face_y(mesh2d.num_faces);

    mesh2d.edge_nodes = edge_nodes.data();
    mesh2d.face_nodes = face_nodes.data();
    mesh2d.nodes_per_face = nodes_per_face.data();
    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_x = edge_x.data();
    mesh2d.edge_y = edge_y.data();
    mesh2d.face_x = face_x.data();
    mesh2d.face_y = face_y.data();
    errorCode = mkernel_mesh2d_get_data(meshKernelId, mesh2d);

    // Assert
    ASSERT_EQ(mesh2d.node_x[0], -0.5);
    ASSERT_EQ(mesh2d.node_y[0], -0.5);
}

TEST_F(CartesianApiTests, MoveNode_OnMesh2DWithInvalidIndex_ShouldReturnAnErrorCode)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    const auto errorCode = meshkernelapi::mkernel_mesh2d_move_node(meshKernelId, -0.5, -0.5, -1);

    // Assert errorCode is not equal to MeshKernelApiErrors::Success
    ASSERT_NE(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
}

TEST_F(CartesianApiTests, GetEdge_OnMesh2D_ShouldGetAnEdgeIndex)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    int edgeIndex;
    const auto errorCode = meshkernelapi::mkernel_mesh2d_get_edge(meshKernelId, 0.5, -0.5, 0.0, 0.0, 3.0, 3.0, edgeIndex);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(edgeIndex, 0);
}

TEST_F(CartesianApiTests, GetNode_OnMesh2D_ShouldGetANodeIndex)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    int nodeIndex;
    const auto errorCode = meshkernelapi::mkernel_mesh2d_get_node_index(meshKernelId, 3.0, 3.0, 10.0, 0.0, 0.0, 3.0, 3.0, nodeIndex);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(nodeIndex, 11);
}

TEST_F(CartesianApiTests, CountSmallFlowEdges_OnMesh2D_ShouldCountSmallFlowEdges)
{
    // Prepare a mesh with two triangles
    meshkernelapi::Mesh2D mesh2d;

    std::vector node_x{0.0, 1.0, 1.0, 1.0};
    std::vector node_y{0.0, 0.0, 0.3, -0.3};
    std::vector edge_nodes{0, 3, 3, 1, 1, 0, 1, 2, 2, 0};

    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_nodes = edge_nodes.data();
    mesh2d.num_edges = static_cast<int>(edge_nodes.size() * 0.5);
    mesh2d.num_nodes = static_cast<int>(node_x.size());

    // Get the meshkernel id
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    auto errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    int numSmallFlowEdges;
    double const smallFlowEdgesThreshold = 100.0;
    errorCode = meshkernelapi::mkernel_mesh2d_count_small_flow_edge_centers(meshKernelId, smallFlowEdgesThreshold, numSmallFlowEdges);
    ASSERT_EQ(1, numSmallFlowEdges);
}

TEST_F(CartesianApiTests, GetSmallFlowEdges_OnMesh2D_ShouldGetSmallFlowEdges)
{
    // Prepare a mesh with two triangles
    meshkernelapi::Mesh2D mesh2d;

    std::vector<double> node_x{0.0, 1.0, 1.0, 1.0};
    std::vector<double> node_y{0.0, 0.0, 0.3, -0.3};
    std::vector<int> edge_nodes{0, 3, 3, 1, 1, 0, 1, 2, 2, 0};

    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_nodes = edge_nodes.data();
    mesh2d.num_edges = static_cast<int>(edge_nodes.size() * 0.5);
    mesh2d.num_nodes = static_cast<int>(node_x.size());

    // Get the meshkernel id
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    auto errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    int numSmallFlowEdges;
    double const smallFlowEdgesThreshold = 100.0;
    errorCode = meshkernelapi::mkernel_mesh2d_count_small_flow_edge_centers(meshKernelId, smallFlowEdgesThreshold, numSmallFlowEdges);

    // Assert
    std::vector<double> coordinates_x(numSmallFlowEdges);
    std::vector<double> coordinates_y(numSmallFlowEdges);
    meshkernelapi::GeometryList result{};
    result.coordinates_x = coordinates_x.data();
    result.coordinates_y = coordinates_y.data();
    result.num_coordinates = numSmallFlowEdges;

    errorCode = mkernel_mesh2d_get_small_flow_edge_centers(meshKernelId, smallFlowEdgesThreshold, result);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    const double tolerance = 1e-6;
    ASSERT_NEAR(result.coordinates_x[0], 0.5, tolerance);
    ASSERT_NEAR(result.coordinates_y[0], 0.0, tolerance);
}

TEST_F(CartesianApiTests, CountObtuseTriangles_OnMesh2DWithOneObtuseTriangle_ShouldCountObtuseTriangles)
{
    // Prepare a mesh with one obtuse triangle
    meshkernelapi::Mesh2D mesh2d;
    std::vector<double> coordinatesX{0.0, 3.0, -1.0, 1.5};
    std::vector<double> coordinatesY{0.0, 0.0, 2.0, -2.0};
    std::vector<int> edge_nodes{0, 1, 1, 2, 2, 0, 0, 3, 3, 1};
    mesh2d.node_x = coordinatesX.data();
    mesh2d.node_y = coordinatesY.data();
    mesh2d.edge_nodes = edge_nodes.data();
    mesh2d.num_edges = static_cast<int>(edge_nodes.size() * 0.5);
    mesh2d.num_nodes = static_cast<int>(coordinatesX.size());
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    auto errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    int numObtuseTriangles;
    errorCode = meshkernelapi::mkernel_mesh2d_count_obtuse_triangles(meshKernelId, numObtuseTriangles);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    ASSERT_EQ(1, numObtuseTriangles);
}

TEST_F(CartesianApiTests, Mesh2DCountObtuseTriangles_OnMesh2DWithOneObtuseTriangle_ShouldGetObtuseTriangle)
{
    // Prepare a mesh with one obtuse triangle
    meshkernelapi::Mesh2D mesh2d;
    std::vector<double> coordinatesX{0.0, 3.0, -1.0, 1.5};
    std::vector<double> coordinatesY{0.0, 0.0, 2.0, -2.0};
    std::vector<int> edge_nodes{0, 1, 1, 2, 2, 0, 0, 3, 3, 1};
    mesh2d.node_x = coordinatesX.data();
    mesh2d.node_y = coordinatesY.data();
    mesh2d.edge_nodes = edge_nodes.data();
    mesh2d.num_edges = static_cast<int>(edge_nodes.size() * 0.5);
    mesh2d.num_nodes = static_cast<int>(coordinatesX.size());
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    auto errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    int numObtuseTriangles;
    errorCode = meshkernelapi::mkernel_mesh2d_count_obtuse_triangles(meshKernelId, numObtuseTriangles);
    meshkernelapi::GeometryList geometryList{};

    std::vector<double> coordinatesObtuseTrianglesX(numObtuseTriangles);
    std::vector<double> coordinatesObtuseTrianglesY(numObtuseTriangles);
    geometryList.coordinates_x = coordinatesObtuseTrianglesX.data();
    geometryList.coordinates_y = coordinatesObtuseTrianglesY.data();
    geometryList.num_coordinates = numObtuseTriangles;
    errorCode = mkernel_mesh2d_get_obtuse_triangles_mass_centers(meshKernelId, geometryList);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(1, numObtuseTriangles);
    const double tolerance = 1e-6;
    std::vector computedCoordinatesX(coordinatesObtuseTrianglesX.data(), coordinatesObtuseTrianglesX.data() + numObtuseTriangles);
    ASSERT_NEAR(computedCoordinatesX[0], 0.66666666666666652, tolerance);
    std::vector computedCoordinatesY(coordinatesObtuseTrianglesY.data(), coordinatesObtuseTrianglesY.data() + numObtuseTriangles);
    ASSERT_NEAR(computedCoordinatesY[0], 0.66666666666666652, tolerance);
}

TEST_F(CartesianApiTests, Mesh2DDeleteSmallFlowEdgesAndSmallTriangles_OnMesh2DWithOneObtuseTriangle_ShouldDeleteOneEdge)
{
    // Prepare a mesh with one obtuse triangle
    meshkernelapi::Mesh2D mesh2d;
    std::vector coordinatesX{0.0, 3.0, -1.0, 1.5};
    std::vector coordinatesY{0.0, 0.0, 2.0, -2.0};
    std::vector edge_nodes{0, 1, 1, 2, 2, 0, 0, 3, 3, 1};

    mesh2d.node_x = coordinatesX.data();
    mesh2d.node_y = coordinatesY.data();
    mesh2d.edge_nodes = edge_nodes.data();
    mesh2d.num_nodes = static_cast<int>(coordinatesX.size());
    mesh2d.num_edges = static_cast<int>(edge_nodes.size() * 0.5);
    auto const meshKernelId = GetMeshKernelId();

    // Execute, with large length threshold
    auto errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    errorCode = meshkernelapi::mkernel_mesh2d_delete_small_flow_edges_and_small_triangles(meshKernelId, 1.0, 0.01);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    meshkernelapi::Mesh2D newMesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, newMesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // One edge is removed
    ASSERT_EQ(4, newMesh2d.num_edges);
}

TEST_F(CartesianApiTests, CurvilinearComputeOrthogonalGridFromSplines_ShouldMakeCurvilinearGrid)
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

    meshkernel::SplinesToCurvilinearParameters splinesToCurvilinearParameters;
    meshkernel::CurvilinearParameters curvilinearParameters{};

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
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert one curvilinear grid is produced
    ASSERT_GT(curvilinearGrid.num_m, 0);
    ASSERT_GT(curvilinearGrid.num_n, 0);
}

TEST_F(CartesianApiTests, CurvilinearSetFrozenLinesOrthogonalize_ShouldSetFrozenLines)
{
    // Setup
    MakeUniformCurvilinearGrid();
    auto const meshKernelId = GetMeshKernelId();
    meshkernel::OrthogonalizationParameters const orthogonalizationParameters{};

    auto errorCode = meshkernelapi::mkernel_curvilinear_initialize_orthogonalize(meshKernelId, orthogonalizationParameters);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Execute
    errorCode = meshkernelapi::mkernel_curvilinear_set_frozen_lines_orthogonalize(meshKernelId, 20.0, 0.0, 20.0, 10.0);

    // Asset
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
}

TEST_F(CartesianApiTests, CurvilinearFinalizeOrthogonalize_ShouldFinalize)
{
    // Setup
    MakeUniformCurvilinearGrid();
    auto const meshKernelId = GetMeshKernelId();
    meshkernel::OrthogonalizationParameters const orthogonalizationParameters{};

    auto errorCode = meshkernelapi::mkernel_curvilinear_initialize_orthogonalize(meshKernelId, orthogonalizationParameters);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Execute
    errorCode = meshkernelapi::mkernel_curvilinear_finalize_orthogonalize(meshKernelId);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
}

TEST_F(CartesianApiTests, CurvilinearInsertFace_ShouldInsertAFace)
{
    // Setup
    MakeUniformCurvilinearGrid();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    auto errorCode = meshkernelapi::mkernel_curvilinear_insert_face(meshKernelId, -5.0, 5.0);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    std::vector<double> node_x(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> node_y(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = node_x.data();
    curvilinearGrid.node_y = node_y.data();
    errorCode = mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert two extra nodes have been inserted (before it was 5 by 5 = 25 nodes, not it is 25 + 2 = 27)
    auto const numValidNodes = CurvilinearGridCountValidNodes(curvilinearGrid);
    ASSERT_EQ(numValidNodes, 27);
}

TEST_F(CartesianApiTests, CurvilinearLineMirror_ShouldInsertANewGridLine)
{
    // Setup
    MakeUniformCurvilinearGrid();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    auto errorCode = meshkernelapi::mkernel_curvilinear_line_mirror(meshKernelId,
                                                                    1.2,
                                                                    0.0,
                                                                    0.0,
                                                                    0.0,
                                                                    50.0);

    // Assert
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert that 5 nodes have been inserted on the bottom boundary, so the total count is 30
    ASSERT_EQ(curvilinearGrid.num_m * curvilinearGrid.num_n, 30);
}

TEST_F(CartesianApiTests, Mesh2dAveragingInterpolation_OnMesh2D_ShouldInterpolateValues)
{
    // Setup
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList samples{};
    std::vector<double> firstGridNodeCoordinateX{1.0, 2.0, 3.0, 1.0};
    std::vector<double> firstGridNodeCoordinateY{1.0, 3.0, 2.0, 4.0};
    std::vector<double> values{3.0, 10, 4.0, 5.0};

    samples.coordinates_x = firstGridNodeCoordinateX.data();
    samples.coordinates_y = firstGridNodeCoordinateY.data();
    samples.values = values.data();
    samples.num_coordinates = static_cast<int>(firstGridNodeCoordinateX.size());

    int const locationType = 1;             // Nodes
    int const averagingMethodType = 1;      // Simple averaging
    double const relativeSearchSize = 1.01; // The relative search size

    meshkernelapi::Mesh2D mesh2d;
    auto errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::GeometryList results{};
    std::vector<double> resultsCoordinateX(mesh2d.num_nodes);
    std::vector<double> resultsCoordinateY(mesh2d.num_nodes);
    std::vector<double> resultsValues(mesh2d.num_nodes);

    results.coordinates_x = resultsCoordinateX.data();
    results.coordinates_y = resultsCoordinateY.data();
    results.values = resultsValues.data();
    results.num_coordinates = static_cast<int>(resultsCoordinateX.size());

    // Execute
    errorCode = mkernel_mesh2d_averaging_interpolation(meshKernelId,
                                                       samples,
                                                       locationType,
                                                       averagingMethodType,
                                                       relativeSearchSize,
                                                       0,
                                                       results);

    // Assert the value has been interpolated
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    const double tolerance = 1e-6;
    std::vector computedResultsValues(resultsValues.data(), resultsValues.data() + mesh2d.num_nodes);
    ASSERT_NEAR(computedResultsValues[4], 3.0, tolerance);
}

TEST_F(CartesianApiTests, Mesh2dTriangulationInterpolation_ShouldInterpolateValues)
{
    // Setup
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::GeometryList samples{};
    std::vector<double> firstGridNodeCoordinateX{1.0, 2.0, 3.0, 1.0};
    std::vector<double> firstGridNodeCoordinateY{1.0, 3.0, 2.0, 4.0};
    std::vector<double> values{3.0, 10, 4.0, 5.0};

    samples.coordinates_x = firstGridNodeCoordinateX.data();
    samples.coordinates_y = firstGridNodeCoordinateY.data();
    samples.values = values.data();
    samples.num_coordinates = static_cast<int>(firstGridNodeCoordinateX.size());

    int const locationType = 1; // Nodes

    meshkernelapi::Mesh2D mesh2d;
    auto errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::GeometryList results{};
    std::vector<double> resultsCoordinateX(mesh2d.num_nodes);
    std::vector<double> resultsCoordinateY(mesh2d.num_nodes);
    std::vector<double> resultsValues(mesh2d.num_nodes);

    results.coordinates_x = resultsCoordinateX.data();
    results.coordinates_y = resultsCoordinateY.data();
    results.values = resultsValues.data();
    results.num_coordinates = static_cast<int>(resultsCoordinateX.size());

    // Execute
    errorCode = mkernel_mesh2d_triangulation_interpolation(meshKernelId,
                                                           samples,
                                                           locationType,
                                                           results);

    // Assert the value has been interpolated
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    const double tolerance = 1e-6;
    std::vector computedResultsValues(resultsValues.data(), resultsValues.data() + mesh2d.num_nodes);
    ASSERT_NEAR(computedResultsValues[8], 5.6666666666666670, tolerance);
}

TEST_F(CartesianApiTests, CurvilinearLineAttractionRepulsion_ShouldAttractGridlines)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();

    MakeUniformCurvilinearGrid(5, 5);

    // Execute
    auto errorCode = meshkernelapi::mkernel_curvilinear_line_attraction_repulsion(meshKernelId,
                                                                                  0.5,
                                                                                  30.0, 0.0, 30.0, 50.0,
                                                                                  10.0, 20.0, 50.0, 20.0);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    std::vector<double> node_x(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> node_y(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = node_x.data();
    curvilinearGrid.node_y = node_y.data();
    errorCode = mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert data
    const double tolerance = 1e-6;
    // Nodes
    ASSERT_NEAR(17.5, curvilinearGrid.node_x[2], tolerance);
    ASSERT_NEAR(0.0, curvilinearGrid.node_y[2], tolerance);
    ASSERT_NEAR(42.5, curvilinearGrid.node_x[4], tolerance);
    ASSERT_NEAR(0.0, curvilinearGrid.node_y[4], tolerance);
}

TEST_F(CartesianApiTests, CurvilinearDeleteNode_ShouldDeleteNode)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();
    MakeUniformCurvilinearGrid(5, 5);

    meshkernelapi::CurvilinearGrid curvilinearGrid{};
    auto errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
    std::vector<double> node_x(curvilinearGrid.num_m * curvilinearGrid.num_n);
    std::vector<double> node_y(curvilinearGrid.num_m * curvilinearGrid.num_n);
    curvilinearGrid.node_x = node_x.data();
    curvilinearGrid.node_y = node_y.data();
    errorCode = mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    auto const numValidNodesBefore = CurvilinearGridCountValidNodes(curvilinearGrid);

    // Execute
    errorCode = meshkernelapi::mkernel_curvilinear_delete_node(meshKernelId, 10.0, 0.0);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Asserts
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    errorCode = mkernel_curvilinear_get_data(meshKernelId, curvilinearGrid);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Two nodes are removed, one by delete node and one by the administration. The node is at the corner and the entire face will be removed
    auto const numValidNodesAfter = CurvilinearGridCountValidNodes(curvilinearGrid);
    ASSERT_EQ(numValidNodesBefore - 2, numValidNodesAfter);
}

TEST_F(CartesianApiTests, Network1DComputeFixedChainages_ShouldGenerateMesh1D)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();

    double separator = meshkernelapi::mkernel_get_separator();
    std::vector<double> polyLineXCoordinate{0.0, 10.0, 20.0, separator, 10.0, 10.0, 10.0};
    std::vector<double> polyLineYCoordinate{0.0, 0.0, 0.0, separator, -10.0, 0.0, 10.0};

    meshkernelapi::GeometryList polylines{};
    polylines.coordinates_x = polyLineXCoordinate.data();
    polylines.coordinates_y = polyLineYCoordinate.data();
    polylines.num_coordinates = static_cast<int>(polyLineXCoordinate.size());
    polylines.geometry_separator = separator;

    // Set the network
    auto errorCode = mkernel_network1d_set(meshKernelId, polylines);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Execute, compute fixed chainages
    int const fixedChainagesSize = 3;
    double const minFaceSize = 0.01;
    double const fixedChainagesOffset = 10.0;
    std::vector<double> fixedChainages{5.0, separator, 5.0};
    errorCode = meshkernelapi::mkernel_network1d_compute_fixed_chainages(meshKernelId, fixedChainages.data(), fixedChainagesSize, minFaceSize, fixedChainagesOffset);
    // Convert network 1d to mesh1d
    errorCode = meshkernelapi::mkernel_network1d_to_mesh1d(meshKernelId, minFaceSize);

    // Asserts
    meshkernelapi::Mesh1D mesh1dResults;
    errorCode = mkernel_mesh1d_get_dimensions(meshKernelId, mesh1dResults);
    ASSERT_EQ(6, mesh1dResults.num_nodes);
    ASSERT_EQ(4, mesh1dResults.num_edges);
}

TEST_F(CartesianApiTests, Network1DToMesh1d_FromPolylines_ShouldGenerateMesh1D)
{
    // Prepare
    auto const meshKernelId = GetMeshKernelId();

    double separator = meshkernelapi::mkernel_get_separator();
    std::vector<double> polyLineXCoordinate{0.0, 10.0, 20.0, separator, 10.0, 10.0, 10.0};
    std::vector<double> polyLineYCoordinate{0.0, 0.0, 0.0, separator, -10.0, 0.0, 10.0};
    meshkernelapi::GeometryList polylines{};
    polylines.coordinates_x = polyLineXCoordinate.data();
    polylines.coordinates_y = polyLineYCoordinate.data();
    polylines.num_coordinates = static_cast<int>(polyLineXCoordinate.size());
    polylines.geometry_separator = separator;

    // Set the network
    auto errorCode = mkernel_network1d_set(meshKernelId, polylines);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Execute, compute offsetted chainages
    double const offset = 1.0;
    errorCode = meshkernelapi::mkernel_network1d_compute_offsetted_chainages(meshKernelId, offset);

    // Convert network 1d to mesh1d
    errorCode = meshkernelapi::mkernel_network1d_to_mesh1d(meshKernelId, 0.01);

    // Asserts
    meshkernelapi::Mesh1D mesh1dResults;
    errorCode = mkernel_mesh1d_get_dimensions(meshKernelId, mesh1dResults);
    ASSERT_EQ(41, mesh1dResults.num_nodes);
    ASSERT_EQ(40, mesh1dResults.num_edges);
}

TEST(Mesh2D, Mesh2DInitializeOrthogonalization_WithHexagon_ShouldOrthogonalize)
{
    // Prepare
    int meshKernelId;
    const int isGeographic = 0;
    meshkernelapi::mkernel_allocate_state(isGeographic, meshKernelId);

    auto [num_nodes, num_edges, node_x, node_y, node_type, edge_nodes, edge_type] = ReadLegacyMeshFile(TEST_FOLDER + "/data/MeshWithHexagon.nc");
    meshkernelapi::Mesh2D mesh2d;
    mesh2d.num_edges = static_cast<int>(num_edges);
    mesh2d.num_nodes = static_cast<int>(num_nodes);
    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_nodes = edge_nodes.data();

    auto errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernel::OrthogonalizationParameters orthogonalizationParameters{};
    orthogonalizationParameters.outer_iterations = 1;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;

    meshkernelapi::GeometryList geometryList{};
    meshkernelapi::GeometryList landBoundaries{};

    // Execute
    errorCode = mkernel_mesh2d_initialize_orthogonalization(meshKernelId,
                                                            1,
                                                            orthogonalizationParameters,
                                                            landBoundaries,
                                                            geometryList);

    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_prepare_outer_iteration_orthogonalization(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_compute_inner_ortogonalization_iteration(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_finalize_inner_ortogonalization_iteration(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_delete_orthogonalization(meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);
}

TEST_F(CartesianApiTests, ContactsComputeSingle_OnMesh2D_ShouldComputeContacts)
{
    auto [nodes_x, nodes_y, edges, face_nodes, num_face_nodes] = MakeMeshWithFaceNodes();
    const auto meshKernelId = GetMeshKernelId();

    meshkernelapi::Mesh2D mesh2d;
    mesh2d.node_x = &nodes_x[0];
    mesh2d.node_y = &nodes_y[0];

    mesh2d.edge_nodes = &edges[0];
    mesh2d.face_nodes = &face_nodes[0];
    mesh2d.nodes_per_face = &num_face_nodes[0];

    mesh2d.num_nodes = static_cast<int>(nodes_x.size());
    mesh2d.num_edges = static_cast<int>(edges.size() / 2);
    mesh2d.num_faces = static_cast<int>(num_face_nodes.size());

    auto errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Init 1d mesh
    meshkernelapi::Mesh1D mesh1d;
    std::vector<double> node_x{
        1.73493900000000,
        2.35659313023165,
        5.38347452702839,
        14.2980910429074,
        22.9324017677239,
        25.3723169493137,
        25.8072280000000};
    std::vector<double> node_y{
        -7.6626510000000,
        1.67281447902331,
        10.3513746546384,
        12.4797224193970,
        15.3007317677239,
        24.1623588554512,
        33.5111870000000};
    mesh1d.node_x = node_x.data();
    mesh1d.node_y = node_y.data();
    mesh1d.num_nodes = 7;

    std::vector<int> edge_nodes{0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6};
    mesh1d.edge_nodes = edge_nodes.data();
    mesh1d.num_edges = 6;

    errorCode = mkernel_mesh1d_set(meshKernelId, mesh1d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Init 1d mask
    std::vector<int> onedNodeMask{1, 1, 1, 1, 1, 1, 1};

    // Init polygon
    meshkernelapi::GeometryList polygon;
    polygon.geometry_separator = meshkernel::constants::missing::doubleValue;

    std::vector<double> xCoordinates{-30, 40, 40, -40, -30};
    std::vector<double> yCoordinates{-20, -20, 50, 50, -20};
    std::vector<double> zCoordinates{0, 0, 0, 0, 0};
    polygon.coordinates_x = xCoordinates.data();
    polygon.coordinates_y = yCoordinates.data();
    polygon.values = zCoordinates.data();
    polygon.num_coordinates = static_cast<int>(xCoordinates.size());

    // Execute
    errorCode = mkernel_contacts_compute_single(meshKernelId, onedNodeMask.data(), polygon, 0.0);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Get the new state
    meshkernelapi::Contacts contacts{};
    errorCode = mkernel_contacts_get_dimensions(meshKernelId, contacts);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    std::vector<int> mesh1d_indices(contacts.num_contacts);
    std::vector<int> mesh2d_indices(contacts.num_contacts);
    contacts.mesh1d_indices = mesh1d_indices.data();
    contacts.mesh2d_indices = mesh2d_indices.data();

    errorCode = mkernel_contacts_get_data(meshKernelId, contacts);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(5, contacts.num_contacts);

    ASSERT_EQ(1, contacts.mesh1d_indices[0]);
    ASSERT_EQ(2, contacts.mesh1d_indices[1]);
    ASSERT_EQ(3, contacts.mesh1d_indices[2]);
    ASSERT_EQ(4, contacts.mesh1d_indices[3]);
    ASSERT_EQ(5, contacts.mesh1d_indices[4]);

    ASSERT_EQ(0, contacts.mesh2d_indices[0]);
    ASSERT_EQ(3, contacts.mesh2d_indices[1]);
    ASSERT_EQ(4, contacts.mesh2d_indices[2]);
    ASSERT_EQ(5, contacts.mesh2d_indices[3]);
    ASSERT_EQ(8, contacts.mesh2d_indices[4]);
}

TEST(Mesh2D, IntersectionsFromPolyline_ShouldIntersectMesh)
{
    // Setup
    int meshKernelId;
    auto errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);

    // Create a curvilinear grid in the back-end and convert to an unstructured grid
    meshkernel::MakeGridParameters makeMeshParameters;
    makeMeshParameters.num_columns = 3;
    makeMeshParameters.num_rows = 3;
    makeMeshParameters.block_size_x = 1.0;
    makeMeshParameters.block_size_y = 1.0;
    makeMeshParameters.origin_x = 0.0;
    makeMeshParameters.origin_y = 0.0;
    makeMeshParameters.angle = 0.0;

    // Empty polygon
    meshkernelapi::GeometryList geometryList;
    geometryList.num_coordinates = 0;

    // Creates an unstructured grid from mesh parameters
    errorCode = mkernel_mesh2d_make_uniform(meshKernelId, makeMeshParameters, geometryList);

    // Get the mesh dimensions
    meshkernelapi::Mesh2D mesh2dDimensions{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2dDimensions);

    // Set the polyLine
    std::vector xCoordinates{0.6, 0.6};
    std::vector yCoordinates{2.5, 0.5};

    meshkernelapi::GeometryList boundaryPolyLine{};
    boundaryPolyLine.geometry_separator = meshkernel::constants::missing::doubleValue;
    boundaryPolyLine.coordinates_x = xCoordinates.data();
    boundaryPolyLine.coordinates_y = yCoordinates.data();
    boundaryPolyLine.values = nullptr;

    boundaryPolyLine.num_coordinates = 2;

    std::vector<int> polylineSegmentIndexes(mesh2dDimensions.num_edges * 2, meshkernel::constants::missing::intValue);

    std::vector<int> edgeNodes(mesh2dDimensions.num_edges * 2, meshkernel::constants::missing::intValue);
    std::vector<int> edgeIndex(mesh2dDimensions.num_edges, meshkernel::constants::missing::intValue);
    std::vector<double> edgeDistances(mesh2dDimensions.num_edges * 2, meshkernel::constants::missing::doubleValue);
    std::vector<double> segmentDistances(mesh2dDimensions.num_edges * 2, meshkernel::constants::missing::doubleValue);
    std::vector<int> segmentIndexes(mesh2dDimensions.num_edges, meshkernel::constants::missing::intValue);

    std::vector<int> faceEdgeIndex(mesh2dDimensions.num_edges, meshkernel::constants::missing::intValue);
    std::vector<int> faceNumEdges(mesh2dDimensions.num_edges, meshkernel::constants::missing::intValue);
    std::vector<int> faceIndexes(mesh2dDimensions.num_edges, meshkernel::constants::missing::intValue);

    errorCode = mkernel_mesh2d_intersections_from_polyline(meshKernelId,
                                                           boundaryPolyLine,
                                                           edgeNodes.data(),
                                                           edgeIndex.data(),
                                                           edgeDistances.data(),
                                                           segmentDistances.data(),
                                                           segmentIndexes.data(),
                                                           faceIndexes.data(),
                                                           faceNumEdges.data(),
                                                           faceEdgeIndex.data());

    /// Assert
    const double tolerance = 1e-6;

    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    ASSERT_EQ(segmentIndexes[0], 0);
    ASSERT_EQ(segmentIndexes[1], 0);
    ASSERT_EQ(segmentIndexes[2], meshkernel::constants::missing::intValue);

    ASSERT_NEAR(segmentDistances[0], 0.25, tolerance);
    ASSERT_NEAR(segmentDistances[1], 0.75, tolerance);
    ASSERT_NEAR(segmentDistances[2], meshkernel::constants::missing::doubleValue, tolerance);

    ASSERT_NEAR(edgeDistances[0], 0.6, tolerance);
    ASSERT_NEAR(edgeDistances[1], 0.6, tolerance);
    ASSERT_NEAR(edgeDistances[2], meshkernel::constants::missing::doubleValue, tolerance);

    ASSERT_EQ(faceIndexes[0], 6);
    ASSERT_EQ(faceIndexes[1], 3);
    ASSERT_EQ(faceIndexes[2], 3);
    ASSERT_EQ(faceIndexes[3], 0);
    ASSERT_EQ(faceIndexes[4], meshkernel::constants::missing::intValue);

    ASSERT_EQ(faceNumEdges[0], 1);
    ASSERT_EQ(faceNumEdges[1], 2);
    ASSERT_EQ(faceNumEdges[2], 1);
    ASSERT_EQ(faceNumEdges[3], meshkernel::constants::missing::intValue);

    ASSERT_EQ(faceEdgeIndex[0], 18);
    ASSERT_EQ(faceEdgeIndex[1], 18);
    ASSERT_EQ(faceEdgeIndex[2], 15);
    ASSERT_EQ(faceEdgeIndex[3], 15);
    ASSERT_EQ(faceEdgeIndex[4], meshkernel::constants::missing::intValue);
}
TEST(Mesh2D, CurvilinearMakeUniformOnExtension_OnSpericalCoordinates_ShouldGenerateCurvilinearMesh)
{
    // Setup
    auto makeGridParameters = meshkernel::MakeGridParameters();
    makeGridParameters.origin_x = -1.0;
    makeGridParameters.origin_y = 49.1;
    makeGridParameters.block_size_x = 0.01;
    makeGridParameters.block_size_y = 0.01;
    makeGridParameters.upper_right_x = -0.2;
    makeGridParameters.upper_right_y = 49.6;

    // Execute
    int meshKernelId;
    const int projectionType = 1;
    auto errorCode = meshkernelapi::mkernel_allocate_state(projectionType, meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_make_uniform_on_extension(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d;
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(mesh2d.num_nodes, 8343);
    ASSERT_EQ(mesh2d.num_edges, 16502);
    ASSERT_EQ(mesh2d.num_faces, 8160);
}

TEST(Mesh2D, Mesh2DRefineBasedOnGriddedSamples_WithGriddedSamples_ShouldRefineMesh)
{
    // Prepare
    int meshKernelId;
    constexpr int isGeographic = 0;
    meshkernelapi::mkernel_allocate_state(isGeographic, meshKernelId);

    auto [num_nodes, num_edges, node_x, node_y, node_type, edge_nodes, edge_type] = ReadLegacyMeshFile(TEST_FOLDER + "/data/MeshRefinementTests/gebco.nc");
    meshkernelapi::Mesh2D mesh2d;
    mesh2d.num_edges = static_cast<int>(num_edges);
    mesh2d.num_nodes = static_cast<int>(num_nodes);
    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_nodes = edge_nodes.data();

    auto errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    auto [numX, numY, xllCenter, yllCenter, cellSize, nodatavalue, values] = ReadAscFile(TEST_FOLDER + "/data/MeshRefinementTests/gebco.asc");
    meshkernelapi::GriddedSamples griddedSamples;
    griddedSamples.num_x = numX;
    griddedSamples.num_y = numY;
    griddedSamples.x_origin = xllCenter;
    griddedSamples.y_origin = yllCenter;
    griddedSamples.cell_size = cellSize;
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
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::Mesh2D mesh2dResults;
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2dResults);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    ASSERT_EQ(417, mesh2dResults.num_nodes);
    ASSERT_EQ(988, mesh2dResults.num_edges);
    ASSERT_EQ(572, mesh2dResults.num_faces);
    ASSERT_EQ(1936, mesh2dResults.num_face_nodes);
}

TEST_F(CartesianApiTests, Mesh2DRefineBasedOnGriddedSamples_WithNotUniformlySpacedSamples_ShouldRefineMesh)
{
    // Prepare
    meshkernel::UInt nRows{5};
    meshkernel::UInt nCols{4};
    MakeMesh(nRows, nCols, 100.0);
    auto const meshKernelId = GetMeshKernelId();

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

    auto errorCode = mkernel_mesh2d_refine_based_on_gridded_samples(meshKernelId, griddedSamples, meshRefinementParameters, true);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::Mesh2D mesh2dResults;
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2dResults);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    ASSERT_EQ(86, mesh2dResults.num_nodes);
    ASSERT_EQ(161, mesh2dResults.num_edges);
    ASSERT_EQ(76, mesh2dResults.num_faces);
}

TEST(Mesh2D, RefineBasedOnGriddedSamples_WithUniformSamplesAndSphericalCoordinates_ShouldRefineMesh2d)
{
    // Prepare
    int meshKernelId;
    constexpr int isSpherical = 1;
    meshkernelapi::mkernel_allocate_state(isSpherical, meshKernelId);

    const auto makeGridParameters = GebcoMakeGridParameters();

    meshkernelapi::GeometryList geometryList{};
    auto errorCode = meshkernelapi::mkernel_mesh2d_make_uniform(meshKernelId, makeGridParameters, geometryList);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    auto [ncols, nrows, xllcenter, yllcenter, cellsize, nodata_value, values] = ReadAscFile(TEST_FOLDER + "/data/MeshRefinementTests/gebco.asc");
    meshkernelapi::GriddedSamples griddedSamples;
    griddedSamples.num_x = ncols;
    griddedSamples.num_y = nrows;
    griddedSamples.x_origin = xllcenter;
    griddedSamples.y_origin = yllcenter;
    griddedSamples.cell_size = cellsize;
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
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::Mesh2D mesh2dResults;
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2dResults);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    ASSERT_EQ(5223, mesh2dResults.num_nodes);
    ASSERT_EQ(10745, mesh2dResults.num_edges);
    ASSERT_EQ(5523, mesh2dResults.num_faces);
    ASSERT_EQ(21212, mesh2dResults.num_face_nodes);
}

TEST(Mesh2D, RefineBasedOnGriddedSamples_WithUniformSamplesAndSphericalCoordinatesAndLargeMinEdgeSize_ShouldNotRefineMesh2d)
{
    // Prepare
    int meshKernelId;
    constexpr int isSpherical = 1;
    meshkernelapi::mkernel_allocate_state(isSpherical, meshKernelId);

    const auto makeGridParameters = GebcoMakeGridParameters();

    meshkernelapi::GeometryList geometryList{};
    auto errorCode = meshkernelapi::mkernel_mesh2d_make_uniform(meshKernelId, makeGridParameters, geometryList);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    auto [ncols, nrows, xllcenter, yllcenter, cellsize, nodata_value, values] = ReadAscFile(TEST_FOLDER + "/data/MeshRefinementTests/gebco.asc");
    meshkernelapi::GriddedSamples griddedSamples;
    griddedSamples.num_x = ncols;
    griddedSamples.num_y = nrows;
    griddedSamples.x_origin = xllcenter;
    griddedSamples.y_origin = yllcenter;
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
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::Mesh2D mesh2dResults;
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2dResults);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert no refinement took place
    ASSERT_EQ(54, (makeGridParameters.num_columns + 1) * (makeGridParameters.num_rows + 1));
}

TEST(CurvilinearGrid, MakeUniform_OnSphericalCoordinates_ShouldMakeCurvilinearGrid)
{
    meshkernel::MakeGridParameters makeGridParameters;
    meshkernelapi::GeometryList geometryList;

    makeGridParameters.origin_x = -1.0;
    makeGridParameters.origin_y = 49.1;
    makeGridParameters.block_size_x = 0.1;
    makeGridParameters.block_size_y = 0.1;
    makeGridParameters.num_columns = 8;
    makeGridParameters.num_rows = 11;

    int meshKernelId = 0;
    const int projectionType = 1;
    auto errorCode = meshkernelapi::mkernel_allocate_state(projectionType, meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    errorCode = mkernel_curvilinear_make_uniform(meshKernelId, makeGridParameters, geometryList);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGridResults;
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGridResults);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    ASSERT_EQ(12, curvilinearGridResults.num_m);
    ASSERT_EQ(9, curvilinearGridResults.num_n);
}
TEST(CurvilinearGrid, MakeUniformOnDefinedExtension_OnSphericalCoordinates_ShouldMakeCurvilinearGrid)
{
    // Setup
    meshkernel::MakeGridParameters makeGridParameters;
    makeGridParameters.origin_x = 0.0;
    makeGridParameters.origin_y = 0.0;
    makeGridParameters.upper_right_x = 10.0;
    makeGridParameters.upper_right_y = 10.0;
    makeGridParameters.block_size_x = 1.0;
    makeGridParameters.block_size_y = 2.0;

    int meshKernelId = 0;
    int projectionType = 0;
    auto errorCode = meshkernelapi::mkernel_allocate_state(projectionType, meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Execute
    errorCode = meshkernelapi::mkernel_curvilinear_make_uniform_on_extension(meshKernelId, makeGridParameters);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    meshkernelapi::CurvilinearGrid curvilinearGridResults;
    errorCode = mkernel_curvilinear_get_dimensions(meshKernelId, curvilinearGridResults);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(6, curvilinearGridResults.num_m);
    ASSERT_EQ(11, curvilinearGridResults.num_n);
}

TEST(MeshState, MKernelGetProjection_ShouldGetProjection)
{
    // Setup
    int meshKernelId = 0;
    int setProjectionType = 0;
    auto errorCode = meshkernelapi::mkernel_allocate_state(setProjectionType, meshKernelId);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    int getProjectionType = 0;
    errorCode = meshkernelapi::mkernel_get_projection(meshKernelId, getProjectionType);
    ASSERT_EQ(meshkernelapi::MeshKernelApiErrors::Success, errorCode);

    // Assert
    ASSERT_EQ(setProjectionType, getProjectionType);
}
