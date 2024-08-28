#include <gtest/gtest.h>

#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/FlipEdges.hpp>
#include <MeshKernel/LandBoundaries.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/MeshRefinement.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/OrthogonalizationAndSmoothing.hpp>
#include <MeshKernel/Orthogonalizer.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/Smoother.hpp>
#include <MeshKernel/UndoActions/UndoAction.hpp>
#include <TestUtils/Definitions.hpp>
#include <TestUtils/MakeMeshes.hpp>

using namespace meshkernel;

TEST(OrthogonalizationAndSmoothing, OrthogonalizationOneQuadOneTriangle)
{
    // Preparation
    std::vector<Point> nodes{{0.0, 0.0},
                             {0.0, 10.0},
                             {10.0, 0.0},
                             {10.0, 10.0},
                             {20.0, 0.0}};

    std::vector<Edge> edges{{1, 0},
                            {0, 2},
                            {2, 4},
                            {4, 3},
                            {3, 2},
                            {3, 1}};

    const auto projectToLandBoundaryOption = LandBoundaries::ProjectToLandBoundaryOption::DoNotProjectToLandBoundary;
    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.inner_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.outer_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor_at_boundary = 0.975;
    orthogonalizationParameters.areal_to_angle_smoothing_factor = 1.0;

    // Execute
    auto mesh = Mesh2D(edges, nodes, Projection::cartesian);
    auto orthogonalizer = std::make_unique<Orthogonalizer>(mesh);
    auto smoother = std::make_unique<Smoother>(mesh);
    auto polygon = std::make_unique<Polygons>();

    std::vector<Point> landBoundary{};
    auto landboundaries = std::make_unique<LandBoundaries>(landBoundary, mesh, *polygon);

    OrthogonalizationAndSmoothing orthogonalization(mesh,
                                                    std::move(smoother),
                                                    std::move(orthogonalizer),
                                                    std::move(polygon),
                                                    std::move(landboundaries),
                                                    projectToLandBoundaryOption,
                                                    orthogonalizationParameters);

    [[maybe_unused]] auto undoAction = orthogonalization.Initialize();
    orthogonalization.Compute();

    // Assert
    constexpr double tolerance = 1e-8;
    ASSERT_NEAR(0.0, mesh.Node(0).x, tolerance);
    ASSERT_NEAR(0.0, mesh.Node(1).x, tolerance);
    ASSERT_NEAR(10.0, mesh.Node(2).x, tolerance);
    ASSERT_NEAR(10.0, mesh.Node(3).x, tolerance);
    ASSERT_NEAR(20.0, mesh.Node(4).x, tolerance);

    ASSERT_NEAR(0.0, mesh.Node(0).y, tolerance);
    ASSERT_NEAR(10.0, mesh.Node(1).y, tolerance);
    ASSERT_NEAR(0.0, mesh.Node(2).y, tolerance);
    ASSERT_NEAR(10.0, mesh.Node(3).y, tolerance);
    ASSERT_NEAR(0.0, mesh.Node(4).y, tolerance);
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizationSmallTriangularGrid)
{

    // now build node-edge mapping
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/SmallTriangularGrid_net.nc");
    const auto projectToLandBoundaryOption = LandBoundaries::ProjectToLandBoundaryOption::DoNotProjectToLandBoundary;
    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor_at_boundary = 1.0;
    orthogonalizationParameters.areal_to_angle_smoothing_factor = 1.0;

    auto orthogonalizer = std::make_unique<Orthogonalizer>(*mesh);
    auto smoother = std::make_unique<Smoother>(*mesh);
    auto polygon = std::make_unique<Polygons>();

    std::vector<Point> landBoundary;
    auto landboundaries = std::make_unique<LandBoundaries>(landBoundary, *mesh, *polygon);

    OrthogonalizationAndSmoothing orthogonalization(*mesh,
                                                    std::move(smoother),
                                                    std::move(orthogonalizer),
                                                    std::move(polygon),
                                                    std::move(landboundaries),
                                                    projectToLandBoundaryOption,
                                                    orthogonalizationParameters);

    [[maybe_unused]] auto undoAction = orthogonalization.Initialize();
    orthogonalization.Compute();

    constexpr double tolerance = 0.8e-2;

    ASSERT_NEAR(325.590101919525, mesh->Node(0).x, tolerance);
    ASSERT_NEAR(229.213730481198, mesh->Node(1).x, tolerance);
    ASSERT_NEAR(263.439319753147, mesh->Node(2).x, tolerance);
    ASSERT_NEAR(429.191105834504, mesh->Node(3).x, tolerance);
    ASSERT_NEAR(535.865215426468, mesh->Node(4).x, tolerance);
    ASSERT_NEAR(503.753784179688, mesh->Node(5).x, tolerance);
    ASSERT_NEAR(354.048340705929, mesh->Node(6).x, tolerance);
    ASSERT_NEAR(346.790050854504, mesh->Node(7).x, tolerance);
    ASSERT_NEAR(315.030130405285, mesh->Node(8).x, tolerance);
    ASSERT_NEAR(424.314957449766, mesh->Node(9).x, tolerance);

    ASSERT_NEAR(455.319334078551, mesh->Node(0).y, tolerance);
    ASSERT_NEAR(362.573521507281, mesh->Node(1).y, tolerance);
    ASSERT_NEAR(241.096458631763, mesh->Node(2).y, tolerance);
    ASSERT_NEAR(211.483073921775, mesh->Node(3).y, tolerance);
    ASSERT_NEAR(311.401495506714, mesh->Node(4).y, tolerance);
    ASSERT_NEAR(432.379974365234, mesh->Node(5).y, tolerance);
    ASSERT_NEAR(458.064836627594, mesh->Node(6).y, tolerance);
    ASSERT_NEAR(405.311585650679, mesh->Node(7).y, tolerance);
    ASSERT_NEAR(319.612138503550, mesh->Node(8).y, tolerance);
    ASSERT_NEAR(327.102805172725, mesh->Node(9).y, tolerance);
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizationSmallTriangularGridAsNcFile)
{

    // now build node-edge mapping
    auto mesh = MakeSmallSizeTriangularMeshForTestingAsNcFile();

    const auto projectToLandBoundaryOption = LandBoundaries::ProjectToLandBoundaryOption::DoNotProjectToLandBoundary;
    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor_at_boundary = 1.0;
    orthogonalizationParameters.areal_to_angle_smoothing_factor = 1.0;

    auto orthogonalizer = std::make_unique<Orthogonalizer>(*mesh);
    auto smoother = std::make_unique<Smoother>(*mesh);
    auto polygon = std::make_unique<Polygons>();

    std::vector<Point> landBoundary;
    auto landboundaries = std::make_unique<LandBoundaries>(landBoundary, *mesh, *polygon);

    OrthogonalizationAndSmoothing orthogonalization(*mesh,
                                                    std::move(smoother),
                                                    std::move(orthogonalizer),
                                                    std::move(polygon),
                                                    std::move(landboundaries),
                                                    projectToLandBoundaryOption,
                                                    orthogonalizationParameters);

    [[maybe_unused]] auto undoAction = orthogonalization.Initialize();
    orthogonalization.Compute();

    constexpr double tolerance = 1e-2;

    ASSERT_NEAR(325.590101919525, mesh->Node(0).x, tolerance);
    ASSERT_NEAR(229.213730481198, mesh->Node(1).x, tolerance);
    ASSERT_NEAR(263.439319753147, mesh->Node(2).x, tolerance);
    ASSERT_NEAR(429.191105834504, mesh->Node(3).x, tolerance);
    ASSERT_NEAR(535.865215426468, mesh->Node(4).x, tolerance);
    ASSERT_NEAR(503.753784179688, mesh->Node(5).x, tolerance);
    ASSERT_NEAR(354.048340705929, mesh->Node(6).x, tolerance);
    ASSERT_NEAR(346.790050854504, mesh->Node(7).x, tolerance);
    ASSERT_NEAR(315.030130405285, mesh->Node(8).x, tolerance);
    ASSERT_NEAR(424.314957449766, mesh->Node(9).x, tolerance);

    ASSERT_NEAR(455.319334078551, mesh->Node(0).y, tolerance);
    ASSERT_NEAR(362.573521507281, mesh->Node(1).y, tolerance);
    ASSERT_NEAR(241.096458631763, mesh->Node(2).y, tolerance);
    ASSERT_NEAR(211.483073921775, mesh->Node(3).y, tolerance);
    ASSERT_NEAR(311.401495506714, mesh->Node(4).y, tolerance);
    ASSERT_NEAR(432.379974365234, mesh->Node(5).y, tolerance);
    ASSERT_NEAR(458.064836627594, mesh->Node(6).y, tolerance);
    ASSERT_NEAR(405.311585650679, mesh->Node(7).y, tolerance);
    ASSERT_NEAR(319.612138503550, mesh->Node(8).y, tolerance);
    ASSERT_NEAR(327.102805172725, mesh->Node(9).y, tolerance);
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizationMediumTriangularGridWithPolygon)
{
    // now build node-edge mapping
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/SmallTriangularGrid_net.nc");

    // The original mesh nodes
    const std::vector<meshkernel::Point> meshNodes(mesh->Nodes());

    const auto projectToLandBoundaryOption = LandBoundaries::ProjectToLandBoundaryOption::DoNotProjectToLandBoundary;
    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor_at_boundary = 1.0;
    orthogonalizationParameters.areal_to_angle_smoothing_factor = 1.0;

    std::vector<Point> nodes{{342.987518, 471.121002},
                             {327.640900, 380.846436},
                             {396.851135, 201.200073},
                             {514.207581, 203.607407},
                             {569.274841, 294.483765},
                             {568.673035, 379.943695},
                             {515.712158, 458.783478},
                             {343.288422, 471.722809},
                             {342.987518, 471.121002}};

    auto orthogonalizer = std::make_unique<Orthogonalizer>(*mesh);
    auto smoother = std::make_unique<Smoother>(*mesh);

    std::vector<Point> landBoundary;
    auto polygon = std::make_unique<Polygons>(nodes, Projection::cartesian);
    auto landboundaries = std::make_unique<LandBoundaries>(landBoundary, *mesh, *polygon);

    OrthogonalizationAndSmoothing orthogonalization(*mesh,
                                                    std::move(smoother),
                                                    std::move(orthogonalizer),
                                                    std::move(polygon),
                                                    std::move(landboundaries),
                                                    projectToLandBoundaryOption,
                                                    orthogonalizationParameters);

    auto undoAction = orthogonalization.Initialize();
    orthogonalization.Compute();

    constexpr double tolerance = 1.8;
    // check the first 10 points

    ASSERT_NEAR(322.252624511719, mesh->Node(0).x, tolerance);
    ASSERT_NEAR(227.002044677734, mesh->Node(1).x, tolerance);
    ASSERT_NEAR(259.252227783203, mesh->Node(2).x, tolerance);
    ASSERT_NEAR(426.680730020925, mesh->Node(3).x, tolerance);
    ASSERT_NEAR(535.742871000462, mesh->Node(4).x, tolerance);
    ASSERT_NEAR(502.614684581472, mesh->Node(5).x, tolerance);
    ASSERT_NEAR(350.424583126127, mesh->Node(6).x, tolerance);
    ASSERT_NEAR(342.832067902036, mesh->Node(7).x, tolerance);
    ASSERT_NEAR(310.300984548069, mesh->Node(8).x, tolerance);
    ASSERT_NEAR(421.791577149564, mesh->Node(9).x, tolerance);

    ASSERT_NEAR(454.880187988281, mesh->Node(0).y, tolerance);
    ASSERT_NEAR(360.379241943359, mesh->Node(1).y, tolerance);
    ASSERT_NEAR(241.878051757812, mesh->Node(2).y, tolerance);
    ASSERT_NEAR(210.624626374983, mesh->Node(3).y, tolerance);
    ASSERT_NEAR(311.862423032534, mesh->Node(4).y, tolerance);
    ASSERT_NEAR(432.575408917275, mesh->Node(5).y, tolerance);
    ASSERT_NEAR(458.587061164954, mesh->Node(6).y, tolerance);
    ASSERT_NEAR(405.390188082498, mesh->Node(7).y, tolerance);
    ASSERT_NEAR(319.410057398020, mesh->Node(8).y, tolerance);
    ASSERT_NEAR(327.001109057344, mesh->Node(9).y, tolerance);

    undoAction->Restore();

    constexpr double nodeTolerance = 1.0e-10;

    for (meshkernel::UInt i = 0; i < meshNodes.size(); ++i)
    {
        EXPECT_NEAR(meshNodes[i].x, mesh->Node(i).x, nodeTolerance);
        EXPECT_NEAR(meshNodes[i].y, mesh->Node(i).y, nodeTolerance);
    }
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizationMediumTriangularGridWithUndo)
{
    // now build node-edge mapping
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/TestOrthogonalizationMediumTriangularGrid_net.nc");

    const auto projectToLandBoundaryOption = LandBoundaries::ProjectToLandBoundaryOption::DoNotProjectToLandBoundary;
    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor_at_boundary = 0.5;
    orthogonalizationParameters.areal_to_angle_smoothing_factor = 1.0;

    auto orthogonalizer = std::make_unique<Orthogonalizer>(*mesh);
    auto smoother = std::make_unique<Smoother>(*mesh);
    auto polygon = std::make_unique<Polygons>();
    std::vector<Point> landBoundary;
    auto landBoundaries = std::make_unique<LandBoundaries>(landBoundary, *mesh, *polygon);

    const std::vector<meshkernel::Point> meshNodes = mesh->Nodes();

    OrthogonalizationAndSmoothing orthogonalization(*mesh,
                                                    std::move(smoother),
                                                    std::move(orthogonalizer),
                                                    std::move(polygon),
                                                    std::move(landBoundaries),
                                                    projectToLandBoundaryOption,
                                                    orthogonalizationParameters);

    std::unique_ptr<meshkernel::UndoAction> undoAction = orthogonalization.Initialize();
    orthogonalization.Compute();

    constexpr double tolerance = 1.2;

    // check the first 10 points
    ASSERT_NEAR(68.771705432835475, mesh->Node(0).x, tolerance);
    ASSERT_NEAR(169.49338272334273, mesh->Node(1).x, tolerance);
    ASSERT_NEAR(262.80128484924921, mesh->Node(2).x, tolerance);
    ASSERT_NEAR(361.60010033352023, mesh->Node(3).x, tolerance);
    ASSERT_NEAR(468.13991812406925, mesh->Node(4).x, tolerance);
    ASSERT_NEAR(549.89461192844624, mesh->Node(5).x, tolerance);
    ASSERT_NEAR(653.02704974527421, mesh->Node(6).x, tolerance);
    ASSERT_NEAR(747.81537706979441, mesh->Node(7).x, tolerance);
    ASSERT_NEAR(853.40641427112951, mesh->Node(8).x, tolerance);
    ASSERT_NEAR(938.69752431820143, mesh->Node(9).x, tolerance);

    ASSERT_NEAR(1399.7751472360221, mesh->Node(0).y, tolerance);
    ASSERT_NEAR(1426.5945287630802, mesh->Node(1).y, tolerance);
    ASSERT_NEAR(1451.4398281457179, mesh->Node(2).y, tolerance);
    ASSERT_NEAR(1477.7472050498141, mesh->Node(3).y, tolerance);
    ASSERT_NEAR(1506.1157955857589, mesh->Node(4).y, tolerance);
    ASSERT_NEAR(1527.8847968946166, mesh->Node(5).y, tolerance);
    ASSERT_NEAR(1555.3460969050145, mesh->Node(6).y, tolerance);
    ASSERT_NEAR(1580.5855923464549, mesh->Node(7).y, tolerance);
    ASSERT_NEAR(1608.7015489976982, mesh->Node(8).y, tolerance);
    ASSERT_NEAR(1631.412199601948, mesh->Node(9).y, tolerance);

    undoAction->Restore();

    EXPECT_EQ(meshNodes.size(), mesh->GetNumNodes());

    for (meshkernel::UInt i = 0; i < meshNodes.size(); ++i)
    {
        EXPECT_EQ(meshNodes[i].x, mesh->Node(i).x);
        EXPECT_EQ(meshNodes[i].y, mesh->Node(i).y);
    }
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizationFourQuads)
{
    auto mesh = MakeRectangularMeshForTesting(3, 3, 1.0, Projection::cartesian);

    const auto projectToLandBoundaryOption = LandBoundaries::ProjectToLandBoundaryOption::DoNotProjectToLandBoundary;
    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.inner_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.outer_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor_at_boundary = 0.975;
    orthogonalizationParameters.areal_to_angle_smoothing_factor = 1.0;

    auto polygon = std::make_unique<Polygons>();
    std::vector<Point> landBoundary{};
    auto landboundaries = std::make_unique<LandBoundaries>(landBoundary, *mesh, *polygon);

    auto orthogonalizer = std::make_unique<Orthogonalizer>(*mesh);
    auto smoother = std::make_unique<Smoother>(*mesh);
    OrthogonalizationAndSmoothing orthogonalization(*mesh,
                                                    std::move(smoother),
                                                    std::move(orthogonalizer),
                                                    std::move(polygon),
                                                    std::move(landboundaries),
                                                    projectToLandBoundaryOption,
                                                    orthogonalizationParameters);

    [[maybe_unused]] auto undoAction = orthogonalization.Initialize();
    orthogonalization.Compute();
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizeAndSnapToLandBoundaries)
{
    using namespace constants;

    // Prepare
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/SmallTriangularGrid_net.nc");

    // the land boundary to use
    std::vector<Point> landBoundary{{235.561218, 290.571899},
                                    {265.953522, 436.515747},
                                    {429.349854, 450.959656},
                                    {535.271545, 386.262909},
                                    {missing::doubleValue, missing::doubleValue},
                                    {246.995941, 262.285858},
                                    {351.112183, 237.309906},
                                    {443.191895, 262.285858},
                                    {553.627319, 327.283539},
                                    {missing::doubleValue, missing::doubleValue}};

    // snap to land boundaries
    const auto projectToLandBoundaryOption = LandBoundaries::ProjectToLandBoundaryOption::OuterMeshBoundaryToLandBoundary;
    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor_at_boundary = 0.975;
    orthogonalizationParameters.areal_to_angle_smoothing_factor = 1.0;

    // no enclosing polygon
    auto polygon = std::make_unique<Polygons>();
    auto landboundaries = std::make_unique<LandBoundaries>(landBoundary, *mesh, *polygon);

    auto orthogonalizer = std::make_unique<Orthogonalizer>(*mesh);
    auto smoother = std::make_unique<Smoother>(*mesh);
    OrthogonalizationAndSmoothing orthogonalization(*mesh,
                                                    std::move(smoother),
                                                    std::move(orthogonalizer),
                                                    std::move(polygon),
                                                    std::move(landboundaries),
                                                    projectToLandBoundaryOption,
                                                    orthogonalizationParameters);

    [[maybe_unused]] auto undoAction = orthogonalization.Initialize();
    orthogonalization.Compute();

    // check the values
    constexpr double tolerance = 0.15;
    ASSERT_NEAR(313.081472564480, mesh->Node(0).x, tolerance);
    ASSERT_NEAR(253.641466857330, mesh->Node(1).x, tolerance);
    ASSERT_NEAR(254.777224294204, mesh->Node(2).x, tolerance);
    ASSERT_NEAR(443.191895000000, mesh->Node(3).x, tolerance);
    ASSERT_NEAR(535.240231516760, mesh->Node(4).x, tolerance);
    ASSERT_NEAR(480.436129612752, mesh->Node(5).x, tolerance);
    ASSERT_NEAR(345.948240805397, mesh->Node(6).x, tolerance);
    ASSERT_NEAR(342.668434889472, mesh->Node(7).x, tolerance);
    ASSERT_NEAR(318.414413615199, mesh->Node(8).x, tolerance);
    ASSERT_NEAR(424.616311031376, mesh->Node(9).x, tolerance);

    ASSERT_NEAR(440.681763586650, mesh->Node(0).y, tolerance);
    ASSERT_NEAR(377.393256506700, mesh->Node(1).y, tolerance);
    ASSERT_NEAR(260.419242817573, mesh->Node(2).y, tolerance);
    ASSERT_NEAR(262.285858000000, mesh->Node(3).y, tolerance);
    ASSERT_NEAR(316.461666783032, mesh->Node(4).y, tolerance);
    ASSERT_NEAR(419.756265860671, mesh->Node(5).y, tolerance);
    ASSERT_NEAR(443.587120174434, mesh->Node(6).y, tolerance);
    ASSERT_NEAR(402.913858250569, mesh->Node(7).y, tolerance);
    ASSERT_NEAR(336.831643075189, mesh->Node(8).y, tolerance);
    ASSERT_NEAR(340.875100904741, mesh->Node(9).y, tolerance);
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizationSphericalRectangular)
{
    // 1 Setup
    auto mesh = MakeRectangularMeshForTesting(4, 4, 0.003, Projection::spherical, {41.1, 41.1});

    const auto projectToLandBoundaryOption = LandBoundaries::ProjectToLandBoundaryOption::DoNotProjectToLandBoundary;
    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor_at_boundary = 1.0;
    orthogonalizationParameters.areal_to_angle_smoothing_factor = 1.0;

    auto orthogonalizer = std::make_unique<Orthogonalizer>(*mesh);
    auto smoother = std::make_unique<Smoother>(*mesh);

    auto polygon = std::make_unique<Polygons>();
    std::vector<Point> landBoundary{};
    auto landboundaries = std::make_unique<LandBoundaries>(landBoundary, *mesh, *polygon);

    OrthogonalizationAndSmoothing orthogonalization(*mesh,
                                                    std::move(smoother),
                                                    std::move(orthogonalizer),
                                                    std::move(polygon),
                                                    std::move(landboundaries),
                                                    projectToLandBoundaryOption,
                                                    orthogonalizationParameters);

    [[maybe_unused]] auto undoAction = orthogonalization.Initialize();
    orthogonalization.Compute();

    // check the values
    constexpr double tolerance = 1e-6;
    ASSERT_NEAR(41.1, mesh->Node(0).x, tolerance);
    ASSERT_NEAR(41.1, mesh->Node(1).x, tolerance);
    ASSERT_NEAR(41.1, mesh->Node(2).x, tolerance);
    ASSERT_NEAR(41.1, mesh->Node(3).x, tolerance);

    ASSERT_NEAR(41.103, mesh->Node(4).x, tolerance);
    ASSERT_NEAR(41.103, mesh->Node(5).x, tolerance);
    ASSERT_NEAR(41.103, mesh->Node(6).x, tolerance);
    ASSERT_NEAR(41.103, mesh->Node(7).x, tolerance);

    ASSERT_NEAR(41.106, mesh->Node(8).x, tolerance);
    ASSERT_NEAR(41.106, mesh->Node(9).x, tolerance);
    ASSERT_NEAR(41.106, mesh->Node(10).x, tolerance);
    ASSERT_NEAR(41.106, mesh->Node(11).x, tolerance);

    ASSERT_NEAR(41.109, mesh->Node(12).x, tolerance);
    ASSERT_NEAR(41.109, mesh->Node(13).x, tolerance);
    ASSERT_NEAR(41.109, mesh->Node(14).x, tolerance);
    ASSERT_NEAR(41.109, mesh->Node(15).x, tolerance);

    ASSERT_NEAR(41.1, mesh->Node(0).y, tolerance);
    ASSERT_NEAR(41.103, mesh->Node(1).y, tolerance);
    ASSERT_NEAR(41.106, mesh->Node(2).y, tolerance);
    ASSERT_NEAR(41.109, mesh->Node(3).y, tolerance);

    ASSERT_NEAR(41.1, mesh->Node(4).y, tolerance);
    ASSERT_NEAR(41.103, mesh->Node(5).y, tolerance);
    ASSERT_NEAR(41.106, mesh->Node(6).y, tolerance);
    ASSERT_NEAR(41.109, mesh->Node(7).y, tolerance);

    ASSERT_NEAR(41.1, mesh->Node(8).y, tolerance);
    ASSERT_NEAR(41.103, mesh->Node(9).y, tolerance);
    ASSERT_NEAR(41.106, mesh->Node(10).y, tolerance);
    ASSERT_NEAR(41.109, mesh->Node(11).y, tolerance);

    ASSERT_NEAR(41.1, mesh->Node(12).y, tolerance);
    ASSERT_NEAR(41.103, mesh->Node(13).y, tolerance);
    ASSERT_NEAR(41.106, mesh->Node(14).y, tolerance);
    ASSERT_NEAR(41.109, mesh->Node(15).y, tolerance);
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizationSmallTriangulargridSpherical)
{
    std::vector<Point> nodes{{41.1019592, 41.1072273},
                             {41.1044655, 41.1043587},
                             {41.1051979, 41.1073151},
                             {41.1080132, 41.1046638},
                             {41.1014137, 41.1039963}};

    std::vector<Edge> edges{{4, 0},
                            {4, 1},
                            {1, 0},
                            {1, 3},
                            {3, 2},
                            {2, 1},
                            {0, 2}};

    auto mesh = std::make_shared<Mesh2D>(edges, nodes, Projection::spherical);

    const auto projectToLandBoundaryOption = LandBoundaries::ProjectToLandBoundaryOption::DoNotProjectToLandBoundary;
    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor_at_boundary = 1.0;
    orthogonalizationParameters.areal_to_angle_smoothing_factor = 1.0;

    // no enclosing polygon
    auto polygon = std::make_unique<Polygons>();
    std::vector<Point> landBoundary{};
    auto landboundaries = std::make_unique<LandBoundaries>(landBoundary, *mesh, *polygon);

    auto orthogonalizer = std::make_unique<Orthogonalizer>(*mesh);
    auto smoother = std::make_unique<Smoother>(*mesh);

    OrthogonalizationAndSmoothing orthogonalization(*mesh,
                                                    std::move(smoother),
                                                    std::move(orthogonalizer),
                                                    std::move(polygon),
                                                    std::move(landboundaries),
                                                    projectToLandBoundaryOption,
                                                    orthogonalizationParameters);

    [[maybe_unused]] auto undoAction = orthogonalization.Initialize();
    orthogonalization.Compute();

    constexpr double tolerance = 1e-3;

    ASSERT_NEAR(41.1019592285156, mesh->Node(0).x, tolerance);
    ASSERT_NEAR(41.1044654597059, mesh->Node(1).x, tolerance);
    ASSERT_NEAR(41.1051978230878, mesh->Node(2).x, tolerance);
    ASSERT_NEAR(41.1080131530762, mesh->Node(3).x, tolerance);
    ASSERT_NEAR(41.1014137268066, mesh->Node(4).x, tolerance);

    ASSERT_NEAR(41.1072273254395, mesh->Node(0).y, tolerance);
    ASSERT_NEAR(41.1043586701373, mesh->Node(1).y, tolerance);
    ASSERT_NEAR(41.1073150612170, mesh->Node(2).y, tolerance);
    ASSERT_NEAR(41.1046638488770, mesh->Node(3).y, tolerance);
    ASSERT_NEAR(41.1039962768555, mesh->Node(4).y, tolerance);
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizationSmallTriangulargridSphericalWithUndo)
{
    std::vector<Point> nodes{{41.1019592, 41.1072273},
                             {41.1044655, 41.1043587},
                             {41.1051979, 41.1073151},
                             {41.1080132, 41.1046638},
                             {41.1014137, 41.1039963}};

    std::vector<Edge> edges{{4, 0},
                            {4, 1},
                            {1, 0},
                            {1, 3},
                            {3, 2},
                            {2, 1},
                            {0, 2}};

    auto mesh = std::make_shared<Mesh2D>(edges, nodes, Projection::spherical);

    const auto projectToLandBoundaryOption = LandBoundaries::ProjectToLandBoundaryOption::DoNotProjectToLandBoundary;
    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor_at_boundary = 1.0;
    orthogonalizationParameters.areal_to_angle_smoothing_factor = 1.0;

    // no enclosing polygon
    auto polygon = std::make_unique<Polygons>();
    std::vector<Point> landBoundary{};
    auto landboundaries = std::make_unique<LandBoundaries>(landBoundary, *mesh, *polygon);

    auto orthogonalizer = std::make_unique<Orthogonalizer>(*mesh);
    auto smoother = std::make_unique<Smoother>(*mesh);

    const std::vector<meshkernel::Point> meshNodes = mesh->Nodes();

    OrthogonalizationAndSmoothing orthogonalization(*mesh,
                                                    std::move(smoother),
                                                    std::move(orthogonalizer),
                                                    std::move(polygon),
                                                    std::move(landboundaries),
                                                    projectToLandBoundaryOption,
                                                    orthogonalizationParameters);

    std::unique_ptr<meshkernel::UndoAction> undoAction = orthogonalization.Initialize();
    orthogonalization.Compute();

    constexpr double tolerance = 1e-3;

    ASSERT_NEAR(41.1019592285156, mesh->Node(0).x, tolerance);
    ASSERT_NEAR(41.1044654597059, mesh->Node(1).x, tolerance);
    ASSERT_NEAR(41.1051978230878, mesh->Node(2).x, tolerance);
    ASSERT_NEAR(41.1080131530762, mesh->Node(3).x, tolerance);
    ASSERT_NEAR(41.1014137268066, mesh->Node(4).x, tolerance);

    ASSERT_NEAR(41.1072273254395, mesh->Node(0).y, tolerance);
    ASSERT_NEAR(41.1043586701373, mesh->Node(1).y, tolerance);
    ASSERT_NEAR(41.1073150612170, mesh->Node(2).y, tolerance);
    ASSERT_NEAR(41.1046638488770, mesh->Node(3).y, tolerance);
    ASSERT_NEAR(41.1039962768555, mesh->Node(4).y, tolerance);

    undoAction->Restore();

    EXPECT_EQ(meshNodes.size(), mesh->GetNumNodes());

    for (meshkernel::UInt i = 0; i < meshNodes.size(); ++i)
    {
        EXPECT_EQ(meshNodes[i].x, mesh->Node(i).x);
        EXPECT_EQ(meshNodes[i].y, mesh->Node(i).y);
    }
}

TEST(OrthogonalizationAndSmoothing, RefineUndoThenOrthogonalise)
{
    auto mesh = MakeRectangularMeshForTestingRand(20, 20, 1.0, Projection::cartesian);

    const std::vector<meshkernel::Point> originalPoints = mesh->Nodes();

    // sample points
    std::vector<Sample> samples{
        {14.7153645, 14.5698833, 1.0},
        {24.7033062, 14.4729137, 1.0},
        {15.5396099, 24.2669525, 1.0},
        {23.8305721, 23.9275551, 1.0}};

    auto interpolator = std::make_unique<AveragingInterpolation>(*mesh,
                                                                 samples,
                                                                 AveragingInterpolation::Method::MinAbsValue,
                                                                 Location::Faces,
                                                                 1.0,
                                                                 false,
                                                                 false,
                                                                 1);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 1;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;
    meshRefinementParameters.min_edge_size = 1.0;
    meshRefinementParameters.account_for_samples_outside = 0;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.refinement_type = 2;
    meshRefinementParameters.smoothing_iterations = 0;

    MeshRefinement meshRefinement(*mesh,
                                  std::move(interpolator),
                                  meshRefinementParameters);

    auto refinementUndoAction = meshRefinement.Compute();
    refinementUndoAction->Restore();

    const auto projectToLandBoundaryOption = LandBoundaries::ProjectToLandBoundaryOption::DoNotProjectToLandBoundary;
    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.inner_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.outer_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor_at_boundary = 0.975;
    orthogonalizationParameters.areal_to_angle_smoothing_factor = 1.0;

    std::vector<meshkernel::Point> polygonNodes{{10.5, 10.5}, {30.5, 10.5}, {30.5, 30.5}, {10.5, 30.5}, {10.5, 10.5}};

    auto polygon = std::make_unique<Polygons>(polygonNodes, Projection::cartesian);
    auto orthogonalisationPolygon = std::make_unique<Polygons>(polygonNodes, Projection::cartesian);
    std::vector<Point> landBoundary{};
    auto landboundaries = std::make_unique<LandBoundaries>(landBoundary, *mesh, *polygon);

    auto orthogonalizer = std::make_unique<Orthogonalizer>(*mesh);
    auto smoother = std::make_unique<Smoother>(*mesh);
    OrthogonalizationAndSmoothing orthogonalization(*mesh,
                                                    std::move(smoother),
                                                    std::move(orthogonalizer),
                                                    std::move(orthogonalisationPolygon),
                                                    std::move(landboundaries),
                                                    projectToLandBoundaryOption,
                                                    orthogonalizationParameters);

    [[maybe_unused]] auto undoAction = orthogonalization.Initialize();

    const std::vector<meshkernel::Point>& nodes = mesh->Nodes();
    const auto& enclosure = polygon->Enclosure(0);

    constexpr double tolerance = 1.0e-12;

    orthogonalization.Compute();

    size_t count = 0;
    size_t orignialCount = 0;

    for (size_t i = 0; i < nodes.size(); ++i)
    {
        if (nodes[i].IsValid())
        {
            if (!enclosure.Contains(nodes[i]))
            {
                EXPECT_NEAR(originalPoints[count].x, nodes[i].x, tolerance) << "x-position: " << i << "  " << count << "  " << orignialCount;
                EXPECT_NEAR(originalPoints[count].y, nodes[i].y, tolerance) << "y-position: " << i << "  " << count << "  " << orignialCount;
                ++orignialCount;
            }
            ++count;
        }
    }

    // There are 400 points in the original mesh
    EXPECT_EQ(count, 400);

    // Only check the number of points not subject to orthogonalisation
    EXPECT_EQ(orignialCount, 319);
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizationWithGapsInNodeAndEdgeLists)
{

    std::vector<Point> polygonPoints{{0.0, 0.0}, {100.0, 20.0}, {150.0, 50.0}, {75.0, 75.0}, {-20.0, 30.0}, {0.0, 0.0}};

    const Polygons polygon(polygonPoints, Projection::cartesian);

    std::vector<Point> refinedPolygonPoints = polygon.RefinePolygon(0, 0, 5, 25.0);

    std::unique_ptr<Polygons> refinedPolygon = std::make_unique<Polygons>(refinedPolygonPoints, Projection::cartesian);

    // generate samples in all polygons
    const std::vector<std::vector<Point>> generatedPoints = refinedPolygon->ComputePointsInPolygons();

    meshkernel::Mesh2D mesh(generatedPoints[0], *refinedPolygon, Projection::cartesian);

    // Create some gaps in the node and edge arrays
    auto [nodeId, nodeInsertUndo] = mesh.InsertNode({0.5, 0.5});
    auto originNodeId = mesh.FindNodeCloseToAPoint({0.0, 0.0}, 0.1);
    [[maybe_unused]] auto [edgeId, edgeInsertUndo] = mesh.ConnectNodes(nodeId, originNodeId);
    [[maybe_unused]] auto nodeRemovaUndo = mesh.DeleteNode(nodeId);
    mesh.Administrate();

    std::unique_ptr<LandBoundaries> boundary = std::make_unique<LandBoundaries>(refinedPolygonPoints, mesh, *refinedPolygon);

    FlipEdges flipEdges1(mesh, *boundary, true, false);

    const auto projectToLandBoundaryOption = LandBoundaries::ProjectToLandBoundaryOption::DoNotProjectToLandBoundary;
    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor_at_boundary = 0.975;
    orthogonalizationParameters.areal_to_angle_smoothing_factor = 1.0;

    FlipEdges flipEdges(mesh, *boundary, true, false);

    OrthogonalizationAndSmoothing orthogonalization(mesh,
                                                    std::make_unique<Smoother>(mesh),
                                                    std::make_unique<Orthogonalizer>(mesh),
                                                    std::move(refinedPolygon),
                                                    std::move(boundary),
                                                    projectToLandBoundaryOption,
                                                    orthogonalizationParameters);

    const std::vector<double> expectedX{28.197593689053271,
                                        45.122094904304319,
                                        67.434041291301682,
                                        91.65019268878028,
                                        121.18012351239659,
                                        150.0,
                                        121.65514995517655,
                                        96.144810647136211,
                                        73.392601318383058,
                                        49.383970314068847,
                                        27.985698958218652,
                                        9.2519807557372626,
                                        -13.113266174727251,
                                        8.0469903108904433,
                                        18.108865376000267,
                                        30.057404995092867,
                                        81.572189207388078,
                                        100.72370484044434,
                                        44.08035834002014,
                                        62.494201377406981,
                                        56.221541072263449,
                                        80.952825997664704,
                                        39.468493116531533,
                                        -999.0};

    const std::vector<double> expectedY{5.6395187378106542,
                                        9.0244189808608635,
                                        13.486808258260337,
                                        18.330038537756057,
                                        32.708074107437959,
                                        50.0,
                                        59.448283348274487,
                                        67.951729784287934,
                                        74.238600624497238,
                                        62.866091201401034,
                                        52.730067927577252,
                                        43.856201410612385,
                                        33.262137075129196,
                                        1.6093980621780888,
                                        22.269932674539941,
                                        35.917300713638298,
                                        51.865232389386854,
                                        45.127820397816109,
                                        43.164232704286427,
                                        48.110596487644287,
                                        28.375434217690007,
                                        34.073083926222459,
                                        23.663968978263402,
                                        -999.0};

    const std::vector<meshkernel::UInt> edgeFirst{13, 14, 12, 13, 0, 22, 0, 0, 20, 18, 22, 15,
                                                  22, 1, 14, 11, 10, 11, 15, 18, 9, 10, 18, 19,
                                                  18, 1, 20, 21, 21, 1, 2, 21, 16, 19, 21, 3, 3,
                                                  2, 21, 17, 4, 16, 8, 16, 7, 8, 4, 6, 4, 5, 6, 7,
                                                  meshkernel::constants::missing::uintValue};

    const std::vector<meshkernel::UInt> edgeSecond{14, 12, 13, 0, 14, 14, 1, 22, 18, 22, 20, 14,
                                                   15, 22, 11, 12, 11, 15, 10, 9, 10, 18, 15, 9,
                                                   19, 20, 2, 20, 19, 2, 21, 16, 19, 20, 3, 4, 17,
                                                   3, 17, 16, 17, 8, 19, 7, 8, 9, 6, 17, 5, 6, 7, 17,
                                                   meshkernel::constants::missing::uintValue};

    [[maybe_unused]] auto undoAction = orthogonalization.Initialize();

    for (int i = 1; i <= 10; ++i)
    {
        auto flipUndoAction = flipEdges.Compute();

        orthogonalization.Compute();
    }

    const double tolerance = 1.0e-8;

    ASSERT_EQ(static_cast<size_t>(mesh.GetNumNodes()), expectedX.size());
    ASSERT_EQ(static_cast<size_t>(mesh.GetNumEdges()), edgeFirst.size());

    for (meshkernel::UInt i = 0; i < mesh.GetNumNodes(); ++i)
    {
        EXPECT_NEAR(expectedX[i], mesh.Node(i).x, tolerance);
        EXPECT_NEAR(expectedY[i], mesh.Node(i).y, tolerance);
    }

    for (meshkernel::UInt i = 0; i < mesh.GetNumEdges(); ++i)
    {
        EXPECT_EQ(edgeFirst[i], mesh.GetEdge(i).first);
        EXPECT_EQ(edgeSecond[i], mesh.GetEdge(i).second);
    }
}
