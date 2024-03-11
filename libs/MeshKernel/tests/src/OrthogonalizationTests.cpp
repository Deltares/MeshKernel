#include <gtest/gtest.h>

#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/LandBoundaries.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/OrthogonalizationAndSmoothing.hpp>
#include <MeshKernel/Orthogonalizer.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/Smoother.hpp>
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

    orthogonalization.Initialize();

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

    orthogonalization.Initialize();

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

    orthogonalization.Initialize();
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

    orthogonalization.Initialize();

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
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizationMediumTriangularGrid)
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

    OrthogonalizationAndSmoothing orthogonalization(*mesh,
                                                    std::move(smoother),
                                                    std::move(orthogonalizer),
                                                    std::move(polygon),
                                                    std::move(landBoundaries),
                                                    projectToLandBoundaryOption,
                                                    orthogonalizationParameters);

    orthogonalization.Initialize();
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

    orthogonalization.Initialize();

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

    orthogonalization.Initialize();

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

    orthogonalization.Initialize();
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

    orthogonalization.Initialize();
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
