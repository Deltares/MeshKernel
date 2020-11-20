#include <MeshKernel/Mesh.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/Constants.hpp>
#include <MeshKernel/OrthogonalizationAndSmoothing.hpp>
#include <MeshKernel/Smoother.hpp>
#include <MeshKernel/LandBoundaries.hpp>
#include <MeshKernel/Orthogonalizer.hpp>
#include <TestUtils/MakeMeshes.hpp>
#include <gtest/gtest.h>
#include <chrono>

#if defined(_WIN32)
#include <Windows.h>
#endif

TEST(OrthogonalizationAndSmoothing, OrthogonalizationOneQuadOneTriangle)
{
    // Preparation
    std::vector<meshkernel::Point> nodes{{0.0, 0.0},
                                         {0.0, 10.0},
                                         {10.0, 0.0},
                                         {10.0, 10.0},
                                         {20.0, 0.0}};

    std::vector<meshkernel::Edge> edges{{1, 0},
                                        {0, 2},
                                        {2, 4},
                                        {4, 3},
                                        {3, 2},
                                        {3, 1}};

    int projectToLandBoundaryOption = 0;
    meshkernelapi::OrthogonalizationParametersNative orthogonalizationParametersNative;
    orthogonalizationParametersNative.InnerIterations = 2;
    orthogonalizationParametersNative.BoundaryIterations = 25;
    orthogonalizationParametersNative.OuterIterations = 25;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactor = 0.975;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactorBoundary = 0.975;
    orthogonalizationParametersNative.Smoothorarea = 1.0;

    // Execute
    auto mesh = std::make_shared<meshkernel::Mesh>();
    mesh->Set(edges, nodes, meshkernel::Projections::cartesian);

    auto orthogonalizer = std::make_shared<meshkernel::Orthogonalizer>(mesh);
    auto smoother = std::make_shared<meshkernel::Smoother>(mesh);
    auto polygon = std::make_shared<meshkernel::Polygons>();

    std::vector<meshkernel::Point> landBoundary{};
    auto landboundaries = std::make_shared<meshkernel::LandBoundaries>(landBoundary, mesh, polygon);

    meshkernel::OrthogonalizationAndSmoothing orthogonalization(mesh,
                                                                smoother,
                                                                orthogonalizer,
                                                                polygon,
                                                                landboundaries,
                                                                projectToLandBoundaryOption,
                                                                orthogonalizationParametersNative);

    orthogonalization.Initialize();

    orthogonalization.Compute();

    // Assert
    constexpr double tolerance = 1e-8;
    ASSERT_NEAR(0.0, mesh->m_nodes[0].x, tolerance);
    ASSERT_NEAR(0.0, mesh->m_nodes[1].x, tolerance);
    ASSERT_NEAR(10.0, mesh->m_nodes[2].x, tolerance);
    ASSERT_NEAR(10.0, mesh->m_nodes[3].x, tolerance);
    ASSERT_NEAR(20.0, mesh->m_nodes[4].x, tolerance);

    ASSERT_NEAR(0.0, mesh->m_nodes[0].y, tolerance);
    ASSERT_NEAR(10.0, mesh->m_nodes[1].y, tolerance);
    ASSERT_NEAR(0.0, mesh->m_nodes[2].y, tolerance);
    ASSERT_NEAR(10.0, mesh->m_nodes[3].y, tolerance);
    ASSERT_NEAR(0.0, mesh->m_nodes[4].y, tolerance);
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizationSmallTriangularGrid)
{

    // now build node-edge mapping
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/SmallTriangularGrid_net.nc");
    int projectToLandBoundaryOption = 0;
    meshkernelapi::OrthogonalizationParametersNative orthogonalizationParametersNative;
    orthogonalizationParametersNative.OuterIterations = 2;
    orthogonalizationParametersNative.BoundaryIterations = 25;
    orthogonalizationParametersNative.InnerIterations = 25;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactor = 0.975;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactorBoundary = 1.0;
    orthogonalizationParametersNative.Smoothorarea = 1.0;

    auto orthogonalizer = std::make_shared<meshkernel::Orthogonalizer>(mesh);
    auto smoother = std::make_shared<meshkernel::Smoother>(mesh);
    auto polygon = std::make_shared<meshkernel::Polygons>();

    std::vector<meshkernel::Point> landBoundary;
    auto landboundaries = std::make_shared<meshkernel::LandBoundaries>(landBoundary, mesh, polygon);

    meshkernel::OrthogonalizationAndSmoothing orthogonalization(mesh,
                                                                smoother,
                                                                orthogonalizer,
                                                                polygon,
                                                                landboundaries,
                                                                projectToLandBoundaryOption,
                                                                orthogonalizationParametersNative);

    orthogonalization.Initialize();

    orthogonalization.Compute();

    constexpr double tolerance = 1e-2;

    ASSERT_NEAR(325.590101919525, mesh->m_nodes[0].x, tolerance);
    ASSERT_NEAR(229.213730481198, mesh->m_nodes[1].x, tolerance);
    ASSERT_NEAR(263.439319753147, mesh->m_nodes[2].x, tolerance);
    ASSERT_NEAR(429.191105834504, mesh->m_nodes[3].x, tolerance);
    ASSERT_NEAR(535.865215426468, mesh->m_nodes[4].x, tolerance);
    ASSERT_NEAR(503.753784179688, mesh->m_nodes[5].x, tolerance);
    ASSERT_NEAR(354.048340705929, mesh->m_nodes[6].x, tolerance);
    ASSERT_NEAR(346.790050854504, mesh->m_nodes[7].x, tolerance);
    ASSERT_NEAR(315.030130405285, mesh->m_nodes[8].x, tolerance);
    ASSERT_NEAR(424.314957449766, mesh->m_nodes[9].x, tolerance);

    ASSERT_NEAR(455.319334078551, mesh->m_nodes[0].y, tolerance);
    ASSERT_NEAR(362.573521507281, mesh->m_nodes[1].y, tolerance);
    ASSERT_NEAR(241.096458631763, mesh->m_nodes[2].y, tolerance);
    ASSERT_NEAR(211.483073921775, mesh->m_nodes[3].y, tolerance);
    ASSERT_NEAR(311.401495506714, mesh->m_nodes[4].y, tolerance);
    ASSERT_NEAR(432.379974365234, mesh->m_nodes[5].y, tolerance);
    ASSERT_NEAR(458.064836627594, mesh->m_nodes[6].y, tolerance);
    ASSERT_NEAR(405.311585650679, mesh->m_nodes[7].y, tolerance);
    ASSERT_NEAR(319.612138503550, mesh->m_nodes[8].y, tolerance);
    ASSERT_NEAR(327.102805172725, mesh->m_nodes[9].y, tolerance);
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizationSmallTriangularGridAsNcFile)
{

    // now build node-edge mapping
    auto mesh = MakeSmallSizeTriangularMeshForTestingAsNcFile();

    int projectToLandBoundaryOption = 0;
    meshkernelapi::OrthogonalizationParametersNative orthogonalizationParametersNative;
    orthogonalizationParametersNative.OuterIterations = 2;
    orthogonalizationParametersNative.BoundaryIterations = 25;
    orthogonalizationParametersNative.InnerIterations = 25;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactor = 0.975;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactorBoundary = 1.0;
    orthogonalizationParametersNative.Smoothorarea = 1.0;

    auto orthogonalizer = std::make_shared<meshkernel::Orthogonalizer>(mesh);
    auto smoother = std::make_shared<meshkernel::Smoother>(mesh);
    auto polygon = std::make_shared<meshkernel::Polygons>();

    std::vector<meshkernel::Point> landBoundary;
    auto landboundaries = std::make_shared<meshkernel::LandBoundaries>(landBoundary, mesh, polygon);

    meshkernel::OrthogonalizationAndSmoothing orthogonalization(mesh,
                                                                smoother,
                                                                orthogonalizer,
                                                                polygon,
                                                                landboundaries,
                                                                projectToLandBoundaryOption,
                                                                orthogonalizationParametersNative);

    orthogonalization.Initialize();
    orthogonalization.Compute();

    constexpr double tolerance = 1e-2;

    ASSERT_NEAR(325.590101919525, mesh->m_nodes[0].x, tolerance);
    ASSERT_NEAR(229.213730481198, mesh->m_nodes[1].x, tolerance);
    ASSERT_NEAR(263.439319753147, mesh->m_nodes[2].x, tolerance);
    ASSERT_NEAR(429.191105834504, mesh->m_nodes[3].x, tolerance);
    ASSERT_NEAR(535.865215426468, mesh->m_nodes[4].x, tolerance);
    ASSERT_NEAR(503.753784179688, mesh->m_nodes[5].x, tolerance);
    ASSERT_NEAR(354.048340705929, mesh->m_nodes[6].x, tolerance);
    ASSERT_NEAR(346.790050854504, mesh->m_nodes[7].x, tolerance);
    ASSERT_NEAR(315.030130405285, mesh->m_nodes[8].x, tolerance);
    ASSERT_NEAR(424.314957449766, mesh->m_nodes[9].x, tolerance);

    ASSERT_NEAR(455.319334078551, mesh->m_nodes[0].y, tolerance);
    ASSERT_NEAR(362.573521507281, mesh->m_nodes[1].y, tolerance);
    ASSERT_NEAR(241.096458631763, mesh->m_nodes[2].y, tolerance);
    ASSERT_NEAR(211.483073921775, mesh->m_nodes[3].y, tolerance);
    ASSERT_NEAR(311.401495506714, mesh->m_nodes[4].y, tolerance);
    ASSERT_NEAR(432.379974365234, mesh->m_nodes[5].y, tolerance);
    ASSERT_NEAR(458.064836627594, mesh->m_nodes[6].y, tolerance);
    ASSERT_NEAR(405.311585650679, mesh->m_nodes[7].y, tolerance);
    ASSERT_NEAR(319.612138503550, mesh->m_nodes[8].y, tolerance);
    ASSERT_NEAR(327.102805172725, mesh->m_nodes[9].y, tolerance);
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizationMediumTriangularGridWithPolygon)
{
    // now build node-edge mapping
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/SmallTriangularGrid_net.nc");

    int projectToLandBoundaryOption = 0;
    meshkernelapi::OrthogonalizationParametersNative orthogonalizationParametersNative;
    orthogonalizationParametersNative.OuterIterations = 2;
    orthogonalizationParametersNative.BoundaryIterations = 25;
    orthogonalizationParametersNative.InnerIterations = 25;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactor = 0.975;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactorBoundary = 1.0;
    orthogonalizationParametersNative.Smoothorarea = 1.0;

    std::vector<meshkernel::Point> nodes{{342.987518, 471.121002},
                                         {327.640900, 380.846436},
                                         {396.851135, 201.200073},
                                         {514.207581, 203.607407},
                                         {569.274841, 294.483765},
                                         {568.673035, 379.943695},
                                         {515.712158, 458.783478},
                                         {343.288422, 471.722809}};

    auto orthogonalizer = std::make_shared<meshkernel::Orthogonalizer>(mesh);
    auto smoother = std::make_shared<meshkernel::Smoother>(mesh);

    std::vector<meshkernel::Point> landBoundary;
    auto polygon = std::make_shared<meshkernel::Polygons>(nodes, meshkernel::Projections::cartesian);
    auto landboundaries = std::make_shared<meshkernel::LandBoundaries>(landBoundary, mesh, polygon);

    meshkernel::OrthogonalizationAndSmoothing orthogonalization(mesh,
                                                                smoother,
                                                                orthogonalizer,
                                                                polygon,
                                                                landboundaries,
                                                                projectToLandBoundaryOption,
                                                                orthogonalizationParametersNative);

    orthogonalization.Initialize();

    orthogonalization.Compute();

    constexpr double tolerance = 1.8;
    // check the first 10 points

    ASSERT_NEAR(322.252624511719, mesh->m_nodes[0].x, tolerance);
    ASSERT_NEAR(227.002044677734, mesh->m_nodes[1].x, tolerance);
    ASSERT_NEAR(259.252227783203, mesh->m_nodes[2].x, tolerance);
    ASSERT_NEAR(426.680730020925, mesh->m_nodes[3].x, tolerance);
    ASSERT_NEAR(535.742871000462, mesh->m_nodes[4].x, tolerance);
    ASSERT_NEAR(502.614684581472, mesh->m_nodes[5].x, tolerance);
    ASSERT_NEAR(350.424583126127, mesh->m_nodes[6].x, tolerance);
    ASSERT_NEAR(342.832067902036, mesh->m_nodes[7].x, tolerance);
    ASSERT_NEAR(310.300984548069, mesh->m_nodes[8].x, tolerance);
    ASSERT_NEAR(421.791577149564, mesh->m_nodes[9].x, tolerance);

    ASSERT_NEAR(454.880187988281, mesh->m_nodes[0].y, tolerance);
    ASSERT_NEAR(360.379241943359, mesh->m_nodes[1].y, tolerance);
    ASSERT_NEAR(241.878051757812, mesh->m_nodes[2].y, tolerance);
    ASSERT_NEAR(210.624626374983, mesh->m_nodes[3].y, tolerance);
    ASSERT_NEAR(311.862423032534, mesh->m_nodes[4].y, tolerance);
    ASSERT_NEAR(432.575408917275, mesh->m_nodes[5].y, tolerance);
    ASSERT_NEAR(458.587061164954, mesh->m_nodes[6].y, tolerance);
    ASSERT_NEAR(405.390188082498, mesh->m_nodes[7].y, tolerance);
    ASSERT_NEAR(319.410057398020, mesh->m_nodes[8].y, tolerance);
    ASSERT_NEAR(327.001109057344, mesh->m_nodes[9].y, tolerance);
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizationMediumTriangularGrid)
{
    // now build node-edge mapping
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/TestOrthogonalizationMediumTriangularGrid_net.nc");

    int projectToLandBoundaryOption = 0;
    meshkernelapi::OrthogonalizationParametersNative orthogonalizationParametersNative;
    orthogonalizationParametersNative.OuterIterations = 2;
    orthogonalizationParametersNative.BoundaryIterations = 25;
    orthogonalizationParametersNative.InnerIterations = 25;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactor = 0.975;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactorBoundary = 0.5;
    orthogonalizationParametersNative.Smoothorarea = 1.0;

    auto orthogonalizer = std::make_shared<meshkernel::Orthogonalizer>(mesh);
    auto smoother = std::make_shared<meshkernel::Smoother>(mesh);
    auto polygon = std::make_shared<meshkernel::Polygons>();
    std::vector<meshkernel::Point> landBoundary;
    auto landBoundaries = std::make_shared<meshkernel::LandBoundaries>(landBoundary, mesh, polygon);

    meshkernel::OrthogonalizationAndSmoothing orthogonalization(mesh,
                                                                smoother,
                                                                orthogonalizer,
                                                                polygon,
                                                                landBoundaries,
                                                                projectToLandBoundaryOption,
                                                                orthogonalizationParametersNative);

    orthogonalization.Initialize();
    orthogonalization.Compute();

    constexpr double tolerance = 1.2;

    // check the first 10 points
    ASSERT_NEAR(68.771705432835475, mesh->m_nodes[0].x, tolerance);
    ASSERT_NEAR(169.49338272334273, mesh->m_nodes[1].x, tolerance);
    ASSERT_NEAR(262.80128484924921, mesh->m_nodes[2].x, tolerance);
    ASSERT_NEAR(361.60010033352023, mesh->m_nodes[3].x, tolerance);
    ASSERT_NEAR(468.13991812406925, mesh->m_nodes[4].x, tolerance);
    ASSERT_NEAR(549.89461192844624, mesh->m_nodes[5].x, tolerance);
    ASSERT_NEAR(653.02704974527421, mesh->m_nodes[6].x, tolerance);
    ASSERT_NEAR(747.81537706979441, mesh->m_nodes[7].x, tolerance);
    ASSERT_NEAR(853.40641427112951, mesh->m_nodes[8].x, tolerance);
    ASSERT_NEAR(938.69752431820143, mesh->m_nodes[9].x, tolerance);

    ASSERT_NEAR(1399.7751472360221, mesh->m_nodes[0].y, tolerance);
    ASSERT_NEAR(1426.5945287630802, mesh->m_nodes[1].y, tolerance);
    ASSERT_NEAR(1451.4398281457179, mesh->m_nodes[2].y, tolerance);
    ASSERT_NEAR(1477.7472050498141, mesh->m_nodes[3].y, tolerance);
    ASSERT_NEAR(1506.1157955857589, mesh->m_nodes[4].y, tolerance);
    ASSERT_NEAR(1527.8847968946166, mesh->m_nodes[5].y, tolerance);
    ASSERT_NEAR(1555.3460969050145, mesh->m_nodes[6].y, tolerance);
    ASSERT_NEAR(1580.5855923464549, mesh->m_nodes[7].y, tolerance);
    ASSERT_NEAR(1608.7015489976982, mesh->m_nodes[8].y, tolerance);
    ASSERT_NEAR(1631.412199601948, mesh->m_nodes[9].y, tolerance);
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizationFourQuads)
{
    auto mesh = MakeRectangularMeshForTesting(3, 3, 1.0, meshkernel::Projections::cartesian);

    int projectToLandBoundaryOption = 0;
    meshkernelapi::OrthogonalizationParametersNative orthogonalizationParametersNative;
    orthogonalizationParametersNative.InnerIterations = 2;
    orthogonalizationParametersNative.BoundaryIterations = 25;
    orthogonalizationParametersNative.OuterIterations = 25;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactor = 0.975;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactorBoundary = 0.975;
    orthogonalizationParametersNative.Smoothorarea = 1.0;

    auto polygon = std::make_shared<meshkernel::Polygons>();
    std::vector<meshkernel::Point> landBoundary{};
    auto landboundaries = std::make_shared<meshkernel::LandBoundaries>(landBoundary, mesh, polygon);

    auto orthogonalizer = std::make_shared<meshkernel::Orthogonalizer>(mesh);
    auto smoother = std::make_shared<meshkernel::Smoother>(mesh);
    meshkernel::OrthogonalizationAndSmoothing orthogonalization(mesh,
                                                                smoother,
                                                                orthogonalizer,
                                                                polygon,
                                                                landboundaries,
                                                                projectToLandBoundaryOption,
                                                                orthogonalizationParametersNative);

    orthogonalization.Initialize();

    orthogonalization.Compute();
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizeAndSnapToLandBoundaries)
{
    // Prepare
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/SmallTriangularGrid_net.nc");

    // the land boundary to use
    std::vector<meshkernel::Point> landBoundary{{235.561218, 290.571899},
                                                {265.953522, 436.515747},
                                                {429.349854, 450.959656},
                                                {535.271545, 386.262909},
                                                {meshkernel::doubleMissingValue, meshkernel::doubleMissingValue},
                                                {246.995941, 262.285858},
                                                {351.112183, 237.309906},
                                                {443.191895, 262.285858},
                                                {553.627319, 327.283539},
                                                {meshkernel::doubleMissingValue, meshkernel::doubleMissingValue}};

    // snap to land boundaries
    int projectToLandBoundaryOption = 2;
    meshkernelapi::OrthogonalizationParametersNative orthogonalizationParametersNative;
    orthogonalizationParametersNative.OuterIterations = 2;
    orthogonalizationParametersNative.BoundaryIterations = 25;
    orthogonalizationParametersNative.InnerIterations = 25;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactor = 0.975;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactorBoundary = 0.975;
    orthogonalizationParametersNative.Smoothorarea = 1.0;

    // no enclosing polygon
    auto polygon = std::make_shared<meshkernel::Polygons>();
    auto landboundaries = std::make_shared<meshkernel::LandBoundaries>(landBoundary, mesh, polygon);

    auto orthogonalizer = std::make_shared<meshkernel::Orthogonalizer>(mesh);
    auto smoother = std::make_shared<meshkernel::Smoother>(mesh);
    meshkernel::OrthogonalizationAndSmoothing orthogonalization(mesh,
                                                                smoother,
                                                                orthogonalizer,
                                                                polygon,
                                                                landboundaries,
                                                                projectToLandBoundaryOption,
                                                                orthogonalizationParametersNative);

    orthogonalization.Initialize();

    orthogonalization.Compute();

    // check the values
    constexpr double tolerance = 0.15;
    ASSERT_NEAR(313.081472564480, mesh->m_nodes[0].x, tolerance);
    ASSERT_NEAR(253.641466857330, mesh->m_nodes[1].x, tolerance);
    ASSERT_NEAR(254.777224294204, mesh->m_nodes[2].x, tolerance);
    ASSERT_NEAR(443.191895000000, mesh->m_nodes[3].x, tolerance);
    ASSERT_NEAR(535.240231516760, mesh->m_nodes[4].x, tolerance);
    ASSERT_NEAR(480.436129612752, mesh->m_nodes[5].x, tolerance);
    ASSERT_NEAR(345.948240805397, mesh->m_nodes[6].x, tolerance);
    ASSERT_NEAR(342.668434889472, mesh->m_nodes[7].x, tolerance);
    ASSERT_NEAR(318.414413615199, mesh->m_nodes[8].x, tolerance);
    ASSERT_NEAR(424.616311031376, mesh->m_nodes[9].x, tolerance);

    ASSERT_NEAR(440.681763586650, mesh->m_nodes[0].y, tolerance);
    ASSERT_NEAR(377.393256506700, mesh->m_nodes[1].y, tolerance);
    ASSERT_NEAR(260.419242817573, mesh->m_nodes[2].y, tolerance);
    ASSERT_NEAR(262.285858000000, mesh->m_nodes[3].y, tolerance);
    ASSERT_NEAR(316.461666783032, mesh->m_nodes[4].y, tolerance);
    ASSERT_NEAR(419.756265860671, mesh->m_nodes[5].y, tolerance);
    ASSERT_NEAR(443.587120174434, mesh->m_nodes[6].y, tolerance);
    ASSERT_NEAR(402.913858250569, mesh->m_nodes[7].y, tolerance);
    ASSERT_NEAR(336.831643075189, mesh->m_nodes[8].y, tolerance);
    ASSERT_NEAR(340.875100904741, mesh->m_nodes[9].y, tolerance);
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizationSphericalRectangular)
{
    //1 Setup
    auto mesh = MakeRectangularMeshForTesting(4, 4, 0.003, meshkernel::Projections::spherical, {41.1, 41.1});

    int projectToLandBoundaryOption = 0;
    meshkernelapi::OrthogonalizationParametersNative orthogonalizationParametersNative;
    orthogonalizationParametersNative.OuterIterations = 2;
    orthogonalizationParametersNative.BoundaryIterations = 25;
    orthogonalizationParametersNative.InnerIterations = 25;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactor = 0.975;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactorBoundary = 1.0;
    orthogonalizationParametersNative.Smoothorarea = 1.0;

    auto orthogonalizer = std::make_shared<meshkernel::Orthogonalizer>(mesh);
    auto smoother = std::make_shared<meshkernel::Smoother>(mesh);

    auto polygon = std::make_shared<meshkernel::Polygons>();
    std::vector<meshkernel::Point> landBoundary{};
    auto landboundaries = std::make_shared<meshkernel::LandBoundaries>(landBoundary, mesh, polygon);

    meshkernel::OrthogonalizationAndSmoothing orthogonalization(mesh,
                                                                smoother,
                                                                orthogonalizer,
                                                                polygon,
                                                                landboundaries,
                                                                projectToLandBoundaryOption,
                                                                orthogonalizationParametersNative);

    orthogonalization.Initialize();
    orthogonalization.Compute();

    // check the values
    constexpr double tolerance = 1e-6;
    ASSERT_NEAR(41.1, mesh->m_nodes[0].x, tolerance);
    ASSERT_NEAR(41.1, mesh->m_nodes[1].x, tolerance);
    ASSERT_NEAR(41.1, mesh->m_nodes[2].x, tolerance);
    ASSERT_NEAR(41.1, mesh->m_nodes[3].x, tolerance);

    ASSERT_NEAR(41.103, mesh->m_nodes[4].x, tolerance);
    ASSERT_NEAR(41.103, mesh->m_nodes[5].x, tolerance);
    ASSERT_NEAR(41.103, mesh->m_nodes[6].x, tolerance);
    ASSERT_NEAR(41.103, mesh->m_nodes[7].x, tolerance);

    ASSERT_NEAR(41.106, mesh->m_nodes[8].x, tolerance);
    ASSERT_NEAR(41.106, mesh->m_nodes[9].x, tolerance);
    ASSERT_NEAR(41.106, mesh->m_nodes[10].x, tolerance);
    ASSERT_NEAR(41.106, mesh->m_nodes[11].x, tolerance);

    ASSERT_NEAR(41.109, mesh->m_nodes[12].x, tolerance);
    ASSERT_NEAR(41.109, mesh->m_nodes[13].x, tolerance);
    ASSERT_NEAR(41.109, mesh->m_nodes[14].x, tolerance);
    ASSERT_NEAR(41.109, mesh->m_nodes[15].x, tolerance);

    ASSERT_NEAR(41.1, mesh->m_nodes[0].y, tolerance);
    ASSERT_NEAR(41.103, mesh->m_nodes[1].y, tolerance);
    ASSERT_NEAR(41.106, mesh->m_nodes[2].y, tolerance);
    ASSERT_NEAR(41.109, mesh->m_nodes[3].y, tolerance);

    ASSERT_NEAR(41.1, mesh->m_nodes[4].y, tolerance);
    ASSERT_NEAR(41.103, mesh->m_nodes[5].y, tolerance);
    ASSERT_NEAR(41.106, mesh->m_nodes[6].y, tolerance);
    ASSERT_NEAR(41.109, mesh->m_nodes[7].y, tolerance);

    ASSERT_NEAR(41.1, mesh->m_nodes[8].y, tolerance);
    ASSERT_NEAR(41.103, mesh->m_nodes[9].y, tolerance);
    ASSERT_NEAR(41.106, mesh->m_nodes[10].y, tolerance);
    ASSERT_NEAR(41.109, mesh->m_nodes[11].y, tolerance);

    ASSERT_NEAR(41.1, mesh->m_nodes[12].y, tolerance);
    ASSERT_NEAR(41.103, mesh->m_nodes[13].y, tolerance);
    ASSERT_NEAR(41.106, mesh->m_nodes[14].y, tolerance);
    ASSERT_NEAR(41.109, mesh->m_nodes[15].y, tolerance);
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizationSmallTriangulargridSpherical)
{
    std::vector<meshkernel::Point> nodes{{41.1019592, 41.1072273},
                                         {41.1044655, 41.1043587},
                                         {41.1051979, 41.1073151},
                                         {41.1080132, 41.1046638},
                                         {41.1014137, 41.1039963}};

    std::vector<meshkernel::Edge> edges{{4, 0},
                                        {4, 1},
                                        {1, 0},
                                        {1, 3},
                                        {3, 2},
                                        {2, 1},
                                        {0, 2}};

    auto mesh = std::make_shared<meshkernel::Mesh>();
    mesh->Set(edges, nodes, meshkernel::Projections::spherical);

    int projectToLandBoundaryOption = 0;
    meshkernelapi::OrthogonalizationParametersNative orthogonalizationParametersNative;
    orthogonalizationParametersNative.OuterIterations = 2;
    orthogonalizationParametersNative.BoundaryIterations = 25;
    orthogonalizationParametersNative.InnerIterations = 25;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactor = 0.975;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactorBoundary = 1.0;
    orthogonalizationParametersNative.Smoothorarea = 1.0;

    // no enclosing polygon
    auto polygon = std::make_shared<meshkernel::Polygons>();
    std::vector<meshkernel::Point> landBoundary{};
    auto landboundaries = std::make_shared<meshkernel::LandBoundaries>(landBoundary, mesh, polygon);

    auto orthogonalizer = std::make_shared<meshkernel::Orthogonalizer>(mesh);
    auto smoother = std::make_shared<meshkernel::Smoother>(mesh);

    meshkernel::OrthogonalizationAndSmoothing orthogonalization(mesh,
                                                                smoother,
                                                                orthogonalizer,
                                                                polygon,
                                                                landboundaries,
                                                                projectToLandBoundaryOption,
                                                                orthogonalizationParametersNative);

    orthogonalization.Initialize();
    orthogonalization.Compute();

    constexpr double tolerance = 1e-3;

    ASSERT_NEAR(41.1019592285156, mesh->m_nodes[0].x, tolerance);
    ASSERT_NEAR(41.1044654597059, mesh->m_nodes[1].x, tolerance);
    ASSERT_NEAR(41.1051978230878, mesh->m_nodes[2].x, tolerance);
    ASSERT_NEAR(41.1080131530762, mesh->m_nodes[3].x, tolerance);
    ASSERT_NEAR(41.1014137268066, mesh->m_nodes[4].x, tolerance);

    ASSERT_NEAR(41.1072273254395, mesh->m_nodes[0].y, tolerance);
    ASSERT_NEAR(41.1043586701373, mesh->m_nodes[1].y, tolerance);
    ASSERT_NEAR(41.1073150612170, mesh->m_nodes[2].y, tolerance);
    ASSERT_NEAR(41.1046638488770, mesh->m_nodes[3].y, tolerance);
    ASSERT_NEAR(41.1039962768555, mesh->m_nodes[4].y, tolerance);
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizationMeshWithEdgeWithNoFaces)
{
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/EdgesWithNoFaces_net.nc");
    ASSERT_EQ(mesh->GetNumNodes(), 376);
    ASSERT_EQ(mesh->GetNumEdges(), 698);
    ASSERT_EQ(mesh->GetNumFaces(), 319);

    int projectToLandBoundaryOption = 0;
    meshkernelapi::OrthogonalizationParametersNative orthogonalizationParametersNative;
    orthogonalizationParametersNative.InnerIterations = 25;
    orthogonalizationParametersNative.BoundaryIterations = 25;
    orthogonalizationParametersNative.OuterIterations = 1;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactor = 0.975;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactorBoundary = 0.975;
    orthogonalizationParametersNative.Smoothorarea = 1.0;

    auto polygon = std::make_shared<meshkernel::Polygons>();
    std::vector<meshkernel::Point> landBoundary{};
    auto landboundaries = std::make_shared<meshkernel::LandBoundaries>(landBoundary, mesh, polygon);

    auto orthogonalizer = std::make_shared<meshkernel::Orthogonalizer>(mesh);
    auto smoother = std::make_shared<meshkernel::Smoother>(mesh);
    meshkernel::OrthogonalizationAndSmoothing orthogonalization(mesh,
                                                                smoother,
                                                                orthogonalizer,
                                                                polygon,
                                                                landboundaries,
                                                                projectToLandBoundaryOption,
                                                                orthogonalizationParametersNative);

    orthogonalization.Initialize();
    orthogonalization.Compute();
}
