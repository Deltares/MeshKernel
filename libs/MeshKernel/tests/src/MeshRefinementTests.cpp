#include "MeshKernel/BilinearInterpolationOnGriddedSamples.hpp"
#include "MeshKernel/CasulliRefinement.hpp"
#include "MeshKernel/SamplesHessianCalculator.hpp"

#include <fstream>

#include <gtest/gtest.h>

#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/MeshRefinement.hpp"
#include "MeshKernel/Parameters.hpp"
#include "MeshKernel/Polygons.hpp"
#include "TestUtils/Definitions.hpp"
#include "TestUtils/MakeMeshes.hpp"
#include "TestUtils/SampleFileReader.hpp"
#include "TestUtils/SampleGenerator.hpp"

#include <MeshKernel/Operations.hpp>

#include <TestUtils/MakeCurvilinearGrids.hpp>

using namespace meshkernel;

TEST(MeshRefinement, MeshRefinementRefinementLevels_OnFourByFourWithFourSamples_ShouldRefinemesh)
{
    auto mesh = MakeRectangularMeshForTesting(5, 5, 10.0, Projection::cartesian);

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

    meshRefinement.Compute();

    // 3 Validation edges connecting hanging nodes

    // bottom side
    ASSERT_EQ(5, mesh->GetEdge(73).first);
    ASSERT_EQ(25, mesh->GetEdge(73).second);

    ASSERT_EQ(10, mesh->GetEdge(72).first);
    ASSERT_EQ(25, mesh->GetEdge(72).second);

    ASSERT_EQ(10, mesh->GetEdge(77).first);
    ASSERT_EQ(28, mesh->GetEdge(77).second);

    ASSERT_EQ(15, mesh->GetEdge(76).first);
    ASSERT_EQ(28, mesh->GetEdge(76).second);

    // right side
    ASSERT_EQ(21, mesh->GetEdge(81).first);
    ASSERT_EQ(35, mesh->GetEdge(81).second);

    ASSERT_EQ(22, mesh->GetEdge(80).first);
    ASSERT_EQ(35, mesh->GetEdge(80).second);

    ASSERT_EQ(22, mesh->GetEdge(83).first);
    ASSERT_EQ(36, mesh->GetEdge(83).second);

    ASSERT_EQ(23, mesh->GetEdge(82).first);
    ASSERT_EQ(36, mesh->GetEdge(82).second);

    // upper side
    ASSERT_EQ(19, mesh->GetEdge(79).first);
    ASSERT_EQ(30, mesh->GetEdge(79).second);

    ASSERT_EQ(14, mesh->GetEdge(78).first);
    ASSERT_EQ(30, mesh->GetEdge(78).second);

    ASSERT_EQ(14, mesh->GetEdge(75).first);
    ASSERT_EQ(27, mesh->GetEdge(75).second);

    ASSERT_EQ(9, mesh->GetEdge(74).first);
    ASSERT_EQ(27, mesh->GetEdge(74).second);

    // left side
    ASSERT_EQ(3, mesh->GetEdge(71).first);
    ASSERT_EQ(32, mesh->GetEdge(71).second);

    ASSERT_EQ(2, mesh->GetEdge(70).first);
    ASSERT_EQ(32, mesh->GetEdge(70).second);

    ASSERT_EQ(2, mesh->GetEdge(69).first);
    ASSERT_EQ(31, mesh->GetEdge(69).second);

    ASSERT_EQ(1, mesh->GetEdge(68).first);
    ASSERT_EQ(31, mesh->GetEdge(68).second);

    // total number of edges
    ASSERT_EQ(84, mesh->GetNumEdges());
}

TEST(MeshRefinement, RefinementOnAFourByFourMeshWithSamplesShouldRefine)
{
    auto mesh = MakeRectangularMeshForTesting(4, 4, 500.0, Projection::cartesian);

    // sample points
    std::vector<Sample> samples{
        {14.7153645, 14.5698833, 0.1},
        {24.7033062, 14.4729137, 0.1},
        {15.5396099, 24.2669525, 0.1},
        {23.8305721, 23.9275551, 0.1}};

    auto interpolator = std::make_unique<AveragingInterpolation>(*mesh,
                                                                 samples,
                                                                 AveragingInterpolation::Method::MinAbsValue,
                                                                 Location::Faces,
                                                                 1.0,
                                                                 false,
                                                                 false,
                                                                 1);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 4;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;
    meshRefinementParameters.min_edge_size = 2.0;
    meshRefinementParameters.account_for_samples_outside = 0;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.refinement_type = 1;
    meshRefinementParameters.smoothing_iterations = 0;

    MeshRefinement meshRefinement(*mesh,
                                  std::move(interpolator),
                                  meshRefinementParameters);
    meshRefinement.Compute();

    // Assert number of edges and nodes
    ASSERT_EQ(60, mesh->GetNumEdges());
    ASSERT_EQ(31, mesh->GetNumNodes());

    // Assert edges
    ASSERT_EQ(0, mesh->GetEdge(0).first);
    ASSERT_EQ(26, mesh->GetEdge(0).second);

    ASSERT_EQ(1, mesh->GetEdge(1).first);
    ASSERT_EQ(17, mesh->GetEdge(1).second);

    ASSERT_EQ(2, mesh->GetEdge(2).first);
    ASSERT_EQ(6, mesh->GetEdge(2).second);

    ASSERT_EQ(3, mesh->GetEdge(3).first);
    ASSERT_EQ(7, mesh->GetEdge(3).second);

    ASSERT_EQ(4, mesh->GetEdge(4).first);
    ASSERT_EQ(8, mesh->GetEdge(4).second);

    ASSERT_EQ(5, mesh->GetEdge(5).first);
    ASSERT_EQ(9, mesh->GetEdge(5).second);

    ASSERT_EQ(6, mesh->GetEdge(6).first);
    ASSERT_EQ(10, mesh->GetEdge(6).second);

    ASSERT_EQ(7, mesh->GetEdge(7).first);
    ASSERT_EQ(11, mesh->GetEdge(7).second);

    ASSERT_EQ(8, mesh->GetEdge(8).first);
    ASSERT_EQ(12, mesh->GetEdge(8).second);

    ASSERT_EQ(9, mesh->GetEdge(9).first);
    ASSERT_EQ(13, mesh->GetEdge(9).second);

    ASSERT_EQ(10, mesh->GetEdge(10).first);
    ASSERT_EQ(14, mesh->GetEdge(10).second);

    ASSERT_EQ(11, mesh->GetEdge(11).first);
    ASSERT_EQ(15, mesh->GetEdge(11).second);

    ASSERT_EQ(1, mesh->GetEdge(12).first);
    ASSERT_EQ(18, mesh->GetEdge(12).second);

    ASSERT_EQ(2, mesh->GetEdge(13).first);
    ASSERT_EQ(1, mesh->GetEdge(13).second);

    ASSERT_EQ(3, mesh->GetEdge(14).first);
    ASSERT_EQ(2, mesh->GetEdge(14).second);

    ASSERT_EQ(5, mesh->GetEdge(15).first);
    ASSERT_EQ(19, mesh->GetEdge(15).second);

    ASSERT_EQ(6, mesh->GetEdge(16).first);
    ASSERT_EQ(5, mesh->GetEdge(16).second);

    ASSERT_EQ(7, mesh->GetEdge(17).first);
    ASSERT_EQ(6, mesh->GetEdge(17).second);

    ASSERT_EQ(9, mesh->GetEdge(18).first);
    ASSERT_EQ(8, mesh->GetEdge(18).second);

    ASSERT_EQ(10, mesh->GetEdge(19).first);
    ASSERT_EQ(9, mesh->GetEdge(19).second);

    ASSERT_EQ(11, mesh->GetEdge(20).first);
    ASSERT_EQ(10, mesh->GetEdge(20).second);
}

TEST(MeshRefinement, MeshRefinementRefinementLevels_SmallTriangualMeshTwoSamples_ShouldRefinemesh)
{
    // Prepare
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/SmallTriangularGrid_net.nc");

    // sample points
    std::vector<Sample> samples{
        {359.8657532, 350.3144836, 1.0},
        {387.5152588, 299.2614746, 1.0}};

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
    meshRefinementParameters.min_edge_size = 50.0;
    meshRefinementParameters.account_for_samples_outside = 0;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.refinement_type = 2;
    meshRefinementParameters.smoothing_iterations = 0;

    MeshRefinement meshRefinement(*mesh,
                                  std::move(interpolator),
                                  meshRefinementParameters);

    meshRefinement.Compute();

    // edges connecting hanging nodes
    ASSERT_EQ(10, mesh->GetEdge(32).first);
    ASSERT_EQ(2, mesh->GetEdge(32).second);

    ASSERT_EQ(14, mesh->GetEdge(33).first);
    ASSERT_EQ(4, mesh->GetEdge(33).second);

    ASSERT_EQ(13, mesh->GetEdge(34).first);
    ASSERT_EQ(5, mesh->GetEdge(34).second);

    ASSERT_EQ(11, mesh->GetEdge(31).first);
    ASSERT_EQ(1, mesh->GetEdge(31).second);

    // total number of edges
    ASSERT_EQ(35, mesh->GetNumEdges());
}

TEST(MeshRefinement, RefineBasedOnPolygonTriangularMesh)
{
    // Prepare
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/SmallTriangularGrid_net.nc");

    // Polygon sample
    std::vector<Point> point{
        {399.638169557229, 504.294564030922},
        {361.827403800769, 129.967983041964},
        {651.709941266965, 113.583317880831},
        {666.834247569549, 411.028008498319},
        {410.981399284167, 505.55492288947},
        {399.638169557229, 504.294564030922}};

    auto polygon = Polygons(point, mesh->m_projection);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 1;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;

    MeshRefinement meshRefinement(*mesh, polygon, meshRefinementParameters);
    meshRefinement.Compute();

    // total number of edges
    ASSERT_EQ(15, mesh->GetNumNodes());
    ASSERT_EQ(33, mesh->GetNumEdges());

    // assert on newly generated edges
    ASSERT_EQ(10, mesh->GetEdge(20).first);
    ASSERT_EQ(11, mesh->GetEdge(20).second);

    ASSERT_EQ(11, mesh->GetEdge(21).first);
    ASSERT_EQ(12, mesh->GetEdge(21).second);

    ASSERT_EQ(12, mesh->GetEdge(22).first);
    ASSERT_EQ(10, mesh->GetEdge(22).second);

    ASSERT_EQ(14, mesh->GetEdge(23).first);
    ASSERT_EQ(13, mesh->GetEdge(23).second);

    ASSERT_EQ(13, mesh->GetEdge(24).first);
    ASSERT_EQ(11, mesh->GetEdge(24).second);

    ASSERT_EQ(11, mesh->GetEdge(25).first);
    ASSERT_EQ(14, mesh->GetEdge(25).second);

    ASSERT_EQ(10, mesh->GetEdge(26).first);
    ASSERT_EQ(4, mesh->GetEdge(26).second);
}

TEST(MeshRefinement, ThreeBythreeWithThreeSamplesPerFace)
{
    // Prepare
    auto mesh = MakeRectangularMeshForTesting(4, 4, 10.0, Projection::cartesian);

    // sample points
    std::vector<Sample> samples{
        {2.7091951, 5.4000854, 0.0000000},
        {6.4910383, 2.4182367, 0.0000000},
        {8.0910482, 6.7091894, 0.0000000},
        {13.2910795, 5.2182646, 0.0000000},
        {16.1274605, 2.0909605, 0.0000000},
        {18.7820244, 7.5091972, 0.0000000},
        {23.5456886, 8.1637497, 0.0000000},
        {24.6366081, 1.5818644, 0.0000000},
        {27.8729897, 6.5273695, 0.0000000},
        {28.0184441, 14.7092705, 0.0000000},
        {23.8366013, 12.4910660, 0.0000000},
        {22.6002312, 17.2183857, 0.0000000},
        {24.0184212, 23.7639065, 0.0000000},
        {27.9457169, 25.9093838, 0.0000000},
        {24.6366081, 27.3275795, 0.0000000},
        {17.4365616, 27.5821285, 0.0000000},
        {16.6001930, 22.8184433, 0.0000000},
        {12.7092581, 27.2548523, 0.0000000},
        {4.7455721, 27.6912193, 0.0000000},
        {2.6728315, 24.7457352, 0.0000000},
        {7.5819540, 22.9638996, 0.0000000},
        {3.5455656, 15.1820030, 0.0000000},
        {4.8546629, 11.8365135, 0.0000000},
        {8.7455969, 17.2183857, 0.0000000},
        {11.8404741, 17.6817989, 3.0000000},
        {13.5837603, 12.1783361, 3.0000000},
        {17.2156067, 16.9106121, 3.0000000}};

    auto interpolator = std::make_unique<AveragingInterpolation>(*mesh,
                                                                 samples,
                                                                 AveragingInterpolation::Method::MinAbsValue,
                                                                 Location::Faces,
                                                                 1.0,
                                                                 false,
                                                                 false,
                                                                 1);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 2;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;
    meshRefinementParameters.min_edge_size = 3.0;
    meshRefinementParameters.account_for_samples_outside = 0;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.refinement_type = 1;

    MeshRefinement meshRefinement(*mesh, std::move(interpolator), meshRefinementParameters);

    meshRefinement.Compute();

    // assert on number of nodes and edges
    ASSERT_EQ(49, mesh->GetNumNodes());
    ASSERT_EQ(84, mesh->GetNumEdges());

    // assert on newly generated edges
    ASSERT_EQ(26, mesh->GetEdge(70).first);
    ASSERT_EQ(14, mesh->GetEdge(70).second);

    ASSERT_EQ(27, mesh->GetEdge(71).first);
    ASSERT_EQ(15, mesh->GetEdge(71).second);

    ASSERT_EQ(28, mesh->GetEdge(72).first);
    ASSERT_EQ(0, mesh->GetEdge(72).second);

    ASSERT_EQ(29, mesh->GetEdge(73).first);
    ASSERT_EQ(1, mesh->GetEdge(73).second);

    ASSERT_EQ(30, mesh->GetEdge(74).first);
    ASSERT_EQ(2, mesh->GetEdge(74).second);

    ASSERT_EQ(31, mesh->GetEdge(75).first);
    ASSERT_EQ(4, mesh->GetEdge(75).second);

    ASSERT_EQ(32, mesh->GetEdge(76).first);
    ASSERT_EQ(5, mesh->GetEdge(76).second);

    ASSERT_EQ(33, mesh->GetEdge(77).first);
    ASSERT_EQ(6, mesh->GetEdge(77).second);

    ASSERT_EQ(34, mesh->GetEdge(78).first);
    ASSERT_EQ(8, mesh->GetEdge(78).second);

    ASSERT_EQ(35, mesh->GetEdge(79).first);
    ASSERT_EQ(9, mesh->GetEdge(79).second);

    ASSERT_EQ(36, mesh->GetEdge(80).first);
    ASSERT_EQ(10, mesh->GetEdge(80).second);
}

TEST(MeshRefinement, WindowOfRefinementFile)
{
    // Prepare
    auto mesh = MakeRectangularMeshForTesting(4, 4, 40.0, Projection::cartesian, {197253.0, 442281.0});

    // Sample points
    std::vector<Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/MeshRefinementTests/WindowOfRefinementFile.xyz");

    auto interpolator = std::make_unique<AveragingInterpolation>(*mesh,
                                                                 samples,
                                                                 AveragingInterpolation::Method::MinAbsValue,
                                                                 Location::Faces,
                                                                 1.0,
                                                                 false,
                                                                 false,
                                                                 1);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 4;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;
    meshRefinementParameters.min_edge_size = 3.0;
    meshRefinementParameters.account_for_samples_outside = 0;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.refinement_type = 1;
    meshRefinementParameters.smoothing_iterations = 0;

    MeshRefinement meshRefinement(*mesh,
                                  std::move(interpolator),
                                  meshRefinementParameters);

    meshRefinement.Compute();

    // total number of edges
    ASSERT_EQ(461, mesh->GetNumNodes());
    ASSERT_EQ(918, mesh->GetNumEdges());

    ASSERT_EQ(233, mesh->GetEdge(906).first);
    ASSERT_EQ(206, mesh->GetEdge(906).second);

    ASSERT_EQ(206, mesh->GetEdge(907).first);
    ASSERT_EQ(70, mesh->GetEdge(907).second);

    ASSERT_EQ(70, mesh->GetEdge(908).first);
    ASSERT_EQ(233, mesh->GetEdge(908).second);

    ASSERT_EQ(179, mesh->GetEdge(909).first);
    ASSERT_EQ(235, mesh->GetEdge(909).second);

    ASSERT_EQ(235, mesh->GetEdge(910).first);
    ASSERT_EQ(62, mesh->GetEdge(910).second);

    ASSERT_EQ(62, mesh->GetEdge(911).first);
    ASSERT_EQ(179, mesh->GetEdge(911).second);

    ASSERT_EQ(249, mesh->GetEdge(912).first);
    ASSERT_EQ(320, mesh->GetEdge(912).second);

    ASSERT_EQ(320, mesh->GetEdge(913).first);
    ASSERT_EQ(84, mesh->GetEdge(913).second);

    ASSERT_EQ(84, mesh->GetEdge(914).first);
    ASSERT_EQ(249, mesh->GetEdge(914).second);

    ASSERT_EQ(327, mesh->GetEdge(915).first);
    ASSERT_EQ(326, mesh->GetEdge(915).second);
}

TEST(MeshRefinement, MeshRefinementRefinementLevels_OnWindowOfRefinementFile_ShouldRefinemesh)
{
    // Prepare
    auto mesh = MakeRectangularMeshForTesting(4, 4, 40.0, Projection::cartesian, {197253.0, 442281.0});

    // Sample points
    std::vector<Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/MeshRefinementTests/WindowOfRefinementFile.xyz");

    auto interpolator = std::make_unique<AveragingInterpolation>(*mesh,
                                                                 samples,
                                                                 AveragingInterpolation::Method::Max,
                                                                 Location::Faces,
                                                                 1.01,
                                                                 false,
                                                                 true,
                                                                 1);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 10;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;
    meshRefinementParameters.min_edge_size = 0.5;
    meshRefinementParameters.account_for_samples_outside = 0;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.refinement_type = 2;
    meshRefinementParameters.smoothing_iterations = 0;

    MeshRefinement meshRefinement(*mesh,
                                  std::move(interpolator),
                                  meshRefinementParameters);

    meshRefinement.Compute();

    // total number of edges
    ASSERT_EQ(413, mesh->GetNumNodes());
    ASSERT_EQ(839, mesh->GetNumEdges());

    ASSERT_EQ(0, mesh->GetEdge(0).first);
    ASSERT_EQ(129, mesh->GetEdge(0).second);

    ASSERT_EQ(1, mesh->GetEdge(1).first);
    ASSERT_EQ(130, mesh->GetEdge(1).second);

    ASSERT_EQ(2, mesh->GetEdge(2).first);
    ASSERT_EQ(131, mesh->GetEdge(2).second);

    ASSERT_EQ(3, mesh->GetEdge(3).first);
    ASSERT_EQ(48, mesh->GetEdge(3).second);

    ASSERT_EQ(4, mesh->GetEdge(4).first);
    ASSERT_EQ(49, mesh->GetEdge(4).second);

    ASSERT_EQ(5, mesh->GetEdge(5).first);
    ASSERT_EQ(132, mesh->GetEdge(5).second);

    ASSERT_EQ(6, mesh->GetEdge(6).first);
    ASSERT_EQ(133, mesh->GetEdge(6).second);

    ASSERT_EQ(7, mesh->GetEdge(7).first);
    ASSERT_EQ(22, mesh->GetEdge(7).second);

    ASSERT_EQ(8, mesh->GetEdge(8).first);
    ASSERT_EQ(23, mesh->GetEdge(8).second);

    ASSERT_EQ(9, mesh->GetEdge(9).first);
    ASSERT_EQ(134, mesh->GetEdge(9).second);

    ASSERT_EQ(10, mesh->GetEdge(10).first);
    ASSERT_EQ(135, mesh->GetEdge(10).second);
}

TEST(MeshRefinement, RefineBasedOnPolygon)
{
    // Prepare
    auto mesh = MakeRectangularMeshForTesting(5, 5, 10.0, Projection::cartesian);

    std::vector<Point> point{
        {25.0, -10.0},
        {25.0, 15.0},
        {45.0, 15.0},
        {45.0, -10.0},
        {25.0, -10.0}};

    const auto polygon = Polygons(point, mesh->m_projection);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 1;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;

    MeshRefinement meshRefinement(*mesh, polygon, meshRefinementParameters);

    meshRefinement.Compute();

    // total number of edges
    ASSERT_EQ(30, mesh->GetNumNodes());
    ASSERT_EQ(52, mesh->GetNumEdges());

    ASSERT_EQ(25, mesh->GetEdge(40).first);
    ASSERT_EQ(29, mesh->GetEdge(40).second);

    ASSERT_EQ(28, mesh->GetEdge(41).first);
    ASSERT_EQ(29, mesh->GetEdge(41).second);

    ASSERT_EQ(26, mesh->GetEdge(42).first);
    ASSERT_EQ(29, mesh->GetEdge(42).second);

    ASSERT_EQ(27, mesh->GetEdge(43).first);
    ASSERT_EQ(29, mesh->GetEdge(43).second);

    ASSERT_EQ(25, mesh->GetEdge(44).first);
    ASSERT_EQ(20, mesh->GetEdge(44).second);

    ASSERT_EQ(26, mesh->GetEdge(45).first);
    ASSERT_EQ(21, mesh->GetEdge(45).second);

    ASSERT_EQ(27, mesh->GetEdge(46).first);
    ASSERT_EQ(15, mesh->GetEdge(46).second);

    ASSERT_EQ(28, mesh->GetEdge(47).first);
    ASSERT_EQ(20, mesh->GetEdge(47).second);

    ASSERT_EQ(10, mesh->GetEdge(48).first);
    ASSERT_EQ(27, mesh->GetEdge(48).second);
}

TEST(MeshRefinement, RefineBasedOnPolygonThreeByThree)
{
    // Prepare
    auto mesh = MakeRectangularMeshForTesting(4, 4, 10.0, Projection::cartesian);

    std::vector<Point> point{
        {9.09836065573771, 34.016393442623},
        {7.18032786885247, 7.75409836065574},
        {34.6229508196721, 6.5},
        {34.4194409808304, 26.6983515050386},
        {34.327868852459, 35.7868852459016},
        {29.0521194370216, 35.4840476661635},
        {9.90983606557378, 34.3852459016394},
        {9.09836065573771, 34.016393442623}};

    const auto polygon = Polygons(point, mesh->m_projection);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 2;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;

    MeshRefinement meshRefinement(*mesh, polygon, meshRefinementParameters);
    meshRefinement.Compute();

    // assert on number of nodes and edges
    ASSERT_EQ(48, mesh->GetNumNodes());
    ASSERT_EQ(96, mesh->GetNumEdges());
    ASSERT_EQ(49, mesh->GetNumFaces());
}

TEST(MeshRefinement, FourByFourWithFourSamplesSpherical)
{

    auto mesh = MakeRectangularMeshForTesting(4, 4, 0.0033, Projection::spherical, {41.1, 41.1});

    // sample points
    std::vector<Sample> samples{
        {41.1050110, 41.1049728, 1.0},
        {41.1084785, 41.1048775, 1.0},
        {41.1085625, 41.1083946, 1.0},
        {41.1052971, 41.1083336, 1.0}};

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
    meshRefinementParameters.min_edge_size = 0.00165;
    meshRefinementParameters.account_for_samples_outside = 0;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.refinement_type = 2;
    meshRefinementParameters.smoothing_iterations = 0;

    MeshRefinement meshRefinement(*mesh,
                                  std::move(interpolator),
                                  meshRefinementParameters);
    meshRefinement.Compute();

    ASSERT_EQ(60, mesh->GetNumEdges());
    ASSERT_EQ(32, mesh->GetNumNodes());

    // sides of the refined part
    ASSERT_EQ(5, mesh->GetEdge(5).first);
    ASSERT_EQ(16, mesh->GetEdge(5).second);

    ASSERT_EQ(16, mesh->GetEdge(40).first);
    ASSERT_EQ(9, mesh->GetEdge(40).second);

    ASSERT_EQ(9, mesh->GetEdge(9).first);
    ASSERT_EQ(19, mesh->GetEdge(9).second);

    ASSERT_EQ(19, mesh->GetEdge(43).first);
    ASSERT_EQ(13, mesh->GetEdge(43).second);

    ASSERT_EQ(6, mesh->GetEdge(16).first);
    ASSERT_EQ(22, mesh->GetEdge(16).second);

    ASSERT_EQ(22, mesh->GetEdge(46).first);
    ASSERT_EQ(5, mesh->GetEdge(46).second);

    ASSERT_EQ(7, mesh->GetEdge(17).first);
    ASSERT_EQ(23, mesh->GetEdge(17).second);

    ASSERT_EQ(23, mesh->GetEdge(47).first);
    ASSERT_EQ(6, mesh->GetEdge(47).second);
}

TEST(MeshRefinement, RefinementFileBasedOnLevels_OnSphericalMesh_ShouldRefine)
{
    // Prepare
    auto mesh = MakeRectangularMeshForTesting(6, 6, 0.0033, Projection::spherical, {41.1, 41.1});

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 1;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;
    meshRefinementParameters.min_edge_size = 0.00165;
    meshRefinementParameters.account_for_samples_outside = 0;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.refinement_type = 2;

    // Execute
    const auto polygon = Polygons();
    MeshRefinement meshRefinement(*mesh, polygon, meshRefinementParameters);
    meshRefinement.Compute();

    // Assert, we passed from 36 to 49 nodes
    ASSERT_EQ(121, mesh->GetNumNodes());

    constexpr double tolerance = 1e-6;
    ASSERT_NEAR(41.100000000000001, mesh->Node(0).x, tolerance);
    ASSERT_NEAR(41.103300000000004, mesh->Node(6).x, tolerance);
    ASSERT_NEAR(41.106600000000000, mesh->Node(12).x, tolerance);
    ASSERT_NEAR(41.109900000000003, mesh->Node(18).x, tolerance);
    ASSERT_NEAR(41.113199999999999, mesh->Node(24).x, tolerance);

    ASSERT_NEAR(41.100000000000001, mesh->Node(0).y, tolerance);
    ASSERT_NEAR(41.100000000000001, mesh->Node(6).y, tolerance);
    ASSERT_NEAR(41.100000000000001, mesh->Node(12).y, tolerance);
    ASSERT_NEAR(41.100000000000001, mesh->Node(18).y, tolerance);
    ASSERT_NEAR(41.100000000000001, mesh->Node(24).y, tolerance);
}

TEST(MeshRefinement, RefineCurvilinearGrid)
{
    auto mesh = MakeCurvilinearGridForTesting();

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 1;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;

    const auto polygon = Polygons();
    MeshRefinement meshRefinement(*mesh, polygon, meshRefinementParameters);
    meshRefinement.Compute();

    mesh->ComputeEdgesLengths();

    // if the circumcenters are wrongly computed, some edges will be smaller than half cell size
    for (size_t i = 0; i < mesh->GetNumEdges(); ++i)
    {
        ASSERT_GT(mesh->m_edgeLengths[i], 0.4);
    }
}

TEST(MeshRefinement, RefineElongatedFaces)
{
    // Prepare
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/MeshRefinementTests/CurvilinearEnlonged.nc");

    std::vector<Point> point{
        {2018.73356016594, 1165.26385465465},
        {1694.83823023783, 1131.75744121381},
        {1708.79923583818, 640.330044081495},
        {2363.7176353495, 645.640193266722},
        {2741.91365026406, 648.706647441705},
        {2722.36824242357, 1173.64045801486},
        {2038.27896800643, 1165.26385465465},
        {2018.73356016594, 1165.26385465465}};

    const auto polygon = Polygons(point, mesh->m_projection);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 2;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;

    MeshRefinement meshRefinement(*mesh, polygon, meshRefinementParameters);

    // Execute
    meshRefinement.Compute();

    // Assert circumcenters are correctly computed
    constexpr double tolerance = 1e-6;

    // Compare the x-location of the circumcentre.
    ASSERT_NEAR(1673.0860169014584, mesh->m_facesCircumcenters[0].x, tolerance);
    ASSERT_NEAR(1660.6851354957175, mesh->m_facesCircumcenters[1].x, tolerance);
    ASSERT_NEAR(1660.5667704694627, mesh->m_facesCircumcenters[2].x, tolerance);
    ASSERT_NEAR(1672.0912775041329, mesh->m_facesCircumcenters[3].x, tolerance);
    ASSERT_NEAR(1659.9354211078053, mesh->m_facesCircumcenters[4].x, tolerance);
    ASSERT_NEAR(1659.8248648846848, mesh->m_facesCircumcenters[5].x, tolerance);
    ASSERT_NEAR(1671.1074693287451, mesh->m_facesCircumcenters[6].x, tolerance);
    ASSERT_NEAR(1659.1859707978906, mesh->m_facesCircumcenters[7].x, tolerance);
    ASSERT_NEAR(1659.0828479935451, mesh->m_facesCircumcenters[8].x, tolerance);
    ASSERT_NEAR(1670.1135380487042, mesh->m_facesCircumcenters[9].x, tolerance);
    ASSERT_NEAR(1658.4375379474418, mesh->m_facesCircumcenters[10].x, tolerance);

    ASSERT_NEAR(645.10565853980427, mesh->m_facesCircumcenters[0].y, tolerance);
    ASSERT_NEAR(646.20173898461292, mesh->m_facesCircumcenters[1].y, tolerance);
    ASSERT_NEAR(654.66978646149676, mesh->m_facesCircumcenters[2].y, tolerance);
    ASSERT_NEAR(662.96936808193914, mesh->m_facesCircumcenters[3].y, tolerance);
    ASSERT_NEAR(664.17198930960728, mesh->m_facesCircumcenters[4].y, tolerance);
    ASSERT_NEAR(672.49588595953537, mesh->m_facesCircumcenters[5].y, tolerance);
    ASSERT_NEAR(680.82566423059006, mesh->m_facesCircumcenters[6].y, tolerance);
    ASSERT_NEAR(682.12720924968505, mesh->m_facesCircumcenters[7].y, tolerance);
    ASSERT_NEAR(690.31918193748425, mesh->m_facesCircumcenters[8].y, tolerance);
    ASSERT_NEAR(698.66471917887850, mesh->m_facesCircumcenters[9].y, tolerance);
    ASSERT_NEAR(700.06356972686194, mesh->m_facesCircumcenters[10].y, tolerance);
}

TEST(MeshRefinement, BilinearInterpolationWithGriddedSamplesOnLandShouldNotRefine)
{
    // Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 10.0, Projection::cartesian);

    std::vector<float> values{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    Point origin{-5.0, -5.0};
    auto interpolator = std::make_unique<BilinearInterpolationOnGriddedSamples<float>>(*mesh, 2, 2, origin, 10.0, values);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 1;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;
    meshRefinementParameters.min_edge_size = 5.0;
    meshRefinementParameters.account_for_samples_outside = 1;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.refinement_type = 1;

    MeshRefinement meshRefinement(*mesh, std::move(interpolator), meshRefinementParameters, true);

    // Execute
    meshRefinement.Compute();

    // Assert: all bathy values are positive and we are in land, so nothing gets refined
    ASSERT_EQ(4, mesh->GetNumEdges());
}

TEST(MeshRefinement, BilinearInterpolationWithGriddedSamplesOnLandAndSeaShouldRefine)
{
    // Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 10.0, Projection::cartesian);

    std::vector<float> values{-1.0, -2.0, 3.0, -4.0, -5.0, 6.0, 7.0, 8.0, 9.0};
    Point origin{-5.0, -5.0};
    auto interpolator = std::make_unique<BilinearInterpolationOnGriddedSamples<float>>(*mesh, 3, 3, origin, 10.0, values);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 1;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;
    meshRefinementParameters.min_edge_size = 5.0;
    meshRefinementParameters.account_for_samples_outside = 1;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.refinement_type = 1;

    MeshRefinement meshRefinement(*mesh, std::move(interpolator), meshRefinementParameters, true);

    // Execute
    meshRefinement.Compute();

    // Assert: all depth values are positive and we are in land, so nothing gets refined
    ASSERT_EQ(12, mesh->GetNumEdges());
}

TEST(MeshRefinement, BilinearInterpolationWithAllGriddedSamplesOnSeaShouldRefine)
{
    // Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 10.0, Projection::cartesian);

    std::vector<float> values{-1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0, -8.0, -9.0};
    Point origin{-5.0, -5.0};
    auto interpolator = std::make_unique<BilinearInterpolationOnGriddedSamples<float>>(*mesh, 2, 2, origin, 10.0, values);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 1;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;
    meshRefinementParameters.min_edge_size = 5.0;
    meshRefinementParameters.account_for_samples_outside = 1;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.refinement_type = 1;

    MeshRefinement meshRefinement(*mesh,
                                  std::move(interpolator),
                                  meshRefinementParameters, true);

    // Execute
    meshRefinement.Compute();

    // Assert: nothing gets refined
    ASSERT_EQ(4, mesh->GetNumEdges());
}

class RidgeRefinementTestCases : public testing::TestWithParam<std::tuple<FunctionTestCase, UInt, UInt>>
{
public:
    [[nodiscard]] static std::vector<std::tuple<FunctionTestCase, UInt, UInt>> GetData()
    {
        return std::vector{
            std::make_tuple<FunctionTestCase, UInt, UInt>(FunctionTestCase::GaussianBump, 1165, 2344),
            std::make_tuple<FunctionTestCase, UInt, UInt>(FunctionTestCase::GaussianWave, 5297, 10784),
            std::make_tuple<FunctionTestCase, UInt, UInt>(FunctionTestCase::RidgeXDirection, 2618, 5694),
            std::make_tuple<FunctionTestCase, UInt, UInt>(FunctionTestCase::ArctanFunction, 2309, 5028)};
    }
};

TEST_P(RidgeRefinementTestCases, expectedResults)
{
    // Prepare
    const auto [testCase, numNodes, numEdges] = GetParam();

    UInt nx = 41;
    UInt ny = 21;

    double deltaX = 10.0;
    double deltaY = 10.0;

    double dimX = (nx - 1) * deltaX;
    double dimY = (ny - 1) * deltaY;

    auto mesh = MakeRectangularMeshForTesting(nx, ny, dimX, dimY, Projection::cartesian);

    UInt superSample = 2;
    UInt sampleNx = (nx - 1) * superSample + 1;
    UInt sampleNy = (ny - 1) * superSample + 1;

    const double sampleDeltaX = deltaX / static_cast<double>(superSample);
    const double sampleDeltaY = deltaY / static_cast<double>(superSample);

    const auto sampleData = generateSampleData(testCase, sampleNx, sampleNy, sampleDeltaX, sampleDeltaY);

    auto samplesHessian = SamplesHessianCalculator::ComputeSamplesHessian(sampleData, mesh->m_projection, 0, sampleNx, sampleNy);

    auto interpolator = std::make_unique<AveragingInterpolation>(*mesh,
                                                                 samplesHessian,
                                                                 AveragingInterpolation::Method::Max,
                                                                 Location::Faces,
                                                                 1.0,
                                                                 false,
                                                                 false,
                                                                 1);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 3;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;
    meshRefinementParameters.min_edge_size = 2.0;
    meshRefinementParameters.account_for_samples_outside = 0;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.refinement_type = 3;
    meshRefinementParameters.smoothing_iterations = 0;

    MeshRefinement meshRefinement(*mesh,
                                  std::move(interpolator),
                                  meshRefinementParameters, false);

    // Execute
    meshRefinement.Compute();

    // Assert
    ASSERT_EQ(numNodes, mesh->GetNumNodes());
    ASSERT_EQ(numEdges, mesh->GetNumEdges());
}

INSTANTIATE_TEST_SUITE_P(RidgeRefinementTestCases,
                         RidgeRefinementTestCases,
                         ::testing::ValuesIn(RidgeRefinementTestCases::GetData()));

TEST(MeshRefinement, CasulliRefinement)
{
    constexpr double tolerance = 1.0e-12;

    auto curviMesh = MakeCurvilinearGrid(0.0, 0.0, 10.0, 10.0, 3, 3);
    Mesh2D mesh(curviMesh->Edges(), curviMesh->Nodes(), Projection::cartesian);

    // Expected values were obtained from a mesh refined using the Casulli refinement algorithm
    std::vector<meshkernel::Point> expectedPoints{{0, 0},
                                                  {20, 0},
                                                  {0, 20},
                                                  {20, 20},
                                                  {2.5, 2.5},
                                                  {7.5, 2.5},
                                                  {7.5, 7.5},
                                                  {2.5, 7.5},
                                                  {12.5, 2.5},
                                                  {17.5, 2.5},
                                                  {17.5, 7.5},
                                                  {12.5, 7.5},
                                                  {2.5, 12.5},
                                                  {7.5, 12.5},
                                                  {7.5, 17.5},
                                                  {2.5, 17.5},
                                                  {12.5, 12.5},
                                                  {17.5, 12.5},
                                                  {17.5, 17.5},
                                                  {12.5, 17.5},
                                                  {2.5, 0},
                                                  {7.5, 0},
                                                  {12.5, 0},
                                                  {17.5, 0},
                                                  {2.5, 20},
                                                  {7.5, 20},
                                                  {12.5, 20},
                                                  {17.5, 20},
                                                  {0, 2.5},
                                                  {0, 7.5},
                                                  {20, 2.5},
                                                  {20, 7.5},
                                                  {0, 12.5},
                                                  {0, 17.5},
                                                  {20, 12.5},
                                                  {20, 17.5}};

    std::vector<int> expectedEdgesStart = {20, 4, 20, 21, 22, 8, 22, 23, 7, 12, 7, 6, 11, 16, 11, 10, 15, 24, 15, 14,
                                           19, 26, 19, 18, 4, 28, 4, 7, 8, 5, 8, 11, 30, 9, 30, 31, 12, 32, 12, 15,
                                           16, 13, 16, 19, 34, 17, 34, 35, 0, 0, 21, 1, 1, 32, 31, 2, 2, 25, 3, 3};

    std::vector<int> expectedEdgesEnd = {21, 5, 4, 5, 23, 9, 8, 9, 6, 13, 12, 13, 10, 17, 16, 17, 14, 25, 24, 25,
                                         18, 27, 26, 27, 7, 29, 28, 29, 11, 6, 5, 6, 31, 10, 9, 10, 15, 33, 32, 33,
                                         19, 14, 13, 14, 35, 18, 17, 18, 20, 28, 22, 23, 30, 29, 34, 24, 33, 26, 27, 35};

    CasulliRefinement meshRefinement;

    meshRefinement.Compute(mesh);

    std::vector<meshkernel::UInt> validNodeMap(mesh.GetValidNodeMapping());
    std::vector<meshkernel::UInt> validEdgeMap(mesh.GetValidEdgeMapping());

    ASSERT_EQ(expectedPoints.size(), mesh.GetNumValidNodes());

    for (size_t i = 0; i < expectedPoints.size(); ++i)
    {
        EXPECT_NEAR(expectedPoints[i].x, mesh.Node(validNodeMap[i]).x, tolerance);
        EXPECT_NEAR(expectedPoints[i].y, mesh.Node(validNodeMap[i]).y, tolerance);
    }

    ASSERT_EQ(expectedEdgesStart.size(), mesh.GetNumValidEdges());

    for (size_t i = 0; i < expectedEdgesStart.size(); ++i)
    {
        EXPECT_EQ(mesh.GetEdge(validEdgeMap[i]).first, validNodeMap[expectedEdgesStart[i]]);
        EXPECT_EQ(mesh.GetEdge(validEdgeMap[i]).second, validNodeMap[expectedEdgesEnd[i]]);
    }
}

void loadCasulliRefinedMeshData(std::vector<Point>& expectedPoints,
                                std::vector<meshkernel::UInt>& expectedEdgeStart,
                                std::vector<meshkernel::UInt>& expectedEdgeEnd)
{

    const std::string fileName = TEST_FOLDER + "/data/CasulliRefinement/casulli_refinement_patch_with_hole.txt";

    std::ifstream asciiFile;
    asciiFile.open(fileName.c_str());

    for (size_t i = 0; i < expectedPoints.size(); ++i)
    {
        asciiFile >> expectedPoints[i].x;
    }

    for (size_t i = 0; i < expectedPoints.size(); ++i)
    {
        asciiFile >> expectedPoints[i].y;
    }

    for (size_t i = 0; i < expectedEdgeStart.size(); ++i)
    {
        asciiFile >> expectedEdgeStart[i];
    }

    for (size_t i = 0; i < expectedEdgeEnd.size(); ++i)
    {
        asciiFile >> expectedEdgeEnd[i];
    }

    asciiFile.close();
}

TEST(MeshRefinement, CasulliPatchRefinement)
{
    const size_t ExpectedNumberOfPoints = 184;
    const size_t ExpectedNumberOfEdges = 360;

    auto curviMesh = MakeCurvilinearGrid(0.0, 0.0, 20.0, 20.0, 11, 11);
    Mesh2D mesh(curviMesh->Edges(), curviMesh->Nodes(), Projection::cartesian);

    std::vector<Point> patch{{45.0, 45.0},
                             {155.0, 45.0},
                             {155.0, 155.0},
                             {45.0, 155.0},
                             {45.0, 45.0},
                             {constants::missing::innerOuterSeparator, constants::missing::innerOuterSeparator},
                             {65.0, 65.0},
                             {115.0, 65.0},
                             {115.0, 115.0},
                             {65.0, 115.0},
                             {65.0, 65.0}};

    std::vector<Point> expectedPoints(ExpectedNumberOfPoints);
    std::vector<meshkernel::UInt> expectedEdgeStart(ExpectedNumberOfEdges);
    std::vector<meshkernel::UInt> expectedEdgeEnd(ExpectedNumberOfEdges);

    Polygons polygon(patch, Projection::cartesian);

    CasulliRefinement meshRefinement;

    meshRefinement.Compute(mesh, polygon);

    constexpr double tolerance = 1.0e-12;

    loadCasulliRefinedMeshData(expectedPoints, expectedEdgeStart, expectedEdgeEnd);

    std::vector<meshkernel::UInt> validNodeMap(mesh.GetValidNodeMapping());
    std::vector<meshkernel::UInt> validEdgeMap(mesh.GetValidEdgeMapping());

    ASSERT_EQ(ExpectedNumberOfPoints, mesh.GetNumValidNodes());
    ASSERT_EQ(ExpectedNumberOfEdges, mesh.GetNumValidEdges());

    for (size_t i = 0; i < expectedPoints.size(); ++i)
    {
        // Map the index i from the nodes array containing only valid points to an array with that may contain in-valid points
        EXPECT_NEAR(expectedPoints[i].x, mesh.Node(validNodeMap[i]).x, tolerance);
        EXPECT_NEAR(expectedPoints[i].y, mesh.Node(validNodeMap[i]).y, tolerance);
    }

    for (size_t i = 0; i < expectedEdgeStart.size(); ++i)
    {
        // Map the index i from the edges array containing only valid edges to an array with that may contain in-valid edges
        EXPECT_EQ(mesh.GetEdge(validEdgeMap[i]).first, validNodeMap[expectedEdgeStart[i]]);
        EXPECT_EQ(mesh.GetEdge(validEdgeMap[i]).second, validNodeMap[expectedEdgeEnd[i]]);
    }
}
