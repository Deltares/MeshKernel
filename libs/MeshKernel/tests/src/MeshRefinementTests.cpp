#include "MeshKernel/BilinearInterpolationOnGriddedSamples.hpp"

#include <gtest/gtest.h>

#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/MeshRefinement.hpp>
#include <MeshKernel/Parameters.hpp>
#include <MeshKernel/Polygons.hpp>
#include <TestUtils/Definitions.hpp>
#include <TestUtils/MakeMeshes.hpp>
#include <TestUtils/SampleFileReader.hpp>

using namespace meshkernel;

TEST(MeshRefinement, FourByFourWithFourSamples)
{
    auto mesh = MakeRectangularMeshForTesting(5, 5, 10.0, Projection::cartesian);

    // sample points
    std::vector<Sample> samples{
        {14.7153645, 14.5698833, 1.0},
        {24.7033062, 14.4729137, 1.0},
        {15.5396099, 24.2669525, 1.0},
        {23.8305721, 23.9275551, 1.0}};

    const auto interpolator = std::make_shared<AveragingInterpolation>(*mesh,
                                                                       samples,
                                                                       AveragingInterpolation::Method::MinAbsValue,
                                                                       Mesh::Location::Faces,
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

    MeshRefinement meshRefinement(mesh, interpolator, meshRefinementParameters);

    meshRefinement.Compute();

    // 3 Validation edges connecting hanging nodes

    // bottom side
    ASSERT_EQ(5, mesh->m_edges[73].first);
    ASSERT_EQ(25, mesh->m_edges[73].second);

    ASSERT_EQ(10, mesh->m_edges[72].first);
    ASSERT_EQ(25, mesh->m_edges[72].second);

    ASSERT_EQ(10, mesh->m_edges[77].first);
    ASSERT_EQ(28, mesh->m_edges[77].second);

    ASSERT_EQ(15, mesh->m_edges[76].first);
    ASSERT_EQ(28, mesh->m_edges[76].second);

    // right side
    ASSERT_EQ(21, mesh->m_edges[81].first);
    ASSERT_EQ(35, mesh->m_edges[81].second);

    ASSERT_EQ(22, mesh->m_edges[80].first);
    ASSERT_EQ(35, mesh->m_edges[80].second);

    ASSERT_EQ(22, mesh->m_edges[83].first);
    ASSERT_EQ(36, mesh->m_edges[83].second);

    ASSERT_EQ(23, mesh->m_edges[82].first);
    ASSERT_EQ(36, mesh->m_edges[82].second);

    // upper side
    ASSERT_EQ(19, mesh->m_edges[79].first);
    ASSERT_EQ(30, mesh->m_edges[79].second);

    ASSERT_EQ(14, mesh->m_edges[78].first);
    ASSERT_EQ(30, mesh->m_edges[78].second);

    ASSERT_EQ(14, mesh->m_edges[75].first);
    ASSERT_EQ(27, mesh->m_edges[75].second);

    ASSERT_EQ(9, mesh->m_edges[74].first);
    ASSERT_EQ(27, mesh->m_edges[74].second);

    // left side
    ASSERT_EQ(3, mesh->m_edges[71].first);
    ASSERT_EQ(32, mesh->m_edges[71].second);

    ASSERT_EQ(2, mesh->m_edges[70].first);
    ASSERT_EQ(32, mesh->m_edges[70].second);

    ASSERT_EQ(2, mesh->m_edges[69].first);
    ASSERT_EQ(31, mesh->m_edges[69].second);

    ASSERT_EQ(1, mesh->m_edges[68].first);
    ASSERT_EQ(31, mesh->m_edges[68].second);

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

    const auto interpolator = std::make_shared<AveragingInterpolation>(*mesh,
                                                                       samples,
                                                                       AveragingInterpolation::Method::MinAbsValue,
                                                                       Mesh::Location::Faces,
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

    MeshRefinement meshRefinement(mesh, interpolator, meshRefinementParameters);
    meshRefinement.Compute();

    // Assert number of edges and nodes
    ASSERT_EQ(60, mesh->GetNumEdges());
    ASSERT_EQ(31, mesh->GetNumNodes());

    // Assert edges
    ASSERT_EQ(0, mesh->m_edges[0].first);
    ASSERT_EQ(26, mesh->m_edges[0].second);

    ASSERT_EQ(1, mesh->m_edges[1].first);
    ASSERT_EQ(17, mesh->m_edges[1].second);

    ASSERT_EQ(2, mesh->m_edges[2].first);
    ASSERT_EQ(6, mesh->m_edges[2].second);

    ASSERT_EQ(3, mesh->m_edges[3].first);
    ASSERT_EQ(7, mesh->m_edges[3].second);

    ASSERT_EQ(4, mesh->m_edges[4].first);
    ASSERT_EQ(8, mesh->m_edges[4].second);

    ASSERT_EQ(5, mesh->m_edges[5].first);
    ASSERT_EQ(9, mesh->m_edges[5].second);

    ASSERT_EQ(6, mesh->m_edges[6].first);
    ASSERT_EQ(10, mesh->m_edges[6].second);

    ASSERT_EQ(7, mesh->m_edges[7].first);
    ASSERT_EQ(11, mesh->m_edges[7].second);

    ASSERT_EQ(8, mesh->m_edges[8].first);
    ASSERT_EQ(12, mesh->m_edges[8].second);

    ASSERT_EQ(9, mesh->m_edges[9].first);
    ASSERT_EQ(13, mesh->m_edges[9].second);

    ASSERT_EQ(10, mesh->m_edges[10].first);
    ASSERT_EQ(14, mesh->m_edges[10].second);

    ASSERT_EQ(11, mesh->m_edges[11].first);
    ASSERT_EQ(15, mesh->m_edges[11].second);

    ASSERT_EQ(1, mesh->m_edges[12].first);
    ASSERT_EQ(18, mesh->m_edges[12].second);

    ASSERT_EQ(2, mesh->m_edges[13].first);
    ASSERT_EQ(1, mesh->m_edges[13].second);

    ASSERT_EQ(3, mesh->m_edges[14].first);
    ASSERT_EQ(2, mesh->m_edges[14].second);

    ASSERT_EQ(5, mesh->m_edges[15].first);
    ASSERT_EQ(19, mesh->m_edges[15].second);

    ASSERT_EQ(6, mesh->m_edges[16].first);
    ASSERT_EQ(5, mesh->m_edges[16].second);

    ASSERT_EQ(7, mesh->m_edges[17].first);
    ASSERT_EQ(6, mesh->m_edges[17].second);

    ASSERT_EQ(9, mesh->m_edges[18].first);
    ASSERT_EQ(8, mesh->m_edges[18].second);

    ASSERT_EQ(10, mesh->m_edges[19].first);
    ASSERT_EQ(9, mesh->m_edges[19].second);

    ASSERT_EQ(11, mesh->m_edges[20].first);
    ASSERT_EQ(10, mesh->m_edges[20].second);
}

TEST(MeshRefinement, SmallTriangualMeshTwoSamples)
{
    // Prepare
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/SmallTriangularGrid_net.nc");

    // sample points
    std::vector<Sample> samples{
        {359.8657532, 350.3144836, 1.0},
        {387.5152588, 299.2614746, 1.0}};

    const auto interpolator = std::make_shared<AveragingInterpolation>(*mesh,
                                                                       samples,
                                                                       AveragingInterpolation::Method::MinAbsValue,
                                                                       Mesh::Location::Faces,
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

    MeshRefinement meshRefinement(mesh, interpolator, meshRefinementParameters);

    meshRefinement.Compute();

    // edges connecting hanging nodes
    ASSERT_EQ(10, mesh->m_edges[32].first);
    ASSERT_EQ(2, mesh->m_edges[32].second);

    ASSERT_EQ(14, mesh->m_edges[33].first);
    ASSERT_EQ(4, mesh->m_edges[33].second);

    ASSERT_EQ(13, mesh->m_edges[34].first);
    ASSERT_EQ(5, mesh->m_edges[34].second);

    ASSERT_EQ(11, mesh->m_edges[31].first);
    ASSERT_EQ(1, mesh->m_edges[31].second);

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

    Polygons polygon(point, mesh->m_projection);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 1;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;

    MeshRefinement meshRefinement(mesh, polygon, meshRefinementParameters);
    meshRefinement.Compute();

    // total number of edges
    ASSERT_EQ(15, mesh->GetNumNodes());
    ASSERT_EQ(33, mesh->GetNumEdges());

    // assert on newly generated edges
    ASSERT_EQ(10, mesh->m_edges[20].first);
    ASSERT_EQ(11, mesh->m_edges[20].second);

    ASSERT_EQ(11, mesh->m_edges[21].first);
    ASSERT_EQ(12, mesh->m_edges[21].second);

    ASSERT_EQ(12, mesh->m_edges[22].first);
    ASSERT_EQ(10, mesh->m_edges[22].second);

    ASSERT_EQ(14, mesh->m_edges[23].first);
    ASSERT_EQ(13, mesh->m_edges[23].second);

    ASSERT_EQ(13, mesh->m_edges[24].first);
    ASSERT_EQ(11, mesh->m_edges[24].second);

    ASSERT_EQ(11, mesh->m_edges[25].first);
    ASSERT_EQ(14, mesh->m_edges[25].second);

    ASSERT_EQ(10, mesh->m_edges[26].first);
    ASSERT_EQ(4, mesh->m_edges[26].second);
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

    const auto interpolator = std::make_shared<AveragingInterpolation>(*mesh,
                                                                       samples,
                                                                       AveragingInterpolation::Method::MinAbsValue,
                                                                       Mesh::Location::Faces,
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

    MeshRefinement meshRefinement(mesh, interpolator, meshRefinementParameters);

    meshRefinement.Compute();

    // assert on number of nodes and edges
    ASSERT_EQ(49, mesh->GetNumNodes());
    ASSERT_EQ(84, mesh->GetNumEdges());

    // assert on newly generated edges
    ASSERT_EQ(26, mesh->m_edges[70].first);
    ASSERT_EQ(14, mesh->m_edges[70].second);

    ASSERT_EQ(27, mesh->m_edges[71].first);
    ASSERT_EQ(15, mesh->m_edges[71].second);

    ASSERT_EQ(28, mesh->m_edges[72].first);
    ASSERT_EQ(0, mesh->m_edges[72].second);

    ASSERT_EQ(29, mesh->m_edges[73].first);
    ASSERT_EQ(1, mesh->m_edges[73].second);

    ASSERT_EQ(30, mesh->m_edges[74].first);
    ASSERT_EQ(2, mesh->m_edges[74].second);

    ASSERT_EQ(31, mesh->m_edges[75].first);
    ASSERT_EQ(4, mesh->m_edges[75].second);

    ASSERT_EQ(32, mesh->m_edges[76].first);
    ASSERT_EQ(5, mesh->m_edges[76].second);

    ASSERT_EQ(33, mesh->m_edges[77].first);
    ASSERT_EQ(6, mesh->m_edges[77].second);

    ASSERT_EQ(34, mesh->m_edges[78].first);
    ASSERT_EQ(8, mesh->m_edges[78].second);

    ASSERT_EQ(35, mesh->m_edges[79].first);
    ASSERT_EQ(9, mesh->m_edges[79].second);

    ASSERT_EQ(36, mesh->m_edges[80].first);
    ASSERT_EQ(10, mesh->m_edges[80].second);
}

TEST(MeshRefinement, WindowOfRefinementFile)
{
    // Prepare
    auto mesh = MakeRectangularMeshForTesting(4, 4, 40.0, Projection::cartesian, {197253.0, 442281.0});

    // Sample points
    std::vector<Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/MeshRefinementTests/WindowOfRefinementFile.xyz");

    const auto interpolator = std::make_shared<AveragingInterpolation>(*mesh,
                                                                       samples,
                                                                       AveragingInterpolation::Method::MinAbsValue,
                                                                       Mesh::Location::Faces,
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

    MeshRefinement meshRefinement(mesh, interpolator, meshRefinementParameters);

    meshRefinement.Compute();

    // total number of edges
    ASSERT_EQ(461, mesh->GetNumNodes());
    ASSERT_EQ(918, mesh->GetNumEdges());

    ASSERT_EQ(233, mesh->m_edges[906].first);
    ASSERT_EQ(206, mesh->m_edges[906].second);

    ASSERT_EQ(206, mesh->m_edges[907].first);
    ASSERT_EQ(70, mesh->m_edges[907].second);

    ASSERT_EQ(70, mesh->m_edges[908].first);
    ASSERT_EQ(233, mesh->m_edges[908].second);

    ASSERT_EQ(179, mesh->m_edges[909].first);
    ASSERT_EQ(235, mesh->m_edges[909].second);

    ASSERT_EQ(235, mesh->m_edges[910].first);
    ASSERT_EQ(62, mesh->m_edges[910].second);

    ASSERT_EQ(62, mesh->m_edges[911].first);
    ASSERT_EQ(179, mesh->m_edges[911].second);

    ASSERT_EQ(249, mesh->m_edges[912].first);
    ASSERT_EQ(320, mesh->m_edges[912].second);

    ASSERT_EQ(320, mesh->m_edges[913].first);
    ASSERT_EQ(84, mesh->m_edges[913].second);

    ASSERT_EQ(84, mesh->m_edges[914].first);
    ASSERT_EQ(249, mesh->m_edges[914].second);

    ASSERT_EQ(327, mesh->m_edges[915].first);
    ASSERT_EQ(326, mesh->m_edges[915].second);
}

TEST(MeshRefinement, WindowOfRefinementFileBasedOnLevels)
{
    // Prepare
    auto mesh = MakeRectangularMeshForTesting(4, 4, 40.0, Projection::cartesian, {197253.0, 442281.0});

    // Sample points
    std::vector<Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/MeshRefinementTests/WindowOfRefinementFile.xyz");

    const auto interpolator = std::make_shared<AveragingInterpolation>(*mesh,
                                                                       samples,
                                                                       AveragingInterpolation::Method::Max,
                                                                       Mesh::Location::Faces,
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

    MeshRefinement meshRefinement(mesh, interpolator, meshRefinementParameters);

    meshRefinement.Compute();

    // total number of edges
    ASSERT_EQ(413, mesh->GetNumNodes());
    ASSERT_EQ(839, mesh->GetNumEdges());

    ASSERT_EQ(0, mesh->m_edges[0].first);
    ASSERT_EQ(129, mesh->m_edges[0].second);

    ASSERT_EQ(1, mesh->m_edges[1].first);
    ASSERT_EQ(130, mesh->m_edges[1].second);

    ASSERT_EQ(2, mesh->m_edges[2].first);
    ASSERT_EQ(131, mesh->m_edges[2].second);

    ASSERT_EQ(3, mesh->m_edges[3].first);
    ASSERT_EQ(48, mesh->m_edges[3].second);

    ASSERT_EQ(4, mesh->m_edges[4].first);
    ASSERT_EQ(49, mesh->m_edges[4].second);

    ASSERT_EQ(5, mesh->m_edges[5].first);
    ASSERT_EQ(132, mesh->m_edges[5].second);

    ASSERT_EQ(6, mesh->m_edges[6].first);
    ASSERT_EQ(133, mesh->m_edges[6].second);

    ASSERT_EQ(7, mesh->m_edges[7].first);
    ASSERT_EQ(22, mesh->m_edges[7].second);

    ASSERT_EQ(8, mesh->m_edges[8].first);
    ASSERT_EQ(23, mesh->m_edges[8].second);

    ASSERT_EQ(9, mesh->m_edges[9].first);
    ASSERT_EQ(134, mesh->m_edges[9].second);

    ASSERT_EQ(10, mesh->m_edges[10].first);
    ASSERT_EQ(135, mesh->m_edges[10].second);
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

    Polygons polygon(point, mesh->m_projection);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 1;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;

    MeshRefinement meshRefinement(mesh, polygon, meshRefinementParameters);

    meshRefinement.Compute();

    // total number of edges
    ASSERT_EQ(30, mesh->GetNumNodes());
    ASSERT_EQ(52, mesh->GetNumEdges());

    ASSERT_EQ(25, mesh->m_edges[40].first);
    ASSERT_EQ(29, mesh->m_edges[40].second);

    ASSERT_EQ(28, mesh->m_edges[41].first);
    ASSERT_EQ(29, mesh->m_edges[41].second);

    ASSERT_EQ(26, mesh->m_edges[42].first);
    ASSERT_EQ(29, mesh->m_edges[42].second);

    ASSERT_EQ(27, mesh->m_edges[43].first);
    ASSERT_EQ(29, mesh->m_edges[43].second);

    ASSERT_EQ(25, mesh->m_edges[44].first);
    ASSERT_EQ(20, mesh->m_edges[44].second);

    ASSERT_EQ(26, mesh->m_edges[45].first);
    ASSERT_EQ(21, mesh->m_edges[45].second);

    ASSERT_EQ(27, mesh->m_edges[46].first);
    ASSERT_EQ(15, mesh->m_edges[46].second);

    ASSERT_EQ(28, mesh->m_edges[47].first);
    ASSERT_EQ(20, mesh->m_edges[47].second);

    ASSERT_EQ(10, mesh->m_edges[48].first);
    ASSERT_EQ(27, mesh->m_edges[48].second);
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

    Polygons polygon(point, mesh->m_projection);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 2;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;

    MeshRefinement meshRefinement(mesh, polygon, meshRefinementParameters);
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

    const auto interpolator = std::make_shared<AveragingInterpolation>(*mesh,
                                                                       samples,
                                                                       AveragingInterpolation::Method::MinAbsValue,
                                                                       Mesh::Location::Faces,
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

    MeshRefinement meshRefinement(mesh, interpolator, meshRefinementParameters);
    meshRefinement.Compute();

    ASSERT_EQ(60, mesh->GetNumEdges());
    ASSERT_EQ(32, mesh->GetNumNodes());

    // sides of the refined part
    ASSERT_EQ(5, mesh->m_edges[5].first);
    ASSERT_EQ(16, mesh->m_edges[5].second);

    ASSERT_EQ(16, mesh->m_edges[40].first);
    ASSERT_EQ(9, mesh->m_edges[40].second);

    ASSERT_EQ(9, mesh->m_edges[9].first);
    ASSERT_EQ(19, mesh->m_edges[9].second);

    ASSERT_EQ(19, mesh->m_edges[43].first);
    ASSERT_EQ(13, mesh->m_edges[43].second);

    ASSERT_EQ(6, mesh->m_edges[16].first);
    ASSERT_EQ(22, mesh->m_edges[16].second);

    ASSERT_EQ(22, mesh->m_edges[46].first);
    ASSERT_EQ(5, mesh->m_edges[46].second);

    ASSERT_EQ(7, mesh->m_edges[17].first);
    ASSERT_EQ(23, mesh->m_edges[17].second);

    ASSERT_EQ(23, mesh->m_edges[47].first);
    ASSERT_EQ(6, mesh->m_edges[47].second);
}

TEST(MeshRefinement, Refine_SphericalMesh_ShouldRefine)
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
    Polygons polygon;
    MeshRefinement meshRefinement(mesh, polygon, meshRefinementParameters);
    meshRefinement.Compute();

    // Assert, we passed from 36 to 49 nodes
    ASSERT_EQ(121, mesh->GetNumNodes());

    constexpr double tolerance = 1e-6;
    ASSERT_NEAR(41.100000000000001, mesh->m_nodes[0].x, tolerance);
    ASSERT_NEAR(41.103300000000004, mesh->m_nodes[6].x, tolerance);
    ASSERT_NEAR(41.106600000000000, mesh->m_nodes[12].x, tolerance);
    ASSERT_NEAR(41.109900000000003, mesh->m_nodes[18].x, tolerance);
    ASSERT_NEAR(41.113199999999999, mesh->m_nodes[24].x, tolerance);

    ASSERT_NEAR(41.100000000000001, mesh->m_nodes[0].y, tolerance);
    ASSERT_NEAR(41.100000000000001, mesh->m_nodes[6].y, tolerance);
    ASSERT_NEAR(41.100000000000001, mesh->m_nodes[12].y, tolerance);
    ASSERT_NEAR(41.100000000000001, mesh->m_nodes[18].y, tolerance);
    ASSERT_NEAR(41.100000000000001, mesh->m_nodes[24].y, tolerance);
}

TEST(MeshRefinement, RefineCurvilinearGrid)
{
    auto mesh = MakeCurvilinearGridForTesting();

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 1;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;

    Polygons polygon;
    MeshRefinement meshRefinement(mesh, polygon, meshRefinementParameters);
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

    Polygons polygon(point, mesh->m_projection);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 2;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;

    MeshRefinement meshRefinement(mesh, polygon, meshRefinementParameters);

    // Execute
    meshRefinement.Compute();

    // The test compares only the curcumcentres of some of the triangles.
    // Since the mesh is dominated by quadrilaterals, these will be sought first in the mesh administration
    // and appear first in the list of elements.
    // So, need to offset by some number (this number is 1173) for the start of the triangles.
    const meshkernel::UInt triangleStart = 1173;
    meshkernel::UInt position = triangleStart;

    // Assert circumcenters are correctly computed
    constexpr double tolerance = 1e-6;

    // Compare the x-location of the circumcentre.
    ASSERT_NEAR(1673.0860169014584, mesh->m_facesCircumcenters[position].x, tolerance);
    ++position;
    ASSERT_NEAR(1660.6851354957175, mesh->m_facesCircumcenters[position].x, tolerance);
    ++position;
    ASSERT_NEAR(1660.5667704694627, mesh->m_facesCircumcenters[position].x, tolerance);
    ++position;
    ASSERT_NEAR(1672.0912775041329, mesh->m_facesCircumcenters[position].x, tolerance);
    ++position;
    ASSERT_NEAR(1659.9354211078053, mesh->m_facesCircumcenters[position].x, tolerance);
    ++position;
    ASSERT_NEAR(1659.8248648846848, mesh->m_facesCircumcenters[position].x, tolerance);
    ++position;
    ASSERT_NEAR(1671.1074693287451, mesh->m_facesCircumcenters[position].x, tolerance);
    ++position;
    ASSERT_NEAR(1659.1859707978906, mesh->m_facesCircumcenters[position].x, tolerance);
    ++position;
    ASSERT_NEAR(1659.0828479935451, mesh->m_facesCircumcenters[position].x, tolerance);
    ++position;
    ASSERT_NEAR(1670.1135380487042, mesh->m_facesCircumcenters[position].x, tolerance);
    ++position;
    ASSERT_NEAR(1658.4375379474418, mesh->m_facesCircumcenters[position].x, tolerance);

    // Reset to the start of the triangles for comparing the y-location of the circumcentre.
    position = triangleStart;

    ASSERT_NEAR(645.10565853980427, mesh->m_facesCircumcenters[position].y, tolerance);
    ++position;
    ASSERT_NEAR(646.20173898461292, mesh->m_facesCircumcenters[position].y, tolerance);
    ++position;
    ASSERT_NEAR(654.66978646149676, mesh->m_facesCircumcenters[position].y, tolerance);
    ++position;
    ASSERT_NEAR(662.96936808193914, mesh->m_facesCircumcenters[position].y, tolerance);
    ++position;
    ASSERT_NEAR(664.17198930960728, mesh->m_facesCircumcenters[position].y, tolerance);
    ++position;
    ASSERT_NEAR(672.49588595953537, mesh->m_facesCircumcenters[position].y, tolerance);
    ++position;
    ASSERT_NEAR(680.82566423059006, mesh->m_facesCircumcenters[position].y, tolerance);
    ++position;
    ASSERT_NEAR(682.12720924968505, mesh->m_facesCircumcenters[position].y, tolerance);
    ++position;
    ASSERT_NEAR(690.31918193748425, mesh->m_facesCircumcenters[position].y, tolerance);
    ++position;
    ASSERT_NEAR(698.66471917887850, mesh->m_facesCircumcenters[position].y, tolerance);
    ++position;
    ASSERT_NEAR(700.06356972686194, mesh->m_facesCircumcenters[position].y, tolerance);
}

TEST(MeshRefinement, BilinearInterpolationWithGriddedSamplesOnLandShouldNotRefine)
{
    // Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 10.0, Projection::cartesian);

    std::vector values{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    Point origin{-5.0, -5.0};
    const auto interpolator = std::make_shared<BilinearInterpolationOnGriddedSamples>(*mesh, 2, 2, origin, 10.0, values);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 1;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;
    meshRefinementParameters.min_edge_size = 5.0;
    meshRefinementParameters.account_for_samples_outside = 1;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.refinement_type = 1;

    MeshRefinement meshRefinement(mesh, interpolator, meshRefinementParameters, true);

    // Execute
    meshRefinement.Compute();

    // Assert: all bathy values are positive and we are in land, so nothing gets refined
    ASSERT_EQ(4, mesh->GetNumEdges());
}

TEST(MeshRefinement, BilinearInterpolationWithGriddedSamplesOnLandAndSeaShouldRefine)
{
    // Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 10.0, Projection::cartesian);

    std::vector values{-1.0, -2.0, 3.0, -4.0, -5.0, 6.0, 7.0, 8.0, 9.0};
    Point origin{-5.0, -5.0};
    const auto interpolator = std::make_shared<BilinearInterpolationOnGriddedSamples>(*mesh, 3, 3, origin, 10.0, values);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 1;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;
    meshRefinementParameters.min_edge_size = 5.0;
    meshRefinementParameters.account_for_samples_outside = 1;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.refinement_type = 1;

    MeshRefinement meshRefinement(mesh, interpolator, meshRefinementParameters, true);

    // Execute
    meshRefinement.Compute();

    // Assert: all depth values are positive and we are in land, so nothing gets refined
    ASSERT_EQ(12, mesh->GetNumEdges());
}

TEST(MeshRefinement, BilinearInterpolationWithAllGriddedSamplesOnSeaShouldRefine)
{
    // Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 10.0, Projection::cartesian);

    std::vector values{-1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0, -8.0, -9.0};
    Point origin{-5.0, -5.0};
    const auto interpolator = std::make_shared<BilinearInterpolationOnGriddedSamples>(*mesh, 2, 2, origin, 10.0, values);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 1;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;
    meshRefinementParameters.min_edge_size = 5.0;
    meshRefinementParameters.account_for_samples_outside = 1;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.refinement_type = 1;

    MeshRefinement meshRefinement(mesh, interpolator, meshRefinementParameters, true);

    // Execute
    meshRefinement.Compute();

    // Assert: nothing gets refined
    ASSERT_EQ(4, mesh->GetNumEdges());
}
