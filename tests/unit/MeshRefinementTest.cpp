#include <MeshKernel/MeshRefinement.hpp>
#include <MeshKernel/Mesh.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/SampleRefineParametersNative.hpp>
#include <MeshKernel/InterpolationParametersNative.hpp>
#include <TestUtils/MakeMeshes.hpp>
#include <TestUtils/SampleFileReader.hpp>
#include <gtest/gtest.h>

TEST(MeshRefinement, FourByFourWithFourSamples)
{
    auto mesh = MakeRectangularMeshForTesting(5, 5, 10.0, meshkernel::Projections::cartesian);

    //sample points
    std::vector<meshkernel::Sample> samples{
        {14.7153645, 14.5698833, 1.0},
        {24.7033062, 14.4729137, 1.0},
        {15.5396099, 24.2669525, 1.0},
        {23.8305721, 23.9275551, 1.0}};

    const auto averaging = std::make_shared<meshkernel::AveragingInterpolation>(mesh,
                                                                                samples,
                                                                                meshkernel::AveragingInterpolation::Method::MinAbsValue,
                                                                                meshkernel::InterpolationLocation::Faces,
                                                                                1.0,
                                                                                false,
                                                                                false);

    meshkernel::MeshRefinement meshRefinement(mesh, averaging);
    meshkernel::Polygons polygon;
    meshkernelapi::SampleRefineParametersNative sampleRefineParametersNative;
    sampleRefineParametersNative.MaximumTimeStepInCourantGrid = 0.32;
    sampleRefineParametersNative.MinimumCellSize = 1.0;
    sampleRefineParametersNative.AccountForSamplesOutside = false;
    sampleRefineParametersNative.ConnectHangingNodes = 1;

    meshkernelapi::InterpolationParametersNative interpolationParametersNative;
    interpolationParametersNative.MaxNumberOfRefinementIterations = 1;
    interpolationParametersNative.RefineIntersected = false;

    meshRefinement.Refine(polygon, sampleRefineParametersNative, interpolationParametersNative);

    // 3 Validation edges connecting hanging nodes

    //bottom side
    ASSERT_EQ(5, mesh->m_edges[73].first);
    ASSERT_EQ(25, mesh->m_edges[73].second);

    ASSERT_EQ(10, mesh->m_edges[72].first);
    ASSERT_EQ(25, mesh->m_edges[72].second);

    ASSERT_EQ(10, mesh->m_edges[77].first);
    ASSERT_EQ(28, mesh->m_edges[77].second);

    ASSERT_EQ(15, mesh->m_edges[76].first);
    ASSERT_EQ(28, mesh->m_edges[76].second);

    //right side
    ASSERT_EQ(21, mesh->m_edges[81].first);
    ASSERT_EQ(35, mesh->m_edges[81].second);

    ASSERT_EQ(22, mesh->m_edges[80].first);
    ASSERT_EQ(35, mesh->m_edges[80].second);

    ASSERT_EQ(22, mesh->m_edges[83].first);
    ASSERT_EQ(36, mesh->m_edges[83].second);

    ASSERT_EQ(23, mesh->m_edges[82].first);
    ASSERT_EQ(36, mesh->m_edges[82].second);

    //upper side
    ASSERT_EQ(19, mesh->m_edges[79].first);
    ASSERT_EQ(30, mesh->m_edges[79].second);

    ASSERT_EQ(14, mesh->m_edges[78].first);
    ASSERT_EQ(30, mesh->m_edges[78].second);

    ASSERT_EQ(14, mesh->m_edges[75].first);
    ASSERT_EQ(27, mesh->m_edges[75].second);

    ASSERT_EQ(9, mesh->m_edges[74].first);
    ASSERT_EQ(27, mesh->m_edges[74].second);

    //left side
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

TEST(MeshRefinement, FourByFourWithFourSamplesEdgeSizeTwo)
{
    auto mesh = MakeRectangularMeshForTesting(4, 4, 10.0, meshkernel::Projections::cartesian);

    //sample points
    std::vector<meshkernel::Sample> samples{
        {14.7153645, 14.5698833, 1.0},
        {24.7033062, 14.4729137, 1.0},
        {15.5396099, 24.2669525, 1.0},
        {23.8305721, 23.9275551, 1.0}};

    const auto averaging = std::make_shared<meshkernel::AveragingInterpolation>(mesh,
                                                                                samples,
                                                                                meshkernel::AveragingInterpolation::Method::MinAbsValue,
                                                                                meshkernel::InterpolationLocation::Faces,
                                                                                1.0,
                                                                                false,
                                                                                false);

    meshkernel::MeshRefinement meshRefinement(mesh, averaging);
    meshkernel::Polygons polygon;
    meshkernelapi::SampleRefineParametersNative sampleRefineParametersNative;
    sampleRefineParametersNative.MaximumTimeStepInCourantGrid = 0.64;
    sampleRefineParametersNative.MinimumCellSize = 2.0;
    sampleRefineParametersNative.AccountForSamplesOutside = false;
    sampleRefineParametersNative.ConnectHangingNodes = 1;

    meshkernelapi::InterpolationParametersNative interpolationParametersNative;
    interpolationParametersNative.MaxNumberOfRefinementIterations = 4;
    interpolationParametersNative.RefineIntersected = false;

    meshRefinement.Refine(polygon, sampleRefineParametersNative, interpolationParametersNative);

    //Assert number of edges and nodes
    ASSERT_EQ(131, mesh->GetNumEdges());
    ASSERT_EQ(62, mesh->GetNumNodes());

    //Assert edges
    ASSERT_EQ(0, mesh->m_edges[0].first);
    ASSERT_EQ(4, mesh->m_edges[0].second);

    ASSERT_EQ(1, mesh->m_edges[1].first);
    ASSERT_EQ(32, mesh->m_edges[1].second);

    ASSERT_EQ(2, mesh->m_edges[2].first);
    ASSERT_EQ(33, mesh->m_edges[2].second);

    ASSERT_EQ(3, mesh->m_edges[3].first);
    ASSERT_EQ(7, mesh->m_edges[3].second);

    ASSERT_EQ(4, mesh->m_edges[4].first);
    ASSERT_EQ(34, mesh->m_edges[4].second);

    ASSERT_EQ(5, mesh->m_edges[5].first);
    ASSERT_EQ(35, mesh->m_edges[5].second);

    ASSERT_EQ(6, mesh->m_edges[6].first);
    ASSERT_EQ(17, mesh->m_edges[6].second);

    ASSERT_EQ(7, mesh->m_edges[7].first);
    ASSERT_EQ(18, mesh->m_edges[7].second);

    ASSERT_EQ(8, mesh->m_edges[8].first);
    ASSERT_EQ(36, mesh->m_edges[8].second);

    ASSERT_EQ(9, mesh->m_edges[9].first);
    ASSERT_EQ(37, mesh->m_edges[9].second);

    ASSERT_EQ(10, mesh->m_edges[10].first);
    ASSERT_EQ(38, mesh->m_edges[10].second);

    ASSERT_EQ(11, mesh->m_edges[11].first);
    ASSERT_EQ(21, mesh->m_edges[11].second);

    ASSERT_EQ(1, mesh->m_edges[12].first);
    ASSERT_EQ(0, mesh->m_edges[12].second);

    ASSERT_EQ(2, mesh->m_edges[13].first);
    ASSERT_EQ(39, mesh->m_edges[13].second);

    ASSERT_EQ(3, mesh->m_edges[14].first);
    ASSERT_EQ(2, mesh->m_edges[14].second);

    ASSERT_EQ(5, mesh->m_edges[15].first);
    ASSERT_EQ(40, mesh->m_edges[15].second);

    ASSERT_EQ(6, mesh->m_edges[16].first);
    ASSERT_EQ(22, mesh->m_edges[16].second);

    ASSERT_EQ(7, mesh->m_edges[17].first);
    ASSERT_EQ(23, mesh->m_edges[17].second);

    ASSERT_EQ(9, mesh->m_edges[18].first);
    ASSERT_EQ(41, mesh->m_edges[18].second);

    ASSERT_EQ(10, mesh->m_edges[19].first);
    ASSERT_EQ(24, mesh->m_edges[19].second);

    ASSERT_EQ(11, mesh->m_edges[20].first);
    ASSERT_EQ(25, mesh->m_edges[20].second);
}

TEST(MeshRefinement, SmallTriangualMeshTwoSamples)
{
    // Prepare
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/SmallTriangularGrid_net.nc");

    //sample points
    std::vector<meshkernel::Sample> samples{
        {359.8657532, 350.3144836, 1.0},
        {387.5152588, 299.2614746, 1.0}};

    const auto averaging = std::make_shared<meshkernel::AveragingInterpolation>(mesh,
                                                                                samples,
                                                                                meshkernel::AveragingInterpolation::Method::MinAbsValue,
                                                                                meshkernel::InterpolationLocation::Faces,
                                                                                1.0,
                                                                                false,
                                                                                false);

    meshkernel::MeshRefinement meshRefinement(mesh, averaging);
    meshkernel::Polygons polygon;
    meshkernelapi::SampleRefineParametersNative sampleRefineParametersNative;
    sampleRefineParametersNative.MaximumTimeStepInCourantGrid = 15.97;
    sampleRefineParametersNative.MinimumCellSize = 50.0;
    sampleRefineParametersNative.AccountForSamplesOutside = false;
    sampleRefineParametersNative.ConnectHangingNodes = 1;

    meshkernelapi::InterpolationParametersNative interpolationParametersNative;
    interpolationParametersNative.MaxNumberOfRefinementIterations = 1;
    interpolationParametersNative.RefineIntersected = false;

    meshRefinement.Refine(polygon, sampleRefineParametersNative, interpolationParametersNative);

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
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/SmallTriangularGrid_net.nc");

    // Polygon sample
    std::vector<meshkernel::Point> point{
        {399.638169557229, 504.294564030922},
        {361.827403800769, 129.967983041964},
        {651.709941266965, 113.583317880831},
        {666.834247569549, 411.028008498319},
        {410.981399284167, 505.55492288947},
        {399.638169557229, 504.294564030922}};

    meshkernel::Polygons polygon(point, mesh->m_projection);

    meshkernel::MeshRefinement meshRefinement(mesh);
    meshkernelapi::SampleRefineParametersNative sampleRefineParametersNative;

    meshkernelapi::InterpolationParametersNative interpolationParametersNative;
    interpolationParametersNative.MaxNumberOfRefinementIterations = 1;
    interpolationParametersNative.RefineIntersected = false;

    std::vector<meshkernel::Sample> samples;
    meshRefinement.Refine(polygon, sampleRefineParametersNative, interpolationParametersNative);

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

TEST(MeshRefinement, ThreeBythreeWithThreeSamplesPerface)
{
    // Prepare

    auto mesh = MakeRectangularMeshForTesting(4, 4, 10.0, meshkernel::Projections::cartesian);

    //sample points
    std::vector<meshkernel::Sample> samples{
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

    const auto averaging = std::make_shared<meshkernel::AveragingInterpolation>(mesh,
                                                                                samples,
                                                                                meshkernel::AveragingInterpolation::Method::MinAbsValue,
                                                                                meshkernel::InterpolationLocation::Faces,
                                                                                1.0,
                                                                                false,
                                                                                false);

    meshkernel::MeshRefinement meshRefinement(mesh, averaging);
    meshkernel::Polygons polygon;
    meshkernelapi::SampleRefineParametersNative sampleRefineParametersNative;
    sampleRefineParametersNative.MaximumTimeStepInCourantGrid = 0.96;
    sampleRefineParametersNative.MinimumCellSize = 3.0;
    sampleRefineParametersNative.AccountForSamplesOutside = false;
    sampleRefineParametersNative.ConnectHangingNodes = 1;

    meshkernelapi::InterpolationParametersNative interpolationParametersNative;
    interpolationParametersNative.MaxNumberOfRefinementIterations = 2;
    interpolationParametersNative.RefineIntersected = false;

    meshRefinement.Refine(polygon, sampleRefineParametersNative, interpolationParametersNative);

    // assert on number of nodes and edges
    ASSERT_EQ(150, mesh->GetNumNodes());
    ASSERT_EQ(293, mesh->GetNumEdges());

    // assert on newly generated edges
    ASSERT_EQ(109, mesh->m_edges[282].first);
    ASSERT_EQ(44, mesh->m_edges[282].second);

    ASSERT_EQ(44, mesh->m_edges[283].first);
    ASSERT_EQ(67, mesh->m_edges[283].second);

    ASSERT_EQ(92, mesh->m_edges[284].first);
    ASSERT_EQ(91, mesh->m_edges[284].second);

    ASSERT_EQ(91, mesh->m_edges[285].first);
    ASSERT_EQ(12, mesh->m_edges[285].second);

    ASSERT_EQ(12, mesh->m_edges[286].first);
    ASSERT_EQ(92, mesh->m_edges[286].second);

    ASSERT_EQ(100, mesh->m_edges[287].first);
    ASSERT_EQ(99, mesh->m_edges[287].second);

    ASSERT_EQ(99, mesh->m_edges[288].first);
    ASSERT_EQ(14, mesh->m_edges[288].second);

    ASSERT_EQ(14, mesh->m_edges[289].first);
    ASSERT_EQ(100, mesh->m_edges[289].second);

    ASSERT_EQ(97, mesh->m_edges[290].first);
    ASSERT_EQ(96, mesh->m_edges[290].second);

    ASSERT_EQ(96, mesh->m_edges[291].first);
    ASSERT_EQ(14, mesh->m_edges[291].second);

    ASSERT_EQ(14, mesh->m_edges[292].first);
    ASSERT_EQ(97, mesh->m_edges[292].second);
}

TEST(MeshRefinement, WindowOfRefinementFile)
{
    // Prepare
    auto mesh = MakeRectangularMeshForTesting(4, 4, 40.0, meshkernel::Projections::cartesian, {197253.0, 442281.0});

    // Sample points
    std::vector<meshkernel::Sample> samples = ReadSampleFile("../../../../tests/data/MeshRefinementTests/WindowOfRefinementFile.xyz");

    const auto averaging = std::make_shared<meshkernel::AveragingInterpolation>(mesh,
                                                                                samples,
                                                                                meshkernel::AveragingInterpolation::Method::MinAbsValue,
                                                                                meshkernel::InterpolationLocation::Faces,
                                                                                1.0,
                                                                                false,
                                                                                false);

    meshkernel::MeshRefinement meshRefinement(mesh, averaging);
    meshkernel::Polygons polygon;
    meshkernelapi::SampleRefineParametersNative sampleRefineParametersNative;
    sampleRefineParametersNative.MaximumTimeStepInCourantGrid = 0.96;
    sampleRefineParametersNative.MinimumCellSize = 3.0;
    sampleRefineParametersNative.AccountForSamplesOutside = false;
    sampleRefineParametersNative.ConnectHangingNodes = 1;

    meshkernelapi::InterpolationParametersNative interpolationParametersNative;
    interpolationParametersNative.MaxNumberOfRefinementIterations = 4;
    interpolationParametersNative.RefineIntersected = false;

    meshRefinement.Refine(polygon, sampleRefineParametersNative, interpolationParametersNative);

    // total number of edges
    ASSERT_EQ(1614, mesh->GetNumNodes());
    ASSERT_EQ(3216, mesh->GetNumEdges());

    ASSERT_EQ(170, mesh->m_edges[3113].first);
    ASSERT_EQ(804, mesh->m_edges[3113].second);

    ASSERT_EQ(1, mesh->m_edges[3114].first);
    ASSERT_EQ(804, mesh->m_edges[3114].second);

    ASSERT_EQ(462, mesh->m_edges[3115].first);
    ASSERT_EQ(856, mesh->m_edges[3115].second);

    ASSERT_EQ(9, mesh->m_edges[3116].first);
    ASSERT_EQ(856, mesh->m_edges[3116].second);

    ASSERT_EQ(211, mesh->m_edges[3117].first);
    ASSERT_EQ(1256, mesh->m_edges[3117].second);

    ASSERT_EQ(538, mesh->m_edges[3118].first);
    ASSERT_EQ(1256, mesh->m_edges[3118].second);

    ASSERT_EQ(77, mesh->m_edges[3119].first);
    ASSERT_EQ(1052, mesh->m_edges[3119].second);

    ASSERT_EQ(258, mesh->m_edges[3120].first);
    ASSERT_EQ(1052, mesh->m_edges[3120].second);

    ASSERT_EQ(290, mesh->m_edges[3121].first);
    ASSERT_EQ(769, mesh->m_edges[3121].second);

    ASSERT_EQ(593, mesh->m_edges[3122].first);
    ASSERT_EQ(769, mesh->m_edges[3122].second);
}

TEST(MeshRefinement, WindowOfRefinementFileBasedOnLevels)
{
    // Prepare
    auto mesh = MakeRectangularMeshForTesting(4, 4, 40.0, meshkernel::Projections::cartesian, {197253.0, 442281.0});

    // Sample points
    std::vector<meshkernel::Sample> samples = ReadSampleFile("../../../../tests/data/MeshRefinementTests/WindowOfRefinementFile.xyz");

    const auto averaging = std::make_shared<meshkernel::AveragingInterpolation>(mesh,
                                                                                samples,
                                                                                meshkernel::AveragingInterpolation::Method::Max,
                                                                                meshkernel::InterpolationLocation::Faces,
                                                                                1.01,
                                                                                false,
                                                                                true);

    meshkernel::MeshRefinement meshRefinement(mesh, averaging);
    meshkernel::Polygons polygon;
    meshkernelapi::SampleRefineParametersNative sampleRefineParametersNative;
    sampleRefineParametersNative.MaximumTimeStepInCourantGrid = 0.96;
    sampleRefineParametersNative.MinimumCellSize = 0.5;
    sampleRefineParametersNative.AccountForSamplesOutside = false;
    sampleRefineParametersNative.ConnectHangingNodes = 1;
    sampleRefineParametersNative.RefinementType = 3;

    meshkernelapi::InterpolationParametersNative interpolationParametersNative;
    interpolationParametersNative.MaxNumberOfRefinementIterations = 10;
    interpolationParametersNative.RefineIntersected = false;

    meshRefinement.Refine(polygon, sampleRefineParametersNative, interpolationParametersNative);

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
    auto mesh = MakeRectangularMeshForTesting(5, 5, 10.0, meshkernel::Projections::cartesian);

    meshkernel::MeshRefinement meshRefinement(mesh);

    std::vector<meshkernel::Point> point{
        {25.0, -10.0},
        {25.0, 15.0},
        {45.0, 15.0},
        {45.0, -10.0},
        {25.0, -10.0}};

    meshkernel::Polygons polygon(point, mesh->m_projection);

    meshkernelapi::SampleRefineParametersNative sampleRefineParametersNative;
    sampleRefineParametersNative.MaximumTimeStepInCourantGrid = 0.96;
    sampleRefineParametersNative.MinimumCellSize = 3.0;
    sampleRefineParametersNative.AccountForSamplesOutside = false;
    sampleRefineParametersNative.ConnectHangingNodes = 1;

    meshkernelapi::InterpolationParametersNative interpolationParametersNative;
    interpolationParametersNative.MaxNumberOfRefinementIterations = 1;
    interpolationParametersNative.RefineIntersected = false;

    meshRefinement.Refine(polygon, sampleRefineParametersNative, interpolationParametersNative);

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
    auto mesh = MakeRectangularMeshForTesting(4, 4, 10.0, meshkernel::Projections::cartesian);

    meshkernel::MeshRefinement meshRefinement(mesh);

    std::vector<meshkernel::Point> point{
        {9.09836065573771, 34.016393442623},
        {7.18032786885247, 7.75409836065574},
        {34.6229508196721, 6.5},
        {34.4194409808304, 26.6983515050386},
        {34.327868852459, 35.7868852459016},
        {29.0521194370216, 35.4840476661635},
        {9.90983606557378, 34.3852459016394},
        {9.09836065573771, 34.016393442623}};

    meshkernel::Polygons polygon(point, mesh->m_projection);

    meshkernelapi::SampleRefineParametersNative sampleRefineParametersNative;
    sampleRefineParametersNative.MaximumTimeStepInCourantGrid = 0.32;
    sampleRefineParametersNative.MinimumCellSize = 1.0;
    sampleRefineParametersNative.AccountForSamplesOutside = false;
    sampleRefineParametersNative.ConnectHangingNodes = 1;
    sampleRefineParametersNative.RefinementType = 3;

    meshkernelapi::InterpolationParametersNative interpolationParametersNative;
    interpolationParametersNative.MaxNumberOfRefinementIterations = 2;
    interpolationParametersNative.RefineIntersected = false;

    meshRefinement.Refine(polygon, sampleRefineParametersNative, interpolationParametersNative);

    // assert on number of nodes and edges
    ASSERT_EQ(48, mesh->GetNumNodes());
    ASSERT_EQ(96, mesh->GetNumEdges());
    ASSERT_EQ(49, mesh->GetNumFaces());
}

TEST(MeshRefinement, FourByFourWithFourSamplesSpherical)
{

    auto mesh = MakeRectangularMeshForTesting(4, 4, 0.0033, meshkernel::Projections::spherical, {41.1, 41.1});

    //sample points
    std::vector<meshkernel::Sample> samples{
        {41.1050110, 41.1049728, 1.0},
        {41.1084785, 41.1048775, 1.0},
        {41.1085625, 41.1083946, 1.0},
        {41.1052971, 41.1083336, 1.0}};

    const auto averaging = std::make_shared<meshkernel::AveragingInterpolation>(mesh,
                                                                                samples,
                                                                                meshkernel::AveragingInterpolation::Method::MinAbsValue,
                                                                                meshkernel::InterpolationLocation::Faces,
                                                                                1.0,
                                                                                false,
                                                                                false);

    meshkernel::MeshRefinement meshRefinement(mesh, averaging);
    meshkernel::Polygons polygon;
    meshkernelapi::SampleRefineParametersNative sampleRefineParametersNative;
    sampleRefineParametersNative.MaximumTimeStepInCourantGrid = 0.000527;
    sampleRefineParametersNative.MinimumCellSize = 0.00165;
    sampleRefineParametersNative.AccountForSamplesOutside = false;
    sampleRefineParametersNative.ConnectHangingNodes = 1;

    meshkernelapi::InterpolationParametersNative interpolationParametersNative;
    interpolationParametersNative.MaxNumberOfRefinementIterations = 1;
    interpolationParametersNative.RefineIntersected = false;
    sampleRefineParametersNative.RefinementType = 2;
    meshRefinement.Refine(polygon, sampleRefineParametersNative, interpolationParametersNative);

    ASSERT_EQ(60, mesh->GetNumEdges());
    ASSERT_EQ(32, mesh->GetNumNodes());

    //sides of the refined part
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

TEST(MeshRefinement, RefineCurvilinearGrid)
{
    auto mesh = MakeCurvilinearGridForTesting();

    meshkernel::MeshRefinement meshRefinement(mesh);
    meshkernel::Polygons polygon;
    meshkernelapi::SampleRefineParametersNative sampleRefineParametersNative;
    sampleRefineParametersNative.MaximumTimeStepInCourantGrid = 0.000527;
    sampleRefineParametersNative.MinimumCellSize = 0.00165;
    sampleRefineParametersNative.AccountForSamplesOutside = false;
    sampleRefineParametersNative.ConnectHangingNodes = 1;

    meshkernelapi::InterpolationParametersNative interpolationParametersNative;
    interpolationParametersNative.MaxNumberOfRefinementIterations = 1;
    interpolationParametersNative.RefineIntersected = false;
    sampleRefineParametersNative.RefinementType = 2;
    meshRefinement.Refine(polygon, sampleRefineParametersNative, interpolationParametersNative);

    mesh->ComputeEdgeLengths();

    // if the circumcenters are wrongly computed, some edges will be smaller than half cell size
    for (int i = 0; i < mesh->GetNumEdges(); ++i)
    {
        ASSERT_GT(mesh->m_edgeLengths[i], 0.4);
    }
}
