
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <MeshKernel/AveragingInterpolation.hpp>
#include <TestUtils/Definitions.hpp>
#include <TestUtils/MakeMeshes.hpp>
#include <TestUtils/SampleFileReader.hpp>

// Simple averaging
TEST(Averaging, AveragingInterpolation_OnNodesWithSphericalCoordinates_Shouldinterpolate)
{
    // Prepare
    auto mesh = MakeRectangularMeshForTesting(5, 5, 1.0, meshkernel::Projection::spherical);

    std::vector<meshkernel::Sample> samples{
        {1.5, 1.5, 2.0},
        {2.5, 1.5, 2.0},
        {3.5, 1.5, 2.0},
        {1.5, 2.5, 2.0},
        {2.5, 2.5, 1.0},
        {3.5, 2.5, 2.0},
        {1.5, 3.5, 2.0},
        {2.5, 3.5, 2.0},
        {3.5, 3.5, 2.0}};

    // Execute
    meshkernel::AveragingInterpolation averaging(*mesh,
                                                 samples,
                                                 meshkernel::AveragingInterpolation::Method::SimpleAveraging,
                                                 meshkernel::Mesh::Location::Nodes,
                                                 1.01,
                                                 false,
                                                 false,
                                                 1);
    averaging.Compute();

    // Asser

    constexpr double tolerance = 1e-6;
    ASSERT_NEAR(-999.0, averaging.GetNodeResult(0), tolerance);
    ASSERT_NEAR(-999.0, averaging.GetNodeResult(1), tolerance);
    ASSERT_NEAR(-999.0, averaging.GetNodeResult(2), tolerance);
    ASSERT_NEAR(-999.0, averaging.GetNodeResult(3), tolerance);
    ASSERT_NEAR(-999.0, averaging.GetNodeResult(4), tolerance);
    ASSERT_NEAR(-999.0, averaging.GetNodeResult(5), tolerance);
    ASSERT_NEAR(2.0000, averaging.GetNodeResult(6), tolerance);
    ASSERT_NEAR(2.0000, averaging.GetNodeResult(7), tolerance);
    ASSERT_NEAR(2.0000, averaging.GetNodeResult(8), tolerance);
    ASSERT_NEAR(2.0000, averaging.GetNodeResult(9), tolerance);
    ASSERT_NEAR(-999.0, averaging.GetNodeResult(10), tolerance);
    ASSERT_NEAR(2.0000, averaging.GetNodeResult(11), tolerance);
    ASSERT_NEAR(1.7500, averaging.GetNodeResult(12), tolerance);
    ASSERT_NEAR(1.7500, averaging.GetNodeResult(13), tolerance);
    ASSERT_NEAR(2.0000, averaging.GetNodeResult(14), tolerance);
    ASSERT_NEAR(-999.0, averaging.GetNodeResult(15), tolerance);
    ASSERT_NEAR(2.0000, averaging.GetNodeResult(16), tolerance);
    ASSERT_NEAR(1.7500, averaging.GetNodeResult(17), tolerance);
    ASSERT_NEAR(1.7500, averaging.GetNodeResult(18), tolerance);
    ASSERT_NEAR(2.0000, averaging.GetNodeResult(19), tolerance);
    ASSERT_NEAR(-999.0, averaging.GetNodeResult(20), tolerance);
    ASSERT_NEAR(2.0000, averaging.GetNodeResult(21), tolerance);
    ASSERT_NEAR(2.0000, averaging.GetNodeResult(22), tolerance);
    ASSERT_NEAR(2.0000, averaging.GetNodeResult(23), tolerance);
    ASSERT_NEAR(2.0000, averaging.GetNodeResult(24), tolerance);
}
TEST(Averaging, InterpolateOnEdgesSimpleAveraging)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/AveragingInterpolationTests/simple_grid_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), static_cast<meshkernel::Index>(0));

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(*mesh,
                                                 samples,
                                                 meshkernel::AveragingInterpolation::Method::SimpleAveraging,
                                                 meshkernel::Mesh::Location::Edges,
                                                 1.01,
                                                 false,
                                                 false,
                                                 0);
    averaging.Compute();

    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(-2.28291315000000, averaging.GetEdgeResult(0), tolerance);
    ASSERT_NEAR(-2.75910365000000, averaging.GetEdgeResult(1), tolerance);
    ASSERT_NEAR(-3.23529410000000, averaging.GetEdgeResult(2), tolerance);
    ASSERT_NEAR(-3.71148460000000, averaging.GetEdgeResult(3), tolerance);
    ASSERT_NEAR(-4.18767505000000, averaging.GetEdgeResult(4), tolerance);
    ASSERT_NEAR(-2.67507005000000, averaging.GetEdgeResult(5), tolerance);
    ASSERT_NEAR(-3.15126050000000, averaging.GetEdgeResult(6), tolerance);
    ASSERT_NEAR(-3.62745095000000, averaging.GetEdgeResult(7), tolerance);
    ASSERT_NEAR(-4.10364145000000, averaging.GetEdgeResult(8), tolerance);
    ASSERT_NEAR(-4.57983195000000, averaging.GetEdgeResult(9), tolerance);
    ASSERT_NEAR(-3.06722690000000, averaging.GetEdgeResult(10), tolerance);
    ASSERT_NEAR(-3.54341735000000, averaging.GetEdgeResult(11), tolerance);
    ASSERT_NEAR(-4.01960785000000, averaging.GetEdgeResult(12), tolerance);
    ASSERT_NEAR(-4.49579835000000, averaging.GetEdgeResult(13), tolerance);
    ASSERT_NEAR(-4.97198880000000, averaging.GetEdgeResult(14), tolerance);
    ASSERT_NEAR(-3.45938375000000, averaging.GetEdgeResult(15), tolerance);
    ASSERT_NEAR(-3.93557425000000, averaging.GetEdgeResult(16), tolerance);
    ASSERT_NEAR(-4.41176470000000, averaging.GetEdgeResult(17), tolerance);
    ASSERT_NEAR(-4.88795520000000, averaging.GetEdgeResult(18), tolerance);
    ASSERT_NEAR(-5.36414565000000, averaging.GetEdgeResult(19), tolerance);
}

TEST(Averaging, InterpolateOnNodesSimpleAveraging)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/AveragingInterpolationTests/simple_grid_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), static_cast<meshkernel::Index>(0));

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(*mesh,
                                                 samples,
                                                 meshkernel::AveragingInterpolation::Method::SimpleAveraging,
                                                 meshkernel::Mesh::Location::Nodes,
                                                 1.01,
                                                 false,
                                                 false,
                                                 0);
    averaging.Compute();

    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(-2.08683470000000, averaging.GetNodeResult(0), tolerance);
    ASSERT_NEAR(-2.47899160000000, averaging.GetNodeResult(1), tolerance);
    ASSERT_NEAR(-2.56302520000000, averaging.GetNodeResult(2), tolerance);
    ASSERT_NEAR(-2.95518210000000, averaging.GetNodeResult(3), tolerance);
    ASSERT_NEAR(-3.03921570000000, averaging.GetNodeResult(4), tolerance);
    ASSERT_NEAR(-3.43137250000000, averaging.GetNodeResult(5), tolerance);
    ASSERT_NEAR(-3.51540620000000, averaging.GetNodeResult(6), tolerance);
    ASSERT_NEAR(-3.90756300000000, averaging.GetNodeResult(7), tolerance);
    ASSERT_NEAR(-3.99159660000000, averaging.GetNodeResult(8), tolerance);
    ASSERT_NEAR(-4.38375350000000, averaging.GetNodeResult(9), tolerance);
    ASSERT_NEAR(-2.87114850000000, averaging.GetNodeResult(10), tolerance);
    ASSERT_NEAR(-3.34733890000000, averaging.GetNodeResult(11), tolerance);
    ASSERT_NEAR(-3.82352940000000, averaging.GetNodeResult(12), tolerance);
    ASSERT_NEAR(-4.29971990000000, averaging.GetNodeResult(13), tolerance);
    ASSERT_NEAR(-4.77591040000000, averaging.GetNodeResult(14), tolerance);
    ASSERT_NEAR(-3.26330530000000, averaging.GetNodeResult(15), tolerance);
    ASSERT_NEAR(-3.73949580000000, averaging.GetNodeResult(16), tolerance);
    ASSERT_NEAR(-4.21568630000000, averaging.GetNodeResult(17), tolerance);
    ASSERT_NEAR(-4.69187680000000, averaging.GetNodeResult(18), tolerance);
    ASSERT_NEAR(-5.16806720000000, averaging.GetNodeResult(19), tolerance);
}

TEST(Averaging, InterpolateOnFacesSimpleAveraging)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/AveragingInterpolationTests/simple_grid_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), static_cast<meshkernel::Index>(0));

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(*mesh,
                                                 samples,
                                                 meshkernel::AveragingInterpolation::Method::SimpleAveraging,
                                                 meshkernel::Mesh::Location::Faces,
                                                 1.01,
                                                 false,
                                                 false,
                                                 1);
    averaging.Compute();

    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(-2.52100840000000, averaging.GetFaceResult(0), tolerance);
    ASSERT_NEAR(-2.91316527500000, averaging.GetFaceResult(1), tolerance);
    ASSERT_NEAR(-2.99719887500000, averaging.GetFaceResult(2), tolerance);
    ASSERT_NEAR(-3.38935572500000, averaging.GetFaceResult(3), tolerance);
    ASSERT_NEAR(-3.47338935000000, averaging.GetFaceResult(4), tolerance);
    ASSERT_NEAR(-3.86554620000000, averaging.GetFaceResult(5), tolerance);
    ASSERT_NEAR(-3.94957982500000, averaging.GetFaceResult(6), tolerance);
    ASSERT_NEAR(-4.34173670000000, averaging.GetFaceResult(7), tolerance);
    ASSERT_NEAR(-3.30532212500000, averaging.GetFaceResult(8), tolerance);
    ASSERT_NEAR(-3.78151260000000, averaging.GetFaceResult(9), tolerance);
    ASSERT_NEAR(-4.25770310000000, averaging.GetFaceResult(10), tolerance);
    ASSERT_NEAR(-4.73389357500000, averaging.GetFaceResult(11), tolerance);
    ASSERT_NEAR(-3.69747900000000, averaging.GetFaceResult(12), tolerance);
    ASSERT_NEAR(-4.17366947500000, averaging.GetFaceResult(13), tolerance);
    ASSERT_NEAR(-4.64985995000000, averaging.GetFaceResult(14), tolerance);
    ASSERT_NEAR(-5.12605042500000, averaging.GetFaceResult(15), tolerance);
    ASSERT_NEAR(-4.08963585000000, averaging.GetFaceResult(16), tolerance);
    ASSERT_NEAR(-4.56582632500000, averaging.GetFaceResult(17), tolerance);
    ASSERT_NEAR(-5.04201680000000, averaging.GetFaceResult(18), tolerance);
    ASSERT_NEAR(-5.51820730000000, averaging.GetFaceResult(19), tolerance);
}

// Closest Point
TEST(Averaging, InterpolateOnEdgesClosestPoint)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/AveragingInterpolationTests/sample_grid_coarse_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), static_cast<meshkernel::Index>(0));

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(*mesh,
                                                 samples,
                                                 meshkernel::AveragingInterpolation::Method::Closest,
                                                 meshkernel::Mesh::Location::Edges,
                                                 1.01,
                                                 false,
                                                 false,
                                                 1);
    averaging.Compute();

    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(-2.47899160000000, averaging.GetEdgeResult(0), tolerance);
    ASSERT_NEAR(-3.43137255000000, averaging.GetEdgeResult(1), tolerance);
    ASSERT_NEAR(-4.38375350000000, averaging.GetEdgeResult(2), tolerance);
    ASSERT_NEAR(-3.26330535000000, averaging.GetEdgeResult(3), tolerance);
    ASSERT_NEAR(-4.21568625000000, averaging.GetEdgeResult(4), tolerance);
    ASSERT_NEAR(-5.16806725000000, averaging.GetEdgeResult(5), tolerance);
    ASSERT_NEAR(-4.04761905000000, averaging.GetEdgeResult(6), tolerance);
    ASSERT_NEAR(-5.00000000000000, averaging.GetEdgeResult(7), tolerance);
    ASSERT_NEAR(-5.95238095000000, averaging.GetEdgeResult(8), tolerance);
    ASSERT_NEAR(-4.83193275000000, averaging.GetEdgeResult(9), tolerance);
    ASSERT_NEAR(-5.78431375000000, averaging.GetEdgeResult(10), tolerance);
    ASSERT_NEAR(-6.73669465000000, averaging.GetEdgeResult(11), tolerance);
    ASSERT_NEAR(-5.61624650000000, averaging.GetEdgeResult(12), tolerance);
    ASSERT_NEAR(-6.56862745000000, averaging.GetEdgeResult(13), tolerance);
    ASSERT_NEAR(-7.52100840000000, averaging.GetEdgeResult(14), tolerance);
    ASSERT_NEAR(-2.56302520000000, averaging.GetEdgeResult(15), tolerance);
    ASSERT_NEAR(-3.51540615000000, averaging.GetEdgeResult(16), tolerance);
    ASSERT_NEAR(-3.34733895000000, averaging.GetEdgeResult(17), tolerance);
    ASSERT_NEAR(-4.29971990000000, averaging.GetEdgeResult(18), tolerance);
    ASSERT_NEAR(-4.13165265000000, averaging.GetEdgeResult(19), tolerance);
}

TEST(Averaging, InterpolateOnNodesClosestPoint)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/AveragingInterpolationTests/sample_grid_coarse_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), static_cast<meshkernel::Index>(0));

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(*mesh,
                                                 samples,
                                                 meshkernel::AveragingInterpolation::Method::Closest,
                                                 meshkernel::Mesh::Location::Nodes,
                                                 1.01,
                                                 false,
                                                 false,
                                                 1);
    averaging.Compute();

    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(-2.0868346999999998, averaging.GetNodeResult(0), tolerance);
    ASSERT_NEAR(-3.0392157000000002, averaging.GetNodeResult(1), tolerance);
    ASSERT_NEAR(-3.9915965999999998, averaging.GetNodeResult(2), tolerance);
    ASSERT_NEAR(-2.8711484999999999, averaging.GetNodeResult(3), tolerance);
    ASSERT_NEAR(-3.8235294000000000, averaging.GetNodeResult(4), tolerance);
    ASSERT_NEAR(-4.7759103999999999, averaging.GetNodeResult(5), tolerance);
    ASSERT_NEAR(-3.6554622000000001, averaging.GetNodeResult(6), tolerance);
    ASSERT_NEAR(-4.6078431000000002, averaging.GetNodeResult(7), tolerance);
    ASSERT_NEAR(-5.5602241000000001, averaging.GetNodeResult(8), tolerance);
    ASSERT_NEAR(-4.4397758999999999, averaging.GetNodeResult(9), tolerance);
    ASSERT_NEAR(-5.3921568999999998, averaging.GetNodeResult(10), tolerance);
    ASSERT_NEAR(-6.3445378000000003, averaging.GetNodeResult(11), tolerance);
    ASSERT_NEAR(-5.2240896000000001, averaging.GetNodeResult(12), tolerance);
    ASSERT_NEAR(-6.1764706000000000, averaging.GetNodeResult(13), tolerance);
    ASSERT_NEAR(-7.1288514999999997, averaging.GetNodeResult(14), tolerance);
    ASSERT_NEAR(-6.0084033999999997, averaging.GetNodeResult(15), tolerance);
    ASSERT_NEAR(-6.9607843000000003, averaging.GetNodeResult(16), tolerance);
    ASSERT_NEAR(-7.9131653000000002, averaging.GetNodeResult(17), tolerance);
}

TEST(Averaging, InterpolateOnFacesClosestPoint)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/AveragingInterpolationTests/sample_grid_coarse_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), static_cast<meshkernel::Index>(0));

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(*mesh,
                                                 samples,
                                                 meshkernel::AveragingInterpolation::Method::Closest,
                                                 meshkernel::Mesh::Location::Faces,
                                                 1.01,
                                                 false,
                                                 false,
                                                 1);
    averaging.Compute();

    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(-2.9551821000000000, averaging.GetFaceResult(0), tolerance);
    ASSERT_NEAR(-3.9075630000000001, averaging.GetFaceResult(1), tolerance);
    ASSERT_NEAR(-3.7394957999999998, averaging.GetFaceResult(2), tolerance);
    ASSERT_NEAR(-4.6918768000000002, averaging.GetFaceResult(3), tolerance);
    ASSERT_NEAR(-4.5238094999999996, averaging.GetFaceResult(4), tolerance);
    ASSERT_NEAR(-5.4761905000000004, averaging.GetFaceResult(5), tolerance);
    ASSERT_NEAR(-5.3081231999999998, averaging.GetFaceResult(6), tolerance);
    ASSERT_NEAR(-6.2605041999999997, averaging.GetFaceResult(7), tolerance);
    ASSERT_NEAR(-6.0924370000000003, averaging.GetFaceResult(8), tolerance);
    ASSERT_NEAR(-7.0448179000000000, averaging.GetFaceResult(9), tolerance);
}

// Max
TEST(Averaging, InterpolateOnEdgesMax)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/AveragingInterpolationTests/sample_grid_coarse_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), static_cast<meshkernel::Index>(0));

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(*mesh,
                                                 samples,
                                                 meshkernel::AveragingInterpolation::Method::Max,
                                                 meshkernel::Mesh::Location::Edges,
                                                 1.01,
                                                 false,
                                                 false,
                                                 1);
    averaging.Compute();

    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(-2.28291315000000, averaging.GetEdgeResult(0), tolerance);
    ASSERT_NEAR(-2.75910365000000, averaging.GetEdgeResult(1), tolerance);
    ASSERT_NEAR(-3.71148460000000, averaging.GetEdgeResult(2), tolerance);
    ASSERT_NEAR(-2.87114845000000, averaging.GetEdgeResult(3), tolerance);
    ASSERT_NEAR(-3.34733895000000, averaging.GetEdgeResult(4), tolerance);
    ASSERT_NEAR(-4.29971990000000, averaging.GetEdgeResult(5), tolerance);
    ASSERT_NEAR(-3.65546215000000, averaging.GetEdgeResult(6), tolerance);
    ASSERT_NEAR(-4.13165265000000, averaging.GetEdgeResult(7), tolerance);
    ASSERT_NEAR(-5.08403365000000, averaging.GetEdgeResult(8), tolerance);
    ASSERT_NEAR(-4.43977590000000, averaging.GetEdgeResult(9), tolerance);
    ASSERT_NEAR(-4.91596635000000, averaging.GetEdgeResult(10), tolerance);
    ASSERT_NEAR(-5.86834735000000, averaging.GetEdgeResult(11), tolerance);
    ASSERT_NEAR(-5.22408965000000, averaging.GetEdgeResult(12), tolerance);
    ASSERT_NEAR(-5.70028010000000, averaging.GetEdgeResult(13), tolerance);
    ASSERT_NEAR(-6.65266105000000, averaging.GetEdgeResult(14), tolerance);
    ASSERT_NEAR(-2.32492995000000, averaging.GetEdgeResult(15), tolerance);
    ASSERT_NEAR(-3.03921570000000, averaging.GetEdgeResult(16), tolerance);
    ASSERT_NEAR(-2.71708685000000, averaging.GetEdgeResult(17), tolerance);
    ASSERT_NEAR(-3.43137255000000, averaging.GetEdgeResult(18), tolerance);
    ASSERT_NEAR(-3.50140055000000, averaging.GetEdgeResult(19), tolerance);
}

TEST(Averaging, InterpolateOnNodesMax)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/AveragingInterpolationTests/sample_grid_coarse_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), static_cast<meshkernel::Index>(0));

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(*mesh,
                                                 samples,
                                                 meshkernel::AveragingInterpolation::Method::Max,
                                                 meshkernel::Mesh::Location::Nodes,
                                                 1.01,
                                                 false,
                                                 false,
                                                 1);
    averaging.Compute();

    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(-2.0868346999999998, averaging.GetNodeResult(0), tolerance);
    ASSERT_NEAR(-2.5630251999999998, averaging.GetNodeResult(1), tolerance);
    ASSERT_NEAR(-3.5154062000000001, averaging.GetNodeResult(2), tolerance);
    ASSERT_NEAR(-2.4789916000000001, averaging.GetNodeResult(3), tolerance);
    ASSERT_NEAR(-2.9551821000000000, averaging.GetNodeResult(4), tolerance);
    ASSERT_NEAR(-3.9075630000000001, averaging.GetNodeResult(5), tolerance);
    ASSERT_NEAR(-3.2633052999999999, averaging.GetNodeResult(6), tolerance);
    ASSERT_NEAR(-3.7394957999999998, averaging.GetNodeResult(7), tolerance);
    ASSERT_NEAR(-4.6918768000000002, averaging.GetNodeResult(8), tolerance);
    ASSERT_NEAR(-4.0476190000000001, averaging.GetNodeResult(9), tolerance);
    ASSERT_NEAR(-4.5238094999999996, averaging.GetNodeResult(10), tolerance);
    ASSERT_NEAR(-5.4761905000000004, averaging.GetNodeResult(11), tolerance);
    ASSERT_NEAR(-4.8319327999999997, averaging.GetNodeResult(12), tolerance);
    ASSERT_NEAR(-5.3081231999999998, averaging.GetNodeResult(13), tolerance);
    ASSERT_NEAR(-6.2605041999999997, averaging.GetNodeResult(14), tolerance);
    ASSERT_NEAR(-5.6162464999999999, averaging.GetNodeResult(15), tolerance);
    ASSERT_NEAR(-6.0924370000000003, averaging.GetNodeResult(16), tolerance);
    ASSERT_NEAR(-7.0448179000000000, averaging.GetNodeResult(17), tolerance);
}

TEST(Averaging, InterpolateOnFacesMax)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/AveragingInterpolationTests/sample_grid_coarse_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), static_cast<meshkernel::Index>(0));

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(*mesh,
                                                 samples,
                                                 meshkernel::AveragingInterpolation::Method::Max,
                                                 meshkernel::Mesh::Location::Faces,
                                                 1.01,
                                                 false,
                                                 false,
                                                 1);
    averaging.Compute();

    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(-2.0868346999999998, averaging.GetFaceResult(0), tolerance);
    ASSERT_NEAR(-3.0392157000000002, averaging.GetFaceResult(1), tolerance);
    ASSERT_NEAR(-2.8711484999999999, averaging.GetFaceResult(2), tolerance);
    ASSERT_NEAR(-3.8235294000000000, averaging.GetFaceResult(3), tolerance);
    ASSERT_NEAR(-3.6554622000000001, averaging.GetFaceResult(4), tolerance);
    ASSERT_NEAR(-4.6078431000000002, averaging.GetFaceResult(5), tolerance);
    ASSERT_NEAR(-4.4397758999999999, averaging.GetFaceResult(6), tolerance);
    ASSERT_NEAR(-5.3921568999999998, averaging.GetFaceResult(7), tolerance);
    ASSERT_NEAR(-5.2240896000000001, averaging.GetFaceResult(8), tolerance);
    ASSERT_NEAR(-6.1764706000000000, averaging.GetFaceResult(9), tolerance);
}

// Min
TEST(Averaging, InterpolateOnEdgesMin)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/AveragingInterpolationTests/sample_grid_coarse_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), static_cast<meshkernel::Index>(0));

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(*mesh,
                                                 samples,
                                                 meshkernel::AveragingInterpolation::Method::Min,
                                                 meshkernel::Mesh::Location::Edges,
                                                 1.01,
                                                 false,
                                                 false,
                                                 1);
    averaging.Compute();

    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(-3.34733895000000, averaging.GetEdgeResult(0), tolerance);
    ASSERT_NEAR(-4.29971990000000, averaging.GetEdgeResult(1), tolerance);
    ASSERT_NEAR(-4.77591035000000, averaging.GetEdgeResult(2), tolerance);
    ASSERT_NEAR(-4.13165265000000, averaging.GetEdgeResult(3), tolerance);
    ASSERT_NEAR(-5.08403365000000, averaging.GetEdgeResult(4), tolerance);
    ASSERT_NEAR(-5.56022410000000, averaging.GetEdgeResult(5), tolerance);
    ASSERT_NEAR(-4.91596635000000, averaging.GetEdgeResult(6), tolerance);
    ASSERT_NEAR(-5.86834735000000, averaging.GetEdgeResult(7), tolerance);
    ASSERT_NEAR(-6.34453785000000, averaging.GetEdgeResult(8), tolerance);
    ASSERT_NEAR(-5.70028010000000, averaging.GetEdgeResult(9), tolerance);
    ASSERT_NEAR(-6.65266105000000, averaging.GetEdgeResult(10), tolerance);
    ASSERT_NEAR(-7.12885155000000, averaging.GetEdgeResult(11), tolerance);
    ASSERT_NEAR(-6.28851540000000, averaging.GetEdgeResult(12), tolerance);
    ASSERT_NEAR(-7.24089635000000, averaging.GetEdgeResult(13), tolerance);
    ASSERT_NEAR(-7.71708685000000, averaging.GetEdgeResult(14), tolerance);
    ASSERT_NEAR(-3.43137255000000, averaging.GetEdgeResult(15), tolerance);
    ASSERT_NEAR(-4.14565825000000, averaging.GetEdgeResult(16), tolerance);
    ASSERT_NEAR(-4.21568630000000, averaging.GetEdgeResult(17), tolerance);
    ASSERT_NEAR(-4.92997200000000, averaging.GetEdgeResult(18), tolerance);
    ASSERT_NEAR(-5.00000000000000, averaging.GetEdgeResult(19), tolerance);
}

TEST(Averaging, InterpolateOnNodesMin)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/AveragingInterpolationTests/sample_grid_coarse_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), static_cast<meshkernel::Index>(0));

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(*mesh,
                                                 samples,
                                                 meshkernel::AveragingInterpolation::Method::Min,
                                                 meshkernel::Mesh::Location::Nodes,
                                                 1.01,
                                                 false,
                                                 false,
                                                 1);
    averaging.Compute();

    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(-2.9551821000000000, averaging.GetNodeResult(0), tolerance);
    ASSERT_NEAR(-3.9075630000000001, averaging.GetNodeResult(1), tolerance);
    ASSERT_NEAR(-4.3837535000000001, averaging.GetNodeResult(2), tolerance);
    ASSERT_NEAR(-3.7394957999999998, averaging.GetNodeResult(3), tolerance);
    ASSERT_NEAR(-4.6918768000000002, averaging.GetNodeResult(4), tolerance);
    ASSERT_NEAR(-5.1680672000000003, averaging.GetNodeResult(5), tolerance);
    ASSERT_NEAR(-4.5238094999999996, averaging.GetNodeResult(6), tolerance);
    ASSERT_NEAR(-5.4761905000000004, averaging.GetNodeResult(7), tolerance);
    ASSERT_NEAR(-5.9523809999999999, averaging.GetNodeResult(8), tolerance);
    ASSERT_NEAR(-5.3081231999999998, averaging.GetNodeResult(9), tolerance);
    ASSERT_NEAR(-6.2605041999999997, averaging.GetNodeResult(10), tolerance);
    ASSERT_NEAR(-6.7366947000000001, averaging.GetNodeResult(11), tolerance);
    ASSERT_NEAR(-6.0924370000000003, averaging.GetNodeResult(12), tolerance);
    ASSERT_NEAR(-7.0448179000000000, averaging.GetNodeResult(13), tolerance);
    ASSERT_NEAR(-7.5210084000000004, averaging.GetNodeResult(14), tolerance);
    ASSERT_NEAR(-6.4845937999999999, averaging.GetNodeResult(15), tolerance);
    ASSERT_NEAR(-7.4369747999999998, averaging.GetNodeResult(16), tolerance);
    ASSERT_NEAR(-7.9131653000000002, averaging.GetNodeResult(17), tolerance);
}

TEST(Averaging, InterpolateOnFacesMin)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/AveragingInterpolationTests/sample_grid_coarse_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), static_cast<meshkernel::Index>(0));

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(*mesh,
                                                 samples,
                                                 meshkernel::AveragingInterpolation::Method::Min,
                                                 meshkernel::Mesh::Location::Faces,
                                                 1.01,
                                                 false,
                                                 false,
                                                 1);
    averaging.Compute();

    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(-3.8235294000000000, averaging.GetFaceResult(0), tolerance);
    ASSERT_NEAR(-4.7759103999999999, averaging.GetFaceResult(1), tolerance);
    ASSERT_NEAR(-4.6078431000000002, averaging.GetFaceResult(2), tolerance);
    ASSERT_NEAR(-5.5602241000000001, averaging.GetFaceResult(3), tolerance);
    ASSERT_NEAR(-5.3921568999999998, averaging.GetFaceResult(4), tolerance);
    ASSERT_NEAR(-6.3445378000000003, averaging.GetFaceResult(5), tolerance);
    ASSERT_NEAR(-6.1764706000000000, averaging.GetFaceResult(6), tolerance);
    ASSERT_NEAR(-7.1288514999999997, averaging.GetFaceResult(7), tolerance);
    ASSERT_NEAR(-6.9607843000000003, averaging.GetFaceResult(8), tolerance);
    ASSERT_NEAR(-7.9131653000000002, averaging.GetFaceResult(9), tolerance);
}

// InverseWeightedDistance
TEST(Averaging, InterpolateOnEdgesInverseWeightedDistance)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/AveragingInterpolationTests/sample_grid_coarse_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), static_cast<meshkernel::Index>(0));

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(*mesh,
                                                 samples,
                                                 meshkernel::AveragingInterpolation::Method::InverseWeightedDistance,
                                                 meshkernel::Mesh::Location::Edges,
                                                 1.01,
                                                 false,
                                                 false,
                                                 1);
    averaging.Compute();

    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(-2.48030306628047, averaging.GetEdgeResult(0), tolerance);
    ASSERT_NEAR(-3.43184384480726, averaging.GetEdgeResult(1), tolerance);
    ASSERT_NEAR(-4.38310967989827, averaging.GetEdgeResult(2), tolerance);
    ASSERT_NEAR(-3.26444992295321, averaging.GetEdgeResult(3), tolerance);
    ASSERT_NEAR(-4.21568625020467, averaging.GetEdgeResult(4), tolerance);
    ASSERT_NEAR(-5.16692267681247, averaging.GetEdgeResult(5), tolerance);
    ASSERT_NEAR(-4.04876362305277, averaging.GetEdgeResult(6), tolerance);
    ASSERT_NEAR(-5.00000000000000, averaging.GetEdgeResult(7), tolerance);
    ASSERT_NEAR(-5.95123637694723, averaging.GetEdgeResult(8), tolerance);
    ASSERT_NEAR(-4.83307732318753, averaging.GetEdgeResult(9), tolerance);
    ASSERT_NEAR(-5.78431374979533, averaging.GetEdgeResult(10), tolerance);
    ASSERT_NEAR(-6.73555007704679, averaging.GetEdgeResult(11), tolerance);
    ASSERT_NEAR(-5.61689032010173, averaging.GetEdgeResult(12), tolerance);
    ASSERT_NEAR(-6.56815615519274, averaging.GetEdgeResult(13), tolerance);
    ASSERT_NEAR(-7.51969693371953, averaging.GetEdgeResult(14), tolerance);
    ASSERT_NEAR(-2.56423567460090, averaging.GetEdgeResult(15), tolerance);
    ASSERT_NEAR(-3.51580591130667, averaging.GetEdgeResult(16), tolerance);
    ASSERT_NEAR(-3.34791123648683, averaging.GetEdgeResult(17), tolerance);
    ASSERT_NEAR(-4.29914761339886, averaging.GetEdgeResult(18), tolerance);
    ASSERT_NEAR(-4.13222493667105, averaging.GetEdgeResult(19), tolerance);
}

TEST(Averaging, InterpolateOnNodesInverseWeightedDistance)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/AveragingInterpolationTests/sample_grid_coarse_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), static_cast<meshkernel::Index>(0));

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(*mesh,
                                                 samples,
                                                 meshkernel::AveragingInterpolation::Method::InverseWeightedDistance,
                                                 meshkernel::Mesh::Location::Nodes,
                                                 1.01,
                                                 false,
                                                 false,
                                                 1);
    averaging.Compute();

    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(-2.0883130596575117, averaging.GetNodeResult(0), tolerance);
    ASSERT_NEAR(-3.0401582895442973, averaging.GetNodeResult(1), tolerance);
    ASSERT_NEAR(-3.9914535330690506, averaging.GetNodeResult(2), tolerance);
    ASSERT_NEAR(-2.8722930729034286, averaging.GetNodeResult(3), tolerance);
    ASSERT_NEAR(-3.8235294000702318, averaging.GetNodeResult(4), tolerance);
    ASSERT_NEAR(-4.7747658267274895, averaging.GetNodeResult(5), tolerance);
    ASSERT_NEAR(-3.6566067730029892, averaging.GetNodeResult(6), tolerance);
    ASSERT_NEAR(-4.6078431003391067, averaging.GetNodeResult(7), tolerance);
    ASSERT_NEAR(-5.5590795268974498, averaging.GetNodeResult(8), tolerance);
    ASSERT_NEAR(-4.4409204731025502, averaging.GetNodeResult(9), tolerance);
    ASSERT_NEAR(-5.3921568996608942, averaging.GetNodeResult(10), tolerance);
    ASSERT_NEAR(-6.3433932269970112, averaging.GetNodeResult(11), tolerance);
    ASSERT_NEAR(-5.2252341732725114, averaging.GetNodeResult(12), tolerance);
    ASSERT_NEAR(-6.1764705999297691, averaging.GetNodeResult(13), tolerance);
    ASSERT_NEAR(-7.1277069270965718, averaging.GetNodeResult(14), tolerance);
    ASSERT_NEAR(-6.0085464669309490, averaging.GetNodeResult(15), tolerance);
    ASSERT_NEAR(-6.9598417104557040, averaging.GetNodeResult(16), tolerance);
    ASSERT_NEAR(-7.9116869403424879, averaging.GetNodeResult(17), tolerance);
}

TEST(Averaging, InterpolateOnFacesInverseWeightedDistance)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/AveragingInterpolationTests/sample_grid_coarse_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), static_cast<meshkernel::Index>(0));

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(*mesh,
                                                 samples,
                                                 meshkernel::AveragingInterpolation::Method::InverseWeightedDistance,
                                                 meshkernel::Mesh::Location::Faces,
                                                 1.01,
                                                 false,
                                                 false,
                                                 1);
    averaging.Compute();

    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(-2.9551820997311253, averaging.GetFaceResult(0), tolerance);
    ASSERT_NEAR(-3.9075630001695529, averaging.GetFaceResult(1), tolerance);
    ASSERT_NEAR(-3.7394957999999998, averaging.GetFaceResult(2), tolerance);
    ASSERT_NEAR(-4.6918767996608945, averaging.GetFaceResult(3), tolerance);
    ASSERT_NEAR(-4.5238095001695529, averaging.GetFaceResult(4), tolerance);
    ASSERT_NEAR(-5.4761904998304463, averaging.GetFaceResult(5), tolerance);
    ASSERT_NEAR(-5.3081232003391055, averaging.GetFaceResult(6), tolerance);
    ASSERT_NEAR(-6.2605042000000006, averaging.GetFaceResult(7), tolerance);
    ASSERT_NEAR(-6.0924369998304471, averaging.GetFaceResult(8), tolerance);
    ASSERT_NEAR(-7.0448179002688773, averaging.GetFaceResult(9), tolerance);
}

// MinAbsValue
TEST(Averaging, InterpolateOnEdgesMinAbsValue)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/AveragingInterpolationTests/sample_grid_coarse_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), static_cast<meshkernel::Index>(0));

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(*mesh,
                                                 samples,
                                                 meshkernel::AveragingInterpolation::Method::MinAbsValue,
                                                 meshkernel::Mesh::Location::Edges,
                                                 1.01,
                                                 false,
                                                 false,
                                                 1);
    averaging.Compute();

    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(2.28291315000000, averaging.GetEdgeResult(0), tolerance);
    ASSERT_NEAR(2.75910365000000, averaging.GetEdgeResult(1), tolerance);
    ASSERT_NEAR(3.71148460000000, averaging.GetEdgeResult(2), tolerance);
    ASSERT_NEAR(2.87114845000000, averaging.GetEdgeResult(3), tolerance);
    ASSERT_NEAR(3.34733895000000, averaging.GetEdgeResult(4), tolerance);
    ASSERT_NEAR(4.29971990000000, averaging.GetEdgeResult(5), tolerance);
    ASSERT_NEAR(3.65546215000000, averaging.GetEdgeResult(6), tolerance);
    ASSERT_NEAR(4.13165265000000, averaging.GetEdgeResult(7), tolerance);
    ASSERT_NEAR(5.08403365000000, averaging.GetEdgeResult(8), tolerance);
    ASSERT_NEAR(4.43977590000000, averaging.GetEdgeResult(9), tolerance);
    ASSERT_NEAR(4.91596635000000, averaging.GetEdgeResult(10), tolerance);
    ASSERT_NEAR(5.86834735000000, averaging.GetEdgeResult(11), tolerance);
    ASSERT_NEAR(5.22408965000000, averaging.GetEdgeResult(12), tolerance);
    ASSERT_NEAR(5.70028010000000, averaging.GetEdgeResult(13), tolerance);
    ASSERT_NEAR(6.65266105000000, averaging.GetEdgeResult(14), tolerance);
    ASSERT_NEAR(2.32492995000000, averaging.GetEdgeResult(15), tolerance);
    ASSERT_NEAR(3.03921570000000, averaging.GetEdgeResult(16), tolerance);
    ASSERT_NEAR(2.71708685000000, averaging.GetEdgeResult(17), tolerance);
    ASSERT_NEAR(3.43137255000000, averaging.GetEdgeResult(18), tolerance);
    ASSERT_NEAR(3.50140055000000, averaging.GetEdgeResult(19), tolerance);
}

TEST(Averaging, InterpolateOnNodesMinAbsValue)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/AveragingInterpolationTests/sample_grid_coarse_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), static_cast<meshkernel::Index>(0));

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(*mesh,
                                                 samples,
                                                 meshkernel::AveragingInterpolation::Method::MinAbsValue,
                                                 meshkernel::Mesh::Location::Nodes,
                                                 1.01,
                                                 false,
                                                 false,
                                                 1);
    averaging.Compute();

    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(2.0868346999999998, averaging.GetNodeResult(0), tolerance);
    ASSERT_NEAR(2.5630251999999998, averaging.GetNodeResult(1), tolerance);
    ASSERT_NEAR(3.5154062000000001, averaging.GetNodeResult(2), tolerance);
    ASSERT_NEAR(2.4789916000000001, averaging.GetNodeResult(3), tolerance);
    ASSERT_NEAR(2.9551821000000000, averaging.GetNodeResult(4), tolerance);
    ASSERT_NEAR(3.9075630000000001, averaging.GetNodeResult(5), tolerance);
    ASSERT_NEAR(3.2633052999999999, averaging.GetNodeResult(6), tolerance);
    ASSERT_NEAR(3.7394957999999998, averaging.GetNodeResult(7), tolerance);
    ASSERT_NEAR(4.6918768000000002, averaging.GetNodeResult(8), tolerance);
    ASSERT_NEAR(4.0476190000000001, averaging.GetNodeResult(9), tolerance);
    ASSERT_NEAR(4.5238094999999996, averaging.GetNodeResult(10), tolerance);
    ASSERT_NEAR(5.4761905000000004, averaging.GetNodeResult(11), tolerance);
    ASSERT_NEAR(4.8319327999999997, averaging.GetNodeResult(12), tolerance);
    ASSERT_NEAR(5.3081231999999998, averaging.GetNodeResult(13), tolerance);
    ASSERT_NEAR(6.2605041999999997, averaging.GetNodeResult(14), tolerance);
    ASSERT_NEAR(5.6162464999999999, averaging.GetNodeResult(15), tolerance);
    ASSERT_NEAR(6.0924370000000003, averaging.GetNodeResult(16), tolerance);
    ASSERT_NEAR(7.0448179000000000, averaging.GetNodeResult(17), tolerance);
}

TEST(Averaging, Interpolate_WithSamplesOnEdges_ShouldNotExtendSampleCoverage)
{
    // 1 Setup
    auto mesh = MakeRectangularMeshForTesting(5, 5, 1.0, meshkernel::Projection::cartesian);

    std::vector<meshkernel::Sample> samples{
        {2.5, 1.0, 1.0},
        {2.5, 2.0, 1.0},
        {2.5, 3.0, 1.0},
        {3.5, 1.0, 1.0},
        {3.5, 2.0, 1.0},
        {3.5, 3.0, 1.0},
        {2.0, 1.5, 1.0},
        {2.0, 2.5, 1.0},
        {3.0, 1.5, 1.0},
        {3.0, 2.5, 1.0},
        {4.0, 1.5, 1.0},
        {4.0, 2.5, 1.0}};

    ASSERT_GT(mesh->GetNumNodes(), static_cast<meshkernel::Index>(0));

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(*mesh,
                                                 samples,
                                                 meshkernel::AveragingInterpolation::Method::SimpleAveraging,
                                                 meshkernel::Mesh::Location::Edges,
                                                 1.01,
                                                 false,
                                                 false,
                                                 1);
    averaging.Compute();

    // Assert: only 12 edges gets a valid value, as the provided sample set
    const auto& interpolationResults = averaging.GetEdgeResults();
    std::vector expectedInterpolationResults(interpolationResults.size(), meshkernel::constants::missing::doubleValue);
    expectedInterpolationResults[11] = 1.0;
    expectedInterpolationResults[12] = 1.0;
    expectedInterpolationResults[13] = 1.0;
    expectedInterpolationResults[16] = 1.0;
    expectedInterpolationResults[17] = 1.0;
    expectedInterpolationResults[18] = 1.0;
    expectedInterpolationResults[29] = 1.0;
    expectedInterpolationResults[30] = 1.0;
    expectedInterpolationResults[33] = 1.0;
    expectedInterpolationResults[34] = 1.0;
    expectedInterpolationResults[37] = 1.0;
    expectedInterpolationResults[38] = 1.0;

    ASSERT_THAT(interpolationResults, ::testing::ContainerEq(expectedInterpolationResults));
}

TEST(Averaging, Interpolate_WithSamplesOnNodes_ShouldNotExtendSampleCoverage)
{
    // 1 Setup
    auto mesh = MakeRectangularMeshForTesting(5, 5, 1.0, meshkernel::Projection::cartesian);

    std::vector<meshkernel::Sample> samples{
        {2.0, 1.0, 1.0},
        {2.0, 2.0, 1.0},
        {2.0, 3.0, 1.0},
        {3.0, 1.0, 1.0},
        {3.0, 2.0, 1.0},
        {3.0, 3.0, 1.0},
        {4.0, 1.0, 1.0},
        {4.0, 2.0, 1.0},
        {4.0, 3.0, 1.0}};

    ASSERT_GT(mesh->GetNumNodes(), static_cast<meshkernel::Index>(0));

    // Execute averaging, this time the minimum number of samples should at least be 2
    meshkernel::AveragingInterpolation averaging(*mesh,
                                                 samples,
                                                 meshkernel::AveragingInterpolation::Method::SimpleAveraging,
                                                 meshkernel::Mesh::Location::Edges,
                                                 1.01,
                                                 false,
                                                 false,
                                                 1);
    averaging.Compute();

    // Assert: only 12 edges gets a valid value
    const auto& interpolationResults = averaging.GetEdgeResults();
    std::vector expectedInterpolationResults(interpolationResults.size(), meshkernel::constants::missing::doubleValue);
    expectedInterpolationResults[11] = 1.0;
    expectedInterpolationResults[12] = 1.0;
    expectedInterpolationResults[13] = 1.0;
    expectedInterpolationResults[16] = 1.0;
    expectedInterpolationResults[17] = 1.0;
    expectedInterpolationResults[18] = 1.0;
    expectedInterpolationResults[29] = 1.0;
    expectedInterpolationResults[30] = 1.0;
    expectedInterpolationResults[33] = 1.0;
    expectedInterpolationResults[34] = 1.0;
    expectedInterpolationResults[37] = 1.0;
    expectedInterpolationResults[38] = 1.0;

    ASSERT_THAT(interpolationResults, ::testing::ContainerEq(expectedInterpolationResults));
}
