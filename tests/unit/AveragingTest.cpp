#include <MeshKernel/AveragingInterpolation.hpp>
#include <TestUtils/SampleFileReader.hpp>
#include <TestUtils/MakeMeshes.hpp>
#include <gtest/gtest.h>

//Simple averaging
TEST(Averaging, InterpolateOnEdgesSimpleAveraging)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile("../../../../tests/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/AveragingInterpolationTests/simple_grid_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), 0);

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(mesh, samples, meshkernel::AveragingInterpolation::Method::SimpleAveraging, meshkernel::InterpolationLocation::Edges, 1.01, false, false);
    averaging.Compute();

    constexpr double tolerance = 1e-6;
    auto interpolationResults = averaging.GetResults();
    ASSERT_NEAR(-2.28291315000000, interpolationResults[0], tolerance);
    ASSERT_NEAR(-2.75910365000000, interpolationResults[1], tolerance);
    ASSERT_NEAR(-3.23529410000000, interpolationResults[2], tolerance);
    ASSERT_NEAR(-3.71148460000000, interpolationResults[3], tolerance);
    ASSERT_NEAR(-4.18767505000000, interpolationResults[4], tolerance);
    ASSERT_NEAR(-2.67507005000000, interpolationResults[5], tolerance);
    ASSERT_NEAR(-3.15126050000000, interpolationResults[6], tolerance);
    ASSERT_NEAR(-3.62745095000000, interpolationResults[7], tolerance);
    ASSERT_NEAR(-4.10364145000000, interpolationResults[8], tolerance);
    ASSERT_NEAR(-4.57983195000000, interpolationResults[9], tolerance);
    ASSERT_NEAR(-3.06722690000000, interpolationResults[10], tolerance);
    ASSERT_NEAR(-3.54341735000000, interpolationResults[11], tolerance);
    ASSERT_NEAR(-4.01960785000000, interpolationResults[12], tolerance);
    ASSERT_NEAR(-4.49579835000000, interpolationResults[13], tolerance);
    ASSERT_NEAR(-4.97198880000000, interpolationResults[14], tolerance);
    ASSERT_NEAR(-3.45938375000000, interpolationResults[15], tolerance);
    ASSERT_NEAR(-3.93557425000000, interpolationResults[16], tolerance);
    ASSERT_NEAR(-4.41176470000000, interpolationResults[17], tolerance);
    ASSERT_NEAR(-4.88795520000000, interpolationResults[18], tolerance);
    ASSERT_NEAR(-5.36414565000000, interpolationResults[19], tolerance);
}

TEST(Averaging, InterpolateOnNodesSimpleAveraging)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile("../../../../tests/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/AveragingInterpolationTests/simple_grid_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), 0);

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(mesh, samples, meshkernel::AveragingInterpolation::Method::SimpleAveraging, meshkernel::InterpolationLocation::Nodes, 1.01, false, false);
    averaging.Compute();

    constexpr double tolerance = 1e-6;
    auto interpolationResults = averaging.GetResults();
    ASSERT_NEAR(-2.08683470000000, interpolationResults[0], tolerance);
    ASSERT_NEAR(-2.47899160000000, interpolationResults[1], tolerance);
    ASSERT_NEAR(-2.56302520000000, interpolationResults[2], tolerance);
    ASSERT_NEAR(-2.95518210000000, interpolationResults[3], tolerance);
    ASSERT_NEAR(-3.03921570000000, interpolationResults[4], tolerance);
    ASSERT_NEAR(-3.43137250000000, interpolationResults[5], tolerance);
    ASSERT_NEAR(-3.51540620000000, interpolationResults[6], tolerance);
    ASSERT_NEAR(-3.90756300000000, interpolationResults[7], tolerance);
    ASSERT_NEAR(-3.99159660000000, interpolationResults[8], tolerance);
    ASSERT_NEAR(-4.38375350000000, interpolationResults[9], tolerance);
    ASSERT_NEAR(-2.87114850000000, interpolationResults[10], tolerance);
    ASSERT_NEAR(-3.34733890000000, interpolationResults[11], tolerance);
    ASSERT_NEAR(-3.82352940000000, interpolationResults[12], tolerance);
    ASSERT_NEAR(-4.29971990000000, interpolationResults[13], tolerance);
    ASSERT_NEAR(-4.77591040000000, interpolationResults[14], tolerance);
    ASSERT_NEAR(-3.26330530000000, interpolationResults[15], tolerance);
    ASSERT_NEAR(-3.73949580000000, interpolationResults[16], tolerance);
    ASSERT_NEAR(-4.21568630000000, interpolationResults[17], tolerance);
    ASSERT_NEAR(-4.69187680000000, interpolationResults[18], tolerance);
    ASSERT_NEAR(-5.16806720000000, interpolationResults[19], tolerance);
}

TEST(Averaging, InterpolateOnFacesSimpleAveraging)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile("../../../../tests/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/AveragingInterpolationTests/simple_grid_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), 0);

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(mesh, samples, meshkernel::AveragingInterpolation::Method::SimpleAveraging, meshkernel::InterpolationLocation::Faces, 1.01, false, false);
    averaging.Compute();

    constexpr double tolerance = 1e-6;
    auto interpolationResults = averaging.GetResults();

    ASSERT_NEAR(-2.52100840000000, interpolationResults[0], tolerance);
    ASSERT_NEAR(-2.91316527500000, interpolationResults[1], tolerance);
    ASSERT_NEAR(-2.99719887500000, interpolationResults[2], tolerance);
    ASSERT_NEAR(-3.38935572500000, interpolationResults[3], tolerance);
    ASSERT_NEAR(-3.47338935000000, interpolationResults[4], tolerance);
    ASSERT_NEAR(-3.86554620000000, interpolationResults[5], tolerance);
    ASSERT_NEAR(-3.94957982500000, interpolationResults[6], tolerance);
    ASSERT_NEAR(-4.34173670000000, interpolationResults[7], tolerance);
    ASSERT_NEAR(-3.30532212500000, interpolationResults[8], tolerance);
    ASSERT_NEAR(-3.78151260000000, interpolationResults[9], tolerance);
    ASSERT_NEAR(-4.25770310000000, interpolationResults[10], tolerance);
    ASSERT_NEAR(-4.73389357500000, interpolationResults[11], tolerance);
    ASSERT_NEAR(-3.69747900000000, interpolationResults[12], tolerance);
    ASSERT_NEAR(-4.17366947500000, interpolationResults[13], tolerance);
    ASSERT_NEAR(-4.64985995000000, interpolationResults[14], tolerance);
    ASSERT_NEAR(-5.12605042500000, interpolationResults[15], tolerance);
    ASSERT_NEAR(-4.08963585000000, interpolationResults[16], tolerance);
    ASSERT_NEAR(-4.56582632500000, interpolationResults[17], tolerance);
    ASSERT_NEAR(-5.04201680000000, interpolationResults[18], tolerance);
    ASSERT_NEAR(-5.51820730000000, interpolationResults[19], tolerance);
}

//Closest Point
TEST(Averaging, InterpolateOnEdgesClosestPoint)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile("../../../../tests/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/AveragingInterpolationTests/sample_grid_coarse_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), 0);

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(mesh, samples, meshkernel::AveragingInterpolation::Method::Closest, meshkernel::InterpolationLocation::Edges, 1.01, false, false);
    averaging.Compute();

    constexpr double tolerance = 1e-6;
    auto interpolationResults = averaging.GetResults();
    ASSERT_NEAR(-2.47899160000000, interpolationResults[0], tolerance);
    ASSERT_NEAR(-3.43137255000000, interpolationResults[1], tolerance);
    ASSERT_NEAR(-4.38375350000000, interpolationResults[2], tolerance);
    ASSERT_NEAR(-3.26330535000000, interpolationResults[3], tolerance);
    ASSERT_NEAR(-4.21568625000000, interpolationResults[4], tolerance);
    ASSERT_NEAR(-5.16806725000000, interpolationResults[5], tolerance);
    ASSERT_NEAR(-4.04761905000000, interpolationResults[6], tolerance);
    ASSERT_NEAR(-5.00000000000000, interpolationResults[7], tolerance);
    ASSERT_NEAR(-5.95238095000000, interpolationResults[8], tolerance);
    ASSERT_NEAR(-4.83193275000000, interpolationResults[9], tolerance);
    ASSERT_NEAR(-5.78431375000000, interpolationResults[10], tolerance);
    ASSERT_NEAR(-6.73669465000000, interpolationResults[11], tolerance);
    ASSERT_NEAR(-5.61624650000000, interpolationResults[12], tolerance);
    ASSERT_NEAR(-6.56862745000000, interpolationResults[13], tolerance);
    ASSERT_NEAR(-7.52100840000000, interpolationResults[14], tolerance);
    ASSERT_NEAR(-2.56302520000000, interpolationResults[15], tolerance);
    ASSERT_NEAR(-3.51540615000000, interpolationResults[16], tolerance);
    ASSERT_NEAR(-3.34733895000000, interpolationResults[17], tolerance);
    ASSERT_NEAR(-4.29971990000000, interpolationResults[18], tolerance);
    ASSERT_NEAR(-4.13165265000000, interpolationResults[19], tolerance);
}

TEST(Averaging, InterpolateOnNodesClosestPoint)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile("../../../../tests/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/AveragingInterpolationTests/sample_grid_coarse_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), 0);

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(mesh, samples, meshkernel::AveragingInterpolation::Method::Closest, meshkernel::InterpolationLocation::Nodes, 1.01, false, false);
    averaging.Compute();

    constexpr double tolerance = 1e-6;
    auto interpolationResults = averaging.GetResults();
    ASSERT_NEAR(-2.0868346999999998, interpolationResults[0], tolerance);
    ASSERT_NEAR(-3.0392157000000002, interpolationResults[1], tolerance);
    ASSERT_NEAR(-3.9915965999999998, interpolationResults[2], tolerance);
    ASSERT_NEAR(-2.8711484999999999, interpolationResults[3], tolerance);
    ASSERT_NEAR(-3.8235294000000000, interpolationResults[4], tolerance);
    ASSERT_NEAR(-4.7759103999999999, interpolationResults[5], tolerance);
    ASSERT_NEAR(-3.6554622000000001, interpolationResults[6], tolerance);
    ASSERT_NEAR(-4.6078431000000002, interpolationResults[7], tolerance);
    ASSERT_NEAR(-5.5602241000000001, interpolationResults[8], tolerance);
    ASSERT_NEAR(-4.4397758999999999, interpolationResults[9], tolerance);
    ASSERT_NEAR(-5.3921568999999998, interpolationResults[10], tolerance);
    ASSERT_NEAR(-6.3445378000000003, interpolationResults[11], tolerance);
    ASSERT_NEAR(-5.2240896000000001, interpolationResults[12], tolerance);
    ASSERT_NEAR(-6.1764706000000000, interpolationResults[13], tolerance);
    ASSERT_NEAR(-7.1288514999999997, interpolationResults[14], tolerance);
    ASSERT_NEAR(-6.0084033999999997, interpolationResults[15], tolerance);
    ASSERT_NEAR(-6.9607843000000003, interpolationResults[16], tolerance);
    ASSERT_NEAR(-7.9131653000000002, interpolationResults[17], tolerance);
}

TEST(Averaging, InterpolateOnFacesClosestPoint)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile("../../../../tests/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/AveragingInterpolationTests/sample_grid_coarse_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), 0);

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(mesh, samples, meshkernel::AveragingInterpolation::Method::Closest, meshkernel::InterpolationLocation::Faces, 1.01, false, false);
    averaging.Compute();

    constexpr double tolerance = 1e-6;
    auto interpolationResults = averaging.GetResults();
    ASSERT_NEAR(-2.9551821000000000, interpolationResults[0], tolerance);
    ASSERT_NEAR(-3.9075630000000001, interpolationResults[1], tolerance);
    ASSERT_NEAR(-3.7394957999999998, interpolationResults[2], tolerance);
    ASSERT_NEAR(-4.6918768000000002, interpolationResults[3], tolerance);
    ASSERT_NEAR(-4.5238094999999996, interpolationResults[4], tolerance);
    ASSERT_NEAR(-5.4761905000000004, interpolationResults[5], tolerance);
    ASSERT_NEAR(-5.3081231999999998, interpolationResults[6], tolerance);
    ASSERT_NEAR(-6.2605041999999997, interpolationResults[7], tolerance);
    ASSERT_NEAR(-6.0924370000000003, interpolationResults[8], tolerance);
    ASSERT_NEAR(-7.0448179000000000, interpolationResults[9], tolerance);
}

// Max
TEST(Averaging, InterpolateOnEdgesMax)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile("../../../../tests/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/AveragingInterpolationTests/sample_grid_coarse_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), 0);

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(mesh, samples, meshkernel::AveragingInterpolation::Method::Max, meshkernel::InterpolationLocation::Edges, 1.01, false, false);
    averaging.Compute();

    constexpr double tolerance = 1e-6;
    auto interpolationResults = averaging.GetResults();
    ASSERT_NEAR(-2.28291315000000, interpolationResults[0], tolerance);
    ASSERT_NEAR(-2.75910365000000, interpolationResults[1], tolerance);
    ASSERT_NEAR(-3.71148460000000, interpolationResults[2], tolerance);
    ASSERT_NEAR(-2.87114845000000, interpolationResults[3], tolerance);
    ASSERT_NEAR(-3.34733895000000, interpolationResults[4], tolerance);
    ASSERT_NEAR(-4.29971990000000, interpolationResults[5], tolerance);
    ASSERT_NEAR(-3.65546215000000, interpolationResults[6], tolerance);
    ASSERT_NEAR(-4.13165265000000, interpolationResults[7], tolerance);
    ASSERT_NEAR(-5.08403365000000, interpolationResults[8], tolerance);
    ASSERT_NEAR(-4.43977590000000, interpolationResults[9], tolerance);
    ASSERT_NEAR(-4.91596635000000, interpolationResults[10], tolerance);
    ASSERT_NEAR(-5.86834735000000, interpolationResults[11], tolerance);
    ASSERT_NEAR(-5.22408965000000, interpolationResults[12], tolerance);
    ASSERT_NEAR(-5.70028010000000, interpolationResults[13], tolerance);
    ASSERT_NEAR(-6.65266105000000, interpolationResults[14], tolerance);
    ASSERT_NEAR(-2.32492995000000, interpolationResults[15], tolerance);
    ASSERT_NEAR(-3.03921570000000, interpolationResults[16], tolerance);
    ASSERT_NEAR(-2.71708685000000, interpolationResults[17], tolerance);
    ASSERT_NEAR(-3.43137255000000, interpolationResults[18], tolerance);
    ASSERT_NEAR(-3.50140055000000, interpolationResults[19], tolerance);
}

TEST(Averaging, InterpolateOnNodesMax)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile("../../../../tests/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/AveragingInterpolationTests/sample_grid_coarse_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), 0);

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(mesh, samples, meshkernel::AveragingInterpolation::Method::Max, meshkernel::InterpolationLocation::Nodes, 1.01, false, false);
    averaging.Compute();

    constexpr double tolerance = 1e-6;
    auto interpolationResults = averaging.GetResults();
    ASSERT_NEAR(-2.0868346999999998, interpolationResults[0], tolerance);
    ASSERT_NEAR(-2.5630251999999998, interpolationResults[1], tolerance);
    ASSERT_NEAR(-3.5154062000000001, interpolationResults[2], tolerance);
    ASSERT_NEAR(-2.4789916000000001, interpolationResults[3], tolerance);
    ASSERT_NEAR(-2.9551821000000000, interpolationResults[4], tolerance);
    ASSERT_NEAR(-3.9075630000000001, interpolationResults[5], tolerance);
    ASSERT_NEAR(-3.2633052999999999, interpolationResults[6], tolerance);
    ASSERT_NEAR(-3.7394957999999998, interpolationResults[7], tolerance);
    ASSERT_NEAR(-4.6918768000000002, interpolationResults[8], tolerance);
    ASSERT_NEAR(-4.0476190000000001, interpolationResults[9], tolerance);
    ASSERT_NEAR(-4.5238094999999996, interpolationResults[10], tolerance);
    ASSERT_NEAR(-5.4761905000000004, interpolationResults[11], tolerance);
    ASSERT_NEAR(-4.8319327999999997, interpolationResults[12], tolerance);
    ASSERT_NEAR(-5.3081231999999998, interpolationResults[13], tolerance);
    ASSERT_NEAR(-6.2605041999999997, interpolationResults[14], tolerance);
    ASSERT_NEAR(-5.6162464999999999, interpolationResults[15], tolerance);
    ASSERT_NEAR(-6.0924370000000003, interpolationResults[16], tolerance);
    ASSERT_NEAR(-7.0448179000000000, interpolationResults[17], tolerance);
}

TEST(Averaging, InterpolateOnFacesMax)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile("../../../../tests/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/AveragingInterpolationTests/sample_grid_coarse_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), 0);

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(mesh, samples, meshkernel::AveragingInterpolation::Method::Max, meshkernel::InterpolationLocation::Faces, 1.01, false, false);
    averaging.Compute();

    constexpr double tolerance = 1e-6;
    auto interpolationResults = averaging.GetResults();
    ASSERT_NEAR(-2.0868346999999998, interpolationResults[0], tolerance);
    ASSERT_NEAR(-3.0392157000000002, interpolationResults[1], tolerance);
    ASSERT_NEAR(-2.8711484999999999, interpolationResults[2], tolerance);
    ASSERT_NEAR(-3.8235294000000000, interpolationResults[3], tolerance);
    ASSERT_NEAR(-3.6554622000000001, interpolationResults[4], tolerance);
    ASSERT_NEAR(-4.6078431000000002, interpolationResults[5], tolerance);
    ASSERT_NEAR(-4.4397758999999999, interpolationResults[6], tolerance);
    ASSERT_NEAR(-5.3921568999999998, interpolationResults[7], tolerance);
    ASSERT_NEAR(-5.2240896000000001, interpolationResults[8], tolerance);
    ASSERT_NEAR(-6.1764706000000000, interpolationResults[9], tolerance);
}

// Min
TEST(Averaging, InterpolateOnEdgesMin)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile("../../../../tests/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/AveragingInterpolationTests/sample_grid_coarse_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), 0);

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(mesh, samples, meshkernel::AveragingInterpolation::Method::Min, meshkernel::InterpolationLocation::Edges, 1.01, false, false);
    averaging.Compute();

    constexpr double tolerance = 1e-6;
    auto interpolationResults = averaging.GetResults();
    ASSERT_NEAR(-3.34733895000000, interpolationResults[0], tolerance);
    ASSERT_NEAR(-4.29971990000000, interpolationResults[1], tolerance);
    ASSERT_NEAR(-4.77591035000000, interpolationResults[2], tolerance);
    ASSERT_NEAR(-4.13165265000000, interpolationResults[3], tolerance);
    ASSERT_NEAR(-5.08403365000000, interpolationResults[4], tolerance);
    ASSERT_NEAR(-5.56022410000000, interpolationResults[5], tolerance);
    ASSERT_NEAR(-4.91596635000000, interpolationResults[6], tolerance);
    ASSERT_NEAR(-5.86834735000000, interpolationResults[7], tolerance);
    ASSERT_NEAR(-6.34453785000000, interpolationResults[8], tolerance);
    ASSERT_NEAR(-5.70028010000000, interpolationResults[9], tolerance);
    ASSERT_NEAR(-6.65266105000000, interpolationResults[10], tolerance);
    ASSERT_NEAR(-7.12885155000000, interpolationResults[11], tolerance);
    ASSERT_NEAR(-6.28851540000000, interpolationResults[12], tolerance);
    ASSERT_NEAR(-7.24089635000000, interpolationResults[13], tolerance);
    ASSERT_NEAR(-7.71708685000000, interpolationResults[14], tolerance);
    ASSERT_NEAR(-3.43137255000000, interpolationResults[15], tolerance);
    ASSERT_NEAR(-4.14565825000000, interpolationResults[16], tolerance);
    ASSERT_NEAR(-4.21568630000000, interpolationResults[17], tolerance);
    ASSERT_NEAR(-4.92997200000000, interpolationResults[18], tolerance);
    ASSERT_NEAR(-5.00000000000000, interpolationResults[19], tolerance);
}

TEST(Averaging, InterpolateOnNodesMin)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile("../../../../tests/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/AveragingInterpolationTests/sample_grid_coarse_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), 0);

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(mesh, samples, meshkernel::AveragingInterpolation::Method::Min, meshkernel::InterpolationLocation::Nodes, 1.01, false, false);
    averaging.Compute();

    constexpr double tolerance = 1e-6;
    auto interpolationResults = averaging.GetResults();
    ASSERT_NEAR(-2.9551821000000000, interpolationResults[0], tolerance);
    ASSERT_NEAR(-3.9075630000000001, interpolationResults[1], tolerance);
    ASSERT_NEAR(-4.3837535000000001, interpolationResults[2], tolerance);
    ASSERT_NEAR(-3.7394957999999998, interpolationResults[3], tolerance);
    ASSERT_NEAR(-4.6918768000000002, interpolationResults[4], tolerance);
    ASSERT_NEAR(-5.1680672000000003, interpolationResults[5], tolerance);
    ASSERT_NEAR(-4.5238094999999996, interpolationResults[6], tolerance);
    ASSERT_NEAR(-5.4761905000000004, interpolationResults[7], tolerance);
    ASSERT_NEAR(-5.9523809999999999, interpolationResults[8], tolerance);
    ASSERT_NEAR(-5.3081231999999998, interpolationResults[9], tolerance);
    ASSERT_NEAR(-6.2605041999999997, interpolationResults[10], tolerance);
    ASSERT_NEAR(-6.7366947000000001, interpolationResults[11], tolerance);
    ASSERT_NEAR(-6.0924370000000003, interpolationResults[12], tolerance);
    ASSERT_NEAR(-7.0448179000000000, interpolationResults[13], tolerance);
    ASSERT_NEAR(-7.5210084000000004, interpolationResults[14], tolerance);
    ASSERT_NEAR(-6.4845937999999999, interpolationResults[15], tolerance);
    ASSERT_NEAR(-7.4369747999999998, interpolationResults[16], tolerance);
    ASSERT_NEAR(-7.9131653000000002, interpolationResults[17], tolerance);
}

TEST(Averaging, InterpolateOnFacesMin)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile("../../../../tests/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/AveragingInterpolationTests/sample_grid_coarse_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), 0);

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(mesh, samples, meshkernel::AveragingInterpolation::Method::Min, meshkernel::InterpolationLocation::Faces, 1.01, false, false);
    averaging.Compute();

    constexpr double tolerance = 1e-6;
    auto interpolationResults = averaging.GetResults();
    ASSERT_NEAR(-3.8235294000000000, interpolationResults[0], tolerance);
    ASSERT_NEAR(-4.7759103999999999, interpolationResults[1], tolerance);
    ASSERT_NEAR(-4.6078431000000002, interpolationResults[2], tolerance);
    ASSERT_NEAR(-5.5602241000000001, interpolationResults[3], tolerance);
    ASSERT_NEAR(-5.3921568999999998, interpolationResults[4], tolerance);
    ASSERT_NEAR(-6.3445378000000003, interpolationResults[5], tolerance);
    ASSERT_NEAR(-6.1764706000000000, interpolationResults[6], tolerance);
    ASSERT_NEAR(-7.1288514999999997, interpolationResults[7], tolerance);
    ASSERT_NEAR(-6.9607843000000003, interpolationResults[8], tolerance);
    ASSERT_NEAR(-7.9131653000000002, interpolationResults[9], tolerance);
}

// InverseWeightedDistance
TEST(Averaging, InterpolateOnEdgesInverseWeightedDistance)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile("../../../../tests/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/AveragingInterpolationTests/sample_grid_coarse_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), 0);

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(mesh, samples, meshkernel::AveragingInterpolation::Method::InverseWeightedDistance, meshkernel::InterpolationLocation::Edges, 1.01, false, false);
    averaging.Compute();

    constexpr double tolerance = 1e-6;
    auto interpolationResults = averaging.GetResults();
    ASSERT_NEAR(-2.48030306628047, interpolationResults[0], tolerance);
    ASSERT_NEAR(-3.43184384480726, interpolationResults[1], tolerance);
    ASSERT_NEAR(-4.38310967989827, interpolationResults[2], tolerance);
    ASSERT_NEAR(-3.26444992295321, interpolationResults[3], tolerance);
    ASSERT_NEAR(-4.21568625020467, interpolationResults[4], tolerance);
    ASSERT_NEAR(-5.16692267681247, interpolationResults[5], tolerance);
    ASSERT_NEAR(-4.04876362305277, interpolationResults[6], tolerance);
    ASSERT_NEAR(-5.00000000000000, interpolationResults[7], tolerance);
    ASSERT_NEAR(-5.95123637694723, interpolationResults[8], tolerance);
    ASSERT_NEAR(-4.83307732318753, interpolationResults[9], tolerance);
    ASSERT_NEAR(-5.78431374979533, interpolationResults[10], tolerance);
    ASSERT_NEAR(-6.73555007704679, interpolationResults[11], tolerance);
    ASSERT_NEAR(-5.61689032010173, interpolationResults[12], tolerance);
    ASSERT_NEAR(-6.56815615519274, interpolationResults[13], tolerance);
    ASSERT_NEAR(-7.51969693371953, interpolationResults[14], tolerance);
    ASSERT_NEAR(-2.56423567460090, interpolationResults[15], tolerance);
    ASSERT_NEAR(-3.51580591130667, interpolationResults[16], tolerance);
    ASSERT_NEAR(-3.34791123648683, interpolationResults[17], tolerance);
    ASSERT_NEAR(-4.29914761339886, interpolationResults[18], tolerance);
    ASSERT_NEAR(-4.13222493667105, interpolationResults[19], tolerance);
}

TEST(Averaging, InterpolateOnNodesInverseWeightedDistance)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile("../../../../tests/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/AveragingInterpolationTests/sample_grid_coarse_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), 0);

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(mesh, samples, meshkernel::AveragingInterpolation::Method::InverseWeightedDistance, meshkernel::InterpolationLocation::Nodes, 1.01, false, false);
    averaging.Compute();

    constexpr double tolerance = 1e-6;
    auto interpolationResults = averaging.GetResults();
    ASSERT_NEAR(-2.0883130596575117, interpolationResults[0], tolerance);
    ASSERT_NEAR(-3.0401582895442973, interpolationResults[1], tolerance);
    ASSERT_NEAR(-3.9914535330690506, interpolationResults[2], tolerance);
    ASSERT_NEAR(-2.8722930729034286, interpolationResults[3], tolerance);
    ASSERT_NEAR(-3.8235294000702318, interpolationResults[4], tolerance);
    ASSERT_NEAR(-4.7747658267274895, interpolationResults[5], tolerance);
    ASSERT_NEAR(-3.6566067730029892, interpolationResults[6], tolerance);
    ASSERT_NEAR(-4.6078431003391067, interpolationResults[7], tolerance);
    ASSERT_NEAR(-5.5590795268974498, interpolationResults[8], tolerance);
    ASSERT_NEAR(-4.4409204731025502, interpolationResults[9], tolerance);
    ASSERT_NEAR(-5.3921568996608942, interpolationResults[10], tolerance);
    ASSERT_NEAR(-6.3433932269970112, interpolationResults[11], tolerance);
    ASSERT_NEAR(-5.2252341732725114, interpolationResults[12], tolerance);
    ASSERT_NEAR(-6.1764705999297691, interpolationResults[13], tolerance);
    ASSERT_NEAR(-7.1277069270965718, interpolationResults[14], tolerance);
    ASSERT_NEAR(-6.0085464669309490, interpolationResults[15], tolerance);
    ASSERT_NEAR(-6.9598417104557040, interpolationResults[16], tolerance);
    ASSERT_NEAR(-7.9116869403424879, interpolationResults[17], tolerance);
}

TEST(Averaging, InterpolateOnFacesInverseWeightedDistance)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile("../../../../tests/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/AveragingInterpolationTests/sample_grid_coarse_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), 0);

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(mesh, samples, meshkernel::AveragingInterpolation::Method::InverseWeightedDistance, meshkernel::InterpolationLocation::Faces, 1.01, false, false);
    averaging.Compute();

    constexpr double tolerance = 1e-6;
    auto interpolationResults = averaging.GetResults();
    ASSERT_NEAR(-2.9551820997311253, interpolationResults[0], tolerance);
    ASSERT_NEAR(-3.9075630001695529, interpolationResults[1], tolerance);
    ASSERT_NEAR(-3.7394957999999998, interpolationResults[2], tolerance);
    ASSERT_NEAR(-4.6918767996608945, interpolationResults[3], tolerance);
    ASSERT_NEAR(-4.5238095001695529, interpolationResults[4], tolerance);
    ASSERT_NEAR(-5.4761904998304463, interpolationResults[5], tolerance);
    ASSERT_NEAR(-5.3081232003391055, interpolationResults[6], tolerance);
    ASSERT_NEAR(-6.2605042000000006, interpolationResults[7], tolerance);
    ASSERT_NEAR(-6.0924369998304471, interpolationResults[8], tolerance);
    ASSERT_NEAR(-7.0448179002688773, interpolationResults[9], tolerance);
}

// MinAbsValue
TEST(Averaging, InterpolateOnEdgesMinAbsValue)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile("../../../../tests/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/AveragingInterpolationTests/sample_grid_coarse_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), 0);

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(mesh, samples, meshkernel::AveragingInterpolation::Method::MinAbsValue, meshkernel::InterpolationLocation::Edges, 1.01, false, false);
    averaging.Compute();

    constexpr double tolerance = 1e-6;
    auto interpolationResults = averaging.GetResults();
    ASSERT_NEAR(2.28291315000000, interpolationResults[0], tolerance);
    ASSERT_NEAR(2.75910365000000, interpolationResults[1], tolerance);
    ASSERT_NEAR(3.71148460000000, interpolationResults[2], tolerance);
    ASSERT_NEAR(2.87114845000000, interpolationResults[3], tolerance);
    ASSERT_NEAR(3.34733895000000, interpolationResults[4], tolerance);
    ASSERT_NEAR(4.29971990000000, interpolationResults[5], tolerance);
    ASSERT_NEAR(3.65546215000000, interpolationResults[6], tolerance);
    ASSERT_NEAR(4.13165265000000, interpolationResults[7], tolerance);
    ASSERT_NEAR(5.08403365000000, interpolationResults[8], tolerance);
    ASSERT_NEAR(4.43977590000000, interpolationResults[9], tolerance);
    ASSERT_NEAR(4.91596635000000, interpolationResults[10], tolerance);
    ASSERT_NEAR(5.86834735000000, interpolationResults[11], tolerance);
    ASSERT_NEAR(5.22408965000000, interpolationResults[12], tolerance);
    ASSERT_NEAR(5.70028010000000, interpolationResults[13], tolerance);
    ASSERT_NEAR(6.65266105000000, interpolationResults[14], tolerance);
    ASSERT_NEAR(2.32492995000000, interpolationResults[15], tolerance);
    ASSERT_NEAR(3.03921570000000, interpolationResults[16], tolerance);
    ASSERT_NEAR(2.71708685000000, interpolationResults[17], tolerance);
    ASSERT_NEAR(3.43137255000000, interpolationResults[18], tolerance);
    ASSERT_NEAR(3.50140055000000, interpolationResults[19], tolerance);
}

TEST(Averaging, InterpolateOnNodesMinAbsValue)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile("../../../../tests/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/AveragingInterpolationTests/sample_grid_coarse_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), 0);

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(mesh, samples, meshkernel::AveragingInterpolation::Method::MinAbsValue, meshkernel::InterpolationLocation::Nodes, 1.01, false, false);
    averaging.Compute();

    constexpr double tolerance = 1e-6;
    auto interpolationResults = averaging.GetResults();
    ASSERT_NEAR(2.0868346999999998, interpolationResults[0], tolerance);
    ASSERT_NEAR(2.5630251999999998, interpolationResults[1], tolerance);
    ASSERT_NEAR(3.5154062000000001, interpolationResults[2], tolerance);
    ASSERT_NEAR(2.4789916000000001, interpolationResults[3], tolerance);
    ASSERT_NEAR(2.9551821000000000, interpolationResults[4], tolerance);
    ASSERT_NEAR(3.9075630000000001, interpolationResults[5], tolerance);
    ASSERT_NEAR(3.2633052999999999, interpolationResults[6], tolerance);
    ASSERT_NEAR(3.7394957999999998, interpolationResults[7], tolerance);
    ASSERT_NEAR(4.6918768000000002, interpolationResults[8], tolerance);
    ASSERT_NEAR(4.0476190000000001, interpolationResults[9], tolerance);
    ASSERT_NEAR(4.5238094999999996, interpolationResults[10], tolerance);
    ASSERT_NEAR(5.4761905000000004, interpolationResults[11], tolerance);
    ASSERT_NEAR(4.8319327999999997, interpolationResults[12], tolerance);
    ASSERT_NEAR(5.3081231999999998, interpolationResults[13], tolerance);
    ASSERT_NEAR(6.2605041999999997, interpolationResults[14], tolerance);
    ASSERT_NEAR(5.6162464999999999, interpolationResults[15], tolerance);
    ASSERT_NEAR(6.0924370000000003, interpolationResults[16], tolerance);
    ASSERT_NEAR(7.0448179000000000, interpolationResults[17], tolerance);
}

TEST(Averaging, InterpolateOnFacesMinAbsValue)
{
    std::vector<meshkernel::Sample> samples = ReadSampleFile("../../../../tests/data/AveragingInterpolationTests/inTestAveragingInterpolation.xyz");
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/AveragingInterpolationTests/sample_grid_coarse_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), 0);

    // Execute averaging
    meshkernel::AveragingInterpolation averaging(mesh, samples, meshkernel::AveragingInterpolation::Method::MinAbsValue, meshkernel::InterpolationLocation::Faces, 1.01, false, false);
    averaging.Compute();

    constexpr double tolerance = 1e-6;
    auto interpolationResults = averaging.GetResults();
    ASSERT_NEAR(2.0868346999999998, interpolationResults[0], tolerance);
    ASSERT_NEAR(3.0392157000000002, interpolationResults[1], tolerance);
    ASSERT_NEAR(2.8711484999999999, interpolationResults[2], tolerance);
    ASSERT_NEAR(3.8235294000000000, interpolationResults[3], tolerance);
    ASSERT_NEAR(3.6554622000000001, interpolationResults[4], tolerance);
    ASSERT_NEAR(4.6078431000000002, interpolationResults[5], tolerance);
    ASSERT_NEAR(4.4397758999999999, interpolationResults[6], tolerance);
    ASSERT_NEAR(5.3921568999999998, interpolationResults[7], tolerance);
    ASSERT_NEAR(5.2240896000000001, interpolationResults[8], tolerance);
    ASSERT_NEAR(6.1764706000000000, interpolationResults[9], tolerance);
}
