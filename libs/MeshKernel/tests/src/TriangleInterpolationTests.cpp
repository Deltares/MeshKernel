#include <gtest/gtest.h>

#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/TriangulationInterpolation.hpp>
#include <TestUtils/Definitions.hpp>
#include <TestUtils/MakeMeshes.hpp>
#include <TestUtils/SampleFileReader.hpp>

TEST(TriangleInterpolation, TriangleInterpolation_OnNodesWithSphericalCoordinates_Shouldinterpolate)
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
    meshkernel::TriangulationInterpolation triangulationInterpolation(mesh->m_nodes, samples, meshkernel::Projection::spherical);
    triangulationInterpolation.Compute();

    // Assert
    auto interpolationResults = triangulationInterpolation.GetResults();
    constexpr double tolerance = 1e-6;
    ASSERT_NEAR(-999.00000000000000, interpolationResults[0], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[1], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[2], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[3], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[4], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[5], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[6], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[7], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[8], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[9], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[10], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[11], tolerance);
    ASSERT_NEAR(2.0000809724570234, interpolationResults[12], tolerance);
    ASSERT_NEAR(1.5001095301135865, interpolationResults[13], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[14], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[15], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[16], tolerance);
    ASSERT_NEAR(1.4999190275429768, interpolationResults[17], tolerance);
    ASSERT_NEAR(2.0000000000000000, interpolationResults[18], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[19], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[20], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[21], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[22], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[23], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[24], tolerance);
}
TEST(TriangleInterpolation, InterpolateOnNodes)
{
    // Set up
    std::vector<meshkernel::Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/TriangleInterpolationTests/inTestTriangleInterpolation.xyz");
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/TriangleInterpolationTests/simple_grid_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), 0);

    meshkernel::TriangulationInterpolation triangulationInterpolation(mesh->m_nodes, samples, meshkernel::Projection::cartesian);
    triangulationInterpolation.Compute();

    const auto results = triangulationInterpolation.GetResults();

    // test internal results
    constexpr double tolerance = 1e-9;
    ASSERT_NEAR(2.0, results[20], tolerance);
    ASSERT_NEAR(2.0, results[21], tolerance);
    ASSERT_NEAR(2.0, results[22], tolerance);
    ASSERT_NEAR(2.0, results[23], tolerance);
    ASSERT_NEAR(2.0, results[24], tolerance);
    ASSERT_NEAR(1.0, results[25], tolerance);
    ASSERT_NEAR(1.0, results[26], tolerance);
    ASSERT_NEAR(1.0, results[27], tolerance);
    ASSERT_NEAR(1.0, results[28], tolerance);
    ASSERT_NEAR(1.0, results[29], tolerance);
}

TEST(TriangleInterpolation, InterpolateOnEdges)
{
    // Set up
    std::vector<meshkernel::Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/TriangleInterpolationTests/inTestTriangleInterpolation.xyz");
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/TriangleInterpolationTests/simple_grid_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), 0);

    mesh->ComputeEdgesCenters();

    meshkernel::TriangulationInterpolation triangulationInterpolation(mesh->m_edgesCenters, samples, meshkernel::Projection::cartesian);
    triangulationInterpolation.Compute();

    const auto results = triangulationInterpolation.GetResults();

    // test internal results
    constexpr double tolerance = 1e-9;
    ASSERT_NEAR(2.0, results[14], tolerance);
    ASSERT_NEAR(2.0, results[15], tolerance);
    ASSERT_NEAR(2.0, results[16], tolerance);
    ASSERT_NEAR(2.0, results[17], tolerance);
    ASSERT_NEAR(2.0, results[18], tolerance);
    ASSERT_NEAR(2.0, results[19], tolerance);
    ASSERT_NEAR(2.0, results[20], tolerance);
    ASSERT_NEAR(2.0, results[21], tolerance);
    ASSERT_NEAR(2.0, results[22], tolerance);
    ASSERT_NEAR(2.0, results[23], tolerance);
}
TEST(TriangleInterpolation, InterpolateOnFaces)
{
    // Set up
    std::vector<meshkernel::Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/TriangleInterpolationTests/inTestTriangleInterpolation.xyz");
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/TriangleInterpolationTests/simple_grid_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), 0);

    meshkernel::TriangulationInterpolation triangulationInterpolation(mesh->m_facesMassCenters, samples, meshkernel::Projection::cartesian);
    triangulationInterpolation.Compute();

    const auto results = triangulationInterpolation.GetResults();

    // test internal results
    constexpr double tolerance = 1e-9;
    ASSERT_NEAR(2.0, results[15], tolerance);
    ASSERT_NEAR(2.0, results[16], tolerance);
    ASSERT_NEAR(2.0, results[17], tolerance);
    ASSERT_NEAR(2.0, results[18], tolerance);
    ASSERT_NEAR(2.0, results[19], tolerance);
    ASSERT_NEAR(0.0, results[20], tolerance);
    ASSERT_NEAR(0.0, results[21], tolerance);
    ASSERT_NEAR(0.0, results[22], tolerance);
    ASSERT_NEAR(0.0, results[23], tolerance);
    ASSERT_NEAR(0.0, results[24], tolerance);
    ASSERT_NEAR(0.0, results[25], tolerance);
}

TEST(TriangleInterpolation, InterpolateOnFacesUsingSphericalAccurateOption)
{
    // Set up
    std::vector<meshkernel::Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/TriangleInterpolationTests/SphericalCutted.xyz");
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/TriangleInterpolationTests/SphericalCutted.nc", meshkernel::Projection::cartesian);
    ASSERT_GT(mesh->GetNumNodes(), 0);
    ASSERT_GT(samples.size(), 0);

    meshkernel::TriangulationInterpolation triangulationInterpolation(mesh->m_facesMassCenters, samples, meshkernel::Projection::sphericalAccurate);
    triangulationInterpolation.Compute();

    const auto results = triangulationInterpolation.GetResults();

    // test internal results
    constexpr double tolerance = 1e-9;
    ASSERT_NEAR(-27.108281995694892, results[0], tolerance);
    ASSERT_NEAR(-26.010901908590927, results[1], tolerance);
    ASSERT_NEAR(-26.761085257511070, results[2], tolerance);
    ASSERT_NEAR(-26.491618738949011, results[3], tolerance);
    ASSERT_NEAR(-26.993955482433620, results[4], tolerance);
    ASSERT_NEAR(-26.897163462761789, results[5], tolerance);
    ASSERT_NEAR(-27.332152155397115, results[6], tolerance);
    ASSERT_NEAR(-28.377394828176936, results[7], tolerance);
    ASSERT_NEAR(-22.334281896214499, results[8], tolerance);
    ASSERT_NEAR(-30.427438741751285, results[9], tolerance);
    ASSERT_NEAR(-22.408269339466585, results[10], tolerance);
    ASSERT_NEAR(-13.373695239170697, results[11], tolerance);
    ASSERT_NEAR(-19.085797819738595, results[12], tolerance);
    ASSERT_NEAR(-33.579059226025457, results[13], tolerance);
    ASSERT_NEAR(-35.306462411394321, results[14], tolerance);
    ASSERT_NEAR(-32.537025752660952, results[15], tolerance);
    ASSERT_NEAR(-28.309418171810119, results[16], tolerance);
    ASSERT_NEAR(-26.425156938934055, results[17], tolerance);
    ASSERT_NEAR(-26.988893382269104, results[18], tolerance);
    ASSERT_NEAR(-29.549320886988440, results[19], tolerance);
}
