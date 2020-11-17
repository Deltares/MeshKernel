#include <MeshKernel/Mesh.hpp>
#include "Meshkernel/Entities.hpp"
#include <MeshKernel/TriangulationInterpolation.hpp>
#include <TestUtils/SampleFileReader.hpp>
#include <TestUtils/MakeMeshes.hpp>
#include <gtest/gtest.h>

TEST(TriangleInterpolation, InterpolateOnNodes)
{
    // Set up
    std::vector<meshkernel::Sample> samples = ReadSampleFile("../../../../tests/data/TriangleInterpolationTests/inTestTriangleInterpolation.xyz");
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/TriangleInterpolationTests/simple_grid_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), 0);

    meshkernel::TriangulationInterpolation triangulationInterpolation(mesh->m_nodes, samples, meshkernel::Projections::cartesian);
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
    std::vector<meshkernel::Sample> samples = ReadSampleFile("../../../../tests/data/TriangleInterpolationTests/inTestTriangleInterpolation.xyz");
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/TriangleInterpolationTests/simple_grid_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), 0);

    mesh->ComputeEdgesCenters();

    meshkernel::TriangulationInterpolation triangulationInterpolation(mesh->m_edgesCenters, samples, meshkernel::Projections::cartesian);
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
    std::vector<meshkernel::Sample> samples = ReadSampleFile("../../../../tests/data/TriangleInterpolationTests/inTestTriangleInterpolation.xyz");
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/TriangleInterpolationTests/simple_grid_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), 0);

    meshkernel::TriangulationInterpolation triangulationInterpolation(mesh->m_facesMassCenters, samples, meshkernel::Projections::cartesian);
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
    std::vector<meshkernel::Sample> samples = ReadSampleFile("../../../../tests/data/TriangleInterpolationTests/SphericalCutted.xyz");
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/TriangleInterpolationTests/SphericalCutted.nc", meshkernel::Projections::cartesian);
    ASSERT_GT(mesh->GetNumNodes(), 0);
    ASSERT_GT(samples.size(), 0);

    meshkernel::TriangulationInterpolation triangulationInterpolation(mesh->m_facesMassCenters, samples, meshkernel::Projections::sphericalAccurate);
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
