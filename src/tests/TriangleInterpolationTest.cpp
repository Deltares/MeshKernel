#include "../Mesh.hpp"
#include "../Entities.hpp"
#include "../TriangulationInterpolation.cpp"
#include "SampleFileReader.cpp"
#include "MakeMeshes.cpp"
#include <gtest/gtest.h>

TEST(TriangleInterpolation, InterpolateOnNodes)
{
    // Set up
    std::vector<meshkernel::Sample> samples = ReadSampleFile("..\\..\\tests\\TriangleInterpolationTests\\inTestTriangleInterpolation.xyz");
    auto mesh = ReadLegacyMeshFromFile("..\\..\\tests\\TriangleInterpolationTests\\simple_grid_net.nc");
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
    std::vector<meshkernel::Sample> samples = ReadSampleFile("..\\..\\tests\\TriangleInterpolationTests\\inTestTriangleInterpolation.xyz");
    auto mesh = ReadLegacyMeshFromFile("..\\..\\tests\\TriangleInterpolationTests\\simple_grid_net.nc");
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
    std::vector<meshkernel::Sample> samples = ReadSampleFile("..\\..\\tests\\TriangleInterpolationTests\\inTestTriangleInterpolation.xyz");
    auto mesh = ReadLegacyMeshFromFile("..\\..\\tests\\TriangleInterpolationTests\\simple_grid_net.nc");
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
    std::vector<meshkernel::Sample> samples = ReadSampleFile("..\\..\\tests\\TriangleInterpolationTests\\SphericalCutted.xyz");
    auto mesh = ReadLegacyMeshFromFile("..\\..\\tests\\TriangleInterpolationTests\\SphericalCutted.nc", meshkernel::Projections::cartesian);
    ASSERT_GT(mesh->GetNumNodes(), 0);
    ASSERT_GT(samples.size(), 0);

    meshkernel::TriangulationInterpolation triangulationInterpolation(mesh->m_facesMassCenters, samples, meshkernel::Projections::sphericalAccurate);
    triangulationInterpolation.Compute();

    const auto results = triangulationInterpolation.GetResults();

    // test internal results
    constexpr double tolerance = 1e-9;
    ASSERT_NEAR(-32.001516910357608, results[0], tolerance);
    ASSERT_NEAR(-32.177930753637561, results[1], tolerance);
    ASSERT_NEAR(-32.387289253798443, results[2], tolerance);
    ASSERT_NEAR(-32.676521526792378, results[3], tolerance);
    ASSERT_NEAR(-32.679307576069114, results[4], tolerance);
    ASSERT_NEAR(-32.580405055756422, results[5], tolerance);
    ASSERT_NEAR(-33.000000000000000, results[6], tolerance);
    ASSERT_NEAR(-33.000000000000000, results[7], tolerance);
    ASSERT_NEAR(-33.000000000000000, results[8], tolerance);
    ASSERT_NEAR(-32.558727666694786, results[9], tolerance);
    ASSERT_NEAR(-32.912803051993762, results[10], tolerance);
    ASSERT_NEAR(-33.000000000000000, results[11], tolerance);
    ASSERT_NEAR(-33.000000000000000, results[12], tolerance);
    ASSERT_NEAR(-33.000000000000000, results[13], tolerance);
    ASSERT_NEAR(-33.000000000000000, results[14], tolerance);
    ASSERT_NEAR(-33.172614414108558, results[15], tolerance);
    ASSERT_NEAR(-33.436592493966508, results[16], tolerance);
    ASSERT_NEAR(-33.443677472308543, results[17], tolerance);
    ASSERT_NEAR(-33.486112239353702, results[18], tolerance);
    ASSERT_NEAR(-33.527992109133827, results[19], tolerance);
    ASSERT_NEAR(-33.569329139845948, results[20], tolerance);
}