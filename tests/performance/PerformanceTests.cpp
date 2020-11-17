#include <chrono>

#include "../unit/MakeMeshes.cpp"
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/FlipEdges.hpp>
#include <MeshKernel/LandBoundaries.hpp>
#include <MeshKernel/Mesh.hpp>
#include <MeshKernel/OrthogonalizationAndSmoothing.hpp>
#include <MeshKernel/OrthogonalizationParametersNative.hpp>
#include <MeshKernel/Orthogonalizer.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/Smoother.hpp>

#include <gtest/gtest.h>

TEST(PerformanceTest, FlipEdgesAndOrthogonalization)
{
    // load half a million mesh
    auto mesh = ReadLegacyMeshFromFile("..\\..\\..\\tests\\PerformanceTests\\SquaredDomainWithTriangles_500by500_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), 0);

    // set landBoundaries
    auto polygon = std::make_shared<meshkernel::Polygons>();
    std::vector<meshkernel::Point> landBoundary;
    auto landBoundaries = std::make_shared<meshkernel::LandBoundaries>(landBoundary, mesh, polygon);

    //execute flip edges
    meshkernel::FlipEdges flipEdges(mesh, landBoundaries, true, false);

    auto successful = flipEdges.Compute();
    ASSERT_TRUE(successful);

    int projectToLandBoundaryOption = 0;
    meshkernelapi::OrthogonalizationParametersNative orthogonalizationParametersNative;
    orthogonalizationParametersNative.OuterIterations = 2;
    orthogonalizationParametersNative.InnerIterations = 25;
    orthogonalizationParametersNative.BoundaryIterations = 25;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactor = 0.975;
    orthogonalizationParametersNative.OrthogonalizationToSmoothingFactorBoundary = 0.975;
    orthogonalizationParametersNative.Smoothorarea = 0.0;

    auto orthogonalizer = std::make_shared<meshkernel::Orthogonalizer>(mesh);
    auto smoother = std::make_shared<meshkernel::Smoother>(mesh);
    auto landboundaries = std::make_shared<meshkernel::LandBoundaries>(landBoundary, mesh, polygon);
    meshkernel::OrthogonalizationAndSmoothing orthogonalization(mesh,
                                                                smoother,
                                                                orthogonalizer,
                                                                polygon,
                                                                landboundaries,
                                                                projectToLandBoundaryOption,
                                                                orthogonalizationParametersNative);
    // execute orthogonalization
    successful = orthogonalization.Initialize();
    ASSERT_TRUE(successful);
    successful = orthogonalization.Compute();
    ASSERT_TRUE(successful);

    // test positions of internal mesh nodes after flipedges and orthogonalization
    const double tolerance = 1e-3;
    ASSERT_NEAR(1060.54825, mesh->m_nodes[10020].x, tolerance);
    ASSERT_NEAR(1055.23662, mesh->m_nodes[10021].x, tolerance);
    ASSERT_NEAR(1058.32669, mesh->m_nodes[10022].x, tolerance);
    ASSERT_NEAR(1091.84353, mesh->m_nodes[10023].x, tolerance);
    ASSERT_NEAR(1089.30253, mesh->m_nodes[10024].x, tolerance);
    ASSERT_NEAR(1110.78646, mesh->m_nodes[10025].x, tolerance);
    ASSERT_NEAR(1117.97993, mesh->m_nodes[10026].x, tolerance);
    ASSERT_NEAR(1119.31475, mesh->m_nodes[10027].x, tolerance);
    ASSERT_NEAR(1120.75278, mesh->m_nodes[10028].x, tolerance);
    ASSERT_NEAR(1133.87846, mesh->m_nodes[10029].x, tolerance);

    ASSERT_NEAR(13.7516096, mesh->m_nodes[10020].y, tolerance);
    ASSERT_NEAR(23.4022358, mesh->m_nodes[10021].y, tolerance);
    ASSERT_NEAR(6.67204648, mesh->m_nodes[10022].y, tolerance);
    ASSERT_NEAR(6.49722192, mesh->m_nodes[10023].y, tolerance);
    ASSERT_NEAR(13.1681871, mesh->m_nodes[10024].y, tolerance);
    ASSERT_NEAR(13.2618403, mesh->m_nodes[10025].y, tolerance);
    ASSERT_NEAR(16.5404488, mesh->m_nodes[10026].y, tolerance);
    ASSERT_NEAR(30.9726827, mesh->m_nodes[10027].y, tolerance);
    ASSERT_NEAR(23.6513297, mesh->m_nodes[10028].y, tolerance);
    ASSERT_NEAR(29.8694322, mesh->m_nodes[10029].y, tolerance);
}
