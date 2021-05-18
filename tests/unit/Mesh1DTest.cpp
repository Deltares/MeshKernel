#include <gtest/gtest.h>

#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh1D.hpp>
#include <TestUtils/MakeMeshes.hpp>

TEST(Mesh1D, GenerateMeshFromPolyLines_WithOverlappingNodes_ShouldRemoveOverlappingNodes)
{
    //1 Setup
    std::vector<std::vector<meshkernel::Point>> polylines{
        {{0.0, 0.0},
         {10.0, 0.0},
         {20.0, 0.0}},
        {{10.0, -10.0},
         {10.0, 0.0},
         {10.0, 10.0}}};

    // 2 Execution
    const auto mesh = meshkernel::Mesh1D(polylines, 5.0, 0.01, meshkernel::Projection::cartesian);

    // 3 Assertion
    const auto tolerance = 1e-6;
    ASSERT_NEAR(0.00, mesh.m_nodes[0].x, tolerance);
    ASSERT_NEAR(5.00, mesh.m_nodes[1].x, tolerance);
    ASSERT_NEAR(15.0, mesh.m_nodes[2].x, tolerance);
    ASSERT_NEAR(20.0, mesh.m_nodes[3].x, tolerance);
    ASSERT_NEAR(10.0, mesh.m_nodes[4].x, tolerance);
    ASSERT_NEAR(10.0, mesh.m_nodes[5].x, tolerance);
    ASSERT_NEAR(10.0, mesh.m_nodes[6].x, tolerance);
    ASSERT_NEAR(10.0, mesh.m_nodes[7].x, tolerance);
    ASSERT_NEAR(10.0, mesh.m_nodes[8].x, tolerance);

    ASSERT_NEAR(0.0, mesh.m_nodes[0].y, tolerance);
    ASSERT_NEAR(0.0, mesh.m_nodes[1].y, tolerance);
    ASSERT_NEAR(0.0, mesh.m_nodes[2].y, tolerance);
    ASSERT_NEAR(0.0, mesh.m_nodes[3].y, tolerance);
    ASSERT_NEAR(-10.0, mesh.m_nodes[4].y, tolerance);
    ASSERT_NEAR(-5.0, mesh.m_nodes[5].y, tolerance);
    ASSERT_NEAR(0.0, mesh.m_nodes[6].y, tolerance);
    ASSERT_NEAR(5.0, mesh.m_nodes[7].y, tolerance);
    ASSERT_NEAR(10.0, mesh.m_nodes[8].y, tolerance);
}