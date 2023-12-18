#include <chrono>
#include <gtest/gtest.h>
#include <random>

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/Mesh2DGenerateGlobalGrid.hpp"
#include "MeshKernel/Polygons.hpp"
#include "TestUtils/MakeMeshes.hpp"

TEST(GlobalGridTest, GenerateGrid)
{
    const std::vector<meshkernel::Point> polygonNodes{};

    const meshkernel::Polygons polygon(polygonNodes, meshkernel::Projection::sphericalAccurate);

    const auto mesh = meshkernel::Mesh2DGenerateGlobalGrid::Compute(192, 250, polygon);

    // Assert data
    const double tolerance = 1e-6;

    const auto numEdges = mesh->GetNumEdges();
    const auto numNodes = mesh->GetNumNodes();
    const auto numFaces = mesh->GetNumFaces();

    ASSERT_EQ(10, mesh->GetNumEdges());
    ASSERT_EQ(10, mesh->GetNumNodes());
    ASSERT_EQ(10, mesh->GetNumFaces());

    // Nodes
    ASSERT_NEAR(0.0, mesh->m_nodes[0].x, tolerance);
    ASSERT_NEAR(0.0, mesh->m_nodes[1].x, tolerance);
    ASSERT_NEAR(0.0, mesh->m_nodes[2].x, tolerance);
    ASSERT_NEAR(0.0, mesh->m_nodes[3].x, tolerance);

    ASSERT_NEAR(0.0, mesh->m_nodes[0].x, tolerance);
    ASSERT_NEAR(0.0, mesh->m_nodes[1].x, tolerance);
    ASSERT_NEAR(0.0, mesh->m_nodes[2].x, tolerance);
    ASSERT_NEAR(0.0, mesh->m_nodes[3].x, tolerance);
}
