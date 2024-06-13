#include <chrono>
#include <gtest/gtest.h>
#include <random>

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/Mesh2DGenerateGlobal.hpp"
#include "MeshKernel/Polygons.hpp"
#include "TestUtils/MakeMeshes.hpp"

TEST(GlobalGridTest, Mesh2DGenerateGlobalCompute_ShouldGenerateMesh)
{
    // Execute
    constexpr meshkernel::UInt numLongitudeNodes = 19;
    constexpr meshkernel::UInt numLatitudeNodes = 25;
    const auto mesh = meshkernel::Mesh2DGenerateGlobal::Compute(numLongitudeNodes, numLatitudeNodes, meshkernel::Projection::spherical);

    // Assert
    ASSERT_EQ(1233, mesh->GetNumValidEdges());
    ASSERT_EQ(629, mesh->GetNumValidNodes());

    const double tolerance = 1e-6;

    ASSERT_NEAR(-161.05263157894737, mesh->Node(0).x, tolerance);
    ASSERT_NEAR(-161.05263157894737, mesh->Node(1).x, tolerance);
    ASSERT_NEAR(-161.05263157894737, mesh->Node(2).x, tolerance);
    ASSERT_NEAR(-142.10526315789474, mesh->Node(3).x, tolerance);

    ASSERT_NEAR(0.0000000000000000, mesh->Node(0).y, tolerance);
    ASSERT_NEAR(18.695753703140564, mesh->Node(1).y, tolerance);
    ASSERT_NEAR(-18.695753703140564, mesh->Node(2).y, tolerance);
    ASSERT_NEAR(0.0000000000000000, mesh->Node(3).y, tolerance);
}

TEST(GlobalGridTest, Mesh2DGenerateGlobalCompute_OnInvalidProjection_ShouldThrowAnException)
{
    // Assert
    const meshkernel::UInt numLongitudeNodes = 19;
    const meshkernel::UInt numLatitudeNodes = 25;
    EXPECT_THROW(meshkernel::Mesh2DGenerateGlobal::Compute(numLongitudeNodes, numLatitudeNodes, meshkernel::Projection::cartesian), meshkernel::MeshKernelError);
}
