#include <chrono>
#include <gtest/gtest.h>
#include <random>

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/Polygons.hpp"
#include "MeshKernel/ProjectionConversions.hpp"

#include "TestUtils/Definitions.hpp"
#include "TestUtils/MakeMeshes.hpp"

namespace mk = meshkernel;

TEST(MeshConversionTests, CartesianToSphericalClassTest)
{
    mk::Point origin{0.0, 0.0};
    mk::ConvertCartesianToSpherical conversion(origin);

    EXPECT_EQ(mk::Projection::cartesian, conversion.SourceProjection());
    EXPECT_EQ(mk::Projection::spherical, conversion.TargetProjection());

#if 0
    mk::Point pnt{111.0, 11.1}; // ≃ 0.001, 0.0001 degrees
    constexpr double tolerance = 1.0e-3;
    mk::Point result = conversion(pnt);
    EXPECT_TRUE(mk::IsEqual(0.001, result.x, tolerance)) << "Expected: 0.001, found: " << result.x;
    EXPECT_TRUE(mk::IsEqual(0.0001, result.y, tolerance)) << "Expected: 0.0001, found: " << result.y;
#endif
}
