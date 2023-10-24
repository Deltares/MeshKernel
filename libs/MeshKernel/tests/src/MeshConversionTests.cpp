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

TEST (MeshConversionTests, SphericalToCartesianClassTest)
{
    constexpr double tolerance = 1.0e-3;

    mk::Point origin{0.0, 0.0};
    mk::ConvertSphericalToCartesian conversion(origin);

    mk::Point pnt{0.001, 0.0001}; // ≃ 111.3 and 11.13 metres

    mk::Point result = conversion (pnt);

    EXPECT_EQ ( mk::Projection::spherical, conversion.SourceProjection ());
    EXPECT_EQ ( mk::Projection::cartesian, conversion.TargetProjection ());
    EXPECT_TRUE ( mk::IsEqual( 111.3, result.x, tolerance)) << "Expected: 111.3, found: " << result.x;
    EXPECT_TRUE ( mk::IsEqual( 11.13, result.y, tolerance)) << "Expected: 11.13, found: " << result.y;

}
