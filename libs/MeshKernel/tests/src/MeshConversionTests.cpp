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

TEST(MeshConversionTests, CartesianToSphericalClassTestFromCode)
{
    constexpr double tolerance = 1.0e-6;

    constexpr int EpsgCode = 3043;

    mk::ConvertCartesianToSphericalEPSG<EpsgCode> conversion;

    mk::Point source{624402.8057, 5808292.0905};
    // coordinate are in (longitude, latitude)
    mk::Point target{4.82898, 52.410695};
    mk::Point result = conversion(source);

    EXPECT_EQ(mk::Projection::cartesian, conversion.SourceProjection());
    EXPECT_EQ(mk::Projection::spherical, conversion.TargetProjection());
    EXPECT_TRUE(mk::IsEqual(target.x, result.x, tolerance));
    EXPECT_TRUE(mk::IsEqual(target.y, result.y, tolerance));

    source = mk::Point(430887.56433058, 4581837.85323929693549871);
    target = mk::Point(2.1734, 41.3851);
    result = conversion(source);

    EXPECT_TRUE(mk::IsEqual(target.x, result.x, tolerance));
    EXPECT_TRUE(mk::IsEqual(target.y, result.y, tolerance));

    // Deltares Delft office
    source = mk::Point(594919.016912, 5760424.72635);
    target = mk::Point(4.382253, 51.986373);
    result = conversion(source);
    EXPECT_TRUE(mk::IsEqual(target.x, result.x, tolerance));
    EXPECT_TRUE(mk::IsEqual(target.y, result.y, tolerance));

    // Brussels
    source = mk::Point(595158.7545, 5633632.17474);
    target = mk::Point(4.351697, 50.846557);
    result = conversion(source);
    EXPECT_TRUE(mk::IsEqual(target.x, result.x, tolerance));
    EXPECT_TRUE(mk::IsEqual(target.y, result.y, tolerance));

    // Paris
    source = mk::Point(452484.160, 5411718.719);
    target = mk::Point(2.352222, 48.856614);
    result = conversion(source);
    EXPECT_TRUE(mk::IsEqual(target.x, result.x, tolerance));
    EXPECT_TRUE(mk::IsEqual(target.y, result.y, tolerance));
}

TEST(MeshConversionTests, CartesianToSphericalClassTestFromStringZone31)
{
    constexpr double tolerance = 1.0e-6;

    mk::ConvertCartesianToSpherical conversion("+proj=utm +lat_1=0.5 +lat_2=2 +n=0.5 +zone=31");

    mk::Point source{624402.8057, 5808292.0905};
    mk::Point target{4.82898, 52.410695};
    mk::Point result = conversion(source);

    EXPECT_EQ(mk::Projection::cartesian, conversion.SourceProjection());
    EXPECT_EQ(mk::Projection::spherical, conversion.TargetProjection());
    EXPECT_TRUE(mk::IsEqual(target.x, result.x, tolerance));
    EXPECT_TRUE(mk::IsEqual(target.y, result.y, tolerance));

    source = mk::Point(430887.56433058, 4581837.85323929693549871);
    target = mk::Point(2.1734, 41.3851);
    result = conversion(source);
    EXPECT_TRUE(mk::IsEqual(target.x, result.x, tolerance));
    EXPECT_TRUE(mk::IsEqual(target.y, result.y, tolerance));

    // Deltares Delft office
    source = mk::Point(594919.016912, 5760424.72635);
    target = mk::Point(4.382253, 51.986373);
    result = conversion(source);
    EXPECT_TRUE(mk::IsEqual(target.x, result.x, tolerance));
    EXPECT_TRUE(mk::IsEqual(target.y, result.y, tolerance));

    // Brussels
    source = mk::Point(595158.7545, 5633632.17474);
    target = mk::Point(4.351697, 50.846557);
    result = conversion(source);
    EXPECT_TRUE(mk::IsEqual(target.x, result.x, tolerance));
    EXPECT_TRUE(mk::IsEqual(target.y, result.y, tolerance));

    // Paris
    source = mk::Point(452484.160, 5411718.719);
    target = mk::Point(2.352222, 48.856614);
    result = conversion(source);
    EXPECT_TRUE(mk::IsEqual(target.x, result.x, tolerance));
    EXPECT_TRUE(mk::IsEqual(target.y, result.y, tolerance));
}

TEST(MeshConversionTests, CartesianToSphericalClassTestFromStringZone30)
{
    constexpr double tolerance = 1.0e-6;

    mk::ConvertCartesianToSpherical conversion("+proj=utm +lat_1=0.5 +lat_2=2 +n=0.5 +zone=30N");

    // Cwmbran
    mk::Point source{498559.553, 5722516.858};
    mk::Point target{-3.020822, 51.653644};
    mk::Point result = conversion(source);

    EXPECT_TRUE(mk::IsEqual(target.x, result.x, tolerance));
    EXPECT_TRUE(mk::IsEqual(target.y, result.y, tolerance));
}
