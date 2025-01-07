#include "MeshKernelApi/MeshKernel.hpp"
#include <gtest/gtest.h>
#include <random>

// namespace aliases
namespace mk = meshkernel;
namespace mkapi = meshkernelapi;

TEST(PolygonTests, PolygonRefinementTests)
{
    int meshKernelId;
    auto errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::GeometryList land;
    land.geometry_separator = meshkernel::constants::missing::doubleValue;
    std::vector xLand{0.0, 10.0};
    std::vector yLand{10.0, 10.0};
    land.coordinates_x = xLand.data();
    land.coordinates_y = yLand.data();
    land.num_coordinates = static_cast<int>(xLand.size());

    meshkernelapi::GeometryList polygon;
    polygon.geometry_separator = meshkernel::constants::missing::doubleValue;
    std::vector xPolygon{1.0, 1.0, 7.0, 7.0, 1.0};
    std::vector yPolygon{9.0, 0.0, 0.0, 9.0, 9.0};
    polygon.coordinates_x = xPolygon.data();
    polygon.coordinates_y = yPolygon.data();
    polygon.num_coordinates = static_cast<int>(xPolygon.size());

    errorCode = mkernel_polygon_snap_to_landboundary(meshKernelId, land, polygon, 3, 0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    constexpr double tolerance = 1e-6;
    ASSERT_NEAR(10.0, yPolygon[0], tolerance);
    ASSERT_NEAR(0.0, yPolygon[1], tolerance);
    ASSERT_NEAR(0.0, yPolygon[2], tolerance);
    ASSERT_NEAR(10.0, yPolygon[3], tolerance);
    ASSERT_NEAR(10.0, yPolygon[4], tolerance);
}

TEST(PolygonTests, PolygonWithHoleSnappingTest)
{
    constexpr double tolerance = 1.0e-8;

    int meshKernelId;
    auto errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::GeometryList land;
    land.geometry_separator = meshkernel::constants::missing::doubleValue;
    std::vector xLand{-1.0, 31.0};
    std::vector yLand{-0.25, -0.25};
    land.coordinates_x = xLand.data();
    land.coordinates_y = yLand.data();
    land.num_coordinates = static_cast<int>(xLand.size());

    meshkernelapi::GeometryList polygon;
    polygon.geometry_separator = meshkernel::constants::missing::doubleValue;
    polygon.inner_outer_separator = meshkernel::constants::missing::innerOuterSeparator;

    std::vector xPolygon{0.0, 10.0, 10.0, 0.0, 0.0, meshkernel::constants::missing::innerOuterSeparator, 2.0, 7.0, 7.0, 2.0, 2.0, meshkernel::constants::missing::doubleValue, 20.0, 30.0, 30.0, 20.0, 20.0};
    std::vector yPolygon{0.0, 0.0, 10.0, 10.0, 0.0, meshkernel::constants::missing::innerOuterSeparator, 2.0, 2.0, 7.0, 7.0, 2.0, meshkernel::constants::missing::doubleValue, 0.0, 0.0, 10.0, 10.0, 0.0};

    polygon.coordinates_x = xPolygon.data();
    polygon.coordinates_y = yPolygon.data();
    polygon.num_coordinates = static_cast<int>(xPolygon.size());

    errorCode = mkernel_polygon_snap_to_landboundary(meshKernelId, land, polygon, 4, 1);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    const std::vector<double> expectedXAfterFirstSnapping{0.0, 10.0, 10.0, 0.0, 0.0, -998.0, 2.0, 7.0, 7.0, 2.0, 2.0, -999.0, 20.0, 30.0, 30.0, 20.0, 20.0};
    const std::vector<double> expectedYAfterFirstSnapping{-0.25, -0.25, 10.0, 10.0, -0.25, -998.0, 2.0, 2.0, 7.0, 7.0, 2.0, -999.0, 0.0, 0.0, 10.0, 10.0, 0.0};

    ASSERT_EQ(polygon.num_coordinates, static_cast<int>(expectedXAfterFirstSnapping.size()));

    for (size_t i = 0; i < expectedXAfterFirstSnapping.size(); ++i)
    {
        EXPECT_NEAR(polygon.coordinates_x[i], expectedXAfterFirstSnapping[i], tolerance);
        EXPECT_NEAR(polygon.coordinates_y[i], expectedYAfterFirstSnapping[i], tolerance);
    }

    errorCode = mkernel_polygon_snap_to_landboundary(meshKernelId, land, polygon, 16, 13);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    const std::vector<double> expectedXAfterSecondSnapping{0.0, 10.0, 10.0, 0.0, 0.0, -998.0, 2.0, 7.0, 7.0, 2.0, 2.0, -999.0, 20.0, 30.0, 30.0, 20.0, 20.0};
    const std::vector<double> expectedYAfterSecondSnapping{-0.25, -0.25, 10.0, 10.0, -0.25, -998.0, 2.0, 2.0, 7.0, 7.0, 2.0, -999.0, -0.25, -0.25, 10.0, 10.0, -0.25};

    for (size_t i = 0; i < expectedXAfterSecondSnapping.size(); ++i)
    {
        EXPECT_NEAR(polygon.coordinates_x[i], expectedXAfterSecondSnapping[i], tolerance);
        EXPECT_NEAR(polygon.coordinates_y[i], expectedYAfterSecondSnapping[i], tolerance);
    }
}
