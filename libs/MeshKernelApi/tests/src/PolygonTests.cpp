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
