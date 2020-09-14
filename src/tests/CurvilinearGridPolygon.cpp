#include "../Mesh.hpp"
#include "../Entities.hpp"
#include "../Polygons.hpp"
#include "../CurvilinearGrid.hpp"
#include "../CurvilinearGridFromPolygon.cpp"
#include <gtest/gtest.h>


TEST(CurvilinearGridInPolygon, ComputeGridInPolygonWithFourthSide)
{
    std::vector<MeshKernel::Point> polygonPoints;
    polygonPoints.push_back(MeshKernel::Point{ 273.502319, 478.880432 });
    polygonPoints.push_back(MeshKernel::Point{ 274.252319, 325.128906 });
    polygonPoints.push_back(MeshKernel::Point{ 275.002350, 172.127350 });
    polygonPoints.push_back(MeshKernel::Point{ 458.003479, 157.127213 });
    polygonPoints.push_back(MeshKernel::Point{ 719.005127, 157.127213 });
    polygonPoints.push_back(MeshKernel::Point{ 741.505249, 328.128937 });
    polygonPoints.push_back(MeshKernel::Point{ 710.755066, 490.880554 });
    polygonPoints.push_back(MeshKernel::Point{ 507.503784, 494.630615 });
    polygonPoints.push_back(MeshKernel::Point{ 305.002533, 493.130615 });

    auto polygon = std::make_shared<MeshKernel::Polygons>();
    polygon->Set(polygonPoints, MeshKernel::Projections::cartesian);

    MeshKernel::CurvilinearGridFromPolygon curvilinearGridFromPolygon(polygon);

    int firstNode = 0;
    int secondNode = 2;
    int thirdNode = 4;
    bool useFourthSide = true;
    MeshKernel::CurvilinearGrid curvilinearGrid;
    bool success = curvilinearGridFromPolygon.Compute(firstNode, secondNode, thirdNode, useFourthSide, curvilinearGrid);

    // check the values
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(273.50231900000000, curvilinearGrid.m_grid[0][0].x, tolerance);
    ASSERT_NEAR(305.00253300000003, curvilinearGrid.m_grid[0][1].x, tolerance);
    ASSERT_NEAR(507.50378400000000, curvilinearGrid.m_grid[0][2].x, tolerance);

    ASSERT_NEAR(478.88043199999998, curvilinearGrid.m_grid[0][0].y, tolerance);
    ASSERT_NEAR(493.13061499999998, curvilinearGrid.m_grid[0][1].y, tolerance);
    ASSERT_NEAR(494.63061499999998, curvilinearGrid.m_grid[0][2].y, tolerance);

    ASSERT_NEAR(274.25231900000000, curvilinearGrid.m_grid[1][0].x, tolerance);
    ASSERT_NEAR(410.51616175207897, curvilinearGrid.m_grid[1][1].x, tolerance);
    ASSERT_NEAR(741.50524900000005, curvilinearGrid.m_grid[1][2].x, tolerance);

    ASSERT_NEAR(325.12890599999997, curvilinearGrid.m_grid[1][0].y, tolerance);
    ASSERT_NEAR(314.33420290273324, curvilinearGrid.m_grid[1][1].y, tolerance);
    ASSERT_NEAR(328.12893700000001, curvilinearGrid.m_grid[1][2].y, tolerance);

    ASSERT_TRUE(success);
}

TEST(CurvilinearGridInPolygon, ComputeGridInPolygonWithoutFourthSide)
{
    std::vector<MeshKernel::Point> polygonPoints;
    polygonPoints.push_back(MeshKernel::Point{ 273.502319, 478.880432 });
    polygonPoints.push_back(MeshKernel::Point{ 274.252319, 325.128906 });
    polygonPoints.push_back(MeshKernel::Point{ 275.002350, 172.127350 });
    polygonPoints.push_back(MeshKernel::Point{ 458.003479, 157.127213 });
    polygonPoints.push_back(MeshKernel::Point{ 719.005127, 157.127213 });
    polygonPoints.push_back(MeshKernel::Point{ 741.505249, 328.128937 });
    polygonPoints.push_back(MeshKernel::Point{ 710.755066, 490.880554 });
    polygonPoints.push_back(MeshKernel::Point{ 507.503784, 494.630615 });
    polygonPoints.push_back(MeshKernel::Point{ 305.002533, 493.130615 });

    auto polygon = std::make_shared<MeshKernel::Polygons>();
    polygon->Set(polygonPoints, MeshKernel::Projections::cartesian);

    MeshKernel::CurvilinearGridFromPolygon curvilinearGridFromPolygon(polygon);

    int firstNode = 0;
    int secondNode = 2;
    int thirdNode = 4;
    bool useFourthSide = false;
    MeshKernel::CurvilinearGrid curvilinearGrid;
    bool success = curvilinearGridFromPolygon.Compute(firstNode, secondNode, thirdNode, useFourthSide, curvilinearGrid);

    // check the values
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(273.50231900000000, curvilinearGrid.m_grid[0][0].x, tolerance);
    ASSERT_NEAR(492.12869250000006, curvilinearGrid.m_grid[0][1].x, tolerance);
    ASSERT_NEAR(710.75506600000006, curvilinearGrid.m_grid[0][2].x, tolerance);

    ASSERT_NEAR(478.88043199999998, curvilinearGrid.m_grid[0][0].y, tolerance);
    ASSERT_NEAR(484.88049300000000, curvilinearGrid.m_grid[0][1].y, tolerance);
    ASSERT_NEAR(490.88055400000002, curvilinearGrid.m_grid[0][2].y, tolerance);

    ASSERT_NEAR(274.25231900000000, curvilinearGrid.m_grid[1][0].x, tolerance);
    ASSERT_NEAR(481.37408996173241, curvilinearGrid.m_grid[1][1].x, tolerance);
    ASSERT_NEAR(741.50524900000005, curvilinearGrid.m_grid[1][2].x, tolerance);

    ASSERT_NEAR(325.12890599999997, curvilinearGrid.m_grid[1][0].y, tolerance);
    ASSERT_NEAR(322.93773596204318, curvilinearGrid.m_grid[1][1].y, tolerance);
    ASSERT_NEAR(328.12893700000001, curvilinearGrid.m_grid[1][2].y, tolerance);

    ASSERT_TRUE(success);
}