#include <MeshKernel/Mesh.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGridFromPolygon.hpp>
#include <gtest/gtest.h>
TEST(CurvilinearGridFromPolygon, ComputeGridInPolygonWithFourthSide)
{
    std::vector<meshkernel::Point> polygonPoints{{273.502319, 478.880432},
                                                 {274.252319, 325.128906},
                                                 {275.002350, 172.127350},
                                                 {458.003479, 157.127213},
                                                 {719.005127, 157.127213},
                                                 {741.505249, 328.128937},
                                                 {710.755066, 490.880554},
                                                 {507.503784, 494.630615},
                                                 {305.002533, 493.130615}};

    auto polygon = std::make_shared<meshkernel::Polygons>(polygonPoints, meshkernel::Projections::cartesian);

    meshkernel::CurvilinearGridFromPolygon curvilinearGridFromPolygon(polygon);

    int firstNode = 0;
    int secondNode = 2;
    int thirdNode = 4;
    bool useFourthSide = true;
    meshkernel::CurvilinearGrid curvilinearGrid;
    curvilinearGridFromPolygon.Compute(firstNode, secondNode, thirdNode, useFourthSide, curvilinearGrid);

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
}

TEST(CurvilinearGridFromPolygon, ComputeGridInPolygonWithoutFourthSide)
{
    std::vector<meshkernel::Point> polygonPoints{{273.502319, 478.880432},
                                                 {274.252319, 325.128906},
                                                 {275.002350, 172.127350},
                                                 {458.003479, 157.127213},
                                                 {719.005127, 157.127213},
                                                 {741.505249, 328.128937},
                                                 {710.755066, 490.880554},
                                                 {507.503784, 494.630615},
                                                 {305.002533, 493.130615}};

    auto polygon = std::make_shared<meshkernel::Polygons>(polygonPoints, meshkernel::Projections::cartesian);

    meshkernel::CurvilinearGridFromPolygon curvilinearGridFromPolygon(polygon);

    int firstNode = 0;
    int secondNode = 2;
    int thirdNode = 4;
    bool useFourthSide = false;
    meshkernel::CurvilinearGrid curvilinearGrid;
    curvilinearGridFromPolygon.Compute(firstNode, secondNode, thirdNode, useFourthSide, curvilinearGrid);

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
}

TEST(CurvilinearGridFromPolygon, ComputeGridTriangle)
{
    std::vector<meshkernel::Point> polygonPoints{{444.504791, 437.155945},
                                                 {427.731781, 382.745758},
                                                 {405.640503, 317.699005},
                                                 {381.094666, 262.470612},
                                                 {451.050354, 262.879700},
                                                 {528.778931, 263.288788},
                                                 {593.416260, 266.561584},
                                                 {558.643005, 324.653687},
                                                 {526.733398, 377.836578},
                                                 {444.095703, 436.746857}};

    auto polygon = std::make_shared<meshkernel::Polygons>(polygonPoints, meshkernel::Projections::cartesian);

    meshkernel::CurvilinearGridFromPolygon curvilinearGridFromPolygon(polygon);

    int firstNode = 0;
    int secondNode = 3;
    int thirdNode = 6;
    meshkernel::CurvilinearGrid curvilinearGrid;
    curvilinearGridFromPolygon.Compute(firstNode, secondNode, thirdNode, curvilinearGrid);

    // check the values
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(444.50479100000001, curvilinearGrid.m_grid[0][0].x, tolerance);
    ASSERT_NEAR(444.09570300000001, curvilinearGrid.m_grid[0][1].x, tolerance);
    ASSERT_NEAR(526.73339799999997, curvilinearGrid.m_grid[0][2].x, tolerance);
    ASSERT_NEAR(558.64300500000002, curvilinearGrid.m_grid[0][3].x, tolerance);
    ASSERT_NEAR(593.41625999999997, curvilinearGrid.m_grid[0][4].x, tolerance);
    ASSERT_NEAR(-999.0000000000000, curvilinearGrid.m_grid[0][5].x, tolerance);

    ASSERT_NEAR(437.15594499999997, curvilinearGrid.m_grid[0][0].y, tolerance);
    ASSERT_NEAR(436.74685699999998, curvilinearGrid.m_grid[0][1].y, tolerance);
    ASSERT_NEAR(377.83657799999997, curvilinearGrid.m_grid[0][2].y, tolerance);
    ASSERT_NEAR(324.65368699999999, curvilinearGrid.m_grid[0][3].y, tolerance);
    ASSERT_NEAR(266.56158399999998, curvilinearGrid.m_grid[0][4].y, tolerance);
    ASSERT_NEAR(-999.0000000000000, curvilinearGrid.m_grid[0][5].y, tolerance);

    ASSERT_NEAR(427.73178100000001, curvilinearGrid.m_grid[1][0].x, tolerance);
    ASSERT_NEAR(455.85723540740742, curvilinearGrid.m_grid[1][1].x, tolerance);
    ASSERT_NEAR(483.98268981481488, curvilinearGrid.m_grid[1][2].x, tolerance);
    ASSERT_NEAR(506.38081040740741, curvilinearGrid.m_grid[1][3].x, tolerance);
    ASSERT_NEAR(528.77893099999994, curvilinearGrid.m_grid[1][4].x, tolerance);
    ASSERT_NEAR(-999.0000000000000, curvilinearGrid.m_grid[1][5].x, tolerance);

    ASSERT_NEAR(382.74575800000002, curvilinearGrid.m_grid[1][0].y, tolerance);
    ASSERT_NEAR(362.14685592592593, curvilinearGrid.m_grid[1][1].y, tolerance);
    ASSERT_NEAR(341.54795385185184, curvilinearGrid.m_grid[1][2].y, tolerance);
    ASSERT_NEAR(302.41837092592596, curvilinearGrid.m_grid[1][3].y, tolerance);
    ASSERT_NEAR(263.28878800000001, curvilinearGrid.m_grid[1][4].y, tolerance);
    ASSERT_NEAR(-999.0000000000000, curvilinearGrid.m_grid[1][5].y, tolerance);
}
