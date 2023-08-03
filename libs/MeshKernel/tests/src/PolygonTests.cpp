#include <gtest/gtest.h>

#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/LandBoundary.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/SplineAlgorithms.hpp>
#include <TestUtils/Definitions.hpp>
#include <TestUtils/MakeMeshes.hpp>

TEST(Polygons, MeshBoundaryToPolygon)
{
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/SmallTriangularGrid_net.nc");

    std::vector<meshkernel::Point> polygonNodes;
    const auto meshBoundaryPolygon = mesh->MeshBoundaryToPolygon(polygonNodes);

    ASSERT_EQ(8, meshBoundaryPolygon.size());

    constexpr double tolerance = 1e-5;

    ASSERT_NEAR(227.00204467773401, meshBoundaryPolygon[0].x, tolerance);
    ASSERT_NEAR(259.25222778320301, meshBoundaryPolygon[1].x, tolerance);
    ASSERT_NEAR(428.00329589843801, meshBoundaryPolygon[2].x, tolerance);
    ASSERT_NEAR(536.00396728515602, meshBoundaryPolygon[3].x, tolerance);
    ASSERT_NEAR(503.75378417968801, meshBoundaryPolygon[4].x, tolerance);
    ASSERT_NEAR(350.75280761718801, meshBoundaryPolygon[5].x, tolerance);
    ASSERT_NEAR(322.25262451171898, meshBoundaryPolygon[6].x, tolerance);
    ASSERT_NEAR(227.00204467773401, meshBoundaryPolygon[7].x, tolerance);

    ASSERT_NEAR(360.37924194335898, meshBoundaryPolygon[0].y, tolerance);
    ASSERT_NEAR(241.87805175781301, meshBoundaryPolygon[1].y, tolerance);
    ASSERT_NEAR(210.37774658203099, meshBoundaryPolygon[2].y, tolerance);
    ASSERT_NEAR(310.87875366210898, meshBoundaryPolygon[3].y, tolerance);
    ASSERT_NEAR(432.37997436523398, meshBoundaryPolygon[4].y, tolerance);
    ASSERT_NEAR(458.63024902343801, meshBoundaryPolygon[5].y, tolerance);
    ASSERT_NEAR(454.88018798828102, meshBoundaryPolygon[6].y, tolerance);
    ASSERT_NEAR(360.37924194335898, meshBoundaryPolygon[7].y, tolerance);
}

TEST(Polygons, CreatePointsInPolygons)
{
    // Prepare

    std::vector<meshkernel::Point> nodes;

    nodes.push_back({302.002502, 472.130371});
    nodes.push_back({144.501526, 253.128174});
    nodes.push_back({368.752930, 112.876755});
    nodes.push_back({707.755005, 358.879242});
    nodes.push_back({301.252502, 471.380371});
    nodes.push_back({302.002502, 472.130371});

    meshkernel::Polygons polygons(nodes, meshkernel::Projection::cartesian);

    // Execute
    const auto generatedPoints = polygons.ComputePointsInPolygons();

    // Assert
    constexpr double tolerance = 1e-5;

    ASSERT_NEAR(302.00250199999999, generatedPoints[0][0].x, tolerance);
    ASSERT_NEAR(472.13037100000003, generatedPoints[0][0].y, tolerance);

    ASSERT_NEAR(144.50152600000001, generatedPoints[0][1].x, tolerance);
    ASSERT_NEAR(253.12817400000000, generatedPoints[0][1].y, tolerance);

    ASSERT_NEAR(368.75292999999999, generatedPoints[0][2].x, tolerance);
    ASSERT_NEAR(112.87675500000000, generatedPoints[0][2].y, tolerance);

    ASSERT_NEAR(707.75500499999998, generatedPoints[0][3].x, tolerance);
    ASSERT_NEAR(358.87924199999998, generatedPoints[0][3].y, tolerance);

    ASSERT_NEAR(301.25250199999999, generatedPoints[0][4].x, tolerance);
    ASSERT_NEAR(471.38037100000003, generatedPoints[0][4].y, tolerance);
}

TEST(Polygons, RefineDefaultPolygon)
{
    meshkernel::Polygons polygon;
    // Should fail, polygon is empty.
    EXPECT_THROW([[maybe_unused]] auto result = polygon.RefineFirstPolygon(0, 0, 1.0), meshkernel::ConstraintError);
    EXPECT_THROW([[maybe_unused]] auto result = polygon.RefineFirstPolygon(0, 1, 1.0), meshkernel::ConstraintError);
}

TEST(Polygons, InvalidRefinePolygonIndex)
{
    std::vector<meshkernel::Point> nodes;

    nodes.push_back({0, 0});
    nodes.push_back({3, 0});
    nodes.push_back({3, 3});
    nodes.push_back({0, 3});
    nodes.push_back({0, 0});

    meshkernel::Polygons polygon(nodes, meshkernel::Projection::cartesian);

    // First is greater than last
    EXPECT_THROW([[maybe_unused]] auto result = polygon.RefineFirstPolygon(10, 8, 1.0), meshkernel::ConstraintError);
    // Last index is out of range
    EXPECT_THROW([[maybe_unused]] auto result = polygon.RefineFirstPolygon(3, 9, 1.0), meshkernel::ConstraintError);
}

TEST(Polygons, RefinePolygon)
{
    // Prepare
    std::vector<meshkernel::Point> nodes;

    nodes.push_back({0, 0});
    nodes.push_back({3, 0});
    nodes.push_back({3, 3});
    nodes.push_back({0, 3});
    nodes.push_back({0, 0});

    meshkernel::Polygons polygons(nodes, meshkernel::Projection::cartesian);

    // Execute
    const auto refinedPolygon = polygons.RefineFirstPolygon(0, 0, 1.0);

    ASSERT_EQ(13, refinedPolygon.size());
    constexpr double tolerance = 1e-5;

    ASSERT_NEAR(0, refinedPolygon[0].x, tolerance);
    ASSERT_NEAR(1, refinedPolygon[1].x, tolerance);
    ASSERT_NEAR(2, refinedPolygon[2].x, tolerance);
    ASSERT_NEAR(3, refinedPolygon[3].x, tolerance);
    ASSERT_NEAR(3, refinedPolygon[4].x, tolerance);
    ASSERT_NEAR(3, refinedPolygon[5].x, tolerance);
    ASSERT_NEAR(3, refinedPolygon[6].x, tolerance);
    ASSERT_NEAR(2, refinedPolygon[7].x, tolerance);
    ASSERT_NEAR(1, refinedPolygon[8].x, tolerance);
    ASSERT_NEAR(0, refinedPolygon[9].x, tolerance);
    ASSERT_NEAR(0, refinedPolygon[10].x, tolerance);
    ASSERT_NEAR(0, refinedPolygon[11].x, tolerance);
    ASSERT_NEAR(0, refinedPolygon[12].x, tolerance);

    ASSERT_NEAR(0, refinedPolygon[0].y, tolerance);
    ASSERT_NEAR(0, refinedPolygon[1].y, tolerance);
    ASSERT_NEAR(0, refinedPolygon[2].y, tolerance);
    ASSERT_NEAR(0, refinedPolygon[3].y, tolerance);
    ASSERT_NEAR(1, refinedPolygon[4].y, tolerance);
    ASSERT_NEAR(2, refinedPolygon[5].y, tolerance);
    ASSERT_NEAR(3, refinedPolygon[6].y, tolerance);
    ASSERT_NEAR(3, refinedPolygon[7].y, tolerance);
    ASSERT_NEAR(3, refinedPolygon[8].y, tolerance);
    ASSERT_NEAR(3, refinedPolygon[9].y, tolerance);
    ASSERT_NEAR(2, refinedPolygon[10].y, tolerance);
    ASSERT_NEAR(1, refinedPolygon[11].y, tolerance);
    ASSERT_NEAR(0, refinedPolygon[12].y, tolerance);
}

TEST(Polygons, RefinePolygonTwiceWithSameRefinement)
{
    // Prepare
    std::vector<meshkernel::Point> nodes;

    nodes.push_back({0, 0});
    nodes.push_back({3, 0});
    nodes.push_back({3, 3});
    nodes.push_back({0, 3});
    nodes.push_back({0, 0});

    meshkernel::Polygons polygons(nodes, meshkernel::Projection::cartesian);

    // Execute
    const auto refinedPolygon = polygons.RefineFirstPolygon(0, 0, 1.0);

    meshkernel::Polygons polygons2(refinedPolygon, meshkernel::Projection::cartesian);
    const auto refinedPolygon2 = polygons2.RefineFirstPolygon(0, 0, 1.0);

    constexpr double tolerance = 1e-5;

    // Only need to check the points from the second refinement, the test above will
    // catch any problems in the first refinement.
    ASSERT_EQ(13, refinedPolygon2.size());

    ASSERT_NEAR(0, refinedPolygon2[0].x, tolerance);
    ASSERT_NEAR(1, refinedPolygon2[1].x, tolerance);
    ASSERT_NEAR(2, refinedPolygon2[2].x, tolerance);
    ASSERT_NEAR(3, refinedPolygon2[3].x, tolerance);
    ASSERT_NEAR(3, refinedPolygon2[4].x, tolerance);
    ASSERT_NEAR(3, refinedPolygon2[5].x, tolerance);
    ASSERT_NEAR(3, refinedPolygon2[6].x, tolerance);
    ASSERT_NEAR(2, refinedPolygon2[7].x, tolerance);
    ASSERT_NEAR(1, refinedPolygon2[8].x, tolerance);
    ASSERT_NEAR(0, refinedPolygon2[9].x, tolerance);
    ASSERT_NEAR(0, refinedPolygon2[10].x, tolerance);
    ASSERT_NEAR(0, refinedPolygon2[11].x, tolerance);
    ASSERT_NEAR(0, refinedPolygon2[12].x, tolerance);

    ASSERT_NEAR(0, refinedPolygon2[0].y, tolerance);
    ASSERT_NEAR(0, refinedPolygon2[1].y, tolerance);
    ASSERT_NEAR(0, refinedPolygon2[2].y, tolerance);
    ASSERT_NEAR(0, refinedPolygon2[3].y, tolerance);
    ASSERT_NEAR(1, refinedPolygon2[4].y, tolerance);
    ASSERT_NEAR(2, refinedPolygon2[5].y, tolerance);
    ASSERT_NEAR(3, refinedPolygon2[6].y, tolerance);
    ASSERT_NEAR(3, refinedPolygon2[7].y, tolerance);
    ASSERT_NEAR(3, refinedPolygon2[8].y, tolerance);
    ASSERT_NEAR(3, refinedPolygon2[9].y, tolerance);
    ASSERT_NEAR(2, refinedPolygon2[10].y, tolerance);
    ASSERT_NEAR(1, refinedPolygon2[11].y, tolerance);
    ASSERT_NEAR(0, refinedPolygon2[12].y, tolerance);
}

TEST(Polygons, RefinePolygonTwiceWithLargerRefinement)
{
    // Prepare
    std::vector<meshkernel::Point> nodes;

    nodes.push_back({0, 0});
    nodes.push_back({3, 0});
    nodes.push_back({3, 3});
    nodes.push_back({0, 3});
    nodes.push_back({0, 0});

    meshkernel::Polygons polygons(nodes, meshkernel::Projection::cartesian);

    // Execute
    const auto refinedPolygon = polygons.RefineFirstPolygon(0, 0, 1.0);

    meshkernel::Polygons polygons2(refinedPolygon, meshkernel::Projection::cartesian);
    const auto refinedPolygon2 = polygons2.RefineFirstPolygon(0, 0, 2.0);

    constexpr double tolerance = 1e-13;

    // Only need to check the points from the second refinement, the test above will
    // catch any problems in the first refinement.
    ASSERT_EQ(13, refinedPolygon.size());
    ASSERT_EQ(refinedPolygon.size(), refinedPolygon2.size());

    for (size_t i = 0; i < refinedPolygon.size(); ++i)
    {
        EXPECT_NEAR(refinedPolygon[i].x, refinedPolygon2[i].x, tolerance);
        EXPECT_NEAR(refinedPolygon[i].y, refinedPolygon2[i].y, tolerance);
    }
}

TEST(Polygons, RefinePolygonTwice)
{
    // Prepare
    std::vector<meshkernel::Point> nodes;

    nodes.push_back({0, 0});
    nodes.push_back({3, 0});
    nodes.push_back({3, 3});
    nodes.push_back({0, 3});
    nodes.push_back({0, 0});

    meshkernel::Polygons polygons(nodes, meshkernel::Projection::cartesian);

    // Execute
    const auto refinedPolygon = polygons.RefineFirstPolygon(0, 0, 1.0);

    meshkernel::Polygons polygons2(refinedPolygon, meshkernel::Projection::cartesian);
    const auto refinedPolygon2 = polygons2.RefineFirstPolygon(0, 0, 0.5);

    constexpr double tolerance = 1e-5;

    // Only need to check the points from the second refinement, the test above will
    // catch any problems in the first refinement.
    ASSERT_EQ(25, refinedPolygon2.size());

    ASSERT_NEAR(0, refinedPolygon2[0].x, tolerance);
    ASSERT_NEAR(0.5, refinedPolygon2[1].x, tolerance);
    ASSERT_NEAR(1, refinedPolygon2[2].x, tolerance);
    ASSERT_NEAR(1.5, refinedPolygon2[3].x, tolerance);
    ASSERT_NEAR(2, refinedPolygon2[4].x, tolerance);
    ASSERT_NEAR(2.5, refinedPolygon2[5].x, tolerance);
    ASSERT_NEAR(3, refinedPolygon2[6].x, tolerance);
    ASSERT_NEAR(3, refinedPolygon2[7].x, tolerance);
    ASSERT_NEAR(3, refinedPolygon2[8].x, tolerance);
    ASSERT_NEAR(3, refinedPolygon2[9].x, tolerance);
    ASSERT_NEAR(3, refinedPolygon2[10].x, tolerance);
    ASSERT_NEAR(3, refinedPolygon2[11].x, tolerance);
    ASSERT_NEAR(3, refinedPolygon2[12].x, tolerance);
    ASSERT_NEAR(2.5, refinedPolygon2[13].x, tolerance);
    ASSERT_NEAR(2, refinedPolygon2[14].x, tolerance);
    ASSERT_NEAR(1.5, refinedPolygon2[15].x, tolerance);
    ASSERT_NEAR(1, refinedPolygon2[16].x, tolerance);
    ASSERT_NEAR(0.5, refinedPolygon2[17].x, tolerance);
    ASSERT_NEAR(0, refinedPolygon2[18].x, tolerance);
    ASSERT_NEAR(0, refinedPolygon2[19].x, tolerance);
    ASSERT_NEAR(0, refinedPolygon2[20].x, tolerance);
    ASSERT_NEAR(0, refinedPolygon2[21].x, tolerance);
    ASSERT_NEAR(0, refinedPolygon2[22].x, tolerance);
    ASSERT_NEAR(0, refinedPolygon2[23].x, tolerance);
    ASSERT_NEAR(0, refinedPolygon2[24].x, tolerance);

    ASSERT_NEAR(0, refinedPolygon2[0].y, tolerance);
    ASSERT_NEAR(0, refinedPolygon2[1].y, tolerance);
    ASSERT_NEAR(0, refinedPolygon2[2].y, tolerance);
    ASSERT_NEAR(0, refinedPolygon2[3].y, tolerance);
    ASSERT_NEAR(0, refinedPolygon2[4].y, tolerance);
    ASSERT_NEAR(0, refinedPolygon2[5].y, tolerance);
    ASSERT_NEAR(0, refinedPolygon2[6].y, tolerance);
    ASSERT_NEAR(0.5, refinedPolygon2[7].y, tolerance);
    ASSERT_NEAR(1, refinedPolygon2[8].y, tolerance);
    ASSERT_NEAR(1.5, refinedPolygon2[9].y, tolerance);
    ASSERT_NEAR(2, refinedPolygon2[10].y, tolerance);
    ASSERT_NEAR(2.5, refinedPolygon2[11].y, tolerance);
    ASSERT_NEAR(3, refinedPolygon2[12].y, tolerance);
    ASSERT_NEAR(3, refinedPolygon2[13].y, tolerance);
    ASSERT_NEAR(3, refinedPolygon2[14].y, tolerance);
    ASSERT_NEAR(3, refinedPolygon2[15].y, tolerance);
    ASSERT_NEAR(3, refinedPolygon2[16].y, tolerance);
    ASSERT_NEAR(3, refinedPolygon2[17].y, tolerance);
    ASSERT_NEAR(3, refinedPolygon2[18].y, tolerance);
    ASSERT_NEAR(2.5, refinedPolygon2[19].y, tolerance);
    ASSERT_NEAR(2, refinedPolygon2[20].y, tolerance);
    ASSERT_NEAR(1.5, refinedPolygon2[21].y, tolerance);
    ASSERT_NEAR(1, refinedPolygon2[22].y, tolerance);
    ASSERT_NEAR(0.5, refinedPolygon2[23].y, tolerance);
    ASSERT_NEAR(0, refinedPolygon2[24].y, tolerance);
}

TEST(Polygons, RefinePolygonOneSide)
{
    // Prepare
    std::vector<meshkernel::Point> nodes;

    nodes.push_back({0, 0});
    nodes.push_back({3, 0});
    nodes.push_back({3, 3});
    nodes.push_back({0, 3});
    nodes.push_back({0, 0});

    meshkernel::Polygons polygons(nodes, meshkernel::Projection::cartesian);

    // Execute
    const auto refinedPolygon = polygons.RefineFirstPolygon(0, 1, 1.0);

    ASSERT_EQ(7, refinedPolygon.size());
    constexpr double tolerance = 1e-5;

    ASSERT_NEAR(0.0, refinedPolygon[0].x, tolerance);
    ASSERT_NEAR(1.0, refinedPolygon[1].x, tolerance);
    ASSERT_NEAR(2.0, refinedPolygon[2].x, tolerance);
    ASSERT_NEAR(3.0, refinedPolygon[3].x, tolerance);
    ASSERT_NEAR(3.0, refinedPolygon[4].x, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[5].x, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[6].x, tolerance);

    ASSERT_NEAR(0.0, refinedPolygon[0].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[1].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[2].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[3].y, tolerance);
    ASSERT_NEAR(3.0, refinedPolygon[4].y, tolerance);
    ASSERT_NEAR(3.0, refinedPolygon[5].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[6].y, tolerance);
}

TEST(Polygons, RefinePolygonTwoTimesOneSideSameRefinement)
{
    // Prepare
    std::vector<meshkernel::Point> nodes;

    nodes.push_back({0, 0});
    nodes.push_back({3, 0});
    nodes.push_back({3, 3});
    nodes.push_back({0, 3});
    nodes.push_back({0, 0});

    meshkernel::Polygons polygons(nodes, meshkernel::Projection::cartesian);

    // Execute
    const auto refinedPolygon = polygons.RefineFirstPolygon(0, 1, 1.0);

    meshkernel::Polygons polygons2(refinedPolygon, meshkernel::Projection::cartesian);
    // Now there are additional segments in the polygon, the segments to refine is different.
    const auto refinedPolygon2 = polygons2.RefineFirstPolygon(0, 3, 1.0);

    ASSERT_EQ(7, refinedPolygon.size());
    constexpr double tolerance = 1e-10;

    ASSERT_NEAR(0.0, refinedPolygon[0].x, tolerance);
    ASSERT_NEAR(1.0, refinedPolygon[1].x, tolerance);
    ASSERT_NEAR(2.0, refinedPolygon[2].x, tolerance);
    ASSERT_NEAR(3.0, refinedPolygon[3].x, tolerance);
    ASSERT_NEAR(3.0, refinedPolygon[4].x, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[5].x, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[6].x, tolerance);

    ASSERT_NEAR(0.0, refinedPolygon[0].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[1].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[2].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[3].y, tolerance);
    ASSERT_NEAR(3.0, refinedPolygon[4].y, tolerance);
    ASSERT_NEAR(3.0, refinedPolygon[5].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[6].y, tolerance);

    ASSERT_EQ(refinedPolygon.size(), refinedPolygon2.size());

    for (size_t i = 0; i < refinedPolygon.size(); ++i)
    {
        ASSERT_NEAR(refinedPolygon[i].x, refinedPolygon2[i].x, tolerance);
        ASSERT_NEAR(refinedPolygon[i].y, refinedPolygon2[i].y, tolerance);
    }
}

TEST(Polygons, RefinePolygonTwoTimesOneSide)
{
    // Prepare
    std::vector<meshkernel::Point> nodes;

    nodes.push_back({0, 0});
    nodes.push_back({4, 0});
    nodes.push_back({4, 4});
    nodes.push_back({0, 4});
    nodes.push_back({0, 0});

    meshkernel::Polygons polygons(nodes, meshkernel::Projection::cartesian);

    // Execute
    const auto refinedPolygon = polygons.RefineFirstPolygon(1, 2, 2.0);

    meshkernel::Polygons polygons2(refinedPolygon, meshkernel::Projection::cartesian);

    // Refine the same edge, this time there should be two segments making up the edge.
    const auto refinedPolygon2 = polygons2.RefineFirstPolygon(1, 3, 1.0);

    constexpr double tolerance = 1.0e-13;

    // Check points from first refinement
    ASSERT_NEAR(0.0, refinedPolygon[0].x, tolerance);
    ASSERT_NEAR(4.0, refinedPolygon[1].x, tolerance);
    ASSERT_NEAR(4.0, refinedPolygon[2].x, tolerance);
    ASSERT_NEAR(4.0, refinedPolygon[3].x, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[4].x, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[5].x, tolerance);

    ASSERT_NEAR(0.0, refinedPolygon[0].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[1].y, tolerance);
    ASSERT_NEAR(2.0, refinedPolygon[2].y, tolerance);
    ASSERT_NEAR(4.0, refinedPolygon[3].y, tolerance);
    ASSERT_NEAR(4.0, refinedPolygon[4].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[5].y, tolerance);

    // Check points from second refinement
    ASSERT_NEAR(0.0, refinedPolygon2[0].x, tolerance);
    ASSERT_NEAR(4.0, refinedPolygon2[1].x, tolerance);
    ASSERT_NEAR(4.0, refinedPolygon2[2].x, tolerance);
    ASSERT_NEAR(4.0, refinedPolygon2[3].x, tolerance);
    ASSERT_NEAR(4.0, refinedPolygon2[4].x, tolerance);
    ASSERT_NEAR(4.0, refinedPolygon2[5].x, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon2[6].x, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon2[7].x, tolerance);

    ASSERT_NEAR(0.0, refinedPolygon2[0].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon2[1].y, tolerance);
    ASSERT_NEAR(1.0, refinedPolygon2[2].y, tolerance);
    ASSERT_NEAR(2.0, refinedPolygon2[3].y, tolerance);
    ASSERT_NEAR(3.0, refinedPolygon2[4].y, tolerance);
    ASSERT_NEAR(4.0, refinedPolygon2[5].y, tolerance);
    ASSERT_NEAR(4.0, refinedPolygon2[6].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon2[7].y, tolerance);
}

TEST(Polygons, RefinePolygonLongerSquare)
{
    // Prepare
    std::vector<meshkernel::Point> nodes;

    nodes.push_back({0, 0});
    nodes.push_back({3, 0});
    nodes.push_back({3, 3});
    nodes.push_back({3.5, 0});
    nodes.push_back({0, 0});

    meshkernel::Polygons polygons(nodes, meshkernel::Projection::cartesian);

    // Execute
    const auto refinedPolygon = polygons.RefineFirstPolygon(0, 0, 1.0);

    ASSERT_EQ(15, refinedPolygon.size());
    constexpr double tolerance = 1e-5;

    ASSERT_NEAR(0.0, refinedPolygon[0].x, tolerance);
    ASSERT_NEAR(1.0, refinedPolygon[1].x, tolerance);
    ASSERT_NEAR(2.0, refinedPolygon[2].x, tolerance);
    ASSERT_NEAR(3.0, refinedPolygon[3].x, tolerance);
    ASSERT_NEAR(3.0, refinedPolygon[4].x, tolerance);
    ASSERT_NEAR(3.0, refinedPolygon[5].x, tolerance);
    ASSERT_NEAR(3.0, refinedPolygon[6].x, tolerance);
    ASSERT_NEAR(3.1643989873053573, refinedPolygon[7].x, tolerance);
    ASSERT_NEAR(3.3287979746107146, refinedPolygon[8].x, tolerance);
    ASSERT_NEAR(3.4931969619160719, refinedPolygon[9].x, tolerance);
    ASSERT_NEAR(3.5, refinedPolygon[10].x, tolerance);
    ASSERT_NEAR(2.5, refinedPolygon[11].x, tolerance);
    ASSERT_NEAR(1.5, refinedPolygon[12].x, tolerance);
    ASSERT_NEAR(0.5, refinedPolygon[13].x, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[14].x, tolerance);

    ASSERT_NEAR(0.0, refinedPolygon[0].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[1].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[2].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[3].y, tolerance);
    ASSERT_NEAR(1.0, refinedPolygon[4].y, tolerance);
    ASSERT_NEAR(2.0, refinedPolygon[5].y, tolerance);
    ASSERT_NEAR(3.0, refinedPolygon[6].y, tolerance);
    ASSERT_NEAR(2.0136060761678563, refinedPolygon[7].y, tolerance);
    ASSERT_NEAR(1.0272121523357125, refinedPolygon[8].y, tolerance);
    ASSERT_NEAR(0.040818228503568754, refinedPolygon[9].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[10].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[11].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[12].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[13].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[14].y, tolerance);
}

TEST(Polygons, OffsetCopy)
{
    std::vector<meshkernel::Point> nodes;
    nodes.push_back({296.752472, 397.879639});
    nodes.push_back({294.502472, 256.128204});
    nodes.push_back({578.754211, 244.128082});
    nodes.push_back({587.754272, 400.129639});
    nodes.push_back({308.002533, 397.879639});
    nodes.push_back({296.752472, 397.879639});

    meshkernel::Polygons polygon(nodes, meshkernel::Projection::cartesian);

    double distance = 10.0;
    bool innerAndOuter = false;
    const auto newPolygon = polygon.OffsetCopy(distance, innerAndOuter);

    constexpr double tolerance = 1e-5;

    ASSERT_NEAR(newPolygon.Node(0).x, 286.75373149966771, tolerance);
    ASSERT_NEAR(newPolygon.Node(1).x, 284.34914611880089, tolerance);
    ASSERT_NEAR(newPolygon.Node(2).x, 588.17047010011993, tolerance);
    ASSERT_NEAR(newPolygon.Node(3).x, 598.35275776004642, tolerance);
    ASSERT_NEAR(newPolygon.Node(4).x, 307.96231942308754, tolerance);
    ASSERT_NEAR(newPolygon.Node(5).x, 296.75247200000001, tolerance);

    ASSERT_NEAR(newPolygon.Node(0).y, 398.03834755999270, tolerance);
    ASSERT_NEAR(newPolygon.Node(1).y, 246.54793497426144, tolerance);
    ASSERT_NEAR(newPolygon.Node(2).y, 233.72165300742589, tolerance);
    ASSERT_NEAR(newPolygon.Node(3).y, 410.21520441451258, tolerance);
    ASSERT_NEAR(newPolygon.Node(4).y, 407.87963900000000, tolerance);
    ASSERT_NEAR(newPolygon.Node(5).y, 407.87963900000000, tolerance);
}

TEST(Polygons, SnapSinglePolygonToSingleLandBoundary)
{
    // Test the algorithm for snapping single polygon to land boundaries.

    constexpr double tolerance = 1.0e-8;

    // The land boundary to which the polygon is to be snapped.
    std::vector<meshkernel::Point> landBoundaryPoints{{139.251465, 497.630615},
                                                      {527.753906, 499.880676},
                                                      {580.254211, 265.878296},
                                                      {194.001801, 212.627762}};

    meshkernel::LandBoundary landBoundary(landBoundaryPoints);

    // The original polygon points.
    std::vector<meshkernel::Point> polygonPoints{{170.001648, 472.880371},
                                                 {263.002228, 472.880371},
                                                 {344.002747, 475.130432},
                                                 {458.753448, 482.630493},
                                                 {515.753845, 487.130554},
                                                 {524.753906, 434.630005},
                                                 {510.503754, 367.129333},
                                                 {557.754089, 297.378601},
                                                 {545.004028, 270.378357},
                                                 {446.003387, 259.128235},
                                                 {340.252716, 244.128067},
                                                 {242.752106, 226.877884}};

    meshkernel::Polygons polygon(polygonPoints, meshkernel::Projection::cartesian);

    // The expected polygon values after snapping to land boundary.
    std::vector<meshkernel::Point> expectedSnappedPoint{{169.8572772242283, 497.8078724305628},
                                                        {262.854737816309, 498.3464789799546},
                                                        {343.8655709877979, 498.8156634613377},
                                                        {458.6558591358565, 499.4804859264834},
                                                        {515.6804060372598, 499.8107507986815},
                                                        {541.5480568270806, 438.3979070214996},
                                                        {555.2836667233159, 377.1760644727631},
                                                        {572.4472626165707, 300.6751319852315},
                                                        {546.2703464583593, 261.1931241088368},
                                                        {447.5942143903486, 247.589178632675},
                                                        {341.7865993173012, 233.0020541046851},
                                                        {243.7707524316129, 219.4891385810638}};

    // Snap the polygon to the land boundary
    polygon.SnapToLandBoundary(landBoundary, 0, static_cast<meshkernel::UInt>(polygonPoints.size() - 1));

    for (meshkernel::UInt i = 0; i < polygonPoints.size(); ++i)
    {
        // Use relative tolerance to compare expected with actual values.
        EXPECT_TRUE(meshkernel::IsEqual(expectedSnappedPoint[i].x, polygon.Node(i).x, tolerance))
            << " Expected x-coordinate: " << expectedSnappedPoint[i].x << ", actual: " << polygon.Node(i).x << ", tolerance: " << tolerance;
        EXPECT_TRUE(meshkernel::IsEqual(expectedSnappedPoint[i].y, polygon.Node(i).y, tolerance))
            << " Expected y-coordinate: " << expectedSnappedPoint[i].y << ", actual: " << polygon.Node(i).y << ", tolerance: " << tolerance;
    }
}

TEST(Polygons, SnapMultiPolygonToSingleLandBoundary)
{
    // Test the algorithm for snapping multi-polygons to land boundaries.

    constexpr double tolerance = 1.0e-8;

    // The land boundary to which the polygon is to be snapped.
    std::vector<meshkernel::Point> landBoundaryPoints{{139.251465, 497.630615},
                                                      {527.753906, 499.880676},
                                                      {580.254211, 265.878296},
                                                      {194.001801, 212.627762}};

    meshkernel::LandBoundary landBoundary(landBoundaryPoints);

    // The original polygon points.
    std::vector<meshkernel::Point> polygonPoints{{170.001648, 472.880371},
                                                 {263.002228, 472.880371},
                                                 {344.002747, 475.130432},
                                                 {458.753448, 482.630493},
                                                 {515.753845, 487.130554},
                                                 {524.753906, 434.630005},
                                                 {meshkernel::constants::missing::doubleValue, meshkernel::constants::missing::doubleValue},
                                                 {510.503754, 367.129333},
                                                 {557.754089, 297.378601},
                                                 {545.004028, 270.378357},
                                                 {446.003387, 259.128235},
                                                 {340.252716, 244.128067},
                                                 {242.752106, 226.877884}};

    meshkernel::Polygons polygon(polygonPoints, meshkernel::Projection::cartesian);

    // The expected polygon values after snapping to land boundary.
    std::vector<meshkernel::Point> expectedSnappedPoint{{169.8572772242283, 497.8078724305628},
                                                        {262.854737816309, 498.3464789799546},
                                                        {343.8655709877979, 498.8156634613377},
                                                        {458.6558591358565, 499.4804859264834},
                                                        {515.6804060372598, 499.8107507986815},
                                                        {541.5480568270806, 438.3979070214996},
                                                        {meshkernel::constants::missing::doubleValue, meshkernel::constants::missing::doubleValue},
                                                        {555.2836667233159, 377.1760644727631},
                                                        {572.4472626165707, 300.6751319852315},
                                                        {546.2703464583593, 261.1931241088368},
                                                        {447.5942143903486, 247.589178632675},
                                                        {341.7865993173012, 233.0020541046851},
                                                        {243.7707524316129, 219.4891385810638}};

    // Snap the polygon to the land boundary
    polygon.SnapToLandBoundary(landBoundary, 0, static_cast<meshkernel::UInt>(polygonPoints.size() - 1));

    for (meshkernel::UInt i = 0; i < polygonPoints.size(); ++i)
    {
        // Use relative tolerance to compare expected with actual values.
        EXPECT_TRUE(meshkernel::IsEqual(expectedSnappedPoint[i].x, polygon.Node(i).x, tolerance))
            << " Expected x-coordinate: " << expectedSnappedPoint[i].x << ", actual: " << polygon.Node(i).x << ", tolerance: " << tolerance;
        EXPECT_TRUE(meshkernel::IsEqual(expectedSnappedPoint[i].y, polygon.Node(i).y, tolerance))
            << " Expected y-coordinate: " << expectedSnappedPoint[i].y << ", actual: " << polygon.Node(i).y << ", tolerance: " << tolerance;
    }
}

TEST(Polygons, SnapMultiPolygonToMultiLandBoundary)
{
    // Test the algorithm for snapping multi-polygons to land boundaries.

    constexpr double tolerance = 1.0e-8;

    // The land boundary to which the polygon is to be snapped.
    // This land boundary is made up of two sections
    std::vector<meshkernel::Point> landBoundaryPoints{{0.0, 1.0},
                                                      {0.0, 0.0},
                                                      {1.0, 0.0},
                                                      {meshkernel::constants::missing::doubleValue, meshkernel::constants::missing::doubleValue},
                                                      {2.0, 0.0},
                                                      {3.0, 0.0},
                                                      {3.0, 1.0}};

    meshkernel::LandBoundary landBoundary(landBoundaryPoints);

    // The original polygon points.
    // The points make up two distinct polygons
    std::vector<meshkernel::Point> polygonPoints{{0.2, 1.1},
                                                 {-0.1, 0.8},
                                                 {0.1, 0.4},
                                                 {0.01, 0.2},
                                                 {0.2, 0.05},
                                                 {0.4, -0.1},
                                                 {0.5, -0.02},
                                                 {1.1, 0.1},
                                                 {meshkernel::constants::missing::doubleValue, meshkernel::constants::missing::doubleValue},
                                                 {1.9, -0.1},
                                                 {2.1, 0.1},
                                                 {2.5, 0.04},
                                                 {3.1, -0.1},
                                                 {2.95, 0.1},
                                                 {2.9, 0.2},
                                                 {3.1, 0.5},
                                                 {2.95, 0.9}};

    meshkernel::Polygons polygon(polygonPoints, meshkernel::Projection::cartesian);

    // The expected polygon values after snapping to land boundary.
    std::vector<meshkernel::Point> expectedSnappedPoint{{0.0, 1.0},
                                                        {0.0, 0.8},
                                                        {0.0, 0.4},
                                                        {0.0, 0.2},
                                                        {0.2, 0.0},
                                                        {0.4, 0.0},
                                                        {0.5, 0.0},
                                                        {1.0, 0.0},
                                                        {meshkernel::constants::missing::doubleValue, meshkernel::constants::missing::doubleValue},
                                                        {2.0, 0.0},
                                                        {2.1, 0.0},
                                                        {2.5, 0.0},
                                                        {3.0, 0.0},
                                                        {3.0, 0.1},
                                                        {3.0, 0.2},
                                                        {3.0, 0.5},
                                                        {3.0, 0.9}};

    // Snap the polygon to the land boundary
    polygon.SnapToLandBoundary(landBoundary, 0, 0);
    // polygon.SnapToLandBoundary(landBoundary, 0, static_cast<meshkernel::UInt>(polygonPoints.size() - 1));

    for (meshkernel::UInt i = 0; i < polygonPoints.size(); ++i)
    {
        // Use relative tolerance to compare expected with actual values.
        EXPECT_TRUE(meshkernel::IsEqual(expectedSnappedPoint[i].x, polygon.Node(i).x, tolerance))
            << " Expected x-coordinate: " << expectedSnappedPoint[i].x << ", actual: " << polygon.Node(i).x << ", tolerance: " << tolerance;
        EXPECT_TRUE(meshkernel::IsEqual(expectedSnappedPoint[i].y, polygon.Node(i).y, tolerance))
            << " Expected y-coordinate: " << expectedSnappedPoint[i].y << ", actual: " << polygon.Node(i).y << ", tolerance: " << tolerance;
    }
}

TEST(Polygons, SnapMultiPolygonPartToSingleLandBoundary)
{
    // Test the algorithm for snapping multi-polygons to land boundaries.
    // Snaps only the first part of the polygon to the land boundary

    constexpr double tolerance = 1.0e-8;

    // The land boundary to which the polygon is to be snapped.
    std::vector<meshkernel::Point> landBoundaryPoints{{139.251465, 497.630615},
                                                      {527.753906, 499.880676},
                                                      {580.254211, 265.878296},
                                                      {194.001801, 212.627762}};

    meshkernel::LandBoundary landBoundary(landBoundaryPoints);

    // The original polygon points.
    std::vector<meshkernel::Point> polygonPoints{{170.001648, 472.880371},
                                                 {263.002228, 472.880371},
                                                 {344.002747, 475.130432},
                                                 {458.753448, 482.630493},
                                                 {515.753845, 487.130554},
                                                 {524.753906, 434.630005},
                                                 {meshkernel::constants::missing::doubleValue, meshkernel::constants::missing::doubleValue},
                                                 {510.503754, 367.129333},
                                                 {557.754089, 297.378601},
                                                 {545.004028, 270.378357},
                                                 {446.003387, 259.128235},
                                                 {340.252716, 244.128067},
                                                 {242.752106, 226.877884}};

    meshkernel::Polygons polygon(polygonPoints, meshkernel::Projection::cartesian);

    // The expected polygon values after snapping to land boundary.
    std::vector<meshkernel::Point> expectedSnappedPoint{{169.8572772242283, 497.8078724305628},
                                                        {262.854737816309, 498.3464789799546},
                                                        {343.8655709877979, 498.8156634613377},
                                                        {458.6558591358565, 499.4804859264834},
                                                        {515.6804060372598, 499.8107507986815},
                                                        {541.5480568270806, 438.3979070214996},
                                                        {meshkernel::constants::missing::doubleValue, meshkernel::constants::missing::doubleValue},
                                                        {510.503754, 367.129333},
                                                        {557.754089, 297.378601},
                                                        {545.004028, 270.378357},
                                                        {446.003387, 259.128235},
                                                        {340.252716, 244.128067},
                                                        {242.752106, 226.877884}};

    // Snap the polygon to the land boundary
    polygon.SnapToLandBoundary(landBoundary, 0, 5);

    for (meshkernel::UInt i = 0; i < polygonPoints.size(); ++i)
    {
        // Use relative tolerance to compare expected with actual values.
        EXPECT_TRUE(meshkernel::IsEqual(expectedSnappedPoint[i].x, polygon.Node(i).x, tolerance))
            << " Expected x-coordinate: " << expectedSnappedPoint[i].x << ", actual: " << polygon.Node(i).x << ", tolerance: " << tolerance;
        EXPECT_TRUE(meshkernel::IsEqual(expectedSnappedPoint[i].y, polygon.Node(i).y, tolerance))
            << " Expected y-coordinate: " << expectedSnappedPoint[i].y << ", actual: " << polygon.Node(i).y << ", tolerance: " << tolerance;
    }
}
