#include <gtest/gtest.h>

#include "MeshKernel/BoundingBox.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/Point.hpp"
#include "MeshKernel/Polygon.hpp"

#include "TestUtils/MakeMeshes.hpp"

namespace mk = meshkernel;

TEST(PolygonTests, BasicConstruction)
{
    std::vector<mk::Point> polygonPoints{{0.0, 0.0}, {10.0, 0.0}, {10.0, 10.0}, {0.0, 10.0}, {0.0, 0.0}};

    mk::Polygon polygon(polygonPoints, mk::Projection::cartesian);

    EXPECT_EQ(polygon.Size(), polygonPoints.size());
    EXPECT_EQ(polygon.Nodes(), polygonPoints);
    EXPECT_EQ(polygon.GetProjection(), mk::Projection::cartesian);

    for (size_t i = 0; i < polygonPoints.size(); ++i)
    {
        EXPECT_TRUE(polygonPoints[i] == polygon.Node(i))
            << "Expected points to be same, expected: " << polygonPoints[i].x << ", " << polygonPoints[i].y
            << ", actual: " << polygon.Node(i).x << ", " << polygon.Node(i).y;
    }

    // Test copy constructor
    mk::Polygon polygonCopy(polygon);

    EXPECT_EQ(polygonCopy.Size(), polygonPoints.size());
    EXPECT_EQ(polygonCopy.Nodes(), polygonPoints);
    EXPECT_EQ(polygonCopy.GetProjection(), mk::Projection::cartesian);

    for (size_t i = 0; i < polygonPoints.size(); ++i)
    {
        EXPECT_TRUE(polygonPoints[i] == polygonCopy.Node(i))
            << "Expected points to be same, expected: " << polygonPoints[i].x << ", " << polygonPoints[i].y
            << ", actual: " << polygonCopy.Node(i).x << ", " << polygonCopy.Node(i).y;
    }

    mk::Polygon polygonAssign;

    EXPECT_EQ(polygonAssign.Size(), 0);

    // Test copy assignment
    polygonAssign = polygon;

    EXPECT_EQ(polygonAssign.Size(), polygonPoints.size());
    EXPECT_EQ(polygonAssign.Nodes(), polygonPoints);
    EXPECT_EQ(polygonAssign.GetProjection(), mk::Projection::cartesian);

    for (size_t i = 0; i < polygonPoints.size(); ++i)
    {
        EXPECT_TRUE(polygonPoints[i] == polygonAssign.Node(i))
            << "Expected points to be same, expected: " << polygonPoints[i].x << ", " << polygonPoints[i].y
            << ", actual: " << polygonAssign.Node(i).x << ", " << polygonAssign.Node(i).y;
    }
}

TEST(PolygonTests, BoundingBoxTest)
{
    std::vector<mk::Point> polygonPointsSquare{{0.0, 0.0}, {10.0, 0.0}, {10.0, 10.0}, {0.0, 10.0}, {0.0, 0.0}};

    mk::Polygon polygon(polygonPointsSquare, mk::Projection::cartesian);
    mk::BoundingBox boundingBox = polygon.GetBoundingBox();

    constexpr double tolerance = 1.0e-10;

    EXPECT_NEAR(boundingBox.lowerLeft().x, 0.0, tolerance);
    EXPECT_NEAR(boundingBox.lowerLeft().y, 0.0, tolerance);
    EXPECT_NEAR(boundingBox.upperRight().x, 10.0, tolerance);
    EXPECT_NEAR(boundingBox.upperRight().y, 10.0, tolerance);

    std::vector<mk::Point> polygonPointsLargerSquare{{-10.0, -5.0}, {20.0, -5.0}, {20.0, 10.0}, {-10.0, 10.0}, {-10.0, -5.0}};
    polygon.Reset(polygonPointsLargerSquare, mk::Projection::cartesian);

    boundingBox = polygon.GetBoundingBox();
    // The boundingbox should be reset too
    EXPECT_NEAR(boundingBox.lowerLeft().x, -10.0, tolerance);
    EXPECT_NEAR(boundingBox.lowerLeft().y, -5.0, tolerance);
    EXPECT_NEAR(boundingBox.upperRight().x, 20.0, tolerance);
    EXPECT_NEAR(boundingBox.upperRight().y, 10.0, tolerance);

    std::vector<mk::Point> polygonPointsPentagon{{-10.0, -5.0}, {20.0, -5.0}, {20.0, 10.0}, {5.0, 15.0}, {-10.0, 10.0}, {-10.0, -5.0}};
    polygon.Reset(polygonPointsPentagon, mk::Projection::cartesian);

    boundingBox = polygon.GetBoundingBox();
    // The boundingbox should be reset too
    EXPECT_NEAR(boundingBox.lowerLeft().x, -10.0, tolerance);
    EXPECT_NEAR(boundingBox.lowerLeft().y, -5.0, tolerance);
    EXPECT_NEAR(boundingBox.upperRight().x, 20.0, tolerance);
    EXPECT_NEAR(boundingBox.upperRight().y, 15.0, tolerance);
}

TEST(PolygonTests, ContainsTest)
{
    constexpr double min = 0.0;
    constexpr double max = 10.0;
    constexpr double extension = 2.0;

    std::vector<mk::Point> polygonPoints{{min, min}, {max, min}, {max, max}, {min, max}, {min, min}};

    mk::Polygon polygon(polygonPoints, mk::Projection::cartesian);

    int numberOfSteps = 100;
    double delta = ((max + extension) - (min - extension)) / static_cast<double>(numberOfSteps - 1);

    double y = min - extension;

    for (int i = 0; i < numberOfSteps; ++i)
    {
        double x = min - extension;

        for (int j = 0; j < numberOfSteps; ++j)
        {
            mk::Point p(x, y);

            // Any point that lies outside of the interval (min..max) for x- or the y-dimension
            // is outside the polygon.
            if (x < min || y < min || x > max || y > max)
            {
                EXPECT_FALSE(polygon.Contains(p));
            }
            else
            {
                EXPECT_TRUE(polygon.Contains(p));
            }

            x += delta;
        }

        y += delta;
    }

    // Check that contains for an invlid point throws an exception.
    mk::Point invalidPoint(mk::constants::missing::doubleValue, mk::constants::missing::doubleValue);
    EXPECT_THROW(polygon.Contains(invalidPoint), mk::ConstraintError);
}

TEST(PolygonTests, AreaCentreAndDirectionTest)
{
    std::vector<mk::Point> polygonPointsLargerSquare{{-10.0, -5.0}, {20.0, -5.0}, {20.0, 10.0}, {-10.0, 10.0}, {-10.0, -5.0}};
    mk::Polygon polygon(polygonPointsLargerSquare, mk::Projection::cartesian);

    double area = 0.0;
    mk::Point centre;
    mk::TraversalDirection direction;

    std::tie(area, centre, direction) = polygon.FaceAreaAndCenterOfMass();

    EXPECT_EQ(area, 450.0);
    EXPECT_EQ(centre.x, 5.0);
    EXPECT_EQ(centre.y, 2.5);
    EXPECT_EQ(direction, mk::TraversalDirection::AntiClockwise);

    // Check the polygon perimeter
    double perimeter = polygon.PerimeterLength();
    EXPECT_EQ(perimeter, 90.0);
}

TEST(PolygonTests, FailureConstructionTests)
{
    std::vector<mk::Point> openPolygon{{0.0, 0.0}, {10.0, 0.0}, {10.0, 10.0}, {0.0, 10.0}};
    EXPECT_THROW([[maybe_unused]] mk::Polygon polygon(openPolygon, mk::Projection::cartesian), mk::ConstraintError);

    std::vector<mk::Point> twoPoints{{0.0, 0.0}, {10.0, 0.0}};
    EXPECT_THROW([[maybe_unused]] mk::Polygon polygon(twoPoints, mk::Projection::cartesian), mk::ConstraintError);

    std::vector<mk::Point> closedPolygonWithNull{{-10.0, -5.0}, {20.0, -5.0}, {20.0, 10.0}, {-10.0, 10.0}, {mk::constants::missing::doubleValue, mk::constants::missing::doubleValue}, {-10.0, -5.0}};
    EXPECT_THROW([[maybe_unused]] mk::Polygon polygon(closedPolygonWithNull, mk::Projection::cartesian), mk::ConstraintError);
}

namespace
{
    void CheckPolygonPointVectors(const std::vector<mk::Point>& actual, const std::vector<mk::Point>& expected)
    {
        ASSERT_EQ(actual.size(), expected.size());

        for (size_t i = 0; i < actual.size(); ++i)
        {
            constexpr double tolerance = 1e-6;
            EXPECT_NEAR(actual[i].x, expected[i].x, tolerance);
            EXPECT_NEAR(actual[i].y, expected[i].y, tolerance);
        }
    }
} // namespace

TEST(PolygonTests, Refine_WhenStartIndexLessThanEndIndex_ThenRefinesBetweenStartIndexAndEndIndex)
{
    // setup
    const std::vector<mk::Point> outer{{0., 0.}, {5., 0.}, {5., 5.}, {0, 5.}, {0., 0.}};
    const std::vector<mk::Point> expected{{0, 0}, {5, 0}, {5, 1.66666666666667}, {5, 3.33333333333333}, {5, 5}, {3.33333333333333, 5}, {1.66666666666667, 5}, {0, 5}, {0, 0}};

    const mk::Polygon polygon(outer, mk::Projection::cartesian);

    // call
    const auto refined = polygon.Refine(1, 3, 2);

    // assert
    SCOPED_TRACE("Refine_WhenStartIndexLessThanEndIndex_ThenRefinesBetweenStartIndexAndEndIndex");
    CheckPolygonPointVectors(refined, expected);
}

TEST(PolygonTests, Refine_WhenStartIndexLargerThanEndIndex_ThenRefinesFromStartIndexToVectorEndAndFromVectorBeginToEndIndex)
{
    // setup
    const std::vector<mk::Point> outer{{0., 0.}, {5., 0.}, {5., 5.}, {0, 5.}, {0., 0.}};
    const std::vector<mk::Point> expected{{0, 0}, {1.66666666666667, 0}, {3.33333333333333, 0}, {5, 0}, {5, 5}, {0, 5}, {0, 3.33333333333333}, {0, 1.66666666666667}, {0, 0}};

    const mk::Polygon polygon(outer, mk::Projection::cartesian);

    // call
    const auto refined = polygon.Refine(3, 1, 2);

    // assert
    SCOPED_TRACE("Refine_WhenStartIndexLargerThanEndIndex_ThenRefinesFromStartIndexToVectorEndAndFromVectorBeginToEndIndex");
    CheckPolygonPointVectors(refined, expected);
}

TEST(PolygonTests, Refine_WhenStartIndexEqualsEndIndex_ThenNoRefinementIsDone)
{
    // setup
    const std::vector<mk::Point> outer{{0., 0.}, {5., 0.}, {5., 5.}, {0, 5.}, {0., 0.}};
    const mk::Polygon polygon(outer, mk::Projection::cartesian);

    // call
    const auto refined = polygon.Refine(2, 2, .5);

    // assert
    SCOPED_TRACE("Refine_WhenStartIndexEqualsEndIndex_ThenNoRefinementIsDone");
    CheckPolygonPointVectors(refined, outer);
}

TEST(PolygonTests, Refine_WhenLastRefinedSegmentSlightlySmallerThanTolerance_AvoidsTinyRefinedSegment)
{
    // setup
    constexpr double d = 2.5 * (1 - .9 * meshkernel::constants::geometric::refinementTolerance);
    const std::vector<mk::Point> outer{{0., 0.}, {5., 0.}, {5., 5.}, {0, 5.}, {0., 0.}};
    const std::vector<mk::Point> expected{{0, 0}, {5, 0}, {5, 2.5}, {5, 5}, {2.5, 5}, {0, 5}, {0, 0}};
    const mk::Polygon polygon(outer, mk::Projection::cartesian);

    // call
    const auto refined = polygon.Refine(1, 3, d);

    // assert
    SCOPED_TRACE("Refine_WhenLastRefinedSegmentSlightlySmallerThanTolerance_AvoidsTinyRefinedSegment");
    CheckPolygonPointVectors(refined, expected);
}

TEST(PolygonTests, Refine_WhenLastRefinementSlightlyLargerThanTolerance_DoesNotOvershootAndUsesOriginalCornerPoint)
{
    // setup
    constexpr double d = 2.5 * (1 + .9 * meshkernel::constants::geometric::refinementTolerance);
    const std::vector<mk::Point> outer{{0., 0.}, {5., 0.}, {5., 5.}, {0, 5.}, {0., 0.}};
    const std::vector<mk::Point> expected{{0, 0}, {5, 0}, {5, 2.5}, {5, 5}, {2.5, 5}, {0, 5}, {0, 0}};

    const mk::Polygon polygon(outer, mk::Projection::cartesian);

    // call
    const auto refined = polygon.Refine(1, 3, d);

    // assert
    SCOPED_TRACE("Refine_WhenLastRefinementSlightlyLargerThanTolerance_DoesNotOvershootAndUsesOriginalCornerPoint");
    CheckPolygonPointVectors(refined, expected);
}

TEST(PolygonTests, Refine_AcceptsRefinedSegmentsLargerThanTheRefinementTolerance)
{
    // setup
    constexpr double d = 2.5 * (1 - 1.1 * meshkernel::constants::geometric::refinementTolerance);
    const std::vector<mk::Point> outer{{0., 0.}, {5., 0.}, {5., 5.}, {0, 5.}, {0., 0.}};
    const std::vector<mk::Point> expected{{0, 0}, {5, 0}, {5, 2.5}, {5, 5}, {2.5, 5}, {0, 5}, {0, 0}};
    const mk::Polygon polygon(outer, mk::Projection::cartesian);

    // call
    const auto refined = polygon.Refine(1, 3, d);

    // assert
    SCOPED_TRACE("Refine_AcceptsRefinedSegmentsLargerThanTheRefinementTolerance");
    CheckPolygonPointVectors(refined, expected);
}

TEST(PolygonTests, LinearRefine_WithEndIndexLargerThanStartIndex_ShouldRefine)
{
    // setup
    const std::vector<mk::Point> outer{
        {1.0704333E+03, 2.3986833E+03},
        {9.6556667E+02, 2.1123167E+03},
        {1.0139667E+03, 1.7533500E+03},
        {1.0946333E+03, 1.5960500E+03},
        {1.1632000E+03, 1.5153833E+03},
        {1.2075667E+03, 1.4710167E+03},
        {1.2801667E+03, 1.4266500E+03},
        {1.7722333E+03, 1.2370833E+03},
        {2.5224333E+03, 1.0233167E+03},
        {3.3492667E+03, 1.0071833E+03},
        {4.0591333E+03, 2.0921500E+03},
        {4.0470333E+03, 2.5196833E+03},
        {2.6394000E+03, 2.8665500E+03},
        {1.6068667E+03, 2.8867167E+03},
        {1.0704333E+03, 2.3986833E+03},
    };

    const mk::Polygon polygon(outer, mk::Projection::cartesian);

    // execute
    std::vector<mk::Point> refinedPolygon = polygon.LinearRefine(0, 8);

    // assert
    constexpr double tolerance = 1.0e-10;
    EXPECT_NEAR(refinedPolygon[0].x, 1070.4332999999999, tolerance);
    EXPECT_NEAR(refinedPolygon[1].x, 980.48220259045979, tolerance);
    EXPECT_NEAR(refinedPolygon[2].x, 1000.6474442197505, tolerance);
    EXPECT_NEAR(refinedPolygon[3].x, 1147.2820369738606, tolerance);
    EXPECT_NEAR(refinedPolygon[4].x, 1509.6856806104856, tolerance);
    EXPECT_NEAR(refinedPolygon[5].x, 1972.1756667754244, tolerance);
    EXPECT_NEAR(refinedPolygon[6].x, 2522.4333000000001, tolerance);
    EXPECT_NEAR(refinedPolygon[7].x, 3349.2667000000001, tolerance);
    EXPECT_NEAR(refinedPolygon[8].x, 4059.1333000000000, tolerance);
    EXPECT_NEAR(refinedPolygon[9].x, 4047.0333000000001, tolerance);
    EXPECT_NEAR(refinedPolygon[10].x, 2639.4000000000001, tolerance);
    EXPECT_NEAR(refinedPolygon[11].x, 1606.8667000000000, tolerance);
    EXPECT_NEAR(refinedPolygon[12].x, 1070.4332999999999, tolerance);

    EXPECT_NEAR(refinedPolygon[0].y, 2398.6833000000001, tolerance);
    EXPECT_NEAR(refinedPolygon[1].y, 2153.0475822179096, tolerance);
    EXPECT_NEAR(refinedPolygon[2].y, 1852.1344283132075, tolerance);
    EXPECT_NEAR(refinedPolygon[3].y, 1534.1103139592642, tolerance);
    EXPECT_NEAR(refinedPolygon[4].y, 1338.2287258438314, tolerance);
    EXPECT_NEAR(refinedPolygon[5].y, 1180.1104928265324, tolerance);
    EXPECT_NEAR(refinedPolygon[6].y, 1023.3167000000000, tolerance);
    EXPECT_NEAR(refinedPolygon[7].y, 1007.1833000000000, tolerance);
    EXPECT_NEAR(refinedPolygon[8].y, 2092.1500000000001, tolerance);
    EXPECT_NEAR(refinedPolygon[9].y, 2519.6833000000001, tolerance);
    EXPECT_NEAR(refinedPolygon[10].y, 2866.5500000000002, tolerance);
    EXPECT_NEAR(refinedPolygon[11].y, 2886.7166999999999, tolerance);
    EXPECT_NEAR(refinedPolygon[12].y, 2398.6833000000001, tolerance);
}

TEST(PolygonTests, LinearRefine_WithStartIndexLargerThanEndIndex_ShouldRefine)
{
    // setup
    const std::vector<mk::Point> outer{
        {1.0704333E+03, 2.3986833E+03},
        {9.6556667E+02, 2.1123167E+03},
        {1.0139667E+03, 1.7533500E+03},
        {1.0946333E+03, 1.5960500E+03},
        {1.1632000E+03, 1.5153833E+03},
        {1.2075667E+03, 1.4710167E+03},
        {1.2801667E+03, 1.4266500E+03},
        {1.7722333E+03, 1.2370833E+03},
        {2.5224333E+03, 1.0233167E+03},
        {3.3492667E+03, 1.0071833E+03},
        {4.0591333E+03, 2.0921500E+03},
        {4.0470333E+03, 2.5196833E+03},
        {2.6394000E+03, 2.8665500E+03},
        {1.6068667E+03, 2.8867167E+03},
        {1.0704333E+03, 2.3986833E+03},
    };

    const mk::Polygon polygon(outer, mk::Projection::cartesian);

    // execute
    std::vector<mk::Point> refinedPolygon = polygon.LinearRefine(8, 0);

    // assert
    constexpr double tolerance = 1.0e-10;
    EXPECT_NEAR(refinedPolygon[0].x, 1070.4332999999999, tolerance);
    EXPECT_NEAR(refinedPolygon[1].x, 965.56667000000004, tolerance);
    EXPECT_NEAR(refinedPolygon[2].x, 1013.9666999999999, tolerance);
    EXPECT_NEAR(refinedPolygon[3].x, 1094.6333000000000, tolerance);
    EXPECT_NEAR(refinedPolygon[4].x, 1163.2000000000000, tolerance);
    EXPECT_NEAR(refinedPolygon[5].x, 1207.5667000000001, tolerance);
    EXPECT_NEAR(refinedPolygon[6].x, 1280.1667000000000, tolerance);
    EXPECT_NEAR(refinedPolygon[7].x, 1772.2333000000001, tolerance);
    EXPECT_NEAR(refinedPolygon[8].x, 2522.4333000000001, tolerance);
    EXPECT_NEAR(refinedPolygon[9].x, 3372.6635744129117, tolerance);
    EXPECT_NEAR(refinedPolygon[10].x, 3839.9884165930689, tolerance);
    EXPECT_NEAR(refinedPolygon[11].x, 4037.5883063359329, tolerance);
    EXPECT_NEAR(refinedPolygon[12].x, 3239.3399285010146, tolerance);
    EXPECT_NEAR(refinedPolygon[13].x, 2450.4753508170843, tolerance);
    EXPECT_NEAR(refinedPolygon[14].x, 1658.7752716626733, tolerance);
    EXPECT_NEAR(refinedPolygon[15].x, 1070.4332999999999, tolerance);

    EXPECT_NEAR(refinedPolygon[0].y, 2398.6833000000001, tolerance);
    EXPECT_NEAR(refinedPolygon[1].y, 2112.3166999999999, tolerance);
    EXPECT_NEAR(refinedPolygon[2].y, 1753.3499999999999, tolerance);
    EXPECT_NEAR(refinedPolygon[3].y, 1596.0500000000000, tolerance);
    EXPECT_NEAR(refinedPolygon[4].y, 1515.3833000000000, tolerance);
    EXPECT_NEAR(refinedPolygon[5].y, 1471.0166999999999, tolerance);
    EXPECT_NEAR(refinedPolygon[6].y, 1426.6500000000001, tolerance);
    EXPECT_NEAR(refinedPolygon[7].y, 1237.0833000000000, tolerance);
    EXPECT_NEAR(refinedPolygon[8].y, 1023.3167000000000, tolerance);
    EXPECT_NEAR(refinedPolygon[9].y, 1042.9433000085805, tolerance);
    EXPECT_NEAR(refinedPolygon[10].y, 1757.2069262282478, tolerance);
    EXPECT_NEAR(refinedPolygon[11].y, 2522.0107199209242, tolerance);
    EXPECT_NEAR(refinedPolygon[12].y, 2718.7137821459733, tolerance);
    EXPECT_NEAR(refinedPolygon[13].y, 2870.2399407725420, tolerance);
    EXPECT_NEAR(refinedPolygon[14].y, 2885.7028590012160, tolerance);
    EXPECT_NEAR(refinedPolygon[15].y, 2398.6833000000001, tolerance);
}
