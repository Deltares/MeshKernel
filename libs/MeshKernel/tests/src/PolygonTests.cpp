#include <gtest/gtest.h>

#include "MeshKernel/BoundingBox.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/LandBoundary.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/Point.hpp"
#include "MeshKernel/Polygon.hpp"

#include "TestUtils/Definitions.hpp"
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
