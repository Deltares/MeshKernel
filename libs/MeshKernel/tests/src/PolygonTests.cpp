//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2025.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation version 3.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// contact: delft3d.support@deltares.nl
// Stichting Deltares
// P.O. Box 177
// 2600 MH Delft, The Netherlands
//
// All indications and logos of, and references to, "Delft3D" and "Deltares"
// are registered trademarks of Stichting Deltares, and remain the property of
// Stichting Deltares. All rights reserved.
//
//------------------------------------------------------------------------------

#include <iomanip> // TODO remove this

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
    const std::vector<mk::Point> expected{{0, 0}, {5, 0}, {5, 2.0}, {5, 4.0}, {4.0, 5}, {2.0, 5}, {0, 5}, {0, 0}};
    // const std::vector<mk::Point> expected{{0, 0}, {5, 0}, {5, 1.66666666666667}, {5, 3.33333333333333}, {5, 5}, {3.33333333333333, 5}, {1.66666666666667, 5}, {0, 5}, {0, 0}};

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
    const std::vector<mk::Point> expected{{5.0, 0}, {5.0, 5.0}, {0, 5.0}, {0, 3.0}, {0, 1.0}, {1.0, 0.0}, {3.0, 0.0}, {5.0, 0.0}};
    // const std::vector<mk::Point> expected{{0, 0}, {1.66666666666667, 0}, {3.33333333333333, 0}, {5, 0}, {5, 5}, {0, 5}, {0, 3.33333333333333}, {0, 1.66666666666667}, {0, 0}};

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
    const std::vector<mk::Point> expected{{0, 0}, {5, 0}, {5, 2}, {5, 4}, {4, 5}, {2, 5}, {0, 5}, {0, 0}};
    // const std::vector<mk::Point> expected{{0, 0}, {5, 0}, {5, 2.5}, {5, 5}, {2.5, 5}, {0, 5}, {0, 0}};
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
    const std::vector<mk::Point> expected{{0, 0}, {5, 0}, {5, 2}, {5, 4}, {4, 5}, {2, 5}, {0, 5}, {0, 0}};
    // const std::vector<mk::Point> expected{{0, 0}, {5, 0}, {5, 2.5}, {5, 5}, {2.5, 5}, {0, 5}, {0, 0}};
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

TEST(PolygonTests, ResampleSmallerThenResampleLarger)
{

    const double tolerance = 1.0e-8;

    // setup
    const std::vector<mk::Point> outer{
        {0.0, 0.0},
        {10.0, 0.0},
        {10.0, 10.0},
        {0.0, 10.0},
        {0.0, 0.0}};

    std::cout.precision(12);
    const mk::Polygon polygon(outer, mk::Projection::cartesian);

    // execute
    std::vector<mk::Point> refinedPolygon = polygon.Refine(0, 4, 2.5);
    const std::vector<mk::Point> expectedFirstRefinement{{0.0, 0.0},
                                                         {2.5, 0.0},
                                                         {5.0, 0.0},
                                                         {7.5, 0.0},
                                                         {10.0, 0.0},
                                                         {10.0, 2.5},
                                                         {10.0, 5.0},
                                                         {10.0, 7.5},
                                                         {10.0, 10.0},
                                                         {7.5, 10.0},
                                                         {5.0, 10.0},
                                                         {2.5, 10.0},
                                                         {0.0, 10.0},
                                                         {0.0, 7.5},
                                                         {0.0, 5.0},
                                                         {0.0, 2.5},
                                                         {0.0, 0.0}};

    ASSERT_EQ(expectedFirstRefinement.size(), refinedPolygon.size());

    for (size_t i = 0; i < refinedPolygon.size(); ++i)
    {
        EXPECT_NEAR(expectedFirstRefinement[i].x, refinedPolygon[i].x, tolerance);
        EXPECT_NEAR(expectedFirstRefinement[i].y, refinedPolygon[i].y, tolerance);
    }

    const mk::Polygon polygon2(refinedPolygon, mk::Projection::cartesian);

    // Resample polygon, with start index larger than end index and resample distance larger than first resampling.
    std::vector<mk::Point> refinedPolygon2 = polygon2.Refine(12, 4, 7.0);

    const std::vector<mk::Point> expectedSecondRefinement{{10.0, 0.0},
                                                          {10.0, 2.5},
                                                          {10.0, 5.0},
                                                          {10.0, 7.5},
                                                          {10.0, 10.0},
                                                          {7.5, 10.0},
                                                          {5.0, 10.0},
                                                          {2.5, 10.0},
                                                          {0.0, 10.0},
                                                          {0.0, 3.33333333333},
                                                          {3.33333333333, 0.0},
                                                          {10.0, 0.0}};

    ASSERT_EQ(expectedSecondRefinement.size(), refinedPolygon2.size());

    for (size_t i = 0; i < refinedPolygon2.size(); ++i)
    {
        EXPECT_NEAR(expectedSecondRefinement[i].x, refinedPolygon2[i].x, tolerance);
        EXPECT_NEAR(expectedSecondRefinement[i].y, refinedPolygon2[i].y, tolerance);
    }
}

TEST(PolygonTests, LinearRefine_WithVaryingSegmentSize_WithStartIndexLargerThanEndIndex_ShouldRefine)
{

    std::vector<meshkernel::Point> outer{{-27.9394258397185, 1665.31356424612},
                                         {76.9523523071914, 1668.31047219317},
                                         {178.847222507047, 1668.31047219317},
                                         {259.763737077521, 1668.31047219317},
                                         {388.630778800867, 1668.31047219317},
                                         {505.510188735996, 1665.31356424612},
                                         {580.43288741236, 1662.31665629906},
                                         {769.238088076798, 1653.3259324579},
                                         {850.154602647272, 1656.32284040495},
                                         {967.0340125824, 1656.32284040495},
                                         {1110.88559404102, 1656.32284040495},
                                         {1281.70934702313, 1653.3259324579},
                                         {1365.62276954066, 1653.3259324579},
                                         {1449.53619205819, 1656.32284040495},
                                         {1614.36612914619, 1653.3259324579},
                                         {1734.24244702837, 1653.3259324579},
                                         {1890.08166027521, 1653.3259324579},
                                         {2048.9177814691, 1650.32902451084},
                                         {2078.88686093965, 556.457623835923},
                                         {-66.899229151428, 589.423611253523},
                                         {-27.9394258397185, 1665.31356424612}};

    std::vector<meshkernel::Point> expected{{259.763737077521, 1668.31047219317},
                                            {388.630778800867, 1668.31047219317},
                                            {505.510188735996, 1665.31356424612},
                                            {580.43288741236, 1662.31665629906},
                                            {769.238088076798, 1653.3259324579},
                                            {850.154602647272, 1656.32284040495},
                                            {967.0340125824, 1656.32284040495},
                                            {1110.88559404102, 1656.32284040495},
                                            {1281.70934702313, 1653.3259324579},
                                            {1365.62276954066, 1653.3259324579},
                                            {1449.53619205819, 1656.32284040495},
                                            {1614.36612914619, 1653.3259324579},
                                            {1734.24244702837, 1653.3259324579},
                                            {1889.0205150208, 1653.3259324579},
                                            {2041.43040301866, 1650.47029580236},
                                            {2052.82432345575, 1507.74024199828},
                                            {2056.87376421424, 1359.93565431344},
                                            {2060.86195131697, 1214.3668250635},
                                            {2064.78981131421, 1070.99993516439},
                                            {2068.65825674078, 929.801677094647},
                                            {2072.4681863281, 790.739247157293},
                                            {2076.220485213, 653.78033785879},
                                            {2041.31270598354, 557.034880406478},
                                            {1908.43169489973, 559.076348453855},
                                            {1777.56070154287, 561.086936340622},
                                            {1648.66932148398, 563.067111173929},
                                            {1521.72761020512, 565.017332995252},
                                            {1396.70607614254, 566.938054887275},
                                            {1273.57567383509, 568.829723079149},
                                            {1152.30779717631, 570.692777050164},
                                            {1032.8742727685, 572.527649631848},
                                            {915.247353377413, 574.334767108527},
                                            {799.399711485897, 576.114549316357},
                                            {685.304432945077, 577.867409740867},
                                            {572.935010721554, 579.593755613016},
                                            {462.265338739193, 581.293988003807},
                                            {353.269705814071, 582.968501917461},
                                            {245.922789681153, 584.61768638319},
                                            {140.199651111336, 586.241924545576},
                                            {36.0757281174708, 587.841593753582},
                                            {-66.4731697519648, 589.417065648224},
                                            {-63.2593151982587, 689.941235037199},
                                            {-59.6592734167226, 789.357773465774},
                                            {-56.1136874806049, 887.270492778563},
                                            {-52.6217336663824, 983.702140417477},
                                            {-49.1826007105429, 1078.67511973643},
                                            {-45.7954896211099, 1172.21149520616},
                                            {-42.4596134920174, 1264.33299754033},
                                            {-39.1741973202934, 1355.06102874409},
                                            {-35.9384778260081, 1444.41666708628},
                                            {-32.7517032749458, 1532.42067199638},
                                            {-29.6131333039599, 1619.09348888745},
                                            {11.2119415203001, 1666.43217474212},
                                            {95.310761132729, 1668.31047219317},
                                            {178.163885853533, 1668.31047219317},
                                            {259.763737077521, 1668.31047219317}};

    const mk::Polygon polygon(outer, mk::Projection::cartesian);

    // execute
    std::vector<mk::Point> refinedPolygon = polygon.LinearRefine(15, 3);

    std::cout.precision(15);

    ASSERT_EQ(refinedPolygon.size(), expected.size());

    const double tolerance = 1.0e-10;

    for (size_t i = 0; i < refinedPolygon.size(); ++i)
    {
        EXPECT_NEAR(expected[i].x, refinedPolygon[i].x, tolerance);
        EXPECT_NEAR(expected[i].y, refinedPolygon[i].y, tolerance);
    }
}
