#include <gtest/gtest.h>

#include "MeshKernel/BoundingBox.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/LandBoundary.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/Point.hpp"
#include "MeshKernel/PolygonalEnclosure.hpp"

#include "TestUtils/Definitions.hpp"
#include "TestUtils/MakeMeshes.hpp"

namespace mk = meshkernel;

TEST(PolygonalEnclosureTests, BasicConstruction)
{

    std::vector<mk::Point> outerPolygon{{0.0, 0.0}, {100.0, 0.0}, {100.0, 100.0}, {0.0, 100.0}, {0.0, 0.0}};
    std::vector<mk::Point> innerPolygon{{40.0, 40.0}, {40.0, 60.0}, {60.0, 60.0}, {40.0, 60.0}, {40.0, 40.0}};

    std::vector<mk::Point> polygonPoints;

    polygonPoints.insert(polygonPoints.end(), outerPolygon.begin(), outerPolygon.end());
    polygonPoints.push_back({mk::constants::missing::innerOuterSeparator, mk::constants::missing::innerOuterSeparator});
    polygonPoints.insert(polygonPoints.end(), innerPolygon.begin(), innerPolygon.end());

    // Start test.

    mk::PolygonalEnclosure enclosure(polygonPoints, mk::Projection::cartesian);

    ASSERT_EQ(enclosure.Outer().Size(), outerPolygon.size());

    for (size_t i = 0; i < outerPolygon.size(); ++i)
    {
        EXPECT_TRUE(outerPolygon[i] == enclosure.Outer().Node(i))
            << "Expected points to be same, expected: " << outerPolygon[i].x << ", " << outerPolygon[i].y
            << ", actual: " << enclosure.Outer().Node(i).x << ", " << enclosure.Outer().Node(i).y;
    }

    ASSERT_EQ(enclosure.NumberOfInner(), 1);
    ASSERT_EQ(enclosure.Inner(0).Size(), innerPolygon.size());

    for (size_t i = 0; i < innerPolygon.size(); ++i)
    {
        EXPECT_TRUE(innerPolygon[i] == enclosure.Inner(0).Node(i))
            << "Expected points to be same, expected: " << innerPolygon[i].x << ", " << innerPolygon[i].y
            << ", actual: " << enclosure.Inner(0).Node(i).x << ", " << enclosure.Inner(0).Node(i).y;
    }
}

TEST(PolygonalEnclosureTests, ContainsTest)
{

    constexpr double outerMin = 0.0;
    constexpr double outerMax = 100.0;

    constexpr double innerMin = 40.0;
    constexpr double innerMax = 60.0;

    constexpr double extension = 2.0;

    std::vector<mk::Point> outerPolygon{{outerMin, outerMin}, {outerMax, outerMin}, {outerMax, outerMax}, {outerMin, outerMax}, {outerMin, outerMin}};
    // Inner polygon needs to be reversed.
    std::vector<mk::Point> innerPolygon{{innerMin, innerMin}, {innerMin, innerMax}, {innerMax, innerMax}, {innerMax, innerMin}, {innerMin, innerMin}};

    std::vector<mk::Point> polygonPoints;

    polygonPoints.insert(polygonPoints.end(), outerPolygon.begin(), outerPolygon.end());
    polygonPoints.push_back({mk::constants::missing::innerOuterSeparator, mk::constants::missing::innerOuterSeparator});
    polygonPoints.insert(polygonPoints.end(), innerPolygon.begin(), innerPolygon.end());

    // Start test.

    mk::PolygonalEnclosure enclosure(polygonPoints, mk::Projection::cartesian);
    ASSERT_EQ(enclosure.NumberOfInner(), 1);

    int numberOfSteps = 100;
    double delta = ((outerMax + extension) - (outerMin - extension)) / static_cast<double>(numberOfSteps - 1);

    double y = outerMin - extension;

    for (int i = 0; i < numberOfSteps; ++i)
    {
        double x = outerMin - extension;
        bool outsideOuterY = (y < outerMin) || (y > outerMax);
        bool insideInnerY = y > innerMin && y < innerMax;

        for (int j = 0; j < numberOfSteps; ++j)
        {
            mk::Point p(x, y);
            bool outsideOuterX = (x < outerMin) || (x > outerMax);
            bool insideInnerX = (x > innerMin) && (x < innerMax);

            if (outsideOuterX || outsideOuterY || (insideInnerX && insideInnerY))
            {
                // If point lies outside the outer polygon or inside the inner polygon then it is considered to be outside the enclosure.
                EXPECT_FALSE(enclosure.Contains(p));
            }
            else
            {
                EXPECT_TRUE(enclosure.Contains(p));
            }

            x += delta;
        }

        y += delta;
    }
}

TEST(PolygonalEnclosureTests, MultipleInnerPolygonsContainsTest)
{

    constexpr double outerMin = 0.0;
    constexpr double outerMax = 100.0;

    constexpr double innerMin1 = 20.0;
    constexpr double innerMax1 = 40.0;
    constexpr double innerMin2 = 60.0;
    constexpr double innerMax2 = 80.0;

    constexpr double extension = 2.0;

    std::vector<mk::Point> outerPolygon{{outerMin, outerMin}, {outerMax, outerMin}, {outerMax, outerMax}, {outerMin, outerMax}, {outerMin, outerMin}};
    // Inner polygon needs to be reversed.
    std::vector<mk::Point> innerPolygon1{{innerMin1, innerMin1}, {innerMin1, innerMax1}, {innerMax1, innerMax1}, {innerMax1, innerMin1}, {innerMin1, innerMin1}};
    std::vector<mk::Point> innerPolygon2{{innerMin2, innerMin2}, {innerMin2, innerMax2}, {innerMax2, innerMax2}, {innerMax2, innerMin2}, {innerMin2, innerMin2}};

    std::vector<mk::Point> polygonPoints;

    polygonPoints.insert(polygonPoints.end(), outerPolygon.begin(), outerPolygon.end());
    polygonPoints.push_back({mk::constants::missing::innerOuterSeparator, mk::constants::missing::innerOuterSeparator});
    polygonPoints.insert(polygonPoints.end(), innerPolygon1.begin(), innerPolygon1.end());
    polygonPoints.push_back({mk::constants::missing::innerOuterSeparator, mk::constants::missing::innerOuterSeparator});
    polygonPoints.insert(polygonPoints.end(), innerPolygon2.begin(), innerPolygon2.end());

    // Start test.

    mk::PolygonalEnclosure enclosure(polygonPoints, mk::Projection::cartesian);

    ASSERT_EQ(enclosure.NumberOfInner(), 2);

    int numberOfSteps = 100;
    double delta = ((outerMax + extension) - (outerMin - extension)) / static_cast<double>(numberOfSteps - 1);

    double y = outerMin - extension;

    for (int i = 0; i < numberOfSteps; ++i)
    {
        double x = outerMin - extension;
        bool outsideOuterY = (y < outerMin) || (y > outerMax);
        bool insideInner1Y = y > innerMin1 && y < innerMax1;
        bool insideInner2Y = y > innerMin2 && y < innerMax2;

        for (int j = 0; j < numberOfSteps; ++j)
        {
            mk::Point p(x, y);
            bool outsideOuterX = (x < outerMin) || (x > outerMax);
            bool insideInner1X = (x > innerMin1) && (x < innerMax1);
            bool insideInner2X = (x > innerMin2) && (x < innerMax2);

            // If point lies outside the outer polygon or inside the inner polygons then it is considered to be outside the enclosure.
            if (outsideOuterX || outsideOuterY || (insideInner1X && insideInner1Y) || (insideInner2X && insideInner2Y))
            {
                EXPECT_FALSE(enclosure.Contains(p));
            }
            else
            {
                // Only if the point lies inside the outer perimeter polygon and outside any inner polygons is it considered to be inside the enclosure.
                EXPECT_TRUE(enclosure.Contains(p));
            }

            x += delta;
        }

        y += delta;
    }
}
