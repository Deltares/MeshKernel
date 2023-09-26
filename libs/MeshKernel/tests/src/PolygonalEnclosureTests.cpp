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

    // Combine the points into a single array
    polygonPoints.insert(polygonPoints.end(), outerPolygon.begin(), outerPolygon.end());
    polygonPoints.push_back({mk::constants::missing::innerOuterSeparator, mk::constants::missing::innerOuterSeparator});
    polygonPoints.insert(polygonPoints.end(), innerPolygon.begin(), innerPolygon.end());

    //--------------------------------
    // Start test.

    mk::PolygonalEnclosure enclosure(polygonPoints, mk::Projection::Cartesian);

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
    std::vector<mk::Point> innerPolygon{{innerMin, innerMin}, {innerMax, innerMin}, {innerMax, innerMax}, {innerMin, innerMax}, {innerMin, innerMin}};

    std::vector<mk::Point> polygonPoints;

    // Combine the points into a single array
    polygonPoints.insert(polygonPoints.end(), outerPolygon.begin(), outerPolygon.end());
    polygonPoints.push_back({mk::constants::missing::innerOuterSeparator, mk::constants::missing::innerOuterSeparator});
    polygonPoints.insert(polygonPoints.end(), innerPolygon.begin(), innerPolygon.end());

    //--------------------------------
    // Start test.

    mk::PolygonalEnclosure enclosure(polygonPoints, mk::Projection::Cartesian);
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
    std::vector<mk::Point> innerPolygon1{{innerMin1, innerMin1}, {innerMax1, innerMin1}, {innerMax1, innerMax1}, {innerMin1, innerMax1}, {innerMin1, innerMin1}};
    std::vector<mk::Point> innerPolygon2{{innerMin2, innerMin2}, {innerMax2, innerMin2}, {innerMax2, innerMax2}, {innerMin2, innerMax2}, {innerMin2, innerMin2}};

    std::vector<mk::Point> polygonPoints;

    // Combine the points into a single array
    polygonPoints.insert(polygonPoints.end(), outerPolygon.begin(), outerPolygon.end());
    polygonPoints.push_back({mk::constants::missing::innerOuterSeparator, mk::constants::missing::innerOuterSeparator});
    polygonPoints.insert(polygonPoints.end(), innerPolygon1.begin(), innerPolygon1.end());
    polygonPoints.push_back({mk::constants::missing::innerOuterSeparator, mk::constants::missing::innerOuterSeparator});
    polygonPoints.insert(polygonPoints.end(), innerPolygon2.begin(), innerPolygon2.end());

    //--------------------------------
    // Start test.

    mk::PolygonalEnclosure enclosure(polygonPoints, mk::Projection::Cartesian);

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

TEST(PolygonalEnclosureTests, RefineTest)
{

    constexpr double outerMin = 0.0;
    constexpr double outerMax = 100.0;

    constexpr double innerMin = 40.0;
    constexpr double innerMax = 60.0;

    std::vector<mk::Point> outerPolygon{{outerMin, outerMin}, {outerMax, outerMin}, {outerMax, outerMax}, {outerMin, outerMax}, {outerMin, outerMin}};
    std::vector<mk::Point> innerPolygon{{innerMin, innerMin}, {innerMax, innerMin}, {innerMax, innerMax}, {innerMin, innerMax}, {innerMin, innerMin}};

    std::vector<mk::Point> polygonPoints;

    // Combine the points into a single array
    polygonPoints.insert(polygonPoints.end(), outerPolygon.begin(), outerPolygon.end());
    polygonPoints.push_back({mk::constants::missing::innerOuterSeparator, mk::constants::missing::innerOuterSeparator});
    polygonPoints.insert(polygonPoints.end(), innerPolygon.begin(), innerPolygon.end());

    //--------------------------------
    // Start test.

    mk::PolygonalEnclosure enclosure(polygonPoints, mk::Projection::Cartesian);

    double delta = 25.0;

    std::vector<mk::Point> refinedPoints = enclosure.Refine(0, 1, delta);
    std::vector<mk::Point> expectedPoints{{outerMin, outerMin}, {outerMin + delta, outerMin}, {outerMin + 2 * delta, outerMin}, {outerMin + 3 * delta, outerMin}, {outerMax, outerMin}, {outerMax, outerMax}, {outerMin, outerMax}, {outerMin, outerMin}};

    ASSERT_EQ(refinedPoints.size(), expectedPoints.size());
    constexpr double tolerance = 1.0e-10;

    for (size_t i = 0; i < refinedPoints.size(); ++i)
    {
        EXPECT_NEAR(refinedPoints[i].x, expectedPoints[i].x, tolerance);
        EXPECT_NEAR(refinedPoints[i].y, expectedPoints[i].y, tolerance);
    }
}

TEST(PolygonalEnclosureTests, RefineFailureTest)
{

    constexpr double outerMin = 0.0;
    constexpr double outerMax = 100.0;

    constexpr double innerMin = 40.0;
    constexpr double innerMax = 60.0;

    std::vector<mk::Point> outerPolygon{{outerMin, outerMin}, {outerMax, outerMin}, {outerMax, outerMax}, {outerMin, outerMax}, {outerMin, outerMin}};
    std::vector<mk::Point> innerPolygon{{innerMin, innerMin}, {innerMax, innerMin}, {innerMax, innerMax}, {innerMin, innerMax}, {innerMin, innerMin}};

    std::vector<mk::Point> polygonPoints;

    // Combine the points into a single array
    polygonPoints.insert(polygonPoints.end(), outerPolygon.begin(), outerPolygon.end());
    polygonPoints.push_back({mk::constants::missing::innerOuterSeparator, mk::constants::missing::innerOuterSeparator});
    polygonPoints.insert(polygonPoints.end(), innerPolygon.begin(), innerPolygon.end());

    mk::PolygonalEnclosure enclosure(polygonPoints, mk::Projection::Cartesian);

    // End point lies in inner polygon, so outside the outer polygon
    EXPECT_THROW([[maybe_unused]] auto refinedPoints = enclosure.Refine(0, 6, 25.0), mk::ConstraintError);
    // All point indices are in inner polygon
    EXPECT_THROW([[maybe_unused]] auto refinedPoints = enclosure.Refine(6, 10, 25.0), mk::ConstraintError);
    // All point indices are in outer polygon, but in wrong order
    EXPECT_THROW([[maybe_unused]] auto refinedPoints = enclosure.Refine(4, 0, 25.0), mk::ConstraintError);
}

TEST(PolygonalEnclosureTests, SnapToLandBoundaryTest)
{
    constexpr double tolerance = 1.0e-10;

    constexpr double outerMin = 0.0;
    constexpr double outerMax = 100.0;

    constexpr double innerMin = 40.0;
    constexpr double innerMax = 60.0;
    constexpr double delta = 10.0;

    std::vector<mk::Point> landBoundaryNodes{{outerMin, outerMax + delta}, {outerMax, outerMax + delta}};
    mk::LandBoundary landBoundary(landBoundaryNodes);

    std::vector<mk::Point> outerPolygon{{outerMin, outerMin}, {outerMax, outerMin}, {outerMax, outerMax}, {outerMin, outerMax}, {outerMin, outerMin}};
    std::vector<mk::Point> innerPolygon{{innerMin, innerMin}, {innerMax, innerMin}, {innerMax, innerMax}, {innerMin, innerMax}, {innerMin, innerMin}};

    std::vector<mk::Point> polygonPoints;

    // Combine the points into a single array
    polygonPoints.insert(polygonPoints.end(), outerPolygon.begin(), outerPolygon.end());
    polygonPoints.push_back({mk::constants::missing::innerOuterSeparator, mk::constants::missing::innerOuterSeparator});
    polygonPoints.insert(polygonPoints.end(), innerPolygon.begin(), innerPolygon.end());

    //--------------------------------
    // Start test.

    mk::PolygonalEnclosure enclosure(polygonPoints, mk::Projection::Cartesian);

    enclosure.SnapToLandBoundary(2, 3, landBoundary);

    ASSERT_EQ(enclosure.Outer().Size(), outerPolygon.size());

    std::vector<mk::Point> expectedOuterPolygon{{outerMin, outerMin}, {outerMax, outerMin}, {outerMax, outerMax + delta}, {outerMin, outerMax + delta}, {outerMin, outerMin}};

    for (size_t i = 0; i < enclosure.Outer().Size(); ++i)
    {
        EXPECT_NEAR(enclosure.Outer().Node(i).x, expectedOuterPolygon[i].x, tolerance);
        EXPECT_NEAR(enclosure.Outer().Node(i).y, expectedOuterPolygon[i].y, tolerance);
    }
}

TEST(PolygonalEnclosureTests, SnappingFailureTest)
{

    constexpr double outerMin = 0.0;
    constexpr double outerMax = 100.0;

    constexpr double innerMin = 40.0;
    constexpr double innerMax = 60.0;
    constexpr double delta = 10.0;

    std::vector<mk::Point> landBoundaryNodes{{outerMin, outerMax + delta}, {outerMax, outerMax + delta}};
    mk::LandBoundary landBoundary(landBoundaryNodes);

    std::vector<mk::Point> outerPolygon{{outerMin, outerMin}, {outerMax, outerMin}, {outerMax, outerMax}, {outerMin, outerMax}, {outerMin, outerMin}};
    std::vector<mk::Point> innerPolygon{{innerMin, innerMin}, {innerMax, innerMin}, {innerMax, innerMax}, {innerMin, innerMax}, {innerMin, innerMin}};

    std::vector<mk::Point> polygonPoints;

    // Combine the points into a single array
    polygonPoints.insert(polygonPoints.end(), outerPolygon.begin(), outerPolygon.end());
    polygonPoints.push_back({mk::constants::missing::innerOuterSeparator, mk::constants::missing::innerOuterSeparator});
    polygonPoints.insert(polygonPoints.end(), innerPolygon.begin(), innerPolygon.end());

    mk::PolygonalEnclosure enclosure(polygonPoints, mk::Projection::Cartesian);

    EXPECT_THROW(enclosure.SnapToLandBoundary(3, 2, landBoundary), mk::ConstraintError);
    EXPECT_THROW(enclosure.SnapToLandBoundary(0, 6, landBoundary), mk::ConstraintError);
}

TEST(PolygonalEnclosureTests, OffsetTest)
{
    // Test the polygon offset.

    constexpr double tolerance = 1.0e-10;

    constexpr double outerMin = 0.0;
    constexpr double outerMax = 100.0;

    constexpr double innerMin = 20.0;
    constexpr double innerMax = 80.0;
    constexpr double distance = 5.0;

    std::vector<mk::Point> outerPolygon{{outerMin, outerMin}, {outerMax, outerMin}, {outerMax, outerMax}, {outerMin, outerMax}, {outerMin, outerMin}};
    std::vector<mk::Point> innerPolygon{{innerMin, innerMin}, {innerMax, innerMin}, {innerMax, innerMax}, {innerMin, innerMax}, {innerMin, innerMin}};

    std::vector<mk::Point> polygonPoints;

    // Combine the points into a single array
    polygonPoints.insert(polygonPoints.end(), outerPolygon.begin(), outerPolygon.end());
    polygonPoints.push_back({mk::constants::missing::innerOuterSeparator, mk::constants::missing::innerOuterSeparator});
    polygonPoints.insert(polygonPoints.end(), innerPolygon.begin(), innerPolygon.end());

    std::vector<mk::Point> expectedOffset{{-distance, -distance}, {distance, -distance}, {distance, distance}, {-distance, distance}, {-distance, -distance}};

    //--------------------------------
    // Start test.

    mk::PolygonalEnclosure enclosure(polygonPoints, mk::Projection::Cartesian);

    auto [outward, inward] = enclosure.OffsetCopy(distance, true);

    // Should not return null pointers, if any pointer is null then abort the test.
    ASSERT_TRUE(outward != nullptr);
    ASSERT_TRUE(inward != nullptr);

    // First check the outer perimeter polygon has been offset correctly
    std::vector<mk::Point> actualPoints = outward->Outer().Nodes();

    for (size_t i = 0; i < actualPoints.size(); ++i)
    {
        EXPECT_NEAR(actualPoints[i].x, outerPolygon[i].x + expectedOffset[i].x, tolerance);
        EXPECT_NEAR(actualPoints[i].y, outerPolygon[i].y + expectedOffset[i].y, tolerance);
    }

    // Now check the inner polygon (there is only 1)
    actualPoints = outward->Inner(0).Nodes();

    for (size_t i = 0; i < actualPoints.size(); ++i)
    {
        EXPECT_NEAR(actualPoints[i].x, innerPolygon[i].x + expectedOffset[i].x, tolerance);
        EXPECT_NEAR(actualPoints[i].y, innerPolygon[i].y + expectedOffset[i].y, tolerance);
    }

    // First check the inner polygon has been offset correctly
    actualPoints = inward->Outer().Nodes();

    for (size_t i = 0; i < actualPoints.size(); ++i)
    {
        EXPECT_NEAR(actualPoints[i].x, outerPolygon[i].x - expectedOffset[i].x, tolerance);
        EXPECT_NEAR(actualPoints[i].y, outerPolygon[i].y - expectedOffset[i].y, tolerance);
    }

    // Now check the inner polygon (there is only 1)
    actualPoints = inward->Inner(0).Nodes();

    for (size_t i = 0; i < actualPoints.size(); ++i)
    {
        EXPECT_NEAR(actualPoints[i].x, innerPolygon[i].x - expectedOffset[i].x, tolerance);
        EXPECT_NEAR(actualPoints[i].y, innerPolygon[i].y - expectedOffset[i].y, tolerance);
    }
}
