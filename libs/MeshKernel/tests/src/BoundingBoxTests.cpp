#include <gtest/gtest.h>

#include "MeshKernel/BoundingBox.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Utilities/RTree.hpp"

TEST(BoundingBox, Contains_WhenPointInside_MustReturnTrue)
{
    // Setup
    const auto boundingBox = meshkernel::BoundingBox({0.0, 0.0}, {10.0, 10.0});

    std::vector<meshkernel::Point> points;
    points.push_back({0.0, 0.0});
    points.push_back({5.5, 5.0});
    points.push_back({10.0, 10.0});
    points.push_back({10.1, 10.1});

    // Assert
    ASSERT_TRUE(boundingBox.Contains(points[0]));
    ASSERT_TRUE(boundingBox.Contains(points[1]));
    ASSERT_TRUE(boundingBox.Contains(points[2]));
    ASSERT_FALSE(boundingBox.Contains(points[3]));
}

TEST(BoundingBox, Contains_WhenPointOutside_MustReturnFalse)
{
    // 1 Setup
    const auto boundingBox = meshkernel::BoundingBox({0.0, 0.0}, {10.0, 10.0});

    std::vector<meshkernel::Point> points;
    points.push_back({-10.0, -10.0});
    points.push_back({20.0, 20.0});

    // Assert
    ASSERT_FALSE(boundingBox.Contains(points[0]));
    ASSERT_FALSE(boundingBox.Contains(points[1]));
}
