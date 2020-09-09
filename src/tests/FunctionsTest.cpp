#include "../Operations.cpp"
#include "../Entities.hpp"
#include <gtest/gtest.h>
#include <chrono>
#include <random>

TEST(FunctionsTest, NormalVectorInsideTestCartesian)
{
    //1 Setup
    MeshKernel::Point firstPoint{0.0,0.0};
    MeshKernel::Point secondPoint{ 1.0,0.0 };
    MeshKernel::Point insidePoint{ 0.0,-1.0 };
    MeshKernel::Point normal;
    bool flippedNormal;

    // 2 Execute
    NormalVectorInside(firstPoint, secondPoint, insidePoint, normal, flippedNormal, MeshKernel::Projections::cartesian);

    // 3 Validation
    ASSERT_EQ(normal.x, 0.0);
    ASSERT_EQ(normal.y, 1.0);
    ASSERT_EQ(flippedNormal, true);
}
