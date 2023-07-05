#include <gtest/gtest.h>

#include <MeshKernel/Entities.hpp>

TEST(PointTests, InplaceAdd)
{
    meshkernel::Point p1(1.0, 2.0);
    meshkernel::Point p2(1.0, 2.0);

    p1 += p2;
    EXPECT_EQ(2.0, p1.x);
    EXPECT_EQ(4.0, p1.y);

    p1 += p2;
    EXPECT_EQ(3.0, p1.x);
    EXPECT_EQ(6.0, p1.y);

    p1 += 1.0;
    EXPECT_EQ(4.0, p1.x);
    EXPECT_EQ(7.0, p1.y);
}

TEST(PointTests, InplaceSubtract)
{
    meshkernel::Point p1(10.0, 20.0);
    meshkernel::Point p2(1.0, 2.0);

    p1 -= p2;
    EXPECT_EQ(9.0, p1.x);
    EXPECT_EQ(18.0, p1.y);

    p1 -= p2;
    EXPECT_EQ(8.0, p1.x);
    EXPECT_EQ(16.0, p1.y);

    p1 -= 1.0;
    EXPECT_EQ(7.0, p1.x);
    EXPECT_EQ(15.0, p1.y);
}

TEST(PointTests, InplaceDivide)
{
    const double tolerance = 1.0e-14;

    meshkernel::Point p1(20.0, 20.0);
    meshkernel::Point p2(2.0, 5.0);

    p1 /= p2;
    EXPECT_NEAR(10.0, p1.x, tolerance);
    EXPECT_NEAR(4.0, p1.y, tolerance);

    p1 /= 2.0;
    EXPECT_NEAR(5.0, p1.x, tolerance);
    EXPECT_NEAR(2.0, p1.y, tolerance);
}

TEST(PointTests, InplaceMultiply)
{
    meshkernel::Point p1(10.0, 20.0);
    meshkernel::Point p2(2.0, 5.0);

    p1 *= p2;
    EXPECT_EQ(20.0, p1.x);
    EXPECT_EQ(100.0, p1.y);

    p1 *= 0.5;
    EXPECT_EQ(10.0, p1.x);
    EXPECT_EQ(50.0, p1.y);
}
