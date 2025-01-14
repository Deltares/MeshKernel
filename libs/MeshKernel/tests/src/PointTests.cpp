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

#include <gtest/gtest.h>

#include <MeshKernel/Entities.hpp>

TEST(PointTests, IsEqualTest)
{
    const double tolerance = 1.0e-14;
    meshkernel::Point p1(1.0, 2.0);
    meshkernel::Point p2(1.0, 2.0);

    EXPECT_TRUE(meshkernel::IsEqual(p1, p2, tolerance));
    p1 += p2;
    EXPECT_FALSE(meshkernel::IsEqual(p1, p2, tolerance));
}

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

TEST(PointTests, NotInplaceAdd)
{
    meshkernel::Point p1(1.0, 2.0);
    meshkernel::Point p2(1.0, 2.0);
    meshkernel::Point p3;

    p3 = p1 + p2;
    EXPECT_EQ(2.0, p3.x);
    EXPECT_EQ(4.0, p3.y);

    p3 = p3 + p2;
    EXPECT_EQ(3.0, p3.x);
    EXPECT_EQ(6.0, p3.y);

    p3 = p3 + 1.0;
    EXPECT_EQ(4.0, p3.x);
    EXPECT_EQ(7.0, p3.y);
}

TEST(PointTests, NotInplaceSubtract)
{
    meshkernel::Point p1(10.0, 20.0);
    meshkernel::Point p2(1.0, 2.0);
    meshkernel::Point p3;

    p3 = p1 - p2;
    EXPECT_EQ(9.0, p3.x);
    EXPECT_EQ(18.0, p3.y);

    p3 = p3 - p2;
    EXPECT_EQ(8.0, p3.x);
    EXPECT_EQ(16.0, p3.y);

    p3 = p3 - 1.0;
    EXPECT_EQ(7.0, p3.x);
    EXPECT_EQ(15.0, p3.y);
}

TEST(PointTests, NotInplaceDivide)
{
    const double tolerance = 1.0e-14;

    meshkernel::Point p1(20.0, 20.0);
    meshkernel::Point p2(2.0, 5.0);
    meshkernel::Point p3;

    p3 = p1 / p2;
    EXPECT_NEAR(10.0, p3.x, tolerance);
    EXPECT_NEAR(4.0, p3.y, tolerance);

    p3 = p3 / 2.0;
    EXPECT_NEAR(5.0, p3.x, tolerance);
    EXPECT_NEAR(2.0, p3.y, tolerance);
}

TEST(PointTests, NotInplaceMultiply)
{
    meshkernel::Point p1(10.0, 20.0);
    meshkernel::Point p2(2.0, 5.0);
    meshkernel::Point p3;

    p3 = p1 * p2;
    EXPECT_EQ(20.0, p3.x);
    EXPECT_EQ(100.0, p3.y);

    p3 = p3 * 0.5;
    EXPECT_EQ(10.0, p3.x);
    EXPECT_EQ(50.0, p3.y);
}

TEST(PointTests, PointAlongLineTest)
{
    // two corners of a 3 by 4 right angled triangle.
    meshkernel::Point start(5.0, 6.0);
    meshkernel::Point end(8.0, 10.0);
    // Should be 5.
    const double hypotenuse = std::hypot(end.x - start.x, end.y - start.y);

    const int numberOfSteps = 10;
    const double h = hypotenuse / static_cast<double>(numberOfSteps - 1);
    double lambda = 0.0;

    const double tolerance = 1.0e-12;

    for (int i = 0; i < numberOfSteps; ++i)
    {
        meshkernel::Point actual = meshkernel::PointAlongLine(start, end, lambda);
        double expectedX = (1.0 - lambda) * start.x + lambda * end.x;
        double expectedY = (1.0 - lambda) * start.y + lambda * end.y;

        EXPECT_NEAR(expectedX, actual.x, tolerance);
        EXPECT_NEAR(expectedY, actual.y, tolerance);

        lambda += h;
    }
}

TEST(PointTests, RotateAboutPointTest)
{
    // A basic rotation test
    // test that a point starting at position (2,0) is rotated correctly with various angles.

    const double tolerance = 1.0e-12;

    const double radius = 2.0;
    const meshkernel::Point centre(0.0, 0.0);

    const int numberOfSteps = 10;
    const double delta = 360.0 / static_cast<double>(numberOfSteps - 1);
    double angle = 0.0;

    meshkernel::Point point(radius, 0.0);

    (void)tolerance;

    for (int i = 0; i < numberOfSteps; ++i)
    {
        meshkernel::Point actual = meshkernel::Rotate(point, angle, centre);
        double expectedX = radius * std::cos(angle * meshkernel::constants::conversion::degToRad);
        double expectedY = radius * std::sin(angle * meshkernel::constants::conversion::degToRad);

        EXPECT_NEAR(expectedX, actual.x, tolerance);
        EXPECT_NEAR(expectedY, actual.y, tolerance);

        angle += delta;
    }
}
