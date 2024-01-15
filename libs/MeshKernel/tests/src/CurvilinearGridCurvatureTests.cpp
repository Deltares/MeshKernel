//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2024.
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

#include <numbers>

#include <gtest/gtest.h>

#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridCurvature.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Operations.hpp"
#include "MeshKernel/Utilities/LinearAlgebra.hpp"
#include "TestUtils/MakeCurvilinearGrids.hpp"

namespace mk = meshkernel;

namespace curvature
{

    std::unique_ptr<meshkernel::CurvilinearGrid> MakeCurvilinearGrid(double originX, double originY, double deltaX, double deltaY, size_t nx, size_t ny)
    {
        double y = originY;

        lin_alg::Matrix<meshkernel::Point> points(nx, ny);

        for (size_t m = 0; m < ny; ++m)
        {
            double x = originX;
            double deltaXValue = deltaX;

            for (size_t n = 0; n < nx; ++n)
            {
                points(n, m) = meshkernel::Point(x, y);
                x += deltaXValue;
            }

            y += deltaY;
        }

        return std::make_unique<meshkernel::CurvilinearGrid>(points, meshkernel::Projection::cartesian);
    }

    std::unique_ptr<meshkernel::CurvilinearGrid> MakeCurvilinearGridCurveX(double originX, double originY, double deltaX, double deltaY, size_t nx, size_t ny)
    {
        double y = originY;
        double thetaDelta = std::numbers::pi_v<double> / static_cast<double>(nx);
        double theta = 0.0;

        lin_alg::Matrix<meshkernel::Point> points(nx, ny);

        for (size_t m = 0; m < ny; ++m)
        {
            double x = originX;
            double yCurve = 0.0;
            theta = 0.0;

            for (size_t n = 0; n < nx; ++n)
            {
                points(n, m) = meshkernel::Point(x, y + yCurve);
                x += deltaX;
                theta += thetaDelta;
                yCurve += 0.5 * deltaY * (std::sin(theta));
            }

            y += deltaY;
        }

        return std::make_unique<meshkernel::CurvilinearGrid>(points, meshkernel::Projection::cartesian);
    }

    std::unique_ptr<meshkernel::CurvilinearGrid> MakeCurvilinearGridCurveY(double originX, double originY, double deltaX, double deltaY, size_t nx, size_t ny)
    {
        double y = originY;
        double theta = 0.0;
        double thetaDelta = 2.0 * std::numbers::pi_v<double> / static_cast<double>(ny);
        double xCurve = 0.0;

        lin_alg::Matrix<meshkernel::Point> points(nx, ny);

        for (size_t m = 0; m < ny; ++m)
        {
            double x = originX;

            for (size_t n = 0; n < nx; ++n)
            {
                points(n, m) = meshkernel::Point(x + xCurve, y);
                x += deltaX;
            }

            theta += thetaDelta;
            xCurve += deltaX * (std::sin(theta));
            y += deltaY;
        }

        return std::make_unique<meshkernel::CurvilinearGrid>(points, meshkernel::Projection::cartesian);
    }

    std::unique_ptr<meshkernel::CurvilinearGrid> MakeCurvilinearGridCurveXandY(double originX, double originY, double deltaX, double deltaY, size_t nx, size_t ny)
    {
        double y = originY;
        double xTheta = 0.0;
        double xThetaDelta = 2.0 * std::numbers::pi_v<double> / static_cast<double>(ny);
        double yThetaDelta = 2.0 * std::numbers::pi_v<double> / static_cast<double>(nx);
        double xCurve = 0.0;

        lin_alg::Matrix<meshkernel::Point> points(nx, ny);

        for (size_t m = 0; m < ny; ++m)
        {
            double x = originX;
            double yTheta = 0.0;
            double yCurve = 0.0;

            for (size_t n = 0; n < nx; ++n)
            {
                points(n, m) = meshkernel::Point(x + xCurve, y + yCurve);
                x += deltaX;
                yTheta += yThetaDelta;
                yCurve += 0.5 * deltaY * (std::sin(yTheta));
            }

            xTheta += xThetaDelta;
            xCurve += 0.5 * deltaX * (std::sin(xTheta));
            y += deltaY;
        }

        return std::make_unique<meshkernel::CurvilinearGrid>(points, meshkernel::Projection::cartesian);
    }

} // namespace curvature

TEST(CurvilinearGridCurvature, Compute_CurvatureOfRegularGrid_ShouldComputeCurvature)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 2.0;

    constexpr size_t nx = 27;
    constexpr size_t ny = 13;

    std::unique_ptr<mk::CurvilinearGrid> grid = curvature::MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);

    lin_alg::Matrix<double> curvature;

    // First test curvature in m-direction

    mk::CurvilinearGridCurvature::Compute(*grid, mk::CurvilinearDirection::M, curvature);

    ASSERT_EQ(nx, curvature.rows());
    ASSERT_EQ(ny, curvature.cols());

    constexpr double tolerance = 1.0e-12;
    // This is the value for a mesh with no curvature
    constexpr double expectedCurvature = 1.0e3 / 999999.0;

    for (long int i = 0; i < curvature.rows(); ++i)
    {
        double expectedCurvatureValue = (i == 0 || i == curvature.rows() - 1 ? mk::constants::missing::doubleValue : expectedCurvature);

        for (long int j = 0; j < curvature.cols(); ++j)
        {
            EXPECT_NEAR(expectedCurvatureValue, curvature(i, j), tolerance);
        }
    }

    // Next test curvature in n-direction

    mk::CurvilinearGridCurvature::Compute(*grid, mk::CurvilinearDirection::N, curvature);

    ASSERT_EQ(nx, curvature.rows());
    ASSERT_EQ(ny, curvature.cols());

    for (long int j = 0; j < curvature.cols(); ++j)
    {
        double expectedCurvatureValue = (j == 0 || j == curvature.cols() - 1 ? mk::constants::missing::doubleValue : expectedCurvature);

        for (long int i = 0; i < curvature.rows(); ++i)
        {
            EXPECT_NEAR(expectedCurvatureValue, curvature(i, j), tolerance);
        }
    }
}

TEST(CurvilinearGridCurvature, Compute_CurvatureOfRegularGrid_ShouldComputeCurvatureInX)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 2.0;

    constexpr size_t nx = 27;
    constexpr size_t ny = 13;

    std::unique_ptr<mk::CurvilinearGrid> grid = curvature::MakeCurvilinearGridCurveX(originX, originY, deltaX, deltaY, nx, ny);

    // Expected values computed by algorithm
    std::vector<double> expectedCurvatureX{-999, 54.692933740071, 49.4363631174515,
                                           42.903008478042, 36.0517223491091, 29.5508089894376,
                                           23.7448255949468, 18.7333942685438, 14.4709747756264,
                                           10.8435956272394, 7.71491613343466, 4.94873108177472,
                                           2.41698338290781, 1.57408203884996e-13, 2.41698357191833,
                                           4.94873423903022, 7.71493326662676, 10.8436549935826,
                                           14.4711362937479, 18.7337700970048, 23.745601428529,
                                           29.552244069831, 36.0540780870016, 42.9063431041705,
                                           49.4402173759351, 54.6961820671308, -999};

    lin_alg::Matrix<double> curvature;

    // First test curvature in m-direction

    mk::CurvilinearGridCurvature::Compute(*grid, mk::CurvilinearDirection::M, curvature);

    ASSERT_EQ(nx, curvature.rows());
    ASSERT_EQ(ny, curvature.cols());

    constexpr double tolerance = 1.0e-11;
    // This is the value for a mesh with no curvature
    constexpr double expectedCurvature = 1.0e3 / 999999.0;

    for (long int i = 0; i < curvature.rows(); ++i)
    {
        double expectedCurvatureValue = (i == 0 || i == curvature.rows() - 1 ? mk::constants::missing::doubleValue : expectedCurvatureX[i]);

        for (long int j = 0; j < curvature.cols(); ++j)
        {
            EXPECT_NEAR(expectedCurvatureValue, curvature(i, j), tolerance);
        }
    }

    // Next test curvature in n-direction

    mk::CurvilinearGridCurvature::Compute(*grid, mk::CurvilinearDirection::N, curvature);

    ASSERT_EQ(nx, curvature.rows());
    ASSERT_EQ(ny, curvature.cols());

    for (long int j = 0; j < curvature.cols(); ++j)
    {
        double expectedCurvatureValue = (j == 0 || j == curvature.cols() - 1 ? mk::constants::missing::doubleValue : expectedCurvature);

        for (long int i = 0; i < curvature.rows(); ++i)
        {
            EXPECT_NEAR(expectedCurvatureValue, curvature(i, j), tolerance);
        }
    }
}

TEST(CurvilinearGridCurvature, Compute_CurvatureOfRegularGrid_ShouldComputeCurvatureInY)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 2.0;

    constexpr size_t nx = 27;
    constexpr size_t ny = 13;

    std::unique_ptr<mk::CurvilinearGrid> grid = curvature::MakeCurvilinearGridCurveY(originX, originY, deltaX, deltaY, nx, ny);

    // Expected values computed by algorithm
    std::vector<double> expectedCurvatureY{-999, 38.5180001187704, 16.0119151992422,
                                           5.27184246903049, 27.1867955829138, 48.9622125450272,
                                           59.4051462180544, 48.9397819786551, 27.1824308973952,
                                           5.27183432158012, 16.012547776799, 38.5311887716354,
                                           -999};

    lin_alg::Matrix<double> curvature;

    // First test curvature in m-direction

    mk::CurvilinearGridCurvature::Compute(*grid, mk::CurvilinearDirection::M, curvature);

    ASSERT_EQ(nx, curvature.rows());
    ASSERT_EQ(ny, curvature.cols());

    constexpr double tolerance = 1.0e-11;
    // This is the value for a mesh with no curvature
    constexpr double expectedCurvature = 1.0e3 / 999999.0;

    for (long int i = 0; i < curvature.rows(); ++i)
    {
        double expectedCurvatureValue = (i == 0 || i == curvature.rows() - 1 ? mk::constants::missing::doubleValue : expectedCurvature);

        for (long int j = 0; j < curvature.cols(); ++j)
        {
            EXPECT_NEAR(expectedCurvatureValue, curvature(i, j), tolerance);
        }
    }

    // Next test curvature in n-direction

    mk::CurvilinearGridCurvature::Compute(*grid, mk::CurvilinearDirection::N, curvature);

    ASSERT_EQ(nx, curvature.rows());
    ASSERT_EQ(ny, curvature.cols());

    for (long int j = 0; j < curvature.cols(); ++j)
    {
        double expectedCurvatureValue = (j == 0 || j == curvature.cols() - 1 ? mk::constants::missing::doubleValue : expectedCurvatureY[j]);

        for (long int i = 0; i < curvature.rows(); ++i)
        {
            EXPECT_NEAR(expectedCurvatureValue, curvature(i, j), tolerance);
        }
    }
}

TEST(CurvilinearGridCurvature, Compute_CurvatureOfRegularGrid_ShouldComputeCurvatureInXandY)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 1.0;

    constexpr size_t nx = 10;
    constexpr size_t ny = 10;

    std::unique_ptr<mk::CurvilinearGrid> grid = curvature::MakeCurvilinearGridCurveXandY(originX, originY, deltaX, deltaY, nx, ny);

    // Expected values computed by algorithm
    std::vector<double> expectedCurvatureX{-999, -999, -999, -999, -999, -999, -999, -999, -999, -999,
                                           73.65778037274387, 73.65778037274379, 73.65778037274383, 73.65778037274389,
                                           73.65778037274389, 73.65778037274389, 73.65778037274387, 73.65778037274387,
                                           73.65778037274387, 73.65778037274387, 4.088622216893358e-14, 1.635448886757343e-13,
                                           4.088622216893357e-14, 0.001000001000001, 0.001000001000001, 0.001000001000001,
                                           0.001000001000001, 0.001000001000001, 0.001000001000001, 0.001000001000001,
                                           73.68572111662814, 73.68572111662822, 73.68572111662814, 73.68572111662814,
                                           73.68572111662814, 73.68572111662814, 73.68572111662814, 73.68572111662814,
                                           73.68572111662814, 73.68572111662814, 140.9838081112912, 140.9838081112912,
                                           140.9838081112912, 140.9838081112912, 140.9838081112912, 140.9838081112912,
                                           140.9838081112912, 140.9838081112912, 140.9838081112912, 140.9838081112912,
                                           140.8629442842, 140.8629442842002, 140.8629442842002, 140.8629442842002,
                                           140.8629442842002, 140.8629442842002, 140.8629442842002, 140.8629442842002,
                                           140.8629442842002, 140.8629442842002, 73.65778037274399, 73.65778037274387,
                                           73.65778037274387, 73.65778037274387, 73.65778037274387, 73.65778037274387,
                                           73.65778037274387, 73.65778037274387, 73.65778037274387, 73.65778037274387,
                                           0.001000001000001, 0.001000001000001, 0.001000001000001, 0.001000001000001,
                                           0.001000001000001, 0.001000001000001, 0.001000001000001, 0.001000001000001,
                                           0.001000001000001, 0.001000001000001, 73.68572111662814, 73.68572111662814,
                                           73.68572111662814, 73.68572111662814, 73.68572111662814, 73.68572111662814,
                                           73.68572111662814, 73.68572111662814, 73.68572111662814, 73.68572111662814,
                                           -999, -999, -999, -999, -999, -999, -999, -999, -999, -999};

    // Expected values computed by algorithm
    std::vector<double> expectedCurvatureY{-999, -999, -999, -999, -999, -999, -999, -999, -999, -999,
                                           73.65778037274387, 73.65778037274379, 73.65778037274383, 73.65778037274389,
                                           73.65778037274389, 73.65778037274389, 73.65778037274387, 73.65778037274387,
                                           73.65778037274387, 73.65778037274387, 4.088622216893358e-14, 1.635448886757343e-13,
                                           4.088622216893357e-14, 0.001000001000001, 0.001000001000001, 0.001000001000001,
                                           0.001000001000001, 0.001000001000001, 0.001000001000001, 0.001000001000001,
                                           73.68572111662814, 73.68572111662822, 73.68572111662814, 73.68572111662814,
                                           73.68572111662814, 73.68572111662814, 73.68572111662814, 73.68572111662814,
                                           73.68572111662814, 73.68572111662814, 140.9838081112912, 140.9838081112912,
                                           140.9838081112912, 140.9838081112912, 140.9838081112912, 140.9838081112912,
                                           140.9838081112912, 140.9838081112912, 140.9838081112912, 140.9838081112912,
                                           140.8629442842, 140.8629442842002, 140.8629442842002, 140.8629442842002,
                                           140.8629442842002, 140.8629442842002, 140.8629442842002, 140.8629442842002,
                                           140.8629442842002, 140.8629442842002, 73.65778037274399, 73.65778037274387,
                                           73.65778037274387, 73.65778037274387, 73.65778037274387, 73.65778037274387,
                                           73.65778037274387, 73.65778037274387, 73.65778037274387, 73.65778037274387,
                                           0.001000001000001, 0.001000001000001, 0.001000001000001, 0.001000001000001,
                                           0.001000001000001, 0.001000001000001, 0.001000001000001, 0.001000001000001,
                                           0.001000001000001, 0.001000001000001, 73.68572111662814, 73.68572111662814,
                                           73.68572111662814, 73.68572111662814, 73.68572111662814, 73.68572111662814,
                                           73.68572111662814, 73.68572111662814, 73.68572111662814, 73.68572111662814,
                                           -999, -999, -999, -999, -999, -999, -999, -999, -999, -999};

    lin_alg::Matrix<double> curvature;

    // First test curvature in m-direction

    mk::CurvilinearGridCurvature::Compute(*grid, mk::CurvilinearDirection::M, curvature);

    ASSERT_EQ(nx, curvature.rows());
    ASSERT_EQ(ny, curvature.cols());

    constexpr double tolerance = 1.0e-11;

    size_t count = 0;

    for (long int i = 0; i < curvature.rows(); ++i)
    {
        for (long int j = 0; j < curvature.cols(); ++j)
        {
            EXPECT_NEAR(expectedCurvatureX[count], curvature(i, j), tolerance);
            ++count;
        }
    }

    // Next test curvature in n-direction

    mk::CurvilinearGridCurvature::Compute(*grid, mk::CurvilinearDirection::N, curvature);

    ASSERT_EQ(nx, curvature.rows());
    ASSERT_EQ(ny, curvature.cols());

    count = 0;

    for (long int j = 0; j < curvature.cols(); ++j)
    {
        for (long int i = 0; i < curvature.rows(); ++i)
        {
            EXPECT_NEAR(expectedCurvatureY[count], curvature(i, j), tolerance);
            ++count;
        }
    }
}
