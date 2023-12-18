#include <gtest/gtest.h>

#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridSmoothness.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Utilities/LinearAlgebra.hpp"
#include "TestUtils/MakeCurvilinearGrids.hpp"

namespace mk = meshkernel;

std::unique_ptr<meshkernel::CurvilinearGrid> MakeMultiplicativeStretchedCurvilinearGrid(double originX, double originY, double deltaX, double deltaXIncrease, double deltaY, double deltaYIncrease, size_t nx, size_t ny)
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
            deltaXValue *= deltaXIncrease;
        }

        y += deltaY;
        deltaY *= deltaYIncrease;
    }

    return std::make_unique<meshkernel::CurvilinearGrid>(points, meshkernel::Projection::cartesian);
}

std::unique_ptr<meshkernel::CurvilinearGrid> MakeAdditiveStretchedCurvilinearGrid(double originX, double originY, double deltaX, double deltaXIncrease, double deltaY, double deltaYIncrease, size_t nx, size_t ny)
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
            deltaXValue += deltaXIncrease;
        }

        y += deltaY;
        deltaY += deltaYIncrease;
    }

    return std::make_unique<meshkernel::CurvilinearGrid>(points, meshkernel::Projection::cartesian);
}

TEST(CurvilinearGridSmoothness, SmoothnessOfRegularGrid)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 2.0;

    constexpr size_t nx = 27;
    constexpr size_t ny = 13;

    std::unique_ptr<mk::CurvilinearGrid> grid = MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);

    lin_alg::Matrix<double> smoothness;

    // First test smoothness in x-direction

    mk::CurvilinearGridSmoothness::Compute(*grid, 1, smoothness);

    ASSERT_EQ(nx, smoothness.rows());
    ASSERT_EQ(ny, smoothness.cols());

    constexpr double tolerance = 1.0e-12;

    for (long int i = 0; i < smoothness.rows(); ++i)
    {
        double expectedSmoothnessValue = (i == 0 || i == smoothness.rows() - 1 ? mk::constants::missing::doubleValue : 1.0);

        for (long int j = 0; j < smoothness.cols(); ++j)
        {
            EXPECT_NEAR(expectedSmoothnessValue, smoothness(i, j), tolerance);
        }
    }

    // Next test smoothness in y-direction

    mk::CurvilinearGridSmoothness::Compute(*grid, 2, smoothness);

    ASSERT_EQ(nx, smoothness.rows());
    ASSERT_EQ(ny, smoothness.cols());

    for (long int j = 0; j < smoothness.cols(); ++j)
    {
        double expectedSmoothnessValue = (j == 0 || j == smoothness.cols() - 1 ? mk::constants::missing::doubleValue : 1.0);

        for (long int i = 0; i < smoothness.rows(); ++i)
        {
            EXPECT_NEAR(expectedSmoothnessValue, smoothness(i, j), tolerance);
        }
    }
}

TEST(CurvilinearGridSmoothness, SmoothnessOfRegularGridWithHole)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 2.0;

    constexpr size_t nx = 10;
    constexpr size_t ny = 20;

    std::unique_ptr<mk::CurvilinearGrid> grid = MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);

    grid->m_gridNodes(4, 4).SetInvalid();
    grid->m_gridNodes(4, 5).SetInvalid();
    grid->m_gridNodes(5, 4).SetInvalid();
    grid->m_gridNodes(5, 5).SetInvalid();

    lin_alg::Matrix<double> smoothness;

    // First test smoothness in x-direction

    mk::CurvilinearGridSmoothness::Compute(*grid, 1, smoothness);

    ASSERT_EQ(nx, smoothness.rows());
    ASSERT_EQ(ny, smoothness.cols());

    constexpr double tolerance = 1.0e-12;

    for (long int i = 0; i < smoothness.rows(); ++i)
    {

        for (long int j = 0; j < smoothness.cols(); ++j)
        {
            double expectedSmoothnessValue = 1.0;

            if (i == 0 || i == smoothness.rows() - 1 || ((i == 3 || i == 4 || i == 5 || i == 6) && (j == 4 || j == 5)))
            {
                expectedSmoothnessValue = mk::constants::missing::doubleValue;
            }

            EXPECT_NEAR(expectedSmoothnessValue, smoothness(i, j), tolerance);
        }
    }

    // Next test smoothness in y-direction

    mk::CurvilinearGridSmoothness::Compute(*grid, 2, smoothness);

    ASSERT_EQ(nx, smoothness.rows());
    ASSERT_EQ(ny, smoothness.cols());

    for (long int j = 0; j < smoothness.cols(); ++j)
    {

        for (long int i = 0; i < smoothness.rows(); ++i)
        {
            double expectedSmoothnessValue = 1.0;

            if (j == 0 || j == smoothness.cols() - 1 || ((j == 3 || j == 4 || j == 5 || j == 6) && (i == 4 || i == 5)))
            {
                expectedSmoothnessValue = mk::constants::missing::doubleValue;
            }

            EXPECT_NEAR(expectedSmoothnessValue, smoothness(i, j), tolerance);
        }
    }
}

TEST(CurvilinearGridSmoothness, SmoothnessOfMultiplicativeIncreasingStretchedGrid)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 2.0;

    constexpr double deltaXIncrease = 1.25;
    constexpr double deltaYIncrease = 1.5;

    constexpr size_t nx = 6;
    constexpr size_t ny = 6;

    std::unique_ptr<mk::CurvilinearGrid> grid = MakeMultiplicativeStretchedCurvilinearGrid(originX, originY, deltaX, deltaXIncrease, deltaY, deltaYIncrease, nx, ny);

    lin_alg::Matrix<double> smoothness;

    // First test smoothness in x-direction

    mk::CurvilinearGridSmoothness::Compute(*grid, 1, smoothness);

    ASSERT_EQ(nx, smoothness.rows());
    ASSERT_EQ(ny, smoothness.cols());

    constexpr double tolerance = 1.0e-12;

    for (long int i = 0; i < smoothness.rows(); ++i)
    {
        double expectedSmoothnessValue = (i == 0 || i == smoothness.rows() - 1 ? mk::constants::missing::doubleValue : deltaXIncrease);

        for (long int j = 0; j < smoothness.cols(); ++j)
        {
            EXPECT_NEAR(expectedSmoothnessValue, smoothness(i, j), tolerance);
        }
    }

    // Next test smoothness in y-direction

    mk::CurvilinearGridSmoothness::Compute(*grid, 2, smoothness);

    ASSERT_EQ(nx, smoothness.rows());
    ASSERT_EQ(ny, smoothness.cols());

    for (long int i = 0; i < smoothness.rows(); ++i)
    {

        for (long int j = 0; j < smoothness.cols(); ++j)
        {
            double expectedSmoothnessValue = (j == 0 || j == smoothness.cols() - 1 ? mk::constants::missing::doubleValue : deltaYIncrease);
            EXPECT_NEAR(expectedSmoothnessValue, smoothness(i, j), tolerance);
        }
    }
}

TEST(CurvilinearGridSmoothness, SmoothnessOfMultiplicativeDecreasingStretchedGrid)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 2.0;

    constexpr double deltaXIncrease = 0.75;
    constexpr double deltaYIncrease = 0.5;

    constexpr size_t nx = 6;
    constexpr size_t ny = 6;

    std::unique_ptr<mk::CurvilinearGrid> grid = MakeMultiplicativeStretchedCurvilinearGrid(originX, originY, deltaX, deltaXIncrease, deltaY, deltaYIncrease, nx, ny);

    lin_alg::Matrix<double> smoothness;

    // First test smoothness in x-direction

    mk::CurvilinearGridSmoothness::Compute(*grid, 1, smoothness);

    ASSERT_EQ(nx, smoothness.rows());
    ASSERT_EQ(ny, smoothness.cols());

    constexpr double tolerance = 1.0e-12;

    for (long int i = 0; i < smoothness.rows(); ++i)
    {
        // 1.0 / deltaXIncrease is used because if the node smoothness value is less than 1.0 then the inverse is used
        double expectedSmoothnessValue = (i == 0 || i == smoothness.rows() - 1 ? mk::constants::missing::doubleValue : 1.0 / deltaXIncrease);

        for (long int j = 0; j < smoothness.cols(); ++j)
        {
            EXPECT_NEAR(expectedSmoothnessValue, smoothness(i, j), tolerance);
        }
    }

    // Next test smoothness in y-direction

    mk::CurvilinearGridSmoothness::Compute(*grid, 2, smoothness);

    ASSERT_EQ(nx, smoothness.rows());
    ASSERT_EQ(ny, smoothness.cols());

    for (long int i = 0; i < smoothness.rows(); ++i)
    {

        for (long int j = 0; j < smoothness.cols(); ++j)
        {
            // 1.0 / deltaYIncrease is used because if the node smoothness value is less than 1.0 then the inverse is used
            double expectedSmoothnessValue = (j == 0 || j == smoothness.cols() - 1 ? mk::constants::missing::doubleValue : 1.0 / deltaYIncrease);
            EXPECT_NEAR(expectedSmoothnessValue, smoothness(i, j), tolerance);
        }
    }
}

TEST(CurvilinearGridSmoothness, SmoothnessOfAdditiveStretchedGrid)
{
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 2.0;

    constexpr double deltaXIncrease = 1.25;
    constexpr double deltaYIncrease = 1.5;

    constexpr size_t nx = 6;
    constexpr size_t ny = 6;

    std::unique_ptr<mk::CurvilinearGrid> grid = MakeAdditiveStretchedCurvilinearGrid(originX, originY, deltaX, deltaXIncrease, deltaY, deltaYIncrease, nx, ny);

    // The expected smoothness for each of the test directions.
    std::vector<double> expectedSmoothnessX({-999.0, 2.25, 1.5555555555555555802, 1.3571428571428572063, 1.2631578947368420351, -999.0});
    std::vector<double> expectedSmoothnessY({-999.0, 1.75, 1.4285714285714286031, 1.3000000000000000444, 1.2307692307692308376, -999.0});

    lin_alg::Matrix<double> smoothness;

    // First test smoothness in x-direction

    mk::CurvilinearGridSmoothness::Compute(*grid, 1, smoothness);

    ASSERT_EQ(nx, smoothness.rows());
    ASSERT_EQ(ny, smoothness.cols());

    constexpr double tolerance = 1.0e-12;

    for (long int i = 0; i < smoothness.rows(); ++i)
    {
        for (long int j = 0; j < smoothness.cols(); ++j)
        {
            EXPECT_NEAR(expectedSmoothnessX[i], smoothness(i, j), tolerance);
        }
    }

    // Next test smoothness in y-direction

    mk::CurvilinearGridSmoothness::Compute(*grid, 2, smoothness);

    ASSERT_EQ(nx, smoothness.rows());
    ASSERT_EQ(ny, smoothness.cols());

    for (long int i = 0; i < smoothness.rows(); ++i)
    {
        for (long int j = 0; j < smoothness.cols(); ++j)
        {
            EXPECT_NEAR(expectedSmoothnessY[j], smoothness(i, j), tolerance);
        }
    }
}
