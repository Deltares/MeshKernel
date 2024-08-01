#include "MeshKernel/CurvilinearGrid/CurvilinearGridGenerateCircularGrid.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Utilities/LinearAlgebra.hpp"

#include <cmath>
#include <numbers>

meshkernel::CurvilinearGrid meshkernel::CurvilinearGridGenerateCircularGrid::Compute(const MakeGridParameters& parameters, const Projection projection)
{
    CheckMakeGridParameters(parameters);
    return CurvilinearGrid(GenerateGridPoints(parameters), projection);
}

std::vector<double> meshkernel::CurvilinearGridGenerateCircularGrid::ComputeXValues(const MakeGridParameters& parameters)
{
    std::vector<double> xValues(static_cast<UInt>(parameters.num_columns + 1));
    const double blockSize = parameters.block_size_x;

    if (parameters.uniform_columns_fraction == 1.0 || parameters.maximum_uniform_size_columns == 1.0)
    {
        // Regular in x-direction

        for (UInt m = 0; m < xValues.size(); ++m)
        {
            xValues[m] = blockSize * static_cast<double>(m);
        }
    }
    else
    {
        // Graded in x-direction

        const int size = parameters.num_columns + 1;
        int uniformElementCount = static_cast<int>(static_cast<double>(size) * parameters.uniform_columns_fraction);
        uniformElementCount += (size - 1 - uniformElementCount) % 2;
        const int uniformCutoffIndex = (size - uniformElementCount) / 2 + 1;
        const double gradingParameter = std::pow(parameters.maximum_uniform_size_columns, 1.0 / static_cast<double>(uniformCutoffIndex));

        double xValue = -blockSize * std::pow(gradingParameter, static_cast<double>(uniformCutoffIndex));

        auto computeDeltaX = [blockSize, gradingParameter, uniformElementCount, uniformCutoffIndex](const int m) -> double
        {
            double deltaX;

            if (m < uniformCutoffIndex)
            {
                deltaX = blockSize * std::pow(gradingParameter, static_cast<double>(uniformCutoffIndex - m));
            }
            else if (m < uniformCutoffIndex + uniformElementCount)
            {
                deltaX = blockSize;
            }
            else
            {
                deltaX = blockSize * std::pow(gradingParameter, static_cast<double>(m + 1 - uniformCutoffIndex - uniformElementCount));
            }

            return deltaX;
        };

        for (int m = 0; m < static_cast<int>(xValues.size()); ++m)
        {
            xValue += computeDeltaX(m);
            xValues[m] = xValue;
        }
    }

    return xValues;
}

std::vector<double> meshkernel::CurvilinearGridGenerateCircularGrid::ComputeYValues(const MakeGridParameters& parameters)
{
    std::vector<double> yValues(static_cast<UInt>(parameters.num_rows + 1));
    double blockSize = parameters.block_size_y;

    if (parameters.uniform_rows_fraction == 1.0 || parameters.maximum_uniform_size_rows == 1.0)
    {
        // Regular in y-direction

        for (UInt n = 0; n < yValues.size(); ++n)
        {
            yValues[n] = blockSize * static_cast<double>(n);
        }
    }
    else
    {
        // Graded in y-direction

        int size = parameters.num_rows + 1;
        int uniformElementCount = static_cast<int>(static_cast<double>(size) * parameters.uniform_rows_fraction);
        uniformElementCount += (size - 1 - uniformElementCount) % 2;
        int uniformCutoffIndex = (size - uniformElementCount) / 2 + 1;
        double gradingParameter = std::pow(parameters.maximum_uniform_size_rows, 1.0 / static_cast<double>(uniformCutoffIndex));
        double yValue = -blockSize * std::pow(gradingParameter, static_cast<double>(1 - uniformElementCount));

        for (int n = 0; n < static_cast<int>(yValues.size()); ++n)
        {
            yValue += (n < uniformElementCount ? blockSize : blockSize * std::pow(gradingParameter, static_cast<double>(n + 1 - uniformElementCount)));
            yValues[static_cast<UInt>(n)] = yValue;
        }
    }

    return yValues;
}

lin_alg::Matrix<meshkernel::Point> meshkernel::CurvilinearGridGenerateCircularGrid::GenerateGradedRectangularGrid(const MakeGridParameters& parameters)
{
    lin_alg::Matrix<Point> gridPoints(parameters.num_rows + 1, parameters.num_columns + 1);

    const double x0 = parameters.origin_x;
    const double y0 = parameters.origin_y;

    const double phi = parameters.angle;
    const double cosTheta = std::cos(phi * constants::conversion::degToRad);
    const double sinTheta = std::sin(phi * constants::conversion::degToRad);

    const std::vector<double> xValues(ComputeXValues(parameters));
    const std::vector<double> yValues(ComputeYValues(parameters));

    for (int n = 0; n < parameters.num_rows + 1; ++n)
    {
        for (int m = 0; m < parameters.num_columns + 1; ++m)
        {
            gridPoints(n, m) = Point(x0 + xValues[m] * cosTheta - yValues[n] * sinTheta,
                                     y0 + xValues[m] * sinTheta + yValues[n] * cosTheta);
        }
    }

    return gridPoints;
}

std::vector<double> meshkernel::CurvilinearGridGenerateCircularGrid::ComputeRadiusValues(const MakeGridParameters& parameters)
{
    const UInt size = static_cast<UInt>(parameters.num_rows + 1);
    const double blockSize = parameters.block_size_y;
    const double r0 = parameters.radius_curvature;

    std::vector<double> radiusValues(size);

    if (parameters.uniform_rows_fraction == 1.0 || parameters.maximum_uniform_size_rows == 1.0)
    {
        double radius = r0;

        for (UInt n = 0; n < size; ++n)
        {
            radiusValues[n] = radius;
            radius += blockSize;
        }
    }
    else
    {
        double radius = r0;

        UInt uniformElementCount = static_cast<UInt>(static_cast<double>(size) * parameters.uniform_rows_fraction);
        uniformElementCount += (size - 1 - uniformElementCount) % 2;
        UInt uniformCutoffIndex = (size - uniformElementCount) / 2 + 1;

        double gradingParameter = std::pow(parameters.maximum_uniform_size_rows, 1.0 / static_cast<double>(uniformCutoffIndex));

        for (UInt n = 0; n < size; ++n)
        {
            radius += (n < uniformElementCount ? blockSize : blockSize * std::pow(gradingParameter, static_cast<double>(n - uniformElementCount)));
            radiusValues[n] = radius;
        }
    }

    return radiusValues;
}

std::vector<double> meshkernel::CurvilinearGridGenerateCircularGrid::ComputeThetaValues(const MakeGridParameters& parameters)
{
    const UInt size = static_cast<UInt>(parameters.num_columns + 1);
    std::vector<double> thetaValues(size);

    double deltaTheta = parameters.block_size_x / parameters.radius_curvature;

    if (deltaTheta * static_cast<double>(size - 1) > 2.0 * std::numbers::pi_v<double>)
    {
        // Log message here: "Increase radius of curvature or decrease M-dimension or delta X"
        deltaTheta = 2.0 * std::numbers::pi_v<double> / static_cast<double>(size - 1);
    }

    double offset = parameters.angle * constants::conversion::degToRad + 0.5 * std::numbers::pi_v<double>;

    for (UInt m = 0; m < size; ++m)
    {
        double theta = offset - static_cast<double>(m) * deltaTheta;
        thetaValues[m] = theta;
    }

    return thetaValues;
}

lin_alg::Matrix<meshkernel::Point> meshkernel::CurvilinearGridGenerateCircularGrid::GenerateCircularGrid(const MakeGridParameters& parameters)
{
    lin_alg::Matrix<Point> gridPoints(parameters.num_rows + 1, parameters.num_columns + 1);

    const double x0 = parameters.origin_x;
    const double y0 = parameters.origin_y;

    const std::vector<double> radiusValues(ComputeRadiusValues(parameters));
    const std::vector<double> thetaValues(ComputeThetaValues(parameters));
    std::vector<double> cosValues(thetaValues.size());
    std::vector<double> sinValues(thetaValues.size());

    for (UInt i = 0; i < thetaValues.size(); ++i)
    {
        cosValues[i] = std::cos(thetaValues[i]);
        sinValues[i] = std::sin(thetaValues[i]);
    }

    for (UInt n = 0; n < static_cast<UInt>(parameters.num_rows + 1); ++n)
    {
        double r = radiusValues[n];

        for (UInt m = 0; m < static_cast<UInt>(parameters.num_columns + 1); ++m)
        {
            gridPoints(n, m) = Point(x0 + r * cosValues[m], y0 + r * sinValues[m]);
        }
    }

    return gridPoints;
}

lin_alg::Matrix<meshkernel::Point> meshkernel::CurvilinearGridGenerateCircularGrid::GenerateGridPoints(const MakeGridParameters& parameters)
{
    if (parameters.radius_curvature == 0.0)
    {
        return GenerateGradedRectangularGrid(parameters);
    }
    else
    {
        return GenerateCircularGrid(parameters);
    }
}
