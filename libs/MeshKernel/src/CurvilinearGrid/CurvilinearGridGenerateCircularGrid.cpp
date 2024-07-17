#include "MeshKernel/CurvilinearGrid/CurvilinearGridGenerateCircularGrid.hpp"
#include "MeshKernel/Utilities/LinearAlgebra.hpp"

#include <cmath>
#include <numbers>

meshkernel::CurvilinearGrid meshkernel::CurvilinearGridGenerateCircularGrid::Compute(const MakeGridParameters& parameters) const
{
    // Check validity of
    // 1 parameters

    Projection projection = Projection::cartesian;

    return CurvilinearGrid(GenerateGridPoints(parameters), projection);
}

lin_alg::Matrix<meshkernel::Point> meshkernel::CurvilinearGridGenerateCircularGrid::GenerateGradedGrid(const MakeGridParameters& parameters) const
{

    lin_alg::Matrix<Point> gridPoints(parameters.num_rows + 1, parameters.num_columns + 1);

    double x0 = parameters.origin_x;
    double y0 = parameters.origin_y;

    double phi = parameters.left_rotation;
    double cs = std::cos(phi * constants::conversion::degToRad);
    double sn = std::sin(phi * constants::conversion::degToRad);
    int mc = parameters.num_rows + 1;
    int nc = parameters.num_columns + 1;

    double fnuni = std::max(0.0, parameters.fraction_rows);
    double fmuni = std::max(0.0, parameters.fraction_columns);
    double fninc = std::max(0.0, parameters.maximum_uniform_rows_size);
    double fminc = std::max(0.0, parameters.maximum_uniform_columns_size);
    double dx = parameters.block_size_x;
    double dy = parameters.block_size_y;
    double dxc;
    double dyc;
    double alfm;
    double alfn;

    int nuni = 1000000000;
    int muni = 1000000000;
    int ninc = 1000000000;
    int minc = 1000000000;

    if (parameters.fraction_columns != 1.0 && parameters.maximum_uniform_columns_size != 1.0)
    {
        muni = static_cast<UInt>(static_cast<double>(mc) * fmuni);
        muni += (mc - 1 - muni) % 2;
        minc = (mc - muni) / 2 + 1;
        alfm = std::pow(fminc, 1.0 / static_cast<double>(minc));
    }

    if (parameters.fraction_rows != 1.0 && parameters.maximum_uniform_rows_size != 1.0)
    {
        nuni = static_cast<UInt>(static_cast<double>(nc) * fnuni);
        nuni += (nc - 1 - nuni) % 2;
        ninc = (nc - nuni) / 2 + 1;
        alfn = std::pow(fninc, 1.0 / static_cast<double>(ninc));
    }

    double yy = 0.0;

    for (int n = 0; n < nc; ++n)
    {
        double xx = 0.0;

        if (fnuni == 1.0 || fninc == 1.0)
        {
            yy = dy * static_cast<double>(n);
        }
        else
        {
            if (n < nuni)
            {
                dyc = dy;
            }
            else
            {
                dyc = dy * std::pow(alfn, static_cast<double>(n + 1 - nuni));
            }

            if (n == 0)
            {
                yy = -dy * std::pow(alfn, static_cast<double>(1 - nuni));
            }

            yy += dyc;
        }

        for (int m = 0; m < mc; ++m)
        {

            if (fmuni == 1.0 || fminc == 1.0)
            {
                xx = dx * static_cast<double>(m);
            }
            else
            {
                if (m < minc)
                {
                    dxc = dx * std::pow(alfm, static_cast<double>(minc - m));
                }
                else if (m < minc + muni)
                {
                    dxc = dx;
                }
                else
                {
                    dxc = dx * std::pow(alfm, static_cast<double>(m + 1 - minc - muni));
                }

                if (m == 0)
                {
                    xx = -dx * std::pow(alfm, static_cast<double>(minc));
                }

                xx += dxc;
            }

            gridPoints(n, m) = Point(x0 + xx * cs - yy * sn, y0 + xx * sn + yy * cs);
        }
    }

    return gridPoints;
}

std::vector<double> meshkernel::CurvilinearGridGenerateCircularGrid::ComputeRadiusValues(const MakeGridParameters& parameters) const
{
    const UInt size = static_cast<UInt>(parameters.num_rows + 1);
    const double deltaY = parameters.block_size_y;
    const double r0 = parameters.column_curvature_radius;

    std::vector<double> radiusValues(size);

    if (parameters.fraction_rows == 1.0 || parameters.maximum_uniform_rows_size == 1.0)
    {
        double yy = r0;

        for (UInt n = 0; n < size; ++n)
        {
            radiusValues[n] = yy;
            yy += deltaY;
        }
    }
    else
    {
        double yy = r0;
        double dyc = 0.0;

        UInt nuni = static_cast<UInt>(static_cast<double>(size) * parameters.fraction_rows);
        nuni += (size - 1 - nuni) % 2;
        UInt ninc = (size - nuni) / 2 + 1;

        double alfn = std::pow(parameters.maximum_uniform_rows_size, 1.0 / static_cast<double>(ninc));

        for (UInt n = 0; n < size; ++n)
        {
            if (n < nuni)
            {
                dyc = deltaY;
            }
            else
            {
                dyc = deltaY * std::pow(alfn, static_cast<double>(n - nuni));
            }

            yy += dyc;
            radiusValues[n] = yy;
        }
    }

    return radiusValues;
}

std::vector<double> meshkernel::CurvilinearGridGenerateCircularGrid::ComputeThetaValues(const MakeGridParameters& parameters) const
{
    const UInt size = static_cast<UInt>(parameters.num_columns + 1);
    std::vector<double> thetaValues(size);

    double deltaTheta = parameters.block_size_x / parameters.column_curvature_radius;

    if (deltaTheta * static_cast<double>(size - 1) > 2.0 * std::numbers::pi_v<double>)
    {
        // call qnmessage('Increase radius of curvature or' // ' decrease M-dimension or delta X')
        deltaTheta = 2.0 * std::numbers::pi_v<double> / static_cast<double>(size - 1);
    }

    double offset = parameters.left_rotation * constants::conversion::degToRad + 0.5 * std::numbers::pi_v<double>;

    for (UInt m = 0; m < size; ++m)
    {
        double theta = offset - static_cast<double>(m) * deltaTheta;
        thetaValues[m] = theta;
    }
    return thetaValues;
}

lin_alg::Matrix<meshkernel::Point> meshkernel::CurvilinearGridGenerateCircularGrid::GenerateRadialGrid(const MakeGridParameters& parameters) const
{
    lin_alg::Matrix<Point> gridPoints(parameters.num_rows + 1, parameters.num_columns + 1);

    double x0 = parameters.origin_x;
    double y0 = parameters.origin_y;

    std::vector<double> radiusValues(ComputeRadiusValues(parameters));
    std::vector<double> thetaValues(ComputeThetaValues(parameters));
    std::vector<double> cosValues(radiusValues.size());
    std::vector<double> sinValues(radiusValues.size());

    for (UInt i = 0; i < radiusValues.size(); ++i)
    {
        cosValues[i] = std::cos(thetaValues[i]);
        sinValues[i] = std::sin(thetaValues[i]);
    }

    for (UInt n = 0; n < static_cast<UInt>(parameters.num_columns + 1); ++n)
    {
        double r = radiusValues[n];

        for (UInt m = 0; m < static_cast<UInt>(parameters.num_rows + 1); ++m)
        {
            gridPoints(n, m) = Point(x0 + r * cosValues[m], y0 + r * sinValues[m]);
        }
    }

    return gridPoints;
}

lin_alg::Matrix<meshkernel::Point> meshkernel::CurvilinearGridGenerateCircularGrid::GenerateGridPoints(const MakeGridParameters& parameters) const
{

    if (parameters.column_curvature_radius == 0.0)
    {
        return GenerateGradedGrid(parameters);
    }
    else
    {
        return GenerateRadialGrid(parameters);
    }
}
