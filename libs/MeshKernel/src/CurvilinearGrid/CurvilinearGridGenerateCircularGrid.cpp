#include "MeshKernel/CurvilinearGrid/CurvilinearGridGenerateCircularGrid.hpp"
#include "MeshKernel/Utilities/LinearAlgebra.hpp"

#include <cmath>

meshkernel::CurvilinearGrid meshkernel::CurvilinearGridGenerateCircularGrid::Compute(const double innerRadius, const double outerRadius, const MakeGridParameters& parameters) const
{
    // Check validity of
    // 1 radius
    // 2 parameters

    Projection projection = Projection::cartesian;

    return CurvilinearGrid(GenerateGridPoints(innerRadius, outerRadius, parameters, projection), projection);
}

std::vector<double> meshkernel::CurvilinearGridGenerateCircularGrid::GenerateRadiusValues(const double innerRadius, const double outerRadius, const UInt nRefinement) const
{
    std::vector<double> radValues(nRefinement);

    double radius = innerRadius;
    double radiusDelta = 2.0 * std::numbers::pi_v<double> * (outerRadius - innerRadius) / static_cast<double>(nRefinement - 1);
    // double radiusDelta = 2.0 * std::numbers::pi_v<double> * outerRadius / static_cast<double>(nRefinement);


    for (UInt i = 0; i < static_cast<UInt>(nRefinement); ++i)
    {
        radValues[i] = radius;
        radius += radiusDelta;
        radiusDelta *= 2.0;
        // radiusDelta += 2.0 * std::numbers::pi_v<double> * innerRadius / static_cast<double>(nRefinement);
    }

    double min = radValues[0];
    double max = radValues[radValues.size () - 1];
    double length = max - min;

    for (UInt i = 0; i < static_cast<UInt>(nRefinement); ++i)
    {
        radValues[i] = (radValues[i] - innerRadius) / length * (outerRadius - innerRadius) + innerRadius;
        std::cout << "radius value: " << radValues[i] << std::endl;
    }


    return radValues;
}

lin_alg::Matrix<meshkernel::Point> meshkernel::CurvilinearGridGenerateCircularGrid::GenerateGridPoints2(const double innerRadius, const double outerRadius, const MakeGridParameters& parameters, const Projection projection [[maybe_unused]]) const
{
    lin_alg::Matrix<Point> gridPoints(parameters.num_rows, parameters.num_columns);

    double deltaTheta = 2.0 * std::numbers::pi_v<double> / static_cast<double>(parameters.num_columns - 0);
    double deltaRadius = (outerRadius - innerRadius) / static_cast<double>(parameters.num_rows - 1);

    std::cout << "radius: " << outerRadius << "  " << innerRadius << "  " << deltaRadius << std::endl;

    std::vector<double> sinValues(parameters.num_columns);
    std::vector<double> cosValues(parameters.num_columns);
    std::vector<double> radValues(GenerateRadiusValues(innerRadius, outerRadius, static_cast<UInt>(parameters.num_rows)));

    double theta = 0.0;

    for (UInt i = 0; i < static_cast<UInt>(parameters.num_columns); ++i)
    {
        sinValues[i] = std::sin(theta);
        cosValues[i] = std::cos(theta);
        theta += deltaTheta;
    }

    double radius = innerRadius;

    // TODO I think the y-direction goes in the wrong direction (need size - n)

    for (UInt n = 0; n < static_cast<UInt>(parameters.num_rows); ++n)
    {

        for (UInt m = 0; m < static_cast<UInt>(parameters.num_columns); ++m)
        {
            gridPoints(n, m) = Point(radValues[n] * cosValues[m], radValues[n] * sinValues[m]);
            // gridPoints(n, m) = Point(radius * cosValues[m], radius * sinValues[m]);

            if (m == 1 && n > 0)
            {
                double xl = std::abs((gridPoints(n - 1, m).x - gridPoints(n, m).x));
                double yl = std::abs((gridPoints(n, m).y - gridPoints(n, m - 1).y));
                std::cout << "lengths: " << xl << "  " << yl << "  " << yl * yl << std::endl;
            }
        }

        radius += deltaRadius;
    }

    return gridPoints;
}

void meshkernel::CurvilinearGridGenerateCircularGrid::GenerateUniformGrid(const MakeGridParameters& parameters, const Projection projection, lin_alg::Matrix<meshkernel::Point>& gridPoints) const
{

    double dxa = 0.5 * parameters.block_size_x;
    double cs = std::cos(fi * constants::conversion::degToRad);
    double sn = std::sin(fi * constants::conversion::degToRad);
    UInt mc = static_cast<UInt>(parameters.num_rows + 1);
    UInt nc = static_cast<UInt>(parameters.num_columns + 1);


        if (parameter.fraction_rows != 1.0 && parameter.maximum_uniform_rows_size != 1.0)
        {
            muni = nint(mc*fmuni);

            if (mod((mc-1)-muni,2)==1)
            {
                muni = muni + 1;// ! at least muni cells, remainder is even
            }

            minc = (mc - muni)/2 + 1;
            dxm = fminc*dx;
            alfm = (dxm/dx)**(1.0_double/dble(minc));
        }


}

void meshkernel::CurvilinearGridGenerateCircularGrid::GenerateGradedGrid(const MakeGridParameters& parameters, const Projection projection, lin_alg::Matrix<meshkernel::Point>& gridPoints) const
{
}

lin_alg::Matrix<meshkernel::Point> meshkernel::CurvilinearGridGenerateCircularGrid::GenerateGridPoints(const double innerRadius, const double outerRadius, const MakeGridParameters& parameters, const Projection projection) const
{
    lin_alg::Matrix<Point> gridPoints(parameters.num_rows, parameters.num_columns);

    if (innerRadius == 0.0)
    {
        GenerateUniformGrid (parameters, projection, gridPoints);
    }
    else
    {
        GenerateGradedGrid (parameters, projection, gridPoints);
    }

    return gridPoints;
}
