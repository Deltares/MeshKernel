#include "MeshKernel/CurvilinearGrid/CurvilinearGridGenerateCircularGrid.hpp"
#include "MeshKernel/Utilities/LinearAlgebra.hpp"

#include <cmath>

meshkernel::CurvilinearGrid meshkernel::CurvilinearGridGenerateCircularGrid::Compute(const double innerRadius, const double outerRadius, const CurvilinearParameters& parameters) const
{
    // Check validity of
    // 1 radius
    // 2 parameters

    return CurvilinearGrid(GenerateGridPoints(innerRadius, outerRadius, parameters), Projection::cartesian);
}

lin_alg::Matrix<meshkernel::Point> meshkernel::CurvilinearGridGenerateCircularGrid::GenerateGridPoints(const double innerRadius, const double outerRadius, const CurvilinearParameters& parameters) const
{
    lin_alg::Matrix<Point> gridPoints(parameters.n_refinement, parameters.m_refinement);

    double deltaTheta = 2.0 * std::numbers::pi_v<double> / static_cast<double>(parameters.m_refinement - 0);
    double deltaRadius = (outerRadius - innerRadius) / static_cast<double>(parameters.n_refinement - 1);

    std::cout << "radius: " << outerRadius << "  " << innerRadius << "  " << deltaRadius << std::endl;

    std::vector<double> sinValues(parameters.m_refinement);
    std::vector<double> cosValues(parameters.m_refinement);

    double theta = 0.0;

    for (UInt i = 0; i < static_cast<UInt>(parameters.m_refinement); ++i)
    {
        sinValues[i] = std::sin(theta);
        cosValues[i] = std::cos(theta);
        theta += deltaTheta;
    }

    double radius = innerRadius;

    // TODO I think the y-direction goes in the wrong direction (need size - n)

    for (UInt n = 0; n < static_cast<UInt>(parameters.n_refinement); ++n)
    {

        for (UInt m = 0; m < static_cast<UInt>(parameters.m_refinement); ++m)
        {
            gridPoints(n, m) = Point(radius * cosValues[m], radius * sinValues[m]);

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
