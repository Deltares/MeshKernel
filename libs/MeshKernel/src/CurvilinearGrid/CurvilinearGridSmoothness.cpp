#include <cmath>

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Vector.hpp"

#include "MeshKernel/CurvilinearGrid/CurvilinearGridSmoothness.hpp"

void meshkernel::CurvilinearGridSmoothness::Compute(const CurvilinearGrid& grid, const int direction, lin_alg::Matrix<double>& smoothness)
{
    smoothness.resize(grid.m_numM, grid.m_numN);
    smoothness.fill(constants::missing::doubleValue);

    if (direction == 1)
    {
        for (UInt i = 1; i < grid.m_numM - 1; ++i)
        {
            for (UInt j = 0; j < grid.m_numN; ++j)
            {
                smoothness(i, j) = ComputeNodeSmoothness(grid.m_gridNodes(i - 1, j), grid.m_gridNodes(i, j), grid.m_gridNodes(i + 1, j));
            }
        }
    }
    else
    {
        for (UInt i = 0; i < grid.m_numM; ++i)
        {
            for (UInt j = 1; j < grid.m_numN - 1; ++j)
            {
                smoothness(i, j) = ComputeNodeSmoothness(grid.m_gridNodes(i, j - 1), grid.m_gridNodes(i, j), grid.m_gridNodes(i, j + 1));
            }
        }
    }
}

double meshkernel::CurvilinearGridSmoothness::ComputeNodeSmoothness(const Point& p0, const Point& p1, const Point& p2)
{
    double nodeSmoothness = constants::missing::doubleValue;

    if (p0.IsValid() && p1.IsValid() && p2.IsValid())
    {
        double diffX21 = p2.x - p1.x;
        double diffY21 = p2.y - p1.y;

        double diffX10 = p1.x - p0.x;
        double diffY10 = p1.y - p0.y;

        double lengthSquared21 = diffX21 * diffX21 + diffY21 * diffY21;
        double lengthSquared10 = diffX10 * diffX10 + diffY10 * diffY10;

        if (lengthSquared10 != 0.0 && lengthSquared21 != 0.0)
        {
            nodeSmoothness = std::sqrt(lengthSquared21 / lengthSquared10);

            if (nodeSmoothness < 1.0)
            {
                nodeSmoothness = 1.0 / nodeSmoothness;
            }
        }
        else
        {
            nodeSmoothness = 1.0;
        }
    }

    return nodeSmoothness;
}
