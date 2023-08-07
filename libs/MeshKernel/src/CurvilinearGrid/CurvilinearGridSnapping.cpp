#include "MeshKernel/CurvilinearGrid/CurvilinearGridSnapping.hpp"

#include "MeshKernel/BoundingBox.hpp"

#include <cmath>
#include <numbers>

meshkernel::CurvilinearGridSnapping::CurvilinearGridSnapping(std::shared_ptr<CurvilinearGrid> grid,
                                                             const LandBoundary& lb) : CurvilinearGridAlgorithm(grid), m_originalGrid(*grid), m_landBoundary(lb) {}

void meshkernel::CurvilinearGridSnapping::ModifyField(const UInt row, const UInt column, const int in, const int jn)
{

    // Value from editgridlineblok.f90
    // TODO Why is this variable called nump?
    constexpr UInt nump = 80;
    constexpr double aspectRatio = 990.0 / 1600.0;

    Point meshDelta = m_grid.m_gridNodes[column][row] - m_originalGrid.m_gridNodes[column][row];

    double dsix = 32.3232336076422; // 30.011936;

    // This is not the real value for dsix.
    double rsx = std::max(dsix, std::hypot(meshDelta.x, meshDelta.y));
    // double rsx = std::max((m_originalGrid.m_gridNodes[0][0].x - m_originalGrid.m_gridNodes[1][0].x) / 6.0, std::hypot(meshDelta.x, meshDelta.y));
    // double rsx = std::max(m_originalGrid.getDeltaX() / 6.0, std::hypot(meshDelta.x, meshDelta.y));

    int m1 = std::max<int>(1, column + 1 - nump * in) - 1;
    int m2 = std::min<int>(m_grid.m_gridNodes.size(), column + 1 + nump * in) - 1;
    int n1 = std::max<int>(1, row + 1 - nump * jn) - 1;
    int n2 = std::min<int>(m_grid.m_gridNodes[0].size(), row + 1 + nump * jn) - 1;

    // UInt m1 = std::max<UInt>(1, row - nump * in) - 1;
    // UInt m2 = std::min<UInt>(m_grid.m_gridNodes.size(), row + nump * in) - 1;
    // UInt n1 = std::max<UInt>(1, column - nump * jn) - 1;
    // UInt n2 = std::min<UInt>(m_grid.m_gridNodes[0].size(), column + nump * jn) - 1;

    BoundingBox gridBb = m_grid.GetBoundingBox();
    BoundingBox lbBb = m_landBoundary.GetBoundingBox();

    [[maybe_unused]] BoundingBox merged = meshkernel::merge(gridBb, lbBb);
    [[maybe_unused]] Point delta = merged.delta();
    [[maybe_unused]] Point centre = merged.centre();

    delta *= 1.2;

    if (delta.y < aspectRatio * delta.x)
    {
        delta.y = aspectRatio * delta.x;
    }

    dsix = delta.y / (6.0 * aspectRatio);

    for (int i = n1; i <= n2; ++i)
    {
        for (int j = m1; j <= m2; ++j)
        {
            Point pointDelta = m_originalGrid.m_gridNodes[j][i] - m_originalGrid.m_gridNodes[column][row];

            if (m_originalGrid.m_gridNodes[j][i].IsValid())
            {
                // m_grid.m_gridNodes[i][j] += ComputeDelta();
                // m_grid.m_gridNodes[i][j] += meshDelta.Compute();

                // if (ismeer)
                // {
                //     auto [scalingH, scalingV, scalingMixed] = m_originalGrid.ComputeDirectionalSmoothingFactors();

                //     dx = dx0 * scalingMixed;
                //     dy = dy0 * scalingMixed;

                //     m_grid.m_gridNodes[i][j] += Point(dx, dy);
                // }
                // else
                {
                    double rn = std::hypot(pointDelta.x, pointDelta.y);

                    if (rn < rsx)
                    {
                        rn = 0.5 * std::numbers::pi * rn / rsx;
                        double fr = 0.5 * (1.0 + std::cos(2.0 * rn));
                        double dx = meshDelta.x * fr;
                        double dy = meshDelta.y * fr;

                        m_grid.m_gridNodes[j][i] = m_originalGrid.m_gridNodes[j][i] + Point(dx, dy);
                    }
                }
            }
        }
    }

    //
}

meshkernel::CurvilinearGrid meshkernel::CurvilinearGridSnapping::Compute()
{
    constexpr double epsilon = 1.0e-5;

    // How to get m_lowerLeft and m_upperRight
    // This is just a grid line

    int in = std::min<int>(1, m_upperRight.m_n - m_lowerLeft.m_n);
    int jn = std::min<int>(1, m_upperRight.m_m - m_lowerLeft.m_m);

    for (UInt i = m_lowerLeft.m_n; i <= m_upperRight.m_n; ++i)
    {
        for (UInt j = m_lowerLeft.m_m; j <= m_upperRight.m_m; ++j)
        {
            Point pnt = m_originalGrid.m_gridNodes[j][i];
            Point snappedPoint = m_landBoundary.FindNearestPoint(pnt, m_grid.m_projection);
            // Point snappedPoint = GetSnappedPoint(m_originalGrid.m_gridNodes[i][j]);
            m_grid.m_gridNodes[j][i] = snappedPoint;

            if (!IsEqual(snappedPoint.x, m_originalGrid.m_gridNodes[j][i].x, epsilon) || !IsEqual(snappedPoint.y, m_originalGrid.m_gridNodes[j][i].y, epsilon))
            {
                ModifyField(i, j, in, jn);
            }
        }
    }

    std::cout.precision(6);
    std::cout.setf(std::ios::scientific);

    for (int i = (int)m_grid.m_gridNodes.size() - 1; i >= 0; --i)
    // for (UInt i = 0; i < m_grid.m_gridNodes.size(); ++i)
    {
        std::cout << "row -> ";
        for (UInt j = 0; j < m_grid.m_gridNodes.size(); ++j)
        {
            std::cout << "{" << m_grid.m_gridNodes[j][i].x << ", " << m_grid.m_gridNodes[j][i].y << "} ";
        }
        std::cout << std::endl;
    }

    return m_grid;
}
