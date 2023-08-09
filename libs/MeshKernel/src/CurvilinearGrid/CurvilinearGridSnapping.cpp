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

    Point meshDelta = m_grid.m_gridNodes[row][column] - m_originalGrid.m_gridNodes[row][column];

    int m1 = std::max<int>(1, row + 1 - nump * in) - 1;
    int m2 = std::min<int>(m_grid.m_gridNodes.size(), row + 1 + nump * in) - 1;
    int n1 = std::max<int>(1, column + 1 - nump * jn) - 1;
    int n2 = std::min<int>(m_grid.m_gridNodes[0].size(), column + 1 + nump * jn) - 1;

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

    double dsix = delta.y / (6.0 * aspectRatio);
    double rsx = std::max(dsix, std::hypot(meshDelta.x, meshDelta.y));

    for (int i = m1; i <= m2; ++i)
    {
        for (int j = n1; j <= n2; ++j)
        {
            Point pointDelta = m_originalGrid.m_gridNodes[i][j] - m_originalGrid.m_gridNodes[row][column];

            if (m_originalGrid.m_gridNodes[i][j].IsValid())
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

                        m_grid.m_gridNodes[i][j] = m_originalGrid.m_gridNodes[i][j] + Point(dx, dy);
                    }
                }
            }
        }
    }
}

meshkernel::CurvilinearGrid meshkernel::CurvilinearGridSnapping::Compute()
{
    constexpr double epsilon = 1.0e-5;

    int in = std::min<int>(1, m_upperRight.m_n - m_lowerLeft.m_n);
    int jn = std::min<int>(1, m_upperRight.m_m - m_lowerLeft.m_m);

    for (UInt i = m_lowerLeft.m_m; i <= m_upperRight.m_m; ++i)
    {
        for (UInt j = m_lowerLeft.m_n; j <= m_upperRight.m_n; ++j)
        {
            Point pnt = m_originalGrid.m_gridNodes[i][j];
            Point snappedPoint = m_landBoundary.FindNearestPoint(pnt, m_grid.m_projection);
            // Point snappedPoint = GetSnappedPoint(m_originalGrid.m_gridNodes[i][j]);
            m_grid.m_gridNodes[i][j] = snappedPoint;

            if (!IsEqual(snappedPoint.x, m_originalGrid.m_gridNodes[i][j].x, epsilon) || !IsEqual(snappedPoint.y, m_originalGrid.m_gridNodes[i][j].y, epsilon))
            {
                ModifyField(i, j, in, jn);
            }
        }
    }

    return m_grid;
}
