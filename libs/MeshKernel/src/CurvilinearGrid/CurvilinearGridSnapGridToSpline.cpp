#include "MeshKernel/CurvilinearGrid/CurvilinearGridSnapGridToSpline.hpp"

meshkernel::CurvilinearGridSnapGridToSpline::CurvilinearGridSnapGridToSpline(CurvilinearGrid& grid,
                                                                             const Splines& spline,
                                                                             const std::vector<Point>& points) : CurvilinearGridSnapping(grid, points), m_spline(spline)
{
}

std::unique_ptr<meshkernel::CurvilinearGridMeshExpansionCalculator> meshkernel::CurvilinearGridSnapGridToSpline::AllocateCurvilinearGridMeshExpansionCalculator(const CurvilinearGrid& originalGrid,
                                                                                                                                                                const CurvilinearGrid& snappedGrid) const
{
    return std::make_unique<DefaultRegionExpasionCalculator>(originalGrid, snappedGrid, m_spline);
}

meshkernel::Point meshkernel::CurvilinearGridSnapGridToSpline::FindNearestPoint(const Point& currentPoint) const
{
    return m_spline.ComputeClosestPoint(0, currentPoint);
}
