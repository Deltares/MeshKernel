#include "MeshKernel/CurvilinearGrid/SnapGridToSpline.hpp"

meshkernel::SnapGridToSpline::SnapGridToSpline(CurvilinearGrid& grid,
                                               const Splines& spline,
                                               const std::vector<Point>& points) : CurvilinearGridSnapping(grid, points), m_spline(spline)
{
}

std::unique_ptr<meshkernel::MeshSmoothingCalculator> meshkernel::SnapGridToSpline::AllocateGridSmoothingCalculator(const CurvilinearGrid& originalGrid,
                                                                                                                   const CurvilinearGrid& snappedGrid) const
{
    return std::make_unique<NonDirectionalSmoothingCalculator>(originalGrid, snappedGrid, m_spline);
}

meshkernel::Point meshkernel::SnapGridToSpline::FindNearestPoint(const Point& currentPoint) const
{
    return m_spline.ComputeClosestPoint(0, currentPoint);
}
