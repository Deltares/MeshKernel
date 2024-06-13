#include "MeshKernel/CurvilinearGrid/CurvilinearGridSnapGridToLandBoundary.hpp"

meshkernel::CurvilinearGridSnapGridToLandBoundary::CurvilinearGridSnapGridToLandBoundary(CurvilinearGrid& grid,
                                                                                         const LandBoundary& landBoundary,
                                                                                         const std::vector<Point>& points) : CurvilinearGridSnapping(grid, points),
                                                                                                                             m_projection(grid.projection()),
                                                                                                                             m_landBoundary(landBoundary)
{
}

std::unique_ptr<meshkernel::CurvilinearGridMeshExpansionCalculator> meshkernel::CurvilinearGridSnapGridToLandBoundary::AllocateCurvilinearGridMeshExpansionCalculator(const CurvilinearGrid& originalGrid,
                                                                                                                                                                      const CurvilinearGrid& snappedGrid) const
{
    return std::make_unique<DefaultRegionExpasionCalculator>(originalGrid, snappedGrid, m_landBoundary);
}

meshkernel::Point meshkernel::CurvilinearGridSnapGridToLandBoundary::FindNearestPoint(const Point& currentPoint) const
{
    return m_landBoundary.FindNearestPoint(currentPoint, m_projection);
}
