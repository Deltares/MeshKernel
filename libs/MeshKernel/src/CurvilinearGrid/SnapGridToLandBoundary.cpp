#include "MeshKernel/CurvilinearGrid/SnapGridToLandBoundary.hpp"

meshkernel::SnapGridToLandBoundary::SnapGridToLandBoundary(CurvilinearGrid& grid,
                                                           const LandBoundary& landBoundary,
                                                           const std::vector<Point>& points) : CurvilinearGridSnapping(grid, points),
                                                                                               m_projection(grid.projection()),
                                                                                               m_landBoundary(landBoundary)
{
}

std::unique_ptr<meshkernel::MeshExpansionCalculator> meshkernel::SnapGridToLandBoundary::AllocateMeshExpansionCalculator(const CurvilinearGrid& originalGrid,
                                                                                                                         const CurvilinearGrid& snappedGrid) const
{
    return std::make_unique<DefaultRegionExpasionCalculator>(originalGrid, snappedGrid, m_landBoundary);
}

meshkernel::Point meshkernel::SnapGridToLandBoundary::FindNearestPoint(const Point& currentPoint) const
{
    return m_landBoundary.FindNearestPoint(currentPoint, m_projection);
}
