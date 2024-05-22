#include "MeshKernel/CurvilinearGrid/MeshSmoothingCalculator.hpp"

#include <numbers>

meshkernel::DirectionalSmoothingCalculator::DirectionalSmoothingCalculator(const CurvilinearGridNodeIndices& lowerLeft,
                                                                           const CurvilinearGridNodeIndices& upperRight,
                                                                           const CurvilinearGridNodeIndices& regionIndicator) : m_indexBoxLowerLeft(lowerLeft),
                                                                                                                                m_indexBoxUpperRight(upperRight),
                                                                                                                                m_smoothingRegionIndicator(regionIndicator) {}

double meshkernel::DirectionalSmoothingCalculator::compute(const CurvilinearGridNodeIndices& snappedNodeIndex,
                                                           const CurvilinearGridNodeIndices& gridLinePointIndex) const
{
    auto [scalingH, scalingV, scalingMixed] = CurvilinearGrid::ComputeDirectionalSmoothingFactors(gridLinePointIndex, snappedNodeIndex, m_indexBoxLowerLeft, m_indexBoxUpperRight);

    double factor = 0.0;

    if (m_smoothingRegionIndicator.m_m == 1 && m_smoothingRegionIndicator.m_n == 1)
    {
        factor = scalingMixed;
    }
    else if (m_smoothingRegionIndicator.m_n == 1)
    {
        factor = scalingV;
    }
    else if (m_smoothingRegionIndicator.m_m == 1)
    {
        factor = scalingH;
    }

    return factor;
}

meshkernel::NonDirectionalSmoothingCalculator::NonDirectionalSmoothingCalculator(const CurvilinearGrid& originalGrid,
                                                                                 const CurvilinearGrid& snappedGrid,
                                                                                 const LandBoundary& landBoundary) : m_originalGrid(originalGrid),
                                                                                                                     m_snappedGrid(snappedGrid),
                                                                                                                     m_smoothingRegionMinimum(CalculateSmoothingRegion(m_originalGrid.GetBoundingBox(),
                                                                                                                                                                       landBoundary.GetBoundingBox())) {}

meshkernel::NonDirectionalSmoothingCalculator::NonDirectionalSmoothingCalculator(const CurvilinearGrid& originalGrid,
                                                                                 const CurvilinearGrid& snappedGrid,
                                                                                 const Splines& spline) : m_originalGrid(originalGrid),
                                                                                                          m_snappedGrid(snappedGrid),
                                                                                                          m_smoothingRegionMinimum(CalculateSmoothingRegion(m_originalGrid.GetBoundingBox(),
                                                                                                                                                            spline.GetBoundingBox(0))) {}

double meshkernel::NonDirectionalSmoothingCalculator::CalculateSmoothingRegion(const BoundingBox& gridBoundingBox,
                                                                               const BoundingBox& entityBoundingBox)
{
    Vector delta = Merge(gridBoundingBox, entityBoundingBox).Delta();

    delta *= smoothingRegionEnlargementFactor;

    delta.y() = std::max(delta.y(), aspectRatio * delta.x());
    return delta.y() / (6.0 * aspectRatio);
}

double meshkernel::NonDirectionalSmoothingCalculator::compute(const CurvilinearGridNodeIndices& snappedNodeIndex,
                                                              const CurvilinearGridNodeIndices& gridLinePointIndex) const
{
    const Point meshDelta = m_snappedGrid.GetNode(snappedNodeIndex) - m_originalGrid.GetNode(snappedNodeIndex);
    const Point pointDelta = m_originalGrid.GetNode(gridLinePointIndex) - m_originalGrid.GetNode(snappedNodeIndex);
    double factor = 0.0;

    double rsx = std::max(m_smoothingRegionMinimum, std::hypot(meshDelta.x, meshDelta.y));

    if (double rn = std::hypot(pointDelta.x, pointDelta.y); rn < rsx)
    {
        factor = 0.5 * (1.0 + std::cos(std::numbers::pi * rn / rsx));
    }

    return factor;
}
