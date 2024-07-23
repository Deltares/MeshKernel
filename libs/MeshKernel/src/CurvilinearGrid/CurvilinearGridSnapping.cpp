#include "MeshKernel/CurvilinearGrid/CurvilinearGridSnapping.hpp"

#include "MeshKernel/BoundingBox.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Vector.hpp"

#include <algorithm>

meshkernel::CurvilinearGridSnapping::CurvilinearGridSnapping(CurvilinearGrid& grid,
                                                             const std::vector<Point>& points) : m_grid(grid),
                                                                                                 m_originalGrid(grid),
                                                                                                 m_controlPoints(points)
{
    if (m_controlPoints.size() <= 1 || m_controlPoints.size() > 4)
    {
        throw ConstraintError("Snapping line or region has not been defined, number of points: {}", m_controlPoints.size());
    }

    if (const bool allValid = std::accumulate(m_controlPoints.begin(), m_controlPoints.end(), true, [](const bool val, const Point& p)
                                              { return val && p.IsValid(); });
        !allValid)
    {
        throw ConstraintError("1 or more of the selected points is invalid");
    }

    Initialise();
}

void meshkernel::CurvilinearGridSnapping::Initialise()
{
    std::tie(m_lineStartIndex, m_lineEndIndex) = m_grid.ComputeBlockFromCornerPoints(m_controlPoints[0], m_controlPoints[1]);

    if (m_controlPoints.size() == 2)
    {
        m_indexBoxLowerLeft = m_lineStartIndex;
        m_indexBoxUpperRight = m_lineEndIndex;
    }
    else
    {

        if (m_controlPoints.size() == 3)
        {

            auto const extentIndexPosition = m_grid.FindLocationIndex(m_controlPoints[2], Location::Nodes);
            CurvilinearGridNodeIndices extentIndex = m_grid.GetCurvilinearGridNodeIndices(extentIndexPosition);

            m_indexBoxLowerLeft = CurvilinearGridNodeIndices(std::min({m_lineStartIndex.m_n, m_lineEndIndex.m_n, extentIndex.m_n}),
                                                             std::min({m_lineStartIndex.m_m, m_lineEndIndex.m_m, extentIndex.m_m}));
            m_indexBoxUpperRight = CurvilinearGridNodeIndices(std::max({m_lineStartIndex.m_n, m_lineEndIndex.m_n, extentIndex.m_n}),
                                                              std::max({m_lineStartIndex.m_m, m_lineEndIndex.m_m, extentIndex.m_m}));
        }
        else
        {
            auto [lowerExtentIndex, upperExtentIndex] = m_grid.ComputeBlockFromCornerPoints(m_controlPoints[2], m_controlPoints[3]);

            m_indexBoxLowerLeft = CurvilinearGridNodeIndices(std::min(m_lineStartIndex.m_n, lowerExtentIndex.m_n), std::min(m_lineStartIndex.m_m, lowerExtentIndex.m_m));
            m_indexBoxUpperRight = CurvilinearGridNodeIndices(std::max(m_lineEndIndex.m_n, upperExtentIndex.m_n), std::max(m_lineEndIndex.m_m, upperExtentIndex.m_m));
        }
    }

    m_expansionRegionIndicator = CurvilinearGridNodeIndices(std::min<UInt>(1, m_lineEndIndex.m_m - m_lineStartIndex.m_m),
                                                            std::min<UInt>(1, m_lineEndIndex.m_n - m_lineStartIndex.m_n));
}

std::tuple<meshkernel::CurvilinearGridNodeIndices, meshkernel::CurvilinearGridNodeIndices>
meshkernel::CurvilinearGridSnapping::ComputeLoopBounds(const CurvilinearGridNodeIndices& snappedNodeIndex) const
{

    if (m_controlPoints.size() == 2)
    {
        const auto m1 = static_cast<UInt>(std::max<int>(1, snappedNodeIndex.m_m + 1 - predefinedExpansionRegionFactor * m_expansionRegionIndicator.m_m) - 1);
        const auto m2 = static_cast<UInt>(std::min<int>(m_grid.NumM(), snappedNodeIndex.m_m + 1 + predefinedExpansionRegionFactor * m_expansionRegionIndicator.m_m) - 1);
        const auto n1 = static_cast<UInt>(std::max<int>(1, snappedNodeIndex.m_n + 1 - predefinedExpansionRegionFactor * m_expansionRegionIndicator.m_n) - 1);
        const auto n2 = static_cast<UInt>(std::min<int>(m_grid.NumN(), snappedNodeIndex.m_n + 1 + predefinedExpansionRegionFactor * m_expansionRegionIndicator.m_n) - 1);
        return {CurvilinearGridNodeIndices(n1, m1), CurvilinearGridNodeIndices(n2, m2)};
    }

    const auto n1 = static_cast<UInt>(std::max<int>(m_indexBoxLowerLeft.m_n, snappedNodeIndex.m_n + 1 - userDefinedExpansionRegionFactor * m_expansionRegionIndicator.m_n - 1));
    const auto m1 = static_cast<UInt>(std::max<int>(m_indexBoxLowerLeft.m_m, snappedNodeIndex.m_m + 1 - userDefinedExpansionRegionFactor * m_expansionRegionIndicator.m_m - 1));

    const auto n2 = static_cast<UInt>(std::min<int>(m_indexBoxUpperRight.m_n, snappedNodeIndex.m_n + 1 + userDefinedExpansionRegionFactor * m_expansionRegionIndicator.m_n - 1));
    const auto m2 = static_cast<UInt>(std::min<int>(m_indexBoxUpperRight.m_m, snappedNodeIndex.m_m + 1 + userDefinedExpansionRegionFactor * m_expansionRegionIndicator.m_m - 1));

    const auto first = CurvilinearGridNodeIndices(n1, m1);
    const auto second = CurvilinearGridNodeIndices(n2, m2);

    return {first, second};
}

std::unique_ptr<meshkernel::CurvilinearGridMeshExpansionCalculator> meshkernel::CurvilinearGridSnapping::AllocateCurvilinearGridMeshExpansionCalculator() const
{
    std::unique_ptr<CurvilinearGridMeshExpansionCalculator> expansionCalculator;

    if (m_controlPoints.size() > 2u)
    {
        // User defined expansion region
        expansionCalculator = std::make_unique<UserDefinedRegionExpasionCalculator>(m_indexBoxLowerLeft, m_indexBoxUpperRight, m_expansionRegionIndicator);
    }
    else
    {
        // Predefined expansion region from mesh entity, e.g. spline or land boundary
        expansionCalculator = AllocateCurvilinearGridMeshExpansionCalculator(m_originalGrid, m_grid);
    }

    return expansionCalculator;
}

void meshkernel::CurvilinearGridSnapping::ApplyExpansionToGrid(const CurvilinearGridNodeIndices& snappedNodeIndex,
                                                               const CurvilinearGridMeshExpansionCalculator& expansionFactor)
{
    auto [lowerBound, upperBound] = ComputeLoopBounds(snappedNodeIndex);

    const Point meshDelta = m_grid.GetNode(snappedNodeIndex) - m_originalGrid.GetNode(snappedNodeIndex);

    CurvilinearGridNodeIndices gridLinePointIndex;

    for (UInt i = lowerBound.m_m; i <= upperBound.m_m; ++i)
    {
        gridLinePointIndex.m_m = i;

        for (UInt j = lowerBound.m_n; j <= upperBound.m_n; ++j)
        {
            gridLinePointIndex.m_n = j;

            if (m_originalGrid.GetNode(gridLinePointIndex).IsValid())
            {
                const double factor = expansionFactor.compute(snappedNodeIndex, gridLinePointIndex);

                if (factor > 0.0)
                {
                    m_grid.GetNode(gridLinePointIndex) = m_originalGrid.GetNode(gridLinePointIndex) + factor * meshDelta;
                }
            }
        }
    }
}

meshkernel::UndoActionPtr meshkernel::CurvilinearGridSnapping::Compute()
{
    // probably can reduce storage required. m_lineStartIndex.m_m and m_lineEndIndex
    std::unique_ptr<CurvilinearGridBlockUndoAction> undoAction = CurvilinearGridBlockUndoAction::Create(m_grid, {0, 0}, {m_grid.NumN(), m_grid.NumM()});
    std::unique_ptr<CurvilinearGridMeshExpansionCalculator> expansionCalculator = AllocateCurvilinearGridMeshExpansionCalculator();

    CurvilinearGridNodeIndices snappedNodeIndex;

    for (UInt i = m_lineStartIndex.m_m; i <= m_lineEndIndex.m_m; ++i)
    {
        snappedNodeIndex.m_m = i;

        for (UInt j = m_lineStartIndex.m_n; j <= m_lineEndIndex.m_n; ++j)
        {
            snappedNodeIndex.m_n = j;

            Point currentPoint = m_originalGrid.GetNode(snappedNodeIndex);

            if (!currentPoint.IsValid())
            {
                continue;
            }

            m_grid.GetNode(snappedNodeIndex) = FindNearestPoint(currentPoint);

            // Only shift the line points in the grid line/region if the current grid point differs from the
            // grid point (at the same index) snapped to the boundary.
            if (!IsEqual(currentPoint, m_grid.GetNode(snappedNodeIndex), epsilon))
            {
                ApplyExpansionToGrid(snappedNodeIndex, *expansionCalculator);
            }
        }
    }

    return undoAction;
}
