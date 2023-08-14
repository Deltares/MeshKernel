#include "MeshKernel/CurvilinearGrid/CurvilinearGridSnapping.hpp"

#include "MeshKernel/BoundingBox.hpp"
#include "MeshKernel/Exceptions.hpp"

#include <cmath>
#include <numbers>

meshkernel::DirectionalSmoothingCalculator::DirectionalSmoothingCalculator(const CurvilinearGrid& grid,
                                                                           const CurvilinearGridNodeIndices& lowerLeft,
                                                                           const CurvilinearGridNodeIndices& upperRight,
                                                                           const CurvilinearGridNodeIndices& regionIndicator) : m_grid(grid),
                                                                                                                                m_indexBoxLowerLeft(lowerLeft),
                                                                                                                                m_indexBoxUpperRight(upperRight),
                                                                                                                                m_smoothingRegionIndicator(regionIndicator) {}

double meshkernel::DirectionalSmoothingCalculator::compute(const CurvilinearGridNodeIndices& snappedNodeIndex,
                                                           const CurvilinearGridNodeIndices& gridLinePointIndex) const
{
    auto [scalingH, scalingV, scalingMixed] = m_grid.ComputeDirectionalSmoothingFactors(gridLinePointIndex, snappedNodeIndex, m_indexBoxLowerLeft, m_indexBoxUpperRight);

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

double meshkernel::NonDirectionalSmoothingCalculator::CalculateSmoothingRegion(const CurvilinearGrid& grid,
                                                                               const LandBoundary& landBoundary)
{
    // Value from editgridlineblok.f90
    constexpr double aspectRatio = 990.0 / 1600.0;

    BoundingBox gridBb = grid.GetBoundingBox();
    BoundingBox landBoundaryBb = landBoundary.GetBoundingBox();

    Point delta = meshkernel::merge(gridBb, landBoundaryBb).delta();

    delta *= 1.2;

    delta.y = std::max(delta.y, aspectRatio * delta.x);
    return delta.y / (6.0 * aspectRatio);
}

meshkernel::NonDirectionalSmoothingCalculator::NonDirectionalSmoothingCalculator(const CurvilinearGrid& originalGrid,
                                                                                 const CurvilinearGrid& snappedGrid,
                                                                                 const LandBoundary& landBoundary) : m_originalGrid(originalGrid),
                                                                                                                     m_snappedGrid(snappedGrid),
                                                                                                                     m_smoothingRegionMinimum(CalculateSmoothingRegion(m_originalGrid, landBoundary)) {}

double meshkernel::NonDirectionalSmoothingCalculator::compute(const CurvilinearGridNodeIndices& snappedNodeIndex,
                                                              const CurvilinearGridNodeIndices& gridLinePointIndex) const
{
    Point meshDelta = m_snappedGrid.GetNode(snappedNodeIndex) - m_originalGrid.GetNode(snappedNodeIndex);
    Point pointDelta = m_originalGrid.GetNode(gridLinePointIndex) - m_originalGrid.GetNode(snappedNodeIndex);
    double factor = 0.0;

    double rsx = std::max(m_smoothingRegionMinimum, std::hypot(meshDelta.x, meshDelta.y));

    if (double rn = std::hypot(pointDelta.x, pointDelta.y); rn < rsx)
    {
        rn = 0.5 * std::numbers::pi * rn / rsx;
        factor = 0.5 * (1.0 + std::cos(2.0 * rn));
    }

    return factor;
}

meshkernel::CurvilinearGridSnapping::CurvilinearGridSnapping(std::shared_ptr<CurvilinearGrid> grid,
                                                             const LandBoundary& lb,
                                                             const std::vector<Point>& points) : CurvilinearGridAlgorithm(grid),
                                                                                                 m_originalGrid(*grid),
                                                                                                 m_landBoundary(lb),
                                                                                                 m_points(points)
{
    if (m_points.size() <= 1 || m_points.size() > 4)
    {
        throw ConstraintError(VariadicErrorMessage("Snapping line or region has not been defined, number of points: {}", m_points.size()));
    }

    if (bool allValid = std::accumulate(m_points.begin(), m_points.end(), true, [](bool val, const Point& p)
                                        { return val && p.IsValid(); });
        !allValid)
    {
        throw ConstraintError("1 or more of the selected points is invalid");
    }

    Initialise();
}

void meshkernel::CurvilinearGridSnapping::Initialise()
{
    std::tie(m_lineStartIndex, m_lineEndIndex) = m_grid.ComputeBlockFromCornerPoints(m_points[0], m_points[1]);

    if (m_points.size() == 2)
    {
        m_indexBoxLowerLeft = m_lineStartIndex;
        m_indexBoxUpperRight = m_lineEndIndex;

        // The intermediate result in the subtraction cannot be less than zero.
        m_smoothingRegionIndicator = CurvilinearGridNodeIndices(std::min<UInt>(1, m_lineEndIndex.m_n - m_lineStartIndex.m_n),
                                                                std::min<UInt>(1, m_lineEndIndex.m_m - m_lineStartIndex.m_m));
    }
    else
    {

        if (m_points.size() == 3)
        {
            CurvilinearGridNodeIndices extentIndex = m_grid.GetNodeIndices(m_points[2]);

            m_indexBoxLowerLeft = CurvilinearGridNodeIndices(std::min({m_lineStartIndex.m_m, m_lineEndIndex.m_m, extentIndex.m_m}),
                                                             std::min({m_lineStartIndex.m_n, m_lineEndIndex.m_n, extentIndex.m_n}));
            m_indexBoxUpperRight = CurvilinearGridNodeIndices(std::max({m_lineStartIndex.m_m, m_lineEndIndex.m_m, extentIndex.m_m}),
                                                              std::max({m_lineStartIndex.m_n, m_lineEndIndex.m_n, extentIndex.m_n}));
        }
        else
        {
            auto [lowerExtentIndex, upperExtentIndex] = m_grid.ComputeBlockFromCornerPoints(m_points[2], m_points[3]);

            m_indexBoxLowerLeft = CurvilinearGridNodeIndices(std::min(m_lineStartIndex.m_m, lowerExtentIndex.m_m), std::min(m_lineStartIndex.m_n, lowerExtentIndex.m_n));
            m_indexBoxUpperRight = CurvilinearGridNodeIndices(std::max(m_lineEndIndex.m_m, upperExtentIndex.m_m), std::max(m_lineEndIndex.m_n, upperExtentIndex.m_n));
        }
    }

    m_smoothingRegionIndicator = CurvilinearGridNodeIndices(std::min<UInt>(1, m_lineEndIndex.m_n - m_lineStartIndex.m_n),
                                                            std::min<UInt>(1, m_lineEndIndex.m_m - m_lineStartIndex.m_m));
}

std::tuple<meshkernel::CurvilinearGridNodeIndices, meshkernel::CurvilinearGridNodeIndices>
meshkernel::CurvilinearGridSnapping::ComputeLoopBounds(const CurvilinearGridNodeIndices& snappedNodeIndex) const
{

    if (m_points.size() == 2)
    {
        // Name from modgr4.f90
        constexpr UInt nump = 80;
        auto m1 = static_cast<UInt>(std::max<int>(1, snappedNodeIndex.m_m + 1 - nump * m_smoothingRegionIndicator.m_m) - 1);
        auto m2 = static_cast<UInt>(std::min<int>(m_grid.m_numM, snappedNodeIndex.m_m + 1 + nump * m_smoothingRegionIndicator.m_m) - 1);
        auto n1 = static_cast<UInt>(std::max<int>(1, snappedNodeIndex.m_n + 1 - nump * m_smoothingRegionIndicator.m_n) - 1);
        auto n2 = static_cast<UInt>(std::min<int>(m_grid.m_numN, snappedNodeIndex.m_n + 1 + nump * m_smoothingRegionIndicator.m_n) - 1);
        return {CurvilinearGridNodeIndices(m1, n1), CurvilinearGridNodeIndices(m2, n2)};
    }
    else
    {
        auto m1 = static_cast<UInt>(std::max<int>(m_indexBoxLowerLeft.m_m, snappedNodeIndex.m_m + 1 - 10000 * m_smoothingRegionIndicator.m_m - 1));
        auto m2 = static_cast<UInt>(std::min<int>(m_indexBoxUpperRight.m_m, snappedNodeIndex.m_m + 1 + 10000 * m_smoothingRegionIndicator.m_m - 1));
        auto n1 = static_cast<UInt>(std::max<int>(m_indexBoxLowerLeft.m_n, snappedNodeIndex.m_n + 1 - 10000 * m_smoothingRegionIndicator.m_n - 1));
        auto n2 = static_cast<UInt>(std::min<int>(m_indexBoxUpperRight.m_n, snappedNodeIndex.m_n + 1 + 10000 * m_smoothingRegionIndicator.m_n - 1));
        return {CurvilinearGridNodeIndices(m1, n1), CurvilinearGridNodeIndices(m2, n2)};
    }
}

void meshkernel::CurvilinearGridSnapping::ApplySmoothingToGrid(const CurvilinearGridNodeIndices& snappedNodeIndex,
                                                               const MeshSmoothingCalculator& smoothingFactor)
{
    auto [lowerBound, upperBound] = ComputeLoopBounds(snappedNodeIndex);

    Point meshDelta = m_grid.GetNode(snappedNodeIndex) - m_originalGrid.GetNode(snappedNodeIndex);

    CurvilinearGridNodeIndices gridLinePointIndex;

    for (UInt i = lowerBound.m_m; i <= upperBound.m_m; ++i)
    {
        gridLinePointIndex.m_m = i;

        for (UInt j = lowerBound.m_n; j <= upperBound.m_n; ++j)
        {
            gridLinePointIndex.m_n = j;

            if (m_originalGrid.GetNode(gridLinePointIndex).IsValid())
            {
                double factor = smoothingFactor.compute(snappedNodeIndex, gridLinePointIndex);

                if (factor > 0.0)
                {
                    m_grid.GetNode(gridLinePointIndex) = m_originalGrid.GetNode(gridLinePointIndex) + factor * meshDelta;
                }
            }
        }
    }
}

meshkernel::CurvilinearGrid meshkernel::CurvilinearGridSnapping::Compute()
{
    constexpr double epsilon = 1.0e-5;

    std::unique_ptr<MeshSmoothingCalculator> smoothingFactorCalculator;

    // Move to allocate function?
    if (m_points.size() > 2)
    {
        smoothingFactorCalculator = std::make_unique<DirectionalSmoothingCalculator>(m_originalGrid, m_indexBoxLowerLeft, m_indexBoxUpperRight, m_smoothingRegionIndicator);
    }
    else
    {
        smoothingFactorCalculator = std::make_unique<NonDirectionalSmoothingCalculator>(m_originalGrid, m_grid, m_landBoundary);
    }

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

            m_grid.GetNode(snappedNodeIndex) = m_landBoundary.FindNearestPoint(currentPoint, m_grid.m_projection);

            // Only shift the line points in the grid line/region if the current grid point differs from the
            // grid point (at the same index) snapped to the boundary.
            if (!IsEqual(currentPoint, m_grid.GetNode(snappedNodeIndex), epsilon))
            {
                ApplySmoothingToGrid(snappedNodeIndex, *smoothingFactorCalculator);
            }
        }
    }

    return m_grid;
}
