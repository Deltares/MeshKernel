#include "MeshKernel/CurvilinearGrid/CurvilinearGridSnapping.hpp"

#include "MeshKernel/BoundingBox.hpp"
#include "MeshKernel/Exceptions.hpp"

#include <cmath>
#include <numbers>

meshkernel::SmeerFunctie::SmeerFunctie(const CurvilinearGrid& grid,
                                       const CurvilinearGridNodeIndices& lowerLeft,
                                       const CurvilinearGridNodeIndices& upperRight,
                                       const CurvilinearGridNodeIndices& regionSize) : m_grid(grid),
                                                                                       m_indexBoxLowerLeft(lowerLeft),
                                                                                       m_indexBoxUpperRight(upperRight),
                                                                                       m_snappedRegionSize(regionSize) {}

double meshkernel::SmeerFunctie::compute(const CurvilinearGridNodeIndices& currentPointIndex,
                                         const CurvilinearGridNodeIndices& gridLineIndex) const
{
    auto [scalingH, scalingV, scalingMixed] = m_grid.ComputeDirectionalSmoothingFactors(gridLineIndex, currentPointIndex, m_indexBoxLowerLeft, m_indexBoxUpperRight);

    double factor = 0.0;

    if (m_snappedRegionSize.m_m == 1 && m_snappedRegionSize.m_n == 1)
    {
        factor = scalingMixed;
    }
    else if (m_snappedRegionSize.m_n == 1)
    {
        factor = scalingV;
    }
    else if (m_snappedRegionSize.m_m == 1)
    {
        factor = scalingH;
    }

    return factor;
}

meshkernel::NotSmeerFunctie::NotSmeerFunctie(const CurvilinearGrid& originalGrid,
                                             const CurvilinearGrid& snappedGrid,
                                             const LandBoundary& landBoundary) : m_originalGrid(originalGrid), m_snappedGrid(snappedGrid)
{
    // Value from editgridlineblok.f90
    constexpr double aspectRatio = 990.0 / 1600.0;

    BoundingBox gridBb = m_originalGrid.GetBoundingBox();
    BoundingBox lbBb = landBoundary.GetBoundingBox();

    Point delta = meshkernel::merge(gridBb, lbBb).delta();

    delta *= 1.2;

    if (delta.y < aspectRatio * delta.x)
    {
        delta.y = aspectRatio * delta.x;
    }

    m_dsix = delta.y / (6.0 * aspectRatio);
}

double meshkernel::NotSmeerFunctie::compute(const CurvilinearGridNodeIndices& currentPointIndex,
                                            const CurvilinearGridNodeIndices& gridLineIndex) const
{
    Point meshDelta = m_snappedGrid.GetNode(currentPointIndex) - m_originalGrid.GetNode(currentPointIndex);
    Point pointDelta = m_originalGrid.GetNode(gridLineIndex) - m_originalGrid.GetNode(currentPointIndex);
    double factor = 0.0;

    double rsx = std::max(m_dsix, std::hypot(meshDelta.x, meshDelta.y));
    double rn = std::hypot(pointDelta.x, pointDelta.y);

    if (rn < rsx)
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

    // Check that all the points are valid.
    struct AllValid
    {
        void operator()(const Point& p) { allValid = allValid && p.IsValid(); }
        operator bool() const { return allValid; }
        bool allValid{true};
    };

    AllValid allValid = std::for_each(m_points.begin(), m_points.end(), AllValid());

    if (!allValid)
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
        m_snappedRegionSize = CurvilinearGridNodeIndices(std::min<UInt>(1, m_lineEndIndex.m_n - m_lineStartIndex.m_n),
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
            m_indexBoxUpperRight = CurvilinearGridNodeIndices(std::max(m_lineStartIndex.m_m, upperExtentIndex.m_m), std::max(m_lineStartIndex.m_n, upperExtentIndex.m_n));
        }
    }

    m_snappedRegionSize = CurvilinearGridNodeIndices(std::min<UInt>(1, m_lineEndIndex.m_n - m_lineStartIndex.m_n),
                                                     std::min<UInt>(1, m_lineEndIndex.m_m - m_lineStartIndex.m_m));
}

std::tuple<meshkernel::CurvilinearGridNodeIndices, meshkernel::CurvilinearGridNodeIndices>
meshkernel::CurvilinearGridSnapping::ComputeLoopBounds(const CurvilinearGridNodeIndices& currentNodeIndex, const CurvilinearGridNodeIndices& snappedRegionSize) const
{

    if (m_points.size() == 2)
    {
        // TODO Add comment
        // TODO Why is this variable called nump?
        constexpr UInt nump = 80;
        UInt m1 = static_cast<UInt>(std::max<int>(1, currentNodeIndex.m_m + 1 - nump * snappedRegionSize.m_m) - 1);
        UInt m2 = static_cast<UInt>(std::min<int>(m_grid.m_numM, currentNodeIndex.m_m + 1 + nump * snappedRegionSize.m_m) - 1);
        UInt n1 = static_cast<UInt>(std::max<int>(1, currentNodeIndex.m_n + 1 - nump * snappedRegionSize.m_n) - 1);
        UInt n2 = static_cast<UInt>(std::min<int>(m_grid.m_numN, currentNodeIndex.m_n + 1 + nump * snappedRegionSize.m_n) - 1);
        return {CurvilinearGridNodeIndices(m1, n1), CurvilinearGridNodeIndices(m2, n2)};
    }
    else
    {
        UInt m1 = static_cast<UInt>(std::max<int>(m_indexBoxLowerLeft.m_m, currentNodeIndex.m_m + 1 - 10000 * snappedRegionSize.m_m - 1));
        UInt m2 = static_cast<UInt>(std::min<int>(m_indexBoxUpperRight.m_m, currentNodeIndex.m_m + 1 + 10000 * snappedRegionSize.m_m - 1));
        UInt n1 = static_cast<UInt>(std::max<int>(m_indexBoxLowerLeft.m_n, currentNodeIndex.m_n + 1 - 10000 * snappedRegionSize.m_n - 1));
        UInt n2 = static_cast<UInt>(std::min<int>(m_indexBoxUpperRight.m_n, currentNodeIndex.m_n + 1 + 10000 * snappedRegionSize.m_n - 1));
        return {CurvilinearGridNodeIndices(m1, n1), CurvilinearGridNodeIndices(m2, n2)};
    }
}

void meshkernel::CurvilinearGridSnapping::ModifyField(const CurvilinearGridNodeIndices& currentNodeIndex,
                                                      const CurvilinearGridNodeIndices& snappedRegionSize,
                                                      const MeshSnappingCalculator& translationFactor)
{
    auto [lowerBound, upperBound] = ComputeLoopBounds(currentNodeIndex, snappedRegionSize);

    Point meshDelta = m_grid.GetNode(currentNodeIndex) - m_originalGrid.GetNode(currentNodeIndex);

    CurvilinearGridNodeIndices gridLineIndex;

    for (UInt i = lowerBound.m_m; i <= upperBound.m_m; ++i)
    {
        gridLineIndex.m_m = i;

        for (UInt j = lowerBound.m_n; j <= upperBound.m_n; ++j)
        {
            gridLineIndex.m_n = j;

            if (m_originalGrid.GetNode(gridLineIndex).IsValid())
            {
                double factor = translationFactor.compute(currentNodeIndex, gridLineIndex);

                if (factor > 0.0)
                {
                    m_grid.GetNode(gridLineIndex) = m_originalGrid.GetNode(gridLineIndex) + factor * meshDelta;
                }
            }
        }
    }
}

meshkernel::CurvilinearGrid meshkernel::CurvilinearGridSnapping::Compute()
{
    constexpr double epsilon = 1.0e-5;

    CurvilinearGridNodeIndices currentNodeIndex;

    std::unique_ptr<MeshSnappingCalculator> translationFactorCalculator;

    // Move to allocate function?
    if (m_points.size() > 2)
    {
        translationFactorCalculator = std::make_unique<SmeerFunctie>(m_originalGrid, m_indexBoxLowerLeft, m_indexBoxUpperRight, m_snappedRegionSize);
    }
    else
    {
        translationFactorCalculator = std::make_unique<NotSmeerFunctie>(m_originalGrid, m_grid, m_landBoundary);
    }

    // TODO try to use range for loop with m_lineStartIndex
    for (UInt i = m_lineStartIndex.m_m; i <= m_lineEndIndex.m_m; ++i)
    {
        currentNodeIndex.m_m = i;

        for (UInt j = m_lineStartIndex.m_n; j <= m_lineEndIndex.m_n; ++j)
        {
            currentNodeIndex.m_n = j;

            Point currentPoint = m_originalGrid.GetNode(currentNodeIndex);

            // // TODO can we check here if the currentPoint is valid
            // if (!currentPoint.IsValid())
            // {
            //     continue;
            // }

            m_grid.GetNode(currentNodeIndex) = m_landBoundary.FindNearestPoint(currentPoint, m_grid.m_projection);

            // Only shift the line points in the grid line/region is the current grid point differs from the
            // same grid point snapped to the boundary.
            if (!IsEqual(currentPoint, m_grid.GetNode(currentNodeIndex), epsilon))
            {
                // TODO now its a member no need to pass parameter
                ModifyField(currentNodeIndex, m_snappedRegionSize, *translationFactorCalculator);
            }
        }
    }

    return m_grid;
}
