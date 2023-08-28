#include "MeshKernel/PolygonalEnclosure.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/LandBoundary.hpp"
#include "MeshKernel/Operations.hpp"

meshkernel::PolygonalEnclosure::PolygonalEnclosure(const std::vector<Point>& points,
                                                   Projection projection)
{
    // The inner polygon indices, the first interval corresponds to the outer polygon
    IndexRangeArray innerIndices = FindIndices(points, 0, points.size() - 1, constants::missing::innerOuterSeparator);

    ConstructOuterPolygon(points, 0, points.size() - 1, innerIndices, projection);

    if (innerIndices.size() >= 1)
    {
        ConstructInnerPolygons(points, innerIndices, projection);
    }
}

meshkernel::Polygon meshkernel::PolygonalEnclosure::ConstructPolygon(const std::vector<Point>& points,
                                                                     size_t start, size_t end,
                                                                     Projection projection)
{
    if (start > end)
    {
        throw ConstraintError(VariadicErrorMessage("The start index is greater than the end index: {} > {}.",
                                                   start, end));
    }

    if (end >= points.size())
    {
        throw ConstraintError(VariadicErrorMessage("The end index is greater than the number of points: {} >= {}.",
                                                   end, points.size()));
    }

    std::vector<Point> polygonPoints;
    polygonPoints.resize(end - start + 1);
    size_t innerI = 0;

    for (size_t i = start; i <= end; ++i)
    {
        polygonPoints[innerI] = points[i];
        ++innerI;
    }

    return Polygon(std::move(polygonPoints), projection);
}

void meshkernel::PolygonalEnclosure::ConstructOuterPolygon(const std::vector<Point>& points,
                                                           size_t start, size_t end,
                                                           const IndexRangeArray& innerIndices,
                                                           Projection projection)
{

    size_t outerStartIndex = start;
    size_t outerEndIndex = end;

    if (innerIndices.size() > 1)
    {
        outerEndIndex = innerIndices[0].second;
    }

    m_outer = ConstructPolygon(points, outerStartIndex, outerEndIndex, projection);
}

void meshkernel::PolygonalEnclosure::ConstructInnerPolygons(const std::vector<Point>& points,
                                                            const IndexRangeArray& innerIndices,
                                                            Projection projection)
{
    if (innerIndices.size() <= 1)
    {
        // Nothing to do
        return;
    }

    m_inner.reserve(innerIndices.size() - 1);

    for (size_t i = 1; i < innerIndices.size(); ++i)
    {
        const IndexRange& innerRange = innerIndices[i];
        m_inner.emplace_back(ConstructPolygon(points, innerRange.first, innerRange.second, projection));
    }
}

bool meshkernel::PolygonalEnclosure::Contains(const Point& pnt) const
{
    // If the point is in one of the inner (island) polygons, then it is considered to be outside the enclosure.
    return ContainsRegion(pnt) == Region::Exterior;
}

meshkernel::PolygonalEnclosure::Region meshkernel::PolygonalEnclosure::ContainsRegion(const Point& pnt) const
{
    Region region = Region::None;

    if (m_outer.Contains(pnt))
    {
        // The point is contained withing the perimeter of the enclosure.
        region = Region::Exterior;

        // Now check if the point is contained within any of the inner polygons.
        for (const Polygon& innerPolygon : m_inner)
        {
            if (innerPolygon.Contains(pnt))
            {
                // If the point is contained in any of the inner polygons then
                // it is considered outside the outer polygon enclosure.
                region = Region::Interior;
                break;
            }
        }
    }

    return region;
}

meshkernel::UInt meshkernel::PolygonalEnclosure::NumberOfPoints(const bool includeInterior) const
{
    UInt pointCount = m_outer.Size();

    if (includeInterior)
    {
        for (size_t i = 0; i < m_inner.size(); ++i)
        {
            pointCount += m_inner[i].Size();
        }
    }

    return pointCount;
}

void meshkernel::PolygonalEnclosure::SnapToLandBoundary(size_t startIndex, size_t endIndex, const LandBoundary& landBoundary)
{
    if (endIndex >= m_outer.Size())
    {
        throw ConstraintError(VariadicErrorMessage("The end index is greater than the number of points in the outer polygon: {} >= {}.",
                                                   endIndex, m_outer.Size()));
    }

    m_outer.SnapToLandBoundary(startIndex, endIndex, landBoundary);
}

std::vector<meshkernel::Point> meshkernel::PolygonalEnclosure::Refine(size_t startIndex, size_t endIndex, double refinementDistance)
{
    if (endIndex >= m_outer.Size())
    {
        throw ConstraintError(VariadicErrorMessage("The end index is greater than the number of points in the outer polygon: {} >= {}.",
                                                   endIndex, m_outer.Size()));
    }

    return m_outer.Refine(startIndex, endIndex, refinementDistance);
}
