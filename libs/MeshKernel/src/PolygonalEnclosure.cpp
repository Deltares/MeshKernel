#include "MeshKernel/PolygonalEnclosure.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/LandBoundary.hpp"
#include "MeshKernel/Operations.hpp"

meshkernel::PolygonalEnclosure::PolygonalEnclosure(const std::vector<Point>& points,
                                                   Projection projection) : PolygonalEnclosure(points, 0, points.size() - 1, projection) {}

meshkernel::PolygonalEnclosure::PolygonalEnclosure(const std::vector<Point>& points,
                                                   size_t start, size_t end,
                                                   Projection projection)
{
    double tolerance = 1.0e-10;

    // TODO check last point is equal to first

    std::cout << "Number of points: " << points.size() << "  " << std::boolalpha << IsEqual(points[0], points[points.size() - 1], tolerance) << std::endl;

    std::cout << points[0].x << "  " << points[0].y << "  " << points[points.size() - 1].x << "  " << points[points.size() - 1].y << std::endl;

    // if (points.size() < 4 || end - start + 1 < 4)
    // {
    //     throw ConstraintError(VariadicErrorMessage("Insufficient points"));
    // }

    // The inner polygon indices, the first interval corresponds to the outer polygon
    IndexRangeArray innerIndices = FindIndices(points, start, end, constants::missing::innerOuterSeparator);

    ConstructOuterPolygon(points, start, end, innerIndices, projection);

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
        // Error
    }

    if (end >= points.size())
    {
        // Error
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
        // Comment on why
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

    bool pointIsContained = m_outer.Contains(pnt);

    for (const Polygon& innerPolygon : m_inner)
    {
        if (innerPolygon.Contains(pnt))
        {
            // If the point is contained in any of the inner polygons then
            // it is considered outside the outer polygon.
            pointIsContained = false;
            break;
        }
    }

    return pointIsContained;
}

int meshkernel::PolygonalEnclosure::ContainsRegion(const Point& pnt) const
{
    int result = 0;

    if (m_outer.Contains(pnt))
    {
        result = 1;
    }

    for (const Polygon& innerPolygon : m_inner)
    {
        if (innerPolygon.Contains(pnt))
        {
            // If the point is contained in any of the inner polygons then
            // it is considered outside the outer polygon.
            result = 2;
            break;
        }
    }

    return result;
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
    m_outer.SnapToLandBoundary(startIndex, endIndex, landBoundary);
}

// meshkernel::PolygonalEnclosure meshkernel::PolygonalEnclosure::Displace(double displacement) const
// {

//     return PolygonalEnclosure();
// }
