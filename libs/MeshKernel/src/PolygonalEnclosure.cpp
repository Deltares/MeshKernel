//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2025.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation version 3.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// contact: delft3d.support@deltares.nl
// Stichting Deltares
// P.O. Box 177
// 2600 MH Delft, The Netherlands
//
// All indications and logos of, and references to, "Delft3D" and "Deltares"
// are registered trademarks of Stichting Deltares, and remain the property of
// Stichting Deltares. All rights reserved.
//
//------------------------------------------------------------------------------

#include <numbers>

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/LandBoundary.hpp"
#include "MeshKernel/Operations.hpp"
#include "MeshKernel/PolygonalEnclosure.hpp"
#include "MeshKernel/TriangulationWrapper.hpp"

meshkernel::PolygonalEnclosure::PolygonalEnclosure(const std::vector<Point>& points,
                                                   Projection projection)
{
    // The inner polygon indices, the first interval corresponds to the outer polygon
    const IndexRangeArray innerIndices = FindIndices(points, 0, points.size() - 1, constants::missing::innerOuterSeparator);

    ConstructOuterPolygon(points, 0, points.size() - 1, innerIndices, projection);

    if (!innerIndices.empty())
    {
        ConstructInnerPolygons(points, innerIndices, projection);
    }
}

meshkernel::Polygon meshkernel::PolygonalEnclosure::ConstructPolygon(const std::vector<Point>& points,
                                                                     size_t start,
                                                                     size_t end,
                                                                     Projection projection)
{
    if (start > end)
    {
        throw ConstraintError("The start index is greater than the end index: {} > {}.", start, end);
    }

    if (end >= points.size())
    {
        throw ConstraintError("The end index is greater than the number of points: {} >= {}.", end, points.size());
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
                                                           size_t start,
                                                           size_t end,
                                                           const IndexRangeArray& innerIndices,
                                                           Projection projection)
{

    const size_t outerStartIndex = start;
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

meshkernel::UInt meshkernel::PolygonalEnclosure::GetNumberOfNodes() const
{
    UInt nodeCount = Outer().Size();

    for (size_t i = 0; i < m_inner.size(); ++i)
    {
        nodeCount += m_inner[i].Size();
    }

    return nodeCount;
}

bool meshkernel::PolygonalEnclosure::Contains(const Point& pnt) const
{
    // If the point is in one of the inner (island) polygons, then it is considered to be outside the enclosure.
    return ContainsRegion(pnt) == Region::Exterior;
}

meshkernel::PolygonalEnclosure::Region meshkernel::PolygonalEnclosure::ContainsRegion(const Point& pnt) const
{
    Region region = Region::None;

    if (pnt.IsValid() && m_outer.Contains(pnt))
    {
        // The point is contained within the perimeter of the enclosure.
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
        throw ConstraintError("The end index is greater than the number of points in the outer polygon: {} >= {}.",
                              endIndex,
                              m_outer.Size());
    }

    m_outer.SnapToLandBoundary(startIndex, endIndex, landBoundary);
}

std::vector<meshkernel::Point> meshkernel::PolygonalEnclosure::Refine(UInt startIndex, UInt endIndex, double refinementDistance) const
{
    if (endIndex >= m_outer.Size())
    {
        throw ConstraintError("The end index is greater than the number of points in the outer polygon: {} >= {}.",
                              endIndex,
                              m_outer.Size());
    }

    return m_outer.Refine(startIndex, endIndex, refinementDistance);
}

std::vector<meshkernel::Point> meshkernel::PolygonalEnclosure::LinearRefine(UInt startIndex, UInt endIndex) const
{
    if (endIndex >= m_outer.Size())
    {
        throw ConstraintError("The end index is greater than the number of points in the outer polygon: {} >= {}.",
                              endIndex,
                              m_outer.Size());
    }

    return m_outer.LinearRefine(startIndex, endIndex);
}

void meshkernel::PolygonalEnclosure::CopyPoints(const std::vector<Point>& source,
                                                const size_t start,
                                                const size_t end,
                                                UInt& count,
                                                std::vector<Point>& target)
{

    for (size_t i = start; i < end; ++i)
    {
        target[count] = source[i];
        ++count;
    }
}

std::tuple<std::unique_ptr<meshkernel::PolygonalEnclosure>, std::unique_ptr<meshkernel::PolygonalEnclosure>>
meshkernel::PolygonalEnclosure::OffsetCopy(const double distance, const bool outwardsAndInwards) const
{
    std::vector<Point> outerOffsetPoints(GetNumberOfNodes() + NumberOfInner(), Point());
    std::vector<Point> innerOffsetPoints;

    // Get offset for the outer perimeter polygon

    std::vector<Point> outerOffsetPolygon(Outer().ComputeOffset(distance, outwardsAndInwards));

    UInt outerCount = 0;
    UInt innerCount = 0;

    CopyPoints(outerOffsetPolygon, 0, Outer().Size(), outerCount, outerOffsetPoints);

    if (outwardsAndInwards)
    {
        innerOffsetPoints.resize(outerOffsetPoints.size(), Point());
        CopyPoints(outerOffsetPolygon, Outer().Size() + 1, 2 * Outer().Size() + 1, innerCount, innerOffsetPoints);
    }

    // Now compute offset for all inner polygons
    for (size_t i = 0; i < m_inner.size(); ++i)
    {
        const Polygon& innerPolygon = Inner(i);

        std::vector<Point> innerOffsetPolygon(innerPolygon.ComputeOffset(distance, outwardsAndInwards));

        outerOffsetPoints[outerCount] = Point(constants::missing::innerOuterSeparator, constants::missing::innerOuterSeparator);
        ++outerCount;

        CopyPoints(innerOffsetPolygon, 0, innerPolygon.Size(), outerCount, outerOffsetPoints);

        if (outwardsAndInwards)
        {
            innerOffsetPoints[innerCount] = Point(constants::missing::innerOuterSeparator, constants::missing::innerOuterSeparator);
            ++innerCount;

            CopyPoints(innerOffsetPolygon, innerPolygon.Size() + 1, 2 * innerPolygon.Size() + 1, innerCount, innerOffsetPoints);
        }
    }

    std::unique_ptr<PolygonalEnclosure> outwardOffset(std::make_unique<PolygonalEnclosure>(outerOffsetPoints, Outer().GetProjection()));
    std::unique_ptr<PolygonalEnclosure> inwardOffset;

    if (outwardsAndInwards)
    {
        inwardOffset = std::make_unique<PolygonalEnclosure>(innerOffsetPoints, Outer().GetProjection());
    }

    return {std::move(outwardOffset), std::move(inwardOffset)};
}

std::vector<meshkernel::Point> meshkernel::PolygonalEnclosure::GeneratePoints(const double scaleFactor) const
{
    TriangulationWrapper triangulationWrapper;

    const auto [localPolygonArea, centerOfMass, direction] = m_outer.FaceAreaAndCenterOfMass();

    // average triangle size
    const auto averageEdgeLength = m_outer.PerimeterLength() / static_cast<double>(m_outer.Size());

    double averageTriangleArea = 0.25 * std::numbers::sqrt3 * averageEdgeLength * averageEdgeLength;

    // estimated number of triangles
    constexpr UInt SafetySize = 11;
    const auto numberOfTriangles = static_cast<UInt>(SafetySize * std::max(1.0, localPolygonArea / averageTriangleArea));

    if (numberOfTriangles == 0)
    {
        throw AlgorithmError("The number of triangles = 0.");
    }

    if (scaleFactor == constants::missing::doubleValue)
    {
        const auto [minimumSegmentLength, maximumSegmentLength] = m_outer.SegmentLengthExtrema();
        const double segmentRatio = (minimumSegmentLength == constants::missing::doubleValue || minimumSegmentLength == 0.0 ? std::numbers::sqrt2 : maximumSegmentLength / minimumSegmentLength);
        averageTriangleArea *= 0.5 * segmentRatio * segmentRatio;
    }
    else
    {
        averageTriangleArea *= scaleFactor;
    }

    triangulationWrapper.Compute(m_outer.Nodes(),
                                 TriangulationWrapper::TriangulationOptions::GeneratePoints,
                                 averageTriangleArea,
                                 numberOfTriangles);

    return triangulationWrapper.SelectNodes(*this);
}
