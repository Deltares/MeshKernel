//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2021.
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
#include <tuple>

#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/LandBoundary.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Polygon.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/TriangulationWrapper.hpp>

using meshkernel::Polygons;

Polygons::Polygons(const std::vector<Point>& polygon, Projection projection) : m_nodes(polygon), m_projection(projection)
{

    std::cout << "Polygons::Polygons: Number of points: " << polygon.size() << std::endl;

    // Find the polygons in the current list of points
    m_outer_polygons_indices = FindIndices(polygon, 0, static_cast<UInt>(polygon.size()), constants::missing::doubleValue);

    for (UInt i = 0; i < m_outer_polygons_indices.size(); ++i)
    {
        const auto& [outer_start, outer_end] = m_outer_polygons_indices[i];

        // The inner polygon indices, the first interval corresponds to the outer polygon
        const auto inner_polygons_indices = FindIndices(polygon, outer_start, outer_end, constants::missing::innerOuterSeparator);

        std::vector<Point> polygonPoints;
        polygonPoints.reserve(outer_end - outer_start + 1);

        for (size_t i = outer_start; i <= outer_end; ++i)
        {
            polygonPoints.emplace_back(m_nodes[i]);
        }

        std::cout << "Constructing polygon " << std::boolalpha << (inner_polygons_indices.size() == 1) << "  "
                  << polygonPoints.size() << "  " << outer_start << "  " << outer_end << "  " << inner_polygons_indices.size() << "  "
                  << m_nodes.size()
                  << std::endl;

        m_enclosures.emplace_back(PolygonalEnclosure(std::move(polygonPoints), m_projection));

        // No inner polygon found
        if (inner_polygons_indices.size() <= 1)
        {
            continue;
        }

        std::cout << " inner indices: "
                  << inner_polygons_indices[0].first << "  " << inner_polygons_indices[0].second << "  "
                  << inner_polygons_indices[1].first << "  " << inner_polygons_indices[1].second << "  "
                  << std::endl;

        // The first inner
        const auto inner_start = inner_polygons_indices[1].first;

        // store inner polygons for this outer polygon
        auto inner_polygons = std::vector<std::pair<UInt, UInt>>{};

        for (UInt j = 1; j < inner_polygons_indices.size(); ++j)
        {
            inner_polygons.emplace_back(inner_polygons_indices[j]);
        }

        // shift the index of the outer polygon, the
        m_outer_polygons_indices[i].second = inner_start - 2;
    }
}

// TOOD what to do with this?
size_t Polygons::GetNumNodes() const
{
    return m_nodes.size();

#if 0
    size_t nodeCount = 0;

    for (size_t i = 0; i < m_enclosures.size(); ++i)
    {
        nodeCount += m_enclosures[i].Outer().Size();
    }

    return nodeCount;
#endif
}

std::vector<std::vector<meshkernel::Point>> Polygons::ComputePointsInPolygons() const
{
    std::vector<std::vector<Point>> generatedPoints(GetNumPolygons(), std::vector<Point>());
    std::vector<Point> localPolygon(GetNumNodes());
    TriangulationWrapper triangulationWrapper;

    for (UInt polygonIndex = 0; polygonIndex < m_enclosures.size(); ++polygonIndex)
    {
        const PolygonalEnclosure& enclosure = m_enclosures[polygonIndex];
        const Polygon& polygon = enclosure.Outer();

        const auto [localPolygonArea, centerOfMass, isCounterClockWise] = polygon.FaceAreaAndCenterOfMass();

        // average triangle size
        const auto averageEdgeLength = polygon.ClosedPerimeterLength() / static_cast<double>(polygon.Size());
        const double averageTriangleArea = 0.25 * std::numbers::sqrt3 * averageEdgeLength * averageEdgeLength;

        // estimated number of triangles
        constexpr UInt SafetySize = 11;
        const auto numberOfTriangles = static_cast<UInt>(SafetySize * localPolygonArea / averageTriangleArea);

        if (numberOfTriangles == 0)
        {
            throw AlgorithmError("Polygons::ComputePointsInPolygons: The number of triangles = 0.");
        }

        triangulationWrapper.Compute(polygon.Nodes(),
                                     TriangulationWrapper::TriangulationOptions::GeneratePoints,
                                     averageTriangleArea,
                                     numberOfTriangles);
        generatedPoints[polygonIndex] = triangulationWrapper.SelectNodes(enclosure);
    }

    return generatedPoints;
}

std::vector<meshkernel::Point> Polygons::RefineFirstPolygon(UInt startIndex,
                                                            UInt endIndex,
                                                            double refinementDistance) const
{

    if (IsEmpty())
    {
        throw ConstraintError("Polygons::RefineFirstPolygon: No nodes in polygon.");
    }

    size_t polygonIndex = 0;
    size_t polygonStartNode = 0;
    size_t polygonEndNode = 0;

    if (startIndex == 0 && endIndex == 0)
    {
        const auto& [outerStart, outerEnd] = m_outer_polygons_indices[0];
        startIndex = outerStart;
        endIndex = outerEnd;
        polygonIndex = 0;
        polygonStartNode = outerStart;
        polygonEndNode = outerEnd;
    }

    if (endIndex <= startIndex)
    {
        throw ConstraintError(VariadicErrorMessage("Polygons::RefineFirstPolygon: The end index is smaller than the start index: {} >= {}.", startIndex, endIndex));
    }

    std::tie(polygonIndex, polygonStartNode, polygonEndNode) = PolygonIndex(startIndex, endIndex);

    return RefinePolygon(polygonIndex, polygonStartNode, polygonEndNode, refinementDistance);
}

std::vector<meshkernel::Point> Polygons::RefinePolygon(UInt polygonIndex, UInt startIndex, UInt endIndex, double refinementDistance) const
{

    if (polygonIndex >= m_enclosures.size())
    {
        throw ConstraintError(VariadicErrorMessage("Invalid polygon index: {} > {}.", polygonIndex, m_enclosures.size() - 1));
    }

    // TODO train wreck-ish.
    return m_enclosures[polygonIndex].Outer().Refine(startIndex, endIndex, refinementDistance);
}

meshkernel::Polygons Polygons::OffsetCopy(double distance, bool innerAndOuter) const
{
    UInt totalNumberOfPoints = GetNumNodes();
    // TODO Should the invalid points between enclosures be added too.
    // TODO use the correct separators between enclosures and inner polygons.
    std::vector<Point> newPolygonPoints(totalNumberOfPoints + (innerAndOuter ? totalNumberOfPoints + 1 : 0), Point());

    UInt innerCount = 0;
    UInt outerCount = totalNumberOfPoints;

    for (const PolygonalEnclosure& enclosure : m_enclosures)
    {
        std::vector<Point> outerOffsetPoints(enclosure.Outer().ComputeOffset(distance, innerAndOuter));

        for (size_t i = 0; i < outerOffsetPoints.size(); ++i)
        {
            newPolygonPoints[innerCount] = outerOffsetPoints[i];
            ++innerCount;
        }

        for (size_t i = 0; i < enclosure.NumberOfInner(); ++i)
        {
            std::vector<Point> innerOffsetPoints(enclosure.Inner(i).ComputeOffset(distance, innerAndOuter));
            // Probably add inner separator here (if not the last)

            for (size_t j = 0; j < outerOffsetPoints.size(); ++j)
            {
                newPolygonPoints[outerCount] = outerOffsetPoints[j];
                ++outerCount;
            }
        }

        // Probably add outer separator here (if not the last)
    }

    // set the new polygon
    Polygons newPolygon{newPolygonPoints, m_projection};
    return newPolygon;
}

void Polygons::SnapToLandBoundary(const LandBoundary& landBoundary [[maybe_unused]], UInt startIndex, UInt endIndex)
{
    if (IsEmpty())
    {
        throw ConstraintError(VariadicErrorMessage("No enclosures."));
    }

    if (startIndex == 0 && endIndex == 0)
    {
        std::cout << "################################" << std::endl;
        // TODO Does this mean all enclosures?
        endIndex = static_cast<UInt>(m_nodes.size()) - 1;
    }

    // TODO is it valid to snap a single point to the land boundary?
    if (startIndex >= endIndex)
    {
        throw ConstraintError(VariadicErrorMessage("The start index is greater than the end index: {} >= {}.", startIndex, endIndex));
    }

    const auto [polygonIndex, polygonStartNode, polygonEndNode] = PolygonIndex(startIndex, endIndex);
    std::cout << " indices: " << polygonIndex << "  " << polygonStartNode << "  " << polygonEndNode << "  " << std::endl;
    m_enclosures[polygonIndex].SnapToLandBoundary(polygonStartNode, polygonEndNode, landBoundary);
}

std::tuple<meshkernel::UInt, meshkernel::UInt, meshkernel::UInt> Polygons::PolygonIndex(UInt startIndex, UInt endIndex) const
{
    if (IsEmpty())
    {
        throw ConstraintError("No enclosures.");
    }

    bool indicesAreValid = false;
    UInt polygonIndex = constants::missing::uintValue;
    UInt polygonStartNode = constants::missing::uintValue;
    UInt polygonEndNode = constants::missing::uintValue;

    if (startIndex == 0 && endIndex == 0)
    {
        const auto& [outerStart, outerEnd] = m_outer_polygons_indices[0];
        startIndex = outerStart;
        endIndex = outerEnd;
        polygonIndex = 0;
        polygonStartNode = outerStart;
        polygonEndNode = outerEnd;
    }

    for (UInt i = 0; i < GetNumPolygons(); ++i)
    {
        const auto& [outerStart, outerEnd] = m_outer_polygons_indices[i];

        if (startIndex >= outerStart && endIndex <= outerEnd)
        {
            indicesAreValid = true;
            polygonIndex = i;
            polygonStartNode = startIndex - outerStart;
            polygonEndNode = endIndex - outerStart;
            break;
        }
    }

    if (!indicesAreValid)
    {
        throw ConstraintError(VariadicErrorMessage("The indices are not valid: {}, {}.", startIndex, endIndex));
    }

    return {polygonIndex, polygonStartNode, polygonEndNode};
}

bool Polygons::IsPointInPolygon(Point const& point, UInt polygonIndex) const
{
    if (IsEmpty())
    {
        return true;
    }

    if (polygonIndex >= GetNumPolygons())
    {
        throw std::invalid_argument("Polygons::IsPointInPolygon: Invalid polygon index.");
    }

    return m_enclosures[polygonIndex].Contains(point);
}

meshkernel::UInt Polygons::GetNumPolygons() const
{
    return static_cast<UInt>(m_enclosures.size());
}

std::tuple<bool, meshkernel::UInt> Polygons::IsPointInPolygons(const Point& point) const
{
    // empty polygon means everything is included
    if (IsEmpty())
    {
        return {true, constants::missing::uintValue};
    }

    for (UInt i = 0; i < m_enclosures.size(); ++i)
    {
        const PolygonalEnclosure& enclosure = m_enclosures[i];

        if (enclosure.ContainsRegion(point) == PolygonalEnclosure::Region::Exterior)
        {
            return {true, i};
        }
        else if (enclosure.ContainsRegion(point) == PolygonalEnclosure::Region::Interior)
        {
            // Point can be found in an hole in the polygon
            break;
        }
    }

    return {false, constants::missing::uintValue};
}

std::vector<bool> Polygons::PointsInPolygons(const std::vector<Point>& points) const
{
    std::vector<bool> result(points.size(), false);

    // TODO if possible improve performance of polygon.Contains, perhaps with
    // multiple points in a single call.
    // Then this loop has to be changed.
    for (UInt i = 0; i < points.size(); ++i)
    {
        const auto [isInPolygon, polygonIndex] = IsPointInPolygons(points[i]);
        result[i] = isInPolygon;
    }

    return result;
}

// TODO put in header.
bool Polygons::IsEmpty() const
{
    return m_enclosures.empty();
}

meshkernel::BoundingBox Polygons::GetBoundingBox(UInt polygonIndex) const
{
    if (IsEmpty())
    {
        throw ConstraintError(VariadicErrorMessage("Enclosures list is empty."));
    }

    if (polygonIndex >= m_enclosures.size())
    {
        throw ConstraintError(VariadicErrorMessage("Invalid enclosure index: {}, maximum index: {}",
                                                   polygonIndex, m_enclosures.size() - 1));
    }

    return m_enclosures[polygonIndex].Outer().GetBoundingBox();
}
