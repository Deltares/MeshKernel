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
    // Find the polygons in the current list of points
    m_outer_polygons_indices = FindIndices(polygon, 0, static_cast<UInt>(polygon.size()), constants::missing::doubleValue);
    for (UInt i = 0; i < m_outer_polygons_indices.size(); ++i)
    {
        m_inner_polygons_indices[i] = std::vector<std::pair<UInt, UInt>>{};

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

        m_polygonGroups.emplace_back(PolygonalEnclosure(std::move(polygonPoints), m_projection));

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

        m_inner_polygons_indices[i] = inner_polygons;

        // shift the index of the outer polygon, the
        m_outer_polygons_indices[i].second = inner_start - 2;
    }
}

std::vector<std::vector<meshkernel::Point>> Polygons::ComputePointsInPolygons() const
{
    std::vector<std::vector<Point>> generatedPoints(GetNumPolygons(), std::vector<Point>());
    std::vector<Point> localPolygon(GetNumNodes());
    TriangulationWrapper triangulationWrapper;

    for (UInt polygonIndex = 0; polygonIndex < m_outer_polygons_indices.size(); ++polygonIndex)
    {
        const auto& [outerStart, outerEnd] = m_outer_polygons_indices[polygonIndex];

        localPolygon.clear();
        for (auto j = outerStart; j <= outerEnd; ++j)
        {
            localPolygon.emplace_back(m_nodes[j]);
        }

        // not a closed polygon
        const auto numLocalPoints = static_cast<UInt>(localPolygon.size());
        if (localPolygon[numLocalPoints - 1] != localPolygon[0] || localPolygon.size() < 4)
        {
            continue;
        }

        const PolygonalEnclosure& group = m_polygonGroups[polygonIndex];
        const Polygon& polygon = group.Outer();

        const auto [localPolygonArea, centerOfMass, isCounterClockWise] = polygon.FaceAreaAndCenterOfMass();
        const auto perimeter = polygon.ClosedPerimeterLength();

        // average triangle size
        const auto averageEdgeLength = perimeter / static_cast<double>(numLocalPoints);
        const double averageTriangleArea = 0.25 * constants::numeric::squareRootOfThree * averageEdgeLength * averageEdgeLength;

        // estimated number of triangles
        constexpr UInt SafetySize = 11;
        const auto numberOfTriangles = static_cast<UInt>(SafetySize * localPolygonArea / averageTriangleArea);

        if (numberOfTriangles == 0)
        {
            throw AlgorithmError("Polygons::ComputePointsInPolygons: The number of triangles is <= 0.");
        }

        triangulationWrapper.Compute(polygon.Points(),
                                     TriangulationWrapper::TriangulationOptions::GeneratePoints,
                                     averageTriangleArea,
                                     numberOfTriangles);

        generatedPoints[polygonIndex].reserve(triangulationWrapper.GetNumNodes());

        for (int i = 0; i < triangulationWrapper.GetNumNodes(); ++i)
        {
            if (Point p(triangulationWrapper.GetCoord(i)); group.Contains(p))
            {
                generatedPoints[polygonIndex].emplace_back(p);
            }
        }
    }

    return generatedPoints;
}

std::vector<meshkernel::Point> Polygons::RefineFirstPolygon(UInt startIndex,
                                                            UInt endIndex,
                                                            double refinementDistance) const
{

    if (m_outer_polygons_indices.empty())
    {
        throw ConstraintError("Polygons::RefineFirstPolygon: No nodes in polygon.");
    }

    if (startIndex == 0 && endIndex == 0)
    {
        const auto& [outerStart, outerEnd] = m_outer_polygons_indices[0];
        startIndex = outerStart;
        endIndex = outerEnd;
    }

    if (endIndex <= startIndex)
    {
        throw ConstraintError(VariadicErrorMessage("Polygons::RefineFirstPolygon: The end index is smaller than the start index: {} >= {}.", startIndex, endIndex));
    }

    //--------------------------------

    bool areIndicesValid = false;
    UInt polygonIndex;
    for (UInt i = 0; i < GetNumPolygons(); ++i)
    {
        const auto& [outerStart, outerEnd] = m_outer_polygons_indices[i];
        if (startIndex >= outerStart && endIndex <= outerEnd)
        {
            areIndicesValid = true;
            polygonIndex = i;
            break;
        }
    }

    if (!areIndicesValid)
    {
        throw ConstraintError(VariadicErrorMessage("Polygons::RefineFirstPolygon: The indices are not valid: {}, {}.", startIndex, endIndex));
    }

    //--------------------------------

    // TODO use Polygon

    const auto& [outerStart, outerEnd] = m_outer_polygons_indices[polygonIndex];

    const auto edgeLengths = PolygonEdgeLengths(m_nodes);
    std::vector<double> nodeLengthCoordinate(edgeLengths.size());
    nodeLengthCoordinate[0] = 0.0;
    for (UInt i = 1; i < edgeLengths.size(); ++i)
    {
        nodeLengthCoordinate[i] = nodeLengthCoordinate[i - 1] + edgeLengths[i - 1];
    }

    // Approximate number of nodes in the refined sections.
    const UInt numNodesRefinedPart = static_cast<UInt>(std::ceil((nodeLengthCoordinate[endIndex] - nodeLengthCoordinate[startIndex]) / refinementDistance)) + endIndex - startIndex;
    UInt numNodesNotRefinedPart = startIndex - outerStart + outerEnd - endIndex;
    // Approximate the number of nodes in the refined polygon.
    UInt totalNumNodes = numNodesRefinedPart + numNodesNotRefinedPart;

    std::vector<Point> refinedPolygon;

    refinedPolygon.reserve(totalNumNodes);

    // Add nodes before the section to be refined
    for (size_t i = outerStart; i <= startIndex; ++i)
    {
        refinedPolygon.emplace_back(m_nodes[i]);
    }

    // Refine each line segment.
    for (size_t i = startIndex; i < endIndex; ++i)
    {
        // Line segment starting point.
        Point p = m_nodes[i];
        const double segmentLength = ComputeDistance(m_nodes[i], m_nodes[i + 1], m_projection);
        // Refined segment step size.
        const Point delta = (m_nodes[i + 1] - m_nodes[i]) * refinementDistance / segmentLength;
        double lengthAlongInterval = refinementDistance;

        // Exit when the lengthAlongInterval is greater or equal than segmentLength
        // To prevent very small refined segment lengths, also exit when lengthAlongInterval is a small fraction less (defined by refinementTolerance)
        while (lengthAlongInterval < segmentLength && !IsEqual(lengthAlongInterval, segmentLength, constants::geometric::refinementTolerance))
        {
            p += delta;
            lengthAlongInterval += refinementDistance;
            refinedPolygon.emplace_back(p);
        }

        // Add last node to the refined polygon point sequence.
        refinedPolygon.emplace_back(m_nodes[i + 1]);
    }

    // Add nodes after the section to be refined
    for (size_t i = endIndex + 1; i <= outerEnd; ++i)
    {
        refinedPolygon.emplace_back(m_nodes[i]);
    }

    return refinedPolygon;
}

Polygons Polygons::OffsetCopy(double distance, bool innerAndOuter) const
{
    auto sizenewPolygon = GetNumNodes();
    if (innerAndOuter)
    {
        sizenewPolygon += GetNumNodes() + 1;
    }

    std::vector<Point> normalVectors(sizenewPolygon);
    double dxNormalPreviousEdge = 0.0;
    double dyNormalPreviousEdge = 0.0;
    double dxNormal = 0.0;
    double dyNormal = 0.0;
    for (UInt n = 0; n < GetNumNodes(); n++)
    {
        if (n < GetNumNodes() - 1)
        {
            const auto dx = GetDx(m_nodes[n], m_nodes[n + 1], m_projection);
            const auto dy = GetDy(m_nodes[n], m_nodes[n + 1], m_projection);
            const auto nodeDistance = std::sqrt(dx * dx + dy * dy);
            dxNormal = -dy / nodeDistance;
            dyNormal = dx / nodeDistance;
        }
        else
        {
            dxNormal = dxNormalPreviousEdge;
            dyNormal = dyNormalPreviousEdge;
        }

        if (n == 0)
        {
            dxNormalPreviousEdge = dxNormal;
            dyNormalPreviousEdge = dyNormal;
        }

        const double factor = 1.0 / (1.0 + dxNormalPreviousEdge * dxNormal + dyNormalPreviousEdge * dyNormal);
        normalVectors[n].x = factor * (dxNormalPreviousEdge + dxNormal);
        normalVectors[n].y = factor * (dyNormalPreviousEdge + dyNormal);

        dxNormalPreviousEdge = dxNormal;
        dyNormalPreviousEdge = dyNormal;
    }

    // negative sign introduced because normal vector pointing inward
    distance = -distance;
    if (m_projection == Projection::spherical)
    {
        distance = distance / (constants::geometric::earth_radius * constants::conversion::degToRad);
    }

    std::vector<Point> newPolygonPoints(sizenewPolygon, {constants::missing::doubleValue, constants::missing::doubleValue});
    for (UInt i = 0; i < GetNumNodes(); ++i)
    {
        auto dx = normalVectors[i].x * distance;
        const auto dy = normalVectors[i].y * distance;
        if (m_projection == Projection::spherical)
        {
            dx = dx / std::cos((m_nodes[i].y + 0.5 * dy) * constants::conversion::degToRad);
        }
        newPolygonPoints[i].x = m_nodes[i].x + dx;
        newPolygonPoints[i].y = m_nodes[i].y + dy;

        if (innerAndOuter)
        {
            newPolygonPoints[i + GetNumNodes() + 1].x = m_nodes[i].x - dx;
            newPolygonPoints[i + GetNumNodes() + 1].y = m_nodes[i].y - dy;
        }
    }

    // set the new polygon
    Polygons newPolygon{newPolygonPoints, m_projection};
    return newPolygon;
}

void Polygons::SnapToLandBoundary(const LandBoundary& landBoundary, UInt startIndex, UInt endIndex)
{
    if (m_nodes.empty())
    {
        throw ConstraintError(VariadicErrorMessage("No nodes in polygon."));
    }

    if (startIndex == 0 && endIndex == 0)
    {
        endIndex = static_cast<UInt>(m_nodes.size()) - 1;
    }

    if (startIndex >= endIndex)
    {
        throw ConstraintError(VariadicErrorMessage("The start index is greater than the end index: {} >= {}.", startIndex, endIndex));
    }

    for (size_t i = startIndex; i <= endIndex; ++i)
    {
        if (m_nodes[i].IsValid())
        {
            m_nodes[i] = landBoundary.FindNearestPoint(m_nodes[i], m_projection);
        }
    }
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

    return m_polygonGroups[polygonIndex].Contains(point);

    // const auto& [outerStart, outerEnd] = m_outer_polygons_indices[polygonIndex];
    // const auto inPolygon = IsPointInPolygonNodes(point, m_nodes, m_projection, Point(), outerStart, outerEnd);
    // return inPolygon;
}

meshkernel::UInt Polygons::GetNumPolygons() const
{
    return static_cast<UInt>(m_outer_polygons_indices.size());
}

std::tuple<bool, meshkernel::UInt> Polygons::IsPointInPolygons(Point point) const
{
    // empty polygon means everything is included
    if (m_outer_polygons_indices.empty())
    {
        return {true, constants::missing::uintValue};
    }

    for (size_t i = 0; i < m_polygonGroups.size(); ++i)
    {
        const PolygonalEnclosure& group = m_polygonGroups[i];

        if (group.ContainsRegion(point) == 1)
        {
            return {true, i};
        }
        else if (group.ContainsRegion(point) == 2)
        {
            // Can we just break here? It feels a bit like a goto
            return {false, constants::missing::uintValue};
        }
    }

    return {false, constants::missing::uintValue};

    bool inPolygon = false;
    for (UInt polygonIndex = 0; polygonIndex < GetNumPolygons(); ++polygonIndex)
    {
        const auto& [polygonStartIndex, polygonEndIndex] = m_outer_polygons_indices[polygonIndex];

        // Calculate the bounding box
        double XMin = std::numeric_limits<double>::max();
        double XMax = std::numeric_limits<double>::lowest();
        double YMin = std::numeric_limits<double>::max();
        double YMax = std::numeric_limits<double>::lowest();

        for (auto n = polygonStartIndex; n <= polygonEndIndex; n++)
        {
            XMin = std::min(XMin, m_nodes[n].x);
            XMax = std::max(XMax, m_nodes[n].x);
            YMin = std::min(YMin, m_nodes[n].y);
            YMax = std::max(YMax, m_nodes[n].y);
        }

        if (point.x >= XMin && point.x <= XMax && (point.y >= YMin && point.y <= YMax))
        {
            inPolygon = IsPointInPolygonNodes(point, m_nodes, m_projection, Point(), polygonStartIndex, polygonEndIndex);
        }

        if (inPolygon)
        {
            for (const auto& [startInner, endInner] : m_inner_polygons_indices.at(polygonIndex))
            {
                if (IsPointInPolygonNodes(point, m_nodes, m_projection, Point(), startInner, endInner))
                {
                    return {false, constants::missing::uintValue};
                }
            }

            return {true, polygonIndex};
        }
    }

    return {false, constants::missing::uintValue};
}

std::vector<bool> Polygons::PointsInPolygons(const std::vector<Point>& points) const
{
    std::vector<bool> result(points.size(), false);
    for (UInt i = 0; i < points.size(); ++i)
    {
        const auto [isInPolygon, polygonIndex] = IsPointInPolygons(points[i]);
        result[i] = isInPolygon;
    }
    return result;
}

bool Polygons::IsEmpty() const
{
    return m_outer_polygons_indices.empty();
}

double Polygons::PerimeterClosedPolygon(const std::vector<Point>& polygonNodes) const
{

    // TODO can this be removed, it is private and is not called in this file.

    if (polygonNodes.front() != polygonNodes.back())
    {
        throw std::invalid_argument("Polygons::PerimeterClosedPolygon: The first and last point of the polygon is not the same.");
    }

    const auto edgeLengths = PolygonEdgeLengths(polygonNodes);
    return std::accumulate(edgeLengths.begin(), edgeLengths.end(), 0.0);
}

std::vector<double> Polygons::PolygonEdgeLengths(const std::vector<Point>& polygonNodes) const
{
    std::vector<double> edgeLengths;
    edgeLengths.reserve(polygonNodes.size());

    for (UInt p = 0; p < polygonNodes.size(); ++p)
    {
        const auto firstNode = p;
        auto secondNode = p + 1;
        if (secondNode == polygonNodes.size())
        {
            secondNode = 0;
        }
        edgeLengths.emplace_back(ComputeDistance(polygonNodes[firstNode], polygonNodes[secondNode], m_projection));
    }
    return edgeLengths;
}

double Polygons::MaximumEdgeLength(const std::vector<Point>& polygonNodes) const
{

    // TODO can this be removed, it is private and is not called in this file.

    if (polygonNodes.front() != polygonNodes.back())
    {
        throw std::invalid_argument("Polygons::MaximumEdgeLength: The first and last point of the polygon is not the same.");
    }

    auto maximumEdgeLength = std::numeric_limits<double>::lowest();
    for (UInt p = 0; p < polygonNodes.size() - 1; ++p)
    {
        double edgeLength = ComputeDistance(m_nodes[p], m_nodes[p + 1], m_projection);
        maximumEdgeLength = std::max(maximumEdgeLength, edgeLength);
    }
    return maximumEdgeLength;
}

meshkernel::BoundingBox Polygons::GetBoundingBox(UInt polygonIndex) const
{
    std::vector<Point> points;
    auto const& [startPolygonIndex, endPolygonIndex] = OuterIndices(polygonIndex);
    for (auto i = startPolygonIndex; i <= endPolygonIndex; ++i)
    {
        auto const& polygonNode = m_nodes[i];
        if (polygonNode.IsValid())
        {
            points.emplace_back(polygonNode);
        }
    }
    return BoundingBox(points);
}
