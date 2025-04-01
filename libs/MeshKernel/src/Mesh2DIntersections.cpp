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

#include <ranges>

#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Mesh2DIntersections.hpp"
#include "MeshKernel/Operations.hpp"

using namespace meshkernel;

Mesh2DIntersections::Mesh2DIntersections(Mesh2D& mesh) : m_mesh(mesh)
{
    if (mesh.GetNumNodes() <= 0)
    {
        throw AlgorithmError("Mesh with no nodes");
    }

    mesh.Administrate();

    if (mesh.GetNumFaces() <= 0)
    {
        throw AlgorithmError("Mesh with no faces");
    }

    // Intersection results
    m_edgesIntersections.resize(mesh.GetNumEdges());
    m_faceIntersections.resize(mesh.GetNumFaces());

    // Declare local caches
    m_edgesIntersectionsCache.resize(mesh.GetNumEdges());
    m_facesIntersectionsCache.resize(mesh.GetNumFaces());

    m_meshBoundingBox = m_mesh.GetBoundingBox();
    m_meshEdgesBoundingBoxes = m_mesh.GetEdgesBoundingBoxes();
}

std::tuple<UInt, UInt> Mesh2DIntersections::GetIntersectionSeed(const Mesh2D& mesh,
                                                                const std::vector<Point>& polyLine,
                                                                const UInt polygonIndexStart,
                                                                const bool checkOnlyBoundarySegments,
                                                                const std::vector<BoundingBox>& polyLineBoundingBoxes,
                                                                const std::vector<bool>& vistedEdges) const
{
    UInt crossedSegmentIndex = constants::missing::uintValue;
    UInt crossedEdgeIndex = constants::missing::uintValue;
    bool isSeedFound = false;

    // Find starting edge and segment
    for (UInt segmentIndex = polygonIndexStart; segmentIndex < polyLine.size() - 1; ++segmentIndex)
    {
        for (UInt edgeIndex = 0; edgeIndex < mesh.GetNumEdges(); ++edgeIndex)
        {
            // edge already crossed, nothing to do
            if (vistedEdges[edgeIndex])
            {
                continue;
            }

            if (checkOnlyBoundarySegments && !mesh.IsEdgeOnBoundary(edgeIndex))
            {
                continue;
            }

            if (!m_meshBoundingBox.Overlaps(polyLineBoundingBoxes[segmentIndex]))
            {
                continue;
            }

            if (!m_meshEdgesBoundingBoxes[edgeIndex].Overlaps(polyLineBoundingBoxes[segmentIndex]))
            {
                continue;
            }

            const auto [isEdgeCrossed,
                        intersectionPoint,
                        crossProductValue,
                        adimensionalPolylineSegmentDistance,
                        adimensionalEdgeDistance] = AreSegmentsCrossing(polyLine[segmentIndex],
                                                                        polyLine[segmentIndex + 1],
                                                                        mesh.Node(mesh.GetEdge(edgeIndex).first),
                                                                        mesh.Node(mesh.GetEdge(edgeIndex).second),
                                                                        false,
                                                                        mesh.m_projection);
            if (!isEdgeCrossed)
            {
                continue;
            }

            crossedSegmentIndex = segmentIndex;
            crossedEdgeIndex = edgeIndex;
            isSeedFound = true;
            break;
        }

        if (isSeedFound)
        {
            break;
        }
    }

    return {crossedEdgeIndex, crossedSegmentIndex};
}

std::tuple<bool, UInt, UInt, double, double, double> Mesh2DIntersections::GetNextEdgeIntersection(const std::vector<Point>& polyLine,
                                                                                                  const std::vector<BoundingBox>& polyLineBoundingBoxes,
                                                                                                  UInt edgeIndex,
                                                                                                  UInt firstIndex,
                                                                                                  UInt secondIndex,
                                                                                                  Direction direction) const
{
    UInt numSteps = 0;
    bool intersectionFound = false;
    double crossProductValue = constants::missing::doubleValue;
    double adimensionalPolylineSegmentDistance = constants::missing::doubleValue;
    double adimensionalEdgeDistance = constants::missing::doubleValue;
    const auto polyLineSize = static_cast<UInt>(polyLine.size());

    const auto checkPolyLineIndex = [&]
    {
        // forward
        if (direction == Direction::Forward)
        {
            return firstIndex < polyLineSize - 2;
        }
        // backward
        if (direction == Direction::Backward)
        {
            return firstIndex >= 1;
        }
        return false;
    };

    while (!intersectionFound && checkPolyLineIndex() && numSteps < maxSearchSegments)
    {
        // forward
        if (direction == Direction::Forward)
        {
            firstIndex = secondIndex;
            secondIndex = firstIndex + 1;
        }

        // backward
        if (direction == Direction::Backward)
        {
            secondIndex = firstIndex;
            firstIndex = firstIndex - 1;
        }

        if (!m_meshBoundingBox.Overlaps(polyLineBoundingBoxes[firstIndex]))
        {
            numSteps++;
            continue;
        }

        if (!m_meshEdgesBoundingBoxes[edgeIndex].Overlaps(polyLineBoundingBoxes[firstIndex]))
        {
            numSteps++;
            continue;
        }

        auto intersection = AreSegmentsCrossing(polyLine[firstIndex],
                                                polyLine[secondIndex],
                                                m_mesh.Node(m_mesh.GetEdge(edgeIndex).first),
                                                m_mesh.Node(m_mesh.GetEdge(edgeIndex).second),
                                                false,
                                                m_mesh.m_projection);

        intersectionFound = std::get<0>(intersection);
        crossProductValue = std::get<2>(intersection);
        adimensionalPolylineSegmentDistance = std::get<3>(intersection);
        adimensionalEdgeDistance = std::get<4>(intersection);

        numSteps++;
    }

    return {intersectionFound, firstIndex, secondIndex, crossProductValue, adimensionalPolylineSegmentDistance, adimensionalEdgeDistance};
}

void Mesh2DIntersections::IntersectFaceEdges(const std::vector<Point>& polyLine,
                                             const std::vector<BoundingBox>& polyLineBoundingBoxes,
                                             const std::vector<double>& cumulativeLength,
                                             UInt currentCrossingEdge,
                                             UInt currentFaceIndex,
                                             UInt segmentIndex,
                                             std::vector<bool>& vistedEdges,
                                             std::vector<bool>& vistedFace,
                                             std::queue<std::array<UInt, 2>>& crossingEdges)
{
    for (UInt e = 0; e < m_mesh.m_facesEdges[currentFaceIndex].size(); ++e)
    {
        const auto edgeIndex = m_mesh.m_facesEdges[currentFaceIndex][e];
        if (vistedEdges[edgeIndex] && vistedFace[currentFaceIndex])
        {
            continue;
        }

        UInt firstIndex = segmentIndex;
        UInt secondIndex = segmentIndex + 1;

        auto intersection = AreSegmentsCrossing(polyLine[firstIndex],
                                                polyLine[secondIndex],
                                                m_mesh.Node(m_mesh.GetEdge(edgeIndex).first),
                                                m_mesh.Node(m_mesh.GetEdge(edgeIndex).second),
                                                false,
                                                m_mesh.m_projection);

        bool intersectionFound = std::get<0>(intersection);
        double crossProductValue = std::get<2>(intersection);
        double adimensionalPolylineSegmentDistance = std::get<3>(intersection);
        double adimensionalEdgeDistance = std::get<4>(intersection);

        if (!intersectionFound)
        {
            std::tie(intersectionFound,
                     firstIndex,
                     secondIndex,
                     crossProductValue,
                     adimensionalPolylineSegmentDistance,
                     adimensionalEdgeDistance) = GetNextEdgeIntersection(polyLine, polyLineBoundingBoxes, edgeIndex, segmentIndex, segmentIndex + 1, Direction::Forward);
        }

        if (!intersectionFound)
        {
            std::tie(intersectionFound,
                     firstIndex,
                     secondIndex,
                     crossProductValue,
                     adimensionalPolylineSegmentDistance,
                     adimensionalEdgeDistance) = GetNextEdgeIntersection(polyLine, polyLineBoundingBoxes, edgeIndex, segmentIndex, segmentIndex + 1, Direction::Backward);
        }

        // none of the polyline intersect the current edge
        if (!intersectionFound)
        {
            continue;
        }

        updateEdgeIntersections(
            firstIndex,
            edgeIndex,
            m_mesh.GetEdge(edgeIndex),
            cumulativeLength,
            crossProductValue,
            adimensionalEdgeDistance,
            adimensionalPolylineSegmentDistance,
            m_edgesIntersectionsCache);

        updateFaceIntersections(currentFaceIndex, edgeIndex, m_facesIntersectionsCache);

        if (edgeIndex != currentCrossingEdge)
        {
            crossingEdges.push({edgeIndex, firstIndex});
        }
        vistedEdges[edgeIndex] = true;
    }
    vistedFace[currentFaceIndex] = true;
}

void Mesh2DIntersections::Compute(const std::vector<Point>& polyLine)
{
    // 1. Find the intersection of any segment of the polyline with the mesh return if nothing is found
    std::ranges::fill(m_edgesIntersectionsCache, EdgeMeshPolyLineIntersection());
    std::ranges::fill(m_facesIntersectionsCache, FaceMeshPolyLineIntersection());

    const auto polyLineSize = static_cast<UInt>(polyLine.size());

    std::vector<double> cumulativeLength(polyLine.size(), 0.0);
    std::vector<BoundingBox> polyLineBoundingBoxes(polyLine.size() - 1);
    for (UInt i = 1; i < polyLineSize; ++i)
    {
        cumulativeLength[i] = cumulativeLength[i - 1] + ComputeDistance(polyLine[i], polyLine[i - 1], m_mesh.m_projection);

        polyLineBoundingBoxes[i - 1] = BoundingBox::CreateBoundingBox(polyLine[i - 1], polyLine[i]);
    }

    std::queue<std::array<UInt, 2>> crossingEdges;
    std::vector<bool> vistedEdges(m_mesh.GetNumEdges(), false);
    std::vector<bool> vistedFaces(m_mesh.GetNumEdges(), false);

    // Keep a track of the last polygon segment that checked
    UInt polygonIndexStart = 0;
    bool checkOnlyBoundarySegments = false;

    // keep traversing the polyline as long crossed edges are found
    while (true)
    {
        // find the seed
        const auto [crossedEdgeIndex, crossedSegmentIndex] = GetIntersectionSeed(
            m_mesh,
            polyLine,
            polygonIndexStart,
            checkOnlyBoundarySegments,
            polyLineBoundingBoxes,
            vistedEdges);

        polygonIndexStart = crossedSegmentIndex;
        checkOnlyBoundarySegments = true;

        // no valid seed found in the entire mesh, we are done
        if (crossedEdgeIndex == constants::missing::uintValue)
        {
            break;
        }

        crossingEdges.push({crossedEdgeIndex, crossedSegmentIndex});

        // use breadth search along the current polyline
        while (!crossingEdges.empty())
        {
            auto [currentCrossingEdge, segmentIndex] = crossingEdges.front();
            crossingEdges.pop();

            for (const auto currentFaceIndex : m_mesh.m_edgesFaces[currentCrossingEdge])
            {
                if (currentFaceIndex == constants::missing::uintValue)
                {
                    continue;
                }
                IntersectFaceEdges(polyLine,
                                   polyLineBoundingBoxes,
                                   cumulativeLength,
                                   currentCrossingEdge,
                                   currentFaceIndex,
                                   segmentIndex,
                                   vistedEdges,
                                   vistedFaces,
                                   crossingEdges);
            }
        }
    }

    // Sort the edges for each face
    for (auto& facesIntersection : m_facesIntersectionsCache)
    {
        if (facesIntersection.edgeIndices.empty())
        {
            continue;
        }

        // sort edge indices based on polyline segment distance
        std::ranges::sort(facesIntersection.edgeIndices, [&](auto first, auto second)
                          { return m_edgesIntersectionsCache[first].adimensionalPolylineSegmentDistance <
                                   m_edgesIntersectionsCache[second].adimensionalPolylineSegmentDistance; });

        // compute the polylineDistance for the face
        double distanceSum = 0.0;
        for (const auto& edgeIndex : facesIntersection.edgeIndices)
        {
            distanceSum += m_edgesIntersectionsCache[edgeIndex].polylineDistance;
        }

        facesIntersection.polylineDistance = distanceSum / static_cast<double>(facesIntersection.edgeIndices.size());

        // push back the face intersection edge nodes
        for (UInt e = 0; e < facesIntersection.edgeIndices.size(); ++e)
        {
            const auto edgeIndex = facesIntersection.edgeIndices[e];
            facesIntersection.edgeNodes.emplace_back(m_edgesIntersectionsCache[edgeIndex].edgeFirstNode);
            facesIntersection.edgeNodes.emplace_back(m_edgesIntersectionsCache[edgeIndex].edgeSecondNode);
        }
    }

    // edge intersections are unique
    for (UInt e = 0; e < m_mesh.GetNumEdges(); ++e)
    {
        if (m_edgesIntersections[e].polylineDistance < 0)
        {
            m_edgesIntersections[e] = m_edgesIntersectionsCache[e];
        }
    }

    // face intersections are not unique and a face could have been intersected already
    for (UInt f = 0; f < m_mesh.GetNumFaces(); ++f)
    {
        if (!m_faceIntersections[f].edgeNodes.empty() &&
            !m_facesIntersectionsCache[f].edgeNodes.empty())
        {
            m_faceIntersections[f].edgeIndices.insert(m_faceIntersections[f].edgeIndices.end(), m_facesIntersectionsCache[f].edgeIndices.begin(), m_facesIntersectionsCache[f].edgeIndices.end());
            m_faceIntersections[f].edgeNodes.insert(m_faceIntersections[f].edgeNodes.end(), m_facesIntersectionsCache[f].edgeNodes.begin(), m_facesIntersectionsCache[f].edgeNodes.end());
            m_faceIntersections[f].polylineDistance = 0.5 * (m_faceIntersections[f].polylineDistance + m_facesIntersectionsCache[f].polylineDistance);
        }
        else if (!m_facesIntersectionsCache[f].edgeNodes.empty())
        {
            m_faceIntersections[f] = m_facesIntersectionsCache[f];
        }
    }
}

void Mesh2DIntersections::Compute(const Polygons& polygon)
{
    // No polygon, nothing is crossed
    if (polygon.IsEmpty())
    {
        return;
    }

    // Multiple polygons here
    for (auto outer = 0u; outer < polygon.GetNumPolygons(); ++outer)
    {
        std::vector<std::vector<Point>> polygonPolylines;

        polygonPolylines.emplace_back(polygon.Enclosure(outer).Outer().Nodes());

        for (auto inner = 0u; inner < polygon.Enclosure(outer).NumberOfInner(); ++inner)
        {
            polygonPolylines.emplace_back(polygon.Enclosure(outer).Inner(inner).Nodes());
        }

        for (const auto& polyLine : polygonPolylines)
        {
            Compute(polyLine);
        }
    }
}
