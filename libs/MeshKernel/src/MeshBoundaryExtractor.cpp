//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2026.
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

#include "MeshKernel/MeshBoundaryExtractor.hpp"

#include <cmath>
#include <numbers>
#include <ranges>
#include <set>

#include "MeshKernel/Operations.hpp"

std::vector<meshkernel::Point> meshkernel::MeshBoundaryExtractor::ExtractConcatenated(const Mesh2D& mesh, BoundarySelection boundaryType)
{

    // The use of std::tie instead of a structured binding is due to limitations in the macos compiler
    // The compiler error "error: capturing a structured binding is not yet supported in OpenMP"
    //
    // auto [boundarySequences, isExterior] = Extract(mesh);
    // Replace the 3 lines below with the line above
    std::vector<std::vector<Point>> boundarySequences;
    std::vector<bool> isExterior;
    std::tie(boundarySequences, isExterior) = Extract(mesh);

    auto isBoundarySelection = [boundaryType, &isExterior](size_t idx)
    {
        using enum BoundarySelection;

        return (boundaryType == All) ||
               (boundaryType == ExteriorOnly && isExterior[idx]) ||
               (boundaryType == InteriorOnly && !isExterior[idx]);
    };

    return ConcatenatePointVectors(boundarySequences, isBoundarySelection);
}

std::tuple<std::vector<std::vector<meshkernel::Point>>, std::vector<bool>> meshkernel::MeshBoundaryExtractor::Extract(const Mesh2D& mesh)
{
    const Polygon emptyPolygon;
    return Extract(mesh, emptyPolygon);
}

std::tuple<std::vector<std::vector<meshkernel::Point>>, std::vector<bool>> meshkernel::MeshBoundaryExtractor::Extract(const Mesh2D& mesh, const Polygon& polygon)
{

    std::vector<std::vector<Point>> allBoundaryPolygons;
    std::vector<std::vector<UInt>> allTouchedFaces;

    const std::vector<Point>& meshNodes(mesh.Nodes());

    FindBoundaryPolygons(meshNodes, mesh.Edges(), mesh.m_edgesFaces, allBoundaryPolygons, allTouchedFaces);
    return SeparateAndDetermineExternality(mesh, polygon, allBoundaryPolygons, allTouchedFaces);
}

bool meshkernel::MeshBoundaryExtractor::IsMultiPolygon(std::span<const Point> boundary)
{

    if (boundary.size() < 4)
    {
        return false;
    }

    std::set<Point> uniquePoints;

    // Since polygons are always closed (first = last) then start from 1 after the first
    for (size_t i = 1; i < boundary.size(); ++i)
    {
        if (!uniquePoints.insert(boundary[i]).second)
        {
            return true;
        }
    }

    return false;
}

std::tuple<std::vector<std::vector<meshkernel::Point>>, std::vector<meshkernel::UInt>> meshkernel::MeshBoundaryExtractor::SplitMultiplePolygons(std::span<const Point> boundaryPoints,
                                                                                                                                                std::span<const UInt> elementIds)
{
    std::vector<std::vector<Point>> completedPolygons;
    std::vector<UInt> firstElementIds;

    std::vector<std::pair<Point, int>> pointFaceStack;
    // Maps a point to its current index in the pointFaceStack
    std::map<Point, size_t> activePoints;

    for (size_t i = 0; i < boundaryPoints.size(); ++i)
    {
        const Point& currentPoint = boundaryPoints[i];

        if (!currentPoint.IsValid())
        {
            continue;
        }

        // Check if this node closes a loop with a previously visited point
        if (auto it = activePoints.find(currentPoint); it != activePoints.end())
        {
            size_t loopStartIndex = it->second;

            std::vector<Point> subPolygon;
            subPolygon.reserve((pointFaceStack.size() - loopStartIndex) + 1);

            int firstEdgeId = pointFaceStack[loopStartIndex].second;

            for (size_t j = loopStartIndex; j < pointFaceStack.size(); ++j)
            {
                subPolygon.push_back(pointFaceStack[j].first);
                activePoints.erase(pointFaceStack[j].first);
            }

            // close the polygon
            if (!subPolygon.empty())
            {
                subPolygon.push_back(subPolygon.front());
            }

            completedPolygons.push_back(subPolygon);
            firstElementIds.push_back(firstEdgeId);

            pointFaceStack.resize(loopStartIndex);
        }

        int currentEdge = (i < elementIds.size()) ? elementIds[i] : -1;

        activePoints[currentPoint] = pointFaceStack.size();
        pointFaceStack.emplace_back(currentPoint, currentEdge);
    }

    return {completedPolygons, firstElementIds};
}

void meshkernel::MeshBoundaryExtractor::ClipToConstrainingPolygon(const Polygon& polygon, std::vector<Point>& nodes)
{
    if (polygon.IsEmpty() || nodes.size() <= MinimumNumberOfPoints)
    {
        return;
    }

    const size_t uniquePointCount = nodes.size() - 1;
    std::vector<Point> clippedPolygon;
    std::vector<bool> nodeIsContained(uniquePointCount);
    UInt numberOfNodesContained = 0;

    for (size_t i = 0; i < uniquePointCount; ++i)
    {
        nodeIsContained[i] = polygon.Contains(nodes[i]);
        numberOfNodesContained += nodeIsContained[i] ? 1 : 0;
    }

    if (numberOfNodesContained < MinimumNumberOfPoints)
    {
        // Cannot make a reasonable polygon with only 2 nodes (and the closing node htat will be added later, making 3 nodes)
        nodes.clear();
        return;
    }

    clippedPolygon.reserve(nodes.size());

    // At this point the boundary polygon is closed, so we have to ignore the last point on the sequence.
    for (size_t i = 0; i < uniquePointCount; ++i)
    {
        size_t nextNodeId = (i + 1) % uniquePointCount;
        size_t previousNodeId = (i + uniquePointCount - 1) % uniquePointCount;

        if (nodeIsContained[i] || nodeIsContained[nextNodeId] || nodeIsContained[previousNodeId])
        {
            clippedPolygon.push_back(nodes[i]);
        }
    }

    // Now re-close the polygon
    clippedPolygon.push_back(clippedPolygon[0]);

    nodes = std::move(clippedPolygon);
}

void meshkernel::MeshBoundaryExtractor::Append(const Polygon& polygon,
                                               const Point& centre,
                                               const Projection projection,
                                               std::vector<Point>& boundaryPolygon,
                                               std::vector<bool>& isExterior,
                                               std::vector<std::vector<Point>>& separatedBoundaryPolygons)
{

    if (auto [area, centreOfMass, orientation] = Polygon::FaceAreaAndCenterOfMass(boundaryPolygon, projection); orientation == TraversalDirection::Clockwise)
    {
        std::ranges::reverse(boundaryPolygon);
    }

    const bool isExteriorBoundary = IsPointInPolygonNodes(centre, boundaryPolygon, projection);

    // Clipping to boundary polygon must be done after determining if the boundary-points form an exterior or interior boundary.
    // Because the element centre (taken from the first element the boundary-polygon touches) may no longer be inside the clipped
    // boundary-polygon
    ClipToConstrainingPolygon(polygon, boundaryPolygon);

    if (boundaryPolygon.size() > MinimumNumberOfPoints)
    {
        // Only add the polygon and the exterior/interior indicator if the polygon has a sufficient number of points
        isExterior.push_back(isExteriorBoundary);
        separatedBoundaryPolygons.push_back(std::move(boundaryPolygon));
    }
}

std::tuple<std::vector<std::vector<meshkernel::Point>>, std::vector<bool>>
meshkernel::MeshBoundaryExtractor::SeparateAndDetermineExternality(const Mesh2D& mesh,
                                                                   const Polygon& polygon,
                                                                   std::vector<std::vector<Point>>& allBoundaryPolygons,
                                                                   const std::vector<std::vector<UInt>>& allTouchedFaces)
{
    std::vector<std::vector<Point>> separatedBoundaryPolygons;
    std::vector<bool> isExterior;

    for (size_t i = 0; i < allBoundaryPolygons.size(); ++i)
    {

        if (IsMultiPolygon(allBoundaryPolygons[i]))
        {
            // It can be that some of the boundaries that are found are composed of multiple sub-boundaries.
            // Some of the sub-boundaries may be combined in a non conforming way.
            // In either case, the boundaries are separated into distinct boundary polygons.
            auto [individualBoundaryPolygons, firstElement] = SplitMultiplePolygons(allBoundaryPolygons[i], allTouchedFaces[i]);

            for (size_t j = 0; j < individualBoundaryPolygons.size(); ++j)
            {
                const Point centre = mesh.m_facesMassCenters[firstElement[j]];

                Append(polygon, centre, mesh.m_projection, individualBoundaryPolygons[j], isExterior, separatedBoundaryPolygons);
            }
        }
        else
        {
            const Point centre = mesh.m_facesMassCenters[allTouchedFaces[i][0]];

            Append(polygon, centre, mesh.m_projection, allBoundaryPolygons[i], isExterior, separatedBoundaryPolygons);
        }
    }

    return {separatedBoundaryPolygons, isExterior};
}

double meshkernel::MeshBoundaryExtractor::NormalizeAngle(double angle)
{
    while (angle < 0)
    {
        angle += 2.0 * std::numbers::pi;
    }

    while (angle >= 2.0 * std::numbers::pi)
    {
        angle -= 2.0 * std::numbers::pi;
    }

    return angle;
}

void meshkernel::MeshBoundaryExtractor::FindAllBoundarEdges(const std::vector<Point>& nodes,
                                                            const std::vector<Edge>& edges,
                                                            const std::vector<std::array<UInt, 2>>& edgesFaces,
                                                            std::unordered_map<UInt, std::vector<BoundaryEdge>>& boundaryAdjacency)
{

    // Collect all boundary edges and compute the angle
    for (UInt count = 0; count < edges.size(); ++count)
    {

        if (edgesFaces[count][1] != constants::missing::uintValue || !IsValidEdge(edges[count]))
        {
            // Edge is either invalid or not on boundary
            continue;
        }

        double dx = nodes[edges[count].second].x - nodes[edges[count].first].x;
        double dy = nodes[edges[count].second].y - nodes[edges[count].first].y;

        // Compute angle for "crossong/pinch point" sorting
        // Pinch points are when the boundary polygon intersects itself, the points will have identical values.
        double angleStartToEnd = NormalizeAngle(std::atan2(dy, dx));
        double angleEndToStart = NormalizeAngle(std::atan2(-dy, -dx));

        // Note: edgesFaces[count][0] is the internal mesh face valid for both directions of boundary traversal
        boundaryAdjacency[edges[count].first].push_back({count, edges[count].second, edgesFaces[count][0], angleStartToEnd});
        boundaryAdjacency[edges[count].second].push_back({count, edges[count].first, edgesFaces[count][0], angleEndToStart});
    }

    auto boundaryAngleLessThan = [](const BoundaryEdge& a, const BoundaryEdge& b)
    { return a.angle < b.angle; };

    // Sort outgoing edges anti-clockwise
    for (auto& [node_id, connected_edges] : boundaryAdjacency)
    {
        std::sort(connected_edges.begin(), connected_edges.end(), boundaryAngleLessThan);
    }
}

meshkernel::UInt meshkernel::MeshBoundaryExtractor::FindEdgeWithMinumumAngle(const std::vector<BoundaryEdge>& boundaryEdges,
                                                                             const std::vector<bool>& edgeVisited,
                                                                             const double incomingAngle)
{

    // Find the smallest visited angle.
    // All edge angles must be in interval [0, 2pi], initialise with number greater than 2pi
    double deltaAngle = 3.0 * std::numbers::pi;
    UInt edgeIndex = constants::missing::uintValue;

    // Find the boundary edge that tracks closest clockwise to keep empty space on the right
    for (UInt e = 0; e < boundaryEdges.size(); ++e)
    {
        if (edgeVisited[boundaryEdges[e].edgeId])
        {
            continue;
        }

        if (double delta = NormalizeAngle(incomingAngle - boundaryEdges[e].angle); delta < deltaAngle)
        {
            deltaAngle = delta;
            edgeIndex = e;
        }
    }

    return edgeIndex;
}

void meshkernel::MeshBoundaryExtractor::FindBoundaryPolygons(const std::vector<Point>& nodes,
                                                             const std::vector<Edge>& edges,
                                                             const std::vector<std::array<UInt, 2>>& edgesFaces,
                                                             std::vector<std::vector<Point>>& allPolygons,
                                                             std::vector<std::vector<UInt>>& allTouchedFaces)
{
    allPolygons.clear();
    allTouchedFaces.clear();

    // Mapping from mesh node-id to a sequence of boundary edges, the boundary edges will be sorted by angle
    std::unordered_map<UInt, std::vector<BoundaryEdge>> boundaryAdjacency;

    FindAllBoundarEdges(nodes, edges, edgesFaces, boundaryAdjacency);

    std::vector<bool> edgeVisited(edges.size(), false);

    // Trace boundary polygons
    for (UInt count = 0; count < edges.size(); ++count)
    {

        if (!IsValidEdge(edges[count]))
        {
            continue;
        }

        if (edgesFaces[count][1] != constants::missing::uintValue || edgeVisited[count])
        {
            continue;
        }

        std::vector<Point> currentNodes;
        std::vector<UInt> currentFaces;

        UInt prevNodeIndex = edges[count].first;
        UInt currentNodeIndex = edges[count].second;

        currentNodes.push_back(nodes[prevNodeIndex]);
        currentFaces.push_back(edgesFaces[count][0]); // First face touched by the initial edge
        edgeVisited[count] = true;                    // Mark edge as having been visited

        double dx = nodes[currentNodeIndex].x - nodes[prevNodeIndex].x;
        double dy = nodes[currentNodeIndex].y - nodes[prevNodeIndex].y;
        double incomingAngle = NormalizeAngle(std::atan2(dy, dx));

        // Loop until we find the start node, making a closed boundary polygon
        while (currentNodeIndex != edges[count].first)
        {
            currentNodes.push_back(nodes[currentNodeIndex]);

            const std::vector<BoundaryEdge>& boundaryEdges = boundaryAdjacency[currentNodeIndex];

            UInt edgeIndex = FindEdgeWithMinumumAngle(boundaryEdges, edgeVisited, incomingAngle);

            // No unvisited edges were found
            // So eigher the boundary polygon was completed or a dead-end reached.
            if (edgeIndex == constants::missing::uintValue)
            {
                break;
            }

            const BoundaryEdge& chosenEdge = boundaryEdges[edgeIndex];
            edgeVisited[chosenEdge.edgeId] = true;

            currentFaces.push_back(chosenEdge.leftFace);

            prevNodeIndex = currentNodeIndex;
            currentNodeIndex = chosenEdge.neighbourNode;
            incomingAngle = chosenEdge.angle;
        }

        if (currentNodes.size() >= MinimumNumberOfPoints)
        {
            // Close the polygon
            // The polygon needs to be closed here. If it is composed of multiple sub-polygons then closing here
            // makes separating them later a much easier task
            currentNodes.push_back(currentNodes.front());
            allPolygons.emplace_back(std::move(currentNodes));
            allTouchedFaces.push_back(std::move(currentFaces));
        }
    }
}
