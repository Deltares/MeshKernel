#include "MeshKernel/MeshBoundaryExtractor.hpp"

#include <cmath>
#include <numbers>
#include <ranges>

#include "MeshKernel/Operations.hpp"

std::vector<meshkernel::Point> meshkernel::MeshBoundaryExtractor::ExtractConcatenated(const Mesh2D& mesh, BoundarySelection boundaryType)
{

    // The use of std::tie instead of a structured binding is due to limitations in the macos compiler
    //
    // auto [boundarySequences, isExterior] = Extract(mesh);
    // Replace the 3 lines below with the line above
    std::vector<std::vector<Point>> boundarySequences;
    std::vector<bool> isExterior;
    std::tie(boundarySequences, isExterior) = Extract(mesh);

    std::vector<Point> allPoints;

    auto isBoundarySelection = [boundaryType, &isExterior](size_t idx)
    {
        using enum BoundarySelection;

        return (boundaryType == All) ||
               (boundaryType == ExteriorOnly && isExterior[idx]) ||
               (boundaryType == InteriorOnly && !isExterior[idx]);
    };

    return ConcatenatePointVectors(boundarySequences, isBoundarySelection);
}

void meshkernel::MeshBoundaryExtractor::Append(const Point& centre,
                                               const Projection projection,
                                               std::vector<Point>& boundaryPolygon,
                                               std::vector<bool>& isExterior,
                                               std::vector<std::vector<Point>>& separatedBoundaryPolygons)
{

    if (auto [area, centreOfMass] = ComputePolygonAreaAndCentre(boundaryPolygon, projection); area < 0.0)
    {
        std::ranges::reverse(boundaryPolygon);
    }

    isExterior.push_back(IsPointInPolygonNodes(centre, boundaryPolygon, projection));
    separatedBoundaryPolygons.push_back(std::move(boundaryPolygon));
}

std::tuple<std::vector<std::vector<meshkernel::Point>>, std::vector<bool>>
meshkernel::MeshBoundaryExtractor::SeparateAndDetermineExternality(const Mesh2D& mesh,
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

            for (size_t i = 0; i < individualBoundaryPolygons.size(); ++i)
            {
                Point centre = mesh.m_facesMassCenters[firstElement[i]];
                Append(centre, mesh.m_projection, individualBoundaryPolygons[i], isExterior, separatedBoundaryPolygons);
            }
        }
        else
        {
            Point centre = mesh.m_facesMassCenters[allTouchedFaces[i][0]];
            Append(centre, mesh.m_projection, allBoundaryPolygons[i], isExterior, separatedBoundaryPolygons);
        }
    }

    return {separatedBoundaryPolygons, isExterior};
}

std::tuple<std::vector<std::vector<meshkernel::Point>>, std::vector<bool>> meshkernel::MeshBoundaryExtractor::Extract(const Mesh2D& mesh)
{

    std::vector<std::vector<Point>> allBoundaryPolygons;
    std::vector<std::vector<UInt>> allTouchedFaces;

    const std::vector<Point>& meshNodes(mesh.Nodes());

    FindBoundaryPolygons(meshNodes, mesh.Edges(), mesh.m_edgesFaces, allBoundaryPolygons, allTouchedFaces);
    return SeparateAndDetermineExternality(mesh, allBoundaryPolygons, allTouchedFaces);
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

        if (!IsValidEdge(edges[count]) || edgesFaces[count][1] != constants::missing::uintValue)
        {
            // Edge is either invalid or not on boundaryy
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
    for (size_t e = 0; e < boundaryEdges.size(); ++e)
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

        // Polygon until we find the start node, making a closed boundary polygon
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

        if (currentNodes.size() >= 3)
        {
            // Close the polygon
            currentNodes.push_back(currentNodes.front());
            allPolygons.emplace_back(std::move(currentNodes));
            allTouchedFaces.push_back(std::move(currentFaces));
        }
    }
}
