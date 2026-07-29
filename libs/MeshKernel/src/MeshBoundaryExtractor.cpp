#include "MeshKernel/MeshBoundaryExtractor.hpp"

#include <numbers>

#include "MeshKernel/Operations.hpp"

std::vector<meshkernel::Point> meshkernel::MeshBoundaryExtractor::ExtractAll(const Mesh2D& mesh)
{

    auto [boundarySequences, isExterior] = Extract(mesh);
    std::vector<Point> allPoints;

    bool isFirst = true;

    for (const std::vector<Point>& loop : boundarySequences)
    {
        if (!isFirst)
        {
            allPoints.push_back({constants::missing::doubleValue, constants::missing::doubleValue});
        }

        allPoints.insert(allPoints.end(), loop.begin(), loop.end());
        isFirst = false;
    }

    return allPoints;
}

std::tuple<std::vector<std::vector<meshkernel::Point>>, std::vector<bool>> meshkernel::MeshBoundaryExtractor::Extract(const Mesh2D& mesh)
{

    std::vector<std::vector<Point>> allBoundaryLoops;
    std::vector<std::vector<UInt>> allTouchedFaces;

    std::vector<std::vector<Point>> separatedBoundaryLoops;
    std::vector<bool> isExterior;

    const std::vector<Point>& meshNodes(mesh.Nodes());

    FindBoundaryLoops(meshNodes, mesh.Edges(), mesh.m_edgesFaces, allBoundaryLoops, allTouchedFaces);

    for (size_t i = 0; i < allBoundaryLoops.size(); ++i)
    {
        if (isMultiPolygon(allBoundaryLoops[i]))
        {
            auto [individualBoundaryPolygons, firstElement] = splitMultiplePolygons(allBoundaryLoops[i], allTouchedFaces[i]);

            for (size_t i = 0; i < individualBoundaryPolygons.size(); ++i)
            {
                Point centre = mesh.m_facesMassCenters[firstElement[i]];
                isExterior.push_back(IsPointInPolygonNodes(centre, individualBoundaryPolygons[i], mesh.m_projection));
                separatedBoundaryLoops.push_back(std::move(individualBoundaryPolygons[i]));
            }
        }
        else
        {
            separatedBoundaryLoops.push_back(std::move(allBoundaryLoops[i]));
        }
    }

    return {separatedBoundaryLoops, isExterior};
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

void meshkernel::MeshBoundaryExtractor::FindBoundaryLoops(const std::vector<Point>& nodes,
                                                          const std::vector<Edge>& edges,
                                                          const std::vector<std::array<UInt, 2>>& edgesFaces,
                                                          std::vector<std::vector<Point>>& allLoops,
                                                          std::vector<std::vector<UInt>>& allTouchedFaces)
{
    allLoops.clear();
    allTouchedFaces.clear();

    // Mapping from mesh node-id to a sequence of boundary edges, the boundary edges will be sorted by angle
    std::unordered_map<UInt, std::vector<BoundaryEdge>> boundaryAdjacency;

    // May be better to use meshkernel::Boolean
    std::vector<bool> edgeVisited(edges.size(), false);

    // Collect all boundary edges and compute the angle
    for (UInt count = 0; count < edges.size(); ++count)
    {

        if (edgesFaces[count][1] == constants::missing::uintValue)
        {
            double dx = nodes[edges[count].second].x - nodes[edges[count].first].x;
            double dy = nodes[edges[count].second].y - nodes[edges[count].first].y;

            // Compute angle for "pinch point" sorting
            // Pinch points are when the boundary polygon intersects itself, the points will have identical values.
            double angleStartToEnd = NormalizeAngle(std::atan2(dy, dx));
            double angleEndToStart = NormalizeAngle(std::atan2(-dy, -dx));

            // Note: edgesFaces[count][0] is the internal mesh face valid for both directions of boundary traversal
            boundaryAdjacency[edges[count].first].push_back({count, edges[count].second, edgesFaces[count][0], angleStartToEnd});
            boundaryAdjacency[edges[count].second].push_back({count, edges[count].first, edgesFaces[count][0], angleEndToStart});
        }
    }

    auto boundaryAngleLessThan = [](const BoundaryEdge& a, const BoundaryEdge& b)
    { return a.angle < b.angle; };

    // Sort outgoing edges anti-clockwise
    for (auto& [node_id, connected_edges] : boundaryAdjacency)
    {
        std::sort(connected_edges.begin(), connected_edges.end(), boundaryAngleLessThan);
    }

    // Trace boundary polygons
    for (UInt count = 0; count < edges.size(); ++count)
    {
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

            // The smallest visited angle.
            // All edge angles must be in interval [-2pi, 2pi],
            double deltaAngle = 3.0 * std::numbers::pi;
            UInt chosenEdgeIndex = constants::missing::uintValue;

            // Find the boundary edge that tracks closest clockwise to keep empty space on the right
            for (size_t e = 0; e < boundaryEdges.size(); ++e)
            {
                if (edgeVisited[boundaryEdges[e].edgeId])
                {
                    continue;
                }

                double delta = NormalizeAngle(incomingAngle - boundaryEdges[e].angle);

                if (delta < deltaAngle)
                {
                    deltaAngle = delta;
                    chosenEdgeIndex = e;
                }
            }

            // No unvisited edges were found
            // So eigher the boundary loop was completed or a dead-end reached.
            if (chosenEdgeIndex == constants::missing::uintValue)
            {
                break;
            }

            const BoundaryEdge& chosenEdge = boundaryEdges[chosenEdgeIndex];
            edgeVisited[chosenEdge.edgeId] = true;

            // Record the face touched by this next step in the loop
            currentFaces.push_back(chosenEdge.leftFace);

            prevNodeIndex = currentNodeIndex;
            currentNodeIndex = chosenEdge.neighbourNode;
            incomingAngle = chosenEdge.angle;
        }

        if (currentNodes.size() >= 3)
        {
            currentNodes.push_back(currentNodes.front());
            allLoops.push_back(std::move(currentNodes));
            allTouchedFaces.push_back(std::move(currentFaces));
        }
    }
}
