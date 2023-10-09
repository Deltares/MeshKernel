#include <ranges>

#include "MeshKernel/Mesh2DIntersections.hpp"

#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Operations.hpp"

using namespace meshkernel;

std::tuple<UInt, UInt> Mesh2DIntersections::GetIntersectionSeed(const Mesh2D& mesh,
                                                                const std::vector<Point>& polyLine,
                                                                const std::vector<bool>& vistedEdges) const
{
    UInt crossedSegmentIndex = constants::missing::uintValue;
    UInt crossedEdgeIndex = constants::missing::uintValue;

    // Find starting edge and segment
    for (UInt segmentIndex = 0; segmentIndex < polyLine.size() - 1; ++segmentIndex)
    {
        for (UInt edgeIndex = 0; edgeIndex < mesh.GetNumEdges(); ++edgeIndex)
        {
            // edge already crossed, nothing to do
            if (vistedEdges[edgeIndex])
            {
                continue;
            }

            const auto [isEdgeCrossed,
                        intersectionPoint,
                        crossProductValue,
                        adimensionalPolylineSegmentDistance,
                        adimensionalEdgeDistance] = AreSegmentsCrossing(polyLine[segmentIndex],
                                                                        polyLine[segmentIndex + 1],
                                                                        mesh.m_nodes[mesh.m_edges[edgeIndex].first],
                                                                        mesh.m_nodes[mesh.m_edges[edgeIndex].second],
                                                                        false,
                                                                        mesh.m_projection);
            if (!isEdgeCrossed)
            {
                continue;
            }

            crossedSegmentIndex = segmentIndex;
            crossedEdgeIndex = edgeIndex;
            break;
        }

        if (crossedSegmentIndex != constants::missing::uintValue)
        {
            break;
        }
    }

    return {crossedEdgeIndex, crossedSegmentIndex};
}

void Mesh2DIntersections::GetPolylineIntersection(const Mesh2D& mesh, const std::vector<Point>& polyLine)
{
    // 1. Find the intersection of any segment of the polyline with the mesh return if nothing is found
    std::ranges::fill(m_edgesIntersectionsCache, EdgeMeshPolylineIntersection());
    std::ranges::fill(m_facesIntersectionsCache, FaceMeshPolylineIntersection());

    const auto polylineSize = static_cast<UInt>(polyLine.size());

    std::vector<double> cumulativeLength(polyLine.size(), 0.0);
    for (UInt i = 1; i < polylineSize; ++i)
    {
        cumulativeLength[i] = cumulativeLength[i - 1] + ComputeDistance(polyLine[i], polyLine[i - 1], mesh.m_projection);
    }

    std::queue<std::array<UInt, 3>> crossingEdges;
    std::vector<bool> vistedEdges(mesh.GetNumEdges(), false);
    std::vector<bool> vistedFace(mesh.GetNumEdges(), false);

    // keep traversing the polyline as long crossed edges are found
    while (true)
    {
        // find a crossed edge on a non-visited segment
        const auto intersectionSeed = GetIntersectionSeed(mesh,
                                                          polyLine,
                                                          vistedEdges);

        const auto crossedEdgeIndex = std::get<0>(intersectionSeed);

        // no valid seed found in the entire mesh, we are done
        if (crossedEdgeIndex == constants::missing::uintValue)
        {
            break;
        }

        const auto crossingSegmentIndex = std::get<1>(intersectionSeed);
        const auto crossingNextSegmentIndex = crossingSegmentIndex + 1;
        crossingEdges.push({crossedEdgeIndex, crossingSegmentIndex, crossingNextSegmentIndex});

        // use breadth search along the current polyline
        while (!crossingEdges.empty())
        {
            auto [currentCrossingEdge, segmentIndex, nextSegmentIndex] = crossingEdges.front();
            crossingEdges.pop();

            for (UInt f = 0; f < mesh.m_edgesFaces[currentCrossingEdge].size(); ++f)
            {
                const auto currentFaceIndex = mesh.m_edgesFaces[currentCrossingEdge][f];

                if (currentFaceIndex == constants::missing::uintValue)
                {
                    continue;
                }

                for (UInt e = 0; e < mesh.m_facesEdges[currentFaceIndex].size(); ++e)
                {
                    const auto edgeIndex = mesh.m_facesEdges[currentFaceIndex][e];
                    if (vistedEdges[edgeIndex] && vistedFace[currentFaceIndex])
                    {
                        continue;
                    }

                    UInt firstIndex = segmentIndex;
                    UInt secondIndex = nextSegmentIndex;

                    auto intersection = AreSegmentsCrossing(polyLine[firstIndex],
                                                            polyLine[secondIndex],
                                                            mesh.m_nodes[mesh.m_edges[edgeIndex].first],
                                                            mesh.m_nodes[mesh.m_edges[edgeIndex].second],
                                                            false,
                                                            mesh.m_projection);

                    auto intersectionFound = std::get<0>(intersection);
                    if (!intersectionFound)
                    {
                        UInt numForwardSteps = 0;
                        while (!intersectionFound && firstIndex < polylineSize - 2 && numForwardSteps < mesh.m_maxSteps)
                        {
                            firstIndex = secondIndex;
                            secondIndex = firstIndex + 1;
                            intersection = AreSegmentsCrossing(polyLine[firstIndex],
                                                               polyLine[secondIndex],
                                                               mesh.m_nodes[mesh.m_edges[edgeIndex].first],
                                                               mesh.m_nodes[mesh.m_edges[edgeIndex].second],
                                                               false,
                                                               mesh.m_projection);

                            intersectionFound = std::get<0>(intersection);
                            numForwardSteps++;
                        }
                    }

                    if (!intersectionFound)
                    {
                        firstIndex = segmentIndex;
                        secondIndex = nextSegmentIndex;
                        UInt numBackwardSteps = 0;
                        while (!intersectionFound && firstIndex >= 1 && numBackwardSteps < mesh.m_maxSteps)
                        {
                            secondIndex = firstIndex;
                            firstIndex = firstIndex - 1;
                            intersection = AreSegmentsCrossing(polyLine[firstIndex],
                                                               polyLine[secondIndex],
                                                               mesh.m_nodes[mesh.m_edges[edgeIndex].first],
                                                               mesh.m_nodes[mesh.m_edges[edgeIndex].second],
                                                               false,
                                                               mesh.m_projection);

                            intersectionFound = std::get<0>(intersection);
                            numBackwardSteps++;
                        }
                    }

                    // none of the polyline intersect the current edge
                    if (!intersectionFound)
                    {
                        continue;
                    }

                    const double crossProductValue = std::get<2>(intersection);
                    const double adimensionalPolylineSegmentDistance = std::get<3>(intersection);
                    const double adimensionalEdgeDistance = std::get<4>(intersection);

                    EdgeMeshPolylineIntersection::updateIntersections(
                        firstIndex,
                        edgeIndex,
                        mesh.m_edges[edgeIndex].first,
                        mesh.m_edges[edgeIndex].second,
                        cumulativeLength,
                        crossProductValue,
                        adimensionalEdgeDistance,
                        adimensionalPolylineSegmentDistance,
                        m_edgesIntersectionsCache);

                    FaceMeshPolylineIntersection::updateIntersections(currentFaceIndex, edgeIndex, m_facesIntersectionsCache);

                    if (edgeIndex != currentCrossingEdge)
                    {
                        crossingEdges.push({edgeIndex, firstIndex, secondIndex});
                    }
                    vistedEdges[edgeIndex] = true;
                }
                vistedFace[currentFaceIndex] = true;
            }
        }
    }

    // compute polylineDistance, sort the edges for each face
    for (auto& facesIntersection : m_facesIntersectionsCache)
    {
        if (facesIntersection.edgeIndexses.empty())
        {
            continue;
        }

        if (facesIntersection.edgeIndexses.size() == 1)
        {
            const auto edgeIndex = facesIntersection.edgeIndexses[0];
            facesIntersection.polylineDistance = m_edgesIntersectionsCache[edgeIndex].polylineDistance;
        }

        if (facesIntersection.edgeIndexses.size() == 2)
        {
            const auto firstEdgeIndex = facesIntersection.edgeIndexses[0];
            const auto secondEdgeIndex = facesIntersection.edgeIndexses[1];

            // swap the edge indexes if needed
            if (m_edgesIntersectionsCache[firstEdgeIndex].adimensionalPolylineSegmentDistance > m_edgesIntersectionsCache[secondEdgeIndex].adimensionalPolylineSegmentDistance)
            {
                std::swap(facesIntersection.edgeIndexses[0], facesIntersection.edgeIndexses[1]);
            }

            // compute the polylineDistance for the face
            facesIntersection.polylineDistance = 0.5 * (m_edgesIntersectionsCache[firstEdgeIndex].polylineDistance + m_edgesIntersectionsCache[secondEdgeIndex].polylineDistance);
        }

        // push back the face intersection edge nodes
        for (UInt e = 0; e < facesIntersection.edgeIndexses.size(); ++e)
        {
            const auto edgeIndex = facesIntersection.edgeIndexses[e];
            facesIntersection.edgeNodes.emplace_back(m_edgesIntersectionsCache[edgeIndex].edgeFirstNode);
            facesIntersection.edgeNodes.emplace_back(m_edgesIntersectionsCache[edgeIndex].edgeSecondNode);
        }
    }

    // edge intersections are unique
    for (UInt e = 0; e < mesh.GetNumEdges(); ++e)
    {
        if (m_edgesIntersectionsResult[e].polylineDistance < 0)
        {
            m_edgesIntersectionsResult[e] = m_edgesIntersectionsCache[e];
        }
    }

    // face intersections are not unique and a face could have been intersected already
    for (UInt f = 0; f < mesh.GetNumFaces(); ++f)
    {
        if (!m_faceIntersectionsResult[f].edgeNodes.empty() &&
            !m_facesIntersectionsCache[f].edgeNodes.empty())
        {
            m_faceIntersectionsResult[f].edgeIndexses.insert(m_faceIntersectionsResult[f].edgeIndexses.end(), m_facesIntersectionsCache[f].edgeIndexses.begin(), m_facesIntersectionsCache[f].edgeIndexses.end());
            m_faceIntersectionsResult[f].edgeNodes.insert(m_faceIntersectionsResult[f].edgeNodes.end(), m_facesIntersectionsCache[f].edgeNodes.begin(), m_facesIntersectionsCache[f].edgeNodes.end());
            m_faceIntersectionsResult[f].polylineDistance = 0.5 * (m_faceIntersectionsResult[f].polylineDistance + m_facesIntersectionsCache[f].polylineDistance);
        }
        else if (!m_facesIntersectionsCache[f].edgeNodes.empty())
        {
            m_faceIntersectionsResult[f] = m_facesIntersectionsCache[f];
        }
    }
}

void Mesh2DIntersections::Compute(Mesh2D& mesh, const Polygons& polygon)
{
    // No polygon, nothing is crossed
    if (polygon.IsEmpty())
    {
        throw AlgorithmError("No polygon provided!");
    }

    // Intersection results
    m_edgesIntersectionsResult.resize(mesh.GetNumEdges());
    m_faceIntersectionsResult.resize(mesh.GetNumFaces());

    // Declare local caches
    m_edgesIntersectionsCache.resize(mesh.GetNumEdges());
    m_facesIntersectionsCache.resize(mesh.GetNumFaces());

    // Make sure face information is available
    mesh.Administrate();

    // Multiple polygons here
    for (auto outer = 0u; outer < polygon.GetNumPolygons(); ++outer)
    {
        std::vector<std::vector<Point>> allPolylines;

        allPolylines.emplace_back(polygon.Enclosure(outer).Outer().Nodes());

        for (auto inner = 0u; inner < polygon.Enclosure(outer).NumberOfInner(); ++inner)
        {
            allPolylines.emplace_back(polygon.Enclosure(outer).Inner(inner).Nodes());
        }

        for (const auto& polyLine : allPolylines)
        {
            GetPolylineIntersection(mesh, polyLine);
        }
    }
}
