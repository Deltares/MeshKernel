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
}

std::tuple<UInt, UInt> Mesh2DIntersections::GetIntersectionSeed(const Mesh2D& mesh,
                                                                const std::vector<Point>& polyLine,
                                                                const std::vector<bool>& vistedEdges) const
{
    UInt crossedSegmentIndex = constants::missing::uintValue;
    UInt crossedEdgeIndex = constants::missing::uintValue;
    bool isSeedFound = false;

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

        auto intersection = AreSegmentsCrossing(polyLine[firstIndex],
                                                polyLine[secondIndex],
                                                m_mesh.m_nodes[m_mesh.m_edges[edgeIndex].first],
                                                m_mesh.m_nodes[m_mesh.m_edges[edgeIndex].second],
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

void Mesh2DIntersections::Compute(const std::vector<Point>& polyLine)
{
    // 1. Find the intersection of any segment of the polyline with the mesh return if nothing is found
    std::ranges::fill(m_edgesIntersectionsCache, EdgeMeshPolyLineIntersection());
    std::ranges::fill(m_facesIntersectionsCache, FaceMeshPolyLineIntersection());

    const auto polyLineSize = static_cast<UInt>(polyLine.size());

    std::vector<double> cumulativeLength(polyLine.size(), 0.0);
    for (UInt i = 1; i < polyLineSize; ++i)
    {
        cumulativeLength[i] = cumulativeLength[i - 1] + ComputeDistance(polyLine[i], polyLine[i - 1], m_mesh.m_projection);
    }

    std::queue<std::array<UInt, 3>> crossingEdges;
    std::vector<bool> vistedEdges(m_mesh.GetNumEdges(), false);
    std::vector<bool> vistedFace(m_mesh.GetNumEdges(), false);

    // keep traversing the polyline as long crossed edges are found
    while (true)
    {
        // find the seed
        const auto [crossedEdgeIndex, crossedSegmentIndex] = GetIntersectionSeed(
            m_mesh,
            polyLine,
            vistedEdges);

        // no valid seed found in the entire mesh, we are done
        if (crossedEdgeIndex == constants::missing::uintValue)
        {
            break;
        }

        const auto crossingNextSegmentIndex = crossedSegmentIndex + 1;
        crossingEdges.push({crossedEdgeIndex, crossedSegmentIndex, crossingNextSegmentIndex});

        // use breadth search along the current polyline
        while (!crossingEdges.empty())
        {
            auto [currentCrossingEdge, segmentIndex, nextSegmentIndex] = crossingEdges.front();
            crossingEdges.pop();

            for (const auto currentFaceIndex : m_mesh.m_edgesFaces[currentCrossingEdge])
            {
                if (currentFaceIndex == constants::missing::uintValue)
                {
                    continue;
                }

                for (UInt e = 0; e < m_mesh.m_facesEdges[currentFaceIndex].size(); ++e)
                {
                    const auto edgeIndex = m_mesh.m_facesEdges[currentFaceIndex][e];
                    if (vistedEdges[edgeIndex] && vistedFace[currentFaceIndex])
                    {
                        continue;
                    }

                    UInt firstIndex = segmentIndex;
                    UInt secondIndex = nextSegmentIndex;

                    auto intersection = AreSegmentsCrossing(polyLine[firstIndex],
                                                            polyLine[secondIndex],
                                                            m_mesh.m_nodes[m_mesh.m_edges[edgeIndex].first],
                                                            m_mesh.m_nodes[m_mesh.m_edges[edgeIndex].second],
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
                                 adimensionalEdgeDistance) = GetNextEdgeIntersection(polyLine, edgeIndex, segmentIndex, nextSegmentIndex, Direction::Forward);
                    }

                    if (!intersectionFound)
                    {
                        std::tie(intersectionFound,
                                 firstIndex,
                                 secondIndex,
                                 crossProductValue,
                                 adimensionalPolylineSegmentDistance,
                                 adimensionalEdgeDistance) = GetNextEdgeIntersection(polyLine, edgeIndex, segmentIndex, nextSegmentIndex, Direction::Backward);
                    }

                    // none of the polyline intersect the current edge
                    if (!intersectionFound)
                    {
                        continue;
                    }

                    updateEdgeIntersections(
                        firstIndex,
                        edgeIndex,
                        m_mesh.m_edges[edgeIndex],
                        cumulativeLength,
                        crossProductValue,
                        adimensionalEdgeDistance,
                        adimensionalPolylineSegmentDistance,
                        m_edgesIntersectionsCache);

                    updateFaceIntersections(currentFaceIndex, edgeIndex, m_facesIntersectionsCache);

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
