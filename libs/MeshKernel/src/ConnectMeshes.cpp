#include "MeshKernel/ConnectMeshes.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Operations.hpp"
#include "MeshKernel/RangeCheck.hpp"

#include <ranges>

void meshkernel::ConnectMeshes::AreEdgesAdjacent(const Mesh2D& mesh,
                                                 const double separationFraction,
                                                 const UInt edge1,
                                                 const UInt edge2,
                                                 bool& areAdjacent,
                                                 UInt& startNode,
                                                 UInt& endNode) const
{
    const Point edge1Start = mesh.m_nodes[mesh.m_edges[edge1].first];
    const Point edge1End = mesh.m_nodes[mesh.m_edges[edge1].second];
    const Point edge2Start = mesh.m_nodes[mesh.m_edges[edge2].first];
    const Point edge2End = mesh.m_nodes[mesh.m_edges[edge2].second];

    areAdjacent = false;
    startNode = constants::missing::uintValue;
    endNode = constants::missing::uintValue;

    if (edge1Start == edge1End || edge2Start == edge2End)
    {
        return;
    }

    const double edge1Length = ComputeDistance(edge1Start, edge1End, mesh.m_projection);
    const double edge2Length = ComputeDistance(edge2Start, edge2End, mesh.m_projection);
    const double minimumLength = separationFraction * std::min(edge1Length, edge2Length);

    if (edge1Length <= edge2Length)
    {
        const Point midPoint = 0.5 * (edge1Start + edge1End);

        const auto [distance, intersection, parameterisedDistance] = DistanceFromLine(midPoint, edge2Start, edge2End, mesh.m_projection);
        areAdjacent = distance != constants::missing::doubleValue && distance < minimumLength;
    }
    else
    {
        const Point midPoint = 0.5 * (edge2Start + edge2End);

        const auto [distance, intersection, parameterisedDistance] = DistanceFromLine(midPoint, edge1Start, edge1End, mesh.m_projection);
        areAdjacent = distance != constants::missing::doubleValue && distance < minimumLength;
    }

    if (areAdjacent)
    {
        if (ComputeDistance(edge1Start, edge2Start, mesh.m_projection) < minimumLength)
        {
            startNode = mesh.m_edges[edge2].first;
        }
        else if (ComputeDistance(edge1Start, edge2End, mesh.m_projection) < minimumLength)
        {
            startNode = mesh.m_edges[edge2].second;
        }

        if (ComputeDistance(edge1End, edge2Start, mesh.m_projection) < minimumLength)
        {
            endNode = mesh.m_edges[edge2].first;
        }
        else if (ComputeDistance(edge1End, edge2End, mesh.m_projection) < minimumLength)
        {
            endNode = mesh.m_edges[edge2].second;
        }

        double absCosPhi = std::abs(NormalizedInnerProductTwoSegments(edge1Start, edge1End, edge2Start, edge2End, mesh.m_projection));
        areAdjacent = (absCosPhi > 0.9 && absCosPhi <= 1.0);
    }
}

void meshkernel::ConnectMeshes::GetQuadrilateralElementsOnDomainBoundary(const Mesh2D& mesh,
                                                                         std::vector<UInt>& elementsOnDomainBoundary,
                                                                         std::vector<UInt>& edgesOnDomainBoundary) const
{
    for (UInt i = 0; i < mesh.m_edges.size(); ++i)
    {
        if (mesh.m_edgesNumFaces[i] == 1)
        {
            const UInt faceId = mesh.m_edgesFaces[i][0];

            // Only store quadrilateral elements
            if (mesh.m_numFacesNodes[faceId] == 4)
            {
                elementsOnDomainBoundary.push_back(faceId);
                edgesOnDomainBoundary.push_back(i);
            }
        }
    }
}

void meshkernel::ConnectMeshes::GatherHangingNodeIds(const Mesh2D& mesh,
                                                     const double separationFraction,
                                                     const std::vector<UInt>& edgesOnDomainBoundary,
                                                     IrregularEdgeInfoArray& irregularEdges) const
{
    for (UInt i = 0; i < edgesOnDomainBoundary.size(); ++i)
    {
        const UInt edgeI = edgesOnDomainBoundary[i];

        for (UInt j = 0; j < edgesOnDomainBoundary.size(); ++j)
        {
            const UInt edgeJ = edgesOnDomainBoundary[j];

            if (i == j)
            {
                continue;
            }

            UInt startNode;
            UInt endNode;
            bool areAdjacent = false;
            IrregularEdgeInfo& edgeInfo = irregularEdges[i];

            AreEdgesAdjacent(mesh, separationFraction, edgeI, edgeJ, areAdjacent, startNode, endNode);

            if (!areAdjacent)
            {
                continue;
            }

            edgeInfo.hangingNodes[edgeInfo.edgeCount] = j;
            ++edgeInfo.edgeCount;

            if (startNode != constants::missing::uintValue)
            {
                edgeInfo.startNode = startNode;
            }

            if (endNode != constants::missing::uintValue)
            {
                edgeInfo.endNode = endNode;
            }
        }
    }
}

void meshkernel::ConnectMeshes::GatherNodesToMerge(const UInt startNode,
                                                   const UInt endNode,
                                                   const Edge& boundaryEdge,
                                                   std::vector<NodesToMerge>& nodesToMerge,
                                                   std::vector<MergeIndicator>& mergeIndicator) const
{
    if (startNode != constants::missing::uintValue && mergeIndicator[boundaryEdge.first] == MergeIndicator::Initial)
    {
        nodesToMerge.emplace_back(boundaryEdge.first, startNode);
        mergeIndicator[boundaryEdge.first] = MergeIndicator::DoMerge;
    }

    if (endNode != constants::missing::uintValue && mergeIndicator[boundaryEdge.second] == MergeIndicator::Initial)
    {
        nodesToMerge.emplace_back(boundaryEdge.second, endNode);
        mergeIndicator[boundaryEdge.second] = MergeIndicator::DoMerge;
    }
}

void meshkernel::ConnectMeshes::GatherHangingNodes(const UInt primaryStartNode,
                                                   const UInt primaryEndNode,
                                                   const Edge& irregularEdge,
                                                   std::vector<UInt>& hangingNodesOnEdge,
                                                   UInt& numberOfHangingNodes,
                                                   std::vector<MergeIndicator>& mergeIndicator) const
{
    const UInt secondaryStartNode = irregularEdge.first;
    const UInt secondaryEndNode = irregularEdge.second;

    if (mergeIndicator[secondaryStartNode] == MergeIndicator::Initial && secondaryStartNode != primaryStartNode && secondaryStartNode != primaryEndNode)
    {
        hangingNodesOnEdge[numberOfHangingNodes] = secondaryStartNode;
        ++numberOfHangingNodes;
        mergeIndicator[secondaryStartNode] = MergeIndicator::DoNotMerge;
    }

    if (mergeIndicator[secondaryEndNode] == MergeIndicator::Initial && secondaryEndNode != primaryStartNode && secondaryEndNode != primaryEndNode)
    {
        hangingNodesOnEdge[numberOfHangingNodes] = secondaryEndNode;
        ++numberOfHangingNodes;
        mergeIndicator[secondaryEndNode] = MergeIndicator::DoNotMerge;
    }
}

void meshkernel::ConnectMeshes::Compute(Mesh2D& mesh, const double separationFraction) const
{
    // Check that the separationFraction is in correct range (0.0, max]
    range_check::CheckInLeftHalfOpenInterval(separationFraction, 0.0, DefaultMaximumSeparationFraction, "Separation Fraction");

    // Only here to shorten several of the initialisations below
    static constexpr UInt missingValue = constants::missing::uintValue;
    const UInt numberOfEdges = mesh.GetNumEdges();

    // Elements with no neighbour
    std::vector<UInt> elementsOnDomainBoundary;
    elementsOnDomainBoundary.reserve(numberOfEdges);

    // Edges with no neighbour
    std::vector<UInt> edgesOnDomainBoundary;
    edgesOnDomainBoundary.reserve(numberOfEdges);

    GetQuadrilateralElementsOnDomainBoundary(mesh, elementsOnDomainBoundary, edgesOnDomainBoundary);

    IrregularEdgeInfoArray irregularEdges(numberOfEdges);
    GatherHangingNodeIds(mesh, separationFraction, edgesOnDomainBoundary, irregularEdges);

    std::vector<bool> adjacentEdgeIndicator(numberOfEdges, true);
    std::vector<NodesToMerge> nodesToMerge;

    // Size of this array needs to be greater than the number of edges
    // The safety margin here is the maximum number of irregular elements along an edge.
    nodesToMerge.reserve(m_maximumNumberOfIrregularElementsAlongEdge * numberOfEdges);

    std::vector<UInt> hangingNodesOnEdge(numberOfEdges);
    std::vector<MergeIndicator> mergeIndicator(numberOfEdges, MergeIndicator::Initial);

    // Free hanging nodes along edges.
    for (UInt i = 0; i < edgesOnDomainBoundary.size(); ++i)
    {
        const IrregularEdgeInfo& irregularEdge = irregularEdges[i];

        if (irregularEdge.edgeCount == 0)
        {
            continue;
        }

        const UInt boundaryEdgeId = edgesOnDomainBoundary[i];
        const Edge boundaryEdge = mesh.m_edges[boundaryEdgeId];

        const UInt primaryStartNode = irregularEdge.startNode;
        const UInt primaryEndNode = irregularEdge.endNode;

        UInt numberOfHangingNodes = 0;
        std::ranges::fill(hangingNodesOnEdge, missingValue);

        // Find hanging nodes along edge and collect nodes to be merged
        for (UInt j = 0; j < irregularEdge.edgeCount; ++j)
        {
            const UInt irregularEdgeIndex = irregularEdge.hangingNodes[j];
            const UInt irregularEdgeId = edgesOnDomainBoundary[irregularEdgeIndex];

            const IrregularEdgeInfo& otherIrregularEdge = irregularEdges[irregularEdgeIndex];

            if (irregularEdge.edgeCount >= otherIrregularEdge.edgeCount)
            {
                GatherNodesToMerge(primaryStartNode, primaryEndNode, boundaryEdge, nodesToMerge, mergeIndicator);

                if (irregularEdge.edgeCount > otherIrregularEdge.edgeCount)
                {
                    GatherHangingNodes(primaryStartNode, primaryEndNode, mesh.m_edges[irregularEdgeId], hangingNodesOnEdge, numberOfHangingNodes, mergeIndicator);
                }

                if (adjacentEdgeIndicator[boundaryEdgeId] && adjacentEdgeIndicator[irregularEdgeId])
                {
                    adjacentEdgeIndicator[boundaryEdgeId] = false;
                    mesh.m_edges[boundaryEdgeId].first = missingValue;
                    mesh.m_edges[boundaryEdgeId].second = missingValue;
                }
            }
        }

        const UInt boundaryFaceId = elementsOnDomainBoundary[i];
        const Point boundaryNode = mesh.m_nodes[boundaryEdge.first];
        FreeHangingNodes(mesh, numberOfHangingNodes, hangingNodesOnEdge, boundaryFaceId, boundaryEdge, boundaryNode, boundaryEdgeId);
    }

    MergeNodes(mesh, nodesToMerge, mergeIndicator);
    mesh.Administrate();
}

void meshkernel::ConnectMeshes::MergeNodes(Mesh2D& mesh, const std::vector<NodesToMerge>& nodesToMerge, std::vector<MergeIndicator>& mergeIndicator) const
{

    for (const auto& [coincidingNodeFirst, coincidingNodeSecond] : nodesToMerge)
    {
        using enum MergeIndicator;

        if (mergeIndicator[coincidingNodeSecond] != DoNotMerge)
        {
            mesh.MergeTwoNodes(coincidingNodeFirst, coincidingNodeSecond);
            // Set to MergeIndicator::DoNotMerge so it will not be processed again.
            mergeIndicator[coincidingNodeFirst] = DoNotMerge;
            mergeIndicator[coincidingNodeSecond] = DoNotMerge;
        }
    }
}

void meshkernel::ConnectMeshes::GetOrderedDistanceFromPoint(const Mesh2D& mesh,
                                                            const std::vector<UInt>& nodeIndices,
                                                            const UInt numberOfNodes,
                                                            const Point& point,
                                                            BoundedIntegerArray& nearestNeighbours) const
{
    std::array<double, m_maximumNumberOfIrregularElementsAlongEdge> distance;
    BoundedIntegerArray distanceIndex;

    distance.fill(0.0);

    // Compute distances of nodes in list from the base point
    for (UInt j = 0; j < numberOfNodes; ++j)
    {
        distance[j] = ComputeDistance(point, mesh.m_nodes[nodeIndices[j]], mesh.m_projection);
    }

    std::iota(distanceIndex.begin(), distanceIndex.begin() + numberOfNodes, 0);
    std::sort(distanceIndex.begin(), distanceIndex.begin() + numberOfNodes, [&distance = std::as_const(distance)](UInt i, UInt j)
              { return distance[i] < distance[j]; });
    nearestNeighbours.fill(constants::missing::uintValue);

    // Get node indices in order of distance from point, closest first.
    for (UInt j = 0; j < numberOfNodes; ++j)
    {
        nearestNeighbours[j] = nodeIndices[distanceIndex[j]];
    }
}

void meshkernel::ConnectMeshes::FreeOneHangingNode(Mesh2D& mesh, const BoundedIntegerArray& hangingNodes, const UInt startNode, const UInt endNode) const
{
    //
    // 2------+------+
    // |      |      |
    // |      x------+
    // |      |      |
    // 1------+------+
    //

    // Connect node marked with 'x' to nodes labeled 1 and 2
    mesh.ConnectNodes(hangingNodes[0], startNode);
    mesh.ConnectNodes(hangingNodes[0], endNode);
}

void meshkernel::ConnectMeshes::FreeTwoHangingNodes(Mesh2D& mesh,
                                                    const UInt faceId,
                                                    const UInt edgeId,
                                                    const BoundedIntegerArray& hangingNodes,
                                                    const UInt startNode,
                                                    const UInt endNode) const
{
    //
    // 2------4------+------+
    // |      |      |      |
    // |      |      x------+
    // |      o      |      |
    // |      |      x------+
    // |      |      |      |
    // 1------3------+------+
    //

    // Compute point labeled with 'o' in ASCII diagram above
    const Point midPoint = PointAlongLine(mesh.m_nodes[startNode], mesh.m_nodes[endNode], 0.5);
    const UInt newNodeIndex = mesh.InsertNode(midPoint);

    // Connect node marked with 'x' to nodes labeled 3 and 'o'
    mesh.ConnectNodes(hangingNodes[0], newNodeIndex);
    mesh.ConnectNodes(hangingNodes[0], startNode);

    // Connect node marked with 'x' to nodes labeled 'o' and 4
    mesh.ConnectNodes(hangingNodes[1], newNodeIndex);
    mesh.ConnectNodes(hangingNodes[1], endNode);

    // Connect node marked with 'o' to nodes labeled 3 and 4
    mesh.ConnectNodes(newNodeIndex, startNode);
    mesh.ConnectNodes(newNodeIndex, endNode);

    const UInt adjacentFaceId = mesh.NextFace(faceId, edgeId);

    mesh.DeleteEdge(edgeId);

    if (adjacentFaceId != constants::missing::uintValue)
    {
        const UInt nextOppositeEdge = mesh.FindOppositeEdge(adjacentFaceId, edgeId);
        // Connect node marked with 'o' to nodes labeled 1 and 2
        mesh.ConnectNodes(newNodeIndex, mesh.m_edges[nextOppositeEdge].first);
        mesh.ConnectNodes(newNodeIndex, mesh.m_edges[nextOppositeEdge].second);
    }
}

void meshkernel::ConnectMeshes::FreeThreeHangingNodes(Mesh2D& mesh,
                                                      const UInt faceId,
                                                      const UInt edgeId,
                                                      const BoundedIntegerArray& hangingNodes,
                                                      const UInt startNode,
                                                      const UInt endNode) const
{
    //
    // 2------4------+------+
    // |      |      |      |
    // |      |      x------+
    // |      |      |      |
    // |      o      x------+
    // |      |      |      |
    // |      |      x------+
    // |      |      |      |
    // 1------3------+------+
    //

    // Compute point labeled with 'o' in ASCII diagram above
    const Point midPoint = PointAlongLine(mesh.m_nodes[startNode], mesh.m_nodes[endNode], 0.5);
    const UInt newNodeIndex = mesh.InsertNode(midPoint);

    mesh.ConnectNodes(hangingNodes[1], newNodeIndex);
    mesh.ConnectNodes(newNodeIndex, endNode);
    mesh.ConnectNodes(newNodeIndex, startNode);

    mesh.ConnectNodes(hangingNodes[0], newNodeIndex);
    mesh.ConnectNodes(hangingNodes[0], startNode);

    mesh.ConnectNodes(hangingNodes[2], newNodeIndex);
    mesh.ConnectNodes(hangingNodes[2], endNode);

    const UInt adjacentFaceId = mesh.NextFace(faceId, edgeId);

    mesh.DeleteEdge(edgeId);

    if (adjacentFaceId != constants::missing::uintValue)
    {
        const UInt nextOppositeEdge = mesh.FindOppositeEdge(adjacentFaceId, edgeId);
        mesh.ConnectNodes(newNodeIndex, mesh.m_edges[nextOppositeEdge].first);
        mesh.ConnectNodes(newNodeIndex, mesh.m_edges[nextOppositeEdge].second);
    }
}

void meshkernel::ConnectMeshes::FreeFourHangingNodes(Mesh2D& mesh,
                                                     const UInt faceId,
                                                     const UInt edgeId,
                                                     const BoundedIntegerArray& hangingNodes,
                                                     const UInt startNode,
                                                     const UInt endNode) const
{
    //
    //  +------+------+------+------+
    //  |      |      |      |      |
    //  +------x      |      |      |
    //  |      |      3      |      |
    //  +------x      |      5      |
    //  |      |      2      |      6
    //  +------x      |      4      |
    //  |      |      1      |      |
    //  +------x      |      |      |
    //  |      |      |      |      |
    //  +------+------+------+------+
    //            ^       ^      ^
    //          first  second  third
    //           next   next    next
    //

    //--------------------------------
    // Create 3 new nodes (labelled 1, 2 and 3 in ASCII diagram)
    // Connect newly created nodes to hanging nodes

    UInt firstNextFace = mesh.NextFace(faceId, edgeId);

    // Compute points labeled 1, 2 or 3 in ASCII diagram above
    const UInt node1 = mesh.InsertNode(PointAlongLine(mesh.m_nodes[startNode], mesh.m_nodes[endNode], 0.25));
    const UInt node2 = mesh.InsertNode(PointAlongLine(mesh.m_nodes[startNode], mesh.m_nodes[endNode], 0.5));
    const UInt node3 = mesh.InsertNode(PointAlongLine(mesh.m_nodes[startNode], mesh.m_nodes[endNode], 0.75));

    // Connect nodes across the face
    mesh.ConnectNodes(hangingNodes[1], node2);
    mesh.ConnectNodes(hangingNodes[2], node2);

    mesh.ConnectNodes(hangingNodes[1], node1);
    mesh.ConnectNodes(hangingNodes[2], node3);

    mesh.ConnectNodes(hangingNodes[0], node1);
    mesh.ConnectNodes(hangingNodes[3], node3);

    // Connect newly created node along the edge
    mesh.ConnectNodes(startNode, node1);
    mesh.ConnectNodes(node1, node2);
    mesh.ConnectNodes(node2, node3);
    mesh.ConnectNodes(node3, endNode);

    // The original edge can now be deleted.
    mesh.DeleteEdge(edgeId);

    // Reached the end of the mesh
    if (firstNextFace == constants::missing::uintValue)
    {
        return;
    }

    const UInt firstNextOppositeEdge = mesh.FindOppositeEdge(firstNextFace, edgeId);
    const UInt firstNextOppositeStartNode = mesh.m_edges[firstNextOppositeEdge].first;
    const UInt firstNextOppositeEndNode = mesh.m_edges[firstNextOppositeEdge].second;

    //--------------------------------
    // Create 2 new nodes (labelled 4 and 5 in ASCII diagram)
    // Connect newly created nodes to (newly created) hanging nodes (1, 2 and 3)

    // Compute points labeled with 4 or 5 in ASCII diagram above
    UInt node4 = mesh.InsertNode(PointAlongLine(mesh.m_nodes[firstNextOppositeStartNode], mesh.m_nodes[firstNextOppositeEndNode], 0.34));
    UInt node5 = mesh.InsertNode(PointAlongLine(mesh.m_nodes[firstNextOppositeStartNode], mesh.m_nodes[firstNextOppositeEndNode], 0.66));

    // Connect nodes across the face
    mesh.ConnectNodes(node1, node4);
    mesh.ConnectNodes(node4, node2);
    mesh.ConnectNodes(node2, node5);
    mesh.ConnectNodes(node5, node3);

    // Connect newly created node along the edge
    mesh.ConnectNodes(firstNextOppositeStartNode, node4);
    mesh.ConnectNodes(node4, node5);
    mesh.ConnectNodes(node5, firstNextOppositeEndNode);

    // The original edge can now be deleted.
    mesh.DeleteEdge(firstNextOppositeEdge);

    const UInt secondNextFaceId = mesh.NextFace(firstNextFace, firstNextOppositeEdge);

    // Reached the end of the mesh
    if (secondNextFaceId == constants::missing::uintValue)
    {
        return;
    }

    const UInt secondNextOppositeEdge = mesh.FindOppositeEdge(secondNextFaceId, firstNextOppositeEdge);
    const UInt secondNextOppositeStartNode = mesh.m_edges[secondNextOppositeEdge].first;
    const UInt secondNextOppositeEndNode = mesh.m_edges[secondNextOppositeEdge].second;

    // Compute point labeled with 6 in ASCII diagram above
    const UInt node6 = mesh.InsertNode(PointAlongLine(mesh.m_nodes[secondNextOppositeStartNode], mesh.m_nodes[secondNextOppositeEndNode], 0.5));

    // Connect nodes across the face1
    mesh.ConnectNodes(node4, node6);
    mesh.ConnectNodes(node5, node6);

    // Connect newly created node along the edge
    mesh.ConnectNodes(secondNextOppositeStartNode, node6);
    mesh.ConnectNodes(secondNextOppositeEndNode, node6);

    // The original edge can now be deleted.
    mesh.DeleteEdge(secondNextOppositeEdge);

    //--------------------------------
    // Create 1 new node (labelled 6 in ASCII diagram)
    // Connect newly created node to existing nodes

    const UInt thirdNextFaceId = mesh.NextFace(secondNextFaceId, secondNextOppositeEdge);

    // Reached the end of the mesh
    if (thirdNextFaceId == constants::missing::uintValue)
    {
        return;
    }

    const UInt thirdNextOppositeEdge = mesh.FindOppositeEdge(thirdNextFaceId, secondNextOppositeEdge);
    mesh.ConnectNodes(node6, mesh.m_edges[thirdNextOppositeEdge].second);
    mesh.ConnectNodes(node6, mesh.m_edges[thirdNextOppositeEdge].first);
}

void meshkernel::ConnectMeshes::FreeHangingNodes(Mesh2D& mesh,
                                                 const UInt numberOfHangingNodes,
                                                 const std::vector<UInt>& hangingNodesOnEdge,
                                                 const UInt faceId,
                                                 const Edge& boundaryEdge,
                                                 const Point& boundaryNode,
                                                 const UInt edgeId) const
{
    if (numberOfHangingNodes == 0)
    {
        return;
    }

    BoundedIntegerArray hangingNodes;
    GetOrderedDistanceFromPoint(mesh, hangingNodesOnEdge, numberOfHangingNodes, boundaryNode, hangingNodes);

    UInt oppositeEdgeId = mesh.FindOppositeEdge(faceId, edgeId);

    if (oppositeEdgeId == constants::missing::uintValue)
    {
        throw ConstraintError("Opposite edge not found for element {} and edge {}.", faceId, edgeId);
    }

    UInt startNode = mesh.m_edges[oppositeEdgeId].first;
    UInt endNode = mesh.m_edges[oppositeEdgeId].second;

    const auto [segmentsCross,
                intersectionPoint,
                crossProduct,
                normalisedPolylineSegmentDistance,
                normalisedEdgeDistance] = AreSegmentsCrossing(mesh.m_nodes[boundaryEdge.first],
                                                              mesh.m_nodes[startNode],
                                                              mesh.m_nodes[boundaryEdge.second],
                                                              mesh.m_nodes[endNode],
                                                              false,
                                                              mesh.m_projection);
    if (segmentsCross)
    {
        std::swap(startNode, endNode);
    }

    switch (numberOfHangingNodes)
    {
    case 1:
        FreeOneHangingNode(mesh, hangingNodes, startNode, endNode);
        break;
    case 2:
        FreeTwoHangingNodes(mesh, faceId, oppositeEdgeId, hangingNodes, startNode, endNode);
        break;
    case 3:
        FreeThreeHangingNodes(mesh, faceId, oppositeEdgeId, hangingNodes, startNode, endNode);
        break;
    case 4:
        FreeFourHangingNodes(mesh, faceId, oppositeEdgeId, hangingNodes, startNode, endNode);
        break;
    default:
        // 0 hanging nodes is handled at the top of this function, so to be here can only be: numberOfHangingNodes > 4
        throw NotImplementedError("Cannot handle more than 4 hanging nodes along irregular edge, number is: {}",
                                  numberOfHangingNodes);
    }
}
