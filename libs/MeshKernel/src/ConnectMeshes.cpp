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
                                                 UInt& endNode)
{
    const Point edge1Start = mesh.Node(mesh.GetEdge(edge1).first);
    const Point edge1End = mesh.Node(mesh.GetEdge(edge1).second);
    const Point edge2Start = mesh.Node(mesh.GetEdge(edge2).first);
    const Point edge2End = mesh.Node(mesh.GetEdge(edge2).second);

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
            startNode = mesh.GetEdge(edge2).first;
        }
        else if (ComputeDistance(edge1Start, edge2End, mesh.m_projection) < minimumLength)
        {
            startNode = mesh.GetEdge(edge2).second;
        }

        if (ComputeDistance(edge1End, edge2Start, mesh.m_projection) < minimumLength)
        {
            endNode = mesh.GetEdge(edge2).first;
        }
        else if (ComputeDistance(edge1End, edge2End, mesh.m_projection) < minimumLength)
        {
            endNode = mesh.GetEdge(edge2).second;
        }

        double absCosPhi = std::abs(NormalizedInnerProductTwoSegments(edge1Start, edge1End, edge2Start, edge2End, mesh.m_projection));
        areAdjacent = (absCosPhi > minimumParallelDeviation);
    }
}

void meshkernel::ConnectMeshes::GetQuadrilateralElementsOnDomainBoundary(const Mesh2D& mesh,
                                                                         std::vector<UInt>& elementsOnDomainBoundary,
                                                                         std::vector<UInt>& edgesOnDomainBoundary)
{
    for (UInt i = 0; i < mesh.GetNumEdges(); ++i)
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
                                                     IrregularEdgeInfoArray& irregularEdges)
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
                                                   std::vector<MergeIndicator>& mergeIndicator)
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
                                                   std::vector<MergeIndicator>& mergeIndicator)
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

std::unique_ptr<meshkernel::UndoAction> meshkernel::ConnectMeshes::Compute(Mesh2D& mesh, const double separationFraction)
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

    std::unique_ptr<meshkernel::CompoundUndoAction> conectMeshesAction = CompoundUndoAction::Create();

    // Free hanging nodes along edges.
    for (UInt i = 0; i < edgesOnDomainBoundary.size(); ++i)
    {
        const IrregularEdgeInfo& irregularEdge = irregularEdges[i];

        if (irregularEdge.edgeCount == 0)
        {
            continue;
        }

        const UInt boundaryEdgeId = edgesOnDomainBoundary[i];
        const Edge boundaryEdge = mesh.GetEdge(boundaryEdgeId);

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
                    GatherHangingNodes(primaryStartNode, primaryEndNode, mesh.GetEdge(irregularEdgeId), hangingNodesOnEdge, numberOfHangingNodes, mergeIndicator);
                }

                if (adjacentEdgeIndicator[boundaryEdgeId] && adjacentEdgeIndicator[irregularEdgeId])
                {
                    adjacentEdgeIndicator[boundaryEdgeId] = false;
                    // TODO should return action
                    mesh.SetEdge(boundaryEdgeId, {missingValue, missingValue});
                }
            }
        }

        const UInt boundaryFaceId = elementsOnDomainBoundary[i];
        const Point boundaryNode = mesh.Node(boundaryEdge.first);
        conectMeshesAction->Add(FreeHangingNodes(mesh, numberOfHangingNodes, hangingNodesOnEdge, boundaryFaceId, boundaryEdge, boundaryNode, boundaryEdgeId));
    }

    conectMeshesAction->Add(MergeNodes(mesh, nodesToMerge, mergeIndicator));
    mesh.Administrate();
    return conectMeshesAction;
}

std::unique_ptr<meshkernel::UndoAction> meshkernel::ConnectMeshes::MergeNodes(Mesh2D& mesh, const std::vector<NodesToMerge>& nodesToMerge, std::vector<MergeIndicator>& mergeIndicator)
{
    std::unique_ptr<meshkernel::CompoundUndoAction> mergeAction = CompoundUndoAction::Create();

    for (const auto& [coincidingNodeFirst, coincidingNodeSecond] : nodesToMerge)
    {
        using enum MergeIndicator;

        if (mergeIndicator[coincidingNodeSecond] != DoNotMerge)
        {
            mergeAction->Add(mesh.MergeTwoNodes(coincidingNodeFirst, coincidingNodeSecond));
            // Set to MergeIndicator::DoNotMerge so it will not be processed again.
            mergeIndicator[coincidingNodeFirst] = DoNotMerge;
            mergeIndicator[coincidingNodeSecond] = DoNotMerge;
        }
    }

    return mergeAction;
}

void meshkernel::ConnectMeshes::GetOrderedDistanceFromPoint(const Mesh2D& mesh,
                                                            const std::vector<UInt>& nodeIndices,
                                                            const UInt numberOfNodes,
                                                            const Point& point,
                                                            BoundedIntegerArray& nearestNeighbours)
{
    std::array<double, m_maximumNumberOfIrregularElementsAlongEdge> distance;
    BoundedIntegerArray distanceIndex;

    distance.fill(0.0);

    // Compute distances of nodes in list from the base point
    for (UInt j = 0; j < numberOfNodes; ++j)
    {
        distance[j] = ComputeDistance(point, mesh.Node(nodeIndices[j]), mesh.m_projection);
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

std::unique_ptr<meshkernel::UndoAction> meshkernel::ConnectMeshes::FreeOneHangingNode(Mesh2D& mesh, const BoundedIntegerArray& hangingNodes, const UInt startNode, const UInt endNode)
{
    //
    // 2------+------+
    // |      |      |
    // |      x------+
    // |      |      |
    // 1------+------+
    //

    std::unique_ptr<CompoundUndoAction> freeHangingNodesAction = CompoundUndoAction::Create();
    UInt edgeId;
    std::unique_ptr<AddEdgeAction> addEdgeAction;

    // Connect node marked with 'x' to nodes labeled 1 and 2
    std::tie(edgeId, addEdgeAction) = mesh.ConnectNodes(hangingNodes[0], startNode);
    freeHangingNodesAction->Add(std::move(addEdgeAction));

    std::tie(edgeId, addEdgeAction) = mesh.ConnectNodes(hangingNodes[0], endNode);
    freeHangingNodesAction->Add(std::move(addEdgeAction));

    return freeHangingNodesAction;
}

std::unique_ptr<meshkernel::UndoAction> meshkernel::ConnectMeshes::FreeTwoHangingNodes(Mesh2D& mesh,
                                                                                       const UInt faceId,
                                                                                       const UInt edgeId,
                                                                                       const BoundedIntegerArray& hangingNodes,
                                                                                       const UInt startNode,
                                                                                       const UInt endNode)
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

    std::unique_ptr<CompoundUndoAction> freeHangingNodesAction = CompoundUndoAction::Create();
    UInt newEdgeId;
    std::unique_ptr<AddEdgeAction> addEdgeAction;

    // Compute point labeled with 'o' in ASCII diagram above
    const Point midPoint = PointAlongLine(mesh.Node(startNode), mesh.Node(endNode), 0.5);
    auto [newNodeIndex, addNodeAction] = mesh.InsertNode(midPoint);
    freeHangingNodesAction->Add(std::move(addNodeAction));

    // Connect node marked with 'x' to nodes labeled 3 and 'o'
    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(hangingNodes[0], newNodeIndex);
    freeHangingNodesAction->Add(std::move(addEdgeAction));

    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(hangingNodes[0], startNode);
    freeHangingNodesAction->Add(std::move(addEdgeAction));

    // Connect node marked with 'x' to nodes labeled 'o' and 4
    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(hangingNodes[1], newNodeIndex);
    freeHangingNodesAction->Add(std::move(addEdgeAction));

    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(hangingNodes[1], endNode);
    freeHangingNodesAction->Add(std::move(addEdgeAction));

    // Connect node marked with 'o' to nodes labeled 3 and 4
    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(newNodeIndex, startNode);
    freeHangingNodesAction->Add(std::move(addEdgeAction));

    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(newNodeIndex, endNode);
    freeHangingNodesAction->Add(std::move(addEdgeAction));

    const UInt adjacentFaceId = mesh.NextFace(faceId, edgeId);

    freeHangingNodesAction->Add(mesh.DeleteEdge(edgeId));

    if (adjacentFaceId != constants::missing::uintValue)
    {
        const UInt nextOppositeEdge = mesh.FindOppositeEdge(adjacentFaceId, edgeId);
        // Connect node marked with 'o' to nodes labeled 1 and 2
        std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(newNodeIndex, mesh.GetEdge(nextOppositeEdge).first);
        freeHangingNodesAction->Add(std::move(addEdgeAction));

        std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(newNodeIndex, mesh.GetEdge(nextOppositeEdge).second);
        freeHangingNodesAction->Add(std::move(addEdgeAction));
    }

    return freeHangingNodesAction;
}

std::unique_ptr<meshkernel::UndoAction> meshkernel::ConnectMeshes::FreeThreeHangingNodes(Mesh2D& mesh,
                                                                                         const UInt faceId,
                                                                                         const UInt edgeId,
                                                                                         const BoundedIntegerArray& hangingNodes,
                                                                                         const UInt startNode,
                                                                                         const UInt endNode)
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

    std::unique_ptr<CompoundUndoAction> freeHangingNodesAction = CompoundUndoAction::Create();
    UInt newEdgeId;
    std::unique_ptr<AddEdgeAction> addEdgeAction;

    // Compute point labeled with 'o' in ASCII diagram above
    const Point midPoint = PointAlongLine(mesh.Node(startNode), mesh.Node(endNode), 0.5);
    auto [newNodeIndex, addNodeAction] = mesh.InsertNode(midPoint);
    freeHangingNodesAction->Add(std::move(addNodeAction));

    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(hangingNodes[1], newNodeIndex);
    freeHangingNodesAction->Add(std::move(addEdgeAction));

    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(newNodeIndex, endNode);
    freeHangingNodesAction->Add(std::move(addEdgeAction));

    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(newNodeIndex, startNode);
    freeHangingNodesAction->Add(std::move(addEdgeAction));

    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(hangingNodes[0], newNodeIndex);
    freeHangingNodesAction->Add(std::move(addEdgeAction));

    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(hangingNodes[0], startNode);
    freeHangingNodesAction->Add(std::move(addEdgeAction));

    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(hangingNodes[2], newNodeIndex);
    freeHangingNodesAction->Add(std::move(addEdgeAction));

    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(hangingNodes[2], endNode);
    freeHangingNodesAction->Add(std::move(addEdgeAction));

    const UInt adjacentFaceId = mesh.NextFace(faceId, edgeId);

    freeHangingNodesAction->Add(mesh.DeleteEdge(edgeId));

    if (adjacentFaceId != constants::missing::uintValue)
    {
        const UInt nextOppositeEdge = mesh.FindOppositeEdge(adjacentFaceId, edgeId);
        std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(newNodeIndex, mesh.GetEdge(nextOppositeEdge).first);
        freeHangingNodesAction->Add(std::move(addEdgeAction));

        std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(newNodeIndex, mesh.GetEdge(nextOppositeEdge).second);
        freeHangingNodesAction->Add(std::move(addEdgeAction));
    }

    return freeHangingNodesAction;
}

std::unique_ptr<meshkernel::UndoAction> meshkernel::ConnectMeshes::FreeFourHangingNodes(Mesh2D& mesh,
                                                                                        const UInt faceId,
                                                                                        const UInt edgeId,
                                                                                        const BoundedIntegerArray& hangingNodes,
                                                                                        const UInt startNode,
                                                                                        const UInt endNode)
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

    std::unique_ptr<CompoundUndoAction> freeHangingNodesAction = CompoundUndoAction::Create();
    UInt newEdgeId;
    std::unique_ptr<AddEdgeAction> addEdgeAction;

    UInt firstNextFace = mesh.NextFace(faceId, edgeId);

    // Compute points labeled 1, 2 or 3 in ASCII diagram above
    auto [node1, node1Insertion] = mesh.InsertNode(PointAlongLine(mesh.Node(startNode), mesh.Node(endNode), 0.25));
    freeHangingNodesAction->Add(std::move(node1Insertion));
    auto [node2, node2Insertion] = mesh.InsertNode(PointAlongLine(mesh.Node(startNode), mesh.Node(endNode), 0.5));
    freeHangingNodesAction->Add(std::move(node2Insertion));
    auto [node3, node3Insertion] = mesh.InsertNode(PointAlongLine(mesh.Node(startNode), mesh.Node(endNode), 0.75));
    freeHangingNodesAction->Add(std::move(node3Insertion));

    // Connect nodes across the face
    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(hangingNodes[1], node2);
    freeHangingNodesAction->Add(std::move(addEdgeAction));
    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(hangingNodes[2], node2);
    freeHangingNodesAction->Add(std::move(addEdgeAction));

    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(hangingNodes[1], node1);
    freeHangingNodesAction->Add(std::move(addEdgeAction));
    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(hangingNodes[2], node3);
    freeHangingNodesAction->Add(std::move(addEdgeAction));

    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(hangingNodes[0], node1);
    freeHangingNodesAction->Add(std::move(addEdgeAction));
    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(hangingNodes[3], node3);
    freeHangingNodesAction->Add(std::move(addEdgeAction));

    // Connect newly created node along the edge
    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(startNode, node1);
    freeHangingNodesAction->Add(std::move(addEdgeAction));
    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(node1, node2);
    freeHangingNodesAction->Add(std::move(addEdgeAction));
    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(node2, node3);
    freeHangingNodesAction->Add(std::move(addEdgeAction));
    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(node3, endNode);
    freeHangingNodesAction->Add(std::move(addEdgeAction));

    // The original edge can now be deleted.
    freeHangingNodesAction->Add(mesh.DeleteEdge(edgeId));

    // Reached the end of the mesh
    if (firstNextFace == constants::missing::uintValue)
    {
        return freeHangingNodesAction;
    }

    const UInt firstNextOppositeEdge = mesh.FindOppositeEdge(firstNextFace, edgeId);
    const UInt firstNextOppositeStartNode = mesh.GetEdge(firstNextOppositeEdge).first;
    const UInt firstNextOppositeEndNode = mesh.GetEdge(firstNextOppositeEdge).second;

    //--------------------------------
    // Create 2 new nodes (labelled 4 and 5 in ASCII diagram)
    // Connect newly created nodes to (newly created) hanging nodes (1, 2 and 3)

    // Compute points labeled with 4 or 5 in ASCII diagram above
    auto [node4, node4Insertion] = mesh.InsertNode(PointAlongLine(mesh.Node(firstNextOppositeStartNode), mesh.Node(firstNextOppositeEndNode), 0.34));
    freeHangingNodesAction->Add(std::move(node4Insertion));
    auto [node5, node5Insertion] = mesh.InsertNode(PointAlongLine(mesh.Node(firstNextOppositeStartNode), mesh.Node(firstNextOppositeEndNode), 0.66));
    freeHangingNodesAction->Add(std::move(node5Insertion));

    // Connect nodes across the face
    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(node1, node4);
    freeHangingNodesAction->Add(std::move(addEdgeAction));
    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(node4, node2);
    freeHangingNodesAction->Add(std::move(addEdgeAction));
    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(node2, node5);
    freeHangingNodesAction->Add(std::move(addEdgeAction));
    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(node5, node3);
    freeHangingNodesAction->Add(std::move(addEdgeAction));

    // Connect newly created node along the edge
    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(firstNextOppositeStartNode, node4);
    freeHangingNodesAction->Add(std::move(addEdgeAction));
    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(node4, node5);
    freeHangingNodesAction->Add(std::move(addEdgeAction));
    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(node5, firstNextOppositeEndNode);
    freeHangingNodesAction->Add(std::move(addEdgeAction));

    // The original edge can now be deleted.
    freeHangingNodesAction->Add(mesh.DeleteEdge(firstNextOppositeEdge));

    const UInt secondNextFaceId = mesh.NextFace(firstNextFace, firstNextOppositeEdge);

    // Reached the end of the mesh
    if (secondNextFaceId == constants::missing::uintValue)
    {
        return freeHangingNodesAction;
    }

    const UInt secondNextOppositeEdge = mesh.FindOppositeEdge(secondNextFaceId, firstNextOppositeEdge);
    const UInt secondNextOppositeStartNode = mesh.GetEdge(secondNextOppositeEdge).first;
    const UInt secondNextOppositeEndNode = mesh.GetEdge(secondNextOppositeEdge).second;

    // Compute point labeled with 6 in ASCII diagram above
    auto [node6, node6Insertion] = mesh.InsertNode(PointAlongLine(mesh.Node(secondNextOppositeStartNode), mesh.Node(secondNextOppositeEndNode), 0.5));
    freeHangingNodesAction->Add(std::move(node6Insertion));

    // Connect nodes across the face1
    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(node4, node6);
    freeHangingNodesAction->Add(std::move(addEdgeAction));
    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(node5, node6);
    freeHangingNodesAction->Add(std::move(addEdgeAction));

    // Connect newly created node along the edge
    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(secondNextOppositeStartNode, node6);
    freeHangingNodesAction->Add(std::move(addEdgeAction));
    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(secondNextOppositeEndNode, node6);
    freeHangingNodesAction->Add(std::move(addEdgeAction));

    // The original edge can now be deleted.
    freeHangingNodesAction->Add(mesh.DeleteEdge(secondNextOppositeEdge));

    //--------------------------------
    // Create 1 new node (labelled 6 in ASCII diagram)
    // Connect newly created node to existing nodes

    const UInt thirdNextFaceId = mesh.NextFace(secondNextFaceId, secondNextOppositeEdge);

    // Reached the end of the mesh
    if (thirdNextFaceId == constants::missing::uintValue)
    {
        return freeHangingNodesAction;
    }

    const UInt thirdNextOppositeEdge = mesh.FindOppositeEdge(thirdNextFaceId, secondNextOppositeEdge);
    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(node6, mesh.GetEdge(thirdNextOppositeEdge).second);
    freeHangingNodesAction->Add(std::move(addEdgeAction));
    std::tie(newEdgeId, addEdgeAction) = mesh.ConnectNodes(node6, mesh.GetEdge(thirdNextOppositeEdge).first);
    freeHangingNodesAction->Add(std::move(addEdgeAction));

    return freeHangingNodesAction;
}

std::unique_ptr<meshkernel::UndoAction> meshkernel::ConnectMeshes::FreeHangingNodes(Mesh2D& mesh,
                                                                                    const UInt numberOfHangingNodes,
                                                                                    const std::vector<UInt>& hangingNodesOnEdge,
                                                                                    const UInt faceId,
                                                                                    const Edge& boundaryEdge,
                                                                                    const Point& boundaryNode,
                                                                                    const UInt edgeId)
{
    if (numberOfHangingNodes == 0)
    {
        return nullptr;
    }

    std::unique_ptr<CompoundUndoAction> freeHangingNodesAction = CompoundUndoAction::Create();

    BoundedIntegerArray hangingNodes;
    GetOrderedDistanceFromPoint(mesh, hangingNodesOnEdge, numberOfHangingNodes, boundaryNode, hangingNodes);

    UInt oppositeEdgeId = mesh.FindOppositeEdge(faceId, edgeId);

    if (oppositeEdgeId == constants::missing::uintValue)
    {
        throw ConstraintError("Opposite edge not found for element {} and edge {}.", faceId, edgeId);
    }

    UInt startNode = mesh.GetEdge(oppositeEdgeId).first;
    UInt endNode = mesh.GetEdge(oppositeEdgeId).second;

    const auto [segmentsCross,
                intersectionPoint,
                crossProduct,
                normalisedPolylineSegmentDistance,
                normalisedEdgeDistance] = AreSegmentsCrossing(mesh.Node(boundaryEdge.first),
                                                              mesh.Node(startNode),
                                                              mesh.Node(boundaryEdge.second),
                                                              mesh.Node(endNode),
                                                              false,
                                                              mesh.m_projection);
    if (segmentsCross)
    {
        std::swap(startNode, endNode);
    }

    switch (numberOfHangingNodes)
    {
    case 1:
        freeHangingNodesAction->Add(FreeOneHangingNode(mesh, hangingNodes, startNode, endNode));
        break;
    case 2:
        freeHangingNodesAction->Add(FreeTwoHangingNodes(mesh, faceId, oppositeEdgeId, hangingNodes, startNode, endNode));
        break;
    case 3:
        freeHangingNodesAction->Add(FreeThreeHangingNodes(mesh, faceId, oppositeEdgeId, hangingNodes, startNode, endNode));
        break;
    case 4:
        freeHangingNodesAction->Add(FreeFourHangingNodes(mesh, faceId, oppositeEdgeId, hangingNodes, startNode, endNode));
        break;
    default:
        // 0 hanging nodes is handled at the top of this function, so to be here can only be: numberOfHangingNodes > 4
        throw NotImplementedError("Cannot handle more than 4 hanging nodes along irregular edge, number is: {}",
                                  numberOfHangingNodes);
    }

    return freeHangingNodesAction;
}
