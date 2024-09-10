#include "MeshKernel/SplitRowColumnOfMesh.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Operations.hpp"

#include <tuple>

std::unique_ptr<meshkernel::UndoAction> meshkernel::SplitRowColumnOfMesh::Compute(Mesh2D& mesh, const UInt edgeId) const
{
    if (!IsValid(edgeId))
    {
        throw ConstraintError("Invalid edge-id");
    }

    if (edgeId >= mesh.GetNumEdges())
    {
        throw ConstraintError("edge-id is greater than the number of edges: {} >= {}", edgeId, mesh.GetNumEdges());
    }

    if (!IsValidEdge(mesh, edgeId))
    {
        throw ConstraintError("Invalid edge-id: {}.", edgeId);
    }

    // Just an estimate on the maximum number number of elements that may be split.
    const UInt initialNumberOfFaces = static_cast<UInt>(std::sqrt(static_cast<double>(mesh.GetNumFaces())));

    std::vector<UInt> elementIds;
    std::vector<UInt> edgeIds;
    std::vector<UInt> edgesToDelete;
    elementIds.reserve(initialNumberOfFaces);
    edgeIds.reserve(initialNumberOfFaces);
    edgesToDelete.reserve(initialNumberOfFaces);

    std::unique_ptr<CompoundUndoAction> undoActions;

    if (CanBeSplit(mesh, edgeId))
    {
        undoActions = CompoundUndoAction::Create();

        CollectElementsToSplit(mesh, edgeId, elementIds, edgeIds);
        SplitAlongRow(mesh, elementIds, edgeIds, *undoActions, edgesToDelete);

        std::ranges::for_each(edgesToDelete, [&](UInt edgeId)
                              { undoActions->Add(mesh.DeleteEdge(edgeId)); });

        mesh.Administrate();
    }

    return undoActions;
}

bool meshkernel::SplitRowColumnOfMesh::IsValidEdge(const Mesh2D& mesh, const UInt edgeId) const
{
    const Edge& edge = mesh.GetEdge(edgeId);

    return meshkernel::IsValidEdge(edge) &&
           mesh.Node(edge.first).IsValid() &&
           mesh.Node(edge.second).IsValid();
}

bool meshkernel::SplitRowColumnOfMesh::CanBeSplit(const Mesh2D& mesh, const UInt edgeId) const
{
    return IsValid(edgeId) && (IsQuadrilateral(mesh, mesh.m_edgesFaces[edgeId][0]) || IsQuadrilateral(mesh, mesh.m_edgesFaces[edgeId][1]));
}

void meshkernel::SplitRowColumnOfMesh::SplitAlongRow(Mesh2D& mesh,
                                                     const std::vector<UInt>& elementIds,
                                                     const std::vector<UInt>& edgeIds,
                                                     CompoundUndoAction& undoActions,
                                                     std::vector<UInt>& edgesToDelete) const
{

    const bool loopDetected = edgeIds.front() == edgeIds.back();
    UInt newNode = constants::missing::uintValue;
    UInt firstNode = constants::missing::uintValue;

    for (UInt i = 0; i < elementIds.size(); ++i)
    {
        if (loopDetected && i == 0)
        {
            SplitFirstLoopElement(mesh, elementIds[i], edgeIds[i], firstNode, newNode, undoActions, edgesToDelete);
        }
        else if (loopDetected && i == elementIds.size() - 1)
        {
            // The last action in the loop is just connecting the two nodes from the neighbouring edges that close the loop of elements
            std::tie(std::ignore, undoActions.Insert()) = mesh.ConnectNodes(newNode, firstNode);
        }
        else
        {
            SplitElement(mesh, elementIds[i], edgeIds[i], newNode, undoActions, edgesToDelete);
        }
    }
}

meshkernel::UInt meshkernel::SplitRowColumnOfMesh::GetNextElement(const Mesh2D& mesh, const UInt elementId, const UInt edgeId) const
{

    if (!IsValid(elementId) || !IsValid(edgeId))
    {
        return constants::missing::uintValue;
    }

    if (!IsQuadrilateral(mesh, elementId))
    {
        return constants::missing::uintValue;
    }

    UInt oppositeEdgeId = OppositeEdgeId(mesh, elementId, edgeId);
    UInt nextElementId = Neighbour(mesh.m_edgesFaces[oppositeEdgeId], elementId);

    return nextElementId;
}

meshkernel::UInt meshkernel::SplitRowColumnOfMesh::OppositeEdgeId(const Mesh2D& mesh, const UInt elementId, const UInt edgeId) const
{
    if (!IsValid(elementId) || !IsValid(edgeId))
    {
        return constants::missing::uintValue;
    }

    UInt edgeIndex = mesh.GetEdgeIndex(elementId, edgeId);

    if (!IsValid(edgeIndex))
    {
        return constants::missing::uintValue;
    }

    UInt oppositeEdgeIndex = (edgeIndex + 2) % constants::geometric::numNodesInQuadrilateral;
    UInt oppositeEdgeId = mesh.m_facesEdges[elementId][oppositeEdgeIndex];

    return oppositeEdgeId;
}

meshkernel::UInt meshkernel::SplitRowColumnOfMesh::SplitEdge(Mesh2D& mesh, const UInt edgeId, std::vector<UInt>& edgesToDelete, CompoundUndoAction& undoActions) const
{
    const Edge edgeNode = mesh.GetEdge(edgeId);

    Point point = 0.5 * (mesh.Node(edgeNode.first) + mesh.Node(edgeNode.second));

    auto [newNodeId, undo] = mesh.InsertNode(point);
    undoActions.Add(std::move(undo));

    std::tie(std::ignore, undoActions.Insert()) = mesh.ConnectNodes(edgeNode.first, newNodeId);
    std::tie(std::ignore, undoActions.Insert()) = mesh.ConnectNodes(newNodeId, edgeNode.second);

    edgesToDelete.push_back(edgeId);
    return newNodeId;
}

void meshkernel::SplitRowColumnOfMesh::SplitFirstLoopElement(Mesh2D& mesh,
                                                             const UInt elementId,
                                                             const UInt edgeId,
                                                             UInt& firstNode,
                                                             UInt& secondNode,
                                                             CompoundUndoAction& undoActions,
                                                             std::vector<UInt>& edgesToDelete) const
{
    UInt oppositeEdgeId = OppositeEdgeId(mesh, elementId, edgeId);
    firstNode = SplitEdge(mesh, edgeId, edgesToDelete, undoActions);
    secondNode = SplitEdge(mesh, oppositeEdgeId, edgesToDelete, undoActions);

    std::tie(std::ignore, undoActions.Insert()) = mesh.ConnectNodes(firstNode, secondNode);
}

void meshkernel::SplitRowColumnOfMesh::SplitElement(Mesh2D& mesh,
                                                    const UInt elementId,
                                                    const UInt edgeId,
                                                    UInt& newNode,
                                                    CompoundUndoAction& undoActions,
                                                    std::vector<UInt>& edgesToDelete) const
{
    const Edge edgeNode = mesh.GetEdge(edgeId);
    const std::array<UInt, 2>& edgeFace = mesh.m_edgesFaces[edgeId];
    const UInt previousElementId = edgeFace[0] + edgeFace[1] - elementId;
    const UInt nextElementId = GetNextElement(mesh, elementId, edgeId);

    if (!IsValid(previousElementId))
    {
        if (!IsValid(nextElementId) || IsQuadrilateral(mesh, nextElementId))
        {
            UInt oppositeEdgeId = OppositeEdgeId(mesh, elementId, edgeId);
            UInt firstNewNodeId = SplitEdge(mesh, edgeId, edgesToDelete, undoActions);
            UInt secondNewNodeId = SplitEdge(mesh, oppositeEdgeId, edgesToDelete, undoActions);

            std::tie(std::ignore, undoActions.Insert()) = mesh.ConnectNodes(firstNewNodeId, secondNewNodeId);

            // Overwrite newNode with the latest new node id.
            newNode = secondNewNodeId;
        }
        else if (!IsQuadrilateral(mesh, nextElementId))
        {
            UInt oppositeEdgeId = OppositeEdgeId(mesh, elementId, edgeId);
            const Edge oppositeEdgeNode = mesh.GetEdge(oppositeEdgeId);
            UInt firstNewNodeId = SplitEdge(mesh, edgeId, edgesToDelete, undoActions);

            std::tie(std::ignore, undoActions.Insert()) = mesh.ConnectNodes(oppositeEdgeNode.first, firstNewNodeId);
            std::tie(std::ignore, undoActions.Insert()) = mesh.ConnectNodes(firstNewNodeId, oppositeEdgeNode.second);

            // Set the newNode with the null id, indicating no new node was created
            newNode = constants::missing::uintValue;
        }
    }
    else if (IsQuadrilateral(mesh, previousElementId))
    {
        if (!IsValid(nextElementId) || IsQuadrilateral(mesh, nextElementId))
        {
            // This is probably the most common case in a sequence of quadrilateral elements
            UInt oppositeEdgeId = OppositeEdgeId(mesh, elementId, edgeId);
            UInt secondNewNodeId = SplitEdge(mesh, oppositeEdgeId, edgesToDelete, undoActions);

            std::tie(std::ignore, undoActions.Insert()) = mesh.ConnectNodes(newNode, secondNewNodeId);

            // Overwrite newNode with the latest new node id.
            newNode = secondNewNodeId;
        }
        else if (!IsQuadrilateral(mesh, nextElementId))
        {
            UInt oppositeEdgeId = OppositeEdgeId(mesh, elementId, edgeId);
            const Edge oppositeEdgeNode = mesh.GetEdge(oppositeEdgeId);

            std::tie(std::ignore, undoActions.Insert()) = mesh.ConnectNodes(oppositeEdgeNode.first, newNode);
            std::tie(std::ignore, undoActions.Insert()) = mesh.ConnectNodes(newNode, oppositeEdgeNode.second);

            // Set the newNode with the null id, indicating no new node was created
            newNode = constants::missing::uintValue;
        }
    }
    else if (!IsQuadrilateral(mesh, previousElementId))
    {
        if (!IsValid(nextElementId) || IsQuadrilateral(mesh, nextElementId))
        {
            UInt oppositeEdgeId = OppositeEdgeId(mesh, elementId, edgeId);
            UInt secondNewNodeId = SplitEdge(mesh, oppositeEdgeId, edgesToDelete, undoActions);

            std::tie(std::ignore, undoActions.Insert()) = mesh.ConnectNodes(edgeNode.first, secondNewNodeId);
            std::tie(std::ignore, undoActions.Insert()) = mesh.ConnectNodes(secondNewNodeId, edgeNode.second);

            // Overwrite newNode with the latest new node id.
            newNode = secondNewNodeId;
        }

        // If neither the previous nor the next elements are quadrilaterals then there is nothing to do.
    }
}

bool meshkernel::SplitRowColumnOfMesh::IsQuadrilateral(const Mesh2D& mesh, const UInt elementId) const
{
    return IsValid(elementId) && mesh.m_numFacesNodes[elementId] == constants::geometric::numNodesInQuadrilateral;
}

void meshkernel::SplitRowColumnOfMesh::CollectElementsOneSideOfEdge(const Mesh2D& mesh,
                                                                    const UInt edgeId,
                                                                    const UInt whichSide,
                                                                    std::vector<UInt>& partialElementIds,
                                                                    std::vector<UInt>& partialEdgeIds,
                                                                    bool& loopDetected) const
{
    const UInt firstElementId = mesh.m_edgesFaces[edgeId][0] != constants::missing::uintValue ? mesh.m_edgesFaces[edgeId][0] : mesh.m_edgesFaces[edgeId][1];

    partialElementIds.clear();
    partialEdgeIds.clear();

    UInt elementId = mesh.m_edgesFaces[edgeId][whichSide];
    UInt currentEdgeId = edgeId;
    bool firstIteration = true;
    loopDetected = false;

    while (IsQuadrilateral(mesh, elementId) && !loopDetected)
    {
        partialElementIds.push_back(elementId);

        // Ensure the seed edge is added only once to the sequence of edges
        if (whichSide == 0 || !firstIteration)
        {
            partialEdgeIds.push_back(currentEdgeId);
        }

        GetNextEdge(mesh, elementId, currentEdgeId);

        if (elementId == firstElementId)
        {
            loopDetected = true;
        }

        firstIteration = false;
    }

    if (IsValid(currentEdgeId))
    {
        partialEdgeIds.push_back(currentEdgeId);
    }
}

void meshkernel::SplitRowColumnOfMesh::CollectElementsToSplit(const Mesh2D& mesh, const UInt edgeId, std::vector<UInt>& elementIds, std::vector<UInt>& edgeIds) const
{
    // Estimate on the number of faces that will be refined.
    const UInt initialNumberOfFaces = static_cast<UInt>(std::sqrt(static_cast<double>(mesh.GetNumFaces())));

    std::vector<UInt> partialElementIds;
    std::vector<UInt> partialEdgeIds;

    partialElementIds.reserve(initialNumberOfFaces);
    partialEdgeIds.reserve(initialNumberOfFaces);

    // Loop over each side of the edge
    for (UInt i = 0; i < 2; ++i)
    {
        if (!IsValid(mesh.m_edgesFaces[edgeId][i]))
        {
            // There is no element on the other side of the edge.
            continue;
        }

        // Indicate if a circular loop of elements has been detected
        bool loopDetected = false;
        CollectElementsOneSideOfEdge(mesh, edgeId, i, partialElementIds, partialEdgeIds, loopDetected);

        if (i == 0)
        {
            elementIds = partialElementIds;
            edgeIds = partialEdgeIds;
        }
        else
        {
            elementIds.insert(elementIds.begin(), partialElementIds.rbegin(), partialElementIds.rend());
            edgeIds.insert(edgeIds.begin(), partialEdgeIds.rbegin(), partialEdgeIds.rend());
        }

        if (loopDetected)
        {
            // If a loop of elements has been detected then do not continue with the neighbouring element of the original edge.
            // It will find the same sequence of elements in reverse.
            break;
        }
    }
}

void meshkernel::SplitRowColumnOfMesh::GetNextEdge(const Mesh2D& mesh, UInt& elementId, UInt& edgeId) const
{
    if (!IsValid(edgeId))
    {
        elementId = constants::missing::uintValue;
        return;
    }

    if (!IsValid(elementId) || !IsQuadrilateral(mesh, elementId))
    {
        edgeId = constants::missing::uintValue;
        elementId = constants::missing::uintValue;
        return;
    }

    edgeId = OppositeEdgeId(mesh, elementId, edgeId);

    if (IsValid(edgeId))
    {
        elementId = Neighbour(mesh.m_edgesFaces[edgeId], elementId);
    }
    else
    {
        elementId = constants::missing::uintValue;
    }
}
