#include "MeshKernel/SplitRowColumnOfMesh3.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Operations.hpp"

std::unique_ptr<meshkernel::UndoAction> meshkernel::SplitRowColumnOfMesh3::Compute(Mesh2D& mesh, const UInt edgeId) const
{
    if (edgeId == constants::missing::uintValue)
    {
        throw ConstraintError("Invalid edge-id");
    }

    if (edgeId >= mesh.GetNumValidEdges())
    {
        throw ConstraintError("edge-id is greater than the number of edges");
    }

    // TODO check elements have been computed.

    // const UInt initialNumberOfElements = static_cast<UInt>(std::sqrt(static_cast<double>(mesh.GetNumFaces())));

    std::vector<UInt> elementIds;
    std::vector<UInt> edgeIds;
    std::vector<UInt> edgesToDelete;
    // elementIds.reserve(initialNumberOfElements);
    // edgeIds.reserve(initialNumberOfElements);

    std::unique_ptr<CompoundUndoAction> undoActions = CompoundUndoAction::Create();

    std::cout << "edge is valid: " << std::boolalpha << IsValidEdge(mesh, edgeId) << "  " << std::boolalpha << MayBeSplit(mesh, edgeId) << std::endl;

    if (IsValidEdge(mesh, edgeId) && MayBeSplit(mesh, edgeId))
    {
        std::cout << "start edge: { "
                  << mesh.Node(mesh.GetEdge(edgeId).first).x << ", " << mesh.Node(mesh.GetEdge(edgeId).first).y << " } -- { "
                  << mesh.Node(mesh.GetEdge(edgeId).second).x << ", " << mesh.Node(mesh.GetEdge(edgeId).second).y << " } -- "
                  << std::endl;

        CollectElementIdsToSplit(mesh, edgeId, elementIds, edgeIds);

        for (UInt i = 0; i < edgeIds.size(); ++i)
        {

            std::cout << "edge: " << edgeIds[i] << "  ";
            if (edgeIds[i] != constants::missing::uintValue)
            {
                std::array<UInt, 2> faces = mesh.m_edgesFaces[edgeIds[i]];
                std::cout << mesh.Node(mesh.GetEdge(edgeIds[i]).first).x << ", " << mesh.Node(mesh.GetEdge(edgeIds[i]).first).y << " -- "
                          << mesh.Node(mesh.GetEdge(edgeIds[i]).second).x << ", " << mesh.Node(mesh.GetEdge(edgeIds[i]).second).y << " -- "
                          << faces[0] << " " << faces[1]
                          << std::endl;
            }
            else
            {
                std::cout << " -- null edge -- " << std::endl;
            }
        }

        for (UInt i = 0; i < elementIds.size(); ++i)
        {
            std::cout << "elementId : " << elementIds[i] << "  " << std::endl;
        }

        SplitAlongRow(mesh, elementIds, edgeIds, *undoActions, edgesToDelete);

        for (UInt i = 0; i < edgesToDelete.size(); ++i)
        {
            undoActions->Add(mesh.DeleteEdge(edgesToDelete[i]));
        }
    }

    mesh.Administrate();
    return undoActions;
}

bool meshkernel::SplitRowColumnOfMesh3::IsValidEdge(const Mesh2D& mesh, const UInt edgeId) const
{
    const Edge& edge = mesh.GetEdge(edgeId);

    return meshkernel::IsValidEdge(edge) &&
           mesh.Node(edge.first).IsValid() &&
           mesh.Node(edge.second).IsValid();
}

bool meshkernel::SplitRowColumnOfMesh3::MayBeSplit(const Mesh2D& mesh, const UInt edgeId) const
{
    if (edgeId == constants::missing::uintValue)
    {
        // Should this be an error condition
        return false;
    }

    const std::array<UInt, 2>& edge = mesh.m_edgesFaces[edgeId];

    if (edge[0] == constants::missing::uintValue && edge[1] == constants::missing::uintValue)
    {
        return false;
    }
    else if (edge[0] != constants::missing::uintValue && edge[1] == constants::missing::uintValue)
    {
        return mesh.m_numFacesNodes[edge[0]] == 4;
    }
    else if (edge[1] == constants::missing::uintValue && edge[1] != constants::missing::uintValue)
    {
        return mesh.m_numFacesNodes[edge[1]] == 4;
    }
    else
    {
        return mesh.m_numFacesNodes[edge[0]] == 4 || mesh.m_numFacesNodes[edge[1]] == 4;
    }
}

void meshkernel::SplitRowColumnOfMesh3::SplitAlongRow(Mesh2D& mesh, const std::vector<UInt>& elementIds, const std::vector<UInt>& edgeIds, CompoundUndoAction& undoActions,
                                                      std::vector<UInt>& edgesToDelete) const
{

    auto elementIter = elementIds.begin();
    UInt newNode = constants::missing::uintValue;

    for (UInt i = 0; i < edgeIds.size() - 1; ++i)
    {
        UInt elementId = elementIter == elementIds.end() ? constants::missing::uintValue : *elementIter;
        SplitEdge5(mesh, elementId, edgeIds[i], newNode, undoActions, edgesToDelete);
        ++elementIter;
    }

    // for (auto edgeIter = edgeIds.begin(); edgeIter != edgeIds.end(); ++edgeIter)
    // {
    //     UInt elementId = elementIter == elementIds.end() ? constants::missing::uintValue : *elementIter;
    //     SplitEdge4(mesh, elementId, *edgeIter, newNode, undoActions, edgesToDelete);
    //     ++elementIter;
    // }
}

meshkernel::UInt meshkernel::SplitRowColumnOfMesh3::GetNextElementId (const Mesh2D& mesh, const UInt elementId, const UInt edgeId) const
{

    if (elementId == constants::missing::uintValue || edgeId == constants::missing::uintValue)
    {
        return constants::missing::uintValue;
    }

    if (elementId != constants::missing::uintValue && mesh.m_numFacesNodes [elementId] != 4)
    {
        return constants::missing::uintValue;
    }

    UInt edgeIndex = mesh.GetEdgeIndex(elementId, edgeId);
    UInt oppositeEdgeIndex = (edgeIndex + 2) % 4;
    UInt oppositeEdgeId = mesh.m_facesEdges[elementId][oppositeEdgeIndex];
    const std::array<UInt, 2>& oppositeEdgeFace = mesh.m_edgesFaces[oppositeEdgeId];
    UInt nextElement = oppositeEdgeFace[0] + oppositeEdgeFace[1] - elementId;

    return nextElement;
}

meshkernel::UInt meshkernel::SplitRowColumnOfMesh3::OppositeEdgeId (const Mesh2D& mesh, const UInt elementId, const UInt edgeId) const
{
    if (elementId == constants::missing::uintValue || edgeId == constants::missing::uintValue)
    {
        return constants::missing::uintValue;
    }

    UInt edgeIndex = mesh.GetEdgeIndex(elementId, edgeId);
    UInt oppositeEdgeIndex = (edgeIndex + 2) % 4;
    UInt oppositeEdgeId = mesh.m_facesEdges[elementId][oppositeEdgeIndex];

    return oppositeEdgeId;
}

void meshkernel::SplitRowColumnOfMesh3::SplitEdge (Mesh2D& mesh, const UInt edgeId, UInt& newNode, std::vector<UInt>& edgesToDelete, CompoundUndoAction& undoActions) const
{
    const Edge& edgeNode = mesh.GetEdge(edgeId);

    Point point = 0.5 * (mesh.Node(edgeNode.first) + mesh.Node(edgeNode.second));
    std::cout << "new node coords: " << point.x << ", " << point.y << std::endl;
    auto [newNodeId, undo] = mesh.InsertNode(point);
    newNode = newNodeId;
    undoActions.Add(std::move(undo));

    auto [newEdgeId1, newEdgeUndo1] = mesh.ConnectNodes(edgeNode.first, newNodeId);
    undoActions.Add(std::move(newEdgeUndo1));

    auto [newEdgeId2, newEdgeUndo2] = mesh.ConnectNodes(newNodeId, edgeNode.second);
    undoActions.Add(std::move(newEdgeUndo2));

    edgesToDelete.push_back (edgeId);
}

void meshkernel::SplitRowColumnOfMesh3::SplitEdge5(Mesh2D& mesh, const UInt elementId, const UInt edgeId, UInt& newNode, CompoundUndoAction& undoActions [[maybe_unused]], std::vector<UInt>& edgesToDelete [[maybe_unused]] ) const
{
    std::cout << " SplitRowColumnOfMesh3::SplitEdge4 " << elementId << "  " << edgeId << "  " << std::endl;

    [[maybe_unused]] const Edge& edgeNode = mesh.GetEdge(edgeId);
    [[maybe_unused]] const std::array<UInt, 2>& edgeFace = mesh.m_edgesFaces[edgeId];
    [[maybe_unused]] UInt previousElementId = edgeFace[0] + edgeFace[1] - elementId;
    [[maybe_unused]] UInt nextElementId = GetNextElementId (mesh, elementId, edgeId);

    std::cout << "edge info: "  << elementId << "  " << edgeId << "  " << edgeFace[0] << "  " << edgeFace[1] << "  " << newNode << "  "  << previousElementId << "  " << nextElementId << std::endl;

    if (previousElementId == constants::missing::uintValue && nextElementId == constants::missing::uintValue)
    {
        std::cout << "case 1" << std::endl;

        UInt oppositeEdgeId = OppositeEdgeId (mesh, elementId, edgeId);

        UInt firstNewNodeId;
        UInt secondNewNodeId;

        SplitEdge (mesh, edgeId, firstNewNodeId, edgesToDelete, undoActions);
        SplitEdge (mesh, oppositeEdgeId, secondNewNodeId, edgesToDelete, undoActions);

        UInt newEdgeId;
        std::tie (newEdgeId, undoActions) = mesh.ConnectNodes(firstNewNodeId, secondNewNodeId);
        // auto [newEdgeId, newEdgeUndo] = mesh.ConnectNodes(firstNewNodeId, secondNewNodeId);
        // undoActions.Add(std::move(newEdgeUndo));

        return;
    }

    if (previousElementId == constants::missing::uintValue && mesh.m_numFacesNodes[nextElementId] == 4)
    {
        UInt oppositeEdgeId = OppositeEdgeId (mesh, elementId, edgeId);

        UInt firstNewNodeId;
        UInt secondNewNodeId;

        SplitEdge (mesh, edgeId, firstNewNodeId, edgesToDelete, undoActions);
        SplitEdge (mesh, oppositeEdgeId, secondNewNodeId, edgesToDelete, undoActions);

        auto [newEdgeId5, newEdgeUndo5] = mesh.ConnectNodes(firstNewNodeId, secondNewNodeId);
        undoActions.Add(std::move(newEdgeUndo5));

        // Overwrite newNode with the latest new node id.
        newNode = secondNewNodeId;

        std::cout << "case 2" << std::endl;
        // Simple case
        return;
    }

    if (previousElementId == constants::missing::uintValue && mesh.m_numFacesNodes[nextElementId] != 4)
    {
        std::cout << "case 3 CHECK the connections are correct" << std::endl;

        UInt oppositeEdgeId = OppositeEdgeId (mesh, elementId, edgeId);
        const Edge& oppositeEdgeNode = mesh.GetEdge(oppositeEdgeId);

        UInt firstNewNodeId;
        SplitEdge (mesh, edgeId, firstNewNodeId, edgesToDelete, undoActions);

        auto [newEdgeId3, newEdgeUndo3] = mesh.ConnectNodes(oppositeEdgeNode.first, firstNewNodeId);
        undoActions.Add(std::move(newEdgeUndo3));

        auto [newEdgeId4, newEdgeUndo4] = mesh.ConnectNodes(firstNewNodeId, oppositeEdgeNode.second);
        undoActions.Add(std::move(newEdgeUndo4));

        // Simple case
        return;
    }

    if (nextElementId == constants::missing::uintValue && mesh.m_numFacesNodes[previousElementId] == 4)
    {
        std::cout << "case 4" << std::endl;

        UInt oppositeEdgeId = OppositeEdgeId (mesh, elementId, edgeId);

        UInt secondNewNodeId;

        SplitEdge (mesh, oppositeEdgeId, secondNewNodeId, edgesToDelete, undoActions);

        auto [newEdgeId5, newEdgeUndo5] = mesh.ConnectNodes(newNode, secondNewNodeId);
        undoActions.Add(std::move(newEdgeUndo5));

        // Overwrite newNode with the latest new node id.
        newNode = secondNewNodeId;

        // Simple case
        return;
    }

    if (nextElementId == constants::missing::uintValue && mesh.m_numFacesNodes[previousElementId] != 4)
    {
        std::cout << "case 5" << std::endl;
        // Should be like case 3
        throw ConstraintError ("case 5");
        // Simple case
        return;
    }

    if (mesh.m_numFacesNodes[previousElementId] != 4 && mesh.m_numFacesNodes[nextElementId] != 4)
    {
        std::cout << "case 6" << std::endl;
        // Simple case, nothing to do
        return;
    }

    if (mesh.m_numFacesNodes[previousElementId] == 4 && mesh.m_numFacesNodes[nextElementId] != 4)
    {
        std::cout << "case 7" << std::endl;

        UInt oppositeEdgeId = OppositeEdgeId (mesh, elementId, edgeId);
        const Edge& oppositeEdgeNode = mesh.GetEdge(oppositeEdgeId);

        std::cout << "connecting nodes: " << mesh.Node (oppositeEdgeNode.first).x << ", " << mesh.Node (oppositeEdgeNode.first).y << " to "
                  << mesh.Node (newNode).x << ", " << mesh.Node (newNode).y << std::endl;

        std::cout << "connecting nodes: " << mesh.Node (oppositeEdgeNode.second).x << ", " << mesh.Node (oppositeEdgeNode.second).y << " to "
                  << mesh.Node (newNode).x << ", " << mesh.Node (newNode).y << std::endl;

        auto [newEdgeId1, newEdgeUndo1] = mesh.ConnectNodes(oppositeEdgeNode.first, newNode);
        undoActions.Add(std::move(newEdgeUndo1));

        auto [newEdgeId2, newEdgeUndo2] = mesh.ConnectNodes(newNode, oppositeEdgeNode.second);
        undoActions.Add(std::move(newEdgeUndo2));


        // Simple case
        return;
    }

    if (mesh.m_numFacesNodes[previousElementId] != 4 && mesh.m_numFacesNodes[nextElementId] == 4)
    {
        std::cout << "case 8" << std::endl;

        UInt oppositeEdgeId = OppositeEdgeId (mesh, elementId, edgeId);

        UInt secondNewNodeId;

        SplitEdge (mesh, oppositeEdgeId, secondNewNodeId, edgesToDelete, undoActions);

        auto [newEdgeId1, newEdgeUndo1] = mesh.ConnectNodes(edgeNode.first, secondNewNodeId);
        undoActions.Add(std::move(newEdgeUndo1));

        auto [newEdgeId2, newEdgeUndo2] = mesh.ConnectNodes(secondNewNodeId, edgeNode.second);
        undoActions.Add(std::move(newEdgeUndo2));

        // Overwrite newNode with the latest new node id.
        newNode = secondNewNodeId;

        // Simple case
        return;
    }

    if (mesh.m_numFacesNodes[previousElementId] == 4 && mesh.m_numFacesNodes[nextElementId] == 4)
    {
        std::cout << "case 9" << std::endl;

        UInt oppositeEdgeId = OppositeEdgeId (mesh, elementId, edgeId);
        UInt secondNewNodeId;

        SplitEdge (mesh, oppositeEdgeId, secondNewNodeId, edgesToDelete, undoActions);

        auto [newEdgeId5, newEdgeUndo5] = mesh.ConnectNodes(newNode, secondNewNodeId);
        undoActions.Add(std::move(newEdgeUndo5));

        // Overwrite newNode with the latest new node id.
        newNode = secondNewNodeId;

        // Simple case
        return;
    }

    std::cout << "no defined case " << std::endl;
    throw ConstraintError ("no defined case");
}

void meshkernel::SplitRowColumnOfMesh3::SplitEdge4(Mesh2D& mesh, const UInt elementId, const UInt edgeId, UInt& newNode, CompoundUndoAction& undoActions, std::vector<UInt>& edgesToDelete) const
{
    std::cout << " SplitRowColumnOfMesh3::SplitEdge4 " << elementId << "  " << edgeId << "  " << std::endl;

    const Edge& edgeNode = mesh.GetEdge(edgeId);
    const std::array<UInt, 2>& edgeFace = mesh.m_edgesFaces[edgeId];
    UInt previousElementId = edgeFace[0] + edgeFace[1] - elementId;

    std::cout << "edge info: "  << elementId << "  " << edgeId << "  " << edgeFace[0] << "  " << edgeFace[1]
              << "  " << newNode << std::endl;

    // previousNode
    // If newNode == null and previousElementId == null
    // then create new previousNode and connect
    // otherwise
    // copy newNode to previousNode

    if (previousElementId == constants::missing::uintValue && elementId != constants::missing::uintValue && mesh.m_numFacesNodes[elementId] == 4)
    {
        std::cout << "case 1" << std::endl;
        // SplitEdge (mesh, edgeNode, undoActions);

        UInt edgeIndex = mesh.GetEdgeIndex(elementId, edgeId);
        UInt oppositeEdgeIndex = (edgeIndex + 2) % 4;
        UInt oppositeEdgeId = mesh.m_facesEdges[elementId][oppositeEdgeIndex];
        const std::array<UInt, 2>& oppositeEdgeFace = mesh.m_edgesFaces[oppositeEdgeId];
        const Edge& oppositeEdgeNode = mesh.GetEdge(oppositeEdgeId);
        UInt nextElement = oppositeEdgeFace[0] + oppositeEdgeFace[1] - elementId;

        Point point = 0.5 * (mesh.Node(edgeNode.first) + mesh.Node(edgeNode.second));
        std::cout << "new node coords: " << point.x << ", " << point.y << std::endl;
        auto [firstNewNodeId, undo] = mesh.InsertNode(point);
        undoActions.Add(std::move(undo));

        auto [newEdgeId1, newEdgeUndo1] = mesh.ConnectNodes(edgeNode.first, firstNewNodeId);
        undoActions.Add(std::move(newEdgeUndo1));

        auto [newEdgeId2, newEdgeUndo2] = mesh.ConnectNodes(firstNewNodeId, edgeNode.second);
        undoActions.Add(std::move(newEdgeUndo2));

        edgesToDelete.push_back(edgeId);

        newNode = firstNewNodeId;

        if (nextElement == constants::missing::uintValue || (nextElement != constants::missing::uintValue && mesh.m_numFacesNodes[elementId] == 4))
        {
            Point point = 0.5 * (mesh.Node(oppositeEdgeNode.first) + mesh.Node(oppositeEdgeNode.second));
            std::cout << "new node coords: " << point.x << ", " << point.y << std::endl;
            auto [secondNewNodeId, undo] = mesh.InsertNode(point);
            undoActions.Add(std::move(undo));

            auto [newEdgeId1, newEdgeUndo1] = mesh.ConnectNodes(oppositeEdgeNode.first, secondNewNodeId);
            undoActions.Add(std::move(newEdgeUndo1));

            auto [newEdgeId2, newEdgeUndo2] = mesh.ConnectNodes(secondNewNodeId, oppositeEdgeNode.second);
            undoActions.Add(std::move(newEdgeUndo2));

            auto [newEdgeId3, newEdgeUndo3] = mesh.ConnectNodes(firstNewNodeId, secondNewNodeId);
            undoActions.Add(std::move(newEdgeUndo1));

            edgesToDelete.push_back(oppositeEdgeId);

            // Overwrite newNode with the latest new node id.
            newNode = secondNewNodeId;
        }
        else
        {
            auto [newEdgeId1, undoAction1] = mesh.ConnectNodes(oppositeEdgeNode.first, firstNewNodeId);
            undoActions.Add(std::move(undoAction1));

            auto [newEdgeId2, undoAction2] = mesh.ConnectNodes(firstNewNodeId, oppositeEdgeNode.second);
            undoActions.Add(std::move(undoAction2));
        }

        return;
    }

    if (previousElementId != constants::missing::uintValue && mesh.m_numFacesNodes[previousElementId] == 4 &&
        elementId != constants::missing::uintValue && mesh.m_numFacesNodes[elementId] == 4)
    {
        std::cout << "case 2 " << mesh.m_numFacesNodes[previousElementId] << std::endl;
        // SplitEdge (mesh, edgeNode, undoActions);

        // UInt edgeIndex = mesh.GetEdgeIndex(elementId, edgeId);
        // UInt oppositeEdgeIndex = (edgeIndex + 2) % 4;
        // UInt oppositeEdgeId = mesh.m_facesEdges[elementId][oppositeEdgeIndex];
        // const std::array<UInt, 2>& oppositeEdgeFace = mesh.m_edgesFaces[oppositeEdgeId];
        // const Edge& oppositeEdgeNode = mesh.GetEdge(oppositeEdgeId);
        // UInt nextElement = oppositeEdgeFace[0] + oppositeEdgeFace[1] - elementId;

        Point point = 0.5 * (mesh.Node(edgeNode.first) + mesh.Node(edgeNode.second));
        std::cout << "new node coords: " << point.x << ", " << point.y << std::endl;
        auto [firstNewNodeId, undo] = mesh.InsertNode(point);
        undoActions.Add(std::move(undo));

        auto [newEdgeId1, newEdgeUndo1] = mesh.ConnectNodes(edgeNode.first, firstNewNodeId);
        undoActions.Add(std::move(newEdgeUndo1));

        auto [newEdgeId2, newEdgeUndo2] = mesh.ConnectNodes(firstNewNodeId, edgeNode.second);
        undoActions.Add(std::move(newEdgeUndo2));

        edgesToDelete.push_back(edgeId);

        auto [newEdgeId3, newEdgeUndo3] = mesh.ConnectNodes(newNode, firstNewNodeId);
        undoActions.Add(std::move(newEdgeUndo1));

        newNode = firstNewNodeId;
        return;
    }

    if (previousElementId != constants::missing::uintValue && elementId == constants::missing::uintValue)
    {
        std::cout << "case 3" << std::endl;
        // SplitEdge (mesh, edgeNode, undoActions);

        Point point = 0.5 * (mesh.Node(edgeNode.first) + mesh.Node(edgeNode.second));
        std::cout << "new node coords: " << point.x << ", " << point.y << std::endl;
        auto [firstNewNodeId, undo] = mesh.InsertNode(point);
        undoActions.Add(std::move(undo));

        auto [newEdgeId1, newEdgeUndo1] = mesh.ConnectNodes(edgeNode.first, firstNewNodeId);
        undoActions.Add(std::move(newEdgeUndo1));

        auto [newEdgeId2, newEdgeUndo2] = mesh.ConnectNodes(firstNewNodeId, edgeNode.second);
        undoActions.Add(std::move(newEdgeUndo2));

        edgesToDelete.push_back(edgeId);

        auto [newEdgeId3, newEdgeUndo3] = mesh.ConnectNodes(newNode, firstNewNodeId);
        undoActions.Add(std::move(newEdgeUndo1));

        newNode = firstNewNodeId;

        return;
    }

    if (elementId != constants::missing::uintValue && mesh.m_numFacesNodes[elementId] == 4)
    {
        std::cout << "case 4" << std::endl;
        // SplitEdge (mesh, edgeNode, undoActions);

        UInt edgeIndex = mesh.GetEdgeIndex(elementId, edgeId);
        UInt oppositeEdgeIndex = (edgeIndex + 2) % 4;
        UInt oppositeEdgeId = mesh.m_facesEdges[elementId][oppositeEdgeIndex];
        const std::array<UInt, 2>& oppositeEdgeFace = mesh.m_edgesFaces[oppositeEdgeId];
        const Edge& oppositeEdgeNode = mesh.GetEdge(oppositeEdgeId);
        UInt nextElement = oppositeEdgeFace[0] + oppositeEdgeFace[1] - elementId;

        if (previousElementId != constants::missing::uintValue && mesh.m_numFacesNodes[previousElementId] == 3)
        {
            Point point = 0.5 * (mesh.Node(oppositeEdgeNode.first) + mesh.Node(oppositeEdgeNode.second));
            std::cout << "new node coords: " << point.x << ", " << point.y << std::endl;
        }
        else if (nextElement != constants::missing::uintValue && mesh.m_numFacesNodes[elementId] == 4)
        {
        }
        else if (nextElement == constants::missing::uintValue)
        {
            Point point = 0.5 * (mesh.Node(edgeNode.first) + mesh.Node(edgeNode.second));
            std::cout << "new node coords: " << point.x << ", " << point.y << std::endl;
            auto [firstNewNodeId, undo] = mesh.InsertNode(point);
            undoActions.Add(std::move(undo));

            auto [newEdgeId1, newEdgeUndo1] = mesh.ConnectNodes(edgeNode.first, firstNewNodeId);
            undoActions.Add(std::move(newEdgeUndo1));

            auto [newEdgeId2, newEdgeUndo2] = mesh.ConnectNodes(firstNewNodeId, edgeNode.second);
            undoActions.Add(std::move(newEdgeUndo2));

            edgesToDelete.push_back(edgeId);

            auto [newEdgeId3, newEdgeUndo3] = mesh.ConnectNodes(newNode, firstNewNodeId);
            undoActions.Add(std::move(newEdgeUndo1));

            // should constants::missing::uintValue be assigned here
            newNode = firstNewNodeId;
        }

        return;
    }

    if (elementId != constants::missing::uintValue && mesh.m_numFacesNodes[elementId] == 3)
    {
        std::cout << "case 5" << std::endl;

        auto [newEdgeId1, newEdgeUndo1] = mesh.ConnectNodes(edgeNode.first, newNode);
        undoActions.Add(std::move(newEdgeUndo1));

        auto [newEdgeId2, newEdgeUndo2] = mesh.ConnectNodes(newNode, edgeNode.second);
        undoActions.Add(std::move(newEdgeUndo2));

        // Must terminate spltting
        newNode = constants::missing::uintValue;

        return;
    }
}

void meshkernel::SplitRowColumnOfMesh3::SplitAlongRow(Mesh2D& mesh, const UInt startEdgeId, CompoundUndoAction& undoActions) const
{

    // TODO need to check for circular connected element loops
    // Also should not do the reverse loop for circular connected element loops

    UInt previousNewNode = constants::missing::uintValue;
    UInt firstPreviousNewNode = constants::missing::uintValue;

    const UInt initialNumberOfElements = static_cast<UInt>(std::sqrt(static_cast<double>(mesh.GetNumFaces())));
    std::vector<UInt> edgesToDelete;
    edgesToDelete.reserve(initialNumberOfElements);
    [[maybe_unused]] bool firstIteration = true;

    // Loop over each side of the edge
    for (UInt i = 0; i < 2; ++i)
    {
        if (mesh.m_edgesFaces[startEdgeId][i] == constants::missing::uintValue)
        {
            continue;
        }

        [[maybe_unused]] UInt firstElementId = mesh.m_edgesFaces[startEdgeId][i];
        UInt elementId = mesh.m_edgesFaces[startEdgeId][i];
        UInt nextElementId = mesh.m_edgesFaces[startEdgeId][i]; // constants::missing::uintValue
        UInt edgeId = startEdgeId;

        UInt iterationCount = 0;

        // Exit also it current is same as first (except when it is the first)
        while ((elementId != constants::missing::uintValue && mesh.m_numFacesNodes[elementId] == 4)) // || (!firstIteration && previousNewNode != constants::missing::uintValue))
        {
            GetNextElement(mesh, edgeId, nextElementId);

            if (nextElementId == firstElementId)
            {
                std::cout << "Element loop found" << std::endl;
                edgeId = constants::missing::uintValue;
                elementId = constants::missing::uintValue;
            }
            else
            {

                // if (elementId == constants::missing::uintValue && edgeId != constants::missing::uintValue)
                {
                    // TODO could use the fact that the edge has been deleted to know when we started in the middle of a row/column
                    SplitEdge3(mesh, elementId, edgeId, previousNewNode, undoActions, edgesToDelete);
                }
            }

            if (firstIteration)
            {
                firstPreviousNewNode = previousNewNode;
            }

            elementId = nextElementId;
            firstIteration = false;
            iterationCount++;
        }

        previousNewNode = firstPreviousNewNode;
        std::cout << "================================" << std::endl;
    }

    for (UInt i = 0; i < edgesToDelete.size(); ++i)
    {
        undoActions.Add(mesh.DeleteEdge(edgesToDelete[i]));
    }
}

void meshkernel::SplitRowColumnOfMesh3::CollectElementIdsToSplit(const Mesh2D& mesh, const UInt edgeId, std::vector<UInt>& elementIds, std::vector<UInt>& edgeIds) const
{

    // TODO need to check for circular connected element loops
    // Also should not do the reverse loop for circular connected element loops

    // Loop over each side of the edge
    for (UInt i = 0; i < 2; ++i)
    {
        if (mesh.m_edgesFaces[edgeId][i] == constants::missing::uintValue)
        {
            continue;
        }

        std::cout << " mesh.m_numFacesNodes[mesh.m_edgesFaces[edgeId][i]] " << mesh.m_numFacesNodes[mesh.m_edgesFaces[edgeId][i]] << std::endl;

        // if (mesh.m_numFacesNodes[mesh.m_edgesFaces[edgeId][i]] != 4)
        // {
        //     continue;
        // }

        [[maybe_unused]] UInt firstElementId = mesh.m_edgesFaces[edgeId][i];
        UInt elementId = mesh.m_edgesFaces[edgeId][i];
        UInt currentEdgeId = edgeId;
        [[maybe_unused]] bool firstIteration = true;

        std::vector<UInt> partialElementIds;
        std::vector<UInt> partialEdgeIds;

        // Exit also it current is same as first (except when it is the first)
        while (elementId != constants::missing::uintValue && mesh.m_numFacesNodes[elementId] == 4)
        {

            // if (elementId != constants::missing::uintValue)
            {
                partialElementIds.push_back(elementId);
            }

            if (i == 0 || !firstIteration)
            {
                partialEdgeIds.push_back(currentEdgeId);
            }

            std::cout << "element: " << elementId << "  " << currentEdgeId << "   ";
            GetNextElement(mesh, currentEdgeId, elementId);
            std::cout << "  next element: " << elementId << "  " << currentEdgeId << "   " << std::endl;
            firstIteration = false;
        }

        if (currentEdgeId != constants::missing::uintValue)
        {
            partialEdgeIds.push_back(currentEdgeId);
        }

        // partialElementIds.push_back(elementId);

        if (elementId != constants::missing::uintValue && mesh.m_numFacesNodes[elementId] == 4)
        {
            partialElementIds.push_back(elementId);
            std::cout << "terminator element: " << elementId << "  " << currentEdgeId << std::endl;
        }

        if (elementId != constants::missing::uintValue)
        {
            std::cout << "terminator element: " << elementId << "  " << currentEdgeId << "  " << mesh.m_numFacesNodes[elementId] << std::endl;
        }

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

    }

    std::cout << "list size: " << elementIds.size() << "  " << edgeIds.size() << std::endl;
}

void meshkernel::SplitRowColumnOfMesh3::GetNextElement(const Mesh2D& mesh, UInt& edgeId, UInt& elementId) const
{
    if (edgeId == constants::missing::uintValue)
    {
        return;
    }

    if (elementId == constants::missing::uintValue || mesh.m_numFacesNodes[elementId] != 4)
    {
        return;
    }

    UInt edgeIndex = mesh.GetEdgeIndex(elementId, edgeId);

    if (edgeIndex == constants::missing::uintValue)
    {
        elementId = constants::missing::uintValue;
        return;
    }

    UInt oppositeEdgeIndex = (edgeIndex + 2) % 4;
    edgeId = mesh.m_facesEdges[elementId][oppositeEdgeIndex];

    const std::array<UInt, 2>& edgesFaces = mesh.m_edgesFaces[edgeId];

    elementId = edgesFaces[0] == elementId ? edgesFaces[1] : edgesFaces[0];
}

void meshkernel::SplitRowColumnOfMesh3::SplitEdge(Mesh2D& mesh, [[maybe_unused]] UInt elementId, UInt edgeId, UInt& previousNewNode, CompoundUndoAction& undoActions) const
{
    const Edge& edge = mesh.GetEdge(edgeId);
    const std::array<UInt, 2>& edgeFaces = mesh.m_edgesFaces[edgeId];
    const UInt startNode = edge.first;
    const UInt endNode = edge.second;
    const UInt nextElement = edgeFaces[0] + edgeFaces[1] - elementId;

    // if (startNode == constants::missing::uintValue || endNode == constants::missing::uintValue || !mesh.Node(startNode).IsValid() || !mesh.Node(endNode).IsValid())
    // {
    //     previousNewNode = constants::missing::uintValue;
    //     return;
    // }

    // if (mesh.NumberOfNodes(elementId) == 3)
    // {
    // }

    // [[maybe_unused]] UInt elementId = constants::missing::uintValue;

    Point point = 0.5 * (mesh.Node(startNode) + mesh.Node(endNode));
    std::cout << "new node coords: " << point.x << ", " << point.y << std::endl;

    auto [newNodeId, undo] = mesh.InsertNode(point);
    undoActions.Add(std::move(undo));

    // undoActions.Add(mesh.DeleteEdge(edgeId));
    // auto [newEdgeId1, newEdgeUndo1] = mesh.ConnectNodes(startNode, newNode);
    // undoActions.Add(std::move(newEdgeUndo1));

    auto [newEdgeId2, newEdgeUndo2] = mesh.ConnectNodes(newNodeId, endNode);
    undoActions.Add(std::move(newEdgeUndo2));

    if (previousNewNode != constants::missing::uintValue)
    {
        auto [newEdgeId, undoAction] = mesh.ConnectNodes(previousNewNode, newNodeId);
        undoActions.Add(std::move(undoAction));
    }

    if (nextElement != constants::missing::uintValue)
    {
        if (mesh.m_numFacesNodes[nextElement] == 3)
        {
            std::cout << "next is triangle" << std::endl;
            UInt edgeIndex = FindIndex(mesh.m_facesEdges[nextElement], edgeId);
            UInt leftEdgeId = mesh.m_facesEdges[nextElement][FindPreviousIndex(mesh.m_facesEdges[nextElement], edgeIndex)];
            UInt rightEdgeId = mesh.m_facesEdges[nextElement][FindNextIndex(mesh.m_facesEdges[nextElement], edgeIndex)];

            const Edge& leftEdge = mesh.GetEdge(leftEdgeId);
            const Edge& rightEdge = mesh.GetEdge(rightEdgeId);

            // TODO in NeighbourElement add checks for bounds and null id's
            [[maybe_unused]] UInt leftElementId = mesh.NeighbourElement(nextElement, leftEdgeId);
            [[maybe_unused]] UInt rightElementId = mesh.NeighbourElement(nextElement, rightEdgeId);

            std::cout << "edge info: " << edgeIndex << "  " << leftEdgeId << "  " << rightEdgeId << "  "
                      << edge.first << "  " << edge.second << "  "
                      << leftEdge.first << "  " << leftEdge.second << "  "
                      << rightEdge.first << "  " << rightEdge.second << "  "
                      << std::endl;

            std::cout << "positions "
                      << mesh.Node(edge.first).x << ", " << mesh.Node(edge.first).y << " -- "
                      << mesh.Node(edge.second).x << ", " << mesh.Node(edge.second).y << " -- "
                      << mesh.Node(rightEdge.second).x << ", " << mesh.Node(rightEdge.second).y << " -- "
                      << std::endl;

            UInt kl;
            [[maybe_unused]] UInt kr;
            UInt k1 = mesh.GetEdge(edgeId).first;
            UInt k2 = mesh.GetEdge(edgeId).second;

            if (mesh.GetEdge(leftEdgeId).first == k1 || mesh.GetEdge(leftEdgeId).second == k1)
            {
                kl = k1;
                kr = k2;
            }
            else
            {
                kl = k2;
                kr = k1;
            }

            std::cout << "long edge nodes: " << kl << "  " << kr << std::endl;

            UInt kll = mesh.GetEdge(leftEdgeId).first + mesh.GetEdge(leftEdgeId).second - kl;
            UInt krr = mesh.GetEdge(rightEdgeId).first + mesh.GetEdge(rightEdgeId).second - kr;

            if (leftElementId != constants::missing::uintValue && mesh.m_numFacesNodes[leftElementId] == 3)
            {
                std::cout << "left element is triangle" << std::endl;

                UInt klll = mesh.m_facesNodes[leftElementId][0] + mesh.m_facesNodes[leftElementId][1] + mesh.m_facesNodes[leftElementId][2] - kl - kll;

                Point n1 = mesh.Node(mesh.m_facesNodes[leftElementId][0]);
                Point n2 = mesh.Node(mesh.m_facesNodes[leftElementId][1]);
                Point n3 = mesh.Node(mesh.m_facesNodes[leftElementId][2]);

                std::cout << "left element: " << n1.x << ", " << n1.y << " -- "
                          << n2.x << ", " << n2.y << " -- "
                          << n3.x << ", " << n3.y << " -- "
                          << std::endl;

                n1 = mesh.Node(klll);
                n2 = mesh.Node(kl);
                n3 = mesh.Node(kll);

                std::cout << "test nodes: " << n1.x << ", " << n1.y << " -- "
                          << n2.x << ", " << n2.y << " -- "
                          << n3.x << ", " << n3.y << " -- "
                          << std::endl;

                [[maybe_unused]] UInt rightLeftEdgeIndex = FindIndex(mesh.m_facesEdges[leftElementId], leftEdgeId);
                [[maybe_unused]] UInt rightLeftEdgeNextIndex = FindPreviousIndex(mesh.m_facesEdges[leftElementId], rightLeftEdgeIndex);
                [[maybe_unused]] UInt rightLeftEdgeId = mesh.m_facesEdges[leftElementId][rightLeftEdgeNextIndex];

                // const Edge& right;

                std::cout << "node ids: " << newNodeId << "  " << kl << "  " << kll << "  " << klll << std::endl;

                double cross1 = NormalizedInnerProductTwoSegments(mesh.Node(klll), mesh.Node(kl), mesh.Node(kll), mesh.Node(newNodeId),
                                                                  mesh.m_projection);
                double cross2 = NormalizedInnerProductTwoSegments(mesh.Node(klll), mesh.Node(kll),
                                                                  mesh.Node(kl), mesh.Node(newNodeId),
                                                                  mesh.m_projection);
                double cross3 = NormalizedInnerProductTwoSegments(mesh.Node(klll), mesh.Node(kll),
                                                                  mesh.Node(kll), mesh.Node(newNodeId),
                                                                  mesh.m_projection);

                // double cross1 = NormalizedInnerProductTwoSegments(mesh.Node(edge.first), mesh.Node(edge.second),
                //                                                   mesh.Node(mesh.GetEdge(rightLeftEdgeId).first), mesh.Node(mesh.GetEdge(rightLeftEdgeId).second),
                //                                                   mesh.m_projection);
                // double cross2 = NormalizedInnerProductTwoSegments(mesh.Node(edge.first), mesh.Node(edge.second),
                //                                                   mesh.Node(mesh.GetEdge(leftEdgeId).first), mesh.Node(mesh.GetEdge(leftEdgeId).second),
                //                                                   mesh.m_projection);
                // double cross3 = NormalizedInnerProductTwoSegments(mesh.Node(edge.first), mesh.Node(edge.second),
                //                                                   mesh.Node(mesh.GetEdge(rightEdgeId).first), mesh.Node(mesh.GetEdge(leftEdgeId).second),
                //                                                   mesh.m_projection);

                std::cout << "cross products: " << cross1 << "  " << cross2 << "  " << cross3 << std::endl;

                UInt nodeToKeep = mesh.FindCommonNode(edgeId, leftEdgeId);
                UInt nodeToConnect1 = mesh.GetEdge(leftEdgeId).first + mesh.GetEdge(leftEdgeId).second - nodeToKeep;
                UInt nodeToConnect2 = mesh.GetEdge(rightLeftEdgeId).first + mesh.GetEdge(rightLeftEdgeId).second - nodeToKeep;

                std::cout << "nodes to connect: " << nodeToConnect1 << "  " << nodeToConnect2 << std::endl;

                undoActions.Add(mesh.DeleteEdge(leftEdgeId));

                auto [newEdgeId1, newEdgeUndo1] = mesh.ConnectNodes(nodeToConnect1, newNodeId);
                undoActions.Add(std::move(newEdgeUndo1));

                auto [newEdgeId2, newEdgeUndo2] = mesh.ConnectNodes(nodeToConnect2, newNodeId);
                undoActions.Add(std::move(newEdgeUndo1));
            }

            if (rightElementId != constants::missing::uintValue && mesh.m_numFacesNodes[rightElementId] == 3)
            {
                std::cout << "right element is triangle" << std::endl;

                UInt krrr = mesh.m_facesNodes[rightElementId][0] + mesh.m_facesNodes[rightElementId][1] + mesh.m_facesNodes[rightElementId][2] - kr - krr;

                Point n1 = mesh.Node(mesh.m_facesNodes[rightElementId][0]);
                Point n2 = mesh.Node(mesh.m_facesNodes[rightElementId][1]);
                Point n3 = mesh.Node(mesh.m_facesNodes[rightElementId][2]);

                std::cout << "right element: " << n1.x << ", " << n1.y << " -- "
                          << n2.x << ", " << n2.y << " -- "
                          << n3.x << ", " << n3.y << " -- "
                          << std::endl;

                n1 = mesh.Node(krrr);
                n2 = mesh.Node(kr);
                n3 = mesh.Node(krr);

                std::cout << "test nodes: " << n1.x << ", " << n1.y << " -- "
                          << n2.x << ", " << n2.y << " -- "
                          << n3.x << ", " << n3.y << " -- "
                          << std::endl;

                [[maybe_unused]] UInt rightRightEdgeIndex = FindIndex(mesh.m_facesEdges[rightElementId], rightEdgeId);
                [[maybe_unused]] UInt rightRightEdgeNextIndex = FindNextIndex(mesh.m_facesEdges[rightElementId], rightRightEdgeIndex);
                [[maybe_unused]] UInt rightRightEdgeId = mesh.m_facesEdges[rightElementId][rightRightEdgeNextIndex];

                // const Edge& right;

                double cross1 = NormalizedInnerProductTwoSegments(mesh.Node(krrr), mesh.Node(kr),
                                                                  mesh.Node(krr), mesh.Node(newNodeId),
                                                                  mesh.m_projection);
                double cross2 = NormalizedInnerProductTwoSegments(mesh.Node(krrr), mesh.Node(krr),
                                                                  mesh.Node(kr), mesh.Node(newNodeId),
                                                                  mesh.m_projection);
                double cross3 = NormalizedInnerProductTwoSegments(mesh.Node(krrr), mesh.Node(krr),
                                                                  mesh.Node(krr), mesh.Node(newNodeId),
                                                                  mesh.m_projection);

                // double cross1 = NormalizedInnerProductTwoSegments(mesh.Node(edge.first), mesh.Node(edge.second),
                //                                                   mesh.Node(mesh.GetEdge(rightRightEdgeId).first), mesh.Node(mesh.GetEdge(rightRightEdgeId).second),
                //                                                   mesh.m_projection);
                // double cross2 = NormalizedInnerProductTwoSegments(mesh.Node(edge.first), mesh.Node(edge.second),
                //                                                   mesh.Node(mesh.GetEdge(leftEdgeId).first), mesh.Node(mesh.GetEdge(leftEdgeId).second),
                //                                                   mesh.m_projection);
                // double cross3 = NormalizedInnerProductTwoSegments(mesh.Node(edge.first), mesh.Node(edge.second),
                //                                                   mesh.Node(mesh.GetEdge(rightEdgeId).first), mesh.Node(mesh.GetEdge(rightEdgeId).second),
                //                                                   mesh.m_projection);

                std::cout << "cross products: " << cross1 << "  " << cross2 << "  " << cross3 << std::endl;

                UInt nodeToKeep = mesh.FindCommonNode(edgeId, rightEdgeId);
                UInt nodeToConnect1 = mesh.GetEdge(rightEdgeId).first + mesh.GetEdge(rightEdgeId).second - nodeToKeep;
                UInt nodeToConnect2 = mesh.GetEdge(rightRightEdgeId).first + mesh.GetEdge(rightRightEdgeId).second - nodeToKeep;

                std::cout << "nodes to connect: " << nodeToConnect1 << "  " << nodeToConnect2 << std::endl;

                undoActions.Add(mesh.DeleteEdge(rightEdgeId));

                auto [newEdgeId1, newEdgeUndo1] = mesh.ConnectNodes(nodeToConnect1, newNodeId);
                undoActions.Add(std::move(newEdgeUndo1));

                auto [newEdgeId2, newEdgeUndo2] = mesh.ConnectNodes(nodeToConnect2, newNodeId);
                undoActions.Add(std::move(newEdgeUndo1));
            }
        }
    }

    undoActions.Add(mesh.DeleteEdge(edgeId));
    auto [newEdgeId1, newEdgeUndo1] = mesh.ConnectNodes(startNode, newNodeId);
    undoActions.Add(std::move(newEdgeUndo1));

    // TODO Set the previous new node to null id if: 1. no next element, 2. merged triangles
    previousNewNode = newNodeId;
}

void meshkernel::SplitRowColumnOfMesh3::SplitEdge2(Mesh2D& mesh, [[maybe_unused]] UInt elementId, UInt edgeId, UInt& previousNewNode, CompoundUndoAction& undoActions) const
{
    std::cout << "SplitRowColumnOfMesh3::SplitEdge2 " << edgeId << std::endl;

    const Edge& edge = mesh.GetEdge(edgeId);

    // if (!meshkernel::IsValidEdge (edge))
    // {
    //     previousNewNode = constants::missing::uintValue;
    //     return;
    // }

    const std::array<UInt, 2>& edgeFaces = mesh.m_edgesFaces[edgeId];
    const UInt startNode = edge.first;
    const UInt endNode = edge.second;
    const UInt nextElement = edgeFaces[0] + edgeFaces[1] - elementId;

    std::cout << "nextElement " << nextElement << std::endl;
    std::cout << "node ids: " << startNode << "  " << endNode << std::endl;

    Point point = 0.5 * (mesh.Node(startNode) + mesh.Node(endNode));
    std::cout << "new node coords: " << point.x << ", " << point.y << std::endl;

    if (nextElement != constants::missing::uintValue && mesh.m_numFacesNodes[nextElement] == 3)
    {
        auto [newEdgeId1, undoAction1] = mesh.ConnectNodes(startNode, previousNewNode);
        undoActions.Add(std::move(undoAction1));

        auto [newEdgeId2, undoAction2] = mesh.ConnectNodes(endNode, previousNewNode);
        undoActions.Add(std::move(undoAction2));

        previousNewNode = constants::missing::uintValue;
    }
    else
    {
        auto [newNodeId, undo] = mesh.InsertNode(point);
        undoActions.Add(std::move(undo));

        std::cout << "new node: " << newNodeId << std::endl;

        std::cout << "deleting edge: " << edgeId << "  "
                  << mesh.Node(edge.first).x << ", " << mesh.Node(edge.first).y << " -- "
                  << mesh.Node(edge.second).x << ", " << mesh.Node(edge.second).y
                  << std::endl;

        undoActions.Add(mesh.DeleteEdge(edgeId));
        auto [newEdgeId1, newEdgeUndo1] = mesh.ConnectNodes(startNode, newNodeId);
        undoActions.Add(std::move(newEdgeUndo1));

        auto [newEdgeId2, newEdgeUndo2] = mesh.ConnectNodes(newNodeId, endNode);
        undoActions.Add(std::move(newEdgeUndo2));

        if (previousNewNode != constants::missing::uintValue)
        {
            // undoActions.Add(mesh.DeleteEdge(edgeId));
            auto [newEdgeId1, newEdgeUndo1] = mesh.ConnectNodes(previousNewNode, newNodeId);
            undoActions.Add(std::move(newEdgeUndo1));
        }

        // TODO Set the previous new node to null id if: 1. no next element, 2. merged triangles
        previousNewNode = newNodeId;
    }

    // if (nextElement != constants::missing::uintValue)
    // {
    //     if (mesh.m_numFacesNodes[nextElement] == 3)
    //     {
    //         auto [newEdgeId1, undoAction1] = mesh.ConnectNodes(startNode, previousNewNode);
    //         undoActions.Add(std::move(undoAction1));

    //         auto [newEdgeId2, undoAction2] = mesh.ConnectNodes(endNode, previousNewNode);
    //         undoActions.Add(std::move(undoAction2));

    //         previousNewNode = constants::missing::uintValue;
    //     }
    //     else
    //     {
    //         auto [newNodeId, undo] = mesh.InsertNode(point);
    //         undoActions.Add(std::move(undo));

    //         undoActions.Add(mesh.DeleteEdge(edgeId));
    //         auto [newEdgeId1, newEdgeUndo1] = mesh.ConnectNodes(startNode, newNodeId);
    //         undoActions.Add(std::move(newEdgeUndo1));

    //         auto [newEdgeId2, newEdgeUndo2] = mesh.ConnectNodes(newNodeId, endNode);
    //         undoActions.Add(std::move(newEdgeUndo2));

    //         if (previousNewNode != constants::missing::uintValue)
    //         {
    //             undoActions.Add(mesh.DeleteEdge(edgeId));
    //             auto [newEdgeId1, newEdgeUndo1] = mesh.ConnectNodes(previousNewNode, newNodeId);
    //             undoActions.Add(std::move(newEdgeUndo1));
    //         }

    //         // TODO Set the previous new node to null id if: 1. no next element, 2. merged triangles
    //         previousNewNode = newNodeId;
    //     }
    // }
}

void meshkernel::SplitRowColumnOfMesh3::SplitEdge3(Mesh2D& mesh, UInt elementId, UInt edgeId, UInt& previousNewNode, CompoundUndoAction& undoActions, std::vector<UInt>& edgesToDelete) const
{
    std::cout << "SplitRowColumnOfMesh3::SplitEdge3 " << elementId << "  " << edgeId << std::endl;

    if (elementId == constants::missing::uintValue || mesh.m_numFacesNodes[elementId] != 4)
    {
        // Could this be handled outside of this function?
        previousNewNode = constants::missing::uintValue;
        return;
    }

    const Edge& edge = mesh.GetEdge(edgeId);

    UInt edgeIndex = mesh.GetEdgeIndex(elementId, edgeId);

    if (edgeIndex == constants::missing::uintValue)
    {
        // Or should this be an error case.
        std::cout << " case 1 " << std::endl;
        previousNewNode = constants::missing::uintValue;
        return;
    }

    UInt oppositeEdgeIndex = (edgeIndex + 2) % 4;
    // This is the next edge
    const UInt previousEdgeId = mesh.m_facesEdges[elementId][oppositeEdgeIndex];
    const std::array<UInt, 2>& previousEdgeFaces = mesh.m_edgesFaces[previousEdgeId];
    const Edge& previousEdge = mesh.GetEdge(previousEdgeId);

    UInt previousElementId = previousEdgeFaces[0] + previousEdgeFaces[1] - elementId;

    const std::array<UInt, 2>& edgeFaces = mesh.m_edgesFaces[edgeId];
    const UInt nextElementId = edgeFaces[0] + edgeFaces[1] - elementId;

    std::cout << "elements " << previousElementId << "  " << elementId << "  " << nextElementId << std::endl;

    std::cout << "number of nodes: "
              << (previousElementId == constants::missing::uintValue ? 0 : mesh.m_numFacesNodes[previousElementId]) << "  "
              << (elementId == constants::missing::uintValue ? 0 : mesh.m_numFacesNodes[elementId]) << "  "
              << (nextElementId == constants::missing::uintValue ? 0 : mesh.m_numFacesNodes[nextElementId]) << "  "
              << std::endl;

    // TODO
    // How best to get started?
    // handle case of elememtId == null,

    if (previousElementId == constants::missing::uintValue && nextElementId == constants::missing::uintValue)
    {
        // TODO
        // Refine previous edge and edge, joining the newly created nodes
        std::cout << " case 2" << std::endl;
        previousNewNode = constants::missing::uintValue;
        return;
    }

    if ((previousElementId == constants::missing::uintValue || mesh.m_numFacesNodes[previousElementId] != 4) && (nextElementId != constants::missing::uintValue && mesh.m_numFacesNodes[nextElementId] != 4))
    {
        Point point = 0.5 * (mesh.Node(edge.first) + mesh.Node(edge.second));
        std::cout << "new node coords: " << point.x << ", " << point.y << std::endl;
        auto [newNodeId, undo] = mesh.InsertNode(point);
        std::cout << "new node: " << newNodeId << std::endl;
        undoActions.Add(std::move(undo));

        edgesToDelete.push_back(edgeId);
        // undoActions.Add(mesh.DeleteEdge(edgeId));

        auto [newEdgeId1, newEdgeUndo1] = mesh.ConnectNodes(previousEdge.first, newNodeId);
        undoActions.Add(std::move(newEdgeUndo1));

        auto [newEdgeId2, newEdgeUndo2] = mesh.ConnectNodes(newNodeId, previousEdge.second);
        undoActions.Add(std::move(newEdgeUndo2));

        previousNewNode = newNodeId;
        std::cout << " case 3" << std::endl;
        return;
    }

    if ((previousElementId == constants::missing::uintValue) && (nextElementId != constants::missing::uintValue && mesh.m_numFacesNodes[nextElementId] == 4))
    {
        Point point = 0.5 * (mesh.Node(previousEdge.first) + mesh.Node(previousEdge.second));
        std::cout << "new node coords: " << point.x << ", " << point.y << std::endl;
        auto [firstNewNodeId, undo1] = mesh.InsertNode(point);
        std::cout << "first new node: " << firstNewNodeId << std::endl;
        undoActions.Add(std::move(undo1));

        std::cout << "this edge "
                  << mesh.Node(edge.first).x << ", " << mesh.Node(edge.first).y << " -- "
                  << mesh.Node(edge.second).x << ", " << mesh.Node(edge.second).y << " -- "
                  << std::endl;

        std::cout << "next edge "
                  << mesh.Node(previousEdge.first).x << ", " << mesh.Node(previousEdge.first).y << " -- "
                  << mesh.Node(previousEdge.second).x << ", " << mesh.Node(previousEdge.second).y << " -- "
                  << std::endl;

        point = 0.5 * (mesh.Node(edge.first) + mesh.Node(edge.second));
        std::cout << "new node coords: " << point.x << ", " << point.y << std::endl;
        auto [secondNewNodeId, undo2] = mesh.InsertNode(point);
        std::cout << "second new node: " << firstNewNodeId << std::endl;
        undoActions.Add(std::move(undo2));

        edgesToDelete.push_back(previousEdgeId);
        edgesToDelete.push_back(edgeId);
        // undoActions.Add(mesh.DeleteEdge(edgeId));

        auto [newEdgeId1, newEdgeUndo1] = mesh.ConnectNodes(previousEdge.first, firstNewNodeId);
        undoActions.Add(std::move(newEdgeUndo1));

        auto [newEdgeId2, newEdgeUndo2] = mesh.ConnectNodes(firstNewNodeId, previousEdge.second);
        undoActions.Add(std::move(newEdgeUndo2));

        auto [newEdgeId3, newEdgeUndo3] = mesh.ConnectNodes(edge.first, secondNewNodeId);
        undoActions.Add(std::move(newEdgeUndo3));

        auto [newEdgeId4, newEdgeUndo4] = mesh.ConnectNodes(secondNewNodeId, edge.second);
        undoActions.Add(std::move(newEdgeUndo4));

        auto [newEdgeId5, newEdgeUndo5] = mesh.ConnectNodes(firstNewNodeId, secondNewNodeId);
        undoActions.Add(std::move(newEdgeUndo5));

        previousNewNode = secondNewNodeId;
        std::cout << " case 3b" << std::endl;
        return;
    }

    if ((previousElementId == constants::missing::uintValue || mesh.m_numFacesNodes[previousElementId] != 4) && (nextElementId != constants::missing::uintValue && mesh.m_numFacesNodes[nextElementId] == 4))
    {
        Point point = 0.5 * (mesh.Node(edge.first) + mesh.Node(edge.second));
        std::cout << "new node coords: " << point.x << ", " << point.y << std::endl;
        auto [newNodeId, undo] = mesh.InsertNode(point);
        std::cout << "new node: " << newNodeId << std::endl;
        undoActions.Add(std::move(undo));

        edgesToDelete.push_back(edgeId);
        // undoActions.Add(mesh.DeleteEdge(edgeId));

        auto [newEdgeId3, newEdgeUndo3] = mesh.ConnectNodes(edge.first, newNodeId);
        undoActions.Add(std::move(newEdgeUndo3));

        auto [newEdgeId4, newEdgeUndo4] = mesh.ConnectNodes(newNodeId, edge.second);
        undoActions.Add(std::move(newEdgeUndo4));

        auto [newEdgeId1, newEdgeUndo1] = mesh.ConnectNodes(previousEdge.first, newNodeId);
        undoActions.Add(std::move(newEdgeUndo1));

        auto [newEdgeId2, newEdgeUndo2] = mesh.ConnectNodes(newNodeId, previousEdge.second);
        undoActions.Add(std::move(newEdgeUndo2));

        previousNewNode = newNodeId;
        std::cout << " case 3a" << std::endl;
        return;
    }

    if (previousElementId == constants::missing::uintValue && nextElementId != constants::missing::uintValue && mesh.m_numFacesNodes[nextElementId] == 4)
    {
        Point point = 0.5 * (mesh.Node(previousEdge.first) + mesh.Node(previousEdge.second));
        std::cout << "new node coords: " << point.x << ", " << point.y << std::endl;
        auto [newNodeId1, undo1] = mesh.InsertNode(point);

        auto [newEdgeId1, newEdgeUndo1] = mesh.ConnectNodes(previousEdge.first, newNodeId1);
        undoActions.Add(std::move(newEdgeUndo1));

        auto [newEdgeId2, newEdgeUndo2] = mesh.ConnectNodes(newNodeId1, previousEdge.second);
        undoActions.Add(std::move(newEdgeUndo2));

        point = 0.5 * (mesh.Node(edge.first) + mesh.Node(edge.second));
        std::cout << "new node coords: " << point.x << ", " << point.y << std::endl;
        auto [newNodeId2, undo2] = mesh.InsertNode(point);

        auto [newEdgeId3, newEdgeUndo3] = mesh.ConnectNodes(edge.first, newNodeId2);
        undoActions.Add(std::move(newEdgeUndo3));

        auto [newEdgeId4, newEdgeUndo4] = mesh.ConnectNodes(newNodeId2, edge.second);
        undoActions.Add(std::move(newEdgeUndo4));

        auto [newEdgeId5, newEdgeUndo5] = mesh.ConnectNodes(newNodeId1, newNodeId2);
        undoActions.Add(std::move(newEdgeUndo5));
        previousNewNode = newNodeId2;
        return;
    }

    if (mesh.m_numFacesNodes[previousElementId] == 4 && (nextElementId != constants::missing::uintValue && mesh.m_numFacesNodes[nextElementId] != 4))
    // if (mesh.m_numFacesNodes[previousElementId] == 4 && (nextElementId == constants::missing::uintValue || mesh.m_numFacesNodes[nextElementId] != 4))
    {

        if (previousNewNode == constants::missing::uintValue)
        {
            Point point = 0.5 * (mesh.Node(previousEdge.first) + mesh.Node(previousEdge.second));
            std::cout << "new node coords: " << point.x << ", " << point.y << std::endl;
            auto [newNodeId1, undo1] = mesh.InsertNode(point);

            auto [newEdgeId1, newEdgeUndo1] = mesh.ConnectNodes(edge.first, newNodeId1);
            undoActions.Add(std::move(newEdgeUndo1));

            auto [newEdgeId2, newEdgeUndo2] = mesh.ConnectNodes(newNodeId1, edge.second);
            undoActions.Add(std::move(newEdgeUndo2));

            previousNewNode = newNodeId1;
        }
        else
        {
            std::cout << "previous node: "
                      << mesh.Node(previousEdge.first).x << ", " << mesh.Node(previousEdge.first).y << " -- "
                      << mesh.Node(previousEdge.second).x << ", " << mesh.Node(previousEdge.second).y << " -- "
                      << mesh.Node(previousNewNode).x << ", " << mesh.Node(previousNewNode).y
                      << std::endl;
        }

        std::cout << " case 4 " << std::endl;

        auto [newEdgeId1, undoAction1] = mesh.ConnectNodes(edge.first, previousNewNode);
        undoActions.Add(std::move(undoAction1));

        auto [newEdgeId2, undoAction2] = mesh.ConnectNodes(previousNewNode, edge.second);
        undoActions.Add(std::move(undoAction2));

        return;
    }

    // if ((previousElementId != constants::missing::uintValue && mesh.m_numFacesNodes[previousElementId] == 4) && (nextElementId != constants::missing::uintValue && mesh.m_numFacesNodes[nextElementId] == 4))
    // // if (mesh.m_numFacesNodes[previousElementId] == 4 && (nextElementId == constants::missing::uintValue || mesh.m_numFacesNodes[nextElementId] != 4))
    // {
    //     // if (previousNewNode == constants::missing::uintValue)
    //     // {
    //     //     throw ConstraintError("previousNewNode is null");
    //     // }

    //     std::cout << " case 5 " << std::endl;
    //     return;
    // }

    if (mesh.m_numFacesNodes[previousElementId] == 4 && (nextElementId == constants::missing::uintValue || (nextElementId != constants::missing::uintValue && mesh.m_numFacesNodes[nextElementId] == 4)))
    // if (mesh.m_numFacesNodes[previousElementId] == 4 && nextElementId != constants::missing::uintValue && mesh.m_numFacesNodes[nextElementId] == 4)
    {
        if (previousNewNode == constants::missing::uintValue)
        {
            Point point = 0.5 * (mesh.Node(previousEdge.first) + mesh.Node(previousEdge.second));
            std::cout << "new node coords: " << point.x << ", " << point.y << std::endl;
            auto [newNodeId, undo] = mesh.InsertNode(point);
            std::cout << "new node: " << newNodeId << std::endl;
            undoActions.Add(std::move(undo));

            edgesToDelete.push_back(previousEdgeId);
            // undoActions.Add(mesh.DeleteEdge(previousEdgeId));

            auto [newEdgeId1, newEdgeUndo1] = mesh.ConnectNodes(previousEdge.first, newNodeId);
            undoActions.Add(std::move(newEdgeUndo1));

            auto [newEdgeId2, newEdgeUndo2] = mesh.ConnectNodes(newNodeId, previousEdge.second);
            undoActions.Add(std::move(newEdgeUndo2));
            previousNewNode = newNodeId;
            // throw ConstraintError("previousNewNode is null");
        }

        Point point = 0.5 * (mesh.Node(edge.first) + mesh.Node(edge.second));
        std::cout << "new node coords: " << point.x << ", " << point.y << std::endl;
        auto [newNodeId, undo] = mesh.InsertNode(point);
        std::cout << "new node: " << newNodeId << std::endl;
        undoActions.Add(std::move(undo));

        edgesToDelete.push_back(edgeId);
        // undoActions.Add(mesh.DeleteEdge(edgeId));

        auto [newEdgeId1, newEdgeUndo1] = mesh.ConnectNodes(edge.first, newNodeId);
        undoActions.Add(std::move(newEdgeUndo1));

        auto [newEdgeId2, newEdgeUndo2] = mesh.ConnectNodes(newNodeId, edge.second);
        undoActions.Add(std::move(newEdgeUndo2));

        auto [newEdgeId3, newEdgeUndo3] = mesh.ConnectNodes(previousNewNode, newNodeId);
        undoActions.Add(std::move(newEdgeUndo2));

        // previousNewNode = newNodeId;

        std::cout << " case 6" << std::endl;
        return;
    }

    // SHould not be here
    throw ConstraintError("Should not be here");
}

void meshkernel::SplitRowColumnOfMesh3::SplitEdges(Mesh2D& mesh, std::vector<UInt>& elementIds, std::vector<UInt>& edgeIds, CompoundUndoAction& undoActions) const
{
    UInt previousNewNode = constants::missing::uintValue;
    UInt currentElement = elementIds[0];
    UInt elementPosition = 0;

    for (UInt currentEdge : edgeIds)
    {
        SplitEdge2(mesh, currentElement, currentEdge, previousNewNode, undoActions);
        currentElement = elementIds[elementPosition];
        ++elementPosition;
    }
}
