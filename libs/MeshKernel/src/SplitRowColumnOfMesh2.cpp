#include "MeshKernel/SplitRowColumnOfMesh2.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Operations.hpp"

std::unique_ptr<meshkernel::UndoAction> meshkernel::SplitRowColumnOfMesh2::Compute(Mesh2D& mesh, const UInt edgeId) const
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

    // std::vector<UInt> elementIds;
    // std::vector<UInt> edgeIds;
    // elementIds.reserve(initialNumberOfElements);
    // edgeIds.reserve(initialNumberOfElements);

    std::unique_ptr<CompoundUndoAction> undoActions = CompoundUndoAction::Create();

    std::cout << "edge is valid: " << std::boolalpha << IsValidEdge(mesh, edgeId) << "  " << std::boolalpha << MayBeSplit(mesh, edgeId) << std::endl;

    if (IsValidEdge(mesh, edgeId) && MayBeSplit(mesh, edgeId))
    {
        SplitEdge(mesh, edgeId, *undoActions);
    }

    mesh.Administrate();
    return undoActions;
}

bool meshkernel::SplitRowColumnOfMesh2::IsValidEdge(const Mesh2D& mesh, const UInt edgeId) const
{
    const Edge& edge = mesh.GetEdge(edgeId);

    return meshkernel::IsValidEdge(edge) &&
           mesh.Node(edge.first).IsValid() &&
           mesh.Node(edge.second).IsValid();
}

bool meshkernel::SplitRowColumnOfMesh2::MayBeSplit(const Mesh2D& mesh, const UInt edgeId) const
{
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

    // return (edge[0] == constants::missing::uintValue ? true : mesh.m_numFacesNodes[edge[0]] == 4) ||
    //        (edge[1] == constants::missing::uintValue ? true : mesh.m_numFacesNodes[edge[1]] == 4);
}

void meshkernel::SplitRowColumnOfMesh2::SplitEdge(Mesh2D& mesh, const UInt startEdgeId, CompoundUndoAction& undoActions [[maybe_unused]]) const
{

    // TODO need to check for circular connected element loops
    // Also should not do the reverse loop for circular connected element loops

    UInt previousNewNode = constants::missing::uintValue;
    UInt firstPreviousNewNode = constants::missing::uintValue;

    const UInt initialNumberOfElements = static_cast<UInt>(std::sqrt(static_cast<double>(mesh.GetNumFaces())));
    std::vector<UInt> edgesToDelete;
    edgesToDelete.reserve(initialNumberOfElements);

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
        [[maybe_unused]] bool firstIteration = true;

        UInt iterationCount = 0;

        // Exit also it current is same as first (except when it is the first)
        while ((elementId != constants::missing::uintValue && mesh.m_numFacesNodes[elementId] == 4) || (!firstIteration && previousNewNode != constants::missing::uintValue))
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
                    SplitEdge3(mesh, elementId, edgeId, previousNewNode, undoActions, edgesToDelete);
                }
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

void meshkernel::SplitRowColumnOfMesh2::CollectElementIdsToSplit(const Mesh2D& mesh, const UInt edgeId, std::vector<UInt>& elementIds, std::vector<UInt>& edgeIds) const
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

        [[maybe_unused]] UInt firstElementId = mesh.m_edgesFaces[edgeId][i];
        UInt elementId = mesh.m_edgesFaces[edgeId][i];
        UInt currentEdgeId = edgeId;
        [[maybe_unused]] bool firstIteration = true;

        // Exit also it current is same as first (except when it is the first)
        while (elementId != constants::missing::uintValue && mesh.m_numFacesNodes[elementId] == 4)
        {

            // if (!firstIteration && (elementId == firstElementId || currentEdgeId == edgeId))
            // {
            //     currentEdgeId = constants::missing::uintValue;
            //     elementId = constants::missing::uintValue;
            //     break;
            // }

            elementIds.push_back(elementId);
            edgeIds.push_back(currentEdgeId);

            std::cout << "element: " << elementId << "  " << currentEdgeId << "   ";
            GetNextElement(mesh, currentEdgeId, elementId);
            std::cout << "  next element: " << elementId << "  " << currentEdgeId << "   " << std::endl;

            // if (elementId == firstElementId)
            // {
            //     currentEdgeId = constants::missing::uintValue;
            //     elementId = constants::missing::uintValue;
            // }

            firstIteration = false;
        }

        if (currentEdgeId != constants::missing::uintValue)
        {
            edgeIds.push_back(currentEdgeId);
        }

        if (elementId != constants::missing::uintValue)
        {
            std::cout << "terminator element: " << elementId << "  " << currentEdgeId << std::endl;
            // terminator element
        }
    }

    std::cout << "list size: " << elementIds.size() << "  " << edgeIds.size() << std::endl;
}

void meshkernel::SplitRowColumnOfMesh2::GetNextElement(const Mesh2D& mesh, UInt& edgeId, UInt& elementId) const
{
    if (elementId == constants::missing::uintValue)
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
    // UInt oppositeEdgeIndex = edgeIndex < 2 ? edgeIndex + 2 : edgeIndex - 2;
    edgeId = mesh.m_facesEdges[elementId][oppositeEdgeIndex];

    const std::array<UInt, 2>& edgesFaces = mesh.m_edgesFaces[edgeId];

    elementId = edgesFaces[0] == elementId ? edgesFaces[1] : edgesFaces[0];
    // std::cout << "edge info :" << edgeId << "  " << edgeIndex << "  " << oppositeEdgeIndex << "  " << elementId << std::endl;
}

void meshkernel::SplitRowColumnOfMesh2::SplitEdge(Mesh2D& mesh, [[maybe_unused]] UInt elementId, UInt edgeId, UInt& previousNewNode, CompoundUndoAction& undoActions) const
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

void meshkernel::SplitRowColumnOfMesh2::SplitEdge2(Mesh2D& mesh, [[maybe_unused]] UInt elementId, UInt edgeId, UInt& previousNewNode, CompoundUndoAction& undoActions) const
{
    std::cout << "SplitRowColumnOfMesh2::SplitEdge2 " << edgeId << std::endl;

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

void meshkernel::SplitRowColumnOfMesh2::SplitEdge3(Mesh2D& mesh, UInt elementId, UInt edgeId, UInt& previousNewNode, CompoundUndoAction& undoActions, std::vector<UInt>& edgesToDelete) const
{
    std::cout << "SplitRowColumnOfMesh2::SplitEdge3 " << elementId << "  " << edgeId << std::endl;

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
    const UInt previousEdgeId = mesh.m_facesEdges[elementId][oppositeEdgeIndex];
    const std::array<UInt, 2>& previousEdgeFaces = mesh.m_edgesFaces[previousEdgeId];
    const Edge& previousEdge = mesh.GetEdge(previousEdgeId);

    UInt previousElementId = previousEdgeFaces[0] + previousEdgeFaces[1] - elementId;

    const std::array<UInt, 2>& edgeFaces = mesh.m_edgesFaces[edgeId];
    const UInt nextElementId = edgeFaces[0] + edgeFaces[1] - elementId;

    std::cout << "elements " << previousElementId << "  " << elementId << "  " << nextElementId << std::endl;

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

        previousNewNode = newNodeId;

        std::cout << " case 6" << std::endl;
        return;
    }

    // SHould not be here
    throw ConstraintError("Should not be here");
}

void meshkernel::SplitRowColumnOfMesh2::SplitEdges(Mesh2D& mesh, std::vector<UInt>& elementIds, std::vector<UInt>& edgeIds, CompoundUndoAction& undoActions) const
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
