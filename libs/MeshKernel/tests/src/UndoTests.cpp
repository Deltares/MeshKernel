#include <gtest/gtest.h>

#include "MeshKernel/AddEdgeAction.hpp"
#include "MeshKernel/AddNodeAction.hpp"
#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Mesh2D.hpp"

#include "TestUtils/Definitions.hpp"
#include "TestUtils/MakeMeshes.hpp"

namespace mk = meshkernel;

TEST(UndoTests, ExpectedEdge)
{
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    mk::UInt edgeId = 10;
    mk::UInt start = 20;
    mk::UInt end = 30;

    std::unique_ptr<mk::AddEdgeAction> action = mk::AddEdgeAction::Create(*mesh, edgeId, start, end);
    EXPECT_EQ(edgeId, action->EdgeId());
    EXPECT_EQ(start, action->GetEdge().first);
    EXPECT_EQ(end, action->GetEdge().second);
}

TEST(UndoTests, ExpectedNode)
{
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    mk::UInt nodeId = 10;
    mk::Point node(0.5, 0.5);

    std::unique_ptr<mk::AddNodeAction> action = mk::AddNodeAction::Create(*mesh, nodeId, node);
    EXPECT_EQ(nodeId, action->NodeId());
    EXPECT_EQ(node.x, action->Node().x);
    EXPECT_EQ(node.y, action->Node().y);
}

TEST(UndoTests, AddNodeToMesh)
{
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    mk::Point node(0.5, 0.5);

    auto [nodeId, action] = mesh->InsertNode2(node);
    EXPECT_EQ(mk::UndoAction::Committed, action->State());
    EXPECT_EQ(4, action->NodeId());
    EXPECT_EQ(node.x, action->Node().x);
    EXPECT_EQ(node.y, action->Node().y);
    EXPECT_EQ(5, mesh->GetNumNodes());
    EXPECT_EQ(5, mesh->GetNumValidNodes());
    EXPECT_EQ(node.x, mesh->Node(action->NodeId()).x);
    EXPECT_EQ(node.y, mesh->Node(action->NodeId()).y);
}

TEST(UndoTests, AddNodeThenUndoInMesh)
{
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    mk::Point node(0.5, 0.5);

    auto [nodeId, action] = mesh->InsertNode2(node);

    EXPECT_EQ(mk::UndoAction::Committed, action->State());
    EXPECT_EQ(4, action->NodeId());
    EXPECT_EQ(node.x, action->Node().x);
    EXPECT_EQ(node.y, action->Node().y);
    EXPECT_EQ(5, mesh->GetNumNodes());
    EXPECT_EQ(5, mesh->GetNumValidNodes());
    EXPECT_EQ(node.x, mesh->Node(action->NodeId()).x);
    EXPECT_EQ(node.y, mesh->Node(action->NodeId()).y);

    action->Restore();
    EXPECT_EQ(mk::UndoAction::Restored, action->State());
    EXPECT_EQ(5, mesh->GetNumNodes());
    EXPECT_EQ(4, mesh->GetNumValidNodes());
    EXPECT_EQ(mk::constants::missing::doubleValue, mesh->Node(action->NodeId()).x);
    EXPECT_EQ(mk::constants::missing::doubleValue, mesh->Node(action->NodeId()).y);
}

TEST(UndoTests, ConnectAddedNodeInMesh)
{
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    mk::Point node(0.5, 0.5);

    auto [nodeId, addNodeAction] = mesh->InsertNode2(node);
    auto [edgeId, addEdgeAction] = mesh->ConnectNodes2(0, nodeId);

    EXPECT_EQ(mk::UndoAction::Committed, addNodeAction->State());
    EXPECT_EQ(mk::UndoAction::Committed, addEdgeAction->State());

    EXPECT_EQ(4, addNodeAction->NodeId());
    EXPECT_EQ(node.x, addNodeAction->Node().x);
    EXPECT_EQ(node.y, addNodeAction->Node().y);
    EXPECT_EQ(5, mesh->GetNumNodes());
    EXPECT_EQ(5, mesh->GetNumValidNodes());
    EXPECT_EQ(node.x, mesh->Node(addNodeAction->NodeId()).x);
    EXPECT_EQ(node.y, mesh->Node(addNodeAction->NodeId()).y);

    EXPECT_EQ(4, addEdgeAction->EdgeId());
    EXPECT_EQ(0, addEdgeAction->GetEdge().first);
    EXPECT_EQ(nodeId, addEdgeAction->GetEdge().second);
    EXPECT_EQ(addEdgeAction->GetEdge().first, mesh->GetEdge(edgeId).first);
    EXPECT_EQ(addEdgeAction->GetEdge().second, mesh->GetEdge(edgeId).second);
}

TEST(UndoTests, ConnectAddedNodeThenUndoInMesh)
{
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    mk::Point node(0.5, 0.5);

    auto [nodeId, addNodeAction] = mesh->InsertNode2(node);
    auto [edgeId, addEdgeAction] = mesh->ConnectNodes2(0, nodeId);

    EXPECT_EQ(mk::UndoAction::Committed, addNodeAction->State());
    EXPECT_EQ(mk::UndoAction::Committed, addEdgeAction->State());

    EXPECT_EQ(4, addNodeAction->NodeId());
    EXPECT_EQ(node.x, addNodeAction->Node().x);
    EXPECT_EQ(node.y, addNodeAction->Node().y);
    EXPECT_EQ(5, mesh->GetNumNodes());
    EXPECT_EQ(5, mesh->GetNumValidNodes());
    EXPECT_EQ(node.x, mesh->Node(addNodeAction->NodeId()).x);
    EXPECT_EQ(node.y, mesh->Node(addNodeAction->NodeId()).y);

    EXPECT_EQ(4, addEdgeAction->EdgeId());
    EXPECT_EQ(0, addEdgeAction->GetEdge().first);
    EXPECT_EQ(nodeId, addEdgeAction->GetEdge().second);
    EXPECT_EQ(addEdgeAction->GetEdge().first, mesh->GetEdge(edgeId).first);
    EXPECT_EQ(addEdgeAction->GetEdge().second, mesh->GetEdge(edgeId).second);

    addEdgeAction->Restore();
    EXPECT_EQ(mk::UndoAction::Restored, addEdgeAction->State());
    EXPECT_EQ(nodeId, addEdgeAction->GetEdge().second);
    EXPECT_EQ(mk::constants::missing::uintValue, mesh->GetEdge(edgeId).first);
    EXPECT_EQ(mk::constants::missing::uintValue, mesh->GetEdge(edgeId).second);
}

TEST(UndoTests, DeleteAddedNodeToMesh)
{
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    mk::Point node(0.5, 0.5);

    ////////////////////////////////
    // First add the node
    ////////////////////////////////

    auto [nodeId, action] = mesh->InsertNode2(node);
    EXPECT_EQ(mk::UndoAction::Committed, action->State());
    EXPECT_EQ(4, action->NodeId());
    EXPECT_EQ(node.x, action->Node().x);
    EXPECT_EQ(node.y, action->Node().y);
    EXPECT_EQ(5, mesh->GetNumNodes());
    EXPECT_EQ(5, mesh->GetNumValidNodes());
    EXPECT_EQ(node.x, mesh->Node(action->NodeId()).x);
    EXPECT_EQ(node.y, mesh->Node(action->NodeId()).y);

    ////////////////////////////////
    // Next delete the node
    ////////////////////////////////

    std::unique_ptr<mk::DeleteNodeAction> deleteAction = mesh->DeleteNode2(nodeId);
    EXPECT_EQ(mk::UndoAction::Committed, deleteAction->State());

    EXPECT_EQ(nodeId, deleteAction->NodeId());
    EXPECT_EQ(node.x, deleteAction->Node().x);
    EXPECT_EQ(node.y, deleteAction->Node().y);
    EXPECT_EQ(5, mesh->GetNumNodes());
    EXPECT_EQ(4, mesh->GetNumValidNodes());
    EXPECT_EQ(mk::constants::missing::doubleValue, mesh->Node(nodeId).x);
    EXPECT_EQ(mk::constants::missing::doubleValue, mesh->Node(nodeId).y);

    ////////////////////////////////
    // Undo the node deletion
    ////////////////////////////////

    deleteAction->Restore();
    EXPECT_EQ(mk::UndoAction::Restored, deleteAction->State());

    EXPECT_EQ(5, mesh->GetNumNodes());
    EXPECT_EQ(5, mesh->GetNumValidNodes());
    EXPECT_EQ(node.x, mesh->Node(nodeId).x);
    EXPECT_EQ(node.y, mesh->Node(nodeId).y);

    ////////////////////////////////
    // Finally redo the node deletion
    ////////////////////////////////

    deleteAction->Commit();
    EXPECT_EQ(mk::UndoAction::Committed, deleteAction->State());

    EXPECT_EQ(5, mesh->GetNumNodes());
    EXPECT_EQ(4, mesh->GetNumValidNodes());
    EXPECT_EQ(mk::constants::missing::doubleValue, mesh->Node(nodeId).x);
    EXPECT_EQ(mk::constants::missing::doubleValue, mesh->Node(nodeId).y);
}

TEST(UndoTests, ConnectNodeToCornersInMesh)
{
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    mk::Point node(0.5, 0.5);

    auto [nodeId, nodeAction] = mesh->InsertNode2(node);
    auto [edgeId1, edgeAction1] = mesh->ConnectNodes2(0, nodeId);
    auto [edgeId2, edgeAction2] = mesh->ConnectNodes2(1, nodeId);
    auto [edgeId3, edgeAction3] = mesh->ConnectNodes2(2, nodeId);
    auto [edgeId4, edgeAction4] = mesh->ConnectNodes2(3, nodeId);

    mesh->Administrate();

    EXPECT_EQ(5, mesh->GetNumValidNodes());
    EXPECT_EQ(8, mesh->GetNumValidEdges());

    // Deleting the node, should also result in deleting of the 4 connecting edges
    std::unique_ptr<mk::DeleteNodeAction> deleteNodeAction = mesh->DeleteNode2(nodeId);
    mesh->Administrate();

    EXPECT_EQ(4, mesh->GetNumValidNodes());
    EXPECT_EQ(4, mesh->GetNumValidEdges());

    EXPECT_EQ(mk::constants::missing::doubleValue, mesh->Node(nodeId).x);
    EXPECT_EQ(mk::constants::missing::doubleValue, mesh->Node(nodeId).y);

    EXPECT_EQ(mk::constants::missing::uintValue, mesh->GetEdge(edgeId1).first);
    EXPECT_EQ(mk::constants::missing::uintValue, mesh->GetEdge(edgeId1).second);

    EXPECT_EQ(mk::constants::missing::uintValue, mesh->GetEdge(edgeId2).first);
    EXPECT_EQ(mk::constants::missing::uintValue, mesh->GetEdge(edgeId2).second);

    EXPECT_EQ(mk::constants::missing::uintValue, mesh->GetEdge(edgeId3).first);
    EXPECT_EQ(mk::constants::missing::uintValue, mesh->GetEdge(edgeId3).second);

    EXPECT_EQ(mk::constants::missing::uintValue, mesh->GetEdge(edgeId4).first);
    EXPECT_EQ(mk::constants::missing::uintValue, mesh->GetEdge(edgeId4).second);

    deleteNodeAction->Restore();

    EXPECT_EQ(5, mesh->GetNumValidNodes());
    EXPECT_EQ(8, mesh->GetNumValidEdges());

    EXPECT_EQ(node.x, mesh->Node(nodeId).x);
    EXPECT_EQ(node.y, mesh->Node(nodeId).y);

    EXPECT_EQ(edgeAction1->GetEdge().first, mesh->GetEdge(edgeId1).first);
    EXPECT_EQ(edgeAction1->GetEdge().second, mesh->GetEdge(edgeId1).second);

    EXPECT_EQ(edgeAction2->GetEdge().first, mesh->GetEdge(edgeId2).first);
    EXPECT_EQ(edgeAction2->GetEdge().second, mesh->GetEdge(edgeId2).second);

    EXPECT_EQ(edgeAction3->GetEdge().first, mesh->GetEdge(edgeId3).first);
    EXPECT_EQ(edgeAction3->GetEdge().second, mesh->GetEdge(edgeId3).second);

    EXPECT_EQ(edgeAction4->GetEdge().first, mesh->GetEdge(edgeId4).first);
    EXPECT_EQ(edgeAction4->GetEdge().second, mesh->GetEdge(edgeId4).second);
}
