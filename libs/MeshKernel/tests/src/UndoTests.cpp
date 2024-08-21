#include <gtest/gtest.h>

#include <iterator>

#include "MeshKernel/UndoActions/AddEdgeAction.hpp"
#include "MeshKernel/UndoActions/AddNodeAction.hpp"
#include "MeshKernel/UndoActions/DeleteEdgeAction.hpp"
#include "MeshKernel/UndoActions/DeleteNodeAction.hpp"
#include "MeshKernel/UndoActions/NoActionUndo.hpp"
#include "MeshKernel/UndoActions/ResetEdgeAction.hpp"
#include "MeshKernel/UndoActions/ResetNodeAction.hpp"
#include "MeshKernel/UndoActions/UndoActionStack.hpp"

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/Operations.hpp" // can delete after test

#include "TestUtils/Definitions.hpp"
#include "TestUtils/MakeMeshes.hpp"

namespace mk = meshkernel;

TEST(UndoTests, AddNodeToMesh)
{
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    mk::Point node(0.5, 0.5);

    auto [nodeId, action] = mesh->InsertNode(node);
    EXPECT_EQ(mk::UndoAction::State::Committed, action->GetState());
    EXPECT_EQ(4, action->NodeId());
    EXPECT_EQ(node.x, action->Node().x);
    EXPECT_EQ(node.y, action->Node().y);
    EXPECT_EQ(5, mesh->GetNumNodes());
    EXPECT_EQ(5, mesh->GetNumValidNodes());
    EXPECT_EQ(node.x, mesh->Node(action->NodeId()).x);
    EXPECT_EQ(node.y, mesh->Node(action->NodeId()).y);
}

TEST(UndoTests, ResetNodeInMesh)
{
    auto mesh = MakeRectangularMeshForTesting(3, 3, 1.0, meshkernel::Projection::cartesian);

    std::vector<mk::Point> originalNodes(mesh->Nodes());

    mk::UInt nodeId = 4;
    mk::Point initialNode = mesh->Node(nodeId);
    mk::Point updatedNode(1.5, 1.5);

    std::unique_ptr<mk::ResetNodeAction> action = mesh->ResetNode(nodeId, updatedNode);
    EXPECT_EQ(mk::UndoAction::State::Committed, action->GetState());
    EXPECT_EQ(nodeId, action->NodeId());
    EXPECT_EQ(updatedNode.x, action->UpdatedNode().x);
    EXPECT_EQ(updatedNode.y, action->UpdatedNode().y);
    EXPECT_EQ(updatedNode.x, mesh->Node(action->NodeId()).x);
    EXPECT_EQ(updatedNode.y, mesh->Node(action->NodeId()).y);
    EXPECT_EQ(initialNode.x, action->InitialNode().x);
    EXPECT_EQ(initialNode.y, action->InitialNode().y);
    EXPECT_EQ(9, mesh->GetNumNodes());
    EXPECT_EQ(9, mesh->GetNumValidNodes());

    action->Restore();
    EXPECT_EQ(mk::UndoAction::State::Restored, action->GetState());

    for (mk::UInt i = 0; i < originalNodes.size(); ++i)
    {
        EXPECT_EQ(originalNodes[i].x, mesh->Node(i).x);
        EXPECT_EQ(originalNodes[i].y, mesh->Node(i).y);
    }

    EXPECT_EQ(initialNode.x, mesh->Node(action->NodeId()).x);
    EXPECT_EQ(initialNode.y, mesh->Node(action->NodeId()).y);
    EXPECT_EQ(9, mesh->GetNumNodes());
    EXPECT_EQ(9, mesh->GetNumValidNodes());
}

TEST(UndoTests, AddNodeThenUndoInMesh)
{
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    mk::Point node(0.5, 0.5);

    auto [nodeId, action] = mesh->InsertNode(node);

    EXPECT_EQ(mk::UndoAction::State::Committed, action->GetState());
    EXPECT_EQ(4, action->NodeId());
    EXPECT_EQ(node.x, action->Node().x);
    EXPECT_EQ(node.y, action->Node().y);
    EXPECT_EQ(5, mesh->GetNumNodes());
    EXPECT_EQ(5, mesh->GetNumValidNodes());
    EXPECT_EQ(node.x, mesh->Node(action->NodeId()).x);
    EXPECT_EQ(node.y, mesh->Node(action->NodeId()).y);

    action->Restore();
    EXPECT_EQ(mk::UndoAction::State::Restored, action->GetState());
    EXPECT_EQ(5, mesh->GetNumNodes());
    EXPECT_EQ(4, mesh->GetNumValidNodes());
    EXPECT_EQ(mk::constants::missing::doubleValue, mesh->Node(action->NodeId()).x);
    EXPECT_EQ(mk::constants::missing::doubleValue, mesh->Node(action->NodeId()).y);
}

TEST(UndoTests, ConnectAddedNodeInMesh)
{
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    mk::Point node(0.5, 0.5);

    auto [nodeId, addNodeAction] = mesh->InsertNode(node);
    auto [edgeId, addEdgeAction] = mesh->ConnectNodes(0, nodeId);

    EXPECT_EQ(mk::UndoAction::State::Committed, addNodeAction->GetState());
    EXPECT_EQ(mk::UndoAction::State::Committed, addEdgeAction->GetState());

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

    auto [nodeId, addNodeAction] = mesh->InsertNode(node);
    auto [edgeId, addEdgeAction] = mesh->ConnectNodes(0, nodeId);

    EXPECT_EQ(mk::UndoAction::State::Committed, addNodeAction->GetState());
    EXPECT_EQ(mk::UndoAction::State::Committed, addEdgeAction->GetState());

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
    EXPECT_EQ(mk::UndoAction::State::Restored, addEdgeAction->GetState());
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

    auto [nodeId, action] = mesh->InsertNode(node);
    EXPECT_EQ(mk::UndoAction::State::Committed, action->GetState());
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

    std::unique_ptr<mk::DeleteNodeAction> deleteAction = mesh->DeleteNode(nodeId);
    EXPECT_EQ(mk::UndoAction::State::Committed, deleteAction->GetState());

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
    EXPECT_EQ(mk::UndoAction::State::Restored, deleteAction->GetState());

    EXPECT_EQ(5, mesh->GetNumNodes());
    EXPECT_EQ(5, mesh->GetNumValidNodes());
    EXPECT_EQ(node.x, mesh->Node(nodeId).x);
    EXPECT_EQ(node.y, mesh->Node(nodeId).y);

    ////////////////////////////////
    // Finally redo the node deletion
    ////////////////////////////////

    deleteAction->Commit();
    EXPECT_EQ(mk::UndoAction::State::Committed, deleteAction->GetState());

    EXPECT_EQ(5, mesh->GetNumNodes());
    EXPECT_EQ(4, mesh->GetNumValidNodes());
    EXPECT_EQ(mk::constants::missing::doubleValue, mesh->Node(nodeId).x);
    EXPECT_EQ(mk::constants::missing::doubleValue, mesh->Node(nodeId).y);
}

TEST(UndoTests, ConnectNodeToCornersInMesh)
{
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    mk::Point node(0.5, 0.5);

    auto [nodeId, nodeAction] = mesh->InsertNode(node);
    auto [edgeId1, edgeAction1] = mesh->ConnectNodes(0, nodeId);
    auto [edgeId2, edgeAction2] = mesh->ConnectNodes(1, nodeId);
    auto [edgeId3, edgeAction3] = mesh->ConnectNodes(2, nodeId);
    auto [edgeId4, edgeAction4] = mesh->ConnectNodes(3, nodeId);

    mesh->Administrate();

    EXPECT_EQ(5, mesh->GetNumValidNodes());
    EXPECT_EQ(8, mesh->GetNumValidEdges());

    // Deleting the node, should also result in deleting of the 4 connecting edges
    std::unique_ptr<mk::DeleteNodeAction> deleteNodeAction = mesh->DeleteNode(nodeId);
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

TEST(UndoTests, MoveNodeMeshTest)
{
    auto mesh = MakeRectangularMeshForTesting(11, 11, 1.0, meshkernel::Projection::cartesian);
    std::vector<mk::Point> originalNodes(mesh->Nodes());

    // Id of node (5.0, 5.0)
    mk::UInt nodeId = 60;
    // Location to where node 60 should be moved.
    mk::Point node(6.5, 6.5);

    std::unique_ptr<mk::UndoAction> action = mesh->MoveNode(node, nodeId);

    std::vector<mk::UInt> expectedNode{38, 48, 49, 50, 58, 59, 60, 61, 62, 70, 71, 72, 82};
    std::vector<double> expectedDispX{0.01207305389, 0.375, 0.8172859213, 0.375, 0.01207305389, 0.8172859213, 1.5, 0.8172859213, 0.01207305389, 0.375, 0.8172859213, 0.375, 0.01207305389};
    std::vector<double> expectedDispY{0.01207305389, 0.375, 0.8172859213, 0.375, 0.01207305389, 0.8172859213, 1.5, 0.8172859213, 0.01207305389, 0.375, 0.8172859213, 0.375, 0.01207305389};

    const double tolerance = 1.0e-8;

    // Check expected nodes have been moved
    for (mk::UInt i = 0; i < expectedNode.size(); ++i)
    {
        EXPECT_NEAR(originalNodes[expectedNode[i]].x + expectedDispX[i], mesh->Node(expectedNode[i]).x, tolerance);
        EXPECT_NEAR(originalNodes[expectedNode[i]].y + expectedDispY[i], mesh->Node(expectedNode[i]).y, tolerance);
    }

    // Check other nodes in mesh have not be touched
    for (mk::UInt i = 0; i < mesh->GetNumNodes(); ++i)
    {
        if (std::ranges::find(expectedNode, i) == expectedNode.end())
        {
            EXPECT_NEAR(originalNodes[i].x, mesh->Node(i).x, tolerance);
            EXPECT_NEAR(originalNodes[i].y, mesh->Node(i).y, tolerance);
        }
    }

    // Restore the original mesh
    action->Restore();

    // All nodes should be as they were before begin moved.
    for (mk::UInt i = 0; i < mesh->GetNumNodes(); ++i)
    {
        EXPECT_EQ(originalNodes[i].x, mesh->Node(i).x);
        EXPECT_EQ(originalNodes[i].y, mesh->Node(i).y);
    }
}

TEST(UndoTests, AddAndRemoveEdgeFromMeshTest)
{
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    mk::UInt start = 0;
    mk::UInt end = 3;

    auto [edgeId, connectNodesAction] = mesh->ConnectNodes(0, 3);

    EXPECT_EQ(edgeId, connectNodesAction->EdgeId());
    EXPECT_EQ(start, connectNodesAction->GetEdge().first);
    EXPECT_EQ(end, connectNodesAction->GetEdge().second);

    EXPECT_EQ(mesh->GetNumEdges(), 5);
    EXPECT_EQ(mesh->GetNumValidEdges(), 5);
    EXPECT_EQ(mesh->GetEdge(4).first, start);
    EXPECT_EQ(mesh->GetEdge(4).second, end);

    std::unique_ptr<mk::DeleteEdgeAction> deleteEdgeAction = mesh->DeleteEdge(connectNodesAction->EdgeId());

    EXPECT_EQ(edgeId, deleteEdgeAction->EdgeId());
    EXPECT_EQ(start, deleteEdgeAction->GetEdge().first);
    EXPECT_EQ(end, deleteEdgeAction->GetEdge().second);

    EXPECT_EQ(mesh->GetNumEdges(), 5);
    EXPECT_EQ(mesh->GetNumValidEdges(), 4);
    EXPECT_EQ(mesh->GetEdge(4).first, mk::constants::missing::uintValue);
    EXPECT_EQ(mesh->GetEdge(4).second, mk::constants::missing::uintValue);

    // Restore the deleted edge
    deleteEdgeAction->Restore();

    // The undo action contents should not change
    EXPECT_EQ(edgeId, deleteEdgeAction->EdgeId());
    EXPECT_EQ(start, deleteEdgeAction->GetEdge().first);
    EXPECT_EQ(end, deleteEdgeAction->GetEdge().second);

    // The edge should be restored
    EXPECT_EQ(mesh->GetNumEdges(), 5);
    EXPECT_EQ(mesh->GetNumValidEdges(), 5);
    EXPECT_EQ(mesh->GetEdge(4).first, start);
    EXPECT_EQ(mesh->GetEdge(4).second, end);

    // Restore the node connection, i.e. remove the added edge
    connectNodesAction->Restore();

    // The contents of the undo action should not change
    EXPECT_EQ(edgeId, connectNodesAction->EdgeId());
    EXPECT_EQ(start, connectNodesAction->GetEdge().first);
    EXPECT_EQ(end, connectNodesAction->GetEdge().second);

    EXPECT_EQ(mesh->GetNumEdges(), 5);
    EXPECT_EQ(mesh->GetNumValidEdges(), 4);
    EXPECT_EQ(mesh->GetEdge(4).first, mk::constants::missing::uintValue);
    EXPECT_EQ(mesh->GetEdge(4).second, mk::constants::missing::uintValue);
}

TEST(UndoTests, AddAndResetEdgeInMeshTest)
{
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    mk::UInt initialStart = 0;
    mk::UInt initialEnd = 3;

    mk::UInt updatedStart = 1;
    mk::UInt updatedEnd = 2;

    auto [edgeId, connectNodesAction] = mesh->ConnectNodes(initialStart, initialEnd);

    EXPECT_EQ(edgeId, connectNodesAction->EdgeId());
    EXPECT_EQ(initialStart, connectNodesAction->GetEdge().first);
    EXPECT_EQ(initialEnd, connectNodesAction->GetEdge().second);

    EXPECT_EQ(mesh->GetNumEdges(), 5);
    EXPECT_EQ(mesh->GetNumValidEdges(), 5);
    EXPECT_EQ(mesh->GetEdge(edgeId).first, initialStart);
    EXPECT_EQ(mesh->GetEdge(edgeId).second, initialEnd);

    std::unique_ptr<mk::ResetEdgeAction> resetEdgeAction = mesh->ResetEdge(edgeId, {updatedStart, updatedEnd});

    EXPECT_EQ(edgeId, resetEdgeAction->EdgeId());
    EXPECT_EQ(initialStart, resetEdgeAction->InitialEdge().first);
    EXPECT_EQ(initialEnd, resetEdgeAction->InitialEdge().second);
    EXPECT_EQ(updatedStart, resetEdgeAction->UpdatedEdge().first);
    EXPECT_EQ(updatedEnd, resetEdgeAction->UpdatedEdge().second);

    EXPECT_EQ(mesh->GetNumEdges(), 5);
    EXPECT_EQ(mesh->GetNumValidEdges(), 5);
    EXPECT_EQ(mesh->GetEdge(edgeId).first, updatedStart);
    EXPECT_EQ(mesh->GetEdge(edgeId).second, updatedEnd);

    ////////////////////////////////
    // Restore the resetting of the edge
    resetEdgeAction->Restore();

    EXPECT_EQ(mesh->GetNumEdges(), 5);
    EXPECT_EQ(mesh->GetNumValidEdges(), 5);
    EXPECT_EQ(mesh->GetEdge(edgeId).first, initialStart);
    EXPECT_EQ(mesh->GetEdge(edgeId).second, initialEnd);

    ////////////////////////////////
    // Remove the added edge
    connectNodesAction->Restore();

    EXPECT_EQ(mesh->GetNumEdges(), 5);
    EXPECT_EQ(mesh->GetNumValidEdges(), 4);
    EXPECT_EQ(mesh->GetEdge(edgeId).first, mk::constants::missing::uintValue);
    EXPECT_EQ(mesh->GetEdge(edgeId).second, mk::constants::missing::uintValue);
}

TEST(UndoTests, AddAndResetEdgeInMeshTestUsingStack)
{
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);
    mk::UndoActionStack undoActionStack;

    mk::UInt initialStart = 0;
    mk::UInt initialEnd = 3;

    mk::UInt updatedStart = 1;
    mk::UInt updatedEnd = 2;

    // Connect diagonally opposite nodes
    auto [edgeId, connectNodesAction] = mesh->ConnectNodes(initialStart, initialEnd);
    undoActionStack.Add(std::move(connectNodesAction));

    EXPECT_EQ(mesh->GetNumEdges(), 5);
    EXPECT_EQ(mesh->GetNumValidEdges(), 5);
    EXPECT_EQ(mesh->GetEdge(edgeId).first, initialStart);
    EXPECT_EQ(mesh->GetEdge(edgeId).second, initialEnd);

    // Reset to other set of diagonnally opposite nodes
    undoActionStack.Add(mesh->ResetEdge(edgeId, {updatedStart, updatedEnd}));

    EXPECT_EQ(mesh->GetNumEdges(), 5);
    EXPECT_EQ(mesh->GetNumValidEdges(), 5);
    EXPECT_EQ(mesh->GetEdge(edgeId).first, updatedStart);
    EXPECT_EQ(mesh->GetEdge(edgeId).second, updatedEnd);

    ////////////////////////////////
    // Restore the resetting of the edge
    EXPECT_TRUE(undoActionStack.Undo());

    EXPECT_EQ(mesh->GetNumEdges(), 5);
    EXPECT_EQ(mesh->GetNumValidEdges(), 5);
    EXPECT_EQ(mesh->GetEdge(edgeId).first, initialStart);
    EXPECT_EQ(mesh->GetEdge(edgeId).second, initialEnd);

    ////////////////////////////////
    // Remove the added edge
    EXPECT_TRUE(undoActionStack.Undo());

    EXPECT_EQ(mesh->GetNumEdges(), 5);
    EXPECT_EQ(mesh->GetNumValidEdges(), 4);
    EXPECT_EQ(mesh->GetEdge(edgeId).first, mk::constants::missing::uintValue);
    EXPECT_EQ(mesh->GetEdge(edgeId).second, mk::constants::missing::uintValue);

    ////////////////////////////////
    // Nothing else to undo
    EXPECT_FALSE(undoActionStack.Undo());
}

TEST(UndoTests, AddNoActionUndoTest)
{
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);
    mk::UndoActionStack undoActionStack;

    mk::UInt initialStart = 0;
    mk::UInt initialEnd = 3;

    // Connect diagonally opposite nodes
    auto [edgeId, connectNodesAction] = mesh->ConnectNodes(initialStart, initialEnd);
    undoActionStack.Add(std::move(connectNodesAction));

    EXPECT_EQ(mesh->GetNumEdges(), 5);
    EXPECT_EQ(mesh->GetNumValidEdges(), 5);
    EXPECT_EQ(mesh->GetEdge(edgeId).first, initialStart);
    EXPECT_EQ(mesh->GetEdge(edgeId).second, initialEnd);

    // Add the no action undo
    undoActionStack.Add(meshkernel::NoActionUndo::Create());

    // Mesh should be same as before adding the action
    EXPECT_EQ(mesh->GetNumEdges(), 5);
    EXPECT_EQ(mesh->GetNumValidEdges(), 5);
    EXPECT_EQ(mesh->GetEdge(edgeId).first, initialStart);
    EXPECT_EQ(mesh->GetEdge(edgeId).second, initialEnd);

    EXPECT_TRUE(undoActionStack.Undo());

    // Mesh should be same as before undo
    EXPECT_EQ(mesh->GetNumEdges(), 5);
    EXPECT_EQ(mesh->GetNumValidEdges(), 5);
    EXPECT_EQ(mesh->GetEdge(edgeId).first, initialStart);
    EXPECT_EQ(mesh->GetEdge(edgeId).second, initialEnd);
}
