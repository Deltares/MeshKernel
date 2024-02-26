#include <gtest/gtest.h>

#include <iterator>

#include "MeshKernel/AddEdgeAction.hpp"
#include "MeshKernel/AddNodeAction.hpp"
#include "MeshKernel/DeleteEdgeAction.hpp"
#include "MeshKernel/DeleteNodeAction.hpp"
#include "MeshKernel/MoveNodeAction.hpp"

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

TEST(UndoTests, MoveNodeActionTest)
{
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    std::unique_ptr<mk::MoveNodeAction> action = mk::MoveNodeAction::Create(*mesh);

    action->AddDisplacement(1, 0.25, 0.25);
    action->AddDisplacement(2, 0.5, 0.25);
    action->AddDisplacement(3, 0.25, 0.5);
    action->AddDisplacement(4, 0.5, 0.5);
    action->AddDisplacement(23, -1.5, -1.5);

    auto iter = action->begin();

    ASSERT_EQ(std::distance(action->begin(), action->end()), 5);

    EXPECT_EQ(iter->m_nodeId, 1);
    EXPECT_EQ(iter->m_displacement.x(), 0.25);
    EXPECT_EQ(iter->m_displacement.y(), 0.25);

    ++iter;
    EXPECT_EQ(iter->m_nodeId, 2);
    EXPECT_EQ(iter->m_displacement.x(), 0.5);
    EXPECT_EQ(iter->m_displacement.y(), 0.25);

    ++iter;
    EXPECT_EQ(iter->m_nodeId, 3);
    EXPECT_EQ(iter->m_displacement.x(), 0.25);
    EXPECT_EQ(iter->m_displacement.y(), 0.5);

    ++iter;
    EXPECT_EQ(iter->m_nodeId, 4);
    EXPECT_EQ(iter->m_displacement.x(), 0.5);
    EXPECT_EQ(iter->m_displacement.y(), 0.5);

    ++iter;
    EXPECT_EQ(iter->m_nodeId, 23);
    EXPECT_EQ(iter->m_displacement.x(), -1.5);
    EXPECT_EQ(iter->m_displacement.y(), -1.5);

    ++iter;
    EXPECT_TRUE(iter == action->end());
}

TEST(UndoTests, AddNodeToMesh)
{
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    mk::Point node(0.5, 0.5);

    auto [nodeId, action] = mesh->InsertNode(node);
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

    auto [nodeId, action] = mesh->InsertNode(node);

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

    auto [nodeId, addNodeAction] = mesh->InsertNode(node);
    auto [edgeId, addEdgeAction] = mesh->ConnectNodes(0, nodeId);

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

    auto [nodeId, addNodeAction] = mesh->InsertNode(node);
    auto [edgeId, addEdgeAction] = mesh->ConnectNodes(0, nodeId);

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

    auto [nodeId, action] = mesh->InsertNode(node);
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

    std::unique_ptr<mk::DeleteNodeAction> deleteAction = mesh->DeleteNode(nodeId);
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
    const size_t size = 13;

    std::vector<mk::Point> originalNodes(mesh->Nodes());

    mk::Point node(6.5, 6.5);
    // Id of node (5.0, 5.0)
    mk::UInt nodeId = 60;

    std::unique_ptr<mk::MoveNodeAction> action = mesh->MoveNode(node, nodeId);

    ASSERT_EQ(static_cast<size_t>(std::distance(action->begin(), action->end())), size);

    std::vector<mk::UInt> expectedNode{38, 48, 49, 50, 58, 59, 60, 61, 62, 70, 71, 72, 82};
    std::vector<double> expectedDispX{0.01207305389, 0.375, 0.8172859213, 0.375, 0.01207305389, 0.8172859213, 1.5, 0.8172859213, 0.01207305389, 0.375, 0.8172859213, 0.375, 0.01207305389};
    std::vector<double> expectedDispY{0.01207305389, 0.375, 0.8172859213, 0.375, 0.01207305389, 0.8172859213, 1.5, 0.8172859213, 0.01207305389, 0.375, 0.8172859213, 0.375, 0.01207305389};

    size_t count = 0;

    const double tolerance = 1.0e-8;

    // Check contents of the move node action
    for (auto iter = action->begin(); iter != action->end(); ++iter)
    {
        EXPECT_EQ(iter->m_nodeId, expectedNode[count]);
        EXPECT_NEAR(iter->m_displacement.x(), expectedDispX[count], tolerance);
        EXPECT_NEAR(iter->m_displacement.y(), expectedDispY[count], tolerance);
        ++count;
    }

    // Check that the action has been applied to the mesh
    for (auto iter = action->begin(); iter != action->end(); ++iter)
    {
        mk::Point expectedPoint = originalNodes[iter->m_nodeId] + iter->m_displacement;
        EXPECT_NEAR(mesh->Node(iter->m_nodeId).x, expectedPoint.x, tolerance);
        EXPECT_NEAR(mesh->Node(iter->m_nodeId).y, expectedPoint.y, tolerance);
    }

    // Restore the original mesh
    action->Restore();

    // Check that the mesh nodes are the same as in the original mesh.
    for (auto iter = action->begin(); iter != action->end(); ++iter)
    {
        mk::Point expectedPoint = originalNodes[iter->m_nodeId];
        EXPECT_NEAR(mesh->Node(iter->m_nodeId).x, expectedPoint.x, tolerance);
        EXPECT_NEAR(mesh->Node(iter->m_nodeId).y, expectedPoint.y, tolerance);
    }
}
