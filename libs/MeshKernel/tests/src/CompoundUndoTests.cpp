#include <gtest/gtest.h>

#include <iterator>

#include "MeshKernel/AddEdgeAction.hpp"
#include "MeshKernel/AddNodeAction.hpp"
#include "MeshKernel/DeleteEdgeAction.hpp"
#include "MeshKernel/DeleteNodeAction.hpp"
#include "MeshKernel/MoveNodeAction.hpp"
#include "MeshKernel/ResetEdgeAction.hpp"
#include "MeshKernel/ResetNodeAction.hpp"
#include "MeshKernel/UndoActionStack.hpp"

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Mesh2D.hpp"

#include "TestUtils/Definitions.hpp"
#include "TestUtils/MakeMeshes.hpp"

namespace mk = meshkernel;

TEST(CompoundUndoTests, ResetNodeInMesh)
{
    auto mesh = MakeRectangularMeshForTesting(3, 3, 1.0, meshkernel::Projection::cartesian);
    mk::UndoActionStack undoActionStack;

    std::vector<mk::Point> originalNodes(mesh->Nodes());

    mk::Point node(1.5, 1.5);

    auto [nodeId, insertNodeAction] = mesh->InsertNode(node);
    undoActionStack.Add(std::move(insertNodeAction));
    EXPECT_EQ(10, mesh->GetNumNodes());
    EXPECT_EQ(10, mesh->GetNumValidNodes());

    // Connect nodes 8 and 4 to newly added node
    auto [edgeId8, connectNodeAction8] = mesh->ConnectNodes(nodeId, 8);
    undoActionStack.Add(std::move(connectNodeAction8));
    auto [edgeId4, connectNodeAction4] = mesh->ConnectNodes(nodeId, 4);
    undoActionStack.Add(std::move(connectNodeAction4));
    mesh->Administrate();

    EXPECT_EQ(14, mesh->GetNumEdges());
    EXPECT_EQ(14, mesh->GetNumValidEdges());

    undoActionStack.Add(mesh->MergeTwoNodes(4, nodeId));

    EXPECT_EQ(14, mesh->GetNumEdges());
    // Edge is invalid because node 4 is invalid
    EXPECT_EQ(13, mesh->GetNumValidEdges());
    EXPECT_EQ(10, mesh->GetNumNodes());
    EXPECT_EQ(9, mesh->GetNumValidNodes());
    EXPECT_EQ(mesh->Node(4).x, mk::constants::missing::doubleValue);
    EXPECT_EQ(mesh->Node(4).y, mk::constants::missing::doubleValue);

    EXPECT_EQ(mesh->Node(9).x, node.x);
    EXPECT_EQ(mesh->Node(9).y, node.y);

    mesh->Administrate();

    ////////////////////////////////
    // Undo node merging
    undoActionStack.Undo();

    EXPECT_EQ(14, mesh->GetNumEdges());
    EXPECT_EQ(14, mesh->GetNumValidEdges());
    EXPECT_EQ(10, mesh->GetNumNodes());
    EXPECT_EQ(10, mesh->GetNumValidNodes());
    EXPECT_EQ(mesh->Node(4).x, 1.0);
    EXPECT_EQ(mesh->Node(4).y, 1.0);

    EXPECT_EQ(mesh->Node(9).x, node.x);
    EXPECT_EQ(mesh->Node(9).y, node.y);
}
