//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2025.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation version 3.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// contact: delft3d.support@deltares.nl
// Stichting Deltares
// P.O. Box 177
// 2600 MH Delft, The Netherlands
//
// All indications and logos of, and references to, "Delft3D" and "Deltares"
// are registered trademarks of Stichting Deltares, and remain the property of
// Stichting Deltares. All rights reserved.
//
//------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <iterator>

#include "MeshKernel/Polygons.hpp"
#include "MeshKernel/UndoActions/AddEdgeAction.hpp"
#include "MeshKernel/UndoActions/AddNodeAction.hpp"
#include "MeshKernel/UndoActions/DeleteEdgeAction.hpp"
#include "MeshKernel/UndoActions/DeleteNodeAction.hpp"
#include "MeshKernel/UndoActions/ResetEdgeAction.hpp"
#include "MeshKernel/UndoActions/ResetNodeAction.hpp"
#include "MeshKernel/UndoActions/UndoActionStack.hpp"
#include "MeshKernel/Utilities/Utilities.hpp"

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Mesh2D.hpp"

#include "TestUtils/Definitions.hpp"
#include "TestUtils/MakeMeshes.hpp"

namespace mk = meshkernel;

TEST(CompoundUndoTests, MergeTwoNodesInMesh)
{
    auto mesh = MakeRectangularMeshForTesting(3, 3, 1.0, meshkernel::Projection::cartesian);
    mk::UndoActionStack undoActionStack;

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

TEST(CompoundUndoTests, MergeNodesInPolygonInMesh)
{
    auto mesh = MakeRectangularMeshForTesting(11, 11, 1.0, meshkernel::Projection::cartesian);
    const std::vector<mk::Point> originalNodes(mesh->Nodes());
    const std::vector<mk::Edge> originalEdges(mesh->Edges());

    mk::UndoActionStack undoActionStack;

    std::vector<mk::Point> polygonNodes{{3.5, 3.5}, {6.5, 3.5}, {6.5, 6.5}, {3.5, 6.5}, {3.5, 3.5}};
    mk::Polygons polygons(polygonNodes, mesh->m_projection);

    EXPECT_EQ(mesh->GetNumNodes(), originalNodes.size());
    EXPECT_EQ(mesh->GetNumValidNodes(), originalNodes.size());
    EXPECT_EQ(mesh->GetNumEdges(), originalEdges.size());
    EXPECT_EQ(mesh->GetNumValidEdges(), originalEdges.size());

    undoActionStack.Add(mesh->MergeNodesInPolygon(polygons, 1.5));

    EXPECT_EQ(mesh->GetNumNodes(), originalNodes.size());
    EXPECT_EQ(mesh->GetNumValidNodes(), 113);
    EXPECT_EQ(mesh->GetNumEdges(), originalEdges.size());
    EXPECT_EQ(mesh->GetNumValidEdges(), 208);

    // Undoing merge should restore the original mesh.
    undoActionStack.Undo();

    EXPECT_EQ(mesh->GetNumValidEdges(), originalEdges.size());

    for (mk::UInt i = 0; i < originalNodes.size(); ++i)
    {
        EXPECT_EQ(originalNodes[i].x, mesh->Node(i).x);
        EXPECT_EQ(originalNodes[i].y, mesh->Node(i).y);
    }

    for (mk::UInt i = 0; i < originalEdges.size(); ++i)
    {
        EXPECT_EQ(originalEdges[i].first, mesh->GetEdge(i).first);
        EXPECT_EQ(originalEdges[i].second, mesh->GetEdge(i).second);
    }
}
