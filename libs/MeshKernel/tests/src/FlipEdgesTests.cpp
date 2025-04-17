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

#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/FlipEdges.hpp>
#include <MeshKernel/LandBoundaries.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Polygons.hpp>
#include <TestUtils/Definitions.hpp>
#include <TestUtils/MakeMeshes.hpp>

TEST(FlipEdges, FlipEdgesWithLandBoundary)
{
    // 1 Setup
    auto mesh = MakeRectangularMeshForTesting(3, 3, 10, meshkernel::Projection::cartesian, {0.0, 0.0});

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    // set landboundaries
    auto polygon = meshkernel::Polygons();
    std::vector<meshkernel::Point> landBoundary{{-1.369282, 21.249086},
                                                {20.885406, 21.539995},
                                                {meshkernel::constants::missing::doubleValue, meshkernel::constants::missing::doubleValue}};

    auto landBoundaries = meshkernel::LandBoundaries(landBoundary, *mesh, polygon);

    // execute flipedges
    meshkernel::FlipEdges flipEdges(*mesh, landBoundaries, true, true);

    auto undoAction = flipEdges.Compute();

    // check the values
    ASSERT_EQ(16, mesh->GetNumEdges());

    // Test the undo action has been computed correctly
    undoAction->Restore();
    // Recompute faces
    mesh->Administrate();

    ASSERT_EQ(originalNodes.size(), mesh->GetNumValidNodes());
    ASSERT_EQ(originalEdges.size(), mesh->GetNumValidEdges());

    meshkernel::UInt count = 0;

    for (meshkernel::UInt i = 0; i < mesh->Nodes().size(); ++i)
    {
        if (mesh->Node(i).IsValid())
        {
            // Check valid nodes
            EXPECT_EQ(originalNodes[count].x, mesh->Node(i).x);
            EXPECT_EQ(originalNodes[count].y, mesh->Node(i).y);
            ++count;
        }
    }

    count = 0;

    for (meshkernel::UInt i = 0; i < mesh->Edges().size(); ++i)
    {
        if (mesh->IsValidEdge(i))
        {
            EXPECT_EQ(originalEdges[count].first, mesh->GetEdge(i).first);
            EXPECT_EQ(originalEdges[count].second, mesh->GetEdge(i).second);
            ++count;
        }
    }
}

TEST(FlipEdges, FlipEdgesMediumTriangularMesh)
{
    // 1 Setup
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/TestOrthogonalizationMediumTriangularGrid_net.nc");

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    // set landboundaries
    auto polygon = meshkernel::Polygons();

    std::vector<meshkernel::Point> landBoundary;
    auto landBoundaries = meshkernel::LandBoundaries(landBoundary, *mesh, polygon);

    // execute flipedges
    meshkernel::FlipEdges flipEdges(*mesh, landBoundaries, true, false);

    auto undoAction = flipEdges.Compute();

    // get the number of edges
    ASSERT_EQ(697, mesh->GetNumEdges());

    // check the values of flipped edges
    ASSERT_EQ(183, mesh->GetEdge(14).first);
    ASSERT_EQ(227, mesh->GetEdge(14).second);

    ASSERT_EQ(58, mesh->GetEdge(33).first);
    ASSERT_EQ(141, mesh->GetEdge(33).second);

    ASSERT_EQ(147, mesh->GetEdge(46).first);
    ASSERT_EQ(145, mesh->GetEdge(46).second);

    ASSERT_EQ(147, mesh->GetEdge(49).first);
    ASSERT_EQ(148, mesh->GetEdge(49).second);

    ASSERT_EQ(242, mesh->GetEdge(68).first);
    ASSERT_EQ(148, mesh->GetEdge(68).second);

    // Test the undo action has been computed correctly
    undoAction->Restore();
    // Recompute faces
    mesh->Administrate();

    ASSERT_EQ(originalNodes.size(), mesh->GetNumValidNodes());
    ASSERT_EQ(originalEdges.size(), mesh->GetNumValidEdges());

    meshkernel::UInt count = 0;

    for (meshkernel::UInt i = 0; i < mesh->Nodes().size(); ++i)
    {
        if (mesh->Node(i).IsValid())
        {
            // Check valid nodes
            EXPECT_EQ(originalNodes[count].x, mesh->Node(i).x);
            EXPECT_EQ(originalNodes[count].y, mesh->Node(i).y);
            ++count;
        }
    }

    count = 0;

    for (meshkernel::UInt i = 0; i < mesh->Edges().size(); ++i)
    {
        if (mesh->IsValidEdge(i))
        {
            EXPECT_EQ(originalEdges[count].first, mesh->GetEdge(i).first);
            EXPECT_EQ(originalEdges[count].second, mesh->GetEdge(i).second);
            ++count;
        }
    }
}
