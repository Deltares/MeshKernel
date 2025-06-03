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

#include <MeshKernel/Utilities/Utilities.hpp>

TEST(FlipEdges, FlipEdgesWithLandBoundary)
{
    // 1 Setup
    auto mesh = MakeRectangularMeshForTesting(3, 3, 4, meshkernel::Projection::cartesian, {0.0, 0.0});

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    // set landboundaries
    auto polygon = meshkernel::Polygons({{-1.0, 19.0}, {12.0, 19.0}, {12.0, 21.0}, {-1.0, 21.0}, {-1.0, 19.0}}, mesh->m_projection);

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

TEST(FlipEdges, FlipEdgesInPolygonMediumQuadrilateralMeshWithTriangulate)
{
    // 1 Setup
    auto mesh = MakeRectangularMeshForTesting(21, 21, 1000.0, 1000.0, meshkernel::Projection::cartesian);

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    // set landboundaries
    auto landBoundarypolygon = meshkernel::Polygons();

    std::vector<meshkernel::Point> landBoundary{{0.0, 1005.0}, {200.0, 1020.0}, {800.0, 995.0}, {1010.0, 1010.0}};
    auto landBoundaries = meshkernel::LandBoundaries(landBoundary, *mesh, landBoundarypolygon);

    // execute flipedges
    meshkernel::FlipEdges flipEdges(*mesh, landBoundaries, true, true);

    meshkernel::Polygons polygon({{201.0, 201.0}, {700.0, 190.0}, {700.0, 700.0}, {190.0, 700.0}, {201.0, 201.0}, {-998.0, -998.0}, {310.0, 310.0}, {600.0, 310.0}, {450.0, 450.0}, {310.0, 310.0}},
                                 mesh->m_projection);

    auto undoAction = flipEdges.Compute(polygon);

    std::vector<meshkernel::UInt> nodesWithSixEdges{111, 112, 113, 114, 115, 116, 117, 118, 131, 134, 135, 136, 137, 138, 139, 152, 156, 157, 158, 159, 160, 173, 178, 179, 180, 181, 194, 200, 201, 202, 215, 221, 222, 223, 236, 241, 242, 243, 244, 257, 258, 261, 262, 263, 264, 265, 278, 279, 280, 281, 282, 283, 284, 285, 286};

    std::vector<meshkernel::UInt> nodesWithFiveEdges{89, 90, 91, 92, 93, 94, 95, 96, 97, 109, 110, 119, 130, 132, 133, 140, 151, 153, 161, 172, 174, 182, 193, 195, 199, 203, 214, 216, 219, 220, 224, 235, 237, 239, 240, 245, 256, 259, 260, 266, 277, 287, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308};

    for (meshkernel::UInt i = 0; i < nodesWithSixEdges.size(); ++i)
    {
        EXPECT_EQ(mesh->GetNumNodesEdges(nodesWithSixEdges[i]), 6);
    }

    for (meshkernel::UInt i = 0; i < nodesWithFiveEdges.size(); ++i)
    {
        EXPECT_EQ(mesh->GetNumNodesEdges(nodesWithFiveEdges[i]), 5);
    }

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

TEST(FlipEdges, FlipEdgesInPolygonMediumTriangularMesh)
{
    // 1 Setup
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/TestOrthogonalizationMediumTriangularGrid_net.nc");

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    // set landboundaries
    std::vector<meshkernel::Point> landBoundary;
    auto landBoundaries = meshkernel::LandBoundaries(landBoundary, *mesh);

    // execute flipedges
    meshkernel::FlipEdges flipEdges(*mesh, landBoundaries, false, false);

    meshkernel::Polygons polygon({{400.0, 900.0}, {1300.0, 900.0}, {1200.0, 1300.0}, {300.0, 1200.0}, {400.0, 900.0}}, mesh->m_projection);

    auto undoAction = flipEdges.Compute(polygon);

    // get the number of edges
    EXPECT_EQ(697, mesh->GetNumEdges());

    // check the values of several flipped edges
    EXPECT_EQ(242, mesh->GetEdge(80).first);
    EXPECT_EQ(148, mesh->GetEdge(80).second);

    EXPECT_EQ(179, mesh->GetEdge(85).first);
    EXPECT_EQ(242, mesh->GetEdge(85).second);

    EXPECT_EQ(203, mesh->GetEdge(155).first);
    EXPECT_EQ(225, mesh->GetEdge(155).second);

    EXPECT_EQ(225, mesh->GetEdge(277).first);
    EXPECT_EQ(221, mesh->GetEdge(277).second);

    EXPECT_EQ(240, mesh->GetEdge(412).first);
    EXPECT_EQ(236, mesh->GetEdge(412).second);

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
