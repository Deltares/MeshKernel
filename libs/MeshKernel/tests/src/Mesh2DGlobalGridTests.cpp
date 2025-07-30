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

#include <chrono>
#include <gtest/gtest.h>
#include <random>

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/Mesh2DGenerateGlobal.hpp"
#include "MeshKernel/Polygons.hpp"
#include "MeshKernel/Utilities/Utilities.hpp"
#include "TestUtils/MakeMeshes.hpp"

TEST(GlobalGridTest, Mesh2DGenerateGlobalCompute_ShouldGenerateMesh)
{
    // Execute
    constexpr meshkernel::UInt numLongitudeNodes = 19;
    constexpr meshkernel::UInt numLatitudeNodes = 25;
    const auto mesh = meshkernel::Mesh2DGenerateGlobal::Compute(numLongitudeNodes, numLatitudeNodes, meshkernel::Projection::spherical);

    // Assert
    ASSERT_EQ(1225, mesh->GetNumValidEdges());
    ASSERT_EQ(621, mesh->GetNumValidNodes());

    const double tolerance = 1e-6;

    // Set of nodes on the left (west) side of the mesh

    EXPECT_NEAR(-180, mesh->Node(0).x, tolerance);
    EXPECT_NEAR(-161.052631578947, mesh->Node(1).x, tolerance);
    EXPECT_NEAR(-161.052631578947, mesh->Node(2).x, tolerance);
    EXPECT_NEAR(-180, mesh->Node(3).x, tolerance);
    EXPECT_NEAR(-161.052631578947, mesh->Node(4).x, tolerance);
    EXPECT_NEAR(-180, mesh->Node(5).x, tolerance);

    EXPECT_NEAR(0, mesh->Node(0).y, tolerance);
    EXPECT_NEAR(0, mesh->Node(1).y, tolerance);
    EXPECT_NEAR(18.6957537031406, mesh->Node(2).y, tolerance);
    EXPECT_NEAR(18.6957537031406, mesh->Node(3).y, tolerance);
    EXPECT_NEAR(-18.6957537031406, mesh->Node(4).y, tolerance);
    EXPECT_NEAR(-18.6957537031406, mesh->Node(5).y, tolerance);

    // Set of nodes on the right (east) side of the mesh

    EXPECT_NEAR(142.105263157895, mesh->Node(51).x, tolerance);
    EXPECT_NEAR(142.105263157895, mesh->Node(52).x, tolerance);
    EXPECT_NEAR(142.105263157895, mesh->Node(53).x, tolerance);
    EXPECT_NEAR(161.052631578947, mesh->Node(54).x, tolerance);
    EXPECT_NEAR(161.052631578947, mesh->Node(55).x, tolerance);
    EXPECT_NEAR(161.052631578947, mesh->Node(56).x, tolerance);

    EXPECT_NEAR(0, mesh->Node(51).y, tolerance);
    EXPECT_NEAR(18.6957537031406, mesh->Node(52).y, tolerance);
    EXPECT_NEAR(-18.6957537031406, mesh->Node(53).y, tolerance);
    EXPECT_NEAR(0, mesh->Node(54).y, tolerance);
    EXPECT_NEAR(18.6957537031406, mesh->Node(55).y, tolerance);
    EXPECT_NEAR(-18.6957537031406, mesh->Node(56).y, tolerance);
}

TEST(GlobalGridTest, Mesh2DGenerateGlobalCompute_OnInvalidProjection_ShouldThrowAnException)
{
    // Assert
    const meshkernel::UInt numLongitudeNodes = 19;
    const meshkernel::UInt numLatitudeNodes = 25;
    EXPECT_THROW(meshkernel::Mesh2DGenerateGlobal::Compute(numLongitudeNodes, numLatitudeNodes, meshkernel::Projection::cartesian), meshkernel::MeshKernelError);
}

TEST(GlobalGridTest, GlobalMeshWithPoles)
{
    // Execute

    constexpr meshkernel::UInt numLongitudeNodes = 21;
    constexpr meshkernel::UInt numLatitudeNodes = 20;
    const auto mesh = meshkernel::Mesh2DGenerateGlobal::Compute(numLongitudeNodes, numLatitudeNodes, meshkernel::Projection::spherical);

    std::vector<meshkernel::Point> polyNodes{{-80, -27}, {80, -27}, {80, 27}, {-80, 27}, {-80, -27}};
    [[maybe_unused]] meshkernel::Polygons polygon(polyNodes, meshkernel::Projection::spherical);

    auto inPoly = mesh->IsLocationInPolygon(polygon, meshkernel::Location::Nodes);

    for (size_t i = 0; i < inPoly.size(); ++i)
    {

        if (inPoly[i])
        {
            [[maybe_unused]] auto undoAction = mesh->DeleteNode(i);
        }
    }

    const meshkernel::UInt northPoleNodeIndex = 727;
    const meshkernel::UInt southPoleNodeIndex = 728;

    ASSERT_EQ(mesh->GetNumNodes(), 729);
    ASSERT_EQ(mesh->GetNumEdges(), 1443);
    EXPECT_EQ(mesh->m_nodesEdges[northPoleNodeIndex].size(), 6);
    EXPECT_EQ(mesh->m_nodesEdges[southPoleNodeIndex].size(), 6);

    const double tolerance = 1.0e-10;

    EXPECT_NEAR(mesh->Node(northPoleNodeIndex).x, -111.428571428571, tolerance);
    EXPECT_NEAR(mesh->Node(northPoleNodeIndex).y, 90.0, tolerance);

    EXPECT_NEAR(mesh->Node(southPoleNodeIndex).x, -111.428571428571, tolerance);
    EXPECT_NEAR(mesh->Node(southPoleNodeIndex).y, -90.0, tolerance);

    // Now check that the nodes are connected "around the back" (the west side should be connected to the east side)

    // Pick a random point along the left boundary north of the equator
    const meshkernel::UInt westBoundaryNodeIndexNorth = 106;
    // This east boundary node should be connected to the west boundary node above
    const meshkernel::UInt eastBoundaryNodeIndexNorth = 145;

    EXPECT_NEAR(mesh->Node(westBoundaryNodeIndexNorth).x, -180.0, tolerance);
    EXPECT_NEAR(mesh->Node(westBoundaryNodeIndexNorth).y, 45.8153394766131, tolerance);

    ASSERT_EQ(mesh->m_nodesEdges[westBoundaryNodeIndexNorth].size(), 4);
    ASSERT_EQ(mesh->GetEdge(mesh->m_nodesEdges[westBoundaryNodeIndexNorth][2]).second, eastBoundaryNodeIndexNorth);

    EXPECT_NEAR(mesh->Node(eastBoundaryNodeIndexNorth).x, 162.857142857143, tolerance);
    EXPECT_NEAR(mesh->Node(eastBoundaryNodeIndexNorth).y, 45.8153394766131, tolerance);

    // Pick a random point along the left boundary south of the equator
    const meshkernel::UInt westBoundaryNodeIndexSouth = 108;
    // This east boundary node should be connected to the west boundary node above
    const meshkernel::UInt eastBoundaryNodeIndexSouth = 146;

    EXPECT_NEAR(mesh->Node(westBoundaryNodeIndexSouth).x, -180.0, tolerance);
    EXPECT_NEAR(mesh->Node(westBoundaryNodeIndexSouth).y, -45.8153394766131, tolerance);

    ASSERT_EQ(mesh->m_nodesEdges[westBoundaryNodeIndexSouth].size(), 4);
    ASSERT_EQ(mesh->GetEdge(mesh->m_nodesEdges[westBoundaryNodeIndexSouth][2]).second, eastBoundaryNodeIndexSouth);

    EXPECT_NEAR(mesh->Node(eastBoundaryNodeIndexSouth).x, 162.857142857143, tolerance);
    EXPECT_NEAR(mesh->Node(eastBoundaryNodeIndexSouth).y, -45.8153394766131, tolerance);
}
