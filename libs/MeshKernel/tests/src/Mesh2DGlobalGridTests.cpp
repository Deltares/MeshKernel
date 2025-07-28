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

    EXPECT_NEAR(-161.05263157894737, mesh->Node(0).x, tolerance);
    EXPECT_NEAR(-161.05263157894737, mesh->Node(1).x, tolerance);
    EXPECT_NEAR(-161.05263157894737, mesh->Node(2).x, tolerance);
    EXPECT_NEAR(-142.10526315789474, mesh->Node(3).x, tolerance);

    EXPECT_NEAR(0.0000000000000000, mesh->Node(0).y, tolerance);
    EXPECT_NEAR(18.695753703140564, mesh->Node(1).y, tolerance);
    EXPECT_NEAR(-18.695753703140564, mesh->Node(2).y, tolerance);
    EXPECT_NEAR(0.0000000000000000, mesh->Node(3).y, tolerance);
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

    auto mesh2 = ReadLegacyMesh2DFromFile("global_net.nc");

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

    // meshkernel::UInt nodeId = mesh->FindNodeCloseToAPoint({0.0, 0.0}, 1.0e-3);
    // [[maybe_unused]] auto mesh->DeleteNode(nodeId);

    // meshkernel::UInt nodeId = mesh->FindNodeCloseToAPoint({0.0, 0.0}, 1.0e-3);
    // [[maybe_unused]] auto mesh->DeleteNode(nodeId);

    meshkernel::Print(mesh->Nodes(), mesh->Edges());
    return;

    meshkernel::Print(mesh2->Nodes(), mesh2->Edges());
    return;

    const meshkernel::UInt northPoleNodeIndex = 690;
    const meshkernel::UInt southPoleNodeIndex = 691;

    ASSERT_EQ(mesh->GetNumNodes(), 692);
    ASSERT_EQ(mesh->GetNumEdges(), 1333);
    ASSERT_EQ(mesh->m_nodesEdges[northPoleNodeIndex].size(), 5);
    ASSERT_EQ(mesh->m_nodesEdges[southPoleNodeIndex].size(), 5);

    const double tolerance = 1.0e-10;

    ASSERT_NEAR(mesh->Node(northPoleNodeIndex).x, 108.0, tolerance);
    ASSERT_NEAR(mesh->Node(northPoleNodeIndex).y, 90.0, tolerance);

    ASSERT_NEAR(mesh->Node(southPoleNodeIndex).x, 108.0, tolerance);
    ASSERT_NEAR(mesh->Node(southPoleNodeIndex).y, -90.0, tolerance);
}
