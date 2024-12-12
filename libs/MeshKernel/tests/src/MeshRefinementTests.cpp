//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2024.
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

#include <iomanip>
#include <iostream>

#include <fstream>

#include <gtest/gtest.h>

#include "MeshKernel/BilinearInterpolationOnGriddedSamples.hpp"
#include "MeshKernel/CasulliDeRefinement.hpp"
#include "MeshKernel/CasulliRefinement.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/MeshRefinement.hpp"
#include "MeshKernel/Operations.hpp"
#include "MeshKernel/Parameters.hpp"
#include "MeshKernel/Polygons.hpp"
#include "MeshKernel/RemoveDisconnectedRegions.hpp"
#include "MeshKernel/SampleAveragingInterpolator.hpp"
#include "MeshKernel/SampleInterpolator.hpp"
#include "MeshKernel/SampleTriangulationInterpolator.hpp"
#include "MeshKernel/SamplesHessianCalculator.hpp"
#include "MeshKernel/SplitRowColumnOfMesh.hpp"
#include "MeshKernel/UndoActions/UndoActionStack.hpp"
#include "MeshKernel/Utilities/Utilities.hpp"

#include "TestUtils/Definitions.hpp"
#include "TestUtils/MakeCurvilinearGrids.hpp"
#include "TestUtils/MakeMeshes.hpp"
#include "TestUtils/MeshReaders.hpp"
#include "TestUtils/SampleFileReader.hpp"
#include "TestUtils/SampleGenerator.hpp"

using namespace meshkernel;

namespace mk = meshkernel;

TEST(MeshRefinement, MeshRefinementRefinementLevels_OnFourByFourWithFourSamples_ShouldRefinemesh)
{
    auto mesh = MakeRectangularMeshForTesting(5, 5, 10.0, Projection::cartesian);

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    // sample points
    std::vector<Sample> samples{
        {14.7153645, 14.5698833, 1.0},
        {24.7033062, 14.4729137, 1.0},
        {15.5396099, 24.2669525, 1.0},
        {23.8305721, 23.9275551, 1.0}};

    auto interpolator = std::make_unique<AveragingInterpolation>(*mesh,
                                                                 samples,
                                                                 AveragingInterpolation::Method::MinAbsValue,
                                                                 Location::Faces,
                                                                 1.0,
                                                                 false,
                                                                 false,
                                                                 1);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 1;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;
    meshRefinementParameters.min_edge_size = 1.0;
    meshRefinementParameters.account_for_samples_outside = 0;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.refinement_type = 2;
    meshRefinementParameters.smoothing_iterations = 0;

    MeshRefinement meshRefinement(*mesh,
                                  std::move(interpolator),
                                  meshRefinementParameters);

    auto undoAction = meshRefinement.Compute();

    // 3 Validation edges connecting hanging nodes

    // bottom side
    ASSERT_EQ(5, mesh->GetEdge(73).first);
    ASSERT_EQ(25, mesh->GetEdge(73).second);

    ASSERT_EQ(10, mesh->GetEdge(72).first);
    ASSERT_EQ(25, mesh->GetEdge(72).second);

    ASSERT_EQ(10, mesh->GetEdge(77).first);
    ASSERT_EQ(28, mesh->GetEdge(77).second);

    ASSERT_EQ(15, mesh->GetEdge(76).first);
    ASSERT_EQ(28, mesh->GetEdge(76).second);

    // right side
    ASSERT_EQ(21, mesh->GetEdge(81).first);
    ASSERT_EQ(35, mesh->GetEdge(81).second);

    ASSERT_EQ(22, mesh->GetEdge(80).first);
    ASSERT_EQ(35, mesh->GetEdge(80).second);

    ASSERT_EQ(22, mesh->GetEdge(83).first);
    ASSERT_EQ(36, mesh->GetEdge(83).second);

    ASSERT_EQ(23, mesh->GetEdge(82).first);
    ASSERT_EQ(36, mesh->GetEdge(82).second);

    // upper side
    ASSERT_EQ(19, mesh->GetEdge(79).first);
    ASSERT_EQ(30, mesh->GetEdge(79).second);

    ASSERT_EQ(14, mesh->GetEdge(78).first);
    ASSERT_EQ(30, mesh->GetEdge(78).second);

    ASSERT_EQ(14, mesh->GetEdge(75).first);
    ASSERT_EQ(27, mesh->GetEdge(75).second);

    ASSERT_EQ(9, mesh->GetEdge(74).first);
    ASSERT_EQ(27, mesh->GetEdge(74).second);

    // left side
    ASSERT_EQ(3, mesh->GetEdge(71).first);
    ASSERT_EQ(32, mesh->GetEdge(71).second);

    ASSERT_EQ(2, mesh->GetEdge(70).first);
    ASSERT_EQ(32, mesh->GetEdge(70).second);

    ASSERT_EQ(2, mesh->GetEdge(69).first);
    ASSERT_EQ(31, mesh->GetEdge(69).second);

    ASSERT_EQ(1, mesh->GetEdge(68).first);
    ASSERT_EQ(31, mesh->GetEdge(68).second);

    // total number of edges
    ASSERT_EQ(84, mesh->GetNumEdges());

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

TEST(MeshRefinement, RefinementOnAFourByFourMeshWithSamplesShouldRefine)
{
    auto mesh = MakeRectangularMeshForTesting(4, 4, 500.0, Projection::cartesian);

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    // sample points
    std::vector<Sample> samples{
        {14.7153645, 14.5698833, 0.1},
        {24.7033062, 14.4729137, 0.1},
        {15.5396099, 24.2669525, 0.1},
        {23.8305721, 23.9275551, 0.1}};

    auto interpolator = std::make_unique<AveragingInterpolation>(*mesh,
                                                                 samples,
                                                                 AveragingInterpolation::Method::MinAbsValue,
                                                                 Location::Faces,
                                                                 1.0,
                                                                 false,
                                                                 false,
                                                                 1);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 4;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;
    meshRefinementParameters.min_edge_size = 2.0;
    meshRefinementParameters.account_for_samples_outside = 0;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.refinement_type = 1;
    meshRefinementParameters.smoothing_iterations = 0;

    MeshRefinement meshRefinement(*mesh,
                                  std::move(interpolator),
                                  meshRefinementParameters);
    auto undoAction = meshRefinement.Compute();

    // Assert number of edges and nodes
    ASSERT_EQ(60, mesh->GetNumEdges());
    ASSERT_EQ(31, mesh->GetNumNodes());

    // Assert edges
    ASSERT_EQ(0, mesh->GetEdge(0).first);
    ASSERT_EQ(26, mesh->GetEdge(0).second);

    ASSERT_EQ(1, mesh->GetEdge(1).first);
    ASSERT_EQ(17, mesh->GetEdge(1).second);

    ASSERT_EQ(2, mesh->GetEdge(2).first);
    ASSERT_EQ(6, mesh->GetEdge(2).second);

    ASSERT_EQ(3, mesh->GetEdge(3).first);
    ASSERT_EQ(7, mesh->GetEdge(3).second);

    ASSERT_EQ(4, mesh->GetEdge(4).first);
    ASSERT_EQ(8, mesh->GetEdge(4).second);

    ASSERT_EQ(5, mesh->GetEdge(5).first);
    ASSERT_EQ(9, mesh->GetEdge(5).second);

    ASSERT_EQ(6, mesh->GetEdge(6).first);
    ASSERT_EQ(10, mesh->GetEdge(6).second);

    ASSERT_EQ(7, mesh->GetEdge(7).first);
    ASSERT_EQ(11, mesh->GetEdge(7).second);

    ASSERT_EQ(8, mesh->GetEdge(8).first);
    ASSERT_EQ(12, mesh->GetEdge(8).second);

    ASSERT_EQ(9, mesh->GetEdge(9).first);
    ASSERT_EQ(13, mesh->GetEdge(9).second);

    ASSERT_EQ(10, mesh->GetEdge(10).first);
    ASSERT_EQ(14, mesh->GetEdge(10).second);

    ASSERT_EQ(11, mesh->GetEdge(11).first);
    ASSERT_EQ(15, mesh->GetEdge(11).second);

    ASSERT_EQ(1, mesh->GetEdge(12).first);
    ASSERT_EQ(18, mesh->GetEdge(12).second);

    ASSERT_EQ(2, mesh->GetEdge(13).first);
    ASSERT_EQ(1, mesh->GetEdge(13).second);

    ASSERT_EQ(3, mesh->GetEdge(14).first);
    ASSERT_EQ(2, mesh->GetEdge(14).second);

    ASSERT_EQ(5, mesh->GetEdge(15).first);
    ASSERT_EQ(19, mesh->GetEdge(15).second);

    ASSERT_EQ(6, mesh->GetEdge(16).first);
    ASSERT_EQ(5, mesh->GetEdge(16).second);

    ASSERT_EQ(7, mesh->GetEdge(17).first);
    ASSERT_EQ(6, mesh->GetEdge(17).second);

    ASSERT_EQ(9, mesh->GetEdge(18).first);
    ASSERT_EQ(8, mesh->GetEdge(18).second);

    ASSERT_EQ(10, mesh->GetEdge(19).first);
    ASSERT_EQ(9, mesh->GetEdge(19).second);

    ASSERT_EQ(11, mesh->GetEdge(20).first);
    ASSERT_EQ(10, mesh->GetEdge(20).second);

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

TEST(MeshRefinement, MeshRefinementRefinementLevels_SmallTriangualMeshTwoSamples_ShouldRefinemesh)
{
    // Prepare
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/SmallTriangularGrid_net.nc");

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    // sample points
    std::vector<Sample> samples{
        {359.8657532, 350.3144836, 1.0},
        {387.5152588, 299.2614746, 1.0}};

    auto interpolator = std::make_unique<AveragingInterpolation>(*mesh,
                                                                 samples,
                                                                 AveragingInterpolation::Method::MinAbsValue,
                                                                 Location::Faces,
                                                                 1.0,
                                                                 false,
                                                                 false,
                                                                 1);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 1;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;
    meshRefinementParameters.min_edge_size = 50.0;
    meshRefinementParameters.account_for_samples_outside = 0;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.refinement_type = 2;
    meshRefinementParameters.smoothing_iterations = 0;

    MeshRefinement meshRefinement(*mesh,
                                  std::move(interpolator),
                                  meshRefinementParameters);

    auto undoAction = meshRefinement.Compute();

    // edges connecting hanging nodes
    ASSERT_EQ(10, mesh->GetEdge(32).first);
    ASSERT_EQ(2, mesh->GetEdge(32).second);

    ASSERT_EQ(14, mesh->GetEdge(33).first);
    ASSERT_EQ(4, mesh->GetEdge(33).second);

    ASSERT_EQ(13, mesh->GetEdge(34).first);
    ASSERT_EQ(5, mesh->GetEdge(34).second);

    ASSERT_EQ(11, mesh->GetEdge(31).first);
    ASSERT_EQ(1, mesh->GetEdge(31).second);

    // total number of edges
    ASSERT_EQ(35, mesh->GetNumEdges());

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

TEST(MeshRefinement, RefineBasedOnPolygonTriangularMesh)
{
    // Prepare
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/SmallTriangularGrid_net.nc");

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    // Polygon sample
    std::vector<Point> point{
        {399.638169557229, 504.294564030922},
        {361.827403800769, 129.967983041964},
        {651.709941266965, 113.583317880831},
        {666.834247569549, 411.028008498319},
        {410.981399284167, 505.55492288947},
        {399.638169557229, 504.294564030922}};

    auto polygon = Polygons(point, mesh->m_projection);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 1;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;

    MeshRefinement meshRefinement(*mesh, polygon, meshRefinementParameters);
    auto undoAction = meshRefinement.Compute();

    // total number of edges
    ASSERT_EQ(15, mesh->GetNumNodes());
    ASSERT_EQ(33, mesh->GetNumEdges());

    // assert on newly generated edges
    ASSERT_EQ(10, mesh->GetEdge(20).first);
    ASSERT_EQ(11, mesh->GetEdge(20).second);

    ASSERT_EQ(11, mesh->GetEdge(21).first);
    ASSERT_EQ(12, mesh->GetEdge(21).second);

    ASSERT_EQ(12, mesh->GetEdge(22).first);
    ASSERT_EQ(10, mesh->GetEdge(22).second);

    ASSERT_EQ(14, mesh->GetEdge(23).first);
    ASSERT_EQ(13, mesh->GetEdge(23).second);

    ASSERT_EQ(13, mesh->GetEdge(24).first);
    ASSERT_EQ(11, mesh->GetEdge(24).second);

    ASSERT_EQ(11, mesh->GetEdge(25).first);
    ASSERT_EQ(14, mesh->GetEdge(25).second);

    ASSERT_EQ(10, mesh->GetEdge(26).first);
    ASSERT_EQ(4, mesh->GetEdge(26).second);

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

TEST(MeshRefinement, ThreeBythreeWithThreeSamplesPerFace)
{
    // Prepare
    auto mesh = MakeRectangularMeshForTesting(4, 4, 10.0, Projection::cartesian);

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    // sample points
    std::vector<Sample> samples{
        {2.7091951, 5.4000854, 0.0000000},
        {6.4910383, 2.4182367, 0.0000000},
        {8.0910482, 6.7091894, 0.0000000},
        {13.2910795, 5.2182646, 0.0000000},
        {16.1274605, 2.0909605, 0.0000000},
        {18.7820244, 7.5091972, 0.0000000},
        {23.5456886, 8.1637497, 0.0000000},
        {24.6366081, 1.5818644, 0.0000000},
        {27.8729897, 6.5273695, 0.0000000},
        {28.0184441, 14.7092705, 0.0000000},
        {23.8366013, 12.4910660, 0.0000000},
        {22.6002312, 17.2183857, 0.0000000},
        {24.0184212, 23.7639065, 0.0000000},
        {27.9457169, 25.9093838, 0.0000000},
        {24.6366081, 27.3275795, 0.0000000},
        {17.4365616, 27.5821285, 0.0000000},
        {16.6001930, 22.8184433, 0.0000000},
        {12.7092581, 27.2548523, 0.0000000},
        {4.7455721, 27.6912193, 0.0000000},
        {2.6728315, 24.7457352, 0.0000000},
        {7.5819540, 22.9638996, 0.0000000},
        {3.5455656, 15.1820030, 0.0000000},
        {4.8546629, 11.8365135, 0.0000000},
        {8.7455969, 17.2183857, 0.0000000},
        {11.8404741, 17.6817989, 3.0000000},
        {13.5837603, 12.1783361, 3.0000000},
        {17.2156067, 16.9106121, 3.0000000}};

    auto interpolator = std::make_unique<AveragingInterpolation>(*mesh,
                                                                 samples,
                                                                 AveragingInterpolation::Method::MinAbsValue,
                                                                 Location::Faces,
                                                                 1.0,
                                                                 false,
                                                                 false,
                                                                 1);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 2;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;
    meshRefinementParameters.min_edge_size = 3.0;
    meshRefinementParameters.account_for_samples_outside = 0;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.refinement_type = 1;

    MeshRefinement meshRefinement(*mesh, std::move(interpolator), meshRefinementParameters);

    auto undoAction = meshRefinement.Compute();

    // assert on number of nodes and edges
    ASSERT_EQ(49, mesh->GetNumNodes());
    ASSERT_EQ(84, mesh->GetNumEdges());

    // assert on newly generated edges
    ASSERT_EQ(26, mesh->GetEdge(70).first);
    ASSERT_EQ(14, mesh->GetEdge(70).second);

    ASSERT_EQ(27, mesh->GetEdge(71).first);
    ASSERT_EQ(15, mesh->GetEdge(71).second);

    ASSERT_EQ(28, mesh->GetEdge(72).first);
    ASSERT_EQ(0, mesh->GetEdge(72).second);

    ASSERT_EQ(29, mesh->GetEdge(73).first);
    ASSERT_EQ(1, mesh->GetEdge(73).second);

    ASSERT_EQ(30, mesh->GetEdge(74).first);
    ASSERT_EQ(2, mesh->GetEdge(74).second);

    ASSERT_EQ(31, mesh->GetEdge(75).first);
    ASSERT_EQ(4, mesh->GetEdge(75).second);

    ASSERT_EQ(32, mesh->GetEdge(76).first);
    ASSERT_EQ(5, mesh->GetEdge(76).second);

    ASSERT_EQ(33, mesh->GetEdge(77).first);
    ASSERT_EQ(6, mesh->GetEdge(77).second);

    ASSERT_EQ(34, mesh->GetEdge(78).first);
    ASSERT_EQ(8, mesh->GetEdge(78).second);

    ASSERT_EQ(35, mesh->GetEdge(79).first);
    ASSERT_EQ(9, mesh->GetEdge(79).second);

    ASSERT_EQ(36, mesh->GetEdge(80).first);
    ASSERT_EQ(10, mesh->GetEdge(80).second);

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

TEST(MeshRefinement, WindowOfRefinementFile)
{
    // Prepare
    auto mesh = MakeRectangularMeshForTesting(4, 4, 40.0, Projection::cartesian, {197253.0, 442281.0});

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    // Sample points
    std::vector<Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/MeshRefinementTests/WindowOfRefinementFile.xyz");

    auto interpolator = std::make_unique<AveragingInterpolation>(*mesh,
                                                                 samples,
                                                                 AveragingInterpolation::Method::MinAbsValue,
                                                                 Location::Faces,
                                                                 1.0,
                                                                 false,
                                                                 false,
                                                                 1);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 4;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;
    meshRefinementParameters.min_edge_size = 3.0;
    meshRefinementParameters.account_for_samples_outside = 0;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.refinement_type = 1;
    meshRefinementParameters.smoothing_iterations = 0;

    MeshRefinement meshRefinement(*mesh,
                                  std::move(interpolator),
                                  meshRefinementParameters);

    auto undoAction = meshRefinement.Compute();

    // total number of edges
    ASSERT_EQ(461, mesh->GetNumNodes());
    ASSERT_EQ(918, mesh->GetNumEdges());

    ASSERT_EQ(233, mesh->GetEdge(906).first);
    ASSERT_EQ(206, mesh->GetEdge(906).second);

    ASSERT_EQ(206, mesh->GetEdge(907).first);
    ASSERT_EQ(70, mesh->GetEdge(907).second);

    ASSERT_EQ(70, mesh->GetEdge(908).first);
    ASSERT_EQ(233, mesh->GetEdge(908).second);

    ASSERT_EQ(179, mesh->GetEdge(909).first);
    ASSERT_EQ(235, mesh->GetEdge(909).second);

    ASSERT_EQ(235, mesh->GetEdge(910).first);
    ASSERT_EQ(62, mesh->GetEdge(910).second);

    ASSERT_EQ(62, mesh->GetEdge(911).first);
    ASSERT_EQ(179, mesh->GetEdge(911).second);

    ASSERT_EQ(249, mesh->GetEdge(912).first);
    ASSERT_EQ(320, mesh->GetEdge(912).second);

    ASSERT_EQ(320, mesh->GetEdge(913).first);
    ASSERT_EQ(84, mesh->GetEdge(913).second);

    ASSERT_EQ(84, mesh->GetEdge(914).first);
    ASSERT_EQ(249, mesh->GetEdge(914).second);

    ASSERT_EQ(327, mesh->GetEdge(915).first);
    ASSERT_EQ(326, mesh->GetEdge(915).second);

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

TEST(MeshRefinement, MeshRefinementRefinementLevels_OnWindowOfRefinementFile_ShouldRefinemesh)
{
    // Prepare
    auto mesh = MakeRectangularMeshForTesting(4, 4, 40.0, Projection::cartesian, {197253.0, 442281.0});

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    // Sample points
    std::vector<Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/MeshRefinementTests/WindowOfRefinementFile.xyz");

    auto interpolator = std::make_unique<AveragingInterpolation>(*mesh,
                                                                 samples,
                                                                 AveragingInterpolation::Method::Max,
                                                                 Location::Faces,
                                                                 1.01,
                                                                 false,
                                                                 true,
                                                                 1);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 10;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;
    meshRefinementParameters.min_edge_size = 0.5;
    meshRefinementParameters.account_for_samples_outside = 0;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.refinement_type = 2;
    meshRefinementParameters.smoothing_iterations = 0;

    MeshRefinement meshRefinement(*mesh,
                                  std::move(interpolator),
                                  meshRefinementParameters);

    auto undoAction = meshRefinement.Compute();

    // total number of edges
    ASSERT_EQ(413, mesh->GetNumNodes());
    ASSERT_EQ(839, mesh->GetNumEdges());

    ASSERT_EQ(0, mesh->GetEdge(0).first);
    ASSERT_EQ(129, mesh->GetEdge(0).second);

    ASSERT_EQ(1, mesh->GetEdge(1).first);
    ASSERT_EQ(130, mesh->GetEdge(1).second);

    ASSERT_EQ(2, mesh->GetEdge(2).first);
    ASSERT_EQ(131, mesh->GetEdge(2).second);

    ASSERT_EQ(3, mesh->GetEdge(3).first);
    ASSERT_EQ(48, mesh->GetEdge(3).second);

    ASSERT_EQ(4, mesh->GetEdge(4).first);
    ASSERT_EQ(49, mesh->GetEdge(4).second);

    ASSERT_EQ(5, mesh->GetEdge(5).first);
    ASSERT_EQ(132, mesh->GetEdge(5).second);

    ASSERT_EQ(6, mesh->GetEdge(6).first);
    ASSERT_EQ(133, mesh->GetEdge(6).second);

    ASSERT_EQ(7, mesh->GetEdge(7).first);
    ASSERT_EQ(22, mesh->GetEdge(7).second);

    ASSERT_EQ(8, mesh->GetEdge(8).first);
    ASSERT_EQ(23, mesh->GetEdge(8).second);

    ASSERT_EQ(9, mesh->GetEdge(9).first);
    ASSERT_EQ(134, mesh->GetEdge(9).second);

    ASSERT_EQ(10, mesh->GetEdge(10).first);
    ASSERT_EQ(135, mesh->GetEdge(10).second);

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

TEST(MeshRefinement, RefineBasedOnPolygon)
{
    // Prepare
    auto mesh = MakeRectangularMeshForTesting(5, 5, 10.0, Projection::cartesian);

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    std::vector<Point> point{
        {25.0, -10.0},
        {25.0, 15.0},
        {45.0, 15.0},
        {45.0, -10.0},
        {25.0, -10.0}};

    const auto polygon = Polygons(point, mesh->m_projection);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 1;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;

    MeshRefinement meshRefinement(*mesh, polygon, meshRefinementParameters);

    auto undoAction = meshRefinement.Compute();

    // total number of edges
    ASSERT_EQ(30, mesh->GetNumNodes());
    ASSERT_EQ(52, mesh->GetNumEdges());

    ASSERT_EQ(25, mesh->GetEdge(40).first);
    ASSERT_EQ(29, mesh->GetEdge(40).second);

    ASSERT_EQ(28, mesh->GetEdge(41).first);
    ASSERT_EQ(29, mesh->GetEdge(41).second);

    ASSERT_EQ(26, mesh->GetEdge(42).first);
    ASSERT_EQ(29, mesh->GetEdge(42).second);

    ASSERT_EQ(27, mesh->GetEdge(43).first);
    ASSERT_EQ(29, mesh->GetEdge(43).second);

    ASSERT_EQ(25, mesh->GetEdge(44).first);
    ASSERT_EQ(20, mesh->GetEdge(44).second);

    ASSERT_EQ(26, mesh->GetEdge(45).first);
    ASSERT_EQ(21, mesh->GetEdge(45).second);

    ASSERT_EQ(27, mesh->GetEdge(46).first);
    ASSERT_EQ(15, mesh->GetEdge(46).second);

    ASSERT_EQ(28, mesh->GetEdge(47).first);
    ASSERT_EQ(20, mesh->GetEdge(47).second);

    ASSERT_EQ(10, mesh->GetEdge(48).first);
    ASSERT_EQ(27, mesh->GetEdge(48).second);

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

TEST(MeshRefinement, RefineBasedOnPolygonThreeByThree)
{
    // Prepare
    auto mesh = MakeRectangularMeshForTesting(4, 4, 10.0, Projection::cartesian);

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    std::vector<Point> point{
        {9.09836065573771, 34.016393442623},
        {7.18032786885247, 7.75409836065574},
        {34.6229508196721, 6.5},
        {34.4194409808304, 26.6983515050386},
        {34.327868852459, 35.7868852459016},
        {29.0521194370216, 35.4840476661635},
        {9.90983606557378, 34.3852459016394},
        {9.09836065573771, 34.016393442623}};

    const auto polygon = Polygons(point, mesh->m_projection);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 2;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;

    MeshRefinement meshRefinement(*mesh, polygon, meshRefinementParameters);
    auto undoAction = meshRefinement.Compute();

    // assert on number of nodes and edges
    ASSERT_EQ(48, mesh->GetNumNodes());
    ASSERT_EQ(96, mesh->GetNumEdges());
    ASSERT_EQ(49, mesh->GetNumFaces());

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

TEST(MeshRefinement, FourByFourWithFourSamplesSpherical)
{

    auto mesh = MakeRectangularMeshForTesting(4, 4, 0.0033, Projection::spherical, {41.1, 41.1});

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    // sample points
    std::vector<Sample> samples{
        {41.1050110, 41.1049728, 1.0},
        {41.1084785, 41.1048775, 1.0},
        {41.1085625, 41.1083946, 1.0},
        {41.1052971, 41.1083336, 1.0}};

    auto interpolator = std::make_unique<AveragingInterpolation>(*mesh,
                                                                 samples,
                                                                 AveragingInterpolation::Method::MinAbsValue,
                                                                 Location::Faces,
                                                                 1.0,
                                                                 false,
                                                                 false,
                                                                 1);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 1;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;
    meshRefinementParameters.min_edge_size = 0.00165;
    meshRefinementParameters.account_for_samples_outside = 0;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.refinement_type = 2;
    meshRefinementParameters.smoothing_iterations = 0;

    MeshRefinement meshRefinement(*mesh,
                                  std::move(interpolator),
                                  meshRefinementParameters);
    auto undoAction = meshRefinement.Compute();

    ASSERT_EQ(60, mesh->GetNumEdges());
    ASSERT_EQ(32, mesh->GetNumNodes());

    // sides of the refined part
    ASSERT_EQ(5, mesh->GetEdge(5).first);
    ASSERT_EQ(16, mesh->GetEdge(5).second);

    ASSERT_EQ(16, mesh->GetEdge(40).first);
    ASSERT_EQ(9, mesh->GetEdge(40).second);

    ASSERT_EQ(9, mesh->GetEdge(9).first);
    ASSERT_EQ(19, mesh->GetEdge(9).second);

    ASSERT_EQ(19, mesh->GetEdge(43).first);
    ASSERT_EQ(13, mesh->GetEdge(43).second);

    ASSERT_EQ(6, mesh->GetEdge(16).first);
    ASSERT_EQ(22, mesh->GetEdge(16).second);

    ASSERT_EQ(22, mesh->GetEdge(46).first);
    ASSERT_EQ(5, mesh->GetEdge(46).second);

    ASSERT_EQ(7, mesh->GetEdge(17).first);
    ASSERT_EQ(23, mesh->GetEdge(17).second);

    ASSERT_EQ(23, mesh->GetEdge(47).first);
    ASSERT_EQ(6, mesh->GetEdge(47).second);

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

TEST(MeshRefinement, RefinementFileBasedOnLevels_OnSphericalMesh_ShouldRefine)
{
    // Prepare
    auto mesh = MakeRectangularMeshForTesting(6, 6, 0.0033, Projection::spherical, {41.1, 41.1});

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 1;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;
    meshRefinementParameters.min_edge_size = 0.00165;
    meshRefinementParameters.account_for_samples_outside = 0;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.refinement_type = 2;

    // Execute
    const auto polygon = Polygons();
    MeshRefinement meshRefinement(*mesh, polygon, meshRefinementParameters);
    auto undoAction = meshRefinement.Compute();

    // Assert, we passed from 36 to 49 nodes
    ASSERT_EQ(121, mesh->GetNumNodes());

    constexpr double tolerance = 1e-6;
    ASSERT_NEAR(41.100000000000001, mesh->Node(0).x, tolerance);
    ASSERT_NEAR(41.103300000000004, mesh->Node(6).x, tolerance);
    ASSERT_NEAR(41.106600000000000, mesh->Node(12).x, tolerance);
    ASSERT_NEAR(41.109900000000003, mesh->Node(18).x, tolerance);
    ASSERT_NEAR(41.113199999999999, mesh->Node(24).x, tolerance);

    ASSERT_NEAR(41.100000000000001, mesh->Node(0).y, tolerance);
    ASSERT_NEAR(41.100000000000001, mesh->Node(6).y, tolerance);
    ASSERT_NEAR(41.100000000000001, mesh->Node(12).y, tolerance);
    ASSERT_NEAR(41.100000000000001, mesh->Node(18).y, tolerance);
    ASSERT_NEAR(41.100000000000001, mesh->Node(24).y, tolerance);

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

TEST(MeshRefinement, RefineCurvilinearGrid)
{
    auto mesh = MakeCurvilinearGridForTesting();

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 1;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;

    const auto polygon = Polygons();
    MeshRefinement meshRefinement(*mesh, polygon, meshRefinementParameters);
    auto undoAction = meshRefinement.Compute();

    mesh->ComputeEdgesLengths();

    // if the circumcenters are wrongly computed, some edges will be smaller than half cell size
    for (meshkernel::UInt i = 0; i < mesh->GetNumEdges(); ++i)
    {
        ASSERT_GT(mesh->m_edgeLengths[i], 0.4);
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

TEST(MeshRefinement, RefineElongatedFaces)
{
    // Prepare
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/MeshRefinementTests/CurvilinearEnlonged.nc");

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    std::vector<Point> point{
        {2018.73356016594, 1165.26385465465},
        {1694.83823023783, 1131.75744121381},
        {1708.79923583818, 640.330044081495},
        {2363.7176353495, 645.640193266722},
        {2741.91365026406, 648.706647441705},
        {2722.36824242357, 1173.64045801486},
        {2038.27896800643, 1165.26385465465},
        {2018.73356016594, 1165.26385465465}};

    const auto polygon = Polygons(point, mesh->m_projection);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 2;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;

    MeshRefinement meshRefinement(*mesh, polygon, meshRefinementParameters);

    // Execute
    auto undoAction = meshRefinement.Compute();

    // Assert circumcenters are correctly computed
    constexpr double tolerance = 1e-6;

    // Compare the x-location of the circumcentre.
    ASSERT_NEAR(1673.0860169014584, mesh->m_facesCircumcenters[0].x, tolerance);
    ASSERT_NEAR(1660.6851354957175, mesh->m_facesCircumcenters[1].x, tolerance);
    ASSERT_NEAR(1660.5667704694627, mesh->m_facesCircumcenters[2].x, tolerance);
    ASSERT_NEAR(1672.0912775041329, mesh->m_facesCircumcenters[3].x, tolerance);
    ASSERT_NEAR(1659.9354211078053, mesh->m_facesCircumcenters[4].x, tolerance);
    ASSERT_NEAR(1659.8248648846848, mesh->m_facesCircumcenters[5].x, tolerance);
    ASSERT_NEAR(1671.1074693287451, mesh->m_facesCircumcenters[6].x, tolerance);
    ASSERT_NEAR(1659.1859707978906, mesh->m_facesCircumcenters[7].x, tolerance);
    ASSERT_NEAR(1659.0828479935451, mesh->m_facesCircumcenters[8].x, tolerance);
    ASSERT_NEAR(1670.1135380487042, mesh->m_facesCircumcenters[9].x, tolerance);
    ASSERT_NEAR(1658.4375379474418, mesh->m_facesCircumcenters[10].x, tolerance);

    ASSERT_NEAR(645.10565853980427, mesh->m_facesCircumcenters[0].y, tolerance);
    ASSERT_NEAR(646.20173898461292, mesh->m_facesCircumcenters[1].y, tolerance);
    ASSERT_NEAR(654.66978646149676, mesh->m_facesCircumcenters[2].y, tolerance);
    ASSERT_NEAR(662.96936808193914, mesh->m_facesCircumcenters[3].y, tolerance);
    ASSERT_NEAR(664.17198930960728, mesh->m_facesCircumcenters[4].y, tolerance);
    ASSERT_NEAR(672.49588595953537, mesh->m_facesCircumcenters[5].y, tolerance);
    ASSERT_NEAR(680.82566423059006, mesh->m_facesCircumcenters[6].y, tolerance);
    ASSERT_NEAR(682.12720924968505, mesh->m_facesCircumcenters[7].y, tolerance);
    ASSERT_NEAR(690.31918193748425, mesh->m_facesCircumcenters[8].y, tolerance);
    ASSERT_NEAR(698.66471917887850, mesh->m_facesCircumcenters[9].y, tolerance);
    ASSERT_NEAR(700.06356972686194, mesh->m_facesCircumcenters[10].y, tolerance);

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

TEST(MeshRefinement, BilinearInterpolationWithGriddedSamplesOnLandShouldNotRefine)
{
    // Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 10.0, Projection::cartesian);

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    std::vector<float> values{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    Point origin{-5.0, -5.0};
    auto interpolator = std::make_unique<BilinearInterpolationOnGriddedSamples<float>>(*mesh, 2, 2, origin, 10.0, values);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 1;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;
    meshRefinementParameters.min_edge_size = 5.0;
    meshRefinementParameters.account_for_samples_outside = 1;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.refinement_type = 1;

    MeshRefinement meshRefinement(*mesh, std::move(interpolator), meshRefinementParameters, true);

    // Execute
    auto undoAction = meshRefinement.Compute();

    // Assert: all bathy values are positive and we are in land, so nothing gets refined
    ASSERT_EQ(4, mesh->GetNumEdges());

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

TEST(MeshRefinement, BilinearInterpolationWithGriddedSamplesOnLandAndSeaShouldRefine)
{
    // Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 10.0, Projection::cartesian);

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    std::vector<float> values{-1.0, -2.0, 3.0, -4.0, -5.0, 6.0, 7.0, 8.0, 9.0};
    Point origin{-5.0, -5.0};
    auto interpolator = std::make_unique<BilinearInterpolationOnGriddedSamples<float>>(*mesh, 3, 3, origin, 10.0, values);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 1;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;
    meshRefinementParameters.min_edge_size = 5.0;
    meshRefinementParameters.account_for_samples_outside = 1;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.refinement_type = 1;

    MeshRefinement meshRefinement(*mesh, std::move(interpolator), meshRefinementParameters, true);

    // Execute
    auto undoAction = meshRefinement.Compute();

    // Assert: all depth values are positive and we are in land, so nothing gets refined
    ASSERT_EQ(12, mesh->GetNumEdges());

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

TEST(MeshRefinement, BilinearInterpolationWithAllGriddedSamplesOnSeaShouldRefine)
{
    // Setup
    auto mesh = MakeRectangularMeshForTesting(2, 2, 10.0, Projection::cartesian);

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    std::vector<float> values{-1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0, -8.0, -9.0};
    Point origin{-5.0, -5.0};
    auto interpolator = std::make_unique<BilinearInterpolationOnGriddedSamples<float>>(*mesh, 2, 2, origin, 10.0, values);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 1;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;
    meshRefinementParameters.min_edge_size = 5.0;
    meshRefinementParameters.account_for_samples_outside = 1;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.refinement_type = 1;

    MeshRefinement meshRefinement(*mesh,
                                  std::move(interpolator),
                                  meshRefinementParameters, true);

    // Execute
    auto undoAction = meshRefinement.Compute();

    // Assert: nothing gets refined
    ASSERT_EQ(4, mesh->GetNumEdges());

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

class RidgeRefinementTestCases : public testing::TestWithParam<std::tuple<FunctionTestCase, UInt, UInt>>
{
public:
    [[nodiscard]] static std::vector<std::tuple<FunctionTestCase, UInt, UInt>> GetData()
    {
        return std::vector{
            std::make_tuple<FunctionTestCase, UInt, UInt>(FunctionTestCase::GaussianBump, 1165, 2344),
            std::make_tuple<FunctionTestCase, UInt, UInt>(FunctionTestCase::GaussianWave, 5297, 10784),
            std::make_tuple<FunctionTestCase, UInt, UInt>(FunctionTestCase::RidgeXDirection, 2618, 5694),
            std::make_tuple<FunctionTestCase, UInt, UInt>(FunctionTestCase::ArctanFunction, 2309, 5028)};
    }
};

TEST_P(RidgeRefinementTestCases, expectedResults)
{
    // Prepare
    const auto [testCase, numNodes, numEdges] = GetParam();

    UInt nx = 41;
    UInt ny = 21;

    double deltaX = 10.0;
    double deltaY = 10.0;

    double dimX = (nx - 1) * deltaX;
    double dimY = (ny - 1) * deltaY;

    auto mesh = MakeRectangularMeshForTesting(nx, ny, dimX, dimY, Projection::cartesian);

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    UInt superSample = 2;
    UInt sampleNx = (nx - 1) * superSample + 1;
    UInt sampleNy = (ny - 1) * superSample + 1;

    const double sampleDeltaX = deltaX / static_cast<double>(superSample);
    const double sampleDeltaY = deltaY / static_cast<double>(superSample);

    const auto sampleData = generateSampleData(testCase, sampleNx, sampleNy, sampleDeltaX, sampleDeltaY);

    auto samplesHessian = SamplesHessianCalculator::ComputeSamplesHessian(sampleData, mesh->m_projection, 0, sampleNx, sampleNy);

    auto interpolator = std::make_unique<AveragingInterpolation>(*mesh,
                                                                 samplesHessian,
                                                                 AveragingInterpolation::Method::Max,
                                                                 Location::Faces,
                                                                 1.0,
                                                                 false,
                                                                 false,
                                                                 1);

    MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 3;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;
    meshRefinementParameters.min_edge_size = 2.0;
    meshRefinementParameters.account_for_samples_outside = 0;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.refinement_type = 3;
    meshRefinementParameters.smoothing_iterations = 0;

    MeshRefinement meshRefinement(*mesh,
                                  std::move(interpolator),
                                  meshRefinementParameters, false);

    // Execute
    auto undoAction = meshRefinement.Compute();

    // Assert
    ASSERT_EQ(numNodes, mesh->GetNumNodes());
    ASSERT_EQ(numEdges, mesh->GetNumEdges());

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

INSTANTIATE_TEST_SUITE_P(RidgeRefinementTestCases,
                         RidgeRefinementTestCases,
                         ::testing::ValuesIn(RidgeRefinementTestCases::GetData()));

TEST(MeshRefinement, CasulliRefinement)
{
    constexpr double tolerance = 1.0e-12;

    auto curviMesh = MakeCurvilinearGrid(0.0, 0.0, 10.0, 10.0, 3, 3);
    const auto edges = curviMesh->ComputeEdges();
    const auto nodes = curviMesh->ComputeNodes();
    Mesh2D mesh(edges, nodes, Projection::cartesian);

    const std::vector<meshkernel::Point> originalNodes(mesh.Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh.Edges());

    // Expected values were obtained from a mesh refined using the Casulli refinement algorithm
    std::vector<meshkernel::Point> expectedPoints{{0, 0},
                                                  {20, 0},
                                                  {0, 20},
                                                  {20, 20},
                                                  {2.5, 2.5},
                                                  {7.5, 2.5},
                                                  {7.5, 7.5},
                                                  {2.5, 7.5},
                                                  {12.5, 2.5},
                                                  {17.5, 2.5},
                                                  {17.5, 7.5},
                                                  {12.5, 7.5},
                                                  {2.5, 12.5},
                                                  {7.5, 12.5},
                                                  {7.5, 17.5},
                                                  {2.5, 17.5},
                                                  {12.5, 12.5},
                                                  {17.5, 12.5},
                                                  {17.5, 17.5},
                                                  {12.5, 17.5},
                                                  {0, 2.5},
                                                  {0, 7.5},
                                                  {20, 2.5},
                                                  {20, 7.5},
                                                  {0, 12.5},
                                                  {0, 17.5},
                                                  {20, 12.5},
                                                  {20, 17.5},
                                                  {2.5, 0},
                                                  {7.5, 0},
                                                  {12.5, 0},
                                                  {17.5, 0},
                                                  {2.5, 20},
                                                  {7.5, 20},
                                                  {12.5, 20},
                                                  {17.5, 20}};

    // std::vector<int> expectedEdgesStart = {9, 25, 9, 12, 13,
    //                                        10, 13, 16, 27, 14,
    //                                        27, 28, 17, 29, 17,
    //                                        20, 21, 18, 21, 24,
    //                                        31, 22, 31, 32, 33,
    //                                        9, 33, 34, 35, 13,
    //                                        35, 36, 12, 17, 12,
    //                                        11, 16, 21, 16, 15,
    //                                        20, 37, 20, 19, 24,
    //                                        39, 24, 23, 0, 0,
    //                                        34, 2, 2, 26, 28,
    //                                        6, 6, 39, 8, 8};

    // std::vector<int> expectedEdgesEnd = {12, 26, 25, 26, 16,
    //                                      11, 10, 11, 28, 15,
    //                                      14, 15, 20, 30, 29,
    //                                      30, 24, 19, 18, 19,
    //                                      32, 23, 22, 23, 34,
    //                                      10, 9, 10, 36, 14,
    //                                      13, 14, 11, 18, 17,
    //                                      18, 15, 22, 21, 22,
    //                                      19, 38, 37, 38, 23,
    //                                      40, 39, 40, 25, 33,
    //                                      35, 27, 36, 29, 31,
    // 30, 37, 38, 32, 40};

    std::vector<int> expectedEdgesStart = {4, 20, 4, 7, 8, 5, 8,
                                           11, 22, 9, 22, 23, 12,
                                           24, 12, 15, 16, 13, 16,
                                           19, 26, 17, 26, 27, 28,
                                           4, 28, 29, 30, 8, 30, 31,
                                           7, 12, 7, 6, 11, 16, 11,
                                           10, 15, 32, 15, 14, 19,
                                           34, 19, 18, 0, 0, 29, 1,
                                           1, 21, 23, 2, 2, 34, 3, 3};

    std::vector<int> expectedEdgesEnd = {7, 21, 20, 21, 11, 6, 5,
                                         6, 23, 10, 9, 10, 15, 25,
                                         24, 25, 19, 14, 13, 14, 27,
                                         18, 17, 18, 29, 5, 4, 5, 31,
                                         9, 8, 9, 6, 13, 12, 13, 10,
                                         17, 16, 17, 14, 33, 32, 33,
                                         18, 35, 34, 35, 20, 28, 30,
                                         22, 31, 24, 26, 25, 32, 33,
                                         27, 35};

    auto undoAction = CasulliRefinement::Compute(mesh);

    std::vector<meshkernel::UInt> validNodeMap(mesh.GetValidNodeMapping());
    std::vector<meshkernel::UInt> validEdgeMap(mesh.GetValidEdgeMapping());

    ASSERT_EQ(expectedPoints.size(), mesh.GetNumValidNodes());

    for (meshkernel::UInt i = 0; i < expectedPoints.size(); ++i)
    {
        EXPECT_NEAR(expectedPoints[i].x, mesh.Node(validNodeMap[i]).x, tolerance);
        EXPECT_NEAR(expectedPoints[i].y, mesh.Node(validNodeMap[i]).y, tolerance);
    }

    ASSERT_EQ(expectedEdgesStart.size(), mesh.GetNumValidEdges());

    for (meshkernel::UInt i = 0; i < expectedEdgesStart.size(); ++i)
    {
        EXPECT_EQ(expectedEdgesStart[i], mesh.GetEdge(validEdgeMap[i]).first);
        EXPECT_EQ(expectedEdgesEnd[i], mesh.GetEdge(validEdgeMap[i]).second);
    }

    // Test the undo action has been computed correctly
    undoAction->Restore();
    // Recompute faces
    mesh.Administrate();

    ASSERT_EQ(originalNodes.size(), mesh.GetNumValidNodes());
    ASSERT_EQ(originalEdges.size(), mesh.GetNumValidEdges());

    meshkernel::UInt count = 0;

    for (meshkernel::UInt i = 0; i < mesh.Nodes().size(); ++i)
    {
        if (mesh.Node(i).IsValid())
        {
            // Check valid nodes
            EXPECT_EQ(originalNodes[count].x, mesh.Node(i).x);
            EXPECT_EQ(originalNodes[count].y, mesh.Node(i).y);
            ++count;
        }
    }

    count = 0;

    for (meshkernel::UInt i = 0; i < mesh.Edges().size(); ++i)
    {
        if (mesh.IsValidEdge(i))
        {
            EXPECT_EQ(originalEdges[count].first, mesh.GetEdge(i).first);
            EXPECT_EQ(originalEdges[count].second, mesh.GetEdge(i).second);
            ++count;
        }
    }
}

void loadCasulliRefinedMeshData(std::vector<Point>& expectedPoints,
                                std::vector<meshkernel::UInt>& expectedEdgeStart,
                                std::vector<meshkernel::UInt>& expectedEdgeEnd)
{

    const std::string fileName = TEST_FOLDER + "/data/CasulliRefinement/casulli_refinement_patch_with_hole.txt";

    std::ifstream asciiFile;
    asciiFile.open(fileName.c_str());

    for (size_t i = 0; i < expectedPoints.size(); ++i)
    {
        asciiFile >> expectedPoints[i].x;
    }

    for (size_t i = 0; i < expectedPoints.size(); ++i)
    {
        asciiFile >> expectedPoints[i].y;
    }

    for (size_t i = 0; i < expectedEdgeStart.size(); ++i)
    {
        asciiFile >> expectedEdgeStart[i];
    }

    for (size_t i = 0; i < expectedEdgeEnd.size(); ++i)
    {
        asciiFile >> expectedEdgeEnd[i];
    }

    asciiFile.close();
}

TEST(MeshRefinement, CasulliPatchRefinement)
{
    const size_t ExpectedNumberOfPoints = 184;
    const size_t ExpectedNumberOfEdges = 360;

    auto curviMesh = MakeCurvilinearGrid(0.0, 0.0, 20.0, 20.0, 11, 11);
    const auto edges = curviMesh->ComputeEdges();
    const auto nodes = curviMesh->ComputeNodes();

    Mesh2D mesh(edges, nodes, Projection::cartesian);

    const std::vector<meshkernel::Point> originalNodes(mesh.Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh.Edges());

    std::vector<Point> patch{{45.0, 45.0},
                             {155.0, 45.0},
                             {155.0, 155.0},
                             {45.0, 155.0},
                             {45.0, 45.0},
                             {constants::missing::innerOuterSeparator, constants::missing::innerOuterSeparator},
                             {65.0, 65.0},
                             {115.0, 65.0},
                             {115.0, 115.0},
                             {65.0, 115.0},
                             {65.0, 65.0}};

    std::vector<Point> expectedPoints(ExpectedNumberOfPoints);
    std::vector<meshkernel::UInt> expectedEdgeStart(ExpectedNumberOfEdges);
    std::vector<meshkernel::UInt> expectedEdgeEnd(ExpectedNumberOfEdges);

    Polygons polygon(patch, Projection::cartesian);

    auto undoAction = CasulliRefinement::Compute(mesh, polygon);

    constexpr double tolerance = 1.0e-12;

    loadCasulliRefinedMeshData(expectedPoints, expectedEdgeStart, expectedEdgeEnd);

    std::vector<meshkernel::UInt> validNodeMap(mesh.GetValidNodeMapping());
    std::vector<meshkernel::UInt> validEdgeMap(mesh.GetValidEdgeMapping());

    ASSERT_EQ(ExpectedNumberOfPoints, mesh.GetNumValidNodes());
    ASSERT_EQ(ExpectedNumberOfEdges, mesh.GetNumValidEdges());

    for (size_t i = 0; i < expectedPoints.size(); ++i)
    {
        // Map the index i from the nodes array containing only valid points to an array with that may contain in-valid points
        EXPECT_NEAR(expectedPoints[i].x, mesh.Node(validNodeMap[i]).x, tolerance);
        EXPECT_NEAR(expectedPoints[i].y, mesh.Node(validNodeMap[i]).y, tolerance);
    }

    for (size_t i = 0; i < expectedEdgeStart.size(); ++i)
    {
        // Map the index i from the edges array containing only valid edges to an array with that may contain in-valid edges
        EXPECT_EQ(mesh.GetEdge(validEdgeMap[i]).first, expectedEdgeStart[i]);
        EXPECT_EQ(mesh.GetEdge(validEdgeMap[i]).second, expectedEdgeEnd[i]);
    }

    // Test the undo action has been computed correctly
    undoAction->Restore();
    // Recompute faces
    mesh.Administrate();

    ASSERT_EQ(originalNodes.size(), mesh.GetNumValidNodes());
    ASSERT_EQ(originalEdges.size(), mesh.GetNumValidEdges());

    meshkernel::UInt count = 0;

    for (meshkernel::UInt i = 0; i < mesh.Nodes().size(); ++i)
    {
        if (mesh.Node(i).IsValid())
        {
            // Check valid nodes
            EXPECT_EQ(originalNodes[count].x, mesh.Node(i).x);
            EXPECT_EQ(originalNodes[count].y, mesh.Node(i).y);
            ++count;
        }
    }

    count = 0;

    for (meshkernel::UInt i = 0; i < mesh.Edges().size(); ++i)
    {
        if (mesh.IsValidEdge(i))
        {
            EXPECT_EQ(originalEdges[count].first, mesh.GetEdge(i).first);
            EXPECT_EQ(originalEdges[count].second, mesh.GetEdge(i).second);
            ++count;
        }
    }
}

void TestDerefinedMesh(const mk::UInt nx, const mk::UInt ny, const std::string& interactorFileName)
{
    constexpr double tolerance = 1.0e-12;

    auto interactorMesh = ReadLegacyMesh2DFromFile(interactorFileName);

    auto curviMesh = MakeRectangularMeshForTesting(nx, ny, 10.0, Projection::cartesian, {0.0, 0.0},
                                                   true /*ewIndexIncreasing*/,
                                                   true /*nsIndexIncreasing*/);
    Mesh2D mesh(curviMesh->Edges(), curviMesh->Nodes(), Projection::cartesian);
    mesh.Administrate();

    const std::vector<mk::Point> originalNodes(mesh.Nodes());
    const std::vector<mk::Edge> originalEdges(mesh.Edges());

    auto undoAction = meshkernel::CasulliDeRefinement::Compute(mesh);

    const std::vector<mk::Point> refinedNodes(mesh.Nodes());
    const std::vector<mk::Edge> refinedEdges(mesh.Edges());

    //--------------------------------
    // Now compare de-refined mesh with one produced by interactor.

    ASSERT_EQ(mesh.GetNumNodes(), interactorMesh->GetNumNodes());
    ASSERT_EQ(mesh.GetNumEdges(), interactorMesh->GetNumEdges());

    for (UInt i = 0u; i < mesh.GetNumNodes(); ++i)
    {
        EXPECT_NEAR(interactorMesh->Node(i).x, mesh.Node(i).x, tolerance);
        EXPECT_NEAR(interactorMesh->Node(i).y, mesh.Node(i).y, tolerance);
    }

    for (UInt i = 0u; i < mesh.GetNumEdges(); ++i)
    {
        EXPECT_EQ(interactorMesh->GetEdge(i).first, mesh.GetEdge(i).first);
        EXPECT_EQ(interactorMesh->GetEdge(i).second, mesh.GetEdge(i).second);
    }

    //--------------------------------
    // Now test undo
    undoAction->Restore();

    ASSERT_EQ(originalNodes.size(), mesh.Nodes().size());
    ASSERT_EQ(originalEdges.size(), mesh.Edges().size());

    for (mk::UInt i = 0; i < originalNodes.size(); ++i)
    {
        EXPECT_NEAR(originalNodes[i].x, mesh.Node(i).x, tolerance);
        EXPECT_NEAR(originalNodes[i].y, mesh.Node(i).y, tolerance);
    }

    for (mk::UInt i = 0; i < originalEdges.size(); ++i)
    {
        EXPECT_EQ(originalEdges[i].first, mesh.GetEdge(i).first);
        EXPECT_EQ(originalEdges[i].second, mesh.GetEdge(i).second);
    }

    //--------------------------------
    // Now test redo
    undoAction->Commit();

    ASSERT_EQ(refinedNodes.size(), mesh.Nodes().size());
    ASSERT_EQ(refinedEdges.size(), mesh.Edges().size());

    for (mk::UInt i = 0; i < refinedNodes.size(); ++i)
    {
        EXPECT_NEAR(refinedNodes[i].x, mesh.Node(i).x, tolerance);
        EXPECT_NEAR(refinedNodes[i].y, mesh.Node(i).y, tolerance);
    }

    for (mk::UInt i = 0; i < refinedEdges.size(); ++i)
    {
        EXPECT_EQ(refinedEdges[i].first, mesh.GetEdge(i).first);
        EXPECT_EQ(refinedEdges[i].second, mesh.GetEdge(i).second);
    }
}

TEST(MeshRefinement, CasulliDeRefinement)
{
    // de-refine the entire mesh (of different sizes) then compare results with precomputed data.

    const std::string prefix(TEST_FOLDER + "/data/CasulliRefinement/");
    TestDerefinedMesh(21, 21, prefix + "casulli_deref_21_21.nc");
    TestDerefinedMesh(21, 22, prefix + "casulli_deref_21_22.nc");
    TestDerefinedMesh(22, 21, prefix + "casulli_deref_22_21.nc");
    TestDerefinedMesh(22, 22, prefix + "casulli_deref_22_22.nc");
}

TEST(MeshRefinement, CasulliDeRefinementPolygon)
{
    // de-refine the mesh inside a polygon then compare results with precomputed data.

    auto interactorMesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/CasulliRefinement/casulli_derefine_polygon.nc");

    std::vector<Point> points{{55.0, 55.0}, {155.0, 105.0}, {175.0, 225.0}, {25.0, 274.0}, {55.0, 55.0}};
    meshkernel::Polygons polygon(points, Projection::cartesian);

    auto curviMesh = MakeRectangularMeshForTesting(20, 31, 10.0, Projection::cartesian, {0.0, 0.0},
                                                   true /*ewIndexIncreasing*/,
                                                   true /*nsIndexIncreasing*/);
    Mesh2D mesh(curviMesh->Edges(), curviMesh->Nodes(), Projection::cartesian);
    mesh.Administrate();

    auto undoAction = meshkernel::CasulliDeRefinement::Compute(mesh, polygon);

    //--------------------------------
    // Now compare de-refined mesh with one produced by interactor.

    ASSERT_EQ(mesh.GetNumNodes(), interactorMesh->GetNumNodes());
    ASSERT_EQ(mesh.GetNumEdges(), interactorMesh->GetNumEdges());

    constexpr double tolerance = 1.0e-10;

    for (UInt i = 0u; i < mesh.GetNumNodes(); ++i)
    {
        EXPECT_NEAR(interactorMesh->Node(i).x, mesh.Node(i).x, tolerance);
        EXPECT_NEAR(interactorMesh->Node(i).y, mesh.Node(i).y, tolerance);
    }

    for (UInt i = 0u; i < mesh.GetNumEdges(); ++i)
    {
        EXPECT_EQ(interactorMesh->GetEdge(i).first, mesh.GetEdge(i).first);
        EXPECT_EQ(interactorMesh->GetEdge(i).second, mesh.GetEdge(i).second);
    }
}

TEST(MeshRefinement, CasulliDeRefinementPolygonThenAll)
{
    // Step 1 de-refine the mesh inside a polygon
    // Step 2 de-refine the entire mesh
    // compare results with precomputed data.

    auto interactorMesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/CasulliRefinement/casulli_polygon_then_all.nc");

    std::vector<Point> points{{55.0, 55.0}, {155.0, 105.0}, {175.0, 225.0}, {25.0, 274.0}, {55.0, 55.0}};
    meshkernel::Polygons polygon(points, Projection::cartesian);

    auto curviMesh = MakeRectangularMeshForTesting(20, 31, 10.0, Projection::cartesian, {0.0, 0.0},
                                                   true /*ewIndexIncreasing*/,
                                                   true /*nsIndexIncreasing*/);

    Mesh2D mesh(curviMesh->Edges(), curviMesh->Nodes(), Projection::cartesian);
    mesh.Administrate();

    // Derefine on polygon
    auto undoAction = meshkernel::CasulliDeRefinement::Compute(mesh, polygon);

    // Derefine on all
    undoAction = meshkernel::CasulliDeRefinement::Compute(mesh);

    //--------------------------------
    // Now compare de-refined mesh with one produced by interactor.

    ASSERT_EQ(mesh.GetNumNodes(), interactorMesh->GetNumNodes());
    ASSERT_EQ(mesh.GetNumEdges(), interactorMesh->GetNumEdges());

    constexpr double tolerance = 1.0e-10;

    for (UInt i = 0u; i < mesh.GetNumNodes(); ++i)
    {
        EXPECT_NEAR(interactorMesh->Node(i).x, mesh.Node(i).x, tolerance);
        EXPECT_NEAR(interactorMesh->Node(i).y, mesh.Node(i).y, tolerance);
    }

    for (UInt i = 0u; i < mesh.GetNumEdges(); ++i)
    {
        EXPECT_EQ(interactorMesh->GetEdge(i).first, mesh.GetEdge(i).first);
        EXPECT_EQ(interactorMesh->GetEdge(i).second, mesh.GetEdge(i).second);
    }
}

TEST(MeshRefinement, CasulliTwoPolygonDeRefinement)
{
    // Step 1 de-refine the mesh inside a polygon
    // Step 2 de-refine the mesh inside an overlapping polygon
    // compare results with precomputed data.

    constexpr double tolerance = 1.0e-10;
    auto interactorMesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/CasulliRefinement/casulli_dref_two_polygon.nc");

    // Centre of element to be deleted
    std::vector<double> elementCentreX{25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 45.0,
                                       45.0, 45.0, 46.66666666666666, 65.0, 65.0, 66.6666666666666,
                                       85.0, 103.3333333333333, 105.0, 125.0, 125.0};
    std::vector<double> elementCentreY{15.0, 35.0, 55.0, 75.0, 95.0, 115.0, 35.0,
                                       55.0, 75.0, 96.6666666666666, 35.0, 55.0, 76.6666666666666,
                                       55.0, 76.6666666666666, 55.0, 98.3333333333333, 75.0};

    std::vector<Point> centrePoints{{55.0, 55.0}, {155.0, 105.0}, {175.0, 225.0}, {25.0, 274.0}, {55.0, 55.0}};
    meshkernel::Polygons centrePolygon(centrePoints, Projection::cartesian);

    std::vector<Point> lowerPoints{{14.5, 125.0}, {15.5, 14.5}, {95.5, 53.5}, {165.5, 115.5}, {14.5, 125.0}};
    meshkernel::Polygons lowerPolygon(lowerPoints, Projection::cartesian);

    // Ensure that the edges are numbered in the correct order.
    auto curviMesh = MakeRectangularMeshForTesting(20, 31, 10.0, Projection::cartesian, {0.0, 0.0},
                                                   true /*ewIndexIncreasing*/,
                                                   true /*nsIndexIncreasing*/);

    Mesh2D mesh(curviMesh->Edges(), curviMesh->Nodes(), Projection::cartesian);
    mesh.Administrate();

    const std::vector<mk::Point> originalNodes(mesh.Nodes());
    const std::vector<mk::Edge> originalEdges(mesh.Edges());

    // Derefine on centre polygon
    auto undoCentrePolygon = meshkernel::CasulliDeRefinement::Compute(mesh, centrePolygon);

    const std::vector<mk::Point> centreRefinedNodes(mesh.Nodes());
    const std::vector<mk::Edge> centreRefinedEdges(mesh.Edges());

    mesh.ComputeCircumcentersMassCentersAndFaceAreas(true);

    // Get the element centres of the elements to be deleted.
    std::vector<meshkernel::Point> toDelete(meshkernel::CasulliDeRefinement::ElementsToDelete(mesh, lowerPolygon));

    // Compare the elements to be deleted by the lowerPolygon with expected data.
    ASSERT_EQ(toDelete.size(), elementCentreX.size());

    for (size_t i = 0u; i < elementCentreX.size(); ++i)
    {
        EXPECT_NEAR(elementCentreX[i], toDelete[i].x, tolerance);
        EXPECT_NEAR(elementCentreY[i], toDelete[i].y, tolerance);
    }

    // Derefine on lower polygon
    auto undoLowerPolygon = meshkernel::CasulliDeRefinement::Compute(mesh, lowerPolygon);

    const std::vector<mk::Point> lowerRefinedNodes(mesh.Nodes());
    const std::vector<mk::Edge> lowerRefinedEdges(mesh.Edges());

    //--------------------------------
    // Now compare de-refined mesh with one produced earlier
    // Originally compared with interactor results.

    ASSERT_EQ(mesh.GetNumNodes(), interactorMesh->GetNumNodes());
    ASSERT_EQ(mesh.GetNumEdges(), interactorMesh->GetNumEdges());

    for (UInt i = 0u; i < mesh.GetNumNodes(); ++i)
    {
        EXPECT_NEAR(interactorMesh->Node(i).x, mesh.Node(i).x, tolerance);
        EXPECT_NEAR(interactorMesh->Node(i).y, mesh.Node(i).y, tolerance);
    }

    for (UInt i = 0u; i < mesh.GetNumEdges(); ++i)
    {
        EXPECT_EQ(interactorMesh->GetEdge(i).first, mesh.GetEdge(i).first);
        EXPECT_EQ(interactorMesh->GetEdge(i).second, mesh.GetEdge(i).second);
    }

    //--------------------------------
    // Now test undo
    // First undo the de-refinement inside the lower polygon
    undoLowerPolygon->Restore();

    ASSERT_EQ(centreRefinedNodes.size(), mesh.Nodes().size());
    ASSERT_EQ(centreRefinedEdges.size(), mesh.Edges().size());

    for (mk::UInt i = 0; i < centreRefinedNodes.size(); ++i)
    {
        EXPECT_NEAR(centreRefinedNodes[i].x, mesh.Node(i).x, tolerance);
        EXPECT_NEAR(centreRefinedNodes[i].y, mesh.Node(i).y, tolerance);
    }

    for (mk::UInt i = 0; i < centreRefinedEdges.size(); ++i)
    {
        EXPECT_EQ(centreRefinedEdges[i].first, mesh.GetEdge(i).first);
        EXPECT_EQ(centreRefinedEdges[i].second, mesh.GetEdge(i).second);
    }

    // Next undo the de-refinement inside the centre polygon
    undoCentrePolygon->Restore();

    ASSERT_EQ(originalNodes.size(), mesh.Nodes().size());
    ASSERT_EQ(originalEdges.size(), mesh.Edges().size());

    for (mk::UInt i = 0; i < originalNodes.size(); ++i)
    {
        EXPECT_NEAR(originalNodes[i].x, mesh.Node(i).x, tolerance);
        EXPECT_NEAR(originalNodes[i].y, mesh.Node(i).y, tolerance);
    }

    for (mk::UInt i = 0; i < originalEdges.size(); ++i)
    {
        EXPECT_EQ(originalEdges[i].first, mesh.GetEdge(i).first);
        EXPECT_EQ(originalEdges[i].second, mesh.GetEdge(i).second);
    }

    //--------------------------------
    // Now test redo

    // First redo the de-refinement inside the centre polygon
    undoCentrePolygon->Commit();

    ASSERT_EQ(centreRefinedNodes.size(), mesh.Nodes().size());
    ASSERT_EQ(centreRefinedEdges.size(), mesh.Edges().size());

    for (mk::UInt i = 0; i < centreRefinedNodes.size(); ++i)
    {
        EXPECT_NEAR(centreRefinedNodes[i].x, mesh.Node(i).x, tolerance);
        EXPECT_NEAR(centreRefinedNodes[i].y, mesh.Node(i).y, tolerance);
    }

    for (mk::UInt i = 0; i < centreRefinedEdges.size(); ++i)
    {
        EXPECT_EQ(centreRefinedEdges[i].first, mesh.GetEdge(i).first);
        EXPECT_EQ(centreRefinedEdges[i].second, mesh.GetEdge(i).second);
    }

    // Next redo the de-refinement inside the lower polygon

    undoLowerPolygon->Commit();

    ASSERT_EQ(lowerRefinedNodes.size(), mesh.Nodes().size());
    ASSERT_EQ(lowerRefinedEdges.size(), mesh.Edges().size());

    for (mk::UInt i = 0; i < lowerRefinedNodes.size(); ++i)
    {
        EXPECT_NEAR(lowerRefinedNodes[i].x, mesh.Node(i).x, tolerance);
        EXPECT_NEAR(lowerRefinedNodes[i].y, mesh.Node(i).y, tolerance);
    }

    for (mk::UInt i = 0; i < lowerRefinedEdges.size(); ++i)
    {
        EXPECT_EQ(lowerRefinedEdges[i].first, mesh.GetEdge(i).first);
        EXPECT_EQ(lowerRefinedEdges[i].second, mesh.GetEdge(i).second);
    }
}

TEST(MeshRefinement, SplitBoundariesOfSmallMesh)
{
    constexpr double tolerance = 1.0e-12;

    auto curviMesh = MakeCurvilinearGrid(0.0, 0.0, 10.0, 10.0, 4, 4);
    Mesh2D mesh(curviMesh->ComputeEdges(), curviMesh->ComputeNodes(), Projection::cartesian);
    mesh.Administrate();

    SplitRowColumnOfMesh splitMesh;
    UndoActionStack undoStack;

    //--------------------------------
    // Split along south boundary of the domain
    UInt node1 = mesh.FindNodeCloseToAPoint({0.0, 0.0}, 1.0e-4);
    UInt node2 = mesh.FindNodeCloseToAPoint({0.0, 10.0}, 1.0e-4);

    UInt edge = mesh.FindEdge(node1, node2);
    undoStack.Add(splitMesh.Compute(mesh, edge));

    //--------------------------------
    // Split along east boundary of the domain
    node1 = mesh.FindNodeCloseToAPoint({20.0, 0.0}, 1.0e-4);
    node2 = mesh.FindNodeCloseToAPoint({30.0, 0.0}, 1.0e-4);

    edge = mesh.FindEdge(node1, node2);
    undoStack.Add(splitMesh.Compute(mesh, edge));

    //--------------------------------
    // Split along north boundary of the domain
    node1 = mesh.FindNodeCloseToAPoint({30.0, 20.0}, 1.0e-4);
    node2 = mesh.FindNodeCloseToAPoint({30.0, 30.0}, 1.0e-4);

    edge = mesh.FindEdge(node1, node2);
    undoStack.Add(splitMesh.Compute(mesh, edge));

    //--------------------------------
    // Split along west boundary of the domain
    node1 = mesh.FindNodeCloseToAPoint({0.0, 30.0}, 1.0e-4);
    node2 = mesh.FindNodeCloseToAPoint({10.0, 30.0}, 1.0e-4);

    edge = mesh.FindEdge(node1, node2);
    undoStack.Add(splitMesh.Compute(mesh, edge));

    //--------------------------------

    std::vector<double> expectedX{0.0, 10.0, 20.0, 30.0, 0.0, 10.0, 20.0, 30.0, 0.0, 10.0,
                                  20.0, 30.0, 0.0, 10.0, 20.0, 30.0, 0.0, 10.0, 20.0, 30.0,
                                  25.0, 25.0, 25.0, 25.0, 25.0, 30.0, 25.0, 20.0, 10.0, 0.0,
                                  5.0, 5.0, 5.0, 5.0, 5.0, 5.0};

    std::vector<double> expectedY{0.0, 0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 10.0, 20.0, 20.0,
                                  20.0, 20.0, 30.0, 30.0, 30.0, 30.0, 5.0, 5.0, 5.0, 5.0,
                                  0.0, 5.0, 10.0, 20.0, 30.0, 25.0, 25.0, 25.0, 25.0, 25.0,
                                  30.0, 25.0, 20.0, 10.0, 5.0, 0.0};

    ASSERT_EQ(static_cast<UInt>(expectedX.size()), mesh.GetNumNodes());

    for (UInt i = 0; i < mesh.GetNumNodes(); ++i)
    {
        EXPECT_NEAR(expectedX[i], mesh.Node(i).x, tolerance);
        EXPECT_NEAR(expectedY[i], mesh.Node(i).y, tolerance);
    }
}

TEST(MeshRefinement, PartialSplittingOfRow)
{
    constexpr double tolerance = 1.0e-12;

    auto curviMesh = MakeCurvilinearGrid(0.0, 0.0, 10.0, 10.0, 5, 5);
    Mesh2D mesh(curviMesh->ComputeEdges(), curviMesh->ComputeNodes(), Projection::cartesian);
    mesh.Administrate();

    SplitRowColumnOfMesh splitMesh;
    UndoActionStack undoStack;

    // Add triangle into mesh at south-east corner

    UInt node1 = mesh.FindNodeCloseToAPoint({30.0, 10.0}, 1.0e-4);
    UInt node2 = mesh.FindNodeCloseToAPoint({40.0, 0.0}, 1.0e-4);

    [[maybe_unused]] auto undo = mesh.ConnectNodes(node1, node2);

    mesh.Administrate();

    //--------------------------------
    // Split along south boundary of the domain from west side
    node1 = mesh.FindNodeCloseToAPoint({0.0, 0.0}, 1.0e-4);
    node2 = mesh.FindNodeCloseToAPoint({0.0, 10.0}, 1.0e-4);

    UInt edge = mesh.FindEdge(node1, node2);
    undoStack.Add(splitMesh.Compute(mesh, edge));

    //--------------------------------
    // Split along east boundary of the domain from north side
    node1 = mesh.FindNodeCloseToAPoint({30.0, 30.0}, 1.0e-4);
    node2 = mesh.FindNodeCloseToAPoint({40.0, 30.0}, 1.0e-4);

    edge = mesh.FindEdge(node1, node2);
    undoStack.Add(splitMesh.Compute(mesh, edge));

    std::vector<double> expectedX{0.0, 10.0, 20.0, 30.0, 40.0, 0.0, 10.0, 20.0, 30.0,
                                  40.0, 0.0, 10.0, 20.0, 30.0, 40.0, 0.0, 10.0, 20.0,
                                  30.0, 40.0, 0.0, 10.0, 20.0, 30.0, 40.0, 0.0, 10.0,
                                  20.0, 35.0, 35.0, 35.0};

    std::vector<double> expectedY{0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 10.0, 10.0,
                                  20.0, 20.0, 20.0, 20.0, 20.0, 30.0, 30.0, 30.0, 30.0,
                                  30.0, 40.0, 40.0, 40.0, 40.0, 40.0, 5.0, 5.0, 5.0, 40.0,
                                  30.0, 20.0};

    ASSERT_EQ(static_cast<UInt>(expectedX.size()), mesh.GetNumNodes());

    for (UInt i = 0; i < mesh.GetNumNodes(); ++i)
    {
        EXPECT_NEAR(expectedX[i], mesh.Node(i).x, tolerance);
        EXPECT_NEAR(expectedY[i], mesh.Node(i).y, tolerance);
    }
}

TEST(MeshRefinement, PartialSplittingOfRowTriangleMidWay)
{
    constexpr double tolerance = 1.0e-12;

    auto curviMesh = MakeCurvilinearGrid(0.0, 0.0, 10.0, 10.0, 5, 5);
    Mesh2D mesh(curviMesh->ComputeEdges(), curviMesh->ComputeNodes(), Projection::cartesian);
    mesh.Administrate();

    SplitRowColumnOfMesh splitMesh;
    UndoActionStack undoStack;

    // Add triangle into mesh at south-east corner

    UInt node1 = mesh.FindNodeCloseToAPoint({20.0, 30.0}, 1.0e-4);
    UInt node2 = mesh.FindNodeCloseToAPoint({30.0, 20.0}, 1.0e-4);

    [[maybe_unused]] auto undo = mesh.ConnectNodes(node1, node2);

    mesh.Administrate();

    //--------------------------------
    // Split along south boundary of the domain from west side
    node1 = mesh.FindNodeCloseToAPoint({0.0, 30.0}, 1.0e-4);
    node2 = mesh.FindNodeCloseToAPoint({0.0, 20.0}, 1.0e-4);

    UInt edge = mesh.FindEdge(node1, node2);
    undoStack.Add(splitMesh.Compute(mesh, edge));

    //--------------------------------
    // Split along east boundary of the domain from south side
    node1 = mesh.FindNodeCloseToAPoint({20.0, 0.0}, 1.0e-4);
    node2 = mesh.FindNodeCloseToAPoint({30.0, 0.0}, 1.0e-4);

    edge = mesh.FindEdge(node1, node2);
    undoStack.Add(splitMesh.Compute(mesh, edge));

    //--------------------------------
    // Split along east boundary of the domain from north side
    node1 = mesh.FindNodeCloseToAPoint({20.0, 40.0}, 1.0e-4);
    node2 = mesh.FindNodeCloseToAPoint({30.0, 40.0}, 1.0e-4);

    edge = mesh.FindEdge(node1, node2);
    undoStack.Add(splitMesh.Compute(mesh, edge));

    //--------------------------------
    // Split along east boundary of the domain from east side
    node1 = mesh.FindNodeCloseToAPoint({40.0, 30.0}, 1.0e-4);
    node2 = mesh.FindNodeCloseToAPoint({40.0, 20.0}, 1.0e-4);

    edge = mesh.FindEdge(node1, node2);
    undoStack.Add(splitMesh.Compute(mesh, edge));

    //--------------------------------

    std::vector<double> expectedX{0.0, 10.0, 20.0, 30.0, 40.0, 0.0, 10.0,
                                  20.0, 30.0, 40.0, 0.0, 10.0, 20.0, 30.0,
                                  40.0, 0.0, 10.0, 20.0, 30.0, 40.0, 0.0,
                                  10.0, 20.0, 30.0, 40.0, 0.0, 10.0, 25.0,
                                  25.0, 25.0, 40.0};

    std::vector<double> expectedY{0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 10.0, 10.0,
                                  10.0, 10.0, 20.0, 20.0, 20.0, 20.0, 20.0,
                                  30.0, 30.0, 30.0, 30.0, 30.0, 40.0, 40.0,
                                  40.0, 40.0, 40.0, 25.0, 25.0, 0.0, 10.0,
                                  40.0, 25.0};

    ASSERT_EQ(static_cast<UInt>(expectedX.size()), mesh.GetNumNodes());

    for (UInt i = 0; i < mesh.GetNumNodes(); ++i)
    {
        EXPECT_NEAR(expectedX[i], mesh.Node(i).x, tolerance);
        EXPECT_NEAR(expectedY[i], mesh.Node(i).y, tolerance);
    }
}

TEST(MeshRefinement, PartialSplittingOfRowBoundedByTriangles)
{
    constexpr double tolerance = 1.0e-12;

    auto curviMesh = MakeCurvilinearGrid(0.0, 0.0, 10.0, 10.0, 7, 4);
    Mesh2D mesh(curviMesh->ComputeEdges(), curviMesh->ComputeNodes(), Projection::cartesian);
    mesh.Administrate();

    SplitRowColumnOfMesh splitMesh;
    UndoActionStack undoStack;

    // Add triangle into mesh at south-east corner

    UInt node1 = mesh.FindNodeCloseToAPoint({0.0, 10.0}, 1.0e-4);
    UInt node2 = mesh.FindNodeCloseToAPoint({10.0, 20.0}, 1.0e-4);

    [[maybe_unused]] auto undo = mesh.ConnectNodes(node1, node2);

    node1 = mesh.FindNodeCloseToAPoint({50.0, 10.0}, 1.0e-4);
    node2 = mesh.FindNodeCloseToAPoint({60.0, 20.0}, 1.0e-4);

    undo = mesh.ConnectNodes(node1, node2);

    mesh.Administrate();

    //--------------------------------
    // Split along south boundary of the domain from middle of the domain
    node1 = mesh.FindNodeCloseToAPoint({30.0, 10.0}, 1.0e-4);
    node2 = mesh.FindNodeCloseToAPoint({30.0, 20.0}, 1.0e-4);

    UInt edge = mesh.FindEdge(node1, node2);
    undoStack.Add(splitMesh.Compute(mesh, edge));

    std::vector<double> expectedX{0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0,
                                  0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0,
                                  0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0,
                                  0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0,
                                  40.0, 30.0, 20.0};

    std::vector<double> expectedY{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0,
                                  10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 20.0,
                                  20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 30.0,
                                  30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 15.0,
                                  15.0, 15.0};

    ASSERT_EQ(static_cast<UInt>(expectedX.size()), mesh.GetNumNodes());

    for (UInt i = 0; i < mesh.GetNumNodes(); ++i)
    {
        EXPECT_NEAR(expectedX[i], mesh.Node(i).x, tolerance);
        EXPECT_NEAR(expectedY[i], mesh.Node(i).y, tolerance);
    }
}

TEST(MeshRefinement, SplitBoundariesOfMesh)
{
    // Refine along all four boundary edges of this mesh
    // Only test for the number of nodes, edges and elements.

    auto curviMesh = MakeCurvilinearGrid(0.0, 0.0, 10.0, 10.0, 11, 11);
    Mesh2D mesh(curviMesh->ComputeEdges(), curviMesh->ComputeNodes(), Projection::cartesian);
    mesh.Administrate();

    ASSERT_EQ(mesh.GetNumNodes(), 121);
    ASSERT_EQ(mesh.GetNumEdges(), 220);
    ASSERT_EQ(mesh.GetNumFaces(), 100);

    SplitRowColumnOfMesh splitMesh;
    UndoActionStack undoStack;

    //--------------------------------
    // Split along south boundary of the domain
    UInt node1 = mesh.FindNodeCloseToAPoint({0.0, 0.0}, 1.0e-4);
    UInt node2 = mesh.FindNodeCloseToAPoint({0.0, 10.0}, 1.0e-4);

    UInt edge = mesh.FindEdge(node1, node2);
    undoStack.Add(splitMesh.Compute(mesh, edge));

    ASSERT_EQ(mesh.GetNumValidNodes(), 132);
    ASSERT_EQ(mesh.GetNumValidEdges(), 241);
    ASSERT_EQ(mesh.GetNumFaces(), 110);

    //--------------------------------
    // Split along east boundary of the domain
    node1 = mesh.FindNodeCloseToAPoint({90.0, 0.0}, 1.0e-4);
    node2 = mesh.FindNodeCloseToAPoint({100.0, 0.0}, 1.0e-4);

    edge = mesh.FindEdge(node1, node2);
    undoStack.Add(splitMesh.Compute(mesh, edge));

    ASSERT_EQ(mesh.GetNumValidNodes(), 144);
    ASSERT_EQ(mesh.GetNumValidEdges(), 264);
    ASSERT_EQ(mesh.GetNumFaces(), 121);

    //--------------------------------
    // Split along north boundary of the domain
    node1 = mesh.FindNodeCloseToAPoint({100.0, 90.0}, 1.0e-4);
    node2 = mesh.FindNodeCloseToAPoint({100.0, 100.0}, 1.0e-4);

    edge = mesh.FindEdge(node1, node2);
    undoStack.Add(splitMesh.Compute(mesh, edge));

    ASSERT_EQ(mesh.GetNumValidNodes(), 156);
    ASSERT_EQ(mesh.GetNumValidEdges(), 287);
    ASSERT_EQ(mesh.GetNumFaces(), 132);

    //--------------------------------
    // Split along west boundary of the domain
    node1 = mesh.FindNodeCloseToAPoint({0.0, 100.0}, 1.0e-4);
    node2 = mesh.FindNodeCloseToAPoint({10.0, 100.0}, 1.0e-4);

    edge = mesh.FindEdge(node1, node2);
    undoStack.Add(splitMesh.Compute(mesh, edge));

    ASSERT_EQ(mesh.GetNumValidNodes(), 169);
    ASSERT_EQ(mesh.GetNumValidEdges(), 312);
    ASSERT_EQ(mesh.GetNumFaces(), 144);
}

void GenerateGridWithLoop(std::vector<Point>& nodes, std::vector<Edge>& edges)
{

    nodes.push_back({-10.0, -10.0});
    nodes.push_back({0.0, -10.0});
    nodes.push_back({10.0, -10.0});

    nodes.push_back({-10.0, 0.0});
    nodes.push_back({0.0, 0.0});
    nodes.push_back({10.0, 0.0});

    nodes.push_back({-10.0, 10.0});
    nodes.push_back({0.0, 10.0});
    nodes.push_back({10.0, 10.0});

    nodes.push_back({-20.0, -20.0});
    nodes.push_back({0.0, -20.0});
    nodes.push_back({20.0, -20.0});

    nodes.push_back({-20.0, 0.0});
    nodes.push_back({20.0, 0.0});

    nodes.push_back({-20.0, 20.0});
    nodes.push_back({0.0, 20.0});
    nodes.push_back({20.0, 20.0});

    edges.push_back({0, 1});
    edges.push_back({1, 2});
    edges.push_back({3, 4});
    edges.push_back({4, 5});
    edges.push_back({6, 7});
    edges.push_back({7, 8});

    edges.push_back({0, 3});
    edges.push_back({1, 4});
    edges.push_back({2, 5});
    edges.push_back({3, 6});
    edges.push_back({4, 7});
    edges.push_back({5, 8});

    edges.push_back({9, 10});
    edges.push_back({10, 11});

    edges.push_back({9, 0});
    edges.push_back({10, 1});
    edges.push_back({11, 2});

    edges.push_back({9, 12});
    edges.push_back({12, 3});
    edges.push_back({11, 13});
    edges.push_back({5, 13});

    edges.push_back({12, 14});
    edges.push_back({14, 6});
    edges.push_back({14, 15});
    edges.push_back({7, 15});
    edges.push_back({15, 16});
    edges.push_back({13, 16});
    edges.push_back({8, 16});
}

TEST(MeshRefinement, SplitElementLoop)
{
    constexpr double tolerance = 1.0e-12;

    std::vector<Point> nodes;
    std::vector<Edge> edges;
    GenerateGridWithLoop(nodes, edges);
    Mesh2D mesh(edges, nodes, Projection::cartesian);

    mesh.Administrate();

    ASSERT_EQ(mesh.GetNumValidNodes(), 17);
    ASSERT_EQ(mesh.GetNumValidEdges(), 28);
    ASSERT_EQ(mesh.GetNumFaces(), 12);

    SplitRowColumnOfMesh splitMesh;
    auto undoSplit = splitMesh.Compute(mesh, 14);

    ASSERT_EQ(mesh.GetNumNodes(), 25);
    ASSERT_EQ(mesh.GetNumEdges(), 52);
    ASSERT_EQ(mesh.GetNumFaces(), 20);

    std::vector<double> expectedX{-10.0, 0.0, 10.0, -10.0, 0.0, 10.0, -10.0, 0.0, 10.0,
                                  -20.0, 0.0, 20.0, -20.0, 20.0, -20.0, 0.0, 20.0, -15.0,
                                  -15.0, -15.0, 0.0, 15.0, 15.0, 15.0, 0.0};

    std::vector<double> expectedY{-10.0, -10.0, -10.0, 0.0, 0.0, 0.0, 10.0, 10.0, 10.0,
                                  -20.0, -20.0, -20.0, 0.0, 0.0, 20.0, 20.0, 20.0, -15.0,
                                  0.0, 15.0, 15.0, 15.0, 0.0, -15.0, -15.0};

    const std::vector<Point>& meshNodes = mesh.Nodes();

    for (UInt i = 0; i < nodes.size(); ++i)
    {
        EXPECT_NEAR(expectedX[i], meshNodes[i].x, tolerance);
        EXPECT_NEAR(expectedY[i], meshNodes[i].y, tolerance);
    }

    const UInt nullValue = constants::missing::uintValue;

    std::vector<UInt> expectedEdgesFirst{0, 1, 3, 4, 6, 7, 0, 1, 2, 3, 4, 5,
                                         9, 10, nullValue, nullValue, nullValue, 9, nullValue, 11, nullValue, 12, nullValue, 14,
                                         nullValue, 15, 13, nullValue, 9, 17, 12, 18, 17, 14, 19, 18,
                                         7, 20, 19, 8, 21, 20, 5, 22};

    std::vector<UInt> expectedEdgesSecond{1, 2, 4, 5, 7, 8, 3, 4, 5, 6, 7, 8, 10,
                                          11, nullValue, nullValue, nullValue, 12, nullValue, 13, nullValue, 14, nullValue, 15, nullValue, 16,
                                          16, nullValue, 17, 0, 18, 3, 18, 19, 6, 19, 20, 15, 20,
                                          21, 16, 21, 22, 13};

    for (UInt i = 0; i < mesh.GetNumValidEdges(); ++i)
    {
        EXPECT_EQ(expectedEdgesFirst[i], mesh.GetEdge(i).first);
        EXPECT_EQ(expectedEdgesSecond[i], mesh.GetEdge(i).second);
    }

    //--------------------------------
    // Undo split

    undoSplit->Restore();
    mesh.Administrate();

    ASSERT_EQ(mesh.GetNumValidNodes(), 17);
    ASSERT_EQ(mesh.GetNumValidEdges(), 28);
    ASSERT_EQ(mesh.GetNumFaces(), 12);
}

TEST(MeshRefinement, RowSplittingFailureTests)
{
    auto curviMesh = MakeCurvilinearGrid(0.0, 0.0, 10.0, 10.0, 11, 11);
    Mesh2D mesh(curviMesh->ComputeEdges(), curviMesh->ComputeNodes(), Projection::cartesian);
    mesh.Administrate();

    mesh.Administrate();
    SplitRowColumnOfMesh splitMeshRow;

    // Edge id is the null value
    EXPECT_THROW([[maybe_unused]] auto undo1 = splitMeshRow.Compute(mesh, constants::missing::uintValue), ConstraintError);
    // Out of bounds edge id
    EXPECT_THROW([[maybe_unused]] auto undo2 = splitMeshRow.Compute(mesh, mesh.GetNumEdges() + 10), ConstraintError);

    auto [newNodeId, undo3] = mesh.InsertNode({110.0, 0.0});
    UInt endNode = mesh.FindNodeCloseToAPoint({100.0, 0.0}, 1.0e-4);

    auto [edgeId, undoInsertEdge] = mesh.ConnectNodes(newNodeId, endNode);
    // Undo  edge insertion, so the edgeId should now be invalid
    undoInsertEdge->Restore();

    // Invalid edge id
    EXPECT_THROW([[maybe_unused]] auto undo4 = splitMeshRow.Compute(mesh, edgeId), ConstraintError);
}

TEST(MeshRefinement, CasulliRefinementBasedOnDepth)
{

    const double delta = 100.0;
    const size_t nodeCount = 21;

    const double limit = delta * static_cast<double>(nodeCount - 1);

    auto curviMesh = MakeCurvilinearGrid(0.0, 0.0, delta, delta, nodeCount, nodeCount);
    const auto edges = curviMesh->ComputeEdges();
    const auto nodes = curviMesh->ComputeNodes();
    Mesh2D mesh(edges, nodes, Projection::cartesian);
    mesh.Administrate();
    mesh.ComputeEdgesCenters();
    mesh.ComputeEdgesLengths();

    std::vector<double> depth(mesh.GetNumNodes());
    mk::Polygons polygon;

    [[maybe_unused]] auto frankesFunction = [limit](const double xp, const double yp)
    {
        const double x = xp / limit;
        const double y = yp / limit;

        double result = 0.75 * std::exp(-0.25 * std::pow(9.0 * x - 2.0, 2) - 0.25 * std::pow(9.0 * y - 2.0, 2)) +
                        0.75 * std::exp(-1.0 / 49.0 * std::pow(9.0 * x - 2.0, 2) - 0.1 * std::pow(9.0 * y - 2.0, 2)) +
                        0.5 * std::exp(-0.25 * std::pow(9.0 * x - 7.0, 2) - 0.25 * std::pow(9.0 * y - 3.0, 2)) -
                        0.2 * std::exp(-std::pow(9.0 * x - 4.0, 2) - std::pow(9.0 * y - 7.0, 2));

        // return result;
        return 100.0 * (1.50053 - result + 0.1);
    };

    [[maybe_unused]] auto sinShore = [limit](const double x, const double y)
    {
        const double seaDepth = 1000.0;
        const double shoreDepth = 2.0;
        double range = limit * 0.3;
        double sx = 0.5 * (std::sin(x / limit * 3.1416 * 4.0) + 1.0);
        double v = 0.6 * limit + range * sx;
        double result;

        if (y > v)
        {
            result = shoreDepth;
        }
        else
        {
            double l = (v - y) / v;

            if (l < 0.4)
            {
                result = (1.0 - l) * shoreDepth + 0.2 * seaDepth * l;
            }
            else
            {
                result = (1.0 - l) * shoreDepth + seaDepth * l;
            }
        }

        return result;
    };

    double max = -1000000.0;
    double min = 1000000.0;

    for (mk::UInt i = 0; i < depth.size(); ++i)
    {
        const auto pnt = mesh.Node(i);
        depth[i] = sinShore(pnt.x, pnt.y);
        max = std::max(max, depth[i]);
        min = std::min(min, depth[i]);
    }

    mk::SampleTriangulationInterpolator depthInterpolator(mesh.Nodes(), mesh.m_projection);
    depthInterpolator.SetData(1, depth);

    MeshRefinementParameters refinementParameters;
    refinementParameters.min_edge_size = 6.5; // 12.5;
    refinementParameters.max_courant_time = 5.0;

    auto undo = mk::CasulliRefinement::Compute(mesh, polygon, depthInterpolator, 1, refinementParameters);
}

TEST(MeshRefinement, CasulliRefinementBasedOnDepthReal)
{
    // return;

    //--------------------------------

    static const size_t MaxSize = 2550923;
    std::vector<double> xNodes(MaxSize);
    std::vector<double> yNodes(MaxSize);
    std::vector<double> depths(MaxSize);

    std::string fileName = TEST_FOLDER + "/data/CasulliRefinement/stpete.xyz";

    std::ifstream asciiFile;
    asciiFile.open(fileName.c_str());

    double maxX = -1.0e30;
    double minX = 1.0e30;
    double maxY = -1.0e30;
    double minY = 1.0e30;
    double maxD = -1.0e30;
    double minD = 1.0e30;

    for (size_t i = 0; i < MaxSize; ++i)
    {
        mk::Point pnt;

        asciiFile >> xNodes[i];
        asciiFile >> yNodes[i];

        // asciiFile >> pnt.x;
        // asciiFile >> pnt.y;

        pnt.x = xNodes[i];
        pnt.y = yNodes[i];
        asciiFile >> depths[i];

        mk::Cartesian3DPoint cpnt = mk::SphericalToCartesian3D(pnt);

        xNodes[i] = cpnt.x;
        yNodes[i] = cpnt.y;

        maxX = std::max(maxX, xNodes[i]);
        minX = std::min(minX, xNodes[i]);

        maxY = std::max(maxY, yNodes[i]);
        minY = std::min(minY, yNodes[i]);

        maxD = std::max(maxD, depths[i]);
        minD = std::min(minD, depths[i]);
    }

    std::cout.precision(16);

    mk::InterpolationParameters interpolationParameters;
    interpolationParameters.m_relativeSearchRadius = 1.0;
    interpolationParameters.m_useClosestIfNoneFound = true;

    mk::SampleAveragingInterpolator depthInterpolator(xNodes, yNodes, mk::Projection::cartesian, interpolationParameters);
    mk::SampleTriangulationInterpolator depthInterpolatorForMesh(xNodes, yNodes, mk::Projection::cartesian);
    depthInterpolator.SetData(1, depths);
    depthInterpolatorForMesh.SetData(1, depths);
    mk::Polygons polygon;

    const size_t nodeCount = 72;
    const double deltaX = (maxX - minX) / static_cast<double>(nodeCount - 1);
    const double deltaY = (maxY - minY) / static_cast<double>(nodeCount - 1);
    const double delta = std::min(deltaX, deltaY);

    // auto curviMesh = MakeCurvilinearGrid(0.0, 0.0, delta, delta, nodeCount, nodeCount);
    auto curviMesh = MakeCurvilinearGrid(minX, minY, delta, delta, nodeCount, nodeCount);

    auto edges = curviMesh->ComputeEdges();
    auto nodes = curviMesh->ComputeNodes();

    for (size_t i = 0; i < nodes.size(); ++i)
    {
        // if (depthInterpolator.InterpolateValue(1, nodes[i]) == mk::constants::missing::doubleValue)
        if (depthInterpolatorForMesh.InterpolateValue(1, nodes[i]) == mk::constants::missing::doubleValue)
        {
            nodes[i].SetInvalid();
        }
    }

    Mesh2D mesh(edges, nodes, Projection::cartesian);
    mesh.Administrate();

    mk::Edge edge{mk::constants::missing::uintValue, mk::constants::missing::uintValue};

    for (mk::UInt i = 0; i < mesh.GetNumEdges(); ++i)
    {
        if (mesh.m_edgesNumFaces[i] == 0)
        {
            [[maybe_unused]] auto undoAction = mesh.ResetEdge(i, edge);
        }
    }

    mesh.Administrate();
    mesh.ComputeEdgesCenters();
    mesh.ComputeEdgesLengths();

    RemoveDisconnectedRegions removeDisconnectedRegions;

    [[maybe_unused]] auto undoAction = removeDisconnectedRegions.Compute(mesh);
    mesh.Administrate();
    mesh.ComputeEdgesCenters();
    mesh.ComputeEdgesLengths();

    for (mk::UInt i = 0; i < mesh.GetNumEdges(); ++i)
    {
        if (mesh.m_edgesNumFaces[i] == 0)
        {
            [[maybe_unused]] auto undoAction = mesh.ResetEdge(i, edge);
        }
    }
    mesh.Administrate();
    mesh.ComputeEdgesCenters();
    mesh.ComputeEdgesLengths();

    std::vector<mk::Edge> usedEdges;
    usedEdges.reserve(mesh.GetNumEdges());

    for (mk::UInt i = 0; i < mesh.GetNumEdges(); ++i)
    {
        if (mesh.m_edgesNumFaces[i] != 0)
        {
            usedEdges.push_back(mesh.GetEdge(i));
        }
    }

    Mesh2D mesh2(usedEdges, mesh.Nodes(), Projection::cartesian);
    mesh2.Administrate();
    mesh2.ComputeEdgesCenters();
    mesh2.ComputeEdgesLengths();

    MeshRefinementParameters refinementParameters;
    refinementParameters.max_num_refinement_iterations = 1;
    refinementParameters.min_edge_size = 6.5; // 12.5;
    refinementParameters.max_courant_time = 100.0;
    refinementParameters.minimum_refinement_depth = -2.0;

    auto undo = mk::CasulliRefinement::Compute(mesh2, polygon, depthInterpolator, 1, refinementParameters);
    mesh2.Administrate();
    mesh2.ComputeEdgesCenters();
    mesh2.ComputeEdgesLengths();

    mk::PrintVtk(mesh2.Nodes(), mesh2.m_facesNodes, "meshdata.vtu");

    // auto ortho = mesh2.GetOrthogonality();

    // double maxOrtho = *std::max_element(ortho.begin(), ortho.end());
    // double minOrtho = 1000.0;

    // for (size_t i = 0; i < ortho.size(); ++i)
    // {

    //     std::cout << std::setw(20) << ortho[i] << "  ";

    //     if ((i + 1) % 10 == 0)
    //     {
    //         std::cout << std::endl;
    //     }

    //     if (ortho[i] != mk::constants::missing::doubleValue && ortho[i] > 0.0)
    //     {
    //         minOrtho = std::min(minOrtho, ortho[i]);
    //     }
    // }

    // std::cout << "orthogonality: " << maxOrtho << "  " << minOrtho << std::endl;
}
