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

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/FlipEdges.hpp"
#include "MeshKernel/LandBoundaries.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/MeshRefinement.hpp"
#include "MeshKernel/Operations.hpp"
#include "MeshKernel/OrthogonalizationAndSmoothing.hpp"
#include "MeshKernel/Polygons.hpp"
#include "MeshKernel/UndoActions/UndoAction.hpp"
#include "MeshKernel/Utilities/Utilities.hpp"
#include "TestUtils/Definitions.hpp"
#include "TestUtils/MakeMeshes.hpp"

using namespace meshkernel;

TEST(OrthogonalizationAndSmoothing, OrthogonalizationOneQuadOneTriangle)
{
    // Preparation
    std::vector<Point> nodes{{0.0, 0.0},
                             {0.0, 10.0},
                             {10.0, 0.0},
                             {10.0, 10.0},
                             {20.0, 0.0}};

    std::vector<Edge> edges{{1, 0},
                            {0, 2},
                            {2, 4},
                            {4, 3},
                            {3, 2},
                            {3, 1}};

    const auto projectToLandBoundaryOption = LandBoundaries::ProjectToLandBoundaryOption::DoNotProjectToLandBoundary;
    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.inner_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.outer_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor_at_boundary = 0.975;
    orthogonalizationParameters.areal_to_angle_smoothing_factor = 1.0;

    // Execute
    auto mesh = Mesh2D(edges, nodes, Projection::cartesian);
    auto polygon = std::make_unique<Polygons>();

    std::vector<Point> landBoundary{};
    auto landboundaries = std::make_unique<LandBoundaries>(landBoundary, mesh, *polygon);

    OrthogonalizationAndSmoothing orthogonalization(mesh,
                                                    std::move(polygon),
                                                    std::move(landboundaries),
                                                    projectToLandBoundaryOption,
                                                    orthogonalizationParameters);

    [[maybe_unused]] auto undoAction = orthogonalization.Initialize();
    orthogonalization.Compute();

    // Assert
    constexpr double tolerance = 1e-8;
    ASSERT_NEAR(0.0, mesh.Node(0).x, tolerance);
    ASSERT_NEAR(0.0, mesh.Node(1).x, tolerance);
    ASSERT_NEAR(10.0, mesh.Node(2).x, tolerance);
    ASSERT_NEAR(10.0, mesh.Node(3).x, tolerance);
    ASSERT_NEAR(20.0, mesh.Node(4).x, tolerance);

    ASSERT_NEAR(0.0, mesh.Node(0).y, tolerance);
    ASSERT_NEAR(10.0, mesh.Node(1).y, tolerance);
    ASSERT_NEAR(0.0, mesh.Node(2).y, tolerance);
    ASSERT_NEAR(10.0, mesh.Node(3).y, tolerance);
    ASSERT_NEAR(0.0, mesh.Node(4).y, tolerance);
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizationSmallTriangularGrid)
{

    // now build node-edge mapping
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/SmallTriangularGrid_net.nc");
    const auto projectToLandBoundaryOption = LandBoundaries::ProjectToLandBoundaryOption::DoNotProjectToLandBoundary;
    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor_at_boundary = 1.0;
    orthogonalizationParameters.areal_to_angle_smoothing_factor = 1.0;

    auto polygon = std::make_unique<Polygons>();

    std::vector<Point> landBoundary;
    auto landboundaries = std::make_unique<LandBoundaries>(landBoundary, *mesh, *polygon);

    OrthogonalizationAndSmoothing orthogonalization(*mesh,
                                                    std::move(polygon),
                                                    std::move(landboundaries),
                                                    projectToLandBoundaryOption,
                                                    orthogonalizationParameters);

    [[maybe_unused]] auto undoAction = orthogonalization.Initialize();
    orthogonalization.Compute();

    constexpr double tolerance = 0.8e-2;

    ASSERT_NEAR(325.590101919525, mesh->Node(0).x, tolerance);
    ASSERT_NEAR(229.213730481198, mesh->Node(1).x, tolerance);
    ASSERT_NEAR(263.439319753147, mesh->Node(2).x, tolerance);
    ASSERT_NEAR(429.191105834504, mesh->Node(3).x, tolerance);
    ASSERT_NEAR(535.865215426468, mesh->Node(4).x, tolerance);
    ASSERT_NEAR(503.753784179688, mesh->Node(5).x, tolerance);
    ASSERT_NEAR(354.048340705929, mesh->Node(6).x, tolerance);
    ASSERT_NEAR(346.790050854504, mesh->Node(7).x, tolerance);
    ASSERT_NEAR(315.030130405285, mesh->Node(8).x, tolerance);
    ASSERT_NEAR(424.314957449766, mesh->Node(9).x, tolerance);

    ASSERT_NEAR(455.319334078551, mesh->Node(0).y, tolerance);
    ASSERT_NEAR(362.573521507281, mesh->Node(1).y, tolerance);
    ASSERT_NEAR(241.096458631763, mesh->Node(2).y, tolerance);
    ASSERT_NEAR(211.483073921775, mesh->Node(3).y, tolerance);
    ASSERT_NEAR(311.401495506714, mesh->Node(4).y, tolerance);
    ASSERT_NEAR(432.379974365234, mesh->Node(5).y, tolerance);
    ASSERT_NEAR(458.064836627594, mesh->Node(6).y, tolerance);
    ASSERT_NEAR(405.311585650679, mesh->Node(7).y, tolerance);
    ASSERT_NEAR(319.612138503550, mesh->Node(8).y, tolerance);
    ASSERT_NEAR(327.102805172725, mesh->Node(9).y, tolerance);
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizationSmallTriangularGridAsNcFile)
{

    // now build node-edge mapping
    auto mesh = MakeSmallSizeTriangularMeshForTestingAsNcFile();

    const auto projectToLandBoundaryOption = LandBoundaries::ProjectToLandBoundaryOption::DoNotProjectToLandBoundary;
    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor_at_boundary = 1.0;
    orthogonalizationParameters.areal_to_angle_smoothing_factor = 1.0;

    auto polygon = std::make_unique<Polygons>();

    std::vector<Point> landBoundary;
    auto landboundaries = std::make_unique<LandBoundaries>(landBoundary, *mesh, *polygon);

    OrthogonalizationAndSmoothing orthogonalization(*mesh,
                                                    std::move(polygon),
                                                    std::move(landboundaries),
                                                    projectToLandBoundaryOption,
                                                    orthogonalizationParameters);

    [[maybe_unused]] auto undoAction = orthogonalization.Initialize();
    orthogonalization.Compute();

    constexpr double tolerance = 1e-2;

    ASSERT_NEAR(325.590101919525, mesh->Node(0).x, tolerance);
    ASSERT_NEAR(229.213730481198, mesh->Node(1).x, tolerance);
    ASSERT_NEAR(263.439319753147, mesh->Node(2).x, tolerance);
    ASSERT_NEAR(429.191105834504, mesh->Node(3).x, tolerance);
    ASSERT_NEAR(535.865215426468, mesh->Node(4).x, tolerance);
    ASSERT_NEAR(503.753784179688, mesh->Node(5).x, tolerance);
    ASSERT_NEAR(354.048340705929, mesh->Node(6).x, tolerance);
    ASSERT_NEAR(346.790050854504, mesh->Node(7).x, tolerance);
    ASSERT_NEAR(315.030130405285, mesh->Node(8).x, tolerance);
    ASSERT_NEAR(424.314957449766, mesh->Node(9).x, tolerance);

    ASSERT_NEAR(455.319334078551, mesh->Node(0).y, tolerance);
    ASSERT_NEAR(362.573521507281, mesh->Node(1).y, tolerance);
    ASSERT_NEAR(241.096458631763, mesh->Node(2).y, tolerance);
    ASSERT_NEAR(211.483073921775, mesh->Node(3).y, tolerance);
    ASSERT_NEAR(311.401495506714, mesh->Node(4).y, tolerance);
    ASSERT_NEAR(432.379974365234, mesh->Node(5).y, tolerance);
    ASSERT_NEAR(458.064836627594, mesh->Node(6).y, tolerance);
    ASSERT_NEAR(405.311585650679, mesh->Node(7).y, tolerance);
    ASSERT_NEAR(319.612138503550, mesh->Node(8).y, tolerance);
    ASSERT_NEAR(327.102805172725, mesh->Node(9).y, tolerance);
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizationMediumTriangularGridWithPolygon)
{
    // now build node-edge mapping
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/SmallTriangularGrid_net.nc");

    // The original mesh nodes
    const std::vector<meshkernel::Point> meshNodes(mesh->Nodes());

    const auto projectToLandBoundaryOption = LandBoundaries::ProjectToLandBoundaryOption::DoNotProjectToLandBoundary;
    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor_at_boundary = 1.0;
    orthogonalizationParameters.areal_to_angle_smoothing_factor = 1.0;

    std::vector<Point> nodes{{342.987518, 471.121002},
                             {327.640900, 380.846436},
                             {396.851135, 201.200073},
                             {514.207581, 203.607407},
                             {569.274841, 294.483765},
                             {568.673035, 379.943695},
                             {515.712158, 458.783478},
                             {343.288422, 471.722809},
                             {342.987518, 471.121002}};

    std::vector<Point> landBoundary;
    auto polygon = std::make_unique<Polygons>(nodes, Projection::cartesian);
    auto landboundaries = std::make_unique<LandBoundaries>(landBoundary, *mesh, *polygon);

    OrthogonalizationAndSmoothing orthogonalization(*mesh,
                                                    std::move(polygon),
                                                    std::move(landboundaries),
                                                    projectToLandBoundaryOption,
                                                    orthogonalizationParameters);

    auto undoAction = orthogonalization.Initialize();
    orthogonalization.Compute();

    constexpr double tolerance = 1.8;
    // check the first 10 points

    ASSERT_NEAR(322.252624511719, mesh->Node(0).x, tolerance);
    ASSERT_NEAR(227.002044677734, mesh->Node(1).x, tolerance);
    ASSERT_NEAR(259.252227783203, mesh->Node(2).x, tolerance);
    ASSERT_NEAR(426.680730020925, mesh->Node(3).x, tolerance);
    ASSERT_NEAR(535.742871000462, mesh->Node(4).x, tolerance);
    ASSERT_NEAR(502.614684581472, mesh->Node(5).x, tolerance);
    ASSERT_NEAR(350.424583126127, mesh->Node(6).x, tolerance);
    ASSERT_NEAR(342.832067902036, mesh->Node(7).x, tolerance);
    ASSERT_NEAR(310.300984548069, mesh->Node(8).x, tolerance);
    ASSERT_NEAR(421.791577149564, mesh->Node(9).x, tolerance);

    ASSERT_NEAR(454.880187988281, mesh->Node(0).y, tolerance);
    ASSERT_NEAR(360.379241943359, mesh->Node(1).y, tolerance);
    ASSERT_NEAR(241.878051757812, mesh->Node(2).y, tolerance);
    ASSERT_NEAR(210.624626374983, mesh->Node(3).y, tolerance);
    ASSERT_NEAR(311.862423032534, mesh->Node(4).y, tolerance);
    ASSERT_NEAR(432.575408917275, mesh->Node(5).y, tolerance);
    ASSERT_NEAR(458.587061164954, mesh->Node(6).y, tolerance);
    ASSERT_NEAR(405.390188082498, mesh->Node(7).y, tolerance);
    ASSERT_NEAR(319.410057398020, mesh->Node(8).y, tolerance);
    ASSERT_NEAR(327.001109057344, mesh->Node(9).y, tolerance);

    undoAction->Restore();

    constexpr double nodeTolerance = 1.0e-10;

    for (meshkernel::UInt i = 0; i < meshNodes.size(); ++i)
    {
        EXPECT_NEAR(meshNodes[i].x, mesh->Node(i).x, nodeTolerance);
        EXPECT_NEAR(meshNodes[i].y, mesh->Node(i).y, nodeTolerance);
    }
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizationMediumTriangularGridWithUndo)
{
    // now build node-edge mapping
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/TestOrthogonalizationMediumTriangularGrid_net.nc");

    const auto projectToLandBoundaryOption = LandBoundaries::ProjectToLandBoundaryOption::DoNotProjectToLandBoundary;
    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor_at_boundary = 0.5;
    orthogonalizationParameters.areal_to_angle_smoothing_factor = 1.0;

    auto polygon = std::make_unique<Polygons>();
    std::vector<Point> landBoundary;
    auto landBoundaries = std::make_unique<LandBoundaries>(landBoundary, *mesh, *polygon);

    const std::vector<meshkernel::Point> meshNodes = mesh->Nodes();

    OrthogonalizationAndSmoothing orthogonalization(*mesh,
                                                    std::move(polygon),
                                                    std::move(landBoundaries),
                                                    projectToLandBoundaryOption,
                                                    orthogonalizationParameters);

    std::unique_ptr<meshkernel::UndoAction> undoAction = orthogonalization.Initialize();
    orthogonalization.Compute();

    constexpr double tolerance = 1.2;

    // check the first 10 points
    ASSERT_NEAR(68.771705432835475, mesh->Node(0).x, tolerance);
    ASSERT_NEAR(169.49338272334273, mesh->Node(1).x, tolerance);
    ASSERT_NEAR(262.80128484924921, mesh->Node(2).x, tolerance);
    ASSERT_NEAR(361.60010033352023, mesh->Node(3).x, tolerance);
    ASSERT_NEAR(468.13991812406925, mesh->Node(4).x, tolerance);
    ASSERT_NEAR(549.89461192844624, mesh->Node(5).x, tolerance);
    ASSERT_NEAR(653.02704974527421, mesh->Node(6).x, tolerance);
    ASSERT_NEAR(747.81537706979441, mesh->Node(7).x, tolerance);
    ASSERT_NEAR(853.40641427112951, mesh->Node(8).x, tolerance);
    ASSERT_NEAR(938.69752431820143, mesh->Node(9).x, tolerance);

    ASSERT_NEAR(1399.7751472360221, mesh->Node(0).y, tolerance);
    ASSERT_NEAR(1426.5945287630802, mesh->Node(1).y, tolerance);
    ASSERT_NEAR(1451.4398281457179, mesh->Node(2).y, tolerance);
    ASSERT_NEAR(1477.7472050498141, mesh->Node(3).y, tolerance);
    ASSERT_NEAR(1506.1157955857589, mesh->Node(4).y, tolerance);
    ASSERT_NEAR(1527.8847968946166, mesh->Node(5).y, tolerance);
    ASSERT_NEAR(1555.3460969050145, mesh->Node(6).y, tolerance);
    ASSERT_NEAR(1580.5855923464549, mesh->Node(7).y, tolerance);
    ASSERT_NEAR(1608.7015489976982, mesh->Node(8).y, tolerance);
    ASSERT_NEAR(1631.412199601948, mesh->Node(9).y, tolerance);

    undoAction->Restore();

    EXPECT_EQ(meshNodes.size(), mesh->GetNumNodes());

    for (meshkernel::UInt i = 0; i < meshNodes.size(); ++i)
    {
        EXPECT_EQ(meshNodes[i].x, mesh->Node(i).x);
        EXPECT_EQ(meshNodes[i].y, mesh->Node(i).y);
    }
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizationFourQuads)
{
    auto mesh = MakeRectangularMeshForTesting(3, 3, 1.0, Projection::cartesian);

    const auto projectToLandBoundaryOption = LandBoundaries::ProjectToLandBoundaryOption::DoNotProjectToLandBoundary;
    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.inner_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.outer_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor_at_boundary = 0.975;
    orthogonalizationParameters.areal_to_angle_smoothing_factor = 1.0;

    auto polygon = std::make_unique<Polygons>();
    std::vector<Point> landBoundary{};
    auto landboundaries = std::make_unique<LandBoundaries>(landBoundary, *mesh, *polygon);

    OrthogonalizationAndSmoothing orthogonalization(*mesh,
                                                    std::move(polygon),
                                                    std::move(landboundaries),
                                                    projectToLandBoundaryOption,
                                                    orthogonalizationParameters);

    [[maybe_unused]] auto undoAction = orthogonalization.Initialize();
    orthogonalization.Compute();
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizeAndSnapToLandBoundaries)
{
    using namespace constants;

    // Prepare
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/SmallTriangularGrid_net.nc");

    // the land boundary to use
    std::vector<Point> landBoundary{{235.561218, 290.571899},
                                    {265.953522, 436.515747},
                                    {429.349854, 450.959656},
                                    {535.271545, 386.262909},
                                    {missing::doubleValue, missing::doubleValue},
                                    {246.995941, 262.285858},
                                    {351.112183, 237.309906},
                                    {443.191895, 262.285858},
                                    {553.627319, 327.283539},
                                    {missing::doubleValue, missing::doubleValue}};

    // snap to land boundaries
    const auto projectToLandBoundaryOption = LandBoundaries::ProjectToLandBoundaryOption::OuterMeshBoundaryToLandBoundary;
    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor_at_boundary = 0.975;
    orthogonalizationParameters.areal_to_angle_smoothing_factor = 1.0;

    // no enclosing polygon
    auto polygon = std::make_unique<Polygons>();
    auto landboundaries = std::make_unique<LandBoundaries>(landBoundary, *mesh, *polygon);

    OrthogonalizationAndSmoothing orthogonalization(*mesh,
                                                    std::move(polygon),
                                                    std::move(landboundaries),
                                                    projectToLandBoundaryOption,
                                                    orthogonalizationParameters);

    [[maybe_unused]] auto undoAction = orthogonalization.Initialize();
    orthogonalization.Compute();

    // check the values
    constexpr double tolerance = 0.15;
    ASSERT_NEAR(313.081472564480, mesh->Node(0).x, tolerance);
    ASSERT_NEAR(253.641466857330, mesh->Node(1).x, tolerance);
    ASSERT_NEAR(254.777224294204, mesh->Node(2).x, tolerance);
    ASSERT_NEAR(443.191895000000, mesh->Node(3).x, tolerance);
    ASSERT_NEAR(535.240231516760, mesh->Node(4).x, tolerance);
    ASSERT_NEAR(480.436129612752, mesh->Node(5).x, tolerance);
    ASSERT_NEAR(345.948240805397, mesh->Node(6).x, tolerance);
    ASSERT_NEAR(342.668434889472, mesh->Node(7).x, tolerance);
    ASSERT_NEAR(318.414413615199, mesh->Node(8).x, tolerance);
    ASSERT_NEAR(424.616311031376, mesh->Node(9).x, tolerance);

    ASSERT_NEAR(440.681763586650, mesh->Node(0).y, tolerance);
    ASSERT_NEAR(377.393256506700, mesh->Node(1).y, tolerance);
    ASSERT_NEAR(260.419242817573, mesh->Node(2).y, tolerance);
    ASSERT_NEAR(262.285858000000, mesh->Node(3).y, tolerance);
    ASSERT_NEAR(316.461666783032, mesh->Node(4).y, tolerance);
    ASSERT_NEAR(419.756265860671, mesh->Node(5).y, tolerance);
    ASSERT_NEAR(443.587120174434, mesh->Node(6).y, tolerance);
    ASSERT_NEAR(402.913858250569, mesh->Node(7).y, tolerance);
    ASSERT_NEAR(336.831643075189, mesh->Node(8).y, tolerance);
    ASSERT_NEAR(340.875100904741, mesh->Node(9).y, tolerance);
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizationSphericalRectangular)
{
    // 1 Setup
    auto mesh = MakeRectangularMeshForTesting(4, 4, 0.003, Projection::spherical, {41.1, 41.1});

    const auto projectToLandBoundaryOption = LandBoundaries::ProjectToLandBoundaryOption::DoNotProjectToLandBoundary;
    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor_at_boundary = 1.0;
    orthogonalizationParameters.areal_to_angle_smoothing_factor = 1.0;

    auto polygon = std::make_unique<Polygons>();
    std::vector<Point> landBoundary{};
    auto landboundaries = std::make_unique<LandBoundaries>(landBoundary, *mesh, *polygon);

    OrthogonalizationAndSmoothing orthogonalization(*mesh,
                                                    std::move(polygon),
                                                    std::move(landboundaries),
                                                    projectToLandBoundaryOption,
                                                    orthogonalizationParameters);

    [[maybe_unused]] auto undoAction = orthogonalization.Initialize();
    orthogonalization.Compute();

    // check the values
    constexpr double tolerance = 1e-6;
    ASSERT_NEAR(41.1, mesh->Node(0).x, tolerance);
    ASSERT_NEAR(41.1, mesh->Node(1).x, tolerance);
    ASSERT_NEAR(41.1, mesh->Node(2).x, tolerance);
    ASSERT_NEAR(41.1, mesh->Node(3).x, tolerance);

    ASSERT_NEAR(41.103, mesh->Node(4).x, tolerance);
    ASSERT_NEAR(41.103, mesh->Node(5).x, tolerance);
    ASSERT_NEAR(41.103, mesh->Node(6).x, tolerance);
    ASSERT_NEAR(41.103, mesh->Node(7).x, tolerance);

    ASSERT_NEAR(41.106, mesh->Node(8).x, tolerance);
    ASSERT_NEAR(41.106, mesh->Node(9).x, tolerance);
    ASSERT_NEAR(41.106, mesh->Node(10).x, tolerance);
    ASSERT_NEAR(41.106, mesh->Node(11).x, tolerance);

    ASSERT_NEAR(41.109, mesh->Node(12).x, tolerance);
    ASSERT_NEAR(41.109, mesh->Node(13).x, tolerance);
    ASSERT_NEAR(41.109, mesh->Node(14).x, tolerance);
    ASSERT_NEAR(41.109, mesh->Node(15).x, tolerance);

    ASSERT_NEAR(41.1, mesh->Node(0).y, tolerance);
    ASSERT_NEAR(41.103, mesh->Node(1).y, tolerance);
    ASSERT_NEAR(41.106, mesh->Node(2).y, tolerance);
    ASSERT_NEAR(41.109, mesh->Node(3).y, tolerance);

    ASSERT_NEAR(41.1, mesh->Node(4).y, tolerance);
    ASSERT_NEAR(41.103, mesh->Node(5).y, tolerance);
    ASSERT_NEAR(41.106, mesh->Node(6).y, tolerance);
    ASSERT_NEAR(41.109, mesh->Node(7).y, tolerance);

    ASSERT_NEAR(41.1, mesh->Node(8).y, tolerance);
    ASSERT_NEAR(41.103, mesh->Node(9).y, tolerance);
    ASSERT_NEAR(41.106, mesh->Node(10).y, tolerance);
    ASSERT_NEAR(41.109, mesh->Node(11).y, tolerance);

    ASSERT_NEAR(41.1, mesh->Node(12).y, tolerance);
    ASSERT_NEAR(41.103, mesh->Node(13).y, tolerance);
    ASSERT_NEAR(41.106, mesh->Node(14).y, tolerance);
    ASSERT_NEAR(41.109, mesh->Node(15).y, tolerance);
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizationSmallTriangulargridSpherical)
{
    std::vector<Point> nodes{{41.1019592, 41.1072273},
                             {41.1044655, 41.1043587},
                             {41.1051979, 41.1073151},
                             {41.1080132, 41.1046638},
                             {41.1014137, 41.1039963}};

    std::vector<Edge> edges{{4, 0},
                            {4, 1},
                            {1, 0},
                            {1, 3},
                            {3, 2},
                            {2, 1},
                            {0, 2}};

    auto mesh = std::make_shared<Mesh2D>(edges, nodes, Projection::spherical);

    const auto projectToLandBoundaryOption = LandBoundaries::ProjectToLandBoundaryOption::DoNotProjectToLandBoundary;
    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor_at_boundary = 1.0;
    orthogonalizationParameters.areal_to_angle_smoothing_factor = 1.0;

    // no enclosing polygon
    auto polygon = std::make_unique<Polygons>();
    std::vector<Point> landBoundary{};
    auto landboundaries = std::make_unique<LandBoundaries>(landBoundary, *mesh, *polygon);

    OrthogonalizationAndSmoothing orthogonalization(*mesh,
                                                    std::move(polygon),
                                                    std::move(landboundaries),
                                                    projectToLandBoundaryOption,
                                                    orthogonalizationParameters);

    [[maybe_unused]] auto undoAction = orthogonalization.Initialize();
    orthogonalization.Compute();

    constexpr double tolerance = 1e-3;

    ASSERT_NEAR(41.1019592285156, mesh->Node(0).x, tolerance);
    ASSERT_NEAR(41.1044654597059, mesh->Node(1).x, tolerance);
    ASSERT_NEAR(41.1051978230878, mesh->Node(2).x, tolerance);
    ASSERT_NEAR(41.1080131530762, mesh->Node(3).x, tolerance);
    ASSERT_NEAR(41.1014137268066, mesh->Node(4).x, tolerance);

    ASSERT_NEAR(41.1072273254395, mesh->Node(0).y, tolerance);
    ASSERT_NEAR(41.1043586701373, mesh->Node(1).y, tolerance);
    ASSERT_NEAR(41.1073150612170, mesh->Node(2).y, tolerance);
    ASSERT_NEAR(41.1046638488770, mesh->Node(3).y, tolerance);
    ASSERT_NEAR(41.1039962768555, mesh->Node(4).y, tolerance);
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizationSmallTriangulargridSphericalWithUndo)
{
    std::vector<Point> nodes{{41.1019592, 41.1072273},
                             {41.1044655, 41.1043587},
                             {41.1051979, 41.1073151},
                             {41.1080132, 41.1046638},
                             {41.1014137, 41.1039963}};

    std::vector<Edge> edges{{4, 0},
                            {4, 1},
                            {1, 0},
                            {1, 3},
                            {3, 2},
                            {2, 1},
                            {0, 2}};

    auto mesh = std::make_shared<Mesh2D>(edges, nodes, Projection::spherical);

    const auto projectToLandBoundaryOption = LandBoundaries::ProjectToLandBoundaryOption::DoNotProjectToLandBoundary;
    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor_at_boundary = 1.0;
    orthogonalizationParameters.areal_to_angle_smoothing_factor = 1.0;

    // no enclosing polygon
    auto polygon = std::make_unique<Polygons>();
    std::vector<Point> landBoundary{};
    auto landboundaries = std::make_unique<LandBoundaries>(landBoundary, *mesh, *polygon);

    const std::vector<meshkernel::Point> meshNodes = mesh->Nodes();

    OrthogonalizationAndSmoothing orthogonalization(*mesh,
                                                    std::move(polygon),
                                                    std::move(landboundaries),
                                                    projectToLandBoundaryOption,
                                                    orthogonalizationParameters);

    std::unique_ptr<meshkernel::UndoAction> undoAction = orthogonalization.Initialize();
    orthogonalization.Compute();

    constexpr double tolerance = 1e-3;

    ASSERT_NEAR(41.1019592285156, mesh->Node(0).x, tolerance);
    ASSERT_NEAR(41.1044654597059, mesh->Node(1).x, tolerance);
    ASSERT_NEAR(41.1051978230878, mesh->Node(2).x, tolerance);
    ASSERT_NEAR(41.1080131530762, mesh->Node(3).x, tolerance);
    ASSERT_NEAR(41.1014137268066, mesh->Node(4).x, tolerance);

    ASSERT_NEAR(41.1072273254395, mesh->Node(0).y, tolerance);
    ASSERT_NEAR(41.1043586701373, mesh->Node(1).y, tolerance);
    ASSERT_NEAR(41.1073150612170, mesh->Node(2).y, tolerance);
    ASSERT_NEAR(41.1046638488770, mesh->Node(3).y, tolerance);
    ASSERT_NEAR(41.1039962768555, mesh->Node(4).y, tolerance);

    undoAction->Restore();

    EXPECT_EQ(meshNodes.size(), mesh->GetNumNodes());

    for (meshkernel::UInt i = 0; i < meshNodes.size(); ++i)
    {
        EXPECT_EQ(meshNodes[i].x, mesh->Node(i).x);
        EXPECT_EQ(meshNodes[i].y, mesh->Node(i).y);
    }
}

TEST(OrthogonalizationAndSmoothing, RefineUndoThenOrthogonalise)
{
    auto mesh = MakeRectangularMeshForTestingRand(20, 20, 1.0, Projection::cartesian);

    const std::vector<meshkernel::Point> originalPoints = mesh->Nodes();

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

    auto refinementUndoAction = meshRefinement.Compute();
    refinementUndoAction->Restore();

    const auto projectToLandBoundaryOption = LandBoundaries::ProjectToLandBoundaryOption::DoNotProjectToLandBoundary;
    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.inner_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.outer_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor_at_boundary = 0.975;
    orthogonalizationParameters.areal_to_angle_smoothing_factor = 1.0;

    std::vector<meshkernel::Point> polygonNodes{{10.5, 10.5}, {30.5, 10.5}, {30.5, 30.5}, {10.5, 30.5}, {10.5, 10.5}};

    auto polygon = std::make_unique<Polygons>(polygonNodes, Projection::cartesian);
    auto orthogonalisationPolygon = std::make_unique<Polygons>(polygonNodes, Projection::cartesian);
    std::vector<Point> landBoundary{};
    auto landboundaries = std::make_unique<LandBoundaries>(landBoundary, *mesh, *polygon);

    OrthogonalizationAndSmoothing orthogonalization(*mesh,
                                                    std::move(orthogonalisationPolygon),
                                                    std::move(landboundaries),
                                                    projectToLandBoundaryOption,
                                                    orthogonalizationParameters);

    [[maybe_unused]] auto undoAction = orthogonalization.Initialize();

    const std::vector<meshkernel::Point>& nodes = mesh->Nodes();
    const auto& enclosure = polygon->Enclosure(0);

    constexpr double tolerance = 1.0e-12;

    orthogonalization.Compute();

    size_t count = 0;
    size_t orignialCount = 0;

    for (size_t i = 0; i < nodes.size(); ++i)
    {
        if (nodes[i].IsValid())
        {
            if (!enclosure.Contains(nodes[i]))
            {
                EXPECT_NEAR(originalPoints[count].x, nodes[i].x, tolerance) << "x-position: " << i << "  " << count << "  " << orignialCount;
                EXPECT_NEAR(originalPoints[count].y, nodes[i].y, tolerance) << "y-position: " << i << "  " << count << "  " << orignialCount;
                ++orignialCount;
            }
            ++count;
        }
    }

    // There are 400 points in the original mesh
    EXPECT_EQ(count, 400);

    // Only check the number of points not subject to orthogonalisation
    EXPECT_EQ(orignialCount, 319);
}

TEST(OrthogonalizationAndSmoothing, OrthogonalizationWithGapsInNodeAndEdgeLists)
{

    std::vector<Point> polygonPoints{{0, 0}, {25, 5}, {50, 10}, {75, 15}, {100, 20}, {125, 35}, {150, 50}, {125, 58.33333333333334}, {100, 66.66666666666667}, {75, 75}, {51.25, 63.75}, {27.5, 52.5}, {3.75, 41.25}, {-20, 30}, {0, 0}};

    std::unique_ptr<Polygons> polygon = std::make_unique<Polygons>(polygonPoints, Projection::cartesian);

    // generate samples in all polygons
    const std::vector<std::vector<Point>> generatedPoints = polygon->ComputePointsInPolygons();

    meshkernel::Mesh2D mesh(generatedPoints[0], *polygon, Projection::cartesian);

    // Create some gaps in the node and edge arrays
    auto [nodeId, nodeInsertUndo] = mesh.InsertNode({0.5, 0.5});
    auto originNodeId = mesh.FindNodeCloseToAPoint({0.0, 0.0}, 0.1);
    [[maybe_unused]] auto [edgeId, edgeInsertUndo] = mesh.ConnectNodes(nodeId, originNodeId);
    [[maybe_unused]] auto nodeRemovaUndo = mesh.DeleteNode(nodeId);

    std::unique_ptr<LandBoundaries> boundary = std::make_unique<LandBoundaries>(polygonPoints, mesh, *polygon);

    FlipEdges flipEdges1(mesh, *boundary, true, false);

    const auto projectToLandBoundaryOption = LandBoundaries::ProjectToLandBoundaryOption::DoNotProjectToLandBoundary;
    OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor_at_boundary = 0.975;
    orthogonalizationParameters.areal_to_angle_smoothing_factor = 1.0;

    FlipEdges flipEdges(mesh, *boundary, true, false);

    OrthogonalizationAndSmoothing orthogonalization(mesh,
                                                    std::move(polygon),
                                                    std::move(boundary),
                                                    projectToLandBoundaryOption,
                                                    orthogonalizationParameters);

    const std::vector<double> expectedX{24.9835753034083, 45.7888024937107, 70.0728593318696,
                                        96.4736020240226, 123.464866047009, 150.0,
                                        123.690307803787, 98.5465046861055, 74.1522156841494,
                                        49.2936207530967, 26.9588231383586, 5.10554425586482,
                                        -20.0, 2.53218892861755, 12.2603085772356,
                                        30.0247966898692, 83.3308053183758, 104.729358926806,
                                        43.9223120527096, 62.551166421051, 57.2356630362903,
                                        87.5391262898118, 41.1303589083358, -999.0};

    const std::vector<double> expectedY{4.99671506068165, 9.15776049874214, 14.0145718663739,
                                        19.2947204048045, 34.0789196282053, 50.0,
                                        58.7698973987377, 67.1511651046315, 74.5984179556497,
                                        62.8232940409406, 52.2436530655383, 41.8920999106728,
                                        30.0, 0.506437785723511, 19.9059354533534,
                                        33.9652316161109, 49.6795938035325, 45.9248709417029,
                                        43.1509423611897, 47.6054287440074, 28.1281033458494,
                                        34.0734735937222, 22.8901803546926, -999.0};

    const std::vector<meshkernel::UInt> edgeFirst{13, 14, 12, 13, 0, 22, 0, 0, 20, 18, 22, 15,
                                                  22, 1, 14, 11, 10, 11, 15, 18, 9, 10, 18, 19,
                                                  18, 1, 20, 21, 21, 1, 2, 21, 16, 19, 21, 3, 3,
                                                  2, 21, 17, 4, 16, 8, 16, 7, 8, 4, 6, 4, 5, 6, 7,
                                                  meshkernel::constants::missing::uintValue};

    const std::vector<meshkernel::UInt> edgeSecond{14, 12, 13, 0, 14, 14, 1, 22, 18, 22, 20, 14,
                                                   15, 22, 11, 12, 11, 15, 10, 9, 10, 18, 15, 9,
                                                   19, 20, 2, 20, 19, 2, 21, 16, 19, 20, 3, 4, 17,
                                                   3, 17, 16, 17, 8, 19, 7, 8, 9, 6, 17, 5, 6, 7, 17,
                                                   meshkernel::constants::missing::uintValue};

    [[maybe_unused]] auto flipUndoAction = flipEdges.Compute();
    [[maybe_unused]] auto orthogUndoAction = orthogonalization.Initialize();

    orthogonalization.Compute();

    const double tolerance = 1.0e-8;

    ASSERT_EQ(static_cast<size_t>(mesh.GetNumNodes()), expectedX.size());
    ASSERT_EQ(static_cast<size_t>(mesh.GetNumEdges()), edgeFirst.size());

    for (meshkernel::UInt i = 0; i < mesh.GetNumNodes(); ++i)
    {
        EXPECT_NEAR(expectedX[i], mesh.Node(i).x, tolerance);
        EXPECT_NEAR(expectedY[i], mesh.Node(i).y, tolerance);
    }

    for (meshkernel::UInt i = 0; i < mesh.GetNumEdges(); ++i)
    {
        EXPECT_EQ(edgeFirst[i], mesh.GetEdge(i).first);
        EXPECT_EQ(edgeSecond[i], mesh.GetEdge(i).second);
    }
}
