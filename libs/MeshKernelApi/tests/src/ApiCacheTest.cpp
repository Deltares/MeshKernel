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

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "MeshKernel/MeshTransformation.hpp"
#include "MeshKernel/Parameters.hpp"

#include "MeshKernelApi/BoundingBox.hpp"
#include "MeshKernelApi/GeometryList.hpp"
#include "MeshKernelApi/Mesh1D.hpp"
#include "MeshKernelApi/Mesh2D.hpp"
#include "MeshKernelApi/MeshKernel.hpp"
#include "MeshKernelApi/SplineIntersections.hpp"
#include "Version/Version.hpp"

#include "TestUtils/Definitions.hpp"
#include "TestUtils/MakeCurvilinearGrids.hpp"
#include "TestUtils/MakeMeshes.hpp"
#include "TestUtils/SampleFileReader.hpp"
#include "TestUtils/SplineReader.hpp"

#include "TestMeshGeneration.hpp"

#include <memory>
#include <numeric>

TEST(ApiCacheTest, GetHangingEdgesMesh2D_WithOneHangingEdges_GetOneHangingEdgesFailures)
{
    int meshKernelId = 0;

    int errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Prepare
    errorCode = GenerateUnstructuredMesh(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // delete an edge at the lower left corner to create a new hanging edge
    errorCode = meshkernelapi::mkernel_mesh2d_delete_edge(meshKernelId, 0.5, 0.0, 0.0, 0.0, 1.0, 1.0);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    int numHangingEdges;
    errorCode = meshkernelapi::mkernel_mesh2d_count_hanging_edges(meshKernelId, numHangingEdges);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Already cached
    errorCode = meshkernelapi::mkernel_mesh2d_count_hanging_edges(meshKernelId, numHangingEdges);
    ASSERT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    // Re-cache
    errorCode = meshkernelapi::mkernel_mesh2d_count_hanging_edges(meshKernelId, numHangingEdges);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<int> hangingEdges(numHangingEdges);
    errorCode = meshkernelapi::mkernel_mesh2d_get_hanging_edges(meshKernelId, hangingEdges.data());
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_get_hanging_edges(meshKernelId, hangingEdges.data());
    ASSERT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    errorCode = meshkernelapi::mkernel_expunge_state(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}

TEST(ApiCacheTest, GetNodesInPolygonMesh2D_OnMesh2D_NodeInPolygonFailures)
{
    int meshKernelId = 0;

    int errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = GenerateUnstructuredMesh(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // By using an empty list, all nodes will be selected
    const meshkernelapi::GeometryList geometryListIn{};

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    int numberOfNodes = -1;
    // Execute
    std::vector<int> selectedNodes(mesh2d.num_nodes);

    // No daata has been cached.
    errorCode = mkernel_mesh2d_get_nodes_in_polygons(meshKernelId, geometryListIn, 0, selectedNodes.data());
    ASSERT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    errorCode = mkernel_mesh2d_count_nodes_in_polygons(meshKernelId, geometryListIn, 1, numberOfNodes);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(numberOfNodes, mesh2d.num_nodes);

    // Values have been cached already
    errorCode = mkernel_mesh2d_count_nodes_in_polygons(meshKernelId, geometryListIn, 1, numberOfNodes);
    ASSERT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    // Re-cache data
    errorCode = mkernel_mesh2d_count_nodes_in_polygons(meshKernelId, geometryListIn, 1, numberOfNodes);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Incorrect parameters, should be the same as the call to count (which does the cahcing)
    errorCode = mkernel_mesh2d_get_nodes_in_polygons(meshKernelId, geometryListIn, 0, selectedNodes.data());
    ASSERT_EQ(meshkernel::ExitCode::ConstraintErrorCode, errorCode);

    // Re-cache data
    errorCode = mkernel_mesh2d_count_nodes_in_polygons(meshKernelId, geometryListIn, 1, numberOfNodes);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Incorrect parameters, should be the same as the call to count (which does the cahcing)
    int* nullArray = nullptr;
    errorCode = mkernel_mesh2d_get_nodes_in_polygons(meshKernelId, geometryListIn, 1, nullArray);
    ASSERT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    errorCode = meshkernelapi::mkernel_expunge_state(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}

TEST(ApiCacheTest, GetSmallFlowEdges_OnMesh2D_GetSmallFlowEdgesFailures)
{
    // Prepare a mesh with two triangles
    meshkernelapi::Mesh2D mesh2d;
    int meshKernelId = 0;

    int errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = GenerateUnstructuredMesh(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> node_x{0.0, 1.0, 1.0, 1.0};
    std::vector<double> node_y{0.0, 0.0, 0.3, -0.3};
    std::vector<int> edge_nodes{0, 3, 3, 1, 1, 0, 1, 2, 2, 0};

    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_nodes = edge_nodes.data();
    mesh2d.num_edges = static_cast<int>(edge_nodes.size() * 0.5);
    mesh2d.num_nodes = static_cast<int>(node_x.size());

    const double smallFlowEdgesThreshold = 100.0;
    meshkernelapi::GeometryList result{};

    errorCode = mkernel_mesh2d_get_small_flow_edge_centers(meshKernelId, smallFlowEdgesThreshold, result);
    ASSERT_EQ(meshkernel::ExitCode::ConstraintErrorCode, errorCode);

    // Execute
    errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    int numSmallFlowEdges;
    errorCode = meshkernelapi::mkernel_mesh2d_count_small_flow_edge_centers(meshKernelId, smallFlowEdgesThreshold, numSmallFlowEdges);

    // Assert
    std::vector<double> coordinates_x(numSmallFlowEdges);
    std::vector<double> coordinates_y(numSmallFlowEdges);
    result.coordinates_x = coordinates_x.data();
    result.coordinates_y = coordinates_y.data();
    result.num_coordinates = numSmallFlowEdges;

    errorCode = mkernel_mesh2d_get_small_flow_edge_centers(meshKernelId, smallFlowEdgesThreshold, result);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Cache has been deleted
    errorCode = mkernel_mesh2d_get_small_flow_edge_centers(meshKernelId, smallFlowEdgesThreshold, result);
    ASSERT_EQ(meshkernel::ExitCode::ConstraintErrorCode, errorCode);

    // Re-cache values
    errorCode = meshkernelapi::mkernel_mesh2d_count_small_flow_edge_centers(meshKernelId, smallFlowEdgesThreshold, numSmallFlowEdges);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Get values with different set of options
    errorCode = mkernel_mesh2d_get_small_flow_edge_centers(meshKernelId, 2.0 * smallFlowEdgesThreshold, result);
    ASSERT_EQ(meshkernel::ExitCode::ConstraintErrorCode, errorCode);

    // Re-cache values
    errorCode = meshkernelapi::mkernel_mesh2d_count_small_flow_edge_centers(meshKernelId, smallFlowEdgesThreshold, numSmallFlowEdges);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    // Attempt at getting the size again will cause an error
    errorCode = meshkernelapi::mkernel_mesh2d_count_small_flow_edge_centers(meshKernelId, smallFlowEdgesThreshold, numSmallFlowEdges);
    ASSERT_EQ(meshkernel::ExitCode::ConstraintErrorCode, errorCode);

    errorCode = meshkernelapi::mkernel_expunge_state(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}

TEST(ApiCacheTest, CountObtuseTriangles_OnMesh2DWithOneObtuseTriangle_ObtuseTrianglesFailures)
{
    // Prepare a mesh with one obtuse triangle
    int meshKernelId = 0;

    int errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d;
    std::vector<double> coordinatesX{0.0, 3.0, -1.0, 1.5};
    std::vector<double> coordinatesY{0.0, 0.0, 2.0, -2.0};
    std::vector<int> edge_nodes{0, 1, 1, 2, 2, 0, 0, 3, 3, 1};
    mesh2d.node_x = coordinatesX.data();
    mesh2d.node_y = coordinatesY.data();
    mesh2d.edge_nodes = edge_nodes.data();
    mesh2d.num_edges = static_cast<int>(edge_nodes.size() * 0.5);
    mesh2d.num_nodes = static_cast<int>(coordinatesX.size());

    // Execute
    errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Need to clear the obtuse triangle cache for the next tests
    meshkernelapi::GeometryList geometryList{};

    // Data has not yet been cached
    errorCode = mkernel_mesh2d_get_obtuse_triangles_mass_centers(meshKernelId, geometryList);
    ASSERT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    int numObtuseTriangles;
    errorCode = meshkernelapi::mkernel_mesh2d_count_obtuse_triangles(meshKernelId, numObtuseTriangles);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Already cached, cached data will be deleted
    errorCode = meshkernelapi::mkernel_mesh2d_count_obtuse_triangles(meshKernelId, numObtuseTriangles);
    ASSERT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    // Re-cache data.
    errorCode = meshkernelapi::mkernel_mesh2d_count_obtuse_triangles(meshKernelId, numObtuseTriangles);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> coordinatesObtuseTrianglesX(numObtuseTriangles);
    std::vector<double> coordinatesObtuseTrianglesY(numObtuseTriangles);
    geometryList.coordinates_x = coordinatesObtuseTrianglesX.data();
    geometryList.coordinates_y = coordinatesObtuseTrianglesY.data();
    geometryList.num_coordinates = numObtuseTriangles;
    errorCode = mkernel_mesh2d_get_obtuse_triangles_mass_centers(meshKernelId, geometryList);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Cache has been deleted in the last call
    errorCode = mkernel_mesh2d_get_obtuse_triangles_mass_centers(meshKernelId, geometryList);
    ASSERT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, errorCode);

    errorCode = meshkernelapi::mkernel_expunge_state(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}

TEST(ApiCacheTest, SplineIntersections)
{
    namespace mk = meshkernel;

    // Prepare a mesh with one obtuse triangle
    int meshKernelId = 0;

    const std::vector<mk::Point> splinePoints (LoadSplinePoints (TEST_FOLDER + "/data/CurvilinearGrids/seventy_splines.spl"));

    std::vector<double> xPoints;
    std::transform (splinePoints.begin (), splinePoints.end (), std::back_inserter (xPoints), [](const mk::Point& p){return p.x;});
    std::vector<double> yPoints;
    std::transform (splinePoints.begin (), splinePoints.end (), std::back_inserter (yPoints), [](const mk::Point& p){return p.y;});



    meshkernelapi::GeometryList splines;
    splines.num_coordinates = static_cast<int> (splinePoints.size ());
    splines.coordinates_x = xPoints.data ();
    splines.coordinates_y = yPoints.data ();

    int errorCode = meshkernelapi::mkernel_allocate_state(0, meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_initialise_spline_intersection (meshKernelId, splines);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> splineCoordX{30825.9994380052, 31503.2022933559, 33217.3720209622, 33153.884253273, 32053.4296133283, 32730.6324686789, 34085.0381793802, 34338.9892501366};
    std::vector<double> splineCoordY{377851.802873428, 379206.208584129, 381915.020005531, 383481.05160853, 384454.530713096, 385576.147942271, 385999.399726865, 387438.455794485};

    meshkernelapi::GeometryList testSpline;
    testSpline.num_coordinates = static_cast<int> (splineCoordX.size ());
    testSpline.coordinates_x = splineCoordX.data ();
    testSpline.coordinates_y = splineCoordY.data ();

    errorCode = meshkernelapi::mkernel_check_spline_intersection(meshKernelId, testSpline);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    int numberOfIntersections = -1;

    errorCode = meshkernelapi::mkernel_get_spline_intersection_dim(meshKernelId, numberOfIntersections);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    ASSERT_EQ (numberOfIntersections, 7);

    std::vector<int> splineIndices (numberOfIntersections);
    std::vector<double> splineIntersectionsAngles (numberOfIntersections);
    std::vector<double> splineIntersectionsCoordX (numberOfIntersections);
    std::vector<double> splineIntersectionsCoordY (numberOfIntersections);

    meshkernelapi::SplineIntersections splineIntersectionData;
    splineIntersectionData.num_intersections = numberOfIntersections;

    splineIntersectionData.spline_index = splineIndices.data ();
    splineIntersectionData.intersection_angle = splineIntersectionsAngles.data ();
    splineIntersectionData.intersection_x = splineIntersectionsCoordX.data ();
    splineIntersectionData.intersection_y = splineIntersectionsCoordY.data ();

    errorCode = meshkernelapi::mkernel_get_spline_intersection_data(meshKernelId, splineIntersectionData);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<int> expectedIntersectedSplines{0, 1, 2, 3, 4, 10, 23};

    std::vector<double> expectedIntersectedAngles{9.369091291645886e+01, 1.128522130510987e+02, 1.133073718558457e+02, 1.146340134566103e+02,
                                                  4.215466692286142e+01, 6.936103716736207e+01, 8.094729661593756e+01};

    std::vector<double> expectedIntersectedCoordX{3.345969673585541e+04,3.291179368046583e+04,3.199660809941513e+04,3.204744958459182e+04,
                                                  3.270633719945683e+04,3.431874691644472e+04,3.108691001430429e+04};

    std::vector<double> expectedIntersectedCoordY{3.826926563132096e+05,3.813522268790855e+05,3.799435990928682e+05,3.848083038591905e+05,
                                                  3.838236381325046e+05,3.865335746218469e+05,3.784979991789844e+05};


    constexpr double tolerance = 1.0e-5;

    for (size_t i = 0; i < splineIndices.size (); ++i)
    {
        EXPECT_EQ(splineIndices[i], expectedIntersectedSplines[i]);
        EXPECT_NEAR(splineIntersectionsAngles[i], expectedIntersectedAngles[i], tolerance);
        EXPECT_NEAR(splineIntersectionsCoordX[i], expectedIntersectedCoordX[i], tolerance);
        EXPECT_NEAR(splineIntersectionsCoordY[i], expectedIntersectedCoordY[i], tolerance);
    }


    errorCode = meshkernelapi::mkernel_finalise_spline_intersection(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}
