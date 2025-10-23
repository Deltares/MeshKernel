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

#include <filesystem>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <MeshKernel/Parameters.hpp>

#include <MeshKernelApi/GeometryList.hpp>
#include <MeshKernelApi/Mesh2D.hpp>
#include <MeshKernelApi/MeshKernel.hpp>

#include <TestUtils/Definitions.hpp>
#include <TestUtils/MakeCurvilinearGrids.hpp>
#include <TestUtils/MakeMeshes.hpp>
#include <TestUtils/SampleFileReader.hpp>

#include "CartesianApiTestFixture.hpp"
#include "MeshKernel/AveragingInterpolation.hpp"
#include "MeshKernel/MeshRefinement.hpp"
#include "SampleGenerator.hpp"

#include "MeshKernel/Utilities/Utilities.hpp"

TEST(InvalidCellsPolygonsTests, OMG)
{
    // Prepare
    int meshKernelId;
    constexpr int isSpherical = 0;
    meshkernelapi::mkernel_allocate_state(isSpherical, meshKernelId);

    meshkernel::MakeGridParameters gridParameters;
    gridParameters.num_columns = 10;
    gridParameters.num_rows = 10;
    gridParameters.block_size_x = 20.0;
    gridParameters.block_size_y = 20.0;
    gridParameters.origin_x = 0.0;
    gridParameters.origin_y = 0.0;
    gridParameters.angle = 0.0;

    auto errorCode = meshkernelapi::mkernel_mesh2d_make_rectangular_mesh(meshKernelId, gridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    //--------------------------------

    std::vector<double> patchX{45.0, 155.0, 155.0, 45.0, 45.0, meshkernel::constants::missing::innerOuterSeparator, 65.0, 115.0, 115.0, 65.0, 65.0};
    std::vector<double> patchY{5.0, 5.0, 155.0, 155.0, 5.0, meshkernel::constants::missing::innerOuterSeparator, 65.0, 65.0, 115.0, 115.0, 65.0};

    meshkernelapi::GeometryList refinementPolygon;
    refinementPolygon.num_coordinates = static_cast<int>(patchX.size());
    refinementPolygon.geometry_separator = meshkernel::constants::missing::doubleValue;
    refinementPolygon.coordinates_x = patchX.data();
    refinementPolygon.coordinates_y = patchY.data();

    errorCode = meshkernelapi::mkernel_mesh2d_casulli_refinement_on_polygon(meshKernelId, refinementPolygon);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    //--------------------------------

    std::vector<double> invalidCellsX{80.0, 85.0, 75.0, 80.0, meshkernel::constants::missing::doubleValue,
                                      180.0, 200.0, 200.0, 180.0, 180.0, meshkernel::constants::missing::doubleValue,
                                      125.0, 135.0, 135.0, 125.0, 125.0, meshkernel::constants::missing::doubleValue,
                                      85.0, 95.0, 95.0, 85.0, 85.0, meshkernel::constants::missing::doubleValue,
                                      100.0, 105.0, 95.0, 100.0, meshkernel::constants::missing::doubleValue,
                                      120.0, 125.0, 115.0, 120.0};

    std::vector<double> invalidCellsY{0.0, 15.0, 15.0, 0.0, meshkernel::constants::missing::doubleValue,
                                      140.0, 140.0, 160.0, 160.0, 140.0, meshkernel::constants::missing::doubleValue,
                                      125.0, 125.0, 135.0, 135.0, 125.0, meshkernel::constants::missing::doubleValue,
                                      15.0, 15.0, 25.0, 25.0, 15.0, meshkernel::constants::missing::doubleValue,
                                      0.0, 15.0, 15.0, 0.0, meshkernel::constants::missing::doubleValue,
                                      0.0, 15.0, 15.0, 0.0};

    meshkernelapi::GeometryList invalidCells;
    invalidCells.num_coordinates = static_cast<int>(invalidCellsX.size());
    invalidCells.geometry_separator = meshkernel::constants::missing::doubleValue;
    invalidCells.coordinates_x = invalidCellsX.data();
    invalidCells.coordinates_y = invalidCellsY.data();

    errorCode = meshkernelapi::mkernel_mesh2d_delete_faces_in_polygons(meshKernelId, invalidCells);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::GeometryList geometryListIn;
    geometryListIn.geometry_separator = meshkernel::constants::missing::doubleValue;
    // geometryListIn.coordinates_x = xCoordinatesIn.data();
    // geometryListIn.coordinates_y = yCoordinatesIn.data();
    // geometryListIn.values = valuesIn.data();
    geometryListIn.num_coordinates = 0; // static_cast<int>(xCoordinatesIn.size());

    meshkernel::MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.max_num_refinement_iterations = 2;
    // meshRefinementParameters.min_edge_size = 200000.0;

    const int isTriangulationRequired = 0;
    const int projectToLandBoundaryOption = 1;
    meshkernelapi::GeometryList selectingPolygon{};
    meshkernelapi::GeometryList landBoundaries{};
    errorCode = mkernel_mesh2d_flip_edges(meshKernelId,
                                          isTriangulationRequired,
                                          projectToLandBoundaryOption,
                                          selectingPolygon,
                                          landBoundaries);

    // Execute
    errorCode = mkernel_mesh2d_refine_based_on_polygon(meshKernelId, geometryListIn, meshRefinementParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    //--------------------------------

    // Get the dimensions
    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Allocate memory and get data
    std::vector<int> edge_faces(mesh2d.num_edges * 2);
    std::vector<int> edge_nodes(mesh2d.num_edges * 2);
    std::vector<int> face_nodes(mesh2d.num_face_nodes);
    std::vector<int> face_edges(mesh2d.num_face_nodes);
    std::vector<int> nodes_per_face(mesh2d.num_faces);
    std::vector<double> node_x(mesh2d.num_nodes);
    std::vector<double> node_y(mesh2d.num_nodes);
    std::vector<double> edge_x(mesh2d.num_edges);
    std::vector<double> edge_y(mesh2d.num_edges);
    std::vector<double> face_x(mesh2d.num_faces);
    std::vector<double> face_y(mesh2d.num_faces);

    mesh2d.edge_faces = edge_faces.data();
    mesh2d.edge_nodes = edge_nodes.data();
    mesh2d.face_nodes = face_nodes.data();
    mesh2d.face_edges = face_edges.data();
    mesh2d.nodes_per_face = nodes_per_face.data();
    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_x = edge_x.data();
    mesh2d.edge_y = edge_y.data();
    mesh2d.face_x = face_x.data();
    mesh2d.face_y = face_y.data();
    errorCode = mkernel_mesh2d_get_data(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // errorCode = meshkernelapi::mkernel_mesh2d_casulli_refinement(meshKernelId);
    // ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // // Get the mesh dimensions
    // meshkernelapi::Mesh2D mesh2d{};
    // errorCode = meshkernelapi::mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    // ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // ASSERT_EQ(mesh2d.num_nodes, static_cast<int>(expectedX.size()));
    // ASSERT_EQ(2 * mesh2d.num_edges, static_cast<int>(expectedEdgesNodes.size()));

    // std::vector<int> edge_faces(mesh2d.num_edges * 2);
    // std::vector<int> edge_nodes(mesh2d.num_edges * 2);
    // std::vector<int> face_nodes(mesh2d.num_face_nodes);
    // std::vector<int> face_edges(mesh2d.num_face_nodes);
    // std::vector<int> nodes_per_face(mesh2d.num_faces);
    // std::vector<double> node_x(mesh2d.num_nodes);
    // std::vector<double> node_y(mesh2d.num_nodes);
    // std::vector<double> edge_x(mesh2d.num_edges);
    // std::vector<double> edge_y(mesh2d.num_edges);
    // std::vector<double> face_x(mesh2d.num_faces);
    // std::vector<double> face_y(mesh2d.num_faces);

    // mesh2d.edge_faces = edge_faces.data();
    // mesh2d.edge_nodes = edge_nodes.data();
    // mesh2d.face_nodes = face_nodes.data();
    // mesh2d.face_edges = face_edges.data();
    // mesh2d.nodes_per_face = nodes_per_face.data();
    // mesh2d.node_x = node_x.data();
    // mesh2d.node_y = node_y.data();
    // mesh2d.edge_x = edge_x.data();
    // mesh2d.edge_y = edge_y.data();
    // mesh2d.face_x = face_x.data();
    // mesh2d.face_y = face_y.data();
    // errorCode = mkernel_mesh2d_get_data(meshKernelId, mesh2d);
    // ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // const double tolerance = 1e-12;

    // for (int i = 0; i < mesh2d.num_nodes; ++i)
    // {
    //     EXPECT_NEAR(expectedX[i], node_x[i], tolerance);
    //     EXPECT_NEAR(expectedY[i], node_y[i], tolerance);
    // }

    // for (int i = 0; i < 2 * mesh2d.num_edges; ++i)
    // {
    //     EXPECT_EQ(expectedEdgesNodes[i], edge_nodes[i]);
    // }
}

TEST(InvalidCellsPolygonsTests, OMG2)
{
    // Prepare
    int meshKernelId;
    constexpr int isSpherical = 1;
    auto errorCode = meshkernelapi::mkernel_allocate_state(isSpherical, meshKernelId);

    meshkernel::MakeGridParameters gridParameters;
    gridParameters.block_size_x = 0.05;
    gridParameters.block_size_y = 0.05;
    gridParameters.origin_x = 147.75;
    gridParameters.origin_y = -40.4;
    gridParameters.upper_right_x = 147.9;
    gridParameters.upper_right_y = -40.25;

    errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid_on_extension(meshKernelId, gridParameters);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_curvilinear_convert_to_mesh2d(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    //--------------------------------

    std::vector<double> gebcoElev{-37, -37, -36, -36, -41, -39, -38, -34, -42, -32, -27, -6, -38, -6, -2, 2, -43, -42, -31, -20};
    std::vector<double> latNp{-40.397917, -40.364583, -40.33125, -40.297917, -40.264583};
    std::vector<double> lonNp{147.76875, 147.802083, 147.835417, 147.86875};

    int interpolationType;
    errorCode = meshkernelapi::mkernel_get_interpolation_type_double(interpolationType);

    meshkernelapi::GriddedSamples griddedSamples;

    griddedSamples.num_x = static_cast<int>(lonNp.size());
    griddedSamples.num_y = static_cast<int>(latNp.size());
    griddedSamples.x_coordinates = lonNp.data();
    griddedSamples.y_coordinates = latNp.data();
    griddedSamples.values = gebcoElev.data();
    griddedSamples.value_type = interpolationType;

    meshkernel::MeshRefinementParameters meshRefinementParameters;
    meshRefinementParameters.refinement_type = 1;
    // meshRefinementParameters.max_num_refinement_iterations = 3;
    meshRefinementParameters.min_edge_size = 500.0;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.smoothing_iterations = 2;
    meshRefinementParameters.max_courant_time = 120;

    meshkernelapi::GeometryList refinementPolygon;

    errorCode = meshkernelapi::mkernel_mesh2d_refine_based_on_gridded_samples(meshKernelId, refinementPolygon, griddedSamples, meshRefinementParameters, true);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    //--------------------------------

    std::vector<double> xx{147.87262097,
                           147.89,
                           147.90366935,
                           147.87262097,
                           147.87262097,
                           -999.0,
                           147.84157258,
                           147.84157258,
                           147.85931452,
                           147.85754032,
                           147.84157258};

    std::vector<double> yy{-40.33993264,
                           -40.33993264,
                           -40.36533877,
                           -40.35490411,
                           -40.33993264,
                           -999.0,
                           -40.30590657,
                           -40.29819399,
                           -40.29819399,
                           -40.30590657,
                           -40.30590657};

    // std::vector<double> xx{147.87262097,
    //                        147.89,
    //                        147.90366935,
    //                        147.87262097,
    //                        147.87262097,
    //                        -999.0,
    //                        147.84157258,
    //                        147.84157258,
    //                        147.85931452,
    //                        147.85754032,
    //                        147.84157258};

    // std::vector<double> yy{-40.33993264,
    //                        -40.33993264,
    //                        -40.36533877,
    //                        -40.35490411,
    //                        -40.33993264,
    //                        -999.0,
    //                        -40.30590657,
    //                        -40.29819399,
    //                        -40.29819399,
    //                        -40.30590657,
    //                        -40.30590657};

    meshkernelapi::GeometryList deletionPolygon;

    deletionPolygon.geometry_separator = meshkernel::constants::missing::doubleValue;
    deletionPolygon.num_coordinates = static_cast<int>(xx.size());
    deletionPolygon.coordinates_x = xx.data();
    deletionPolygon.coordinates_y = yy.data();

    errorCode = meshkernelapi::mkernel_mesh2d_delete_faces_in_polygons(meshKernelId, deletionPolygon);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    //--------------------------------

    // meshkernelapi::GeometryList refinementPolygon;
    // refinementPolygon.num_coordinates = static_cast<int>(patchX.size());
    // refinementPolygon.geometry_separator = meshkernel::constants::missing::doubleValue;
    // refinementPolygon.coordinates_x = patchX.data();
    // refinementPolygon.coordinates_y = patchY.data();

    // errorCode = meshkernelapi::mkernel_mesh2d_casulli_refinement_on_polygon(meshKernelId, refinementPolygon);
    // ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // //--------------------------------

    // std::vector<double> invalidCellsX{80.0, 85.0, 75.0, 80.0, meshkernel::constants::missing::doubleValue,
    //                                   180.0, 200.0, 200.0, 180.0, 180.0, meshkernel::constants::missing::doubleValue,
    //                                   125.0, 135.0, 135.0, 125.0, 125.0, meshkernel::constants::missing::doubleValue,
    //                                   85.0, 95.0, 95.0, 85.0, 85.0, meshkernel::constants::missing::doubleValue,
    //                                   100.0, 105.0, 95.0, 100.0, meshkernel::constants::missing::doubleValue,
    //                                   120.0, 125.0, 115.0, 120.0};

    // std::vector<double> invalidCellsY{0.0, 15.0, 15.0, 0.0, meshkernel::constants::missing::doubleValue,
    //                                   140.0, 140.0, 160.0, 160.0, 140.0, meshkernel::constants::missing::doubleValue,
    //                                   125.0, 125.0, 135.0, 135.0, 125.0, meshkernel::constants::missing::doubleValue,
    //                                   15.0, 15.0, 25.0, 25.0, 15.0, meshkernel::constants::missing::doubleValue,
    //                                   0.0, 15.0, 15.0, 0.0, meshkernel::constants::missing::doubleValue,
    //                                   0.0, 15.0, 15.0, 0.0};

    // meshkernelapi::GeometryList invalidCells;
    // invalidCells.num_coordinates = static_cast<int>(invalidCellsX.size());
    // invalidCells.geometry_separator = meshkernel::constants::missing::doubleValue;
    // invalidCells.coordinates_x = invalidCellsX.data();
    // invalidCells.coordinates_y = invalidCellsY.data();

    // errorCode = meshkernelapi::mkernel_mesh2d_delete_faces_in_polygons(meshKernelId, invalidCells);
    // ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // meshkernelapi::GeometryList geometryListIn;
    // geometryListIn.geometry_separator = meshkernel::constants::missing::doubleValue;
    // // geometryListIn.coordinates_x = xCoordinatesIn.data();
    // // geometryListIn.coordinates_y = yCoordinatesIn.data();
    // // geometryListIn.values = valuesIn.data();
    // geometryListIn.num_coordinates = 0; // static_cast<int>(xCoordinatesIn.size());

    // meshkernel::MeshRefinementParameters meshRefinementParameters2;
    // meshRefinementParameters2.max_num_refinement_iterations = 2;
    // meshRefinementParameters.min_edge_size = 200000.0;

    // const int isTriangulationRequired = 0;
    // const int projectToLandBoundaryOption = 1;
    // meshkernelapi::GeometryList selectingPolygon{};
    // meshkernelapi::GeometryList landBoundaries{};
    // errorCode = mkernel_mesh2d_flip_edges(meshKernelId,
    //                                       isTriangulationRequired,
    //                                       projectToLandBoundaryOption,
    //                                       selectingPolygon,
    //                                       landBoundaries);

    // // Execute
    // errorCode = mkernel_mesh2d_refine_based_on_polygon(meshKernelId, geometryListIn, meshRefinementParameters2);
    // ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    //--------------------------------

    // errorCode = meshkernelapi::mkernel_mesh2d_casulli_refinement(meshKernelId);
    // ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    //--------------------------------

    // Get the dimensions
    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Allocate memory and get data
    std::vector<int> edge_faces(mesh2d.num_edges * 2);
    std::vector<int> edge_nodes(mesh2d.num_edges * 2);
    std::vector<int> face_nodes(mesh2d.num_face_nodes);
    std::vector<int> face_edges(mesh2d.num_face_nodes);
    std::vector<int> nodes_per_face(mesh2d.num_faces);
    std::vector<double> node_x(mesh2d.num_nodes);
    std::vector<double> node_y(mesh2d.num_nodes);
    std::vector<double> edge_x(mesh2d.num_edges);
    std::vector<double> edge_y(mesh2d.num_edges);
    std::vector<double> face_x(mesh2d.num_faces);
    std::vector<double> face_y(mesh2d.num_faces);

    mesh2d.edge_faces = edge_faces.data();
    mesh2d.edge_nodes = edge_nodes.data();
    mesh2d.face_nodes = face_nodes.data();
    mesh2d.face_edges = face_edges.data();
    mesh2d.nodes_per_face = nodes_per_face.data();
    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_x = edge_x.data();
    mesh2d.edge_y = edge_y.data();
    mesh2d.face_x = face_x.data();
    mesh2d.face_y = face_y.data();
    errorCode = mkernel_mesh2d_get_data(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::Print(node_x, node_y, edge_nodes);
    meshkernel::SaveVtk(node_x, node_y, face_nodes, nodes_per_face, "mesh3.vtu");
}
