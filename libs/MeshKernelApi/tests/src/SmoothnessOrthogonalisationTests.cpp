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

#include <memory>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "MeshKernel/MeshTransformation.hpp"
#include "MeshKernel/Parameters.hpp"

#include "MeshKernelApi/BoundingBox.hpp"
#include "MeshKernelApi/GeometryList.hpp"
#include "MeshKernelApi/Mesh2D.hpp"
#include "MeshKernelApi/MeshKernel.hpp"

#include "CartesianApiTestFixture.hpp"
#include "TestUtils/Definitions.hpp"
#include "TestUtils/MakeCurvilinearGrids.hpp"
#include "TestUtils/MakeMeshes.hpp"
#include "TestUtils/SampleFileReader.hpp"

// namespace aliases
namespace mk = meshkernel;
namespace mkapi = meshkernelapi;

TEST(SmoothnessOrthogonalisationTests, DISABLED_Orthogonalize_OnInvalidMesh_ShouldThrowAMeshGeometryError)
{
    // Prepare
    int meshKernelId;
    int isGeographic = 0;
    meshkernelapi::mkernel_allocate_state(isGeographic, meshKernelId);

    auto [num_nodes, num_edges, node_x, node_y, node_type, edge_nodes, edge_type] =
        ReadLegacyMeshFile(TEST_FOLDER + "/data/InvalidMeshes/invalid_orthogonalization_net.nc");
    meshkernelapi::Mesh2D mesh2d;
    mesh2d.num_edges = static_cast<int>(num_edges);
    mesh2d.num_nodes = static_cast<int>(num_nodes);
    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_nodes = edge_nodes.data();

    auto errorCode = mkernel_mesh2d_set(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernel::OrthogonalizationParameters orthogonalizationParameters{};
    orthogonalizationParameters.outer_iterations = 1;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;

    meshkernelapi::GeometryList geometryList{};
    meshkernelapi::GeometryList landBoundaries{};

    // Execute
    errorCode = mkernel_mesh2d_initialize_orthogonalization(meshKernelId,
                                                            1,
                                                            orthogonalizationParameters,
                                                            landBoundaries,
                                                            geometryList);

    // Assert
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Assert there is a geometry error
    errorCode = meshkernelapi::mkernel_mesh2d_prepare_outer_iteration_orthogonalization(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::MeshGeometryErrorCode, errorCode);

    // Delete orthogonalization instance
    errorCode = meshkernelapi::mkernel_mesh2d_delete_orthogonalization(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Get the message
    auto exceptionMessage = std::make_unique<char[]>(512);
    errorCode = meshkernelapi::mkernel_get_error(exceptionMessage.get());
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    // Get the index of the invalid location
    int invalidIndex;
    int type;
    errorCode = meshkernelapi::mkernel_get_geometry_error(invalidIndex, type);
    ASSERT_EQ(static_cast<int>(meshkernel::Location::Nodes), type);
    ASSERT_EQ(478, invalidIndex);
}

TEST_F(CartesianApiTestFixture, OrthogonalizationThroughApi)
{
    // Set a new mesh in mesh
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Prepare
    meshkernel::OrthogonalizationParameters orthogonalizationParameters{};
    orthogonalizationParameters.outer_iterations = 1;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;

    meshkernelapi::GeometryList geometryList{};
    meshkernelapi::GeometryList landBoundaries{};

    // Execute
    auto errorCode = mkernel_mesh2d_initialize_orthogonalization(meshKernelId,
                                                                 1,
                                                                 orthogonalizationParameters,
                                                                 landBoundaries,
                                                                 geometryList);

    // Assert (nothing is done, just check that the api communication works)
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_prepare_outer_iteration_orthogonalization(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_compute_inner_ortogonalization_iteration(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_finalize_inner_ortogonalization_iteration(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    errorCode = meshkernelapi::mkernel_mesh2d_delete_orthogonalization(meshKernelId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::Mesh2D mesh2d{};
    errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);

    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(12, mesh2d.num_nodes);
    ASSERT_EQ(17, mesh2d.num_edges);
}

TEST_F(CartesianApiTestFixture, ComputeOrthogonalizationMesh2D_WithOrthogonalMesh2D_ShouldOrthogonalize)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    // Execute
    meshkernel::OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = 2;
    orthogonalizationParameters.boundary_iterations = 25;
    orthogonalizationParameters.inner_iterations = 25;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
    orthogonalizationParameters.orthogonalization_to_smoothing_factor_at_boundary = 1.0;
    orthogonalizationParameters.areal_to_angle_smoothing_factor = 1.0;
    // By using an empty polygon the entire mesh will be orthogonalized
    meshkernelapi::GeometryList polygons{};
    // No land boundaries accounted
    meshkernelapi::GeometryList landBoundaries{};

    auto errorCode = mkernel_mesh2d_compute_orthogonalization(meshKernelId, 1, orthogonalizationParameters, polygons, landBoundaries);

    // Assert
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}

TEST_F(CartesianApiTestFixture, GetOrthogonalityMesh2D_OnMesh2D_ShouldGetOrthogonality)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::Mesh2D mesh2d{};
    auto errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::GeometryList edgeOrthogonality;
    std::vector<double> values(mesh2d.num_edges);
    edgeOrthogonality.values = values.data();
    edgeOrthogonality.num_coordinates = mesh2d.num_edges;

    // Execute
    errorCode = mkernel_mesh2d_get_orthogonality(meshKernelId, edgeOrthogonality);

    // Assert
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}

TEST_F(CartesianApiTestFixture, GetSmoothnessMesh2D_OnMesh2D_ShouldGetSmoothness)
{
    // Prepare
    MakeMesh();
    auto const meshKernelId = GetMeshKernelId();

    meshkernelapi::Mesh2D mesh2d{};
    auto errorCode = mkernel_mesh2d_get_dimensions(meshKernelId, mesh2d);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::GeometryList edgeSmoothness;
    std::vector<double> values(mesh2d.num_edges);
    edgeSmoothness.values = values.data();
    edgeSmoothness.num_coordinates = mesh2d.num_edges;

    // Execute
    errorCode = mkernel_mesh2d_get_smoothness(meshKernelId, edgeSmoothness);

    // Assert
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
}

TEST_F(CartesianApiTestFixture, CurvilinearSetFrozenLinesOrthogonalize_ShouldSet2FrozenLines)
{
    // Setup
    const double deltaX = 10.0;
    const double deltaY = 10.0;

    int sizeX = 15;
    int sizeY = 15;

    MakeRectangularCurvilinearGrid(sizeX, sizeY, deltaX, deltaY);
    auto const meshKernelId = GetMeshKernelId();

    const std::vector<double> randomValues{-0.368462, -0.0413499, -0.281041, 0.178865, 0.434693, 0.0194164,
                                           -0.465428, 0.0297002, -0.492302, -0.433158, 0.186773, 0.430436,
                                           0.0269288, 0.153919, 0.201191, 0.262198, -0.452535, -0.171766,
                                           0.25641, -0.134661, 0.48255, 0.253356, -0.427314, 0.384707,
                                           -0.0635886, -0.0222682, -0.225093, -0.333493, 0.397656, -0.439436,
                                           0.00452289, -0.180967, -0.00602331, -0.409267, -0.426251, -0.115858,
                                           0.413817, -0.0355542, -0.449916, 0.270205, -0.374635, 0.188455, 0.129543};

    // displace grid
    auto random = [randomValues]() mutable
    {
        static size_t randomCounter = 0;

        if (randomCounter == randomValues.size() - 1)
        {
            randomCounter = 0;
        }
        else
        {
            ++randomCounter;
        }

        return randomValues[randomCounter];
    };

    // Use same indices for x and y for the block
    const int blockStartIndex = 2;
    const int blockEndIndex = 12;

    const int line1IndexX = 10;
    const int line1StartIndex = 0;
    const int line1EndIndex = 14;

    const int line2IndexX = 6;
    const int line2StartIndex = 4;
    const int line2EndIndex = 7;

    double line1StartPointX = -999.0;
    double line1StartPointY = -999.0;
    double line1EndPointX = -999.0;
    double line1EndPointY = -999.0;

    double line2StartPointX = -999.0;
    double line2StartPointY = -999.0;
    double line2EndPointX = -999.0;
    double line2EndPointY = -999.0;

    double blockStartPointX = -999.0;
    double blockStartPointY = -999.0;
    double blockEndPointX = -999.0;
    double blockEndPointY = -999.0;

    for (int i = 0; i <= sizeX; ++i)
    {
        double xPoint = static_cast<double>(i) * deltaX;

        for (int j = 0; j <= sizeY; ++j)
        {
            double yPoint = static_cast<double>(j) * deltaY;
            double xDisplaced = xPoint + 0.6 * random() * deltaX;
            double yDisplaced = yPoint + 0.6 * random() * deltaY;

            auto errorCode = meshkernelapi::mkernel_curvilinear_move_node(meshKernelId, xPoint, yPoint, xDisplaced, yDisplaced);
            ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

            if (i == j && i == blockStartIndex)
            {
                blockStartPointX = xDisplaced;
                blockStartPointY = yDisplaced;
            }

            if (i == j && i == blockEndIndex)
            {
                blockEndPointX = xDisplaced;
                blockEndPointY = yDisplaced;
            }

            // While jiggling the mesh, collect the start and end points of the frozen lines
            if (i == line1IndexX && j == line1StartIndex)
            {
                line1StartPointX = xDisplaced;
                line1StartPointY = yDisplaced;
            }

            if (i == line1IndexX && j == line1EndIndex)
            {
                line1EndPointX = xDisplaced;
                line1EndPointY = yDisplaced;
            }

            if (i == line2IndexX && j == line2StartIndex)
            {
                line2StartPointX = xDisplaced;
                line2StartPointY = yDisplaced;
            }

            if (i == line2IndexX && j == line2EndIndex)
            {
                line2EndPointX = xDisplaced;
                line2EndPointY = yDisplaced;
            }
        }
    }

    // Execute
    meshkernel::OrthogonalizationParameters const orthogonalizationParameters{};
    int frozenLinesId = -1;

    auto errorCode = meshkernelapi::mkernel_curvilinear_frozen_line_add(meshKernelId,
                                                                        line1StartPointX,
                                                                        line1StartPointY,
                                                                        line1EndPointX,
                                                                        line1EndPointY,
                                                                        frozenLinesId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(0, frozenLinesId);

    errorCode = meshkernelapi::mkernel_curvilinear_frozen_line_add(meshKernelId,
                                                                   line2StartPointX,
                                                                   line2StartPointY,
                                                                   line2EndPointX,
                                                                   line2EndPointY,
                                                                   frozenLinesId);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);
    ASSERT_EQ(1, frozenLinesId);

    errorCode = meshkernelapi::mkernel_curvilinear_orthogonalize(meshKernelId,
                                                                 orthogonalizationParameters,
                                                                 blockStartPointX,
                                                                 blockStartPointY,
                                                                 blockEndPointX,
                                                                 blockEndPointY);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    meshkernelapi::CurvilinearGrid clgData;
    errorCode = meshkernelapi::mkernel_curvilinear_get_dimensions(meshKernelId, clgData);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    std::vector<double> clgPointsX(clgData.num_m * clgData.num_n);
    std::vector<double> clgPointsY(clgData.num_m * clgData.num_n);
    clgData.node_x = clgPointsX.data();
    clgData.node_y = clgPointsY.data();

    errorCode = meshkernelapi::mkernel_curvilinear_get_data(meshKernelId, clgData);
    ASSERT_EQ(meshkernel::ExitCode::Success, errorCode);

    const double tolerance = 1.0e-10;
    // The start of the frozen line section for line 1
    size_t count = 10;

    std::vector<double> expectedLine1X{1.0289530000e+02, 9.7436116000e+01, 9.9618468400e+01, 9.8649442e+01,
                                       1.0238593600e+02, 1.0002713734e+02, 9.996386014e+01, 9.7442494000e+01,
                                       1.0248290200e+02, 9.7300504000e+01, 9.7752190000e+01, 1.0077725800e+02,
                                       9.9751900600e+01, 1.0107319000e+02, 1.0011649840e+02, 1.0017820120e+02};
    std::vector<double> expectedLine1Y{1.5201360000e+00, 1.2308242000e+01, 1.9866390800e+01, 2.7999042000e+01,
                                       3.7363384000e+01, 4.8914198000e+01, 5.7544398000e+01, 6.9304852000e+01,
                                       7.9786674800e+01, 9.1621230000e+01, 1.0113073000e+02, 1.0778922800e+02,
                                       1.1831375400e+02, 1.3260815800e+02, 1.3720743200e+02, 1.4704618800e+02};
    std::vector<double> expectedLine2X{5.7363384000e+01, 5.8914198000e+01, 5.7544398000e+01, 5.9304852000e+01};
    std::vector<double> expectedLine2Y{4.0027137340e+01, 4.9963860140e+01, 5.7442494000e+01, 7.2482902000e+01};

    for (size_t i = 0; i < expectedLine1X.size(); ++i)
    {
        EXPECT_NEAR(expectedLine1X[i], clgPointsX[count], tolerance);
        EXPECT_NEAR(expectedLine1Y[i], clgPointsY[count], tolerance);
        count += sizeY + 1;
    }

    // The start of the frozen line section for line 2.
    count = 6 + line2StartIndex * (sizeY + 1);

    for (size_t i = 0; i < expectedLine2X.size(); ++i)
    {
        EXPECT_NEAR(expectedLine2X[i], clgPointsX[count], tolerance);
        EXPECT_NEAR(expectedLine2Y[i], clgPointsY[count], tolerance);
        count += sizeY + 1;
    }
}
