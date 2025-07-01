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

#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/MeshEdgeCenters.hpp>
#include <MeshKernel/TriangulationInterpolation.hpp>
#include <TestUtils/Definitions.hpp>
#include <TestUtils/MakeMeshes.hpp>
#include <TestUtils/SampleFileReader.hpp>

TEST(TriangleInterpolation, TriangleInterpolation_OnNodesWithSphericalCoordinates_Shouldinterpolate)
{
    // Prepare
    auto mesh = MakeRectangularMeshForTesting(5, 5, 1.0, meshkernel::Projection::spherical);

    std::vector<meshkernel::Sample> samples{
        {1.5, 1.5, 2.0},
        {2.5, 1.5, 2.0},
        {3.5, 1.5, 2.0},
        {1.5, 2.5, 2.0},
        {2.5, 2.5, 1.0},
        {3.5, 2.5, 2.0},
        {1.5, 3.5, 2.0},
        {2.5, 3.5, 2.0},
        {3.5, 3.5, 2.0}};

    // Execute
    meshkernel::TriangulationInterpolation triangulationInterpolation(mesh->Nodes(), samples, meshkernel::Projection::spherical);
    triangulationInterpolation.Compute();

    // Assert
    auto interpolationResults = triangulationInterpolation.GetResults();
    constexpr double tolerance = 1e-6;
    ASSERT_NEAR(-999.00000000000000, interpolationResults[0], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[1], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[2], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[3], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[4], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[5], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[6], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[7], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[8], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[9], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[10], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[11], tolerance);
    ASSERT_NEAR(2.0000809724570234, interpolationResults[12], tolerance);
    ASSERT_NEAR(1.5000000000000002, interpolationResults[13], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[14], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[15], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[16], tolerance);
    ASSERT_NEAR(1.5000000000000000, interpolationResults[17], tolerance);
    ASSERT_NEAR(2.0000000000000000, interpolationResults[18], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[19], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[20], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[21], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[22], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[23], tolerance);
    ASSERT_NEAR(-999.00000000000000, interpolationResults[24], tolerance);
}
TEST(TriangleInterpolation, InterpolateOnNodes)
{
    // Set up
    std::vector<meshkernel::Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/TriangleInterpolationTests/inTestTriangleInterpolation.xyz");
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/TriangleInterpolationTests/simple_grid_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), static_cast<meshkernel::UInt>(0));

    meshkernel::TriangulationInterpolation triangulationInterpolation(mesh->Nodes(), samples, meshkernel::Projection::cartesian);
    triangulationInterpolation.Compute();

    const auto results = triangulationInterpolation.GetResults();

    // test internal results
    constexpr double tolerance = 1e-9;
    ASSERT_NEAR(2.0, results[20], tolerance);
    ASSERT_NEAR(2.0, results[21], tolerance);
    ASSERT_NEAR(2.0, results[22], tolerance);
    ASSERT_NEAR(2.0, results[23], tolerance);
    ASSERT_NEAR(2.0, results[24], tolerance);
    ASSERT_NEAR(1.0, results[25], tolerance);
    ASSERT_NEAR(1.0, results[26], tolerance);
    ASSERT_NEAR(1.0, results[27], tolerance);
    ASSERT_NEAR(1.0, results[28], tolerance);
    ASSERT_NEAR(1.0, results[29], tolerance);
}

TEST(TriangleInterpolation, InterpolateOnEdges)
{
    // Set up
    std::vector<meshkernel::Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/TriangleInterpolationTests/inTestTriangleInterpolation.xyz");
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/TriangleInterpolationTests/simple_grid_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), static_cast<meshkernel::UInt>(0));

    std::vector<meshkernel::Point> edgeCentres = meshkernel::algo::ComputeEdgeCentres(*mesh);
    meshkernel::TriangulationInterpolation triangulationInterpolation(edgeCentres, samples, meshkernel::Projection::cartesian);
    triangulationInterpolation.Compute();

    const auto results = triangulationInterpolation.GetResults();

    // test internal results
    constexpr double tolerance = 1e-9;
    ASSERT_NEAR(2.0, results[14], tolerance);
    ASSERT_NEAR(2.0, results[15], tolerance);
    ASSERT_NEAR(2.0, results[16], tolerance);
    ASSERT_NEAR(2.0, results[17], tolerance);
    ASSERT_NEAR(2.0, results[18], tolerance);
    ASSERT_NEAR(2.0, results[19], tolerance);
    ASSERT_NEAR(2.0, results[20], tolerance);
    ASSERT_NEAR(2.0, results[21], tolerance);
    ASSERT_NEAR(2.0, results[22], tolerance);
    ASSERT_NEAR(2.0, results[23], tolerance);
}
TEST(TriangleInterpolation, InterpolateOnFaces)
{
    // Set up
    std::vector<meshkernel::Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/TriangleInterpolationTests/inTestTriangleInterpolation.xyz");
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/TriangleInterpolationTests/simple_grid_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), static_cast<meshkernel::UInt>(0));

    meshkernel::TriangulationInterpolation triangulationInterpolation(mesh->m_facesMassCenters, samples, meshkernel::Projection::cartesian);
    triangulationInterpolation.Compute();

    const auto results = triangulationInterpolation.GetResults();

    // test internal results
    constexpr double tolerance = 1e-9;
    ASSERT_NEAR(2.0, results[15], tolerance);
    ASSERT_NEAR(2.0, results[16], tolerance);
    ASSERT_NEAR(2.0, results[17], tolerance);
    ASSERT_NEAR(2.0, results[18], tolerance);
    ASSERT_NEAR(2.0, results[19], tolerance);
    ASSERT_NEAR(0.0, results[20], tolerance);
    ASSERT_NEAR(0.0, results[21], tolerance);
    ASSERT_NEAR(0.0, results[22], tolerance);
    ASSERT_NEAR(0.0, results[23], tolerance);
    ASSERT_NEAR(0.0, results[24], tolerance);
    ASSERT_NEAR(0.0, results[25], tolerance);
}

TEST(TriangleInterpolation, InterpolateOnFacesUsingSphericalAccurateOption)
{
    // Set up
    std::vector<meshkernel::Sample> samples = ReadSampleFile(TEST_FOLDER + "/data/TriangleInterpolationTests/SphericalCutted.xyz");
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/TriangleInterpolationTests/SphericalCutted.nc", meshkernel::Projection::cartesian);
    ASSERT_GT(mesh->GetNumNodes(), static_cast<meshkernel::UInt>(0));
    ASSERT_GT(samples.size(), 0);

    meshkernel::TriangulationInterpolation triangulationInterpolation(mesh->m_facesMassCenters, samples, meshkernel::Projection::sphericalAccurate);
    triangulationInterpolation.Compute();

    const auto results = triangulationInterpolation.GetResults();

    // test internal results
    constexpr double tolerance = 1e-9;

    EXPECT_NEAR(-27.10868779318686, results[0], tolerance);
    EXPECT_NEAR(-26.00930207321026, results[1], tolerance);
    EXPECT_NEAR(-26.76102146173949, results[2], tolerance);
    EXPECT_NEAR(-26.49275221859572, results[3], tolerance);
    EXPECT_NEAR(-26.99345542316746, results[4], tolerance);
    EXPECT_NEAR(-26.89671692900715, results[5], tolerance);
    EXPECT_NEAR(-27.33232890626062, results[6], tolerance);
    EXPECT_NEAR(-28.37412297707877, results[7], tolerance);
    EXPECT_NEAR(-22.33073000076475, results[8], tolerance);
    EXPECT_NEAR(-30.42755591027593, results[9], tolerance);
    EXPECT_NEAR(-22.40438283953654, results[10], tolerance);
    EXPECT_NEAR(-13.36606610319540, results[11], tolerance);
    EXPECT_NEAR(-19.08344680981395, results[12], tolerance);
    EXPECT_NEAR(-33.57513770674175, results[13], tolerance);
    EXPECT_NEAR(-35.30789016071103, results[14], tolerance);
    EXPECT_NEAR(-32.53606308736682, results[15], tolerance);
    EXPECT_NEAR(-28.30846634340706, results[16], tolerance);
    EXPECT_NEAR(-26.42547910691146, results[17], tolerance);
    EXPECT_NEAR(-26.98947539456887, results[18], tolerance);
    EXPECT_NEAR(-29.54900509035342, results[19], tolerance);
}
