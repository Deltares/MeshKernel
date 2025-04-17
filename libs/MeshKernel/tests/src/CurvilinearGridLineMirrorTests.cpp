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

#include "MeshKernel/ConnectMeshes.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridLineMirror.hpp"
#include "MeshKernel/Mesh2DToCurvilinear.hpp"
#include "MeshKernel/UndoActions/UndoActionStack.hpp"
#include "MeshKernel/Utilities/Utilities.hpp"

#include "TestUtils/MakeCurvilinearGrids.hpp"
#include "TestUtils/MakeMeshes.hpp"

TEST(CurvilinearLineMirror, Compute_LineMirrorOnLeftBoundary_ShouldCorrectlySumContributionsFromSubsequentColumns)
{
    // Set-up
    const auto curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 2.0, 3, 2);
    EXPECT_EQ(2, curvilinearGrid->NumN());
    EXPECT_EQ(3, curvilinearGrid->NumM());

    constexpr double mirroringFactor = 1.2;
    constexpr int numRowsToMirror = 1;
    meshkernel::CurvilinearGridLineMirror curvilinearLineMirror(*curvilinearGrid, mirroringFactor, numRowsToMirror);
    curvilinearLineMirror.SetLine({0, 0}, {0, 2});

    const auto p0 = curvilinearGrid->GetNode(1, 0);
    const auto p1 = curvilinearGrid->GetNode(1, 1);

    // Execute
    [[maybe_unused]] auto dummyUndoAction = curvilinearLineMirror.Compute();

    EXPECT_EQ(4, curvilinearGrid->NumM());
    EXPECT_EQ(2, curvilinearGrid->NumN());

    // Asserts
    constexpr double tolerance = 1e-6;
    const auto p_expected = (1 + mirroringFactor) * p0 + (-mirroringFactor) * p1;
    const auto p_actual = curvilinearGrid->GetNode(1, 0);
    ASSERT_TRUE(meshkernel::IsEqual(p_expected, p_actual, tolerance));
}

TEST(CurvilinearLineMirror, Compute_LineMirrorOnRightBoundary_ShouldCorrectlySumContributionsFromPrecedingColumns)
{
    // Set-up
    const auto curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 2.0, 3, 2);
    EXPECT_EQ(3, curvilinearGrid->NumM());
    EXPECT_EQ(2, curvilinearGrid->NumN());

    constexpr double f = 1.2;
    constexpr int numRowsToMirror = 1;
    meshkernel::CurvilinearGridLineMirror curvilinearLineMirror(*curvilinearGrid, f, numRowsToMirror);
    curvilinearLineMirror.SetLine({2, 0}, {2, 2});

    const auto p0 = curvilinearGrid->GetNode(1, 2);
    const auto p1 = curvilinearGrid->GetNode(1, 1);

    // Execute
    [[maybe_unused]] auto dummyUndoAction = curvilinearLineMirror.Compute();

    EXPECT_EQ(4, curvilinearGrid->NumM());
    EXPECT_EQ(2, curvilinearGrid->NumN());

    // Asserts
    constexpr double tolerance = 1e-6;
    const auto p_expected = (1 + f) * p0 + (-f) * p1;
    const auto p_actual = curvilinearGrid->GetNode(1, 3);
    ASSERT_TRUE(meshkernel::IsEqual(p_expected, p_actual, tolerance));
}

TEST(CurvilinearLineMirror, Compute_LineMirrorOnBottomBoundary_ShouldAddFacesOnBottomBoundary)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();
    const double mirroringFactor = 1.2;
    const int numRowsToMirror = 1;
    meshkernel::CurvilinearGridLineMirror curvilinearLineMirror(*curvilinearGrid, mirroringFactor, numRowsToMirror);
    curvilinearLineMirror.SetLine({79983.0, 366936.2}, {80155.8, 366529.5});

    // Execute
    [[maybe_unused]] auto dummyUndoAction = curvilinearLineMirror.Compute();

    // Asserts
    constexpr double tolerance = 1e-6;
    ASSERT_NEAR(79885.972404917018, curvilinearGrid->GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(79945.113707304932, curvilinearGrid->GetNode(1, 0).x, tolerance);
    ASSERT_NEAR(79997.648681153471, curvilinearGrid->GetNode(2, 0).x, tolerance);
    ASSERT_NEAR(80022.752290748060, curvilinearGrid->GetNode(3, 0).x, tolerance);
    ASSERT_NEAR(80049.721398047535, curvilinearGrid->GetNode(4, 0).x, tolerance);

    ASSERT_NEAR(366871.50371491740, curvilinearGrid->GetNode(0, 0).y, tolerance);
    ASSERT_NEAR(366772.69280839822, curvilinearGrid->GetNode(1, 0).y, tolerance);
    ASSERT_NEAR(366659.31789138837, curvilinearGrid->GetNode(2, 0).y, tolerance);
    ASSERT_NEAR(366598.79661874950, curvilinearGrid->GetNode(3, 0).y, tolerance);
    ASSERT_NEAR(366516.53233619139, curvilinearGrid->GetNode(4, 0).y, tolerance);
}

TEST(CurvilinearLineMirror, Compute_LineMirrorOnBottomBoundaryWithZeroMirrowingFactor_ShouldNotAddFacesOnBottomBoundary)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();

    // Assert
    ASSERT_THROW(meshkernel::CurvilinearGridLineMirror(*curvilinearGrid, 0.0, 1), std::invalid_argument);
}

TEST(CurvilinearLineMirror, Compute_LineMirrorOnUpperBoundary_ShouldAddFacesOnUpperBoundary)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();
    ASSERT_EQ(9, curvilinearGrid->NumM());
    const int numRowsToMirror = 1;
    const double mirroringFactor = 1.2;
    meshkernel::CurvilinearGridLineMirror curvilinearLineMirror(*curvilinearGrid, mirroringFactor, numRowsToMirror);
    curvilinearLineMirror.SetLine({80960.2, 366520.7}, {80609.8, 367406.0});

    // Execute
    [[maybe_unused]] auto dummyUndoAction = curvilinearLineMirror.Compute();

    // Asserts
    ASSERT_EQ(10, curvilinearGrid->NumM());

    constexpr double tolerance = 1e-6;
    const auto last = (meshkernel::UInt)curvilinearGrid->NumM() - 1;
    ASSERT_NEAR(80703.065731618568, curvilinearGrid->GetNode(0, last).x, tolerance);
    ASSERT_NEAR(80878.447265919545, curvilinearGrid->GetNode(1, last).x, tolerance);
    ASSERT_NEAR(81010.674000571220, curvilinearGrid->GetNode(2, last).x, tolerance);
    ASSERT_NEAR(81097.939900259138, curvilinearGrid->GetNode(3, last).x, tolerance);
    ASSERT_NEAR(81096.681464918671, curvilinearGrid->GetNode(4, last).x, tolerance);

    ASSERT_NEAR(367480.69596951915, curvilinearGrid->GetNode(0, last).y, tolerance);
    ASSERT_NEAR(367242.36746145069, curvilinearGrid->GetNode(1, last).y, tolerance);
    ASSERT_NEAR(367002.07900554762, curvilinearGrid->GetNode(2, last).y, tolerance);
    ASSERT_NEAR(366719.73488287395, curvilinearGrid->GetNode(3, last).y, tolerance);
    ASSERT_NEAR(366511.72792605805, curvilinearGrid->GetNode(4, last).y, tolerance);
}

TEST(CurvilinearLineMirror, Compute_LineMirrorOnLeftBoundary_ShouldAddFacesOnLeftBoundary)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();
    const int numRowsToMirror = 1;
    const double mirroringFactor = 1.2;
    meshkernel::CurvilinearGridLineMirror curvilinearLineMirror(*curvilinearGrid, mirroringFactor, numRowsToMirror);
    curvilinearLineMirror.SetLine({79983.0, 366936.2}, {80609.8, 367406.0});

    // Execute
    [[maybe_unused]] auto dummyUndoAction = curvilinearLineMirror.Compute();

    // Asserts
    constexpr double tolerance = 1e-6;
    ASSERT_NEAR(79899.555713630645, curvilinearGrid->GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(79970.016439946601, curvilinearGrid->GetNode(0, 1).x, tolerance);
    ASSERT_NEAR(80039.226795186711, curvilinearGrid->GetNode(0, 2).x, tolerance);
    ASSERT_NEAR(80107.373974174203, curvilinearGrid->GetNode(0, 3).x, tolerance);
    ASSERT_NEAR(80174.667394645265, curvilinearGrid->GetNode(0, 4).x, tolerance);
    ASSERT_NEAR(80233.729789995661, curvilinearGrid->GetNode(0, 5).x, tolerance);
    ASSERT_NEAR(80293.934851283411, curvilinearGrid->GetNode(0, 6).x, tolerance);
    ASSERT_NEAR(80355.263363179940, curvilinearGrid->GetNode(0, 7).x, tolerance);
    ASSERT_NEAR(80417.692693760589, curvilinearGrid->GetNode(0, 8).x, tolerance);

    ASSERT_NEAR(367068.55540497036, curvilinearGrid->GetNode(0, 0).y, tolerance);
    ASSERT_NEAR(367133.95424016198, curvilinearGrid->GetNode(0, 1).y, tolerance);
    ASSERT_NEAR(367201.17228871223, curvilinearGrid->GetNode(0, 2).y, tolerance);
    ASSERT_NEAR(367270.06081416988, curvilinearGrid->GetNode(0, 3).y, tolerance);
    ASSERT_NEAR(367340.45903014857, curvilinearGrid->GetNode(0, 4).y, tolerance);
    ASSERT_NEAR(367425.25311140029, curvilinearGrid->GetNode(0, 5).y, tolerance);
    ASSERT_NEAR(367508.49440735363, curvilinearGrid->GetNode(0, 6).y, tolerance);
    ASSERT_NEAR(367590.20902999706, curvilinearGrid->GetNode(0, 7).y, tolerance);
    ASSERT_NEAR(367670.42773418076, curvilinearGrid->GetNode(0, 8).y, tolerance);
}

TEST(CurvilinearLineMirror, Compute_LineMirrorOnRightBoundary_ShouldAddFacesOnRightBoundary)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();
    ASSERT_EQ(5, curvilinearGrid->NumN());
    const int numRowsToMirror = 1;
    const double mirroringFactor = 1.2;
    meshkernel::CurvilinearGridLineMirror curvilinearLineMirror(*curvilinearGrid, mirroringFactor, numRowsToMirror);
    curvilinearLineMirror.SetLine({80155.8, 366529.5}, {80960.2, 366520.72});

    // Execute
    [[maybe_unused]] auto dummyUndoAction = curvilinearLineMirror.Compute();

    // Asserts
    ASSERT_EQ(6, curvilinearGrid->NumN());

    constexpr double tolerance = 1e-6;
    const std::vector<meshkernel::Point> expected{
        {80180.468545087249, 366414.38889474329},
        {80262.455393654236, 366411.53364276129},
        {80344.165583941241, 366400.34976007964},
        {80425.29492251901, 366381.16727329313},
        {80505.648164048485, 366354.53710556292},
        {80619.791294589813, 366341.94004566048},
        {80733.43286732884, 366326.7255580065},
        {80846.495013789012, 366309.39825786441},
        {80959.52970866223, 366287.90794878011}};
    const auto last = curvilinearGrid->NumN() - 1;
    for (int i = 0; i < 9; ++i)
    {
        EXPECT_TRUE(meshkernel::IsEqual(expected[i], curvilinearGrid->GetNode(last, i), tolerance));
    }
}

class CurvilinearLineMirrorInternalBoundaryTest : public ::testing::TestWithParam<std::tuple<std::pair<int, int>,
                                                                                             std::pair<int, int>,
                                                                                             int,
                                                                                             std::vector<std::pair<int, int>>,
                                                                                             std::vector<meshkernel::Point>>>
{
protected:
    void RunTest()
    {
        // Set-up
        const auto curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, 15, 15);

        // Delete internal nodes for creating a hole
        std::vector<meshkernel::Point> holeNodes = {
            {7, 6},
            {7, 7},
            {7, 8},
            {7, 9},

            {8, 6},
            {8, 7},
            {8, 8},
            {8, 9},

            {9, 6},
            {9, 7},
            {9, 8},
            {9, 9},

            {10, 6},
            {10, 7},
            {10, 8},
            {10, 9}};

        for (const auto& node : holeNodes)
        {
            [[maybe_unused]] auto dummyUndoAction = curvilinearGrid->DeleteNode(node);
        }

        // Get test parameters
        auto [start, end, numRowsToMirror, mirroredCoords, expectedCoordinates] = GetParam();
        constexpr double mirroringFactor = 1.2;
        meshkernel::CurvilinearGridLineMirror curvilinearLineMirror(*curvilinearGrid, mirroringFactor, numRowsToMirror);

        // Apply line mirroring
        meshkernel::Point startLinePoint{static_cast<double>(start.first), static_cast<double>(start.second)};
        meshkernel::Point secondLinePoint{static_cast<double>(end.first), static_cast<double>(end.second)};

        curvilinearLineMirror.SetLine(startLinePoint, secondLinePoint);
        [[maybe_unused]] auto dummyUndoAction = curvilinearLineMirror.Compute();

        const double tolerance = 1e-6;
        for (auto i = 0u; i < mirroredCoords.size(); ++i)
        {
            const auto [nIndex, mIndex] = mirroredCoords[i];
            auto node = curvilinearGrid->GetNode(nIndex, mIndex);
            ASSERT_NEAR(node.x, expectedCoordinates[i].x, tolerance);
            ASSERT_NEAR(node.y, expectedCoordinates[i].y, tolerance);
        }
    }
};

INSTANTIATE_TEST_SUITE_P(
    CurvilinearLineMirrorTests,
    CurvilinearLineMirrorInternalBoundaryTest,
    ::testing::Values(
        std::make_tuple(std::make_pair(6, 6),
                        std::make_pair(6, 9),
                        2,
                        std::vector<std::pair<int, int>>{{6, 7},
                                                         {7, 7},
                                                         {8, 7},
                                                         {9, 7},

                                                         {6, 8},
                                                         {7, 8},
                                                         {8, 8},
                                                         {9, 8}},

                        std::vector<meshkernel::Point>{{7.2, 6},
                                                       {7.2, 7},
                                                       {7.2, 8},
                                                       {7.2, 9},

                                                       {8.64, 6},
                                                       {8.64, 7},
                                                       {8.64, 8},
                                                       {8.64, 9}}), // Right-side test
        std::make_tuple(std::make_pair(11, 6),
                        std::make_pair(11, 9),
                        2,
                        std::vector<std::pair<int, int>>{
                            {6, 10},
                            {7, 10},
                            {8, 10},
                            {9, 10},

                            {6, 9},
                            {7, 9},
                            {8, 9},
                            {9, 9}},
                        std::vector<meshkernel::Point>{
                            {9.8, 6},
                            {9.8, 7},
                            {9.8, 8},
                            {9.8, 9},

                            {8.36, 6},
                            {8.36, 7},
                            {8.36, 8},
                            {8.36, 9}}), // Left-side test

        std::make_tuple(std::make_pair(7, 5),
                        std::make_pair(10, 5),
                        2,
                        std::vector<std::pair<int, int>>{
                            {6, 7},
                            {6, 8},
                            {6, 9},
                            {6, 10},

                            {7, 7},
                            {7, 8},
                            {7, 9},
                            {7, 10}},
                        std::vector<meshkernel::Point>{
                            {7.0, 6.2},
                            {8.0, 6.2},
                            {9.0, 6.2},
                            {10.0, 6.2},

                            {7.0, 7.64},
                            {8.0, 7.64},
                            {9.0, 7.64},
                            {10.0, 7.64}}), // Bottom-side test

        std::make_tuple(std::make_pair(7, 10),
                        std::make_pair(10, 10),
                        2,
                        std::vector<std::pair<int, int>>{
                            {9, 7},
                            {9, 8},
                            {9, 9},
                            {9, 10},

                            {8, 7},
                            {8, 8},
                            {8, 9},
                            {8, 10}},
                        std::vector<meshkernel::Point>{
                            {7.0, 8.8},
                            {8.0, 8.8},
                            {9.0, 8.8},
                            {10.0, 8.8},

                            {7.0, 7.36},
                            {8.0, 7.36},
                            {9.0, 7.36},
                            {10.0, 7.36}}) // Top-side test
        ));

TEST_P(CurvilinearLineMirrorInternalBoundaryTest, Compute_LineMirrorInsideHole_ShouldCorrectlyAddLine)
{
    RunTest();
}

namespace CurvilinearLineMirrorTest
{
    size_t CountValidNodes(const meshkernel::CurvilinearGrid& grid)
    {

        size_t count = 0;

        for (meshkernel::UInt n = 0; n < grid.NumN(); ++n)
        {
            for (meshkernel::UInt m = 0; m < grid.NumM(); ++m)
            {
                if (grid.GetNode(n, m).IsValid())
                {
                    ++count;
                }
            }
        }

        return count;
    }
} // namespace CurvilinearLineMirrorTest

TEST(CurvilinearLineMirror, MultiStageMirrorCells)
{

    // This tests consists of several stages
    // 1. generate two unstructured grid, s1 slightly offset from the other
    // 2. combine these two grids into a single grid
    // 3. convert the merged grid into a single curvilinear grid.
    // 4. add 3 layers of mirror cells
    // 5. move the middle node of the mirrored cells
    // 6. check the nodes of the mirrored cells are correct
    // 7. undo node move and mirrored cells
    // 8. check the nodes are same as the original mesh

    // Prepare
    const auto mesh1 = MakeRectangularMeshForTesting(6,
                                                     6,
                                                     50.0,
                                                     50.0,
                                                     meshkernel::Projection::cartesian);
    // Prepare
    const auto mesh2 = MakeRectangularMeshForTesting(6,
                                                     6,
                                                     50.0,
                                                     50.0,
                                                     meshkernel::Projection::cartesian,
                                                     {50.0, 20.0});

    auto mergedMesh = meshkernel::Mesh2D::Merge(*mesh1, *mesh2);

    meshkernel::UndoActionStack undoStack;

    undoStack.Add(meshkernel::ConnectMeshes::Compute(*mergedMesh));

    meshkernel::Mesh2DToCurvilinear toClg(*mergedMesh);

    auto convertedClg = toClg.Compute({2.5, 2.5});

    ASSERT_EQ(convertedClg->NumN(), 11);
    ASSERT_EQ(convertedClg->NumM(), 8);

    std::vector<meshkernel::Point> originalMeshPoints(convertedClg->ComputeNodes());

    meshkernel::CurvilinearGridLineMirror lineMirror(*convertedClg, 1.0, 3);

    lineMirror.SetLine({10.0, 50.0}, {30.0, 50.0});

    undoStack.Add(lineMirror.Compute());

    undoStack.Add(convertedClg->MoveNode(meshkernel::Point(20.0, 80.0), meshkernel::Point(20.0, 85.0)));

    std::vector<meshkernel::UInt> expectedIndices{87, 78, 69, 88, 79, 70, 89, 80, 71};
    std::vector<meshkernel::Point> mirroredPoints{{10.0, 60.0}, {20.0, 60.0}, {30.0, 60.0}, {10.0, 70.0}, {20.0, 70.0}, {30.0, 70.0}, {10.0, 80.0}, {20.0, 85.0}, {30.0, 80.0}};

    ASSERT_EQ(convertedClg->NumN(), 11);
    ASSERT_EQ(convertedClg->NumM(), 9);

    for (size_t i = 0; i < mirroredPoints.size(); ++i)
    {
        EXPECT_EQ(expectedIndices[i], convertedClg->FindLocationIndex(mirroredPoints[i], meshkernel::Location::Nodes));
    }

    undoStack.Undo(); // move node
    undoStack.Undo(); // line mirror

    std::vector<meshkernel::Point> finalMeshPoints(convertedClg->ComputeNodes());

    ASSERT_EQ(finalMeshPoints.size(), originalMeshPoints.size());

    ASSERT_EQ(convertedClg->NumN(), 11);
    ASSERT_EQ(convertedClg->NumM(), 8);

    for (size_t i = 0; i < originalMeshPoints.size(); ++i)
    {
        EXPECT_EQ(originalMeshPoints[i].x, finalMeshPoints[i].x);
        EXPECT_EQ(originalMeshPoints[i].y, finalMeshPoints[i].y);
    }
}

TEST(CurvilinearLineMirror, UndoAndRecomputeMirrorCellsRightwards)
{

    // Prepare
    const auto mesh1 = MakeRectangularMeshForTesting(11,
                                                     11,
                                                     100.0,
                                                     100.0,
                                                     meshkernel::Projection::cartesian);
    // Prepare
    const auto mesh2 = MakeRectangularMeshForTesting(11,
                                                     11,
                                                     100.0,
                                                     100.0,
                                                     meshkernel::Projection::cartesian,
                                                     {70.0, 100.0});

    auto mergedMesh = meshkernel::Mesh2D::Merge(*mesh1, *mesh2);

    meshkernel::UndoActionStack undoStack;

    undoStack.Add(meshkernel::ConnectMeshes::Compute(*mergedMesh));

    meshkernel::Mesh2DToCurvilinear toClg(*mergedMesh);

    auto convertedClg = toClg.Compute({2.5, 2.5});

    ASSERT_EQ(convertedClg->NumN(), 18);
    ASSERT_EQ(convertedClg->NumM(), 21);

    meshkernel::CurvilinearGridLineMirror lineMirror(*convertedClg, 1.0, 9);
    lineMirror.SetLine({100.0, 30.0}, {100.0, 50.0});
    undoStack.Add(lineMirror.Compute());
    EXPECT_EQ(CurvilinearLineMirrorTest::CountValidNodes(*convertedClg), 265);

    ASSERT_EQ(convertedClg->NumN(), 20);
    ASSERT_EQ(convertedClg->NumM(), 21);

    undoStack.Undo(); // line mirror 9

    ASSERT_EQ(convertedClg->NumN(), 18);
    ASSERT_EQ(convertedClg->NumM(), 21);

    meshkernel::CurvilinearGridLineMirror lineMirror2(*convertedClg, 1.0, 10);
    lineMirror2.SetLine({100.0, 30.0}, {100.0, 50.0});
    undoStack.Add(lineMirror2.Compute());
    EXPECT_EQ(CurvilinearLineMirrorTest::CountValidNodes(*convertedClg), 268);

    ASSERT_EQ(convertedClg->NumN(), 21);
    ASSERT_EQ(convertedClg->NumM(), 21);

    undoStack.Undo(); // line mirror 10

    EXPECT_EQ(CurvilinearLineMirrorTest::CountValidNodes(*convertedClg), 238);
    ASSERT_EQ(convertedClg->NumN(), 18);
    ASSERT_EQ(convertedClg->NumM(), 21);
}

TEST(CurvilinearLineMirror, UndoAndRecomputeMirrorCellsLeftwards)
{

    // Prepare
    const auto mesh1 = MakeRectangularMeshForTesting(11,
                                                     11,
                                                     100.0,
                                                     100.0,
                                                     meshkernel::Projection::cartesian);
    // Prepare
    const auto mesh2 = MakeRectangularMeshForTesting(11,
                                                     11,
                                                     100.0,
                                                     100.0,
                                                     meshkernel::Projection::cartesian,
                                                     {70.0, 100.0});

    auto mergedMesh = meshkernel::Mesh2D::Merge(*mesh1, *mesh2);

    meshkernel::UndoActionStack undoStack;

    undoStack.Add(meshkernel::ConnectMeshes::Compute(*mergedMesh));

    meshkernel::Mesh2DToCurvilinear toClg(*mergedMesh);

    auto convertedClg = toClg.Compute({2.5, 2.5});

    ASSERT_EQ(convertedClg->NumN(), 18);
    ASSERT_EQ(convertedClg->NumM(), 21);
    EXPECT_EQ(CurvilinearLineMirrorTest::CountValidNodes(*convertedClg), 238);

    meshkernel::CurvilinearGridLineMirror lineMirror(*convertedClg, 1.0, 9);
    lineMirror.SetLine({70.0, 130.0}, {70.0, 150.0});
    undoStack.Add(lineMirror.Compute());

    ASSERT_EQ(convertedClg->NumN(), 20);
    ASSERT_EQ(convertedClg->NumM(), 21);
    EXPECT_EQ(CurvilinearLineMirrorTest::CountValidNodes(*convertedClg), 265);

    undoStack.Undo(); // line mirror 9

    ASSERT_EQ(convertedClg->NumN(), 18);
    ASSERT_EQ(convertedClg->NumM(), 21);
    EXPECT_EQ(CurvilinearLineMirrorTest::CountValidNodes(*convertedClg), 238);

    meshkernel::CurvilinearGridLineMirror lineMirror2(*convertedClg, 1.0, 10);
    lineMirror2.SetLine({70.0, 140.0}, {70.0, 160.0});
    undoStack.Add(lineMirror2.Compute());
    EXPECT_EQ(CurvilinearLineMirrorTest::CountValidNodes(*convertedClg), 268);

    ASSERT_EQ(convertedClg->NumN(), 21);
    ASSERT_EQ(convertedClg->NumM(), 21);

    undoStack.Undo(); // line mirror 10
    EXPECT_EQ(CurvilinearLineMirrorTest::CountValidNodes(*convertedClg), 238);

    ASSERT_EQ(convertedClg->NumN(), 18);
    ASSERT_EQ(convertedClg->NumM(), 21);
}

TEST(CurvilinearLineMirror, UndoAndRecomputeMirrorCellsDownwards)
{

    // Prepare
    const auto mesh1 = MakeRectangularMeshForTesting(11,
                                                     11,
                                                     100.0,
                                                     100.0,
                                                     meshkernel::Projection::cartesian);
    // Prepare
    const auto mesh2 = MakeRectangularMeshForTesting(11,
                                                     11,
                                                     100.0,
                                                     100.0,
                                                     meshkernel::Projection::cartesian,
                                                     {100.0, 70.0});

    auto mergedMesh = meshkernel::Mesh2D::Merge(*mesh1, *mesh2);

    meshkernel::UndoActionStack undoStack;

    undoStack.Add(meshkernel::ConnectMeshes::Compute(*mergedMesh));

    meshkernel::Mesh2DToCurvilinear toClg(*mergedMesh);

    auto convertedClg = toClg.Compute({2.5, 2.5});

    ASSERT_EQ(convertedClg->NumN(), 21);
    ASSERT_EQ(convertedClg->NumM(), 18);
    EXPECT_EQ(CurvilinearLineMirrorTest::CountValidNodes(*convertedClg), 238);

    meshkernel::CurvilinearGridLineMirror lineMirror(*convertedClg, 1.0, 5);
    lineMirror.SetLine({120.0, 70.0}, {140.0, 70.0});
    undoStack.Add(lineMirror.Compute());
    EXPECT_EQ(CurvilinearLineMirrorTest::CountValidNodes(*convertedClg), 253);

    undoStack.Undo(); // line mirror

    ASSERT_EQ(convertedClg->NumN(), 21);
    ASSERT_EQ(convertedClg->NumM(), 18);

    meshkernel::CurvilinearGridLineMirror lineMirror2(*convertedClg, 1.0, 8);
    lineMirror2.SetLine({120.0, 70.0}, {140.0, 70.0});
    undoStack.Add(lineMirror2.Compute());

    ASSERT_EQ(convertedClg->NumN(), 21);
    ASSERT_EQ(convertedClg->NumM(), 19);
    EXPECT_EQ(CurvilinearLineMirrorTest::CountValidNodes(*convertedClg), 262);

    undoStack.Undo(); // line mirror

    meshkernel::CurvilinearGridLineMirror lineMirror3(*convertedClg, 1.0, 5);
    lineMirror3.SetLine({120.0, 70.0}, {140.0, 70.0});
    undoStack.Add(lineMirror3.Compute());
    EXPECT_EQ(CurvilinearLineMirrorTest::CountValidNodes(*convertedClg), 253);
    undoStack.Undo(); // line mirror

    meshkernel::CurvilinearGridLineMirror lineMirror4(*convertedClg, 1.0, 9);
    lineMirror4.SetLine({120.0, 70.0}, {140.0, 70.0});
    undoStack.Add(lineMirror4.Compute());

    ASSERT_EQ(convertedClg->NumN(), 21);
    ASSERT_EQ(convertedClg->NumM(), 20);
    EXPECT_EQ(CurvilinearLineMirrorTest::CountValidNodes(*convertedClg), 265);
}

TEST(CurvilinearLineMirror, UndoAndExtendMirrorCells)
{

    // Prepare
    const auto mesh1 = MakeCurvilinearGrid(0.0, 0.0, 10.0, 10.0, 11, 11);
    meshkernel::UndoActionStack undoStack;

    meshkernel::CurvilinearGridLineMirror lineMirror(*mesh1, 1.0, 5);

    lineMirror.SetLine({20.0, 0.0}, {40.0, 0.0}); // bottom

    undoStack.Add(lineMirror.Compute());

    EXPECT_EQ(CurvilinearLineMirrorTest::CountValidNodes(*mesh1), 136);

    undoStack.Undo();

    meshkernel::CurvilinearGridLineMirror lineMirror2(*mesh1, 1.0, 3);

    lineMirror2.SetLine({20.0, 0.0}, {40.0, 0.0}); // bottom
    undoStack.Add(lineMirror2.Compute());
    EXPECT_EQ(CurvilinearLineMirrorTest::CountValidNodes(*mesh1), 130);

    undoStack.Undo();

    meshkernel::CurvilinearGridLineMirror lineMirror3(*mesh1, 1.0, 6);

    lineMirror3.SetLine({20.0, 0.0}, {40.0, 0.0}); // bottom
    undoStack.Add(lineMirror3.Compute());
    EXPECT_EQ(CurvilinearLineMirrorTest::CountValidNodes(*mesh1), 139);

    meshkernel::CurvilinearGridLineMirror lineMirror4(*mesh1, 1.0, 2);

    lineMirror4.SetLine({20.0, -60.0}, {40.0, -60.0}); // bottom
    undoStack.Add(lineMirror4.Compute());
    EXPECT_EQ(CurvilinearLineMirrorTest::CountValidNodes(*mesh1), 145);
}
