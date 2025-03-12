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

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLineMirror.hpp>
#include <TestUtils/MakeCurvilinearGrids.hpp>

TEST(CurvilinearLineMirror, Compute_LineMirrorOnLeftBoundary_ShouldCorrectlySumContributionsFromSubsequentColumns)
{
    // Set-up
    const auto curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 2.0, 3, 2);
    EXPECT_EQ(2, curvilinearGrid->NumN());
    EXPECT_EQ(3, curvilinearGrid->NumM());

    constexpr double f = 1.2;
    meshkernel::CurvilinearGridLineMirror curvilinearLineMirror(*curvilinearGrid, f);
    curvilinearLineMirror.SetLine({0, 0}, {0, 2});

    const auto p0 = curvilinearGrid->GetNode(1, 0);
    const auto p1 = curvilinearGrid->GetNode(1, 1);

    // Execute
    [[maybe_unused]] auto dummyUndoAction = curvilinearLineMirror.Compute();

    EXPECT_EQ(4, curvilinearGrid->NumM());
    EXPECT_EQ(2, curvilinearGrid->NumN());

    // Asserts
    constexpr double tolerance = 1e-6;
    const auto p_expected = (1 + f) * p0 + (-f) * p1;
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
    meshkernel::CurvilinearGridLineMirror curvilinearLineMirror(*curvilinearGrid, f);
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
    meshkernel::CurvilinearGridLineMirror curvilinearLineMirror(*curvilinearGrid, 1.2);
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
    ASSERT_THROW(meshkernel::CurvilinearGridLineMirror(*curvilinearGrid, 0.0), std::invalid_argument);
}

TEST(CurvilinearLineMirror, Compute_LineMirrorOnUpperBoundary_ShouldAddFacesOnUpperBoundary)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();
    ASSERT_EQ(9, curvilinearGrid->NumM());
    meshkernel::CurvilinearGridLineMirror curvilinearLineMirror(*curvilinearGrid, 1.2);
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
    meshkernel::CurvilinearGridLineMirror curvilinearLineMirror(*curvilinearGrid, 1.2);
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
    meshkernel::CurvilinearGridLineMirror curvilinearLineMirror(*curvilinearGrid, 1.2);
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

class CurvilinearLineMirrorInternalBoundaryTest : public ::testing::TestWithParam<std::tuple<std::pair<int, int>, std::pair<int, int>, std::vector<std::pair<int, int>>, double>>
{
protected:
    void RunTest()
    {
        // Set-up
        const auto curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 1.0, 10, 10);

        // Delete internal nodes for creating a hole
        std::vector<meshkernel::Point> holeNodes = {
            {5, 3}, {5, 4}, {5, 5}, {6, 3}, {6, 4}, {6, 5}, {7, 3}, {7, 4}, {7, 5}};

        for (const auto& node : holeNodes)
        {
            [[maybe_unused]] auto dummyUndoAction = curvilinearGrid->DeleteNode(node);
        }

        // Get test parameters
        auto [start, end, mirroredCoords, expectedX] = GetParam();
        constexpr double f = 1.2;
        meshkernel::CurvilinearGridLineMirror curvilinearLineMirror(*curvilinearGrid, f);

        // Apply line mirroring
        meshkernel::Point startLinePoint{static_cast<double>(start.first), static_cast<double>(start.second)};
        meshkernel::Point secondLinePoint{static_cast<double>(end.first), static_cast<double>(end.second)};

        curvilinearLineMirror.SetLine(startLinePoint, secondLinePoint);
        [[maybe_unused]] auto dummyUndoAction = curvilinearLineMirror.Compute();

        const double tolerance = 1e-6;
        for (const auto& coord : mirroredCoords)
        {
            auto node = curvilinearGrid->GetNode(coord.first, coord.second);
            EXPECT_NEAR(expectedX, node.x, tolerance);
        }
    }
};

INSTANTIATE_TEST_SUITE_P(
    CurvilinearLineMirrorTests,
    CurvilinearLineMirrorInternalBoundaryTest,
    ::testing::Values(
        std::make_tuple(std::make_pair(4, 3), std::make_pair(4, 5),
                        std::vector<std::pair<int, int>>{{3, 5}, {4, 5}, {5, 5}},
                        5.2), // Right-side test
        std::make_tuple(std::make_pair(8, 3), std::make_pair(8, 5),
                        std::vector<std::pair<int, int>>{{3, 7}, {4, 7}, {5, 7}},
                        6.8) // Left-side test
        ));

TEST_P(CurvilinearLineMirrorInternalBoundaryTest, Compute_LineMirrorInsideHole_ShouldCorrectlyAddLine)
{
    RunTest();
}
