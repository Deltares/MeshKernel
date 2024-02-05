#include <gtest/gtest.h>

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLineMirror.hpp>
#include <TestUtils/MakeCurvilinearGrids.hpp>

TEST(CurvilinearLineMirror, Compute_LineMirrorOnLeftBoundary_ShouldCorrectlySumContributionsFromSubsequentColumns)
{
    // Set-up
    const auto curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 2.0, 3, 2);
    EXPECT_EQ(3, curvilinearGrid->m_gridNodes.rows());
    EXPECT_EQ(2, curvilinearGrid->m_gridNodes.cols());

    constexpr double f = 1.2;
    meshkernel::CurvilinearGridLineMirror curvilinearLineMirror(*curvilinearGrid, f);
    curvilinearLineMirror.SetLine({0, 0}, {0, 2});

    const auto p0 = curvilinearGrid->m_gridNodes(0, 1);
    const auto p1 = curvilinearGrid->m_gridNodes(1, 1);

    // Execute
    curvilinearLineMirror.Compute();

    EXPECT_EQ(4, curvilinearGrid->m_gridNodes.rows());
    EXPECT_EQ(2, curvilinearGrid->m_gridNodes.cols());

    // Asserts
    constexpr double tolerance = 1e-6;
    const auto p_expected = (1 + f) * p0 + (-f) * p1;
    const auto p_actual = curvilinearGrid->m_gridNodes(0, 1);
    ASSERT_TRUE(meshkernel::IsEqual(p_expected, p_actual, tolerance));
}

TEST(CurvilinearLineMirror, Compute_LineMirrorOnRightBoundary_ShouldCorrectlySumContributionsFromPrecedingColumns)
{
    // Set-up
    const auto curvilinearGrid = MakeCurvilinearGrid(0.0, 0.0, 1.0, 2.0, 3, 2);
    EXPECT_EQ(3, curvilinearGrid->m_gridNodes.rows());
    EXPECT_EQ(2, curvilinearGrid->m_gridNodes.cols());

    constexpr double f = 1.2;
    meshkernel::CurvilinearGridLineMirror curvilinearLineMirror(*curvilinearGrid, f);
    curvilinearLineMirror.SetLine({2, 0}, {2, 2});

    const auto p0 = curvilinearGrid->m_gridNodes(2, 1);
    const auto p1 = curvilinearGrid->m_gridNodes(1, 1);

    // Execute
    curvilinearLineMirror.Compute();

    EXPECT_EQ(4, curvilinearGrid->m_gridNodes.rows());
    EXPECT_EQ(2, curvilinearGrid->m_gridNodes.cols());

    // Asserts
    constexpr double tolerance = 1e-6;
    const auto p_expected = (1 + f) * p0 + (-f) * p1;
    const auto p_actual = curvilinearGrid->m_gridNodes(3, 1);
    ASSERT_TRUE(meshkernel::IsEqual(p_expected, p_actual, tolerance));
}

TEST(CurvilinearLineMirror, Compute_LineMirrorOnBottomBoundary_ShouldAddFacesOnBottomBoundary)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();
    meshkernel::CurvilinearGridLineMirror curvilinearLineMirror(*curvilinearGrid, 1.2);
    curvilinearLineMirror.SetLine({79983.0, 366936.2}, {80155.8, 366529.5});

    // Execute
    curvilinearLineMirror.Compute();

    // Asserts
    constexpr double tolerance = 1e-6;
    ASSERT_NEAR(79885.972404917018, curvilinearGrid->m_gridNodes(0, 0).x, tolerance);
    ASSERT_NEAR(79945.113707304932, curvilinearGrid->m_gridNodes(1, 0).x, tolerance);
    ASSERT_NEAR(79997.648681153471, curvilinearGrid->m_gridNodes(2, 0).x, tolerance);
    ASSERT_NEAR(80022.752290748060, curvilinearGrid->m_gridNodes(3, 0).x, tolerance);
    ASSERT_NEAR(80049.721398047535, curvilinearGrid->m_gridNodes(4, 0).x, tolerance);

    ASSERT_NEAR(366871.50371491740, curvilinearGrid->m_gridNodes(0, 0).y, tolerance);
    ASSERT_NEAR(366772.69280839822, curvilinearGrid->m_gridNodes(1, 0).y, tolerance);
    ASSERT_NEAR(366659.31789138837, curvilinearGrid->m_gridNodes(2, 0).y, tolerance);
    ASSERT_NEAR(366598.79661874950, curvilinearGrid->m_gridNodes(3, 0).y, tolerance);
    ASSERT_NEAR(366516.53233619139, curvilinearGrid->m_gridNodes(4, 0).y, tolerance);
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
    meshkernel::CurvilinearGridLineMirror curvilinearLineMirror(*curvilinearGrid, 1.2);
    curvilinearLineMirror.SetLine({80960.2, 366520.7}, {80609.8, 367406.0});

    // Execute
    curvilinearLineMirror.Compute();

    // Asserts
    constexpr double tolerance = 1e-6;
    Eigen::Index const last = curvilinearGrid->m_gridNodes.cols() - 1;
    ASSERT_NEAR(80703.065731618568, curvilinearGrid->m_gridNodes(0, last).x, tolerance);
    ASSERT_NEAR(80878.447265919545, curvilinearGrid->m_gridNodes(1, last).x, tolerance);
    ASSERT_NEAR(81010.674000571220, curvilinearGrid->m_gridNodes(2, last).x, tolerance);
    ASSERT_NEAR(81097.939900259138, curvilinearGrid->m_gridNodes(3, last).x, tolerance);
    ASSERT_NEAR(81096.681464918671, curvilinearGrid->m_gridNodes(4, last).x, tolerance);

    ASSERT_NEAR(367480.69596951915, curvilinearGrid->m_gridNodes(0, last).y, tolerance);
    ASSERT_NEAR(367242.36746145069, curvilinearGrid->m_gridNodes(1, last).y, tolerance);
    ASSERT_NEAR(367002.07900554762, curvilinearGrid->m_gridNodes(2, last).y, tolerance);
    ASSERT_NEAR(366719.73488287395, curvilinearGrid->m_gridNodes(3, last).y, tolerance);
    ASSERT_NEAR(366511.72792605805, curvilinearGrid->m_gridNodes(4, last).y, tolerance);
}

TEST(CurvilinearLineMirror, Compute_LineMirrorOnLeftBoundary_ShouldAddFacesOnLeftBoundary)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();
    meshkernel::CurvilinearGridLineMirror curvilinearLineMirror(*curvilinearGrid, 1.2);
    curvilinearLineMirror.SetLine({79983.0, 366936.2}, {80609.8, 367406.0});

    // Execute
    curvilinearLineMirror.Compute();

    // Asserts
    constexpr double tolerance = 1e-6;
    ASSERT_NEAR(79899.555713630645, curvilinearGrid->m_gridNodes(0, 0).x, tolerance);
    ASSERT_NEAR(79970.016439946601, curvilinearGrid->m_gridNodes(0, 1).x, tolerance);
    ASSERT_NEAR(80039.226795186711, curvilinearGrid->m_gridNodes(0, 2).x, tolerance);
    ASSERT_NEAR(80107.373974174203, curvilinearGrid->m_gridNodes(0, 3).x, tolerance);
    ASSERT_NEAR(80174.667394645265, curvilinearGrid->m_gridNodes(0, 4).x, tolerance);
    ASSERT_NEAR(80233.729789995661, curvilinearGrid->m_gridNodes(0, 5).x, tolerance);
    ASSERT_NEAR(80293.934851283411, curvilinearGrid->m_gridNodes(0, 6).x, tolerance);
    ASSERT_NEAR(80355.263363179940, curvilinearGrid->m_gridNodes(0, 7).x, tolerance);
    ASSERT_NEAR(80417.692693760589, curvilinearGrid->m_gridNodes(0, 8).x, tolerance);

    ASSERT_NEAR(367068.55540497036, curvilinearGrid->m_gridNodes(0, 0).y, tolerance);
    ASSERT_NEAR(367133.95424016198, curvilinearGrid->m_gridNodes(0, 1).y, tolerance);
    ASSERT_NEAR(367201.17228871223, curvilinearGrid->m_gridNodes(0, 2).y, tolerance);
    ASSERT_NEAR(367270.06081416988, curvilinearGrid->m_gridNodes(0, 3).y, tolerance);
    ASSERT_NEAR(367340.45903014857, curvilinearGrid->m_gridNodes(0, 4).y, tolerance);
    ASSERT_NEAR(367425.25311140029, curvilinearGrid->m_gridNodes(0, 5).y, tolerance);
    ASSERT_NEAR(367508.49440735363, curvilinearGrid->m_gridNodes(0, 6).y, tolerance);
    ASSERT_NEAR(367590.20902999706, curvilinearGrid->m_gridNodes(0, 7).y, tolerance);
    ASSERT_NEAR(367670.42773418076, curvilinearGrid->m_gridNodes(0, 8).y, tolerance);
}

TEST(CurvilinearLineMirror, Compute_LineMirrorOnRightBoundary_ShouldAddFacesOnRightBoundary)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();
    meshkernel::CurvilinearGridLineMirror curvilinearLineMirror(*curvilinearGrid, 1.2);
    curvilinearLineMirror.SetLine({80155.8, 366529.5}, {80960.2, 366520.72});

    // Execute
    curvilinearLineMirror.Compute();

    // Asserts
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
    Eigen::Index const last = curvilinearGrid->m_gridNodes.rows() - 1;
    for (int i = 0; i < 9; ++i)
    {
        EXPECT_TRUE(meshkernel::IsEqual(expected[i], curvilinearGrid->m_gridNodes(last, i), tolerance));
    }
}
