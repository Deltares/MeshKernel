#include <gtest/gtest.h>

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLineMirror.hpp>
#include <TestUtils/MakeCurvilinearGrids.hpp>

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
    Eigen::Index const last = curvilinearGrid->m_gridNodes.rows() - 1;
    ASSERT_NEAR(272501.90233055683, curvilinearGrid->m_gridNodes(last, 0).x, tolerance);
    ASSERT_NEAR(272806.24608551903, curvilinearGrid->m_gridNodes(last, 1).x, tolerance);
    ASSERT_NEAR(273113.11689828755, curvilinearGrid->m_gridNodes(last, 2).x, tolerance);
    ASSERT_NEAR(273421.82178223983, curvilinearGrid->m_gridNodes(last, 3).x, tolerance);
    ASSERT_NEAR(273731.95226017234, curvilinearGrid->m_gridNodes(last, 4).x, tolerance);
    ASSERT_NEAR(274115.35433400975, curvilinearGrid->m_gridNodes(last, 5).x, tolerance);
    ASSERT_NEAR(274499.35511354846, curvilinearGrid->m_gridNodes(last, 6).x, tolerance);
    ASSERT_NEAR(274883.99770888803, curvilinearGrid->m_gridNodes(last, 7).x, tolerance);
    ASSERT_NEAR(275268.64743354439, curvilinearGrid->m_gridNodes(last, 8).x, tolerance);

    ASSERT_NEAR(1246326.3868120508, curvilinearGrid->m_gridNodes(last, 0).y, tolerance);
    ASSERT_NEAR(1246385.9365869928, curvilinearGrid->m_gridNodes(last, 1).y, tolerance);
    ASSERT_NEAR(1246430.2519299081, curvilinearGrid->m_gridNodes(last, 2).y, tolerance);
    ASSERT_NEAR(1246461.7442429555, curvilinearGrid->m_gridNodes(last, 3).y, tolerance);
    ASSERT_NEAR(1246482.6999190229, curvilinearGrid->m_gridNodes(last, 4).y, tolerance);
    ASSERT_NEAR(1246466.1900003948, curvilinearGrid->m_gridNodes(last, 5).y, tolerance);
    ASSERT_NEAR(1246448.5627491081, curvilinearGrid->m_gridNodes(last, 6).y, tolerance);
    ASSERT_NEAR(1246429.5894027406, curvilinearGrid->m_gridNodes(last, 7).y, tolerance);
    ASSERT_NEAR(1246411.3593545749, curvilinearGrid->m_gridNodes(last, 8).y, tolerance);
}
