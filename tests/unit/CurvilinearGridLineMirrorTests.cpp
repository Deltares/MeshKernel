#include <gtest/gtest.h>

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLineMirror.hpp>
#include <TestUtils/MakeCurvilinearGrids.hpp>

TEST(CurvilinearLineMirror, Compute_LineMirrorOnBottomBoundary_ShouldAddFacesOnBottomBoundary)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();
    meshkernel::CurvilinearGridLineMirror curvilinearLineMirror(curvilinearGrid, 1.2);
    curvilinearLineMirror.SetLine({79983.0, 366936.2}, {80155.8, 366529.5});

    // Execute
    const auto mirroredGrid = curvilinearLineMirror.Compute();

    // Asserts
    const double tolerance = 1e-6;
    ASSERT_NEAR(79885.972404917018, mirroredGrid.m_gridNodes[0].front().x, tolerance);
    ASSERT_NEAR(79945.113707304932, mirroredGrid.m_gridNodes[1].front().x, tolerance);
    ASSERT_NEAR(79997.648681153471, mirroredGrid.m_gridNodes[2].front().x, tolerance);
    ASSERT_NEAR(80022.752290748060, mirroredGrid.m_gridNodes[3].front().x, tolerance);
    ASSERT_NEAR(80049.721398047535, mirroredGrid.m_gridNodes[4].front().x, tolerance);

    ASSERT_NEAR(366871.50371491740, mirroredGrid.m_gridNodes[0].front().y, tolerance);
    ASSERT_NEAR(366772.69280839822, mirroredGrid.m_gridNodes[1].front().y, tolerance);
    ASSERT_NEAR(366659.31789138837, mirroredGrid.m_gridNodes[2].front().y, tolerance);
    ASSERT_NEAR(366598.79661874950, mirroredGrid.m_gridNodes[3].front().y, tolerance);
    ASSERT_NEAR(366516.53233619139, mirroredGrid.m_gridNodes[4].front().y, tolerance);
}

TEST(CurvilinearLineMirror, Compute_LineMirrorOnUpperBoundary_ShouldAddFacesOnUpperBoundary)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();
    meshkernel::CurvilinearGridLineMirror curvilinearLineMirror(curvilinearGrid, 1.2);
    curvilinearLineMirror.SetLine({80960.2, 366520.7}, {80609.8, 367406.0});

    // Execute
    const auto mirroredGrid = curvilinearLineMirror.Compute();

    // Asserts
    const double tolerance = 1e-6;
    ASSERT_NEAR(80703.065731618568, mirroredGrid.m_gridNodes[0].back().x, tolerance);
    ASSERT_NEAR(80878.447265919545, mirroredGrid.m_gridNodes[1].back().x, tolerance);
    ASSERT_NEAR(81010.674000571220, mirroredGrid.m_gridNodes[2].back().x, tolerance);
    ASSERT_NEAR(81097.939900259138, mirroredGrid.m_gridNodes[3].back().x, tolerance);
    ASSERT_NEAR(81096.681464918671, mirroredGrid.m_gridNodes[4].back().x, tolerance);

    ASSERT_NEAR(367480.69596951915, mirroredGrid.m_gridNodes[0].back().y, tolerance);
    ASSERT_NEAR(367242.36746145069, mirroredGrid.m_gridNodes[1].back().y, tolerance);
    ASSERT_NEAR(367002.07900554762, mirroredGrid.m_gridNodes[2].back().y, tolerance);
    ASSERT_NEAR(366719.73488287395, mirroredGrid.m_gridNodes[3].back().y, tolerance);
    ASSERT_NEAR(366511.72792605805, mirroredGrid.m_gridNodes[4].back().y, tolerance);
}

TEST(CurvilinearLineMirror, Compute_LineMirrorOnLeftBoundary_ShouldAddFacesOnLeftBoundary)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();
    meshkernel::CurvilinearGridLineMirror curvilinearLineMirror(curvilinearGrid, 1.2);
    curvilinearLineMirror.SetLine({79983.0, 366936.2}, {80609.8, 367406.0});

    // Execute
    const auto mirroredGrid = curvilinearLineMirror.Compute();

    // Asserts
    const double tolerance = 1e-6;
    ASSERT_NEAR(79899.555713630645, mirroredGrid.m_gridNodes.front()[0].x, tolerance);
    ASSERT_NEAR(79970.016439946601, mirroredGrid.m_gridNodes.front()[1].x, tolerance);
    ASSERT_NEAR(80039.226795186711, mirroredGrid.m_gridNodes.front()[2].x, tolerance);
    ASSERT_NEAR(80107.373974174203, mirroredGrid.m_gridNodes.front()[3].x, tolerance);
    ASSERT_NEAR(80174.667394645265, mirroredGrid.m_gridNodes.front()[4].x, tolerance);
    ASSERT_NEAR(80233.729789995661, mirroredGrid.m_gridNodes.front()[5].x, tolerance);
    ASSERT_NEAR(80293.934851283411, mirroredGrid.m_gridNodes.front()[6].x, tolerance);
    ASSERT_NEAR(80355.263363179940, mirroredGrid.m_gridNodes.front()[7].x, tolerance);
    ASSERT_NEAR(80417.692693760589, mirroredGrid.m_gridNodes.front()[8].x, tolerance);

    ASSERT_NEAR(367068.55540497036, mirroredGrid.m_gridNodes.front()[0].y, tolerance);
    ASSERT_NEAR(367133.95424016198, mirroredGrid.m_gridNodes.front()[1].y, tolerance);
    ASSERT_NEAR(367201.17228871223, mirroredGrid.m_gridNodes.front()[2].y, tolerance);
    ASSERT_NEAR(367270.06081416988, mirroredGrid.m_gridNodes.front()[3].y, tolerance);
    ASSERT_NEAR(367340.45903014857, mirroredGrid.m_gridNodes.front()[4].y, tolerance);
    ASSERT_NEAR(367425.25311140029, mirroredGrid.m_gridNodes.front()[5].y, tolerance);
    ASSERT_NEAR(367508.49440735363, mirroredGrid.m_gridNodes.front()[6].y, tolerance);
    ASSERT_NEAR(367590.20902999706, mirroredGrid.m_gridNodes.front()[7].y, tolerance);
    ASSERT_NEAR(367670.42773418076, mirroredGrid.m_gridNodes.front()[8].y, tolerance);
}

TEST(CurvilinearLineMirror, Compute_LineMirrorOnRightBoundary_ShouldAddFacesOnRightBoundary)
{
    // Set-up
    const auto curvilinearGrid = MakeSmallCurvilinearGrid();
    meshkernel::CurvilinearGridLineMirror curvilinearLineMirror(curvilinearGrid, 1.2);
    curvilinearLineMirror.SetLine({80155.8, 366529.5}, {80960.2, 366520.72});

    // Execute
    const auto mirroredGrid = curvilinearLineMirror.Compute();

    // Asserts
    const double tolerance = 1e-6;
    ASSERT_NEAR(272501.90233055683, mirroredGrid.m_gridNodes.back()[0].x, tolerance);
    ASSERT_NEAR(272806.24608551903, mirroredGrid.m_gridNodes.back()[1].x, tolerance);
    ASSERT_NEAR(273113.11689828755, mirroredGrid.m_gridNodes.back()[2].x, tolerance);
    ASSERT_NEAR(273421.82178223983, mirroredGrid.m_gridNodes.back()[3].x, tolerance);
    ASSERT_NEAR(273731.95226017234, mirroredGrid.m_gridNodes.back()[4].x, tolerance);
    ASSERT_NEAR(274115.35433400975, mirroredGrid.m_gridNodes.back()[5].x, tolerance);
    ASSERT_NEAR(274499.35511354846, mirroredGrid.m_gridNodes.back()[6].x, tolerance);
    ASSERT_NEAR(274883.99770888803, mirroredGrid.m_gridNodes.back()[7].x, tolerance);
    ASSERT_NEAR(275268.64743354439, mirroredGrid.m_gridNodes.back()[8].x, tolerance);

    ASSERT_NEAR(1246326.3868120508, mirroredGrid.m_gridNodes.back()[0].y, tolerance);
    ASSERT_NEAR(1246385.9365869928, mirroredGrid.m_gridNodes.back()[1].y, tolerance);
    ASSERT_NEAR(1246430.2519299081, mirroredGrid.m_gridNodes.back()[2].y, tolerance);
    ASSERT_NEAR(1246461.7442429555, mirroredGrid.m_gridNodes.back()[3].y, tolerance);
    ASSERT_NEAR(1246482.6999190229, mirroredGrid.m_gridNodes.back()[4].y, tolerance);
    ASSERT_NEAR(1246466.1900003948, mirroredGrid.m_gridNodes.back()[5].y, tolerance);
    ASSERT_NEAR(1246448.5627491081, mirroredGrid.m_gridNodes.back()[6].y, tolerance);
    ASSERT_NEAR(1246429.5894027406, mirroredGrid.m_gridNodes.back()[7].y, tolerance);
    ASSERT_NEAR(1246411.3593545749, mirroredGrid.m_gridNodes.back()[8].y, tolerance);
}