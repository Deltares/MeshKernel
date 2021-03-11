#include <gtest/gtest.h>

#include <MeshKernel/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGridOrthogonalization.hpp>
#include <MeshKernel/Entities.hpp>

std::shared_ptr<meshkernel::CurvilinearGrid> SmallRealWorldCurvilinearGrid()
{
    std::vector<std::vector<meshkernel::Point>> grid{

        {{7.998379637459554942E+04, 3.669368953805413912E+05},
         {8.006531634932766610E+04, 3.669913884352280875E+05},
         {8.014597083144872158E+04, 3.670473627646105597E+05},
         {8.022590004201814008E+04, 3.671046293496827129E+05},
         {8.030524375682926620E+04, 3.671630169196527568E+05},
         {8.038174798275028297E+04, 3.672241090446270537E+05},
         {8.045825220867129974E+04, 3.672852011696013506E+05},
         {8.053475643459231651E+04, 3.673462932945756475E+05},
         {8.061126066051333328E+04, 3.674073854195499443E+05}},

        {{8.005399692539963871E+04, 3.668271786935172859E+05},
         {8.014473294047857053E+04, 3.668725835977831739E+05},
         {8.023492419500040705E+04, 3.669191881611925783E+05},
         {8.032467176522142836E+04, 3.669667697959434590E+05},
         {8.041405739198261290E+04, 3.670151484942396637E+05},
         {8.050509647671248240E+04, 3.670564889889827464E+05},
         {8.059518333982788317E+04, 3.670991234714745078E+05},
         {8.068433399410264974E+04, 3.671430301817245199E+05},
         {8.077256729947395797E+04, 3.671881834906909498E+05}},

        {{8.010795455132638745E+04, 3.667094073585610604E+05},
         {8.019987610980382306E+04, 3.667511485812049941E+05},
         {8.029244910801970400E+04, 3.667925081235445105E+05},
         {8.038589000320057676E+04, 3.668319564031905611E+05},
         {8.047986028559294937E+04, 3.668702388785204384E+05},
         {8.058312060736966669E+04, 3.668916236319428426E+05},
         {8.068569967994136096E+04, 3.669153340029007522E+05},
         {8.078764917724052793E+04, 3.669411437571864226E+05},
         {8.088902409693629306E+04, 3.669688415973505471E+05}},

        {{8.013393074394566065E+04, 3.666299991322114947E+05},
         {8.022657945494366868E+04, 3.666560012267631828E+05},
         {8.032039638097764691E+04, 3.666791259040951263E+05},
         {8.041521952488367970E+04, 3.667002404040259426E+05},
         {8.051096004005162104E+04, 3.667200678389416425E+05},
         {8.062315126642497489E+04, 3.667184374811393791E+05},
         {8.073580093592485355E+04, 3.667174321629589540E+05},
         {8.084895945629126800E+04, 3.667167463103650953E+05},
         {8.096213238536757126E+04, 3.667181047524145106E+05}},

        {{8.015508428991910478E+04, 3.665319944788168068E+05},
         {8.024288669981209387E+04, 3.665448795976713882E+05},
         {8.033120056414291321E+04, 3.665524094749971409E+05},
         {8.041979925108155294E+04, 3.665552071627837140E+05},
         {8.050854555095927208E+04, 3.665539175055876258E+05},
         {8.062162400649990013E+04, 3.665473022831943817E+05},
         {8.073472454110847320E+04, 3.665398382516169222E+05},
         {8.084783925515387091E+04, 3.665315881046829163E+05},
         {8.096094935050149797E+04, 3.665225607507624663E+05}}};

    return std::make_shared<meshkernel::CurvilinearGrid>(std::move(grid), meshkernel::Projection::cartesian);
}

TEST(CurvilinearGridOrthogonalization, Compute_OnOrthogonalCurvilinearGrid_ShouldNotModifyGrid)
{
    // Set-up
    std::vector<std::vector<meshkernel::Point>> grid{
        {{0, 0}, {0, 10}, {0, 20}, {0, 30}},
        {{10, 0}, {10, 10}, {10, 20}, {10, 30}},
        {{20, 0}, {20, 10}, {20, 20}, {20, 30}},
        {{30, 0}, {30, 10}, {30, 20}, {30, 30}}};

    const auto curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(std::move(grid), meshkernel::Projection::cartesian);

    meshkernelapi::OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.OuterIterations = 1;
    orthogonalizationParameters.BoundaryIterations = 25;
    orthogonalizationParameters.InnerIterations = 25;
    orthogonalizationParameters.OrthogonalizationToSmoothingFactor = 0.975;
    meshkernel::CurvilinearGridOrthogonalization curvilinearGridOrthogonalization(curvilinearGrid, orthogonalizationParameters, {0, 0}, {30, 30});

    // Execute
    curvilinearGridOrthogonalization.Compute();

    // Assert nodes are on the same location because the grid is already orthogonal
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(0.0, curvilinearGrid->m_gridNodes[0][0].x, tolerance);
    ASSERT_NEAR(0.0, curvilinearGrid->m_gridNodes[0][1].x, tolerance);
    ASSERT_NEAR(0.0, curvilinearGrid->m_gridNodes[0][2].x, tolerance);
    ASSERT_NEAR(0.0, curvilinearGrid->m_gridNodes[0][3].x, tolerance);

    ASSERT_NEAR(10.0, curvilinearGrid->m_gridNodes[1][0].x, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid->m_gridNodes[1][1].x, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid->m_gridNodes[1][2].x, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid->m_gridNodes[1][3].x, tolerance);

    ASSERT_NEAR(20.0, curvilinearGrid->m_gridNodes[2][0].x, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid->m_gridNodes[2][1].x, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid->m_gridNodes[2][2].x, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid->m_gridNodes[2][3].x, tolerance);

    ASSERT_NEAR(30.0, curvilinearGrid->m_gridNodes[3][0].x, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid->m_gridNodes[3][1].x, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid->m_gridNodes[3][2].x, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid->m_gridNodes[3][3].x, tolerance);

    ASSERT_NEAR(0.0, curvilinearGrid->m_gridNodes[0][0].y, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid->m_gridNodes[0][1].y, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid->m_gridNodes[0][2].y, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid->m_gridNodes[0][3].y, tolerance);

    ASSERT_NEAR(0.0, curvilinearGrid->m_gridNodes[1][0].y, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid->m_gridNodes[1][1].y, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid->m_gridNodes[1][2].y, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid->m_gridNodes[1][3].y, tolerance);

    ASSERT_NEAR(0.0, curvilinearGrid->m_gridNodes[2][0].y, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid->m_gridNodes[2][1].y, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid->m_gridNodes[2][2].y, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid->m_gridNodes[2][3].y, tolerance);

    ASSERT_NEAR(0.0, curvilinearGrid->m_gridNodes[3][0].y, tolerance);
    ASSERT_NEAR(10.0, curvilinearGrid->m_gridNodes[3][1].y, tolerance);
    ASSERT_NEAR(20.0, curvilinearGrid->m_gridNodes[3][2].y, tolerance);
    ASSERT_NEAR(30.0, curvilinearGrid->m_gridNodes[3][3].y, tolerance);
}

TEST(CurvilinearGridOrthogonalization, Compute_OnONonOrthogonalCurvilinearGrid_ShouldOrthogonalizeGrid)
{
    // Set-up
    const auto curvilinearGrid = SmallRealWorldCurvilinearGrid();

    meshkernelapi::OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.OuterIterations = 2;
    orthogonalizationParameters.BoundaryIterations = 25;
    orthogonalizationParameters.InnerIterations = 25;
    orthogonalizationParameters.OrthogonalizationToSmoothingFactor = 0.975;
    meshkernel::CurvilinearGridOrthogonalization curvilinearGridOrthogonalization(curvilinearGrid, orthogonalizationParameters, {80154, 366530}, {80610, 367407});

    // Execute
    curvilinearGridOrthogonalization.Compute();

    // Assert
    constexpr double tolerance = 1e-6;

    ASSERT_NEAR(79983.796374595549, curvilinearGrid->m_gridNodes[0][0].x, tolerance);
    ASSERT_NEAR(80062.458148124875, curvilinearGrid->m_gridNodes[0][1].x, tolerance);
    ASSERT_NEAR(80139.231311212163, curvilinearGrid->m_gridNodes[0][2].x, tolerance);
    ASSERT_NEAR(80215.764666262330, curvilinearGrid->m_gridNodes[0][3].x, tolerance);
    ASSERT_NEAR(80293.275800678923, curvilinearGrid->m_gridNodes[0][4].x, tolerance);
    ASSERT_NEAR(80369.452240930041, curvilinearGrid->m_gridNodes[0][5].x, tolerance);
    ASSERT_NEAR(80446.051457860507, curvilinearGrid->m_gridNodes[0][6].x, tolerance);
    ASSERT_NEAR(80525.600108469254, curvilinearGrid->m_gridNodes[0][7].x, tolerance);
    ASSERT_NEAR(80611.260660513333, curvilinearGrid->m_gridNodes[0][8].x, tolerance);

    ASSERT_NEAR(80046.074998067590, curvilinearGrid->m_gridNodes[1][0].x, tolerance);
    ASSERT_NEAR(80126.368056777312, curvilinearGrid->m_gridNodes[1][1].x, tolerance);
    ASSERT_NEAR(80206.601210925364, curvilinearGrid->m_gridNodes[1][2].x, tolerance);
    ASSERT_NEAR(80286.287381013084, curvilinearGrid->m_gridNodes[1][3].x, tolerance);
    ASSERT_NEAR(80365.556974053165, curvilinearGrid->m_gridNodes[1][4].x, tolerance);
    ASSERT_NEAR(80443.587651289563, curvilinearGrid->m_gridNodes[1][5].x, tolerance);
    ASSERT_NEAR(80518.402282663534, curvilinearGrid->m_gridNodes[1][6].x, tolerance);
    ASSERT_NEAR(80591.265795857922, curvilinearGrid->m_gridNodes[1][7].x, tolerance);
    ASSERT_NEAR(80662.822394333387, curvilinearGrid->m_gridNodes[1][8].x, tolerance);

    ASSERT_NEAR(366936.89538054139, curvilinearGrid->m_gridNodes[0][0].y, tolerance);
    ASSERT_NEAR(366989.44915708830, curvilinearGrid->m_gridNodes[0][1].y, tolerance);
    ASSERT_NEAR(367042.59840712498, curvilinearGrid->m_gridNodes[0][2].y, tolerance);
    ASSERT_NEAR(367097.32852206373, curvilinearGrid->m_gridNodes[0][3].y, tolerance);
    ASSERT_NEAR(367153.90433931275, curvilinearGrid->m_gridNodes[0][4].y, tolerance);
    ASSERT_NEAR(367214.16882124962, curvilinearGrid->m_gridNodes[0][5].y, tolerance);
    ASSERT_NEAR(367275.49042299920, curvilinearGrid->m_gridNodes[0][6].y, tolerance);
    ASSERT_NEAR(367338.97499635734, curvilinearGrid->m_gridNodes[0][7].y, tolerance);
    ASSERT_NEAR(367407.38541954994, curvilinearGrid->m_gridNodes[0][8].y, tolerance);

    ASSERT_NEAR(366841.35146471433, curvilinearGrid->m_gridNodes[1][0].y, tolerance);
    ASSERT_NEAR(366887.54534671927, curvilinearGrid->m_gridNodes[1][1].y, tolerance);
    ASSERT_NEAR(366936.88064741943, curvilinearGrid->m_gridNodes[1][2].y, tolerance);
    ASSERT_NEAR(366990.14406133827, curvilinearGrid->m_gridNodes[1][3].y, tolerance);
    ASSERT_NEAR(367047.87322551024, curvilinearGrid->m_gridNodes[1][4].y, tolerance);
    ASSERT_NEAR(367110.98782724276, curvilinearGrid->m_gridNodes[1][5].y, tolerance);
    ASSERT_NEAR(367179.02693904703, curvilinearGrid->m_gridNodes[1][6].y, tolerance);
    ASSERT_NEAR(367254.16935723333, curvilinearGrid->m_gridNodes[1][7].y, tolerance);
    ASSERT_NEAR(367340.26390535629, curvilinearGrid->m_gridNodes[1][8].y, tolerance);
}