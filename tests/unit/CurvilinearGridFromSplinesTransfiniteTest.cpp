#include <MeshKernel/Mesh.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Splines.hpp>
#include <MeshKernel/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGridFromSplinesTransfinite.hpp>
#include <gtest/gtest.h>

TEST(CurvilinearGridFromSplinesTransfinite, FourSplines)
{
    std::vector<meshkernel::Point> firstSpline{{2.172341E+02, -2.415445E+01},
                                               {4.314185E+02, 1.947381E+02},
                                               {8.064374E+02, 3.987241E+02}};

    auto splines = std::make_shared<meshkernel::Splines>(meshkernel::Projections::cartesian);
    splines->AddSpline(firstSpline, 0, firstSpline.size());

    std::vector<meshkernel::Point> secondSpline{{2.894012E+01, 2.010146E+02},
                                                {2.344944E+02, 3.720490E+02},
                                                {6.424647E+02, 5.917262E+02}};

    splines->AddSpline(secondSpline, 0, secondSpline.size());

    std::vector<meshkernel::Point> thirdSpline{{2.265137E+00, 2.802553E+02},
                                               {2.799988E+02, -2.807726E+01}};
    splines->AddSpline(thirdSpline, 0, thirdSpline.size());

    std::vector<meshkernel::Point> fourthSpline{{5.067361E+02, 6.034946E+02},
                                                {7.475956E+02, 3.336055E+02}};
    splines->AddSpline(fourthSpline, 0, fourthSpline.size());

    meshkernelapi::CurvilinearParametersNative curvilinearParametersNative;
    curvilinearParametersNative.NRefinement = 40;
    curvilinearParametersNative.MRefinement = 20;
    meshkernel::CurvilinearGridFromSplinesTransfinite curvilinearGridFromSplinesTransfinite(splines, curvilinearParametersNative);

    meshkernel::CurvilinearGrid curvilinearGrid;
    curvilinearGridFromSplinesTransfinite.Compute(curvilinearGrid);

    // check the values
    constexpr double tolerance = 1e-6;
    ASSERT_NEAR(244.84733455150598, curvilinearGrid.m_grid[0][0].x, tolerance);
    ASSERT_NEAR(240.03223719861575, curvilinearGrid.m_grid[0][1].x, tolerance);
    ASSERT_NEAR(235.21721587684686, curvilinearGrid.m_grid[0][2].x, tolerance);
    ASSERT_NEAR(230.40187707543339, curvilinearGrid.m_grid[0][3].x, tolerance);
    ASSERT_NEAR(225.58666038317327, curvilinearGrid.m_grid[0][4].x, tolerance);
    ASSERT_NEAR(220.77175891290770, curvilinearGrid.m_grid[0][5].x, tolerance);
    ASSERT_NEAR(215.95654192442103, curvilinearGrid.m_grid[0][6].x, tolerance);
    ASSERT_NEAR(211.14151904110099, curvilinearGrid.m_grid[0][7].x, tolerance);
    ASSERT_NEAR(206.32630377949152, curvilinearGrid.m_grid[0][8].x, tolerance);
    ASSERT_NEAR(201.51108480926104, curvilinearGrid.m_grid[0][9].x, tolerance);
    ASSERT_NEAR(196.69606411139034, curvilinearGrid.m_grid[0][10].x, tolerance);

    ASSERT_NEAR(10.946966348412502, curvilinearGrid.m_grid[0][0].y, tolerance);
    ASSERT_NEAR(16.292559955716278, curvilinearGrid.m_grid[0][1].y, tolerance);
    ASSERT_NEAR(21.638069155281958, curvilinearGrid.m_grid[0][2].y, tolerance);
    ASSERT_NEAR(26.983930812344244, curvilinearGrid.m_grid[0][3].y, tolerance);
    ASSERT_NEAR(32.329656907057121, curvilinearGrid.m_grid[0][4].y, tolerance);
    ASSERT_NEAR(37.675033050656793, curvilinearGrid.m_grid[0][5].y, tolerance);
    ASSERT_NEAR(43.020759474232527, curvilinearGrid.m_grid[0][6].y, tolerance);
    ASSERT_NEAR(48.366270407390992, curvilinearGrid.m_grid[0][7].y, tolerance);
    ASSERT_NEAR(53.711994913833436, curvilinearGrid.m_grid[0][8].y, tolerance);
    ASSERT_NEAR(59.057723537488684, curvilinearGrid.m_grid[0][9].y, tolerance);
    ASSERT_NEAR(64.403232044419184, curvilinearGrid.m_grid[0][10].y, tolerance);

    ASSERT_NEAR(263.67028430842242, curvilinearGrid.m_grid[1][0].x, tolerance);
    ASSERT_NEAR(259.11363739326902, curvilinearGrid.m_grid[1][1].x, tolerance);
    ASSERT_NEAR(254.53691267796933, curvilinearGrid.m_grid[1][2].x, tolerance);
    ASSERT_NEAR(249.93698634609487, curvilinearGrid.m_grid[1][3].x, tolerance);
    ASSERT_NEAR(245.31456069699095, curvilinearGrid.m_grid[1][4].x, tolerance);
    ASSERT_NEAR(240.66785332275725, curvilinearGrid.m_grid[1][5].x, tolerance);
    ASSERT_NEAR(235.99933187522288, curvilinearGrid.m_grid[1][6].x, tolerance);
    ASSERT_NEAR(231.30940727936030, curvilinearGrid.m_grid[1][7].x, tolerance);
    ASSERT_NEAR(226.60252865287427, curvilinearGrid.m_grid[1][8].x, tolerance);
    ASSERT_NEAR(221.88022520931327, curvilinearGrid.m_grid[1][9].x, tolerance);
    ASSERT_NEAR(217.14743651601677, curvilinearGrid.m_grid[1][10].x, tolerance);

    ASSERT_NEAR(34.264668045745267, curvilinearGrid.m_grid[1][0].y, tolerance);
    ASSERT_NEAR(39.307546170495868, curvilinearGrid.m_grid[1][1].y, tolerance);
    ASSERT_NEAR(44.379080332661857, curvilinearGrid.m_grid[1][2].y, tolerance);
    ASSERT_NEAR(49.481517460105827, curvilinearGrid.m_grid[1][3].y, tolerance);
    ASSERT_NEAR(54.613111211730796, curvilinearGrid.m_grid[1][4].y, tolerance);
    ASSERT_NEAR(59.775023214376127, curvilinearGrid.m_grid[1][5].y, tolerance);
    ASSERT_NEAR(64.963841851929189, curvilinearGrid.m_grid[1][6].y, tolerance);
    ASSERT_NEAR(70.178519042215470, curvilinearGrid.m_grid[1][7].y, tolerance);
    ASSERT_NEAR(75.413628528250186, curvilinearGrid.m_grid[1][8].y, tolerance);
    ASSERT_NEAR(80.667056521594716, curvilinearGrid.m_grid[1][9].y, tolerance);
    ASSERT_NEAR(85.932983124208747, curvilinearGrid.m_grid[1][10].y, tolerance);
}

TEST(CurvilinearGridFromSplinesTransfinite, FourSplinesOneNSwapped)
{
    std::vector<meshkernel::Point> firstSpline{{2.172341E+02, -2.415445E+01},
                                               {4.314185E+02, 1.947381E+02},
                                               {8.064374E+02, 3.987241E+02}};

    auto splines = std::make_shared<meshkernel::Splines>(meshkernel::Projections::cartesian);
    splines->AddSpline(firstSpline, 0, firstSpline.size());

    std::vector<meshkernel::Point> secondSpline{{2.894012E+01, 2.010146E+02},
                                                {2.344944E+02, 3.720490E+02},
                                                {6.424647E+02, 5.917262E+02}};
    splines->AddSpline(secondSpline, 0, secondSpline.size());

    std::vector<meshkernel::Point> fourthSpline{{5.067361E+02, 6.034946E+02},
                                                {7.475956E+02, 3.336055E+02}};
    splines->AddSpline(fourthSpline, 0, fourthSpline.size());

    std::vector<meshkernel::Point> thirdSpline{{2.265137E+00, 2.802553E+02},
                                               {2.799988E+02, -2.807726E+01}};
    splines->AddSpline(thirdSpline, 0, thirdSpline.size());

    meshkernelapi::CurvilinearParametersNative curvilinearParametersNative;
    curvilinearParametersNative.NRefinement = 40;
    curvilinearParametersNative.MRefinement = 20;
    meshkernel::CurvilinearGridFromSplinesTransfinite curvilinearGridFromSplinesTransfinite(splines, curvilinearParametersNative);

    meshkernel::CurvilinearGrid curvilinearGrid;
    curvilinearGridFromSplinesTransfinite.Compute(curvilinearGrid);

    // check the values
    constexpr double tolerance = 1e-6;
    ASSERT_NEAR(244.84733455150598, curvilinearGrid.m_grid[0][0].x, tolerance);
    ASSERT_NEAR(240.03223719861575, curvilinearGrid.m_grid[0][1].x, tolerance);
    ASSERT_NEAR(235.21721587684686, curvilinearGrid.m_grid[0][2].x, tolerance);
    ASSERT_NEAR(230.40187707543339, curvilinearGrid.m_grid[0][3].x, tolerance);
    ASSERT_NEAR(225.58666038317327, curvilinearGrid.m_grid[0][4].x, tolerance);
    ASSERT_NEAR(220.77175891290770, curvilinearGrid.m_grid[0][5].x, tolerance);
    ASSERT_NEAR(215.95654192442103, curvilinearGrid.m_grid[0][6].x, tolerance);
    ASSERT_NEAR(211.14151904110099, curvilinearGrid.m_grid[0][7].x, tolerance);
    ASSERT_NEAR(206.32630377949152, curvilinearGrid.m_grid[0][8].x, tolerance);
    ASSERT_NEAR(201.51108480926104, curvilinearGrid.m_grid[0][9].x, tolerance);
    ASSERT_NEAR(196.69606411139034, curvilinearGrid.m_grid[0][10].x, tolerance);

    ASSERT_NEAR(10.946966348412502, curvilinearGrid.m_grid[0][0].y, tolerance);
    ASSERT_NEAR(16.292559955716278, curvilinearGrid.m_grid[0][1].y, tolerance);
    ASSERT_NEAR(21.638069155281958, curvilinearGrid.m_grid[0][2].y, tolerance);
    ASSERT_NEAR(26.983930812344244, curvilinearGrid.m_grid[0][3].y, tolerance);
    ASSERT_NEAR(32.329656907057121, curvilinearGrid.m_grid[0][4].y, tolerance);
    ASSERT_NEAR(37.675033050656793, curvilinearGrid.m_grid[0][5].y, tolerance);
    ASSERT_NEAR(43.020759474232527, curvilinearGrid.m_grid[0][6].y, tolerance);
    ASSERT_NEAR(48.366270407390992, curvilinearGrid.m_grid[0][7].y, tolerance);
    ASSERT_NEAR(53.711994913833436, curvilinearGrid.m_grid[0][8].y, tolerance);
    ASSERT_NEAR(59.057723537488684, curvilinearGrid.m_grid[0][9].y, tolerance);
    ASSERT_NEAR(64.403232044419184, curvilinearGrid.m_grid[0][10].y, tolerance);

    ASSERT_NEAR(263.67028430842242, curvilinearGrid.m_grid[1][0].x, tolerance);
    ASSERT_NEAR(259.11363739326902, curvilinearGrid.m_grid[1][1].x, tolerance);
    ASSERT_NEAR(254.53691267796933, curvilinearGrid.m_grid[1][2].x, tolerance);
    ASSERT_NEAR(249.93698634609487, curvilinearGrid.m_grid[1][3].x, tolerance);
    ASSERT_NEAR(245.31456069699095, curvilinearGrid.m_grid[1][4].x, tolerance);
    ASSERT_NEAR(240.66785332275725, curvilinearGrid.m_grid[1][5].x, tolerance);
    ASSERT_NEAR(235.99933187522288, curvilinearGrid.m_grid[1][6].x, tolerance);
    ASSERT_NEAR(231.30940727936030, curvilinearGrid.m_grid[1][7].x, tolerance);
    ASSERT_NEAR(226.60252865287427, curvilinearGrid.m_grid[1][8].x, tolerance);
    ASSERT_NEAR(221.88022520931327, curvilinearGrid.m_grid[1][9].x, tolerance);
    ASSERT_NEAR(217.14743651601677, curvilinearGrid.m_grid[1][10].x, tolerance);

    ASSERT_NEAR(34.264668045745267, curvilinearGrid.m_grid[1][0].y, tolerance);
    ASSERT_NEAR(39.307546170495868, curvilinearGrid.m_grid[1][1].y, tolerance);
    ASSERT_NEAR(44.379080332661857, curvilinearGrid.m_grid[1][2].y, tolerance);
    ASSERT_NEAR(49.481517460105827, curvilinearGrid.m_grid[1][3].y, tolerance);
    ASSERT_NEAR(54.613111211730796, curvilinearGrid.m_grid[1][4].y, tolerance);
    ASSERT_NEAR(59.775023214376127, curvilinearGrid.m_grid[1][5].y, tolerance);
    ASSERT_NEAR(64.963841851929189, curvilinearGrid.m_grid[1][6].y, tolerance);
    ASSERT_NEAR(70.178519042215470, curvilinearGrid.m_grid[1][7].y, tolerance);
    ASSERT_NEAR(75.413628528250186, curvilinearGrid.m_grid[1][8].y, tolerance);
    ASSERT_NEAR(80.667056521594716, curvilinearGrid.m_grid[1][9].y, tolerance);
    ASSERT_NEAR(85.932983124208747, curvilinearGrid.m_grid[1][10].y, tolerance);
}

TEST(CurvilinearGridFromSplinesTransfinite, FiveSplines)
{
    std::vector<meshkernel::Point> firstSpline{{2.172341E+02, -2.415445E+01},
                                               {4.314185E+02, 1.947381E+02},
                                               {8.064374E+02, 3.987241E+02}};

    auto splines = std::make_shared<meshkernel::Splines>(meshkernel::Projections::cartesian);
    splines->AddSpline(firstSpline, 0, firstSpline.size());

    std::vector<meshkernel::Point> secondSpline{{2.894012E+01, 2.010146E+02},
                                                {2.344944E+02, 3.720490E+02},
                                                {6.424647E+02, 5.917262E+02}};
    splines->AddSpline(secondSpline, 0, secondSpline.size());

    std::vector<meshkernel::Point> thirdSpline{{2.265137E+00, 2.802553E+02},
                                               {2.799988E+02, -2.807726E+01}};
    splines->AddSpline(thirdSpline, 0, thirdSpline.size());

    std::vector<meshkernel::Point> fourthSpline{{5.067361E+02, 6.034946E+02},
                                                {7.475956E+02, 3.336055E+02}};
    splines->AddSpline(fourthSpline, 0, fourthSpline.size());

    std::vector<meshkernel::Point> fifthSpline{{2.673223E+02, 4.706788E+02},
                                               {5.513401E+02, 1.545069E+02}};
    splines->AddSpline(fifthSpline, 0, fifthSpline.size());

    meshkernelapi::CurvilinearParametersNative curvilinearParametersNative;
    curvilinearParametersNative.NRefinement = 40;
    curvilinearParametersNative.MRefinement = 20;
    meshkernel::CurvilinearGridFromSplinesTransfinite curvilinearGridFromSplinesTransfinite(splines, curvilinearParametersNative);

    meshkernel::CurvilinearGrid curvilinearGrid;
    curvilinearGridFromSplinesTransfinite.Compute(curvilinearGrid);

    constexpr double tolerance = 1e-6;
    ASSERT_NEAR(244.84733455150598, curvilinearGrid.m_grid[0][0].x, tolerance);
    ASSERT_NEAR(240.03223719861575, curvilinearGrid.m_grid[0][1].x, tolerance);
    ASSERT_NEAR(235.21721587684686, curvilinearGrid.m_grid[0][2].x, tolerance);
    ASSERT_NEAR(230.40187707543339, curvilinearGrid.m_grid[0][3].x, tolerance);
    ASSERT_NEAR(225.58666038317327, curvilinearGrid.m_grid[0][4].x, tolerance);
    ASSERT_NEAR(220.77175891290770, curvilinearGrid.m_grid[0][5].x, tolerance);
    ASSERT_NEAR(215.95654192442103, curvilinearGrid.m_grid[0][6].x, tolerance);
    ASSERT_NEAR(211.14151904110099, curvilinearGrid.m_grid[0][7].x, tolerance);
    ASSERT_NEAR(206.32630377949152, curvilinearGrid.m_grid[0][8].x, tolerance);
    ASSERT_NEAR(201.51108480926104, curvilinearGrid.m_grid[0][9].x, tolerance);
    ASSERT_NEAR(196.69606411139034, curvilinearGrid.m_grid[0][10].x, tolerance);

    ASSERT_NEAR(10.946966348412502, curvilinearGrid.m_grid[0][0].y, tolerance);
    ASSERT_NEAR(16.292559955716278, curvilinearGrid.m_grid[0][1].y, tolerance);
    ASSERT_NEAR(21.638069155281958, curvilinearGrid.m_grid[0][2].y, tolerance);
    ASSERT_NEAR(26.983930812344244, curvilinearGrid.m_grid[0][3].y, tolerance);
    ASSERT_NEAR(32.329656907057121, curvilinearGrid.m_grid[0][4].y, tolerance);
    ASSERT_NEAR(37.675033050656793, curvilinearGrid.m_grid[0][5].y, tolerance);
    ASSERT_NEAR(43.020759474232527, curvilinearGrid.m_grid[0][6].y, tolerance);
    ASSERT_NEAR(48.366270407390992, curvilinearGrid.m_grid[0][7].y, tolerance);
    ASSERT_NEAR(53.711994913833436, curvilinearGrid.m_grid[0][8].y, tolerance);
    ASSERT_NEAR(59.057723537488684, curvilinearGrid.m_grid[0][9].y, tolerance);
    ASSERT_NEAR(64.403232044419184, curvilinearGrid.m_grid[0][10].y, tolerance);

    ASSERT_NEAR(255.89614293923407, curvilinearGrid.m_grid[1][0].x, tolerance);
    ASSERT_NEAR(251.26839070344425, curvilinearGrid.m_grid[1][1].x, tolerance);
    ASSERT_NEAR(246.62717589518911, curvilinearGrid.m_grid[1][2].x, tolerance);
    ASSERT_NEAR(241.96945582856105, curvilinearGrid.m_grid[1][3].x, tolerance);
    ASSERT_NEAR(237.29374836322307, curvilinearGrid.m_grid[1][4].x, tolerance);
    ASSERT_NEAR(232.59945837385263, curvilinearGrid.m_grid[1][5].x, tolerance);
    ASSERT_NEAR(227.88656387177011, curvilinearGrid.m_grid[1][6].x, tolerance);
    ASSERT_NEAR(223.15709488341233, curvilinearGrid.m_grid[1][7].x, tolerance);
    ASSERT_NEAR(218.41314240105709, curvilinearGrid.m_grid[1][8].x, tolerance);
    ASSERT_NEAR(213.65762819876193, curvilinearGrid.m_grid[1][9].x, tolerance);
    ASSERT_NEAR(208.89353710816445, curvilinearGrid.m_grid[1][10].x, tolerance);

    ASSERT_NEAR(24.731736741118521, curvilinearGrid.m_grid[1][0].y, tolerance);
    ASSERT_NEAR(29.842940652626876, curvilinearGrid.m_grid[1][1].y, tolerance);
    ASSERT_NEAR(34.982267945763468, curvilinearGrid.m_grid[1][2].y, tolerance);
    ASSERT_NEAR(40.148526703963910, curvilinearGrid.m_grid[1][3].y, tolerance);
    ASSERT_NEAR(45.340177298582923, curvilinearGrid.m_grid[1][4].y, tolerance);
    ASSERT_NEAR(50.555639961868344, curvilinearGrid.m_grid[1][5].y, tolerance);
    ASSERT_NEAR(55.793467784299104, curvilinearGrid.m_grid[1][6].y, tolerance);
    ASSERT_NEAR(61.050433278839293, curvilinearGrid.m_grid[1][7].y, tolerance);
    ASSERT_NEAR(66.323655962424397, curvilinearGrid.m_grid[1][8].y, tolerance);
    ASSERT_NEAR(71.609593896396262, curvilinearGrid.m_grid[1][9].y, tolerance);
    ASSERT_NEAR(76.904826220873304, curvilinearGrid.m_grid[1][10].y, tolerance);
}
