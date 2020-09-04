#include "../Mesh.hpp"
#include "../Entities.hpp"
#include "../Polygons.hpp"
#include "../Constants.cpp"
#include "../CurvilinearGrid.hpp"
#include "../CurvilinearGridFromSplinesTransfinite.cpp"
#include <gtest/gtest.h>

TEST(CurvilinearGridFromSplinesTransfinite, FourSplines)
{
    std::vector<GridGeom::Point> firstSpline;
    firstSpline.push_back(GridGeom::Point{ 2.172341E+02, -2.415445E+01 });
    firstSpline.push_back(GridGeom::Point{ 4.314185E+02, 1.947381E+02 });
    firstSpline.push_back(GridGeom::Point{ 8.064374E+02, 3.987241E+02 });

    auto splines = std::make_shared<GridGeom::Splines>(GridGeom::Projections::cartesian);
    bool success = splines->AddSpline(firstSpline, 0, firstSpline.size());
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> secondSpline;
    secondSpline.push_back(GridGeom::Point{ 2.894012E+01, 2.010146E+02 });
    secondSpline.push_back(GridGeom::Point{ 2.344944E+02, 3.720490E+02 });
    secondSpline.push_back(GridGeom::Point{ 6.424647E+02, 5.917262E+02 });
    success = splines->AddSpline(secondSpline, 0, secondSpline.size());
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> thirdSpline;
    thirdSpline.push_back(GridGeom::Point{ 2.265137E+00, 2.802553E+02 });
    thirdSpline.push_back(GridGeom::Point{ 2.799988E+02, -2.807726E+01 });
    success = splines->AddSpline(thirdSpline, 0, thirdSpline.size());
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> fourthSpline;
    fourthSpline.push_back(GridGeom::Point{ 5.067361E+02, 6.034946E+02 });
    fourthSpline.push_back(GridGeom::Point{ 7.475956E+02, 3.336055E+02 });
    success = splines->AddSpline(fourthSpline, 0, fourthSpline.size());
    ASSERT_TRUE(success);

    GridGeom::CurvilinearGridFromSplinesTransfinite curvilinearGridFromSplinesTransfinite(splines);

    GridGeomApi::CurvilinearParametersNative curvilinearParametersNative;
    curvilinearParametersNative.NRefinement = 40;
    curvilinearParametersNative.MRefinement = 20;
    curvilinearGridFromSplinesTransfinite.Set(curvilinearParametersNative);

    GridGeom::CurvilinearGrid curvilinearGrid;
    success = curvilinearGridFromSplinesTransfinite.Compute(curvilinearGrid);
    ASSERT_TRUE(success);

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
    ASSERT_NEAR(259.87111655945489, curvilinearGrid.m_grid[1][1].x, tolerance);
    ASSERT_NEAR(255.11348812342322, curvilinearGrid.m_grid[1][2].x, tolerance);
    ASSERT_NEAR(250.35555048158136, curvilinearGrid.m_grid[1][3].x, tolerance);
    ASSERT_NEAR(245.59773444154465, curvilinearGrid.m_grid[1][4].x, tolerance);
    ASSERT_NEAR(240.84023176016686, curvilinearGrid.m_grid[1][5].x, tolerance);
    ASSERT_NEAR(236.08241927564706, curvilinearGrid.m_grid[1][6].x, tolerance);
    ASSERT_NEAR(231.32479892771738, curvilinearGrid.m_grid[1][7].x, tolerance);
    ASSERT_NEAR(226.56699197029900, curvilinearGrid.m_grid[1][8].x, tolerance);
    ASSERT_NEAR(221.80918330170789, curvilinearGrid.m_grid[1][9].x, tolerance);
    ASSERT_NEAR(217.05157087339163, curvilinearGrid.m_grid[1][10].x, tolerance);

    ASSERT_NEAR(34.264668045745267, curvilinearGrid.m_grid[1][0].y, tolerance);
    ASSERT_NEAR(38.464566052057620, curvilinearGrid.m_grid[1][1].y, tolerance);
    ASSERT_NEAR(43.760683877699215, curvilinearGrid.m_grid[1][2].y, tolerance);
    ASSERT_NEAR(49.057146729216853, curvilinearGrid.m_grid[1][3].y, tolerance);
    ASSERT_NEAR(54.353475063988888, curvilinearGrid.m_grid[1][4].y, tolerance);
    ASSERT_NEAR(59.649455413876531, curvilinearGrid.m_grid[1][5].y, tolerance);
    ASSERT_NEAR(64.945781472508884, curvilinearGrid.m_grid[1][6].y, tolerance);
    ASSERT_NEAR(70.241894499612627, curvilinearGrid.m_grid[1][7].y, tolerance);
    ASSERT_NEAR(75.538216086511994, curvilinearGrid.m_grid[1][8].y, tolerance);
    ASSERT_NEAR(80.834540418041755, curvilinearGrid.m_grid[1][9].y, tolerance);
    ASSERT_NEAR(86.130647149070327, curvilinearGrid.m_grid[1][10].y, tolerance);
}

TEST(CurvilinearGridFromSplinesTransfinite, FourSplinesOneNSwapped)
{
    std::vector<GridGeom::Point> firstSpline;
    firstSpline.push_back(GridGeom::Point{ 2.172341E+02, -2.415445E+01 });
    firstSpline.push_back(GridGeom::Point{ 4.314185E+02, 1.947381E+02 });
    firstSpline.push_back(GridGeom::Point{ 8.064374E+02, 3.987241E+02 });

    auto splines = std::make_shared<GridGeom::Splines>(GridGeom::Projections::cartesian);
    bool success = splines->AddSpline(firstSpline, 0, firstSpline.size());
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> secondSpline;
    secondSpline.push_back(GridGeom::Point{ 2.894012E+01, 2.010146E+02 });
    secondSpline.push_back(GridGeom::Point{ 2.344944E+02, 3.720490E+02 });
    secondSpline.push_back(GridGeom::Point{ 6.424647E+02, 5.917262E+02 });
    success = splines->AddSpline(secondSpline, 0, secondSpline.size());
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> fourthSpline;
    fourthSpline.push_back(GridGeom::Point{ 5.067361E+02, 6.034946E+02 });
    fourthSpline.push_back(GridGeom::Point{ 7.475956E+02, 3.336055E+02 });
    success = splines->AddSpline(fourthSpline, 0, fourthSpline.size());
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> thirdSpline;
    thirdSpline.push_back(GridGeom::Point{ 2.265137E+00, 2.802553E+02 });
    thirdSpline.push_back(GridGeom::Point{ 2.799988E+02, -2.807726E+01 });
    success = splines->AddSpline(thirdSpline, 0, thirdSpline.size());
    ASSERT_TRUE(success);

    GridGeom::CurvilinearGridFromSplinesTransfinite curvilinearGridFromSplinesTransfinite(splines);
    
    GridGeomApi::CurvilinearParametersNative curvilinearParametersNative;
    curvilinearParametersNative.NRefinement = 40;
    curvilinearParametersNative.MRefinement = 20;
    curvilinearGridFromSplinesTransfinite.Set(curvilinearParametersNative);
    
    GridGeom::CurvilinearGrid curvilinearGrid;
    success = curvilinearGridFromSplinesTransfinite.Compute(curvilinearGrid);
    ASSERT_TRUE(success);

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
    ASSERT_NEAR(259.87111655945489, curvilinearGrid.m_grid[1][1].x, tolerance);
    ASSERT_NEAR(255.11348812342322, curvilinearGrid.m_grid[1][2].x, tolerance);
    ASSERT_NEAR(250.35555048158136, curvilinearGrid.m_grid[1][3].x, tolerance);
    ASSERT_NEAR(245.59773444154465, curvilinearGrid.m_grid[1][4].x, tolerance);
    ASSERT_NEAR(240.84023176016686, curvilinearGrid.m_grid[1][5].x, tolerance);
    ASSERT_NEAR(236.08241927564706, curvilinearGrid.m_grid[1][6].x, tolerance);
    ASSERT_NEAR(231.32479892771738, curvilinearGrid.m_grid[1][7].x, tolerance);
    ASSERT_NEAR(226.56699197029900, curvilinearGrid.m_grid[1][8].x, tolerance);
    ASSERT_NEAR(221.80918330170789, curvilinearGrid.m_grid[1][9].x, tolerance);
    ASSERT_NEAR(217.05157087339163, curvilinearGrid.m_grid[1][10].x, tolerance);

    ASSERT_NEAR(34.264668045745267, curvilinearGrid.m_grid[1][0].y, tolerance);
    ASSERT_NEAR(38.464566052057620, curvilinearGrid.m_grid[1][1].y, tolerance);
    ASSERT_NEAR(43.760683877699215, curvilinearGrid.m_grid[1][2].y, tolerance);
    ASSERT_NEAR(49.057146729216853, curvilinearGrid.m_grid[1][3].y, tolerance);
    ASSERT_NEAR(54.353475063988888, curvilinearGrid.m_grid[1][4].y, tolerance);
    ASSERT_NEAR(59.649455413876531, curvilinearGrid.m_grid[1][5].y, tolerance);
    ASSERT_NEAR(64.945781472508884, curvilinearGrid.m_grid[1][6].y, tolerance);
    ASSERT_NEAR(70.241894499612627, curvilinearGrid.m_grid[1][7].y, tolerance);
    ASSERT_NEAR(75.538216086511994, curvilinearGrid.m_grid[1][8].y, tolerance);
    ASSERT_NEAR(80.834540418041755, curvilinearGrid.m_grid[1][9].y, tolerance);
    ASSERT_NEAR(86.130647149070327, curvilinearGrid.m_grid[1][10].y, tolerance);
}

TEST(CurvilinearGridFromSplinesTransfinite, FiveSplines)
{
    std::vector<GridGeom::Point> firstSpline;
    firstSpline.push_back(GridGeom::Point{ 2.172341E+02, -2.415445E+01 });
    firstSpline.push_back(GridGeom::Point{ 4.314185E+02, 1.947381E+02 });
    firstSpline.push_back(GridGeom::Point{ 8.064374E+02, 3.987241E+02 });

    auto splines = std::make_shared<GridGeom::Splines>(GridGeom::Projections::cartesian);
    bool success = splines->AddSpline(firstSpline, 0, firstSpline.size());
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> secondSpline;
    secondSpline.push_back(GridGeom::Point{ 2.894012E+01, 2.010146E+02 });
    secondSpline.push_back(GridGeom::Point{ 2.344944E+02, 3.720490E+02 });
    secondSpline.push_back(GridGeom::Point{ 6.424647E+02, 5.917262E+02 });
    success = splines->AddSpline(secondSpline, 0, secondSpline.size());
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> thirdSpline;
    thirdSpline.push_back(GridGeom::Point{ 2.265137E+00, 2.802553E+02 });
    thirdSpline.push_back(GridGeom::Point{ 2.799988E+02, -2.807726E+01 });
    success = splines->AddSpline(thirdSpline, 0, thirdSpline.size());
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> fourthSpline;
    fourthSpline.push_back(GridGeom::Point{ 5.067361E+02, 6.034946E+02 });
    fourthSpline.push_back(GridGeom::Point{ 7.475956E+02, 3.336055E+02 });
    success = splines->AddSpline(fourthSpline, 0, fourthSpline.size());
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> fifthSpline;
    fifthSpline.push_back(GridGeom::Point{ 2.673223E+02, 4.706788E+02 });
    fifthSpline.push_back(GridGeom::Point{ 5.513401E+02, 1.545069E+02 });
    success = splines->AddSpline(fifthSpline, 0, fifthSpline.size());
    ASSERT_TRUE(success);

    GridGeom::CurvilinearGridFromSplinesTransfinite curvilinearGridFromSplinesTransfinite(splines);
    
    GridGeomApi::CurvilinearParametersNative curvilinearParametersNative;
    curvilinearParametersNative.NRefinement = 40;
    curvilinearParametersNative.MRefinement = 20;
    curvilinearGridFromSplinesTransfinite.Set(curvilinearParametersNative);
    
    GridGeom::CurvilinearGrid curvilinearGrid;
    success = curvilinearGridFromSplinesTransfinite.Compute(curvilinearGrid);
    ASSERT_TRUE(success);

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
    ASSERT_NEAR(251.94499680206604, curvilinearGrid.m_grid[1][1].x, tolerance);
    ASSERT_NEAR(247.16412810980611, curvilinearGrid.m_grid[1][2].x, tolerance);
    ASSERT_NEAR(242.38294374839489, curvilinearGrid.m_grid[1][3].x, tolerance);
    ASSERT_NEAR(237.60187401525590, curvilinearGrid.m_grid[1][4].x, tolerance);
    ASSERT_NEAR(232.82111604316972, curvilinearGrid.m_grid[1][5].x, tolerance);
    ASSERT_NEAR(228.04005951332920, curvilinearGrid.m_grid[1][6].x, tolerance);
    ASSERT_NEAR(223.25918473247705, curvilinearGrid.m_grid[1][7].x, tolerance);
    ASSERT_NEAR(218.47812352940775, curvilinearGrid.m_grid[1][8].x, tolerance);
    ASSERT_NEAR(213.69705929422133, curvilinearGrid.m_grid[1][9].x, tolerance);
    ASSERT_NEAR(208.91618836124883, curvilinearGrid.m_grid[1][10].x, tolerance);

    ASSERT_NEAR(24.731736741118521, curvilinearGrid.m_grid[1][0].y, tolerance);
    ASSERT_NEAR(29.209237191214346, curvilinearGrid.m_grid[1][1].y, tolerance);
    ASSERT_NEAR(34.522745499288575, curvilinearGrid.m_grid[1][2].y, tolerance);
    ASSERT_NEAR(39.836605428332383, curvilinearGrid.m_grid[1][3].y, tolerance);
    ASSERT_NEAR(45.150338752180389, curvilinearGrid.m_grid[1][4].y, tolerance);
    ASSERT_NEAR(50.463726365361609, curvilinearGrid.m_grid[1][5].y, tolerance);
    ASSERT_NEAR(55.777446554148980, curvilinearGrid.m_grid[1][6].y, tolerance);
    ASSERT_NEAR(61.090965546380446, curvilinearGrid.m_grid[1][7].y, tolerance);
    ASSERT_NEAR(66.404692503592742, curvilinearGrid.m_grid[1][8].y, tolerance);
    ASSERT_NEAR(71.718423610914982, curvilinearGrid.m_grid[1][9].y, tolerance);
    ASSERT_NEAR(77.031940667468888, curvilinearGrid.m_grid[1][10].y, tolerance);
}