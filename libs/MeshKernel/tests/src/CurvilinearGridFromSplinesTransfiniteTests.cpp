#include <gtest/gtest.h>

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridFromSplinesTransfinite.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Splines.hpp>

using namespace meshkernel;

namespace testmk
{
    void TestCurvilinearGridFromSplines(const std::vector<meshkernel::Point>& firstSpline,
                                        const std::vector<meshkernel::Point>& secondSpline,
                                        const std::vector<meshkernel::Point>& thirdSpline,
                                        const std::vector<meshkernel::Point>& fourthSpline,
                                        const std::vector<meshkernel::Point>& expectedPoints)
    {
        auto splines = std::make_shared<Splines>(Projection::cartesian);

        splines->AddSpline(firstSpline);
        splines->AddSpline(secondSpline);
        splines->AddSpline(thirdSpline);
        splines->AddSpline(fourthSpline);

        CurvilinearParameters curvilinearParameters;
        curvilinearParameters.n_refinement = 4;
        curvilinearParameters.m_refinement = 4;
        CurvilinearGridFromSplinesTransfinite curvilinearGridFromSplinesTransfinite(splines, curvilinearParameters);

        const auto curvilinearGrid = curvilinearGridFromSplinesTransfinite.Compute();

        ASSERT_EQ(curvilinearGrid->m_nodes.size(), expectedPoints.size());

        constexpr double tolerance = 1e-4;

        for (size_t i = 0; i < curvilinearGrid->m_nodes.size(); ++i)
        {
            EXPECT_NEAR(expectedPoints[i].x, curvilinearGrid->m_nodes[i].x, tolerance);
            EXPECT_NEAR(expectedPoints[i].y, curvilinearGrid->m_nodes[i].y, tolerance);
        }
    }
} // namespace testmk

TEST(CurvilinearGridFromSplinesTransfinite, FourSplines)
{
    std::vector<Point> firstSpline{{2.172341E+02, -2.415445E+01},
                                   {4.314185E+02, 1.947381E+02},
                                   {8.064374E+02, 3.987241E+02}};

    auto splines = std::make_shared<Splines>(Projection::cartesian);
    splines->AddSpline(firstSpline);

    std::vector<Point> secondSpline{{2.894012E+01, 2.010146E+02},
                                    {2.344944E+02, 3.720490E+02},
                                    {6.424647E+02, 5.917262E+02}};
    splines->AddSpline(secondSpline);

    std::vector<Point> thirdSpline{{2.265137E+00, 2.802553E+02},
                                   {2.799988E+02, -2.807726E+01}};
    splines->AddSpline(thirdSpline);

    std::vector<Point> fourthSpline{{5.067361E+02, 6.034946E+02},
                                    {7.475956E+02, 3.336055E+02}};
    splines->AddSpline(fourthSpline);

    CurvilinearParameters curvilinearParameters;
    curvilinearParameters.n_refinement = 40;
    curvilinearParameters.m_refinement = 20;
    CurvilinearGridFromSplinesTransfinite curvilinearGridFromSplinesTransfinite(splines, curvilinearParameters);

    const auto curvilinearGrid = curvilinearGridFromSplinesTransfinite.Compute();

    // check the values
    constexpr double tolerance = 1e-6;
    ASSERT_NEAR(244.84733455150598, curvilinearGrid->GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(240.03223719861575, curvilinearGrid->GetNode(0, 1).x, tolerance);
    ASSERT_NEAR(235.21721587684686, curvilinearGrid->GetNode(0, 2).x, tolerance);
    ASSERT_NEAR(230.40187707543339, curvilinearGrid->GetNode(0, 3).x, tolerance);
    ASSERT_NEAR(225.58666038317327, curvilinearGrid->GetNode(0, 4).x, tolerance);
    ASSERT_NEAR(220.77175891290770, curvilinearGrid->GetNode(0, 5).x, tolerance);
    ASSERT_NEAR(215.95654192442103, curvilinearGrid->GetNode(0, 6).x, tolerance);
    ASSERT_NEAR(211.14151904110099, curvilinearGrid->GetNode(0, 7).x, tolerance);
    ASSERT_NEAR(206.32630377949152, curvilinearGrid->GetNode(0, 8).x, tolerance);
    ASSERT_NEAR(201.51108480926104, curvilinearGrid->GetNode(0, 9).x, tolerance);
    ASSERT_NEAR(196.69606411139034, curvilinearGrid->GetNode(0, 10).x, tolerance);

    ASSERT_NEAR(10.946966348412502, curvilinearGrid->GetNode(0, 0).y, tolerance);
    ASSERT_NEAR(16.292559955716278, curvilinearGrid->GetNode(0, 1).y, tolerance);
    ASSERT_NEAR(21.638069155281958, curvilinearGrid->GetNode(0, 2).y, tolerance);
    ASSERT_NEAR(26.983930812344244, curvilinearGrid->GetNode(0, 3).y, tolerance);
    ASSERT_NEAR(32.329656907057121, curvilinearGrid->GetNode(0, 4).y, tolerance);
    ASSERT_NEAR(37.675033050656793, curvilinearGrid->GetNode(0, 5).y, tolerance);
    ASSERT_NEAR(43.020759474232527, curvilinearGrid->GetNode(0, 6).y, tolerance);
    ASSERT_NEAR(48.366270407390992, curvilinearGrid->GetNode(0, 7).y, tolerance);
    ASSERT_NEAR(53.711994913833436, curvilinearGrid->GetNode(0, 8).y, tolerance);
    ASSERT_NEAR(59.057723537488684, curvilinearGrid->GetNode(0, 9).y, tolerance);
    ASSERT_NEAR(64.403232044419184, curvilinearGrid->GetNode(0, 10).y, tolerance);

    ASSERT_NEAR(263.67028430842242, curvilinearGrid->GetNode(1, 0).x, tolerance);
    ASSERT_NEAR(259.11363739326902, curvilinearGrid->GetNode(1, 1).x, tolerance);
    ASSERT_NEAR(254.53691267796933, curvilinearGrid->GetNode(1, 2).x, tolerance);
    ASSERT_NEAR(249.93698634609487, curvilinearGrid->GetNode(1, 3).x, tolerance);
    ASSERT_NEAR(245.31456069699095, curvilinearGrid->GetNode(1, 4).x, tolerance);
    ASSERT_NEAR(240.66785332275725, curvilinearGrid->GetNode(1, 5).x, tolerance);
    ASSERT_NEAR(235.99933187522288, curvilinearGrid->GetNode(1, 6).x, tolerance);
    ASSERT_NEAR(231.30940727936030, curvilinearGrid->GetNode(1, 7).x, tolerance);
    ASSERT_NEAR(226.60252865287427, curvilinearGrid->GetNode(1, 8).x, tolerance);
    ASSERT_NEAR(221.88022520931327, curvilinearGrid->GetNode(1, 9).x, tolerance);
    ASSERT_NEAR(217.14743651601677, curvilinearGrid->GetNode(1, 10).x, tolerance);

    ASSERT_NEAR(34.264668045745267, curvilinearGrid->GetNode(1, 0).y, tolerance);
    ASSERT_NEAR(39.307546170495868, curvilinearGrid->GetNode(1, 1).y, tolerance);
    ASSERT_NEAR(44.379080332661857, curvilinearGrid->GetNode(1, 2).y, tolerance);
    ASSERT_NEAR(49.481517460105827, curvilinearGrid->GetNode(1, 3).y, tolerance);
    ASSERT_NEAR(54.613111211730796, curvilinearGrid->GetNode(1, 4).y, tolerance);
    ASSERT_NEAR(59.775023214376127, curvilinearGrid->GetNode(1, 5).y, tolerance);
    ASSERT_NEAR(64.963841851929189, curvilinearGrid->GetNode(1, 6).y, tolerance);
    ASSERT_NEAR(70.178519042215470, curvilinearGrid->GetNode(1, 7).y, tolerance);
    ASSERT_NEAR(75.413628528250186, curvilinearGrid->GetNode(1, 8).y, tolerance);
    ASSERT_NEAR(80.667056521594716, curvilinearGrid->GetNode(1, 9).y, tolerance);
    ASSERT_NEAR(85.932983124208747, curvilinearGrid->GetNode(1, 10).y, tolerance);
}

TEST(CurvilinearGridFromSplinesTransfinite, FourSplinesOneNSwapped)
{
    std::vector<Point> firstSpline{{2.172341E+02, -2.415445E+01},
                                   {4.314185E+02, 1.947381E+02},
                                   {8.064374E+02, 3.987241E+02}};

    auto splines = std::make_shared<Splines>(Projection::cartesian);
    splines->AddSpline(firstSpline);

    std::vector<Point> secondSpline{{2.894012E+01, 2.010146E+02},
                                    {2.344944E+02, 3.720490E+02},
                                    {6.424647E+02, 5.917262E+02}};
    splines->AddSpline(secondSpline);

    std::vector<Point> fourthSpline{{5.067361E+02, 6.034946E+02},
                                    {7.475956E+02, 3.336055E+02}};
    splines->AddSpline(fourthSpline);

    std::vector<Point> thirdSpline{{2.265137E+00, 2.802553E+02},
                                   {2.799988E+02, -2.807726E+01}};
    splines->AddSpline(thirdSpline);

    CurvilinearParameters curvilinearParameters;
    curvilinearParameters.n_refinement = 40;
    curvilinearParameters.m_refinement = 20;
    CurvilinearGridFromSplinesTransfinite curvilinearGridFromSplinesTransfinite(splines, curvilinearParameters);

    const auto curvilinearGrid = curvilinearGridFromSplinesTransfinite.Compute();

    // check the values
    constexpr double tolerance = 1e-6;
    ASSERT_NEAR(244.84733455150598, curvilinearGrid->GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(240.03223719861575, curvilinearGrid->GetNode(0, 1).x, tolerance);
    ASSERT_NEAR(235.21721587684686, curvilinearGrid->GetNode(0, 2).x, tolerance);
    ASSERT_NEAR(230.40187707543339, curvilinearGrid->GetNode(0, 3).x, tolerance);
    ASSERT_NEAR(225.58666038317327, curvilinearGrid->GetNode(0, 4).x, tolerance);
    ASSERT_NEAR(220.77175891290770, curvilinearGrid->GetNode(0, 5).x, tolerance);
    ASSERT_NEAR(215.95654192442103, curvilinearGrid->GetNode(0, 6).x, tolerance);
    ASSERT_NEAR(211.14151904110099, curvilinearGrid->GetNode(0, 7).x, tolerance);
    ASSERT_NEAR(206.32630377949152, curvilinearGrid->GetNode(0, 8).x, tolerance);
    ASSERT_NEAR(201.51108480926104, curvilinearGrid->GetNode(0, 9).x, tolerance);
    ASSERT_NEAR(196.69606411139034, curvilinearGrid->GetNode(0, 10).x, tolerance);

    ASSERT_NEAR(10.946966348412502, curvilinearGrid->GetNode(0, 0).y, tolerance);
    ASSERT_NEAR(16.292559955716278, curvilinearGrid->GetNode(0, 1).y, tolerance);
    ASSERT_NEAR(21.638069155281958, curvilinearGrid->GetNode(0, 2).y, tolerance);
    ASSERT_NEAR(26.983930812344244, curvilinearGrid->GetNode(0, 3).y, tolerance);
    ASSERT_NEAR(32.329656907057121, curvilinearGrid->GetNode(0, 4).y, tolerance);
    ASSERT_NEAR(37.675033050656793, curvilinearGrid->GetNode(0, 5).y, tolerance);
    ASSERT_NEAR(43.020759474232527, curvilinearGrid->GetNode(0, 6).y, tolerance);
    ASSERT_NEAR(48.366270407390992, curvilinearGrid->GetNode(0, 7).y, tolerance);
    ASSERT_NEAR(53.711994913833436, curvilinearGrid->GetNode(0, 8).y, tolerance);
    ASSERT_NEAR(59.057723537488684, curvilinearGrid->GetNode(0, 9).y, tolerance);
    ASSERT_NEAR(64.403232044419184, curvilinearGrid->GetNode(0, 10).y, tolerance);

    ASSERT_NEAR(263.67028430842242, curvilinearGrid->GetNode(1, 0).x, tolerance);
    ASSERT_NEAR(259.11363739326902, curvilinearGrid->GetNode(1, 1).x, tolerance);
    ASSERT_NEAR(254.53691267796933, curvilinearGrid->GetNode(1, 2).x, tolerance);
    ASSERT_NEAR(249.93698634609487, curvilinearGrid->GetNode(1, 3).x, tolerance);
    ASSERT_NEAR(245.31456069699095, curvilinearGrid->GetNode(1, 4).x, tolerance);
    ASSERT_NEAR(240.66785332275725, curvilinearGrid->GetNode(1, 5).x, tolerance);
    ASSERT_NEAR(235.99933187522288, curvilinearGrid->GetNode(1, 6).x, tolerance);
    ASSERT_NEAR(231.30940727936030, curvilinearGrid->GetNode(1, 7).x, tolerance);
    ASSERT_NEAR(226.60252865287427, curvilinearGrid->GetNode(1, 8).x, tolerance);
    ASSERT_NEAR(221.88022520931327, curvilinearGrid->GetNode(1, 9).x, tolerance);
    ASSERT_NEAR(217.14743651601677, curvilinearGrid->GetNode(1, 10).x, tolerance);

    ASSERT_NEAR(34.264668045745267, curvilinearGrid->GetNode(1, 0).y, tolerance);
    ASSERT_NEAR(39.307546170495868, curvilinearGrid->GetNode(1, 1).y, tolerance);
    ASSERT_NEAR(44.379080332661857, curvilinearGrid->GetNode(1, 2).y, tolerance);
    ASSERT_NEAR(49.481517460105827, curvilinearGrid->GetNode(1, 3).y, tolerance);
    ASSERT_NEAR(54.613111211730796, curvilinearGrid->GetNode(1, 4).y, tolerance);
    ASSERT_NEAR(59.775023214376127, curvilinearGrid->GetNode(1, 5).y, tolerance);
    ASSERT_NEAR(64.963841851929189, curvilinearGrid->GetNode(1, 6).y, tolerance);
    ASSERT_NEAR(70.178519042215470, curvilinearGrid->GetNode(1, 7).y, tolerance);
    ASSERT_NEAR(75.413628528250186, curvilinearGrid->GetNode(1, 8).y, tolerance);
    ASSERT_NEAR(80.667056521594716, curvilinearGrid->GetNode(1, 9).y, tolerance);
    ASSERT_NEAR(85.932983124208747, curvilinearGrid->GetNode(1, 10).y, tolerance);
}

TEST(CurvilinearGridFromSplinesTransfinite, FiveSplines)
{
    std::vector<Point> firstSpline{{2.172341E+02, -2.415445E+01},
                                   {4.314185E+02, 1.947381E+02},
                                   {8.064374E+02, 3.987241E+02}};

    auto splines = std::make_shared<Splines>(Projection::cartesian);
    splines->AddSpline(firstSpline);

    std::vector<Point> secondSpline{{2.894012E+01, 2.010146E+02},
                                    {2.344944E+02, 3.720490E+02},
                                    {6.424647E+02, 5.917262E+02}};
    splines->AddSpline(secondSpline);

    std::vector<Point> thirdSpline{{2.265137E+00, 2.802553E+02},
                                   {2.799988E+02, -2.807726E+01}};
    splines->AddSpline(thirdSpline);

    std::vector<Point> fourthSpline{{5.067361E+02, 6.034946E+02},
                                    {7.475956E+02, 3.336055E+02}};
    splines->AddSpline(fourthSpline);

    std::vector<Point> fifthSpline{{2.673223E+02, 4.706788E+02},
                                   {5.513401E+02, 1.545069E+02}};
    splines->AddSpline(fifthSpline);

    CurvilinearParameters curvilinearParameters;
    curvilinearParameters.n_refinement = 40;
    curvilinearParameters.m_refinement = 20;
    CurvilinearGridFromSplinesTransfinite curvilinearGridFromSplinesTransfinite(splines, curvilinearParameters);

    const auto curvilinearGrid = curvilinearGridFromSplinesTransfinite.Compute();

    constexpr double tolerance = 1e-6;
    ASSERT_NEAR(244.84733455150598, curvilinearGrid->GetNode(0, 0).x, tolerance);
    ASSERT_NEAR(240.03223719861575, curvilinearGrid->GetNode(0, 1).x, tolerance);
    ASSERT_NEAR(235.21721587684686, curvilinearGrid->GetNode(0, 2).x, tolerance);
    ASSERT_NEAR(230.40187707543339, curvilinearGrid->GetNode(0, 3).x, tolerance);
    ASSERT_NEAR(225.58666038317327, curvilinearGrid->GetNode(0, 4).x, tolerance);
    ASSERT_NEAR(220.77175891290770, curvilinearGrid->GetNode(0, 5).x, tolerance);
    ASSERT_NEAR(215.95654192442103, curvilinearGrid->GetNode(0, 6).x, tolerance);
    ASSERT_NEAR(211.14151904110099, curvilinearGrid->GetNode(0, 7).x, tolerance);
    ASSERT_NEAR(206.32630377949152, curvilinearGrid->GetNode(0, 8).x, tolerance);
    ASSERT_NEAR(201.51108480926104, curvilinearGrid->GetNode(0, 9).x, tolerance);
    ASSERT_NEAR(196.69606411139034, curvilinearGrid->GetNode(0, 10).x, tolerance);

    ASSERT_NEAR(10.946966348412502, curvilinearGrid->GetNode(0, 0).y, tolerance);
    ASSERT_NEAR(16.292559955716278, curvilinearGrid->GetNode(0, 1).y, tolerance);
    ASSERT_NEAR(21.638069155281958, curvilinearGrid->GetNode(0, 2).y, tolerance);
    ASSERT_NEAR(26.983930812344244, curvilinearGrid->GetNode(0, 3).y, tolerance);
    ASSERT_NEAR(32.329656907057121, curvilinearGrid->GetNode(0, 4).y, tolerance);
    ASSERT_NEAR(37.675033050656793, curvilinearGrid->GetNode(0, 5).y, tolerance);
    ASSERT_NEAR(43.020759474232527, curvilinearGrid->GetNode(0, 6).y, tolerance);
    ASSERT_NEAR(48.366270407390992, curvilinearGrid->GetNode(0, 7).y, tolerance);
    ASSERT_NEAR(53.711994913833436, curvilinearGrid->GetNode(0, 8).y, tolerance);
    ASSERT_NEAR(59.057723537488684, curvilinearGrid->GetNode(0, 9).y, tolerance);
    ASSERT_NEAR(64.403232044419184, curvilinearGrid->GetNode(0, 10).y, tolerance);

    ASSERT_NEAR(255.89614293923407, curvilinearGrid->GetNode(1, 0).x, tolerance);
    ASSERT_NEAR(251.26839070344425, curvilinearGrid->GetNode(1, 1).x, tolerance);
    ASSERT_NEAR(246.62717589518911, curvilinearGrid->GetNode(1, 2).x, tolerance);
    ASSERT_NEAR(241.96945582856105, curvilinearGrid->GetNode(1, 3).x, tolerance);
    ASSERT_NEAR(237.29374836322307, curvilinearGrid->GetNode(1, 4).x, tolerance);
    ASSERT_NEAR(232.59945837385263, curvilinearGrid->GetNode(1, 5).x, tolerance);
    ASSERT_NEAR(227.88656387177011, curvilinearGrid->GetNode(1, 6).x, tolerance);
    ASSERT_NEAR(223.15709488341233, curvilinearGrid->GetNode(1, 7).x, tolerance);
    ASSERT_NEAR(218.41314240105709, curvilinearGrid->GetNode(1, 8).x, tolerance);
    ASSERT_NEAR(213.65762819876193, curvilinearGrid->GetNode(1, 9).x, tolerance);
    ASSERT_NEAR(208.89353710816445, curvilinearGrid->GetNode(1, 10).x, tolerance);

    ASSERT_NEAR(24.731736741118521, curvilinearGrid->GetNode(1, 0).y, tolerance);
    ASSERT_NEAR(29.842940652626876, curvilinearGrid->GetNode(1, 1).y, tolerance);
    ASSERT_NEAR(34.982267945763468, curvilinearGrid->GetNode(1, 2).y, tolerance);
    ASSERT_NEAR(40.148526703963910, curvilinearGrid->GetNode(1, 3).y, tolerance);
    ASSERT_NEAR(45.340177298582923, curvilinearGrid->GetNode(1, 4).y, tolerance);
    ASSERT_NEAR(50.555639961868344, curvilinearGrid->GetNode(1, 5).y, tolerance);
    ASSERT_NEAR(55.793467784299104, curvilinearGrid->GetNode(1, 6).y, tolerance);
    ASSERT_NEAR(61.050433278839293, curvilinearGrid->GetNode(1, 7).y, tolerance);
    ASSERT_NEAR(66.323655962424397, curvilinearGrid->GetNode(1, 8).y, tolerance);
    ASSERT_NEAR(71.609593896396262, curvilinearGrid->GetNode(1, 9).y, tolerance);
    ASSERT_NEAR(76.904826220873304, curvilinearGrid->GetNode(1, 10).y, tolerance);
}

void TestCurvilinearGridFromSplines(const std::vector<meshkernel::Point>& firstSpline,
                                    const std::vector<meshkernel::Point>& secondSpline,
                                    const std::vector<meshkernel::Point>& thirdSpline,
                                    const std::vector<meshkernel::Point>& fourthSpline,
                                    const std::vector<meshkernel::Point>& expectedPoints)
{
    auto splines = std::make_shared<Splines>(Projection::cartesian);

    splines->AddSpline(firstSpline);
    splines->AddSpline(secondSpline);
    splines->AddSpline(thirdSpline);
    splines->AddSpline(fourthSpline);

    CurvilinearParameters curvilinearParameters;
    curvilinearParameters.n_refinement = 4;
    curvilinearParameters.m_refinement = 4;
    CurvilinearGridFromSplinesTransfinite curvilinearGridFromSplinesTransfinite(splines, curvilinearParameters);

    const auto curvilinearGrid = curvilinearGridFromSplinesTransfinite.Compute();

    ASSERT_EQ(curvilinearGrid->m_nodes.size(), expectedPoints.size());

    constexpr double tolerance = 1e-4;

    for (size_t i = 0; i < curvilinearGrid->m_nodes.size(); ++i)
    {
        EXPECT_NEAR(expectedPoints[i].x, curvilinearGrid->m_nodes[i].x, tolerance);
        EXPECT_NEAR(expectedPoints[i].y, curvilinearGrid->m_nodes[i].y, tolerance);
    }
}

TEST(CurvilinearGridFromSplinesTransfinite, FourSplinesSimpleRectangleShortFirstSpline)
{
    std::vector<Point> firstSpline{{0.0, 2.0}, {1.0, 2.0}, {2.0, 2.0}, {3.0, 2.0}, {4.0, 2.0}, {5.0, 2.0}};
    std::vector<Point> secondSpline{{1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0}, {5.0, 0.0}};
    std::vector<Point> thirdSpline{{2.5, -4.0}, {2.5, 4.0}};
    std::vector<Point> fourthSpline{{4.5, 4.0}, {4.5, -4.0}};

    std::vector<Point> expectedPoints{{4.5, 0.0}, {4.0, 0.0}, {3.5, 0.0}, {3.0, 0.0}, {2.5, 0.0}, {4.5, 0.5}, {4.0, 0.5}, {3.5, 0.5}, {3.0, 0.5}, {2.5, 0.5}, {4.5, 1.0}, {4.0, 1.0}, {3.5, 1.0}, {3.0, 1.0}, {2.5, 1.0}, {4.5, 1.5}, {4.0, 1.5}, {3.5, 1.5}, {3.0, 1.5}, {2.5, 1.5}, {4.5, 2.0}, {4.0, 2.0}, {3.5, 2.0}, {3.0, 2.0}, {2.5, 2.0}};

    SCOPED_TRACE("FourSplinesSimpleRectangleShortFirstSpline");
    TestCurvilinearGridFromSplines(thirdSpline, firstSpline, secondSpline, fourthSpline, expectedPoints);
}

TEST(CurvilinearGridFromSplinesTransfinite, FourSplinesSimpleRectangleShortSecondSpline)
{
    std::vector<Point> firstSpline{{0.0, 2.0}, {1.0, 2.0}, {2.0, 2.0}, {3.0, 2.0}, {4.0, 2.0}, {5.0, 2.0}};
    std::vector<Point> secondSpline{{1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0}, {5.0, 0.0}};
    std::vector<Point> thirdSpline{{2.5, -4.0}, {2.5, 4.0}};
    std::vector<Point> fourthSpline{{4.5, 4.0}, {4.5, -4.0}};

    std::vector<Point> expectedPoints{{2.5, 2.0}, {3.0, 2.0}, {3.5, 2.0}, {4.0, 2.0}, {4.5, 2.0}, {2.5, 1.5}, {3.0, 1.5}, {3.5, 1.5}, {4.0, 1.5}, {4.5, 1.5}, {2.5, 1.0}, {3.0, 1.0}, {3.5, 1.0}, {4.0, 1.0}, {4.5, 1.0}, {2.5, 0.5}, {3.0, 0.5}, {3.5, 0.5}, {4.0, 0.5}, {4.5, 0.5}, {2.5, 0.0}, {3.0, 0.0}, {3.5, 0.0}, {4.0, 0.0}, {4.5, 0.0}};

    SCOPED_TRACE("FourSplinesSimpleRectangleShortSecondSpline");
    TestCurvilinearGridFromSplines(fourthSpline, secondSpline, thirdSpline, firstSpline, expectedPoints);
}

TEST(CurvilinearGridFromSplinesTransfinite, FourSplinesSimpleRectangleLongFirstSpline)
{
    std::vector<Point> firstSpline{{0.0, 2.0}, {1.0, 2.0}, {2.0, 2.0}, {3.0, 2.0}, {4.0, 2.0}, {5.0, 2.0}};
    std::vector<Point> secondSpline{{1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0}, {5.0, 0.0}};
    std::vector<Point> thirdSpline{{2.5, -4.0}, {2.5, 4.0}};
    std::vector<Point> fourthSpline{{4.5, 4.0}, {4.5, -4.0}};

    std::vector<Point> expectedPoints{{2.5, 0.0}, {2.5, 0.5}, {2.5, 1.0}, {2.5, 1.5}, {2.5, 2.0}, {3.0, 0.0}, {3.0, 0.5}, {3.0, 1.0}, {3.0, 1.5}, {3.0, 2.0}, {3.5, 0.0}, {3.5, 0.5}, {3.5, 1.0}, {3.5, 1.5}, {3.5, 2.0}, {4.0, 0.0}, {4.0, 0.5}, {4.0, 1.0}, {4.0, 1.5}, {4.0, 2.0}, {4.5, 0.0}, {4.5, 0.5}, {4.5, 1.0}, {4.5, 1.5}, {4.5, 2.0}};

    SCOPED_TRACE("FourSplinesSimpleRectangleLongFirstSpline");
    TestCurvilinearGridFromSplines(firstSpline, secondSpline, thirdSpline, fourthSpline, expectedPoints);
}

TEST(CurvilinearGridFromSplinesTransfinite, FourSplinesSimpleRectangleLongSecondSpline)
{
    std::vector<Point> firstSpline{{1.0, 2.0}, {2.0, 2.0}, {3.0, 2.0}, {4.0, 2.0}, {5.0, 2.0}};
    std::vector<Point> secondSpline{{0.0, 0.0}, {1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0}, {5.0, 0.0}};
    std::vector<Point> thirdSpline{{2.5, -4.0}, {2.5, 4.0}};
    std::vector<Point> fourthSpline{{4.5, 4.0}, {4.5, -4.0}};

    std::vector<Point> expectedPoints{{2.5, 0.0}, {2.5, 0.5}, {2.5, 1.0}, {2.5, 1.5}, {2.5, 2.0}, {3.0, 0.0}, {3.0, 0.5}, {3.0, 1.0}, {3.0, 1.5}, {3.0, 2.0}, {3.5, 0.0}, {3.5, 0.5}, {3.5, 1.0}, {3.5, 1.5}, {3.5, 2.0}, {4.0, 0.0}, {4.0, 0.5}, {4.0, 1.0}, {4.0, 1.5}, {4.0, 2.0}, {4.5, 0.0}, {4.5, 0.5}, {4.5, 1.0}, {4.5, 1.5}, {4.5, 2.0}};

    SCOPED_TRACE("FourSplinesSimpleRectangleLongSecondSpline");
    TestCurvilinearGridFromSplines(firstSpline, secondSpline, thirdSpline, fourthSpline, expectedPoints);
}
