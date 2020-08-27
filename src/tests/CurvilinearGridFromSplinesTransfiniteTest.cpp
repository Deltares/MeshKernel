#include "../Mesh.hpp"
#include "../Entities.hpp"
#include "../Polygons.hpp"
#include "../Constants.cpp"
#include "../CurvilinearGrid.hpp"
#include "../CurvilinearGridFromSplinesTransfinite.cpp"
#include <gtest/gtest.h>

TEST(CurvilinearGridFromSplinesTransfinite, ComputeSplineIntersections)
{
    std::vector<GridGeom::Point> firstSpline;
    firstSpline.push_back(GridGeom::Point{ 2.172341E+02, -2.415445E+01 });
    firstSpline.push_back(GridGeom::Point{ 4.314185E+02, 1.947381E+02 });
    firstSpline.push_back(GridGeom::Point{ 8.064374E+02, 3.987241E+02 });

    GridGeom::Splines splines(GridGeom::Projections::cartesian);
    bool success = splines.AddSpline(firstSpline, 0, firstSpline.size());
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> secondSpline;
    secondSpline.push_back(GridGeom::Point{ 2.894012E+01, 2.010146E+02 });
    secondSpline.push_back(GridGeom::Point{ 2.344944E+02, 3.720490E+02 });
    secondSpline.push_back(GridGeom::Point{ 6.424647E+02, 5.917262E+02 });
    success = splines.AddSpline(secondSpline, 0, secondSpline.size());
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> thirdSpline;
    thirdSpline.push_back(GridGeom::Point{ 2.265137E+00, 2.802553E+02 });
    thirdSpline.push_back(GridGeom::Point{ 2.799988E+02, -2.807726E+01 });
    success = splines.AddSpline(thirdSpline, 0, thirdSpline.size());
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> fourthSpline;
    fourthSpline.push_back(GridGeom::Point{ 5.067361E+02, 6.034946E+02 });
    fourthSpline.push_back(GridGeom::Point{ 7.475956E+02, 3.336055E+02 });
    success = splines.AddSpline(fourthSpline, 0, fourthSpline.size());
    ASSERT_TRUE(success);

    GridGeom::CurvilinearGridFromSplinesTransfinite curvilinearGridFromSplinesTransfinite(&splines);
    GridGeom::CurvilinearGrid curvilinearGrid;
    success = curvilinearGridFromSplinesTransfinite.Compute(curvilinearGrid);
    ASSERT_TRUE(success);
}

TEST(CurvilinearGridFromSplinesTransfinite, ComputeSplineTypeSwappedN)
{
    std::vector<GridGeom::Point> firstSpline;
    firstSpline.push_back(GridGeom::Point{ 2.172341E+02, -2.415445E+01 });
    firstSpline.push_back(GridGeom::Point{ 4.314185E+02, 1.947381E+02 });
    firstSpline.push_back(GridGeom::Point{ 8.064374E+02, 3.987241E+02 });

    GridGeom::Splines splines(GridGeom::Projections::cartesian);
    bool success = splines.AddSpline(firstSpline, 0, firstSpline.size());
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> secondSpline;
    secondSpline.push_back(GridGeom::Point{ 2.894012E+01, 2.010146E+02 });
    secondSpline.push_back(GridGeom::Point{ 2.344944E+02, 3.720490E+02 });
    secondSpline.push_back(GridGeom::Point{ 6.424647E+02, 5.917262E+02 });
    success = splines.AddSpline(secondSpline, 0, secondSpline.size());
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> fourthSpline;
    fourthSpline.push_back(GridGeom::Point{ 5.067361E+02, 6.034946E+02 });
    fourthSpline.push_back(GridGeom::Point{ 7.475956E+02, 3.336055E+02 });
    success = splines.AddSpline(fourthSpline, 0, fourthSpline.size());
    ASSERT_TRUE(success);

    std::vector<GridGeom::Point> thirdSpline;
    thirdSpline.push_back(GridGeom::Point{ 2.265137E+00, 2.802553E+02 });
    thirdSpline.push_back(GridGeom::Point{ 2.799988E+02, -2.807726E+01 });
    success = splines.AddSpline(thirdSpline, 0, thirdSpline.size());
    ASSERT_TRUE(success);

    GridGeom::CurvilinearGridFromSplinesTransfinite curvilinearGridFromSplinesTransfinite(&splines);
    GridGeom::CurvilinearGrid curvilinearGrid;
    success = curvilinearGridFromSplinesTransfinite.Compute(curvilinearGrid);
    ASSERT_TRUE(success);
}