#include "GridGeomTest.hpp"
#include "../Splines.hpp"

TEST(SplineTests, SetSpline)
{
    //One gets the edges
    std::vector<GridGeom::Point> splineNodes;

    splineNodes.push_back(GridGeom::Point{ 322.252624511719,454.880187988281 });
    splineNodes.push_back(GridGeom::Point{ 227.002044677734,360.379241943359 });
    splineNodes.push_back(GridGeom::Point{ 259.252227783203,241.878051757813 });
    splineNodes.push_back(GridGeom::Point{ 428.003295898438,210.377746582031 });
    splineNodes.push_back(GridGeom::Point{ 536.003967285156,310.878753662109 });
    splineNodes.push_back(GridGeom::Point{ 503.753784179688,432.379974365234 });
    splineNodes.push_back(GridGeom::Point{ 350.752807617188,458.630249023438 });
    splineNodes.push_back(GridGeom::Point{ 343.15053976393,406.232256102912 });

    GridGeom::Splines splines;
    bool success = splines.Set(splineNodes);

    ASSERT_TRUE(success);
    EXPECT_EQ(splines.m_numSplines, 1);
    EXPECT_EQ(splines.m_splines[0].size(), splineNodes.size());
    EXPECT_EQ(splines.m_numAllocatedSplines, 5);
    EXPECT_EQ(splines.m_numAllocatedSplineNodes[0], 10);
}

TEST(SplineTests, CubicSplineInterpolation)
{
    //One gets the edges
    std::vector<GridGeom::Point> splineNodes;

    splineNodes.push_back(GridGeom::Point{ 322.252624511719,454.880187988281 });
    splineNodes.push_back(GridGeom::Point{ 227.002044677734,360.379241943359 });
    splineNodes.push_back(GridGeom::Point{ 259.252227783203,241.878051757813 });
    splineNodes.push_back(GridGeom::Point{ 428.003295898438,210.377746582031 });
    splineNodes.push_back(GridGeom::Point{ 536.003967285156,310.878753662109 });
    splineNodes.push_back(GridGeom::Point{ 503.753784179688,432.379974365234 });
    splineNodes.push_back(GridGeom::Point{ 350.752807617188,458.630249023438 });
    splineNodes.push_back(GridGeom::Point{ 343.15053976393,406.232256102912 });

    int number_of_points_between_vertices = 20;
    std::vector<GridGeom::Point> coordinatesDerivatives(splineNodes.size());
    GridGeom::Splines::Derivative(splineNodes, coordinatesDerivatives);
    std::vector<GridGeom::Point> splineCoordinates;

    for (int n = 0; n < splineNodes.size(); n++)
    {
        for (int p = 0; p < number_of_points_between_vertices; p++)
        {
            const double pointAdimensionalCoordinate = n + p / number_of_points_between_vertices;
            GridGeom::Point pointCoordinate;
            GridGeom::Splines::Interpolate(splineNodes, coordinatesDerivatives, pointAdimensionalCoordinate, pointCoordinate);
            splineCoordinates.push_back({ pointCoordinate.x, pointCoordinate.y });
        }
    }

}