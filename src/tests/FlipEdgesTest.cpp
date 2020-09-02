#include "MakeMeshes.cpp"
#include "../Entities.hpp"
#include "../Polygons.hpp"
#include "../Constants.cpp"
#include "../Mesh.hpp"
#include "../FlipEdges.cpp"

#include <gtest/gtest.h>

TEST(FlipEdges, FlipEdgesWithLandBoundary)
{
    //1 Setup
    auto mesh = MakeRectangularMeshForTesting(3, 3, 10, GridGeom::Projections::cartesian, { 0.0,0.0 });
    
    //set landboundaries
    GridGeom::Polygons polygon;
    
    std::vector<GridGeom::Point> landBoundary{  {-1.369282,21.249086},
                                                {20.885406,21.539995},
                                                {GridGeom::doubleMissingValue,GridGeom::doubleMissingValue} };

    GridGeom::LandBoundaries  landBoundaries;
    landBoundaries.Set(landBoundary, &mesh, &polygon);

    //execute flipedges
    GridGeom::FlipEdges flipEdges(&mesh, &landBoundaries, true, true);

    auto success = flipEdges.Compute();
    ASSERT_TRUE(success);

    // check the values
    //constexpr double tolerance = 1e-6;
    //ASSERT_NEAR(244.84733455150598, curvilinearGrid.m_grid[0][0].x, tolerance);

}