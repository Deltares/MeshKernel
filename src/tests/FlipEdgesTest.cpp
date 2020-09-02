#include "../Entities.hpp"
#include "../Polygons.hpp"
#include "../Constants.cpp"
#include "../Mesh.hpp"
#include "../FlipEdges.cpp"
#include <gtest/gtest.h>

TEST(FlipEdges, FlipEdgesWithLandBoundary)
{
    //1 Setup
    std::vector<GridGeom::Point> nodes;
    nodes.push_back({ 0.0,0.0 });
    nodes.push_back({ 0.0,10.0 });
    nodes.push_back({ 10.0,0.0 });
    nodes.push_back({ 10.0,10.0 });
    std::vector<GridGeom::Edge> edges;
    edges.push_back({ 0, 2 });
    edges.push_back({ 1, 3 });
    edges.push_back({ 0, 1 });
    edges.push_back({ 2, 3 });
    
    //set mesh
    GridGeom::Mesh mesh;
    mesh.Set(edges, nodes, GridGeom::Projections::cartesian);

    //set landboundaries
    GridGeom::Polygons polygon;
    GridGeom::LandBoundaries  landBoundaries;
    std::vector<GridGeom::Point> landBoundary;
    landBoundaries.Set(landBoundary, &mesh, &polygon);

    //execute flipedges
    GridGeom::FlipEdges flipEdges(&mesh, &landBoundaries);

    //auto success = flipEdges.Compute();
    //ASSERT_TRUE(success);

    // check the values
    //constexpr double tolerance = 1e-6;
    //ASSERT_NEAR(244.84733455150598, curvilinearGrid.m_grid[0][0].x, tolerance);

}