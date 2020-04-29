# pragma once
#include "../Entities.hpp"
#include "../SpatialTrees.hpp"
#include <gtest/gtest.h>
#include <chrono>
#include <random>

TEST(SpatialTrees, RTreeOnePoint)
{
    const int n = 4; // x
    const int m = 4; // y

    std::vector<GridGeom::Point> nodes(n * m);
    std::size_t nodeIndex = 0;
    for (int j = 0; j < m; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            nodes[nodeIndex] = { (double)i, (double)j };
            nodeIndex++;
        }
    }

    GridGeom::SpatialTrees::RTree rtree;
    rtree.BuildTree(nodes, GridGeom::Projections::cartesian);
    std::vector<GridGeom::Point> pointToSearch(1, { (n-1.0)/2.0, (n-1.0)/2.0 });
    auto successful = rtree.NearestNeighbours(pointToSearch[0], 0.708);
    ASSERT_EQ(true, successful);
    ASSERT_EQ(rtree.GetQueryResultSize(),4);
    successful  = rtree.NearestNeighbours(pointToSearch[0], 0.700);
    ASSERT_EQ(rtree.GetQueryResultSize(), 0);
}

TEST(SpatialTrees, RTreeRemovePoint)
{
    const int n = 4; // x
    const int m = 4; // y

    std::vector<GridGeom::Point> nodes(n * m);
    std::size_t nodeIndex = 0;
    for (int j = 0; j < m; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            nodes[nodeIndex] = { (double)i, (double)j };
            nodeIndex++;
        }
    }

    GridGeom::SpatialTrees::RTree rtree;
    rtree.BuildTree(nodes, GridGeom::Projections::cartesian);

    rtree.RemoveNode(0);

    ASSERT_EQ(rtree.Size(), 15);
}

TEST(SpatialTrees, RTreeManyPoints)
{
    const int n = 10; // x
    const int m = 10; // y
    std::vector<GridGeom::Point> nodes(n * m );
    std::size_t nodeIndex = 0;
    for (int j = 0; j < m; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            nodes[nodeIndex] = { (double)i, (double)j };
            nodeIndex++;
        }
    }

    auto start(std::chrono::steady_clock::now());
    GridGeom::SpatialTrees::RTree rtree;
    rtree.BuildTree(nodes, GridGeom::Projections::cartesian);
    auto end = std::chrono::steady_clock::now();
    double elapsedTime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    std::cout << "Elapsed time build " << elapsedTime << " s " << std::endl;

    //generate random points
    //std::uniform_real_distribution<double>  xDistrution(0.0, double(n - 1));
    //std::uniform_real_distribution<double>  yDistrution(0.0, double(m - 1));
    //std::random_device                  rand_dev;
    //std::mt19937                        generator(rand_dev());

    //const int numPoints = m*n;
    //std::vector<GridGeom::Point> pointsToSearch(numPoints);
    //for (int i = 0; i < numPoints; ++i)
    //{
    //    pointsToSearch[i].x = xDistrution(generator);
    //    pointsToSearch[i].y = yDistrution(generator);
    //}

    start = std::chrono::steady_clock::now();
    for (int i = 0; i < nodes.size(); ++i)
    {
        auto successful = rtree.NearestNeighbours(nodes[i], 1e-4);
        ASSERT_EQ(successful, true);
        ASSERT_EQ(rtree.GetQueryResultSize(), 1);
    }
    end = std::chrono::steady_clock::now();
    elapsedTime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    std::cout << "Elapsed time search " << elapsedTime << " s " << std::endl;
}
