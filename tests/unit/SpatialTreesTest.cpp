#pragma once
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/SpatialTrees.hpp>
#include <gtest/gtest.h>
#include <chrono>
#include <random>

TEST(SpatialTrees, RTreeRemovePoint)
{
    const int n = 4; // x
    const int m = 4; // y

    std::vector<meshkernel::Point> nodes(n * m);
    std::size_t nodeIndex = 0;
    for (int j = 0; j < m; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            nodes[nodeIndex] = {(double)i, (double)j};
            nodeIndex++;
        }
    }

    meshkernel::SpatialTrees::RTree rtree;
    rtree.BuildTree(nodes);

    rtree.RemoveNode(0);

    ASSERT_EQ(rtree.Size(), 15);
}

TEST(SpatialTrees, PerformanceTestBuildAndSearchRTree)
{
    const int n = 10; // x
    const int m = 10; // y
    std::vector<meshkernel::Point> nodes(n * m);
    std::size_t nodeIndex = 0;
    for (int j = 0; j < m; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            nodes[nodeIndex] = {(double)i, (double)j};
            nodeIndex++;
        }
    }

    auto start(std::chrono::steady_clock::now());
    meshkernel::SpatialTrees::RTree rtree;
    rtree.BuildTree(nodes);
    auto end = std::chrono::steady_clock::now();
    double elapsedTime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    std::cout << "Elapsed time build " << nodes.size() << " mesh nodes RTree: " << elapsedTime << " s " << std::endl;

    start = std::chrono::steady_clock::now();
    for (int i = 0; i < nodes.size(); ++i)
    {
        rtree.NearestNeighboursOnSquaredDistance(nodes[i], 1e-8);
        ASSERT_EQ(rtree.GetQueryResultSize(), 1);
    }
    end = std::chrono::steady_clock::now();
    elapsedTime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    std::cout << "Elapsed time search " << nodes.size() << " nodes in RTree: " << elapsedTime << " s " << std::endl;
}

TEST(SpatialTrees, FindNodesInSquare)
{
    const int n = 4; // x
    const int m = 4; // y

    std::vector<meshkernel::Point> nodes(n * m);
    std::size_t nodeIndex = 0;
    for (int j = 0; j < m; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            nodes[nodeIndex] = {(double)i, (double)j};
            nodeIndex++;
        }
    }

    meshkernel::SpatialTrees::RTree rtree;
    rtree.BuildTree(nodes);

    // large search size, node found
    std::vector<meshkernel::Point> pointToSearch(1, {(n - 1.0) / 2.0, (n - 1.0) / 2.0});
    double squaredDistance = 0.708 * 0.708;
    rtree.NearestNeighboursOnSquaredDistance(pointToSearch[0], squaredDistance);
    ASSERT_EQ(rtree.GetQueryResultSize(), 4);

    // smaller search size, node not found
    squaredDistance = 0.700 * 0.700;
    rtree.NearestNeighboursOnSquaredDistance(pointToSearch[0], squaredDistance);
    ASSERT_EQ(rtree.GetQueryResultSize(), 0);
}
