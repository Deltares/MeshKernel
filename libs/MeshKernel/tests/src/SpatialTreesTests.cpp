#include <chrono>
#include <gtest/gtest.h>
#include <random>

#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Utilities/RTreeFactory.hpp>

TEST(RTree, RTreeRemovePoint)
{
    const int n = 4; // x
    const int m = 4; // y

    std::vector<meshkernel::Point> nodes(n * m);
    std::size_t nodeIndex = 0;
    for (auto j = 0; j < m; ++j)
    {
        for (auto i = 0; i < n; ++i)
        {
            nodes[nodeIndex] = {(double)i, (double)j};
            nodeIndex++;
        }
    }

    const auto rtree = meshkernel::RTreeFactory::create(meshkernel::Projection::cartesian);
    rtree->BuildTree(nodes);
    rtree->DeleteNode(0);
    ASSERT_EQ(rtree->Size(), 15);
}

TEST(RTree, PerformanceTestBuildAndSearchRTree)
{
    const int n = 10; // x
    const int m = 10; // y
    std::vector<meshkernel::Point> nodes(n * m);
    std::size_t nodeIndex = 0;
    for (auto j = 0; j < m; ++j)
    {
        for (auto i = 0; i < n; ++i)
        {
            nodes[nodeIndex] = {(double)i, (double)j};
            nodeIndex++;
        }
    }

    auto start(std::chrono::steady_clock::now());
    const auto rtree = meshkernel::RTreeFactory::create(meshkernel::Projection::cartesian);
    rtree->BuildTree(nodes);

    auto end = std::chrono::steady_clock::now();
    double elapsedTime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    std::cout << "Elapsed time build " << nodes.size() << " mesh nodes RTree: " << elapsedTime << " s " << std::endl;

    start = std::chrono::steady_clock::now();
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        rtree->SearchPoints(nodes[i], 1e-8);
        ASSERT_EQ(rtree->GetQueryResultSize(), 1);
    }
    end = std::chrono::steady_clock::now();
    elapsedTime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    std::cout << "Elapsed time search " << nodes.size() << " nodes in RTree: " << elapsedTime << " s " << std::endl;
}

TEST(RTree, SearchPoints_MustComputeCorrectQuerySize)
{
    const int n = 4; // x
    const int m = 4; // y

    std::vector<meshkernel::Point> nodes(n * m);
    std::size_t nodeIndex = 0;
    for (auto j = 0; j < m; ++j)
    {
        for (auto i = 0; i < n; ++i)
        {
            nodes[nodeIndex] = {static_cast<double>(i), static_cast<double>(j)};
            nodeIndex++;
        }
    }

    const auto rtree = meshkernel::RTreeFactory::create(meshkernel::Projection::cartesian);
    rtree->BuildTree(nodes);

    // large search size, node found
    meshkernel::Point const pointToSearch{(n - 1.0) * 0.5, (n - 1.0) * 0.5};
    double squaredDistance = 0.708 * 0.708;
    rtree->SearchPoints(pointToSearch, squaredDistance);
    ASSERT_EQ(rtree->GetQueryResultSize(), 4);

    // smaller search size, node not found
    squaredDistance = 0.700 * 0.700;
    rtree->SearchPoints(pointToSearch, squaredDistance);
    ASSERT_EQ(rtree->GetQueryResultSize(), 0);
}

TEST(RTree, BuildTree_WithBoundigBox_MustCreateaSmallerTree)
{
    // 1 Setup
    const auto boundingBox = meshkernel::BoundingBox({0.0, 0.0}, {10.0, 10.0});

    std::vector<meshkernel::Point> points;
    points.push_back({0.0, 0.0});
    points.push_back({5.5, 5.0});
    points.push_back({10.0, 10.0});
    points.push_back({-10.0, -10.0});
    points.push_back({20.0, 20.0});
    points.push_back({-20.0, 10.0});

    const auto rtree = meshkernel::RTreeFactory::create(meshkernel::Projection::cartesian);

    // 2 Execute
    rtree->BuildTree(points, boundingBox);

    // 3 Assert
    ASSERT_EQ(3, rtree->Size());
}
