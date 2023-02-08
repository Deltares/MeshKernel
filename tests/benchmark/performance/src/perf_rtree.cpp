#include "memory_management.hpp"

#include <MeshKernel/Entities.hpp>
#include <MeshKernel/RTree.hpp>

#include <benchmark/benchmark.h>

static void BM_RTree(benchmark::State& state)
{
    for (auto _ : state)
    {
        size_t const n = state.range(0); // x
        size_t const m = state.range(0); // y
        // std::cout << "N = " << n << ", M = " << m << '\n';
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

        meshkernel::RTree rtree;

        rtree.BuildTree(nodes);

        for (size_t i = 0; i < nodes.size(); ++i)
        {
            rtree.SearchPoints(nodes[i], 1e-8);
        }
    }
}
BENCHMARK(BM_RTree)->RangeMultiplier(2) //
    ->Args({500, 500})
    ->Args({1000, 1000})
    ->Args({2000, 2000})
    ->Args({4000, 4000})
    ->Args({5000, 5000});
//->Range(1024, 8192);
// BENCHMARK(BM_RTree)->DenseRange(100, 1024, 100);