#include <MeshKernel/Entities.hpp>
#include <MeshKernel/RTree.hpp>

#include <TestUtils/MakeMeshes.hpp>

#include <benchmark/benchmark.h>

using namespace meshkernel;

template <typename T>
static void BM_RTree(benchmark::State& state)
{
    for (auto _ : state)
    {
        // pause the timers to prepare the benchmark (excludes operation
        // that are irrelevant to the benchmark and should not be measured)
        state.PauseTiming();

        size_t const n = static_cast<size_t>(state.range(0)); // number of nodes in x-dir
        size_t const m = static_cast<size_t>(state.range(1)); // number of nodes in y-dir
        double const delta = 1.0;
        auto const mesh = MakeRectangularMeshForTesting(m, n, delta, Projection::cartesian);
        std::vector<Point> const& nodes = mesh->m_nodes;

        // resume the timers to begin benchmarking
        state.ResumeTiming();

        RTree<T> rtree;

        rtree.BuildTree(nodes);

        {
            double const search_radius = 2.0 * delta;
            double const search_radius_squared = search_radius * search_radius;
            std::size_t nodeIndex = 0;
            for (size_t j = 0; j < m; ++j)
            {
                for (size_t i = 0; i < n; ++i)
                {
                    rtree.SearchPoints(nodes[nodeIndex], search_radius_squared);
                    nodeIndex++;
                }
            }
        }
    }
}

// linear tests

BENCHMARK(BM_RTree<bgi::linear<16>>)
    ->ArgNames({"x-nodes", "y-nodes"})
    ->Args({500, 500})
    ->Args({1000, 1000})
    ->Args({2000, 2000})
    ->Args({4000, 4000})
    ->Args({5000, 5000});

BENCHMARK(BM_RTree<bgi::linear<32>>)
    ->ArgNames({"x-nodes", "y-nodes"})
    ->Args({500, 500})
    ->Args({1000, 1000})
    ->Args({2000, 2000})
    ->Args({4000, 4000})
    ->Args({5000, 5000});

BENCHMARK(BM_RTree<bgi::linear<64>>)
    ->ArgNames({"x-nodes", "y-nodes"})
    ->Args({500, 500})
    ->Args({1000, 1000})
    ->Args({2000, 2000})
    ->Args({4000, 4000})
    ->Args({5000, 5000});

// quadratic tests

BENCHMARK(BM_RTree<bgi::quadratic<16>>)
    ->ArgNames({"x-nodes", "y-nodes"})
    ->Args({500, 500})
    ->Args({1000, 1000})
    ->Args({2000, 2000})
    ->Args({4000, 4000})
    ->Args({5000, 5000});

BENCHMARK(BM_RTree<bgi::quadratic<32>>)
    ->ArgNames({"x-nodes", "y-nodes"})
    ->Args({500, 500})
    ->Args({1000, 1000})
    ->Args({2000, 2000})
    ->Args({4000, 4000})
    ->Args({5000, 5000});

BENCHMARK(BM_RTree<bgi::quadratic<64>>)
    ->ArgNames({"x-nodes", "y-nodes"})
    ->Args({500, 500})
    ->Args({1000, 1000})
    ->Args({2000, 2000})
    ->Args({4000, 4000})
    ->Args({5000, 5000});

// rstar tets

BENCHMARK(BM_RTree<bgi::rstar<16>>)
    ->ArgNames({"x-nodes", "y-nodes"})
    ->Args({500, 500})
    ->Args({1000, 1000})
    ->Args({2000, 2000})
    ->Args({4000, 4000})
    ->Args({5000, 5000});

BENCHMARK(BM_RTree<bgi::rstar<32>>)
    ->ArgNames({"x-nodes", "y-nodes"})
    ->Args({500, 500})
    ->Args({1000, 1000})
    ->Args({2000, 2000})
    ->Args({4000, 4000})
    ->Args({5000, 5000});

BENCHMARK(BM_RTree<bgi::rstar<64>>)
    ->ArgNames({"x-nodes", "y-nodes"})
    ->Args({500, 500})
    ->Args({1000, 1000})
    ->Args({2000, 2000})
    ->Args({4000, 4000})
    ->Args({5000, 5000});