#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridCreateUniform.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Parameters.hpp>
#include <MeshKernel/Polygons.hpp>

#include <benchmark/benchmark.h>

using namespace meshkernel;

static void BM_CurvilinearUniform(benchmark::State& state)
{
    for (auto _ : state)
    {
        // pause the timers to prepare the benchmark (excludes operation
        // that are irrelevant to the benchmark and should not be measured)
        state.PauseTiming();

        std::vector<meshkernel::Point> polygon_points{
            {2.21527777777778, 5.08143004115226},
            {3.88194444444444, 6.5397633744856},
            {6.27777777777778, 5.8522633744856},
            {6.21527777777778, 4.14393004115226},
            {3.98611111111111, 3.0397633744856},
            {2.38194444444444, 3.39393004115226},
            {2.04861111111111, 4.5397633744856},
            {2.21527777777778, 5.08143004115226}};

        auto const polygons = std::make_shared<Polygons>(polygon_points, Projection::cartesian);

        double const dim_x = 10.0;
        double const dim_y = 15.0;
        double const delta_x = dim_x / static_cast<double>(state.range(0) - 1);
        double const delta_y = dim_y / static_cast<double>(state.range(1) - 1);

        MakeGridParameters make_grid_arameters;
        make_grid_arameters.angle = 0.0;
        make_grid_arameters.origin_x = 0.0;
        make_grid_arameters.origin_y = 0.0;
        make_grid_arameters.num_columns = 3;
        make_grid_arameters.num_rows = 3;
        make_grid_arameters.block_size_x = delta_x;
        make_grid_arameters.block_size_y = delta_y;

        // resume the timers to begin benchmarking
        state.ResumeTiming();

        CurvilinearGridCreateUniform const curvilinear_grid_create_uniform(Projection::cartesian);
        const auto curvilinearGrid = std::make_shared<CurvilinearGrid>(curvilinear_grid_create_uniform.Compute(make_grid_arameters.angle, make_grid_arameters.block_size_x, make_grid_arameters.block_size_y, polygons, 0));
        // const auto curvilinearGrid = std::make_shared<CurvilinearGrid>(curvilinear_grid_create_uniform.Compute(make_grid_arameters, polygons, 0));
    }
}
BENCHMARK(BM_CurvilinearUniform)
    ->ArgNames({"x-nodes", "y-nodes"})
    ->Args({500, 500})
    ->Args({1000, 1000})
    ->Args({2000, 2000})
    ->Args({4000, 4000})
    ->Args({5000, 5000});
