#include "MeshKernel/CurvilinearGrid/CurvilinearGridFromSplinesTransfinite.hpp"

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridCreateUniform.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridFromSplines.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Parameters.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/Splines.hpp>

#include <benchmark/benchmark.h>

using namespace meshkernel;

static void BM_CurvilinearUniform(benchmark::State& state)
{
    for (auto _ : state)
    {
        // pause the timers to prepare the benchmark (excludes operation
        // that are irrelevant to the benchmark and should not be measured)
        state.PauseTiming();

        const std::vector<meshkernel::Point> polygon_points{
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

        MakeGridParameters make_grid_parameters;
        make_grid_parameters.angle = 0.0;
        make_grid_parameters.origin_x = 0.0;
        make_grid_parameters.origin_y = 0.0;
        make_grid_parameters.num_columns = 3;
        make_grid_parameters.num_rows = 3;
        make_grid_parameters.block_size_x = delta_x;
        make_grid_parameters.block_size_y = delta_y;

        // resume the timers to begin benchmarking
        state.ResumeTiming();

        CurvilinearGridCreateUniform const curvilinear_grid_create_uniform(Projection::cartesian);
        const auto curvilinearGrid = std::make_shared<CurvilinearGrid>(
            curvilinear_grid_create_uniform.Compute(make_grid_parameters.angle,
                                                    make_grid_parameters.block_size_x,
                                                    make_grid_parameters.block_size_y,
                                                    polygons,
                                                    0));
    }
}
BENCHMARK(BM_CurvilinearUniform)
    ->ArgNames({"x-nodes", "y-nodes"})
    ->Args({500, 500})
    ->Args({1000, 1000})
    ->Args({2000, 2000})
    ->Args({4000, 4000})
    ->Args({5000, 5000});

static void BM_CurvilinearUniform_add_faces_to_left_boundary(benchmark::State& state)
{
    for (auto _ : state)
    {
        // pause the timers to prepare the benchmark (excludes operation
        // that are irrelevant to the benchmark and should not be measured)
        state.PauseTiming();

        const int numColumns = static_cast<int>(state.range(0));
        const int numRows = static_cast<int>(state.range(1));
        const double block_size = 10.0;

        const double origin_x = 0.0;
        const double origin_y = 0.0;
        const double angle = 0.0;
        const double blockSizeX = block_size;
        const double blockSizeY = block_size;

        CurvilinearGridCreateUniform const curvilinear_grid_create_uniform(Projection::cartesian);
        const auto curvilinearGrid = std::make_shared<CurvilinearGrid>(
            curvilinear_grid_create_uniform.Compute(numColumns,
                                                    numRows,
                                                    origin_x,
                                                    origin_y,
                                                    angle,
                                                    blockSizeX,
                                                    blockSizeY));

        int const faces_to_add = static_cast<int>(state.range(2));

        // resume the timers to begin benchmarking
        state.ResumeTiming();

        Point point{-5.0, 5.0};
        for (int i = 0; i < faces_to_add; ++i)
        {
            curvilinearGrid->InsertFace(point);
            point.x -= block_size;
        }
    }
}
BENCHMARK(BM_CurvilinearUniform_add_faces_to_left_boundary)
    ->ArgNames({"x-nodes", "y-nodes", "faces_to_add"})
    ->Args({500, 500, 10})
    ->Args({1000, 1000, 100});

static void BM_CurvilinearFromSplines(benchmark::State& state)
{
    for (auto _ : state)
    {
        // pause the timers to prepare the benchmark (excludes operation
        // that are irrelevant to the benchmark and should not be measured)
        state.PauseTiming();

        const auto splines = std::make_shared<Splines>(Projection::cartesian);

        const std::vector<Point>
            firstSpline{
                {-429.606973475825, 246.105697581108},
                {-47.0754200262393, 684.481273432093},
                {949.740379838741, 1156.36326272392},
                {2479.86659363708, 1753.89430241889},
                {2870.77475044688, 1977.27039202449}};

        const std::vector<Point> secondSpline{{-549.671621638834, 701.234480152513},
                                              {-83.3740345871489, 248.897898701178}};

        const std::vector<Point> thirdSpline{{422.014368145516, 1265.25910640665},
                                             {731.948692473283, 690.065675672233}};

        const std::vector<Point> fourthSpline{{1502.59620161259, 1801.36172146008},
                                              {1848.82914050127, 1086.55823472217}};

        const std::vector<Point> fifthSpline{{2494.60941555105, 2140.44662548138},
                                             {2802.19829093796, 1489.08194819146}};

        splines->AddSpline(firstSpline);
        splines->AddSpline(secondSpline);
        splines->AddSpline(thirdSpline);
        splines->AddSpline(fourthSpline);
        splines->AddSpline(fifthSpline);

        SplinesToCurvilinearParameters splinesToCurvilinearParameters;
        splinesToCurvilinearParameters.aspect_ratio = 0.1;
        splinesToCurvilinearParameters.aspect_ratio_grow_factor = 1.1;
        splinesToCurvilinearParameters.average_width = 50.0;
        splinesToCurvilinearParameters.nodes_on_top_of_each_other_tolerance = 1e-4;
        splinesToCurvilinearParameters.min_cosine_crossing_angles = 0.95;
        splinesToCurvilinearParameters.check_front_collisions = false;
        splinesToCurvilinearParameters.curvature_adapted_grid_spacing = true;
        splinesToCurvilinearParameters.remove_skinny_triangles = 0;

        CurvilinearParameters curvilinearParameters;
        curvilinearParameters.m_refinement = static_cast<int>(state.range(0));
        curvilinearParameters.n_refinement = static_cast<int>(state.range(1));

        CurvilinearGridFromSplines curvilinearGridFromSplines(splines,
                                                              curvilinearParameters,
                                                              splinesToCurvilinearParameters);

        // resume the timers to begin benchmarking
        state.ResumeTiming();

        const auto curvilinearGrid = curvilinearGridFromSplines.Compute();
    }
}
BENCHMARK(BM_CurvilinearFromSplines)
    ->ArgNames({"m_refinement", "n_refinement"})
    ->Args({20, 40})
    ->Args({200, 400});

static void BM_CurvilinearFromSplinesTransfinite(benchmark::State& state)
{
    for (auto _ : state)
    {
        // pause the timers to prepare the benchmark (excludes operation
        // that are irrelevant to the benchmark and should not be measured)
        state.PauseTiming();

        const std::vector<Point> firstSpline{{2.172341E+02, -2.415445E+01},
                                             {4.314185E+02, 1.947381E+02},
                                             {8.064374E+02, 3.987241E+02}};

        const std::vector<Point> secondSpline{{2.894012E+01, 2.010146E+02},
                                              {2.344944E+02, 3.720490E+02},
                                              {6.424647E+02, 5.917262E+02}};

        const std::vector<Point> thirdSpline{{2.265137E+00, 2.802553E+02},
                                             {2.799988E+02, -2.807726E+01}};

        const std::vector<Point> fourthSpline{{5.067361E+02, 6.034946E+02},
                                              {7.475956E+02, 3.336055E+02}};

        const auto splines = std::make_shared<Splines>(Projection::cartesian);
        splines->AddSpline(firstSpline);
        splines->AddSpline(secondSpline);
        splines->AddSpline(thirdSpline);
        splines->AddSpline(fourthSpline);

        CurvilinearParameters curvilinearParameters;
        curvilinearParameters.m_refinement = static_cast<int>(state.range(0));
        curvilinearParameters.n_refinement = static_cast<int>(state.range(1));

        CurvilinearGridFromSplinesTransfinite curvilinearGridFromSplinesTransfinite(splines, curvilinearParameters);

        // resume the timers to begin benchmarking
        state.ResumeTiming();

        const auto curvilinearGrid = curvilinearGridFromSplinesTransfinite.Compute();
    }
}
BENCHMARK(BM_CurvilinearFromSplinesTransfinite)
    ->ArgNames({"m_refinement", "n_refinement"})
    ->Args({200, 400})
    ->Args({2000, 4000});
