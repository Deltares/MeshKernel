#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/MeshRefinement.hpp>
#include <MeshKernel/Parameters.hpp>

#include <TestUtils/Definitions.hpp>
#include <TestUtils/MakeMeshes.hpp>

#include <benchmark/benchmark.h>

using namespace meshkernel;

static void BM_MeshRefinementBasedOnSamples(benchmark::State& state)
{
    for (auto _ : state)
    {
        // pause the timers to prepare the benchmark (excludes operation
        // that are irrelevant to the benchmark and should not be measured)
        state.PauseTiming();

        std::shared_ptr<meshkernel::Mesh2D> mesh =
            MakeRectangularMeshForTesting(static_cast<UInt>(state.range(0)),
                                          static_cast<UInt>(state.range(1)),
                                          10.0,
                                          15.0,
                                          Projection::Cartesian);

        // sample points
        std::vector<Sample> samples;
        for (UInt i = 0; i < mesh->GetNumNodes(); i++)
        {
            if (mesh->m_nodes[i].y > 7.0 && mesh->m_nodes[i].y < 8.0)
            {
                samples.push_back({mesh->m_nodes[i].x, mesh->m_nodes[i].y, 20.0});
            }
        }

        auto const interpolator = std::make_shared<AveragingInterpolation>(
            *mesh,
            samples,
            AveragingInterpolationMethod::Method::MinAbsValue,
            MeshLocation::Faces,
            1.0,
            false,
            false,
            1);

        MeshRefinementParameters mesh_refinement_parameters;
        mesh_refinement_parameters.max_num_refinement_iterations = 1;
        mesh_refinement_parameters.refine_intersected = 0;
        mesh_refinement_parameters.use_mass_center_when_refining = 0;
        mesh_refinement_parameters.min_edge_size = 1.e-5;
        mesh_refinement_parameters.account_for_samples_outside = 0;
        mesh_refinement_parameters.connect_hanging_nodes = 1;
        mesh_refinement_parameters.refinement_type = static_cast<int>(state.range(2));

        // resume the timers to begin benchmarking
        state.ResumeTiming();

        MeshRefinement meshRefinement(mesh, interpolator, mesh_refinement_parameters);

        meshRefinement.Compute();
    }
}

BENCHMARK(BM_MeshRefinementBasedOnSamples)
    ->ArgNames({"x-nodes", "y-nodes", "refinement type"})
    ->Args({500, 500, 1})
    ->Args({500, 500, 2})
    ->Args({1000, 1000, 1})
    ->Args({1000, 1000, 2})
    ->Args({2000, 2000, 1})
    ->Args({2000, 2000, 2});

static void BM_MeshRefinementBasedOnPolygons(benchmark::State& state)
{
    for (auto _ : state)
    {
        // pause the timers to prepare the benchmark (excludes operation
        // that are irrelevant to the benchmark and should not be measured)
        state.PauseTiming();

        std::shared_ptr<meshkernel::Mesh2D> mesh =
            MakeRectangularMeshForTesting(static_cast<UInt>(state.range(0)),
                                          static_cast<UInt>(state.range(1)),
                                          10.0,
                                          15.0,
                                          Projection::Cartesian);

        std::vector<meshkernel::Point> polygon_points{
            {2.21527777777778, 5.08143004115226},
            {3.88194444444444, 6.5397633744856},
            {6.27777777777778, 5.8522633744856},
            {6.21527777777778, 4.14393004115226},
            {3.98611111111111, 3.0397633744856},
            {2.38194444444444, 3.39393004115226},
            {2.04861111111111, 4.5397633744856},
            {2.21527777777778, 5.08143004115226}};

        meshkernel::Polygons polygon(polygon_points, mesh->m_projection);

        MeshRefinementParameters mesh_refinement_parameters;
        mesh_refinement_parameters.max_num_refinement_iterations = 1;
        mesh_refinement_parameters.refine_intersected = 0;
        mesh_refinement_parameters.use_mass_center_when_refining = 0;
        mesh_refinement_parameters.min_edge_size = 1.e-5;
        mesh_refinement_parameters.account_for_samples_outside = 0;
        mesh_refinement_parameters.connect_hanging_nodes = 1;
        mesh_refinement_parameters.refinement_type = 2;

        // resume the timers to begin benchmarking
        state.ResumeTiming();

        MeshRefinement meshRefinement(mesh, polygon, mesh_refinement_parameters);

        meshRefinement.Compute();
    }
}
BENCHMARK(BM_MeshRefinementBasedOnPolygons)
    ->ArgNames({"x-nodes", "y-nodes", "refinement type"})
    ->Args({500, 500, 1})
    ->Args({500, 500, 2})
    ->Args({1000, 1000, 1})
    ->Args({1000, 1000, 2})
    ->Args({2000, 2000, 1})
    ->Args({2000, 2000, 2});