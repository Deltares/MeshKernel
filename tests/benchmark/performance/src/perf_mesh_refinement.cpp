#include "custom_memory_manager.hpp"

#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/MeshRefinement.hpp>
#include <MeshKernelApi/MeshRefinementParameters.hpp>
#include <TestUtils/Definitions.hpp>
#include <TestUtils/MakeMeshes.hpp>

#include <benchmark/benchmark.h>

using namespace meshkernel;
namespace mkapi = meshkernelapi;

static void BM_MeshRefinement(benchmark::State& state)
{
    for (auto _ : state)
    {
        CUSTOM_MEMORY_MANAGER.ResetStatistics();
        std::shared_ptr<meshkernel::Mesh2D> mesh =
            MakeRectangularMeshForTesting(static_cast<size_t>(state.range(0)),
                                          static_cast<size_t>(state.range(1)),
                                          10.0,
                                          15.0,
                                          Projection::cartesian);

        // sample points
        double const dummy_sample_value = 0.0;
        std::vector<Sample> samples{
            {14.7153645, 14.5698833, dummy_sample_value},
            {24.7033062, 14.4729137, dummy_sample_value},
            {15.5396099, 24.2669525, dummy_sample_value},
            {23.8305721, 23.9275551, dummy_sample_value}};

        auto const interpolator = std::make_shared<AveragingInterpolation>(
            mesh,
            samples,
            AveragingInterpolation::Method::MinAbsValue,
            Mesh::Location::Faces,
            1.0,
            false,
            false,
            1);

        mkapi::MeshRefinementParameters mesh_refinement_parameters;
        mesh_refinement_parameters.max_num_refinement_iterations = 1;
        mesh_refinement_parameters.refine_intersected = 0;
        mesh_refinement_parameters.use_mass_center_when_refining = 0;
        mesh_refinement_parameters.min_face_size = 1.e-5;
        mesh_refinement_parameters.account_for_samples_outside = 0;
        mesh_refinement_parameters.connect_hanging_nodes = 1;
        mesh_refinement_parameters.refinement_type = 2;

        MeshRefinement meshRefinement(mesh, interpolator, mesh_refinement_parameters);

        meshRefinement.Compute();
    }
}
BENCHMARK(BM_MeshRefinement)
    ->ArgNames({"x-nodes", "y-nodes"})
    ->Args({500, 500})
    ->Args({1000, 1000})
    ->Args({2000, 2000})
    ->Args({4000, 4000})
    ->Args({5000, 5000});
