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
            MakeRectangularMeshForTesting(static_cast<int>(state.range(0)),
                                          static_cast<int>(state.range(1)),
                                          10.0,
                                          Projection::cartesian);

        // sample points
        std::vector<Sample> samples{
            {14.7153645, 14.5698833, 1.0},
            {24.7033062, 14.4729137, 1.0},
            {15.5396099, 24.2669525, 1.0},
            {23.8305721, 23.9275551, 1.0}};

        const auto interpolator = std::make_shared<AveragingInterpolation>(mesh,
                                                                           samples,
                                                                           AveragingInterpolation::Method::MinAbsValue,
                                                                           Mesh::Location::Faces,
                                                                           1.0,
                                                                           false,
                                                                           false,
                                                                           1);

        mkapi::MeshRefinementParameters meshRefinementParameters;
        meshRefinementParameters.max_num_refinement_iterations = 1;
        meshRefinementParameters.refine_intersected = 0;
        meshRefinementParameters.use_mass_center_when_refining = 0;
        meshRefinementParameters.min_face_size = 1.e-5;
        meshRefinementParameters.account_for_samples_outside = 0;
        meshRefinementParameters.connect_hanging_nodes = 1;
        meshRefinementParameters.refinement_type = 2;

        MeshRefinement meshRefinement(mesh, interpolator, meshRefinementParameters);

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
