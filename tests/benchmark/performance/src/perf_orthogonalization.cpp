#include "custom_memory_manager.hpp"

#include <MeshKernel/LandBoundaries.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/MeshRefinement.hpp>
#include <MeshKernel/OrthogonalizationAndSmoothing.hpp>
#include <MeshKernel/Orthogonalizer.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/Smoother.hpp>
#include <MeshKernelApi/MeshRefinementParameters.hpp>
#include <TestUtils/Definitions.hpp>
#include <TestUtils/MakeMeshes.hpp>

#include <benchmark/benchmark.h>

#include <cmath>

using namespace meshkernel;
namespace mkapi = meshkernelapi;

static void BM_Orthogonalization(benchmark::State& state)
{
    for (auto _ : state)
    {
        CUSTOM_MEMORY_MANAGER.ResetStatistics();

        // create a rectangular grid
        double const delta = 10.0;
        std::shared_ptr<Mesh2D> mesh =
            MakeRectangularMeshForTesting(static_cast<int>(state.range(0)),
                                          static_cast<int>(state.range(1)),
                                          delta,
                                          Projection::cartesian);

        // move nodes to skew the mesh
        for (size_t i = 0; i < mesh->m_nodes.size(); ++i)
        {
            // only move inetrnal nodes
            if (!mesh->IsNodeOnBoundary(i))
            {
                Point& node = mesh->m_nodes[i];
                double delta_x;
                double delta_y;
                if (i % 2 == 0)
                {
                    delta_x = delta / 3.0;
                    delta_y = -delta / 4.0;
                }
                else
                {
                    delta_x = -delta / 2.0;
                    delta_y = delta / 3.0;
                }
                node.x += delta_x;
                node.y += delta_y;
            }
        }

        mesh->AdministrateNodesEdges();

        auto const projectToLandBoundaryOption = LandBoundaries::ProjectToLandBoundaryOption::DoNotProjectToLandBoundary;
        mkapi::OrthogonalizationParameters orthogonalizationParameters;
        orthogonalizationParameters.inner_iterations = 25;
        orthogonalizationParameters.boundary_iterations = 25;
        orthogonalizationParameters.outer_iterations = 1;
        orthogonalizationParameters.orthogonalization_to_smoothing_factor = 0.975;
        orthogonalizationParameters.orthogonalization_to_smoothing_factor_at_boundary = 0.975;
        orthogonalizationParameters.areal_to_angle_smoothing_factor = 1.0;

        auto polygon = std::make_shared<Polygons>();
        std::vector<Point> landBoundary{};
        auto landboundaries = std::make_shared<LandBoundaries>(landBoundary, mesh, polygon);

        auto orthogonalizer = std::make_shared<Orthogonalizer>(mesh);
        auto smoother = std::make_shared<Smoother>(mesh);
        OrthogonalizationAndSmoothing orthogonalization(mesh,
                                                        smoother,
                                                        orthogonalizer,
                                                        polygon,
                                                        landboundaries,
                                                        projectToLandBoundaryOption,
                                                        orthogonalizationParameters);

        orthogonalization.Initialize();

        orthogonalization.Compute();
    }
}
BENCHMARK(BM_Orthogonalization)
    ->ArgNames({"x-nodes", "y-nodes"})
    ->Args({500, 500})
    ->Args({1000, 1000})
    ->Args({2000, 2000})
    ->Args({4000, 4000})
    ->Args({5000, 5000});
