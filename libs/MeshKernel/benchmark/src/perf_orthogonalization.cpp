#include <MeshKernel/Constants.hpp>
#include <MeshKernel/LandBoundaries.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/MeshRefinement.hpp>
#include <MeshKernel/OrthogonalizationAndSmoothing.hpp>
#include <MeshKernel/Orthogonalizer.hpp>
#include <MeshKernel/Parameters.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/Smoother.hpp>
#include <TestUtils/Definitions.hpp>
#include <TestUtils/MakeMeshes.hpp>

#include <benchmark/benchmark.h>

#include <cmath>

using namespace meshkernel;

static void BM_Orthogonalization(benchmark::State& state)
{
    for (auto _ : state)
    {
        // pause the timers to prepare the benchmark (excludes operation
        // that are irrelevant to the benchmark and should not be measured)
        state.PauseTiming();

        // create a rectangular grid
        UInt const n = static_cast<UInt>(state.range(0));
        UInt const m = static_cast<UInt>(state.range(1));
        double const dim_x = 10.0;
        double const dim_y = 12.0;
        std::shared_ptr<Mesh2D> mesh =
            MakeRectangularMeshForTesting(n,
                                          m,
                                          dim_x,
                                          dim_y,
                                          Projection::cartesian);

        // move nodes to skew the mesh
        double const delta_x = dim_x / static_cast<double>(n - 1);
        double const delta_y = dim_y / static_cast<double>(m - 1);
        for (UInt i = 0; i < mesh->m_nodes.size(); ++i)
        {
            // only move inetrnal nodes
            if (!mesh->IsNodeOnBoundary(i))
            {
                Point& node = mesh->m_nodes[i];
                double trans_x;
                double trans_y;
                if (i % 2 == 0)
                {
                    trans_x = delta_x / 3.0;
                    trans_y = -delta_y / 4.0;
                }
                else
                {
                    trans_x = -delta_x / 2.0;
                    trans_y = delta_y / 3.0;
                }
                node.x += trans_x;
                node.y += trans_y;
            }
        }

        mesh->AdministrateNodesEdges();

        auto const project_to_land_Boundary = LandBoundaries::ProjectToLandBoundaryOption::DoNotProjectToLandBoundary;
        OrthogonalizationParameters orthogonalization_parameters;
        orthogonalization_parameters.inner_iterations = 25;
        orthogonalization_parameters.boundary_iterations = 25;
        orthogonalization_parameters.outer_iterations = 1;
        orthogonalization_parameters.orthogonalization_to_smoothing_factor = 0.975;
        orthogonalization_parameters.orthogonalization_to_smoothing_factor_at_boundary = 0.975;
        orthogonalization_parameters.areal_to_angle_smoothing_factor = 1.0;

        auto polygon = std::make_shared<Polygons>();
        std::vector<Point> land_boundary{};
        auto landboundaries = std::make_shared<LandBoundaries>(land_boundary, mesh, polygon);

        // resume the timers to begin benchmarking
        state.ResumeTiming();

        auto orthogonalizer = std::make_shared<Orthogonalizer>(mesh);
        auto smoother = std::make_shared<Smoother>(mesh);
        OrthogonalizationAndSmoothing orthogonalization(
            mesh,
            smoother,
            orthogonalizer,
            polygon,
            landboundaries,
            project_to_land_Boundary,
            orthogonalization_parameters);

        orthogonalization.Initialize();

        orthogonalization.Compute();
    }
}
BENCHMARK(BM_Orthogonalization)
    ->ArgNames({"x-nodes", "y-nodes"})
    ->Args({500, 500})
    ->Args({1000, 1000});
