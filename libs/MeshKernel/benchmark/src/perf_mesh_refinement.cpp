//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2025.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation version 3.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// contact: delft3d.support@deltares.nl
// Stichting Deltares
// P.O. Box 177
// 2600 MH Delft, The Netherlands
//
// All indications and logos of, and references to, "Delft3D" and "Deltares"
// are registered trademarks of Stichting Deltares, and remain the property of
// Stichting Deltares. All rights reserved.
//
//------------------------------------------------------------------------------

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

        auto mesh = MakeRectangularMeshForTesting(static_cast<UInt>(state.range(0)),
                                                  static_cast<UInt>(state.range(1)),
                                                  10.0,
                                                  15.0,
                                                  Projection::cartesian);

        // sample points
        std::vector<Sample> samples;
        for (UInt i = 0; i < mesh->GetNumNodes(); i++)
        {
            if (mesh->Node(i).y > 7.0 && mesh->Node(i).y < 8.0)
            {
                samples.push_back({mesh->Node(i).x, mesh->Node(i).y, 20.0});
            }
        }

        auto interpolator = std::make_unique<AveragingInterpolation>(*mesh,
                                                                     samples,
                                                                     AveragingInterpolation::Method::MinAbsValue,
                                                                     Location::Faces,
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

        MeshRefinement meshRefinement(*mesh,
                                      std::move(interpolator),
                                      mesh_refinement_parameters);

        [[maybe_unused]] auto dummyUndoAction = meshRefinement.Compute();
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
                                          Projection::cartesian);

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

        MeshRefinement meshRefinement(*mesh, polygon, mesh_refinement_parameters);

        [[maybe_unused]] auto dummyUndoAction = meshRefinement.Compute();
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
