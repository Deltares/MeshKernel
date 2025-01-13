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

#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Utilities/RTree.hpp>

#include <benchmark/benchmark.h>

static void BM_RTree(benchmark::State& state)
{
    for (auto _ : state)
    {
        // pause the timers to prepare the benchmark (excludes operation
        // that are irrelevant to the benchmark and should not be measured)
        state.PauseTiming();

        int64_t const n = state.range(0); // number of nodes in x-dir
        int64_t const m = state.range(1); // number of nodes in y-dir
        std::vector<meshkernel::Point> nodes(n * m);
        std::size_t nodeIndex = 0;
        for (int64_t j = 0; j < m; ++j)
        {
            for (int64_t i = 0; i < n; ++i)
            {
                nodes[nodeIndex] = {static_cast<double>(i), static_cast<double>(j)};
                nodeIndex++;
            }
        }

        // resume the timers to begin benchmarking
        state.ResumeTiming();

        meshkernel::RTree rtree;

        rtree.BuildTree(nodes);

        rtree.SearchPoints(nodes[0], 1.0e-8);
    }
}
BENCHMARK(BM_RTree)
    ->ArgNames({"x-nodes", "y-nodes"})
    ->Args({500, 500})
    ->Args({1000, 1000})
    ->Args({2000, 2000})
    ->Args({4000, 4000})
    ->Args({5000, 5000});
