//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2024.
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

#include <numbers>
#include <random>
#include <utility>

#include <gtest/gtest.h>

#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridGenerateCircularGrid.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Operations.hpp"
#include "MeshKernel/Parameters.hpp"
#include "MeshKernel/UndoActions/UndoActionStack.hpp"
#include "MeshKernel/Utilities/LinearAlgebra.hpp"
#include "TestUtils/MakeCurvilinearGrids.hpp"

namespace mk = meshkernel;

TEST(CurvilinearGridGenerationTests, Basic)
{

    mk::CurvilinearGridGenerateCircularGrid generate;

    mk::MakeGridParameters parameters = {
        .num_columns = 10,
        .num_rows = 10,
        .left_rotation = 0.0,
        .column_curvature_radius = 10.0,
        .fraction_columns = 0.5,
        .fraction_rows = 0.5,
        .maximum_uniform_columns_size = 10.0,
        .maximum_uniform_rows_size = 10.0};

    auto grid = generate.Compute(parameters);

    mk::Print(grid.ComputeNodes(), grid.ComputeEdges());
}
