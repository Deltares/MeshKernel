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

#include "CartesianApiTestFixture.hpp"

#include "MeshKernelApi/MeshKernel.hpp"

#include <memory>

#include <gtest/gtest.h>

TEST_F(CartesianApiTestFixture, RangeError)
{
    this->MakeRectangularCurvilinearGrid(3, 3);
    int const meshKernelId = this->GetMeshKernelId();
    meshkernel::OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = -1; // invalid, must be > 0
    int const errorCode = meshkernelapi::mkernel_curvilinear_orthogonalize(meshKernelId, orthogonalizationParameters, 0.0, 0.0, 10.0, 10.0);
    EXPECT_EQ(meshkernel::ExitCode::RangeErrorCode, errorCode);
}

TEST(ErrorHandling, ExitCodes)
{
    int exitCode;

    meshkernelapi::mkernel_get_exit_code_success(exitCode);
    EXPECT_EQ(meshkernel::ExitCode::Success, exitCode);

    meshkernelapi::mkernel_get_exit_code_meshkernel_error(exitCode);
    EXPECT_EQ(meshkernel::ExitCode::MeshKernelErrorCode, exitCode);

    meshkernelapi::mkernel_get_exit_code_not_implemented_error(exitCode);
    EXPECT_EQ(meshkernel::ExitCode::NotImplementedErrorCode, exitCode);

    meshkernelapi::mkernel_get_exit_code_algorithm_error(exitCode);
    EXPECT_EQ(meshkernel::ExitCode::AlgorithmErrorCode, exitCode);

    meshkernelapi::mkernel_get_exit_code_constraint_error(exitCode);
    EXPECT_EQ(meshkernel::ExitCode::ConstraintErrorCode, exitCode);

    meshkernelapi::mkernel_get_exit_code_mesh_geometry_error(exitCode);
    EXPECT_EQ(meshkernel::ExitCode::MeshGeometryErrorCode, exitCode);

    meshkernelapi::mkernel_get_exit_code_linear_algebra_error(exitCode);
    EXPECT_EQ(meshkernel::ExitCode::LinearAlgebraErrorCode, exitCode);

    meshkernelapi::mkernel_get_exit_code_range_error(exitCode);
    EXPECT_EQ(meshkernel::ExitCode::RangeErrorCode, exitCode);

    meshkernelapi::mkernel_get_exit_code_stdlib_exception(exitCode);
    EXPECT_EQ(meshkernel::ExitCode::StdLibExceptionCode, exitCode);

    meshkernelapi::mkernel_get_exit_code_unknown_exception(exitCode);
    EXPECT_EQ(meshkernel::ExitCode::UnknownExceptionCode, exitCode);
}
