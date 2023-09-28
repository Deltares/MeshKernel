#include "CartesianApiTestFixture.hpp"

#include "MeshKernelApi/MeshKernel.hpp"

#include <memory>

#include <gtest/gtest.h>

TEST_F(CartesianApiTestFixture, RangeError)
{
    this->MakeUniformCurvilinearGrid(3, 3);
    int const meshKernelId = this->GetMeshKernelId();
    meshkernel::OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.outer_iterations = -1; // invalid, must be > 0
    int const errorCode = meshkernelapi::mkernel_curvilinear_initialize_orthogonalize(meshKernelId,
                                                                                      orthogonalizationParameters);
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
