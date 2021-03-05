#include <gtest/gtest.h>

#include <MeshKernel/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGridOrthogonalization.hpp>
#include <MeshKernel/Entities.hpp>

TEST(CurvilinearGridOrthogonalization, Compute_Orthogonalization)
{
    // Set-up
    std::vector<std::vector<meshkernel::Point>> grid{
        {{0, 0}, {0, 10}, {0, 20}, {0, 30}},
        {{10, 0}, {10, 10}, {10, 20}, {10, 30}},
        {{20, 0}, {20, 10}, {20, 20}, {20, 30}},
        {{30, 0}, {30, 10}, {30, 20}, {30, 30}}};

    const auto curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(std::move(grid), meshkernel::Projection::cartesian);

    meshkernelapi::OrthogonalizationParameters orthogonalizationParameters;
    orthogonalizationParameters.OuterIterations = 1;
    orthogonalizationParameters.BoundaryIterations = 25;
    orthogonalizationParameters.InnerIterations = 25;
    orthogonalizationParameters.OrthogonalizationToSmoothingFactor = 0.975;
    meshkernel::CurvilinearGridOrthogonalization curvilinearGridOrthogonalization(curvilinearGrid, orthogonalizationParameters, {0, 0}, {30, 30});

    // Execute
    // curvilinearGridOrthogonalization.Compute();
}