#include "MakeCurvilinearGrids.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"

#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>

using namespace meshkernel;

TEST(CurvilinearGridBoundaryPolygon, ComputeBoundaryToPolygon_OnValidCurvilinearGrid_ShouldComputeBoundaryPolygon)

{
    // Prepare
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 1.0;

    constexpr size_t nx = 3;
    constexpr size_t ny = 3;

    std::unique_ptr<CurvilinearGrid> grid = MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);

    // Execute
    const auto boundaryPolygon = grid->ComputeBoundaryToPolygon();

    // Assert
    std::vector<Point> expectedBoundaryPolygon{{0.0, 0.0},
                                               {1.0, 0.0},
                                               {2.0, 0.0},
                                               {2.0, 1.0},
                                               {2.0, 2.0},
                                               {1.0, 2.0},
                                               {0.0, 2.0},
                                               {0.0, 1.0},
                                               {0.0, 0.0}};
    ASSERT_THAT(boundaryPolygon, ::testing::ContainerEq(expectedBoundaryPolygon));
}

TEST(CurvilinearGridBoundaryPolygon, ComputeBoundaryToPolygon_OnValidCurvilinearGridWithInvalidNodes_ShouldComputeBoundaryPolygon)

{
    // Prepare
    constexpr double originX = 0.0;
    constexpr double originY = 0.0;

    constexpr double deltaX = 1.0;
    constexpr double deltaY = 1.0;

    constexpr size_t nx = 4;
    constexpr size_t ny = 4;

    std::unique_ptr<CurvilinearGrid> grid = MakeCurvilinearGrid(originX, originY, deltaX, deltaY, nx, ny);
    auto undoAction = grid->DeleteNode({3.0, 0.0});
    undoAction = grid->DeleteNode({3.0, 3.0});

    // Execute
    const auto boundaryPolygon = grid->ComputeBoundaryToPolygon();

    // Assert
    std::vector<Point> expectedBoundaryPolygon{{0.0, 0.0},
                                               {1.0, 0.0},
                                               {2.0, 0.0},
                                               {2.0, 1.0},
                                               {3.0, 1.0},
                                               {3.0, 2.0},
                                               {2.0, 2.0},
                                               {2.0, 3.0},
                                               {1.0, 3.0},
                                               {0.0, 3.0},
                                               {0.0, 2.0},
                                               {0.0, 1.0},
                                               {0.0, 0.0}};

    ASSERT_THAT(boundaryPolygon, ::testing::ContainerEq(expectedBoundaryPolygon));
}
