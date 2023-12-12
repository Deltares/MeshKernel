#include <chrono>
#include <gtest/gtest.h>
#include <random>

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/Mesh2DGenerateGlobalGrid.hpp"
#include "MeshKernel/Polygons.hpp"
#include "TestUtils/MakeMeshes.hpp"

TEST(GlobalGridTest, BasicTest)
{
    const std::vector<meshkernel::Point> polygonNodes{};

    const meshkernel::Polygons polygon(polygonNodes, meshkernel::Projection::sphericalAccurate);

    const auto mesh = meshkernel::Mesh2DGenerateGlobalGrid::Compute(192, 250, polygon);
}
