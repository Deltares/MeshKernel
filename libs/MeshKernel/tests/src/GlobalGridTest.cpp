#include <chrono>
#include <gtest/gtest.h>
#include <random>

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/GenerateGlobalGrid.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/Polygons.hpp"
#include "TestUtils/Definitions.hpp"
#include "TestUtils/MakeMeshes.hpp"

namespace mk = meshkernel;

TEST(GlobalGridTest, BasicTest)
{

    mk::Mesh2D mesh;

    mk::GenerateGlobalGrid::Compute(192, 250, mesh);
}
