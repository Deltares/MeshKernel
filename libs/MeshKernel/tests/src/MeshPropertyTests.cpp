#include <chrono>
#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>
#include <random>
#include <span>

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Mesh.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/Mesh2DIntersections.hpp"
#include "MeshKernel/Mesh2DToCurvilinear.hpp"
#include "MeshKernel/MeshTriangulation.hpp"
#include "MeshKernel/Operations.hpp"
#include "MeshKernel/Polygons.hpp"
#include "MeshKernel/RemoveDisconnectedRegions.hpp"
#include "MeshKernel/SampleInterpolator.hpp"

#include "TestUtils/Definitions.hpp"
#include "TestUtils/MakeMeshes.hpp"

namespace mk = meshkernel;

TEST(MeshPropertyTests, SimpleSampleInterpolationTest)
{
    const double tolerance = 1.0e-13;

    std::vector<double> xValues{0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0};
    std::vector<double> yValues{0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0};

    std::vector<mk::Point> nodes{{0.25, 0.25}, {1.25, 0.25}, {2.25, 0.25}, {0.25, 1.25}, {1.25, 1.25}, {2.25, 1.25}, {0.25, 2.25}, {1.25, 2.25}, {2.25, 2.25}};
    std::vector<double> interpolated(nodes.size(), 0.0);

    mk::SampleInterpolator properties(xValues, yValues, mk::Projection::cartesian);

    std::vector<double> expectedXDepths{0.25, 1.25, 2.25, 0.25, 1.25, 2.25, 0.25, 1.25, 2.25};
    std::vector<double> expectedYDepths{0.25, 0.25, 0.25, 1.25, 1.25, 1.25, 2.25, 2.25, 2.25};

    properties.SetData("xdepth", std::vector{0.0, 1.0, 2.0, 3.0,
                                             0.0, 1.0, 2.0, 3.0,
                                             0.0, 1.0, 2.0, 3.0,
                                             0.0, 1.0, 2.0, 3.0});

    properties.SetData("ydepth", std::vector{0.0, 0.0, 0.0, 0.0,
                                             1.0, 1.0, 1.0, 1.0,
                                             2.0, 2.0, 2.0, 2.0,
                                             3.0, 3.0, 3.0, 3.0});

    properties.Interpolate("xdepth", nodes, interpolated);

    for (size_t i = 0; i < interpolated.size(); ++i)
    {
        EXPECT_NEAR(interpolated[i], expectedXDepths[i], tolerance);
    }

    std::span interpolatedData(interpolated.data(), interpolated.size());

    properties.Interpolate("ydepth", nodes, interpolatedData);

    for (size_t i = 0; i < interpolated.size(); ++i)
    {
        EXPECT_NEAR(interpolatedData[i], expectedYDepths[i], tolerance);
    }
}

TEST(MeshPropertyTests, AveragePointTest)
{
    const double tolerance = 1.0e-13;

    std::vector<mk::Point> pnts{{1.0, 2.0}, {-999.0, -999.0}, {1.0, 2.0}, {1.0, 2.0}};
    mk::Point average = mk::ComputeAverageCoordinate(pnts, mk::Projection::sphericalAccurate);

    EXPECT_NEAR(average.x, 1.0, tolerance);
    EXPECT_NEAR(average.y, 2.0, tolerance);
}
