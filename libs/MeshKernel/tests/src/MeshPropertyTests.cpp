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

TEST(MeshPropertyTests, TriangulationTest)
{

    std::vector<double> xValues{0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0};
    std::vector<double> yValues{0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0};
    // std::vector<double> xValues{0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0};
    // std::vector<double> yValues{0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0};

    std::vector<mk::Point> nodes{{0.25, 0.25}, {1.25, 0.25}, {2.25, 0.25}, {0.25, 1.25}, {1.25, 1.25}, {2.25, 1.25}, {0.25, 2.25}, {1.25, 2.25}, {2.25, 2.25}};
    std::vector<double> interpolated(nodes.size(), 0.0);

    mk::MeshTriangulation triangulation(xValues, yValues);

    // std::cout << "nullId = 4294967295;" << std::endl;
    // std::cout << "nullValue = -999.0;" << std::endl;
    // std::cout << "nodex = zeros (" << triangulation.NumberOfNodes() << ", 1);" << std::endl;
    // std::cout << "nodey = zeros (" << triangulation.NumberOfNodes() << ", 1);" << std::endl;
    // std::cout << "edges = zeros (" << triangulation.NumberOfEdges() << ", 2);" << std::endl;

    // for (mk::UInt i = 0; i < triangulation.NumberOfNodes(); ++i)
    // {
    //     std::cout << "nodex (" << i + 1 << ") = " << triangulation.GetNode(i).x << std::endl;
    // }

    // for (mk::UInt i = 0; i < triangulation.NumberOfNodes(); ++i)
    // {
    //     std::cout << "nodey (" << i + 1 << ") = " << triangulation.GetNode(i).y << std::endl;
    // }

    // for (mk::UInt i = 0; i < triangulation.NumberOfEdges(); ++i)
    // {
    //     std::cout << "edge (" << i + 1 << ", 1) = " << triangulation.GetEdge(i).first << std::endl;
    //     std::cout << "edge (" << i + 1 << ", 2) = " << triangulation.GetEdge(i).second << std::endl;
    // }

    // return;

    mk::SampleInterpolator properties(xValues, yValues);

    properties.SetData("xdepth", std::vector{0.0, 1.0, 2.0, 3.0,
                                             0.0, 1.0, 2.0, 3.0,
                                             0.0, 1.0, 2.0, 3.0,
                                             0.0, 1.0, 2.0, 3.0});

    properties.SetData("ydepth", std::vector{0.0, 0.0, 0.0, 0.0,
                                             1.0, 1.0, 1.0, 1.0,
                                             2.0, 2.0, 2.0, 2.0,
                                             3.0, 3.0, 3.0, 3.0});

    // properties.SetData("xdepth", std::vector{0.0, 1.0, 2.0,
    //                                          0.0, 1.0, 2.0,
    //                                          0.0, 1.0, 2.0});

    // properties.SetData("ydepth", std::vector{0.0, 0.0, 0.0,
    //                                          1.0, 1.0, 1.0,
    //                                          2.0, 2.0, 2.0});

    properties.Interpolate("xdepth", nodes, interpolated);

    for (size_t i = 0; i < nodes.size(); ++i)
    {
        std::cout << "xdepth interpolated " << interpolated[i] << std::endl;
    }

    std::span interpolatedData(interpolated.data(), interpolated.size());

    properties.Interpolate("ydepth", nodes, interpolatedData);

    for (size_t i = 0; i < nodes.size(); ++i)
    {
        std::cout << "ydepth interpolated " << interpolatedData[i] << std::endl;
    }
}

TEST(MeshPropertyTests, AveragePointTest)
{
    std::vector<mk::Point> pnts{{1.0, 1.0}, {-999.0, -999.0}, {1.0, 1.0}, {1.0, 1.0}};
    mk::Point average = mk::ComputeAverageCoordinate(pnts, mk::Projection::sphericalAccurate);

    std::cout << "average: " << average.x << ", " << average.y << std::endl;
}
