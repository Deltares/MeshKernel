#include <gtest/gtest.h>

#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Splines.hpp>

TEST(Splines, SetSpline)
{
    //a valid spline
    std::vector<meshkernel::Point> splineNodes({{212.001953125000, 155.627197265625},
                                                {529.253906250000, 432.379974365234},
                                                {930.506469726562, 453.380187988281},
                                                {1030.506469726562, 653.380187988281}});

    meshkernel::Splines splines(meshkernel::Projection::cartesian);
    splines.AddSpline(splineNodes, 0, splineNodes.size());

    ASSERT_EQ(1, splines.GetNumSplines());
    ASSERT_EQ(4, splines.m_splineNodes[0].size());
}

TEST(Splines, CubicSplineInterpolation)
{
    //One gets the edges
    std::vector<meshkernel::Point> splineNodes;
    splineNodes.push_back(meshkernel::Point{212.001953125000, 155.627197265625});
    splineNodes.push_back(meshkernel::Point{529.253906250000, 432.379974365234});
    splineNodes.push_back(meshkernel::Point{930.506469726562, 453.380187988281});

    int pointsBetweenNodes = 20;
    auto coordinatesDerivatives = meshkernel::Splines::SecondOrderDerivative(splineNodes, 0, splineNodes.size() - 1);
    std::vector<meshkernel::Point> splineCoordinates;

    for (size_t n = 0; n < splineNodes.size() - 1; n++)
    {
        for (auto p = 0; p <= pointsBetweenNodes; p++)
        {
            const double pointAdimensionalCoordinate = n + double(p) / double(pointsBetweenNodes);
            auto pointCoordinate = ComputePointOnSplineAtAdimensionalDistance(splineNodes, coordinatesDerivatives, pointAdimensionalCoordinate);
            ASSERT_TRUE(pointCoordinate.IsValid());

            splineCoordinates.push_back({pointCoordinate.x, pointCoordinate.y});
        }
    }

    const double tolerance = 1e-6;
    ASSERT_NEAR(226.817168170929, splineCoordinates[1].x, tolerance);
    ASSERT_NEAR(241.648133331299, splineCoordinates[2].x, tolerance);
    ASSERT_NEAR(256.510598720551, splineCoordinates[3].x, tolerance);
    ASSERT_NEAR(271.420314453125, splineCoordinates[4].x, tolerance);
    ASSERT_NEAR(286.393030643463, splineCoordinates[5].x, tolerance);
    ASSERT_NEAR(930.506469726562, splineCoordinates.back().x, tolerance);

    ASSERT_NEAR(172.653750896454, splineCoordinates[1].y, tolerance);
    ASSERT_NEAR(189.632350921631, splineCoordinates[2].y, tolerance);
    ASSERT_NEAR(206.515043735504, splineCoordinates[3].y, tolerance);
    ASSERT_NEAR(223.253875732422, splineCoordinates[4].y, tolerance);
    ASSERT_NEAR(453.380187988281, splineCoordinates.back().y, tolerance);
}

TEST(Splines, SplineIntersection)
{
    std::vector<meshkernel::Point> firstSpline;
    firstSpline.push_back(meshkernel::Point{152.001571655273, 86.6264953613281});
    firstSpline.push_back(meshkernel::Point{374.752960205078, 336.378997802734});
    firstSpline.push_back(meshkernel::Point{850.255920410156, 499.130676269531});

    meshkernel::Splines splines(meshkernel::Projection::cartesian);

    splines.AddSpline(firstSpline, 0, firstSpline.size());

    std::vector<meshkernel::Point> secondSpline;
    secondSpline.push_back(meshkernel::Point{72.5010681152344, 391.129577636719});
    secondSpline.push_back(meshkernel::Point{462.503479003906, 90.3765411376953});
    splines.AddSpline(secondSpline, 0, secondSpline.size());

    double crossProductIntersection;
    meshkernel::Point dimensionalIntersection;
    double firstSplineRatio;
    double secondSplineRatio;

    splines.GetSplinesIntersection(0, 1, crossProductIntersection, dimensionalIntersection, firstSplineRatio, secondSplineRatio);

    const double tolerance = 1e-6;
    ASSERT_NEAR(261.736770097059, dimensionalIntersection.x, tolerance);
    ASSERT_NEAR(245.199166962145, dimensionalIntersection.y, tolerance);
    ASSERT_NEAR(0.601498208554790, firstSplineRatio, tolerance);
    ASSERT_NEAR(0.485216749175026, secondSplineRatio, tolerance);
    ASSERT_NEAR(-0.996215079635043, crossProductIntersection, tolerance);
}
