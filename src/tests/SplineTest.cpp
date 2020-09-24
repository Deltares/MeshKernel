#include "../Mesh.hpp"
#include "../Entities.hpp"
#include "../Polygons.hpp"
#include "../Constants.cpp"
#include "../CurvilinearGrid.hpp"
#include "../Splines.hpp"
#include <gtest/gtest.h>

TEST(Splines, SetSpline)
{
    //One gets the edges
    std::vector<MeshKernel::Point> splineNodes;
    splineNodes.push_back(MeshKernel::Point{ MeshKernel::doubleMissingValue, MeshKernel::doubleMissingValue });
    splineNodes.push_back(MeshKernel::Point{ MeshKernel::doubleMissingValue, MeshKernel::doubleMissingValue });
    splineNodes.push_back(MeshKernel::Point{ MeshKernel::doubleMissingValue, MeshKernel::doubleMissingValue });
    splineNodes.push_back(MeshKernel::Point{ MeshKernel::doubleMissingValue, MeshKernel::doubleMissingValue });

    MeshKernel::Splines splines(MeshKernel::Projections::cartesian);
    bool successful = splines.AddSpline(splineNodes, 0, int(splineNodes.size()));
    ASSERT_TRUE(successful);

    ASSERT_EQ(1, splines.m_numSplines);
    ASSERT_EQ(10, splines.m_splineNodes[0].size());
    ASSERT_EQ(5, splines.m_numAllocatedSplines);
    ASSERT_EQ(10, splines.m_numAllocatedSplineNodes[0]);
}

TEST(Splines, CubicSplineInterpolation)
{
    //One gets the edges
    std::vector<MeshKernel::Point> splineNodes;

    splineNodes.push_back(MeshKernel::Point{ 212.001953125000, 155.627197265625 });
    splineNodes.push_back(MeshKernel::Point{ 529.253906250000, 432.379974365234 });
    splineNodes.push_back(MeshKernel::Point{ 930.506469726562, 453.380187988281 });

    int pointsBetweenVertices = 20;
    std::vector<MeshKernel::Point> coordinatesDerivatives(splineNodes.size());
    MeshKernel::Splines::SecondOrderDerivative(splineNodes, splineNodes.size(), coordinatesDerivatives);
    std::vector<MeshKernel::Point> splineCoordinates;

    for (int n = 0; n < splineNodes.size() - 1; n++)
    {
        for (int p = 0; p <= pointsBetweenVertices; p++)
        {
            const double pointAdimensionalCoordinate = n + double(p) / double(pointsBetweenVertices);
            MeshKernel::Point pointCoordinate;
            InterpolateSplinePoint(splineNodes, coordinatesDerivatives, pointAdimensionalCoordinate, pointCoordinate);
            splineCoordinates.push_back({ pointCoordinate.x, pointCoordinate.y });
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
    std::vector<MeshKernel::Point> firstSpline;
    firstSpline.push_back(MeshKernel::Point{ 152.001571655273, 86.6264953613281 });
    firstSpline.push_back(MeshKernel::Point{ 374.752960205078, 336.378997802734 });
    firstSpline.push_back(MeshKernel::Point{ 850.255920410156, 499.130676269531 });

    MeshKernel::Splines splines(MeshKernel::Projections::cartesian);

    bool successful = splines.AddSpline(firstSpline, 0, firstSpline.size());
    ASSERT_TRUE(successful);

    std::vector<MeshKernel::Point> secondSpline;
    secondSpline.push_back(MeshKernel::Point{ 72.5010681152344,391.129577636719 });
    secondSpline.push_back(MeshKernel::Point{ 462.503479003906, 90.3765411376953 });
    successful = splines.AddSpline(secondSpline, 0, secondSpline.size());

    double crossProductIntersection;
    MeshKernel::Point dimensionalIntersection;
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