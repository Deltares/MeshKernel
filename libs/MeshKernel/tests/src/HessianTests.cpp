#include <gtest/gtest.h>

#include <algorithm>
#include <iomanip>
#include <random>
#include <vector>

#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Point.hpp"
#include "MeshKernel/RidgeRefinement.hpp"

namespace mk = meshkernel;

std::vector<std::vector<mk::Point>> GenerateGridPoints(const mk::UInt rows, const mk::UInt cols)
{

    std::vector<std::vector<mk::Point>> meshPoints(cols, std::vector<mk::Point>(rows));

    mk::Point origin = mk::Point(0.0, 0.0);

    double hx = 1.0;
    double hy = 1.0;

    mk::Point current = origin;

    for (mk::UInt ix = 0; ix < cols; ++ix)
    {
        current.y = origin.y;

        for (mk::UInt jy = 0; jy < rows; ++jy)
        {
            meshPoints[ix][jy] = current;
            current.y += hy;
        }

        current.x += hx;
    }

    return meshPoints;
}

TEST(HessianTests, TidyPointsNoDuplicatesTest)
{
    // Test that the TidySamples does not change a data set that has no duplicates.

    const mk::UInt size = 100;

    std::vector<mk::Point> samplePoints(size);
    std::vector<mk::Point> generated;
    std::vector<double> sampleData(size);

    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_real_distribution<double> distribution(0.0, 100.0);

    for (mk::UInt i = 0; i < size; ++i)
    {
        samplePoints[i] = mk::Point(distribution(generator), distribution(generator));
        sampleData[i] = distribution(generator);
    }

    std::vector<mk::Point> expectedSamplePoints = samplePoints;
    std::vector<double> expectedSampleData = sampleData;

    mk::RidgeRefinement ridgeRefinement;
    ridgeRefinement.TidySamples(samplePoints, sampleData);

    for (mk::UInt i = 0; i < size; ++i)
    {
        EXPECT_EQ(expectedSamplePoints[i].x, samplePoints[i].x);
        EXPECT_EQ(expectedSamplePoints[i].y, samplePoints[i].y);
        EXPECT_EQ(expectedSampleData[i], sampleData[i]);
    }
}

TEST(HessianTests, TidyPointsTest)
{

    const mk::UInt size = 100;

    std::vector<mk::Point> samplePoints(size);
    std::vector<mk::Point> expectedSamplePoints;
    std::vector<mk::Point> generated;
    std::vector<double> sampleData(size);

    const double invMax = 1.0 / static_cast<double>(RAND_MAX);

    for (mk::UInt i = 0; i < size; ++i)
    {
        // Use the simple pseudo random number generator to have repeatable values between tests when debugging.
        samplePoints[i] = mk::Point(rand() * invMax, rand() * invMax);
        sampleData[i] = rand() * invMax;

        if (i > 0 && i % 10 == 0)
        {
            samplePoints[i - 1] = samplePoints[i];
            sampleData[i - 1] = sampleData[i];
        }
        else if (i > 10 && i % 21 == 0)
        {
            samplePoints[i - 10] = samplePoints[i];
            sampleData[i - 10] = sampleData[i];
        }
        else
        {
            expectedSamplePoints.push_back(samplePoints[i]);
        }
    }

    mk::RidgeRefinement ridgeRefinement;
    ridgeRefinement.TidySamples(samplePoints, sampleData);

    EXPECT_EQ(expectedSamplePoints.size(), samplePoints.size());

#if 0
    for (const auto& point : expectedSamplePoints)
    {
        EXPECT_TRUE(std::find(samplePoints.begin(), samplePoints.end(), point) != samplePoints.end());
    }
#endif
}
