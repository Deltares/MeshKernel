#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <vector>

#include "MeshKernel/SampleAveragingInterpolator.hpp"
#include "TestUtils/Definitions.hpp"
#include "TestUtils/MakeMeshes.hpp"
#include "TestUtils/SampleFileReader.hpp"

namespace mk = meshkernel;

TEST(SampleInterpolationTests, AveragingInterpolation)
{
    const mk::UInt numberOfPointsX = 11;
    const mk::UInt numberOfPointsY = 11;
    const mk::UInt numberOfPoints = numberOfPointsX * numberOfPointsY;

    std::vector<double> xPoints(numberOfPoints);
    std::vector<double> yPoints(numberOfPoints);
    std::vector<double> data(numberOfPoints);

    // Generate sample data points
    double delta = 1000.0;
    double x = 0.0;
    double y = 0.0;
    size_t count = 0;

    for (size_t i = 0; i < numberOfPointsY; ++i)
    {
        x = 0.0;

        for (size_t j = 0; j < numberOfPointsX; ++j)
        {
            xPoints[count] = x;
            yPoints[count] = y;
            data[count] = x;
            x += delta;
            ++count;
        }

        y += delta;
    }

    mk::InterpolationParameters params{.m_absoluteSearchRadius = 3.0 * delta};
    mk::SampleAveragingInterpolator interpolator(xPoints, yPoints, mk::Projection::cartesian, params);

    int propertyId = 1;
    interpolator.SetData(propertyId, data);

    // Execute
    ASSERT_EQ(interpolator.Size(), numberOfPoints);

    std::vector<mk::Point> interpolationPoints{{0.5 * delta, 0.5 * delta}, {9.5 * delta, 9.5 * delta}};
    const double initialValue = -1.0e20;
    std::vector<double> interpolationResult(interpolationPoints.size(), initialValue);
    std::vector<double> expectedResult{1400.0, 8600.0};

    interpolator.Interpolate(propertyId, interpolationPoints, interpolationResult);

    const double tolerance = 1.0e-8;

    for (size_t i = 0; i < expectedResult.size(); ++i)
    {
        EXPECT_NEAR(expectedResult[i], interpolationResult[i], tolerance);
    }
}
