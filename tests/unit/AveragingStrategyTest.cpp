#include "MeshKernel/AveragingInterpolation.hpp"

#include <gtest/gtest.h>

#include <memory>

#include "MeshKernel/ClosestAveragingStrategy.hpp"
#include "MeshKernel/InverseWeightedAveragingStrategy.hpp"
#include "MeshKernel/MaxAveragingStrategy.hpp"
#include "MeshKernel/MinAbsAveragingStrategy.hpp"
#include "MeshKernel/MinAveragingStrategy.hpp"
#include "MeshKernel/SimpleAveragingStrategy.hpp"

namespace meshkernel::averaging
{

    std::unique_ptr<AveragingStrategy> GetAveragingStrategy(
        AveragingInterpolation::Method const averagingMethod,
        double const missingValue,
        Point const& interpolationPoint,
        Projection const projection)
    {
        switch (averagingMethod)
        {
        case AveragingInterpolation::Method::SimpleAveraging:
            return std::make_unique<SimpleAveragingStrategy>(missingValue);
        case AveragingInterpolation::Method::Closest:
            return std::make_unique<ClosestAveragingStrategy>(missingValue, interpolationPoint, projection);
        case AveragingInterpolation::Method::Max:
            return std::make_unique<MaxAveragingStrategy>(missingValue);
        case AveragingInterpolation::Method::Min:
            return std::make_unique<MinAveragingStrategy>(missingValue);
        case AveragingInterpolation::Method::InverseWeightedDistance:
            return std::make_unique<InverseWeightedAveragingStrategy>(missingValue, interpolationPoint, projection);
        case AveragingInterpolation::Method::MinAbsValue:
            return std::make_unique<MinAbsAveragingStrategy>(missingValue);
        default:
            throw std::invalid_argument("Unsupported averagingMethod");
        }
    }

    class MissingValueIsReturnedWhenNoValueIsAddedTest : public ::testing::TestWithParam<std::pair<AveragingInterpolation::Method, double>>
    {
    public:
        [[nodiscard]] static std::vector<std::pair<AveragingInterpolation::Method, double>> get_data()
        {
            return {
                std::pair(AveragingInterpolation::Method::SimpleAveraging, doubleMissingValue),
                std::pair(AveragingInterpolation::Method::Closest, doubleMissingValue),
                std::pair(AveragingInterpolation::Method::Max, doubleMissingValue),
                std::pair(AveragingInterpolation::Method::Min, doubleMissingValue),
                std::pair(AveragingInterpolation::Method::InverseWeightedDistance, doubleMissingValue),
                std::pair(AveragingInterpolation::Method::MinAbsValue, doubleMissingValue),
            };
        }
    };

    TEST_P(MissingValueIsReturnedWhenNoValueIsAddedTest, expected_results)
    {
        // Setup
        Point p = Point(0.0, 0.0);
        std::unique_ptr<AveragingStrategy> pStrategy = GetAveragingStrategy(GetParam().first,
                                                                            GetParam().second,
                                                                            p, Projection::cartesian);

        // Call
        double result = pStrategy->Calculate();

        // Assert
        ASSERT_EQ(result, GetParam().second);
    }

    INSTANTIATE_TEST_SUITE_P(AveragingStrategyTest,
                             MissingValueIsReturnedWhenNoValueIsAddedTest,
                             ::testing::ValuesIn(MissingValueIsReturnedWhenNoValueIsAddedTest::get_data()));

} // namespace meshkernel::averaging
