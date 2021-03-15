#include <memory>
#include <utility>
#include <gtest/gtest.h>

#include "MeshKernel/AveragingInterpolation.hpp"
#include "MeshKernel/AveragingStrategies/ClosestAveragingStrategy.hpp"
#include "MeshKernel/AveragingStrategies/InverseWeightedAveragingStrategy.hpp"
#include "MeshKernel/AveragingStrategies/MaxAveragingStrategy.hpp"
#include "MeshKernel/AveragingStrategies/MinAbsAveragingStrategy.hpp"
#include "MeshKernel/AveragingStrategies/MinAveragingStrategy.hpp"
#include "MeshKernel/AveragingStrategies/SimpleAveragingStrategy.hpp"

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
            return std::make_unique<SimpleAveragingStrategy>();
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
        [[nodiscard]] static std::vector<std::pair<AveragingInterpolation::Method, double>> GetData()
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
        ASSERT_EQ(GetParam().second, result);
    }

    INSTANTIATE_TEST_SUITE_P(AveragingStrategyTest,
                             MissingValueIsReturnedWhenNoValueIsAddedTest,
                             ::testing::ValuesIn(MissingValueIsReturnedWhenNoValueIsAddedTest::GetData()));

    class CalculateWithAddedValuesData
    {
    public:
        CalculateWithAddedValuesData(AveragingInterpolation::Method const method, std::vector<std::pair<Point, double>> const& addData, double const expectedResult) : method_(method), addData_(addData), expectedResult(expectedResult) {}

        AveragingInterpolation::Method const method_;
        std::vector<std::pair<Point, double>> const addData_;
        double const expectedResult;
    };

    class CalculateWithAddedValuesTest : public ::testing::TestWithParam<CalculateWithAddedValuesData>
    {
    public:
        [[nodiscard]] static std::vector<CalculateWithAddedValuesData> GetData()
        {
            Point p = Point(0, 0);
            auto simpleData = {
                std::make_pair(p, -1.0),
                std::make_pair(p, -2.0),
                std::make_pair(p, -3.0),
                std::make_pair(p, -4.0),
            };

            auto closestData = {
                std::make_pair(Point(4.0, 0.0), 4.0),
                std::make_pair(Point(-3.0, 0.0), -3.0),
                std::make_pair(Point(2.0, 0.0), 2.0),
                std::make_pair(Point(-1.0, 0.0), -1.0),
            };

            auto inverseWeigthedData = {
                std::make_pair(Point(4.0, 2.0), 5.0),
                std::make_pair(Point(-3.0, 1.0), -3.0),
                std::make_pair(Point(2.0, -3.0), 2.0),
                std::make_pair(Point(-1.0, -2.0), -1.0),
            };

            return {
                CalculateWithAddedValuesData(AveragingInterpolation::Method::SimpleAveraging, simpleData, -2.5),
                CalculateWithAddedValuesData(AveragingInterpolation::Method::Max, simpleData, -1.0),
                CalculateWithAddedValuesData(AveragingInterpolation::Method::Min, simpleData, -4.0),
                CalculateWithAddedValuesData(AveragingInterpolation::Method::MinAbsValue, simpleData, 1.0),
                CalculateWithAddedValuesData(AveragingInterpolation::Method::Closest, closestData, -1.0),
                CalculateWithAddedValuesData(AveragingInterpolation::Method::InverseWeightedDistance, inverseWeigthedData, 0.218948),
            };
        }
    };

    TEST_P(CalculateWithAddedValuesTest, expected_ersults)
    {
        // Setup
        Point p = Point(0.0, 0.0);
        std::unique_ptr<AveragingStrategy> pStrategy = GetAveragingStrategy(GetParam().method_,
                                                                            doubleMissingValue,
                                                                            p, Projection::cartesian);

        for (auto const p : GetParam().addData_)
        {
            pStrategy->Add(p.first, p.second);
        }

        // Call
        double result = pStrategy->Calculate();

        // Assert
        constexpr double tolerance = 1e-6;
        ASSERT_NEAR(GetParam().expectedResult, result, tolerance);
    }

    INSTANTIATE_TEST_SUITE_P(AveragingStrategyTest,
                             CalculateWithAddedValuesTest,
                             ::testing::ValuesIn(CalculateWithAddedValuesTest::GetData()));

} // namespace meshkernel::averaging
