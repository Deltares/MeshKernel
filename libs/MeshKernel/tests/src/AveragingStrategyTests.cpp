#include <gtest/gtest.h>
#include <memory>
#include <utility>

#include <MeshKernel/AveragingInterpolation.hpp>
#include <MeshKernel/AveragingStrategies/AveragingStrategyFactory.hpp>
#include <MeshKernel/AveragingStrategies/ClosestAveragingStrategy.hpp>

namespace meshkernel::averaging
{

    class MissingValueIsReturnedWhenNoValueIsAddedTest : public ::testing::TestWithParam<std::pair<AveragingInterpolationMethod::Method, double>>
    {
    public:
        [[nodiscard]] static std::vector<std::pair<AveragingInterpolationMethod::Method, double>> GetData()
        {
            return {
                {AveragingInterpolationMethod::Method::SimpleAveraging, constants::missing::doubleValue},
                {AveragingInterpolationMethod::Method::Closest, constants::missing::doubleValue},
                {AveragingInterpolationMethod::Method::Max, constants::missing::doubleValue},
                {AveragingInterpolationMethod::Method::Min, constants::missing::doubleValue},
                {AveragingInterpolationMethod::Method::InverseWeightedDistance, constants::missing::doubleValue},
                {AveragingInterpolationMethod::Method::MinAbsValue, constants::missing::doubleValue},
            };
        }
    };

    TEST_P(MissingValueIsReturnedWhenNoValueIsAddedTest, expected_results)
    {
        // Setup
        Point p = Point(0.0, 0.0);
        std::unique_ptr<AveragingStrategy> pStrategy = AveragingStrategyFactory::GetAveragingStrategy(GetParam().first,
                                                                                                      1,
                                                                                                      p, Projection::Cartesian);

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
        CalculateWithAddedValuesData(AveragingInterpolationMethod::Method const method, std::vector<std::pair<Point, double>> const& addData, double const expectedResult) : method_(method), addData_(addData), expectedResult(expectedResult) {}

        AveragingInterpolationMethod::Method const method_;
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
                CalculateWithAddedValuesData(AveragingInterpolationMethod::Method::SimpleAveraging, simpleData, -2.5),
                CalculateWithAddedValuesData(AveragingInterpolationMethod::Method::Max, simpleData, -1.0),
                CalculateWithAddedValuesData(AveragingInterpolationMethod::Method::Min, simpleData, -4.0),
                CalculateWithAddedValuesData(AveragingInterpolationMethod::Method::MinAbsValue, simpleData, 1.0),
                CalculateWithAddedValuesData(AveragingInterpolationMethod::Method::Closest, closestData, -1.0),
                CalculateWithAddedValuesData(AveragingInterpolationMethod::Method::InverseWeightedDistance, inverseWeigthedData, 0.218948),
            };
        }
    };

    TEST_P(CalculateWithAddedValuesTest, expected_ersults)
    {
        // Setup
        Point p = Point(0.0, 0.0);
        std::unique_ptr<AveragingStrategy> pStrategy = AveragingStrategyFactory::GetAveragingStrategy(GetParam().method_, 1,
                                                                                                      p, Projection::Cartesian);

        for (auto const& p : GetParam().addData_)
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
