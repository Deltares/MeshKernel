#pragma once
#include "AveragingStrategy.hpp"

namespace meshkernel::averaging
{
    /// @brief SimpleAveragingStrategy implements the averaging strategy which takes the average of the values of all added points.
    class SimpleAveragingStrategy final : public AveragingStrategy
    {
    public:
        /// @brief Construct a new SimpleAveragingStrategy implements the simple averaging strategy.
        explicit SimpleAveragingStrategy(double missingValue);

        void Add(Point const& samplePoint, double sampleValue) override;
        [[nodiscard]] double Calculate() const override;

    private:
        /// @brief The current result from which Calculate calculates the final value.
        double result_ = 0.0;

        /// @brief The number of times a value has been added to this strategy.
        size_t nAdds_ = 0;

        /// @brief The current result returned in Calculate
        double const missingValue_;
    };
} // namespace meshkernel::averaging