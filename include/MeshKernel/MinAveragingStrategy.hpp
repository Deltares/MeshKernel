#pragma once
#include "AveragingStrategy.hpp"

namespace meshkernel::averaging
{
    /// @brief MinAbsAveragingStrategy implements the AveragingStrategy which takes the minimum value of the added points.
    class MinAveragingStrategy final : public AveragingStrategy
    {
    public:
        /// @brief Construct a new MinAveragingStrategy.
        /// @param[in] missingValue The value used to indicate a missing value.
        explicit MinAveragingStrategy(double missingValue);

        void Add(Point const& samplePoint, double sampleValue) override;

        [[nodiscard]] double Calculate() const override;

    private:
        double result_ = std::numeric_limits<double>::max();
        double const missingValue_;
    };
}
