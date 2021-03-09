#pragma once
#include "AveragingStrategy.hpp"

namespace meshkernel::averaging
{
    /// @brief MinAbsAveragingStrategy implements the AveragingStrategy which takes the absolute minimum value of the added points.
    class MinAbsAveragingStrategy final : public AveragingStrategy
    {
    public:
        /// @brief Construct a new MinAbsAveragingStrategy.
        /// @param[in] missingValue The value used to indicate a missing value.
        explicit MinAbsAveragingStrategy(double missingValue);

        void Add(Point const& samplePoint, double sampleValue) override;

        [[nodiscard]] double Calculate() const override;

    private:
        /// @brief The current result returned in Calculate
        double m_result = std::numeric_limits<double>::max();

        /// @brief The value returned when no valid value can be returned.
        double const m_missingValue;
    };
} // namespace meshkernel::averaging
