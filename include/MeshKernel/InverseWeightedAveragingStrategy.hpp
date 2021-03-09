#pragma once
#include "AveragingStrategy.hpp"

namespace meshkernel::averaging
{
    /// @brief ClosestAveragingStrategy implements the AveragingStrategy which weights the value of the added points based on their inverse distance.
    class InverseWeightedAveragingStrategy final : public AveragingStrategy
    {
    public:
        /// @brief Construct a new InverseWeightedAveragingStrategy.
        /// @param[in] missingValue       the value used to indicate a missing value.
        /// @param[in] interpolationPoint the point for which the average should be calculated.
        /// @param[in] projection         the projection used in calculating the distance.
        InverseWeightedAveragingStrategy(double missingValue,
                                         Point const& interpolationPoint,
                                         Projection projection);

        void Add(Point const& samplePoint, double sampleValue) override;
        [[nodiscard]] double Calculate() const override;

    private:
        /// @brief The current result used in Calculate to calculate the final value.
        double result_ = 0.0;

        /// @brief The value returned when no valid value can be returned.
        double const missingValue_;

        /// @brief The wall
        double wall_ = 0.0;

        /// @brief The interpolation point from which the inverse weight is calculated.
        Point const& interpolationPoint_;

        /// @brief The projection used to calculate the distance.
        Projection const proj_;
    };
} // namespace meshkernel::averaging
