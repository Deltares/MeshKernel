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
        double result_ = 0.0;
        double const missingValue_;
        double wall_ = 0.0;

        Point const& interpolationPoint_;
        Projection const proj_;
    };
}
