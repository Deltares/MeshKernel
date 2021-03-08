#pragma once
#include "AveragingStrategy.hpp"

namespace meshkernel::averaging
{    
    /// @brief ClosestAveragingStrategy implements the AveragingStrategy which selects the value of the closest added point.
    class ClosestAveragingStrategy final : public AveragingStrategy
    {
    public:
        /// @brief Construct a new ClosestAveragingStrategy.
        /// @param[in] missingValue the value used to indicate a missing value.
        /// @param[in] interpolationPoint the point for which the average should be calculated.
        /// @param[in] the projection.
        ClosestAveragingStrategy(double missingValue,
                                 Point const& interpolationPoint,
                                 Projection projection);

        void Add(Point const& samplePoint, double sampleValue) override;
        [[nodiscard]] double Calculate() const override;
    private:
        double result_;
        double closestSquaredValue_ = std::numeric_limits<double>::max();
        Point const& interpolationPoint_;
        Projection const proj_;
    };
}