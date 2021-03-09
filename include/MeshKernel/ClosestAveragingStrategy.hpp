#pragma once
#include "AveragingStrategy.hpp"

namespace meshkernel::averaging
{
    /// @brief ClosestAveragingStrategy implements the AveragingStrategy which selects the value of the closest added point.
    class ClosestAveragingStrategy final : public AveragingStrategy
    {
    public:
        /// @brief Construct a new ClosestAveragingStrategy.
        /// @param[in] missingValue       The value used to indicate a missing value.
        /// @param[in] interpolationPoint The point for which the average should be calculated.
        /// @param[in] projection         The projection used to calculate distances with.
        ClosestAveragingStrategy(double missingValue,
                                 Point const& interpolationPoint,
                                 Projection projection);

        void Add(Point const& samplePoint, double sampleValue) override;
        [[nodiscard]] double Calculate() const override;

    private:
        /// @brief The result used to calculate the final value in Calculate.
        double m_result;

        /// @brief The closest squared value currently found.
        double m_closestSquaredValue = std::numeric_limits<double>::max();

        /// @brief The interpolation point from which the closest value is calculated.
        Point const& m_interpolationPoint;

        /// @brief The projection used to calculate the squared distance.
        Projection const m_projection;
    };
} // namespace meshkernel::averaging