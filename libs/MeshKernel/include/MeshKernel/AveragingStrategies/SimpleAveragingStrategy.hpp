﻿#pragma once
#include <MeshKernel/AveragingStrategies/AveragingStrategy.hpp>
#include <MeshKernel/AveragingStrategies/SimpleAveragingStrategy.hpp>
#include <MeshKernel/Entities.hpp>

namespace meshkernel::averaging
{
    /// @brief SimpleAveragingStrategy implements the averaging strategy which takes the average of the values of all added points.
    class SimpleAveragingStrategy final : public AveragingStrategy
    {
    public:
        /// @brief Constructor taking the minimum amount of points parameter
        /// @param minNumSamples[in] The minimum amount of samples for a valid interpolation
        SimpleAveragingStrategy(size_t minNumSamples);

        /// @brief Reset the state of the simple averaging-strategy.
        void Reset(const Point& interpolationPoint) override;

        void Add(Point const& samplePoint, double sampleValue) override;
        [[nodiscard]] double Calculate() const override;

        /// @brief Calculates the average value based on the sample values.
        /// @param[in] interpolationPoint The point for which the average should be calculated.
        /// @param[in] samplePoints The sample points to used by this strategy.
        /// @param[in] sampleValues The sample values  associated with each sample point.
        /// @return The calculated average
        double Calculate (const Point& interpolationPoint,
                          const std::vector<Point>& samplePoints,
                          const std::vector<double>& sampleValues) const override;

    private:
        /// @brief The minimum number of points for a valid interpolation.
        const size_t m_minNumPoints;

        /// @brief The current result from which Calculate calculates the final value.
        double m_result = 0.0;

        /// @brief The number of times a value has been added to this strategy.
        size_t m_nAdds = 0;
    };
} // namespace meshkernel::averaging
