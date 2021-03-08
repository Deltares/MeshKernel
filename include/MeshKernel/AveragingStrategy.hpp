#pragma once
#include "Entities.hpp"

namespace meshkernel::averaging
{    
    /// @brief AveragingStrategy defines the averaging strategy to use in AveragingInterpolation
    class AveragingStrategy
    {
    public:
        virtual ~AveragingStrategy() = default;
        
        /// @brief Adds the specified sample point.
        /// @param[in] samplePoint The sample point to add to this strategy.
        /// @param[in] sampleValue The value associated with the sample point.
        virtual void Add(Point const& samplePoint, double sampleValue) = 0;
        
        /// @brief Calculates the average value based on the values added.
        /// @return The calculated average
        [[nodiscard]] virtual double Calculate() const = 0;
    };
}
