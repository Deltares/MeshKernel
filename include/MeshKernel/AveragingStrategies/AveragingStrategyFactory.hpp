#pragma once

#include <MeshKernel/AveragingInterpolation.hpp>
#include <MeshKernel/Entities.hpp>

namespace meshkernel::averaging
{
    /// @brief A factory struct for getting the averaging strategies
    struct AveragingStrategyFactory
    {
        /// @brief The static method returning the strategy
        /// @param averagingMethod The averaging method enumeration value
        /// @param interpolationPoint The interpolation point
        /// @param projection  The projection to use
        /// @return The interpolation strategy to use
        [[nodiscard]] std::unique_ptr<AveragingStrategy> static GetAveragingStrategy(AveragingInterpolation::Method averagingMethod,
                                                                                     Point const& interpolationPoint,
                                                                                     Projection projection);
    };
} // namespace meshkernel::averaging