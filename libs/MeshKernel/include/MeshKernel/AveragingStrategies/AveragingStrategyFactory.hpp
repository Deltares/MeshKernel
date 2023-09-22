//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2021.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation version 3.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// contact: delft3d.support@deltares.nl
// Stichting Deltares
// P.O. Box 177
// 2600 MH Delft, The Netherlands
//
// All indications and logos of, and references to, "Delft3D" and "Deltares"
// are registered trademarks of Stichting Deltares, and remain the property of
// Stichting Deltares. All rights reserved.
//
//------------------------------------------------------------------------------

#pragma once

#include <MeshKernel/AveragingInterpolation.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Projection.hpp>

namespace meshkernel::averaging
{
    /// @brief A factory struct for getting the averaging strategies
    struct AveragingStrategyFactory
    {
        /// @brief The static method returning the strategy
        /// @param[in] averagingMethod The averaging method enumeration value
        /// @param[in] minNumSamples The minimum a of samples used for certain interpolation algorithms
        /// @param[in] interpolationPoint The interpolation point
        /// @param[in] projection  The projection to use
        /// @return The interpolation strategy to use
        [[nodiscard]] std::unique_ptr<AveragingStrategy> static GetAveragingStrategy(AveragingInterpolation::Method averagingMethod,
                                                                                     size_t minNumSamples,
                                                                                     Point const& interpolationPoint,
                                                                                     Projection::Type projection);
    };
} // namespace meshkernel::averaging
