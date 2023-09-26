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

#include <MeshKernel/AveragingInterpolation.hpp>
#include <MeshKernel/AveragingInterpolationMethod.hpp>
#include <MeshKernel/AveragingStrategies/AveragingStrategyFactory.hpp>
#include <MeshKernel/AveragingStrategies/ClosestAveragingStrategy.hpp>
#include <MeshKernel/AveragingStrategies/InverseWeightedAveragingStrategy.hpp>
#include <MeshKernel/AveragingStrategies/MaxAveragingStrategy.hpp>
#include <MeshKernel/AveragingStrategies/MinAbsAveragingStrategy.hpp>
#include <MeshKernel/AveragingStrategies/MinAveragingStrategy.hpp>
#include <MeshKernel/AveragingStrategies/SimpleAveragingStrategy.hpp>
#include <MeshKernel/Exceptions.hpp>

namespace meshkernel
{

    std::unique_ptr<averaging::AveragingStrategy> averaging::AveragingStrategyFactory::GetAveragingStrategy(AveragingInterpolationMethod::Method averagingMethod,
                                                                                                            size_t minNumSamples,
                                                                                                            Point const& interpolationPoint,
                                                                                                            Projection::Type projection)
    {
        switch (averagingMethod)
        {
        case AveragingInterpolationMethod::Method::SimpleAveraging:
            return std::make_unique<SimpleAveragingStrategy>(minNumSamples);
        case AveragingInterpolationMethod::Method::Closest:
            return std::make_unique<ClosestAveragingStrategy>(interpolationPoint, projection);
        case AveragingInterpolationMethod::Method::Max:
            return std::make_unique<MaxAveragingStrategy>();
        case AveragingInterpolationMethod::Method::Min:
            return std::make_unique<MinAveragingStrategy>();
        case AveragingInterpolationMethod::Method::InverseWeightedDistance:
            return std::make_unique<InverseWeightedAveragingStrategy>(interpolationPoint, minNumSamples, projection);
        case AveragingInterpolationMethod::Method::MinAbsValue:
            return std::make_unique<MinAbsAveragingStrategy>();
        default:
            throw MeshKernelError("Unsupported averaging method. Valid methods are {}",
                                  AveragingInterpolationMethod::ValidValuesString());
        }
    }
} // namespace meshkernel
