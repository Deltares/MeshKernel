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
#include <MeshKernel/Entities.hpp>

namespace meshkernel::averaging
{
    /// @brief AveragingStrategy defines the averaging strategy to use in AveragingInterpolation
    class AveragingStrategy
    {
    public:
        virtual ~AveragingStrategy() = default;

        /// @brief Reset the state of the averaging strategy, ready for the next calculation.
        virtual void Reset(const Point& interpolationPoint) = 0;

        /// @brief Adds the specified sample point.
        /// @param[in] samplePoint The sample point to add to this strategy.
        /// @param[in] sampleValue The value associated with the sample point.
        virtual void Add(Point const& samplePoint, double sampleValue) = 0;

        /// @brief Calculates the average value based on the values added.
        /// @return The calculated average
        [[nodiscard]] virtual double Calculate() const = 0;

        virtual double Calculate (const Point& interpolationPoint,
                                  const std::vector<Point>& samplePoints,
                                  const std::vector<double>& sampleValues) {


            Reset (interpolationPoint);

            for (UInt i = 0; i < samplePoints.size (); ++i)
            {
                Add(samplePoints[i], sampleValues[i]);
            }

            return Calculate ();
        }
    };
} // namespace meshkernel::averaging
