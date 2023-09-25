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

#include <string>
#include <vector>

namespace meshkernel
{
    /// @brief Projection::Type class
    class AveragingInterpolationMethod
    {
    public:
        /// @brief Class constructor
        AveragingInterpolationMethod() = delete;

        /// @enum Method
        /// @brief Enumerator describing the averaging interpolation methods
        enum Method
        {
            SimpleAveraging = 1,         ///< Computes a simple mean
            Closest = 2,                 ///< Takes the value of the closest sample to the interpolation location
            Max = 3,                     ///< Takes the maximum sample value
            Min = 4,                     ///< Takes the minimum sample value
            InverseWeightedDistance = 5, ///< Computes the inverse weighted sample mean
            MinAbsValue = 6,             ///< Computes the minimum absolute value
            Unknown = 7                  ///> Unknown
        };

        /// @brief Returns the mesh location string
        /// @param type The mesh location enumeration
        /// @return The mesh lcoations string
        inline static std::string ToString(Method const type)
        {
            switch (type)
            {
            case SimpleAveraging:
                return "Simple averaging (mean)";
            case Closest:
                return "Closest";
            case Max:
                return "Max";
            case Min:
                return "Min";
            case InverseWeightedDistance:
                return "Inverse weighted distance";
            case MinAbsValue:
                return "Minimum absolute value";
            case Unknown:
            default:
                return "Unknown";
            }
        }

        /// @brief Vector of valid projections stored as integers
        inline static std::vector<int> const ValidValues = {
            SimpleAveraging,
            Closest,
            Max,
            Min,
            InverseWeightedDistance,
            MinAbsValue //
        };
    };

} // namespace meshkernel
