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

#include "MeshKernel/RangeCheck.hpp"
#include "MeshKernel/Utilities/UnderlyingType.hpp"

#include <string>
#include <vector>

namespace meshkernel
{
    /// @brief Projection::Type class
    class Projection
    {
    public:
        /// @brief Enumerator describing the projections
        enum class Type
        {
            Cartesian = 0,         ///> Cartesian
            Spherical = 1,         ///> Spherical
            SphericalAccurate = 2, ///> Spherical accurate
            Unknown = 3            ///> Unknown
        };

        /// @brief Returns the projection as a string
        /// @param type Type of projection
        /// @return Type of the projection as a string
        inline static std::string ToString(Type const type)
        {
            switch (type)
            {
            case Type::Cartesian:
                return "Cartesian";
            case Type::Spherical:
                return "Spherical";
            case Type::SphericalAccurate:
                return "Spherical accurate";
            case Type::Unknown:
            default:
                return "Unknown";
            }
        }

        /// @brief Checks the validity of the projection
        /// @param type Type of projection
        inline static void CheckValidity(Type const type)
        {
            range_check::CheckOneOf(to_underlying(type),
                                    m_valid_projections,
                                    "Projection type");
        }

    private:
        /// @brief Vector of valid projections stored as integers
        inline static std::vector<int> const m_valid_projections = {
            to_underlying(Type::Cartesian),
            to_underlying(Type::Spherical),
            to_underlying(Type::SphericalAccurate) //
        };
    };

} // namespace meshkernel
