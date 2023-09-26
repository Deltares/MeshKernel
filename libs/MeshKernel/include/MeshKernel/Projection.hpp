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

#include "MeshKernel/Definitions.hpp"

namespace meshkernel
{
    /// @brief Projection::Type class
    class Projection
    {
    public:
        /// @brief Class constructor
        Projection() = delete;

        /// @enum Type
        /// @brief Enumerator describing the projections
        enum Type
        {
            Cartesian = 0,         ///> Cartesian
            Spherical = 1,         ///> Spherical
            SphericalAccurate = 2, ///> Spherical accurate
            Unknown = 3            ///> Unknown
        };

        /// @brief Returns the projection as a string
        /// @param type The projection enumeration
        /// @return The projection string
        inline static std::string ToString(Type const type)
        {
            switch (type)
            {
            case Cartesian:
                return "Cartesian";
            case Spherical:
                return "Spherical";
            case SphericalAccurate:
                return "Spherical accurate";
            case Unknown:
            default:
                return "Unknown";
            }
        }

        /// @brief Returns the mesh location string
        /// @param type The mesh location enumeration as an integer
        /// @return The mesh lcoations string
        inline static std::string ToString(int const type)
        {
            return ToString(static_cast<Type>(type));
        }

        /// @brief Gets vector of valid values (integers)
        /// @return Vector of valid values
        inline static std::vector<int> const& ValidValues() { return m_ValidValues; }

        /// @brief Gets valid values as a string
        /// @return String of valid values
        inline static std::string const& ValidValuesString() { return m_ValidValuesString; }

    private:
        /// @brief Vector of valid values
        inline static std::vector<int> const m_ValidValues = {Cartesian, Spherical, SphericalAccurate};

        /// @brief String of valid values
        inline static std::string const m_ValidValuesString = MakeValidValuesString<Projection>();
    };

} // namespace meshkernel
