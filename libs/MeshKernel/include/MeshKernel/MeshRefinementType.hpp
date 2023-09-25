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
    class MeshRefinementType
    {
    public:
        /// @brief Class constructor
        MeshRefinementType() = delete;

        /// @enum Type
        /// @brief Enumerator describing the mesh refinements types
        enum Type
        {
            WaveCourant = 1,      ///> Wave Couranr
            RefinementLevels = 2, ///> Refinement levels
            Unknown = 3           ///> Unknown
        };

        /// @brief Returns the mesh location string
        /// @param type The mesh location enumeration
        /// @return The mesh lcoations string
        inline static std::string ToString(Type const type)
        {
            switch (type)
            {
            case WaveCourant:
                return "Wave Courant";
            case RefinementLevels:
                return "Refinement levels";
            case Unknown:
            default:
                return "Unknown";
            }
        }

        /// @brief Vector of valid projections stored as integers
        inline static std::vector<int> const ValidValues = {WaveCourant, RefinementLevels};
    };

} // namespace meshkernel
