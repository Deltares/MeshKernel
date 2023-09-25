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
    class ProjectionToLandBoundary
    {
    public:
        /// @brief Class constructor
        ProjectionToLandBoundary() = delete;

        /// @enum Type
        /// @brief Enumerator describing the projections
        enum Type
        {
            None = 0,
            ToOriginalNetBoundary = 1,
            OuterMeshBoundaryToLandBoundary = 2,
            InnerAndOuterMeshBoundaryToLandBoundary = 3,
            WholeMesh = 4
        };

        /// @brief Returns the projection as a string
        /// @param type The projection enumeration
        /// @return The projection string
        inline static std::string ToString(Type const type)
        {
            switch (type)
            {
            case None:
                return "None";
            case ToOriginalNetBoundary:
                return "To original net boundary";
            case OuterMeshBoundaryToLandBoundary:
                return "Outer mesh boundary to land boundary";
            case InnerAndOuterMeshBoundaryToLandBoundary:
                return "Inner and outer mesh boundary to land boundary";
            case WholeMesh:
                return "Whole mesh";
            default:
                return "Unknown";
            }
        }

        /// @brief Vector of valid projections stored as integers
        inline static std::vector<int> const ValidValues = {
            None,
            ToOriginalNetBoundary,
            OuterMeshBoundaryToLandBoundary,
            InnerAndOuterMeshBoundaryToLandBoundary,
            WholeMesh //
        };
    };

} // namespace meshkernel
