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

namespace meshkernelapi
{
    /// @brief A struct used to describe a list of geometries in a C-compatible manner
    struct GeometryList
    {
        /// @brief The value used as separator in coordinates_x, coordinates_y and values
        double geometry_separator = -999.0;

        /// @brief The value used to separate the inner part of a polygon from its outer part
        double inner_outer_separator = -998.0;

        /// @brief The number of coordinate values present
        int num_coordinates = 0;

        /// @brief The x coordinate values
        double* coordinates_x = nullptr;

        /// @brief The y coordinate values
        double* coordinates_y = nullptr;

        /// @brief The values at the point
        double* values = nullptr;
    };
} // namespace meshkernelapi
