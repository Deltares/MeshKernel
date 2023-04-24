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
    /// @brief A struct describing gridded samples
    struct GriddedSamples
    {
        /// @brief Number of grid columns
        int n_cols = 0;

        /// @brief Number of grid rows
        int n_rows = 0;

        /// @brief X coordinate of the grid origin
        double x_origin = 0.0;

        /// @brief Y coordinate of the grid origin
        double y_origin = 0.0;

        /// @brief Type of the origin (centre=0 / corner=1)
        int origin_location_type = 0;

        /// @brief Constant grid cell size
        double cell_size = 0.0;

        /// @brief If not nullptr, coordinates for non-uniform grid spacing in x direction
        double* x_coordinates = nullptr;

        /// @brief If not nullptr, coordinates for non-uniform grid spacing in y direction
        double* y_coordinates = nullptr;

        /// @brief Value for missing data
        double missing_value = 1e10;

        /// @brief Sample values
        double* values = nullptr;
    };
} // namespace meshkernelapi
