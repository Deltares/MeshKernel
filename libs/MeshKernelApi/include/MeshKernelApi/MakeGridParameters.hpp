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
    /// @brief This struct describes the necessary parameters to create a new curvilinear grid in a C-compatible manner
    ///
    /// @see mkernel_make_mesh
    struct MakeGridParameters
    {
        /// @brief The number of columns in x direction
        int num_columns = 3;

        /// @brief The number of columns in y direction
        int num_rows = 3;

        /// @brief The grid angle
        double angle = 0.0;

        /// @brief The grid block size, used in x and y direction
        double block_size = 10.0;

        /// @brief The x coordinate of the origin, located at the bottom left corner
        double origin_x = 0.0;

        /// @brief The y coordinate of the origin, located at the bottom left corner
        double origin_y = 0.0;

        /// @brief The grid block size in x dimension, used only for squared grids
        double block_size_x = 10.0;

        /// @brief The grid block size in y dimension, used only for squared grids
        double block_size_y = 10.0;
    };
} // namespace meshkernelapi
