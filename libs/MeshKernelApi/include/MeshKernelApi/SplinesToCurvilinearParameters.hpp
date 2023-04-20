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
    /// @brief A struct used to describe the spline to curvilinear grid parameters in a C-compatible manner
    struct SplinesToCurvilinearParameters
    {
        /// @brief Aspect ratio (mfacmax)
        double aspect_ratio = 0.1;

        /// @brief Grow factor of aspect ratio
        double aspect_ratio_grow_factor = 1.1;

        /// @brief Average mesh width on center spline
        double average_width = 0.005;

        /// @brief Curvature adapted grid spacing, 1 or not 0
        int curvature_adapted_grid_spacing = 1;

        /// @brief Grow the grid outside the prescribed grid height
        int grow_grid_outside = 0;

        /// @brief Maximum number of layers in the uniform part
        int maximum_num_faces_in_uniform_part = 5;

        /// @brief On-top-of-each-other tolerance
        double nodes_on_top_of_each_other_tolerance = 0.0001;

        /// @brief Minimum allowed absolute value of crossing-angle cosine
        double min_cosine_crossing_angles = 0.95;

        /// @brief Check for collisions with other parts of the front, 1 or not 0
        int check_front_collisions = 0;

        /// @brief Check for collisions with other parts of the front
        int remove_skinny_triangles = 1;
    };
} // namespace meshkernelapi
