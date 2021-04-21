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
        /// @brief Aspect ratio (mfacmax, 0.1)
        double aspect_ratio;

        /// @brief Grow factor of aspect ratio (1.1)
        double aspect_ratio_grow_factor;

        /// @brief Average mesh width on center spline (0.005)
        double average_width;

        /// @brief Curvature adapted grid spacing, 1 or not 0 (1)
        int curvature_adapted_grid_spacing;

        /// @brief Grow the grid outside the prescribed grid height (1)
        int grow_grid_outside = 0;

        /// @brief Maximum number of layers in the uniform part (5)
        int maximum_num_faces_in_uniform_part;

        /// @brief On-top-of-each-other tolerance (0.0001)
        double nodes_on_top_of_each_other_tolerance;

        /// @brief Minimum allowed absolute value of crossing-angle cosine (0.95)
        double min_cosine_crossing_angles;

        /// @brief Check for collisions with other parts of the front, 1 or not 0 (0)
        int check_front_collisions;

        /// @brief Remove skinny triangles (1)
        int remove_skinny_triangles;
    };
} // namespace meshkernelapi
