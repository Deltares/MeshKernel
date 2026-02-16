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
    /// @brief A struct used to describe the intersection points of a spline with a number of other splines
    struct SplineIntersections
    {
        /// @brief The number of spline intersections
        int num_intersections = 0;

        /// @brief The index of the intersected spline
        int* spline_index = nullptr;

        /// @brief The angle of the intersection
        double* intersection_angle = nullptr;

        /// @brief The x coordinate of the intersection point
        double* intersection_x = nullptr;

        /// @brief The y coordinate of the intersection point
        double* intersection_y = nullptr;
    };

} // namespace meshkernelapi
