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
    /// @brief A struct used to describe the interpolation parameters in a C-compatible manner
    struct InterpolationParameters
    {
        /// Maximum number of refinement iterations, set to 1 if only one refinement is wanted (10)
        int max_num_refinement_iterations;

        /// Averaging method : 1 = simple averaging, 2 = closest point, 3 = max, 4 = min, 5 = inverse weighted distance, 6 = minabs, 7 = kdtree (1)
        int averaging_method;

        /// Minimum number of points needed inside cell to handle the cell (1)
        int minimum_num_points;

        /// Relative search cell size, default 1= actual cell size, 2= twice as large, search radius can be larger than cell so more sample are included. (1.01)
        double relative_search_radius;

        /// Interpolation settings, 1=bathy, 2=zk, 3=s1, 4=Zc (2)
        int interpolate_to;

        /// Whether to compute faces intersected by polygon (yes=1/no=0)
        int refine_intersected;

        /// Whether to use the mass center when splitting a face in the refinement process (yes=1/no=0)
        int use_mass_center_when_refining;
    };
} // namespace meshkernelapi
