//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2020.
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
    struct InterpolationParametersNative
    {
        /// Actual interpolation type (1)
        int InterpolationType;

        /// Variable related to interactor behaviour (0)
        int DisplayInterpolationProcess;

        /// Maximum number of refinement iterations, set to 1 if only one refinement is wanted (10)
        int MaxNumberOfRefinementIterations;

        /// Averaging method : 1 = simple averaging, 2 = closest point, 3 = max, 4 = min, 5 = inverse weighted distance, 6 = minabs, 7 = kdtree (1)
        int AveragingMethod;

        /// Minimum number of points needed inside cell to handle the cell (1)
        int MinimumNumberOfPoints;

        /// Relative search cell size, default 1= actual cell size, 2= twice as large, search radius can be larger than cell so more sample are included. (1.01)
        double RelativeSearchRadius;

        /// Interpolation settings, 1=bathy, 2=zk, 3=s1, 4=Zc (2)
        int InterpolateTo;

        /// Refine faces intersected by polygon
        bool RefineIntersected;
    };
} // namespace meshkernelapi
