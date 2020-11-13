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
    struct SplinesToCurvilinearParametersNative
    {
        /// <summary>
        /// * Aspect ratio (mfacmax, 0.1)
        /// </summary>
        double AspectRatio;

        /// <summary>
        /// * Grow factor of aspect ratio (1.1)
        /// </summary>
        double AspectRatioGrowFactor;

        /// <summary>
        /// * Average mesh width on center spline (0.005)
        /// </summary>
        double AverageWidth;

        /// <summary>
        /// * Curvature adapted grid spacing, 1 or not 0 (1)
        /// </summary>
        int CurvatureAdapetedGridSpacing;

        /// <summary>
        /// * Grow the grid outside the prescribed grid height (1)
        /// </summary>
        int GrowGridOutside;

        /// <summary>
        /// * Maximum number of layers in the uniform part (5)
        /// </summary>
        int MaximumNumberOfGridCellsInTheUniformPart;

        /// <summary>
        /// * On-top-of-each-other tolerance (0.0001)
        /// </summary>
        double GridsOnTopOfEachOtherTolerance;

        /// <summary>
        /// * Minimum allowed absolute value of crossing-angle cosine (0.95)
        /// </summary>
        double MinimumCosineOfCrossingAngles;

        /// <summary>
        /// * Check for collisions with other parts of the front, 1 or not 0 (0)
        /// </summary>
        int CheckFrontCollisions;

        /// <summary>
        /// * Uniform grid size, netboundary to grid only (0.0)
        /// </summary>
        double UniformGridSize;

        /// <summary>
        /// * Remove skinny triangles (1)
        /// </summary>
        int RemoveSkinnyTriangles;
    };
} // namespace meshkernelapi
