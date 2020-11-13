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
    struct SampleRefineParametersNative
    {
        /// Sample vector dimension
        int SampleVectorDimension;

        /// Maximum number of refinement iterations, set to 1 if only one refinement is wanted
        int MaxNumberOfRefinementIterations;

        /// Minimum cell size
        double MinimumCellSize;

        /// Directional refinement, 1 yes 0 no
        int DirectionalRefinement;

        /// Refinement criterion type
        int RefinementType;

        /// Connect hanging nodes at the end of the iteration, 1 yes 0 no
        int ConnectHangingNodes;

        /// Maximum time-step in courant grid
        double MaximumTimeStepInCourantGrid;

        /// Take samples outside face into account , 1 yes 0 no
        int AccountForSamplesOutside;
    };
} // namespace meshkernelapi
