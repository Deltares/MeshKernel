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
    /// @brief A struct used to describe the sample refine parameters in a C-compatible manner
    struct SampleRefineParameters
    {
        /// @brief Sample vector dimension
        int SampleVectorDimension;

        /// @brief Maximum number of refinement iterations, set to 1 if only one refinement is wanted
        int MaxNumberOfRefinementIterations;

        /// @brief Minimum cell size
        double MinimumCellSize;

        /// @brief Directional refinement, 1 yes 0 no
        int DirectionalRefinement;

        /// @brief Refinement criterion type
        int RefinementType;

        /// @brief Connect hanging nodes at the end of the iteration, yes or no
        int ConnectHangingNodes;

        /// @brief Maximum time-step in courant grid
        double MaximumTimeStepInCourantGrid;

        /// @brief Take samples outside face into account , 1 yes 0 no
        int AccountForSamplesOutside;
    };
} // namespace meshkernelapi
