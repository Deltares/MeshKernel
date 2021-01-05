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
    /// @brief A struct used to describe parameters for generating a curvilinear grid in a C-compatible manner
    struct CurvilinearParameters
    {
        /// @brief M-refinement factor for regular grid generation (mfacmax, 2000)
        int MRefinement;

        /// @brief N-refinement factor for regular grid generation (nfacmax, 40)
        int NRefinement;

        /// @brief Nr. of inner iterations in regular grid smoothing (10).
        int SmoothingIterations;

        /// @brief Smoothing parameter (0.5).
        double SmoothingParameter;

        /// @brief Attraction/repulsion parameter (0).
        double AttractionParameter;
    };
} // namespace meshkernelapi
