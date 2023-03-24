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
    /// @brief A struct used to describe the orthogonalization parameters in a C-compatible manner
    struct OrthogonalizationParameters
    {
        /// @brief Number of outer iterations in orthogonalization. Increase this parameter for complex grids (2)
        int outer_iterations;

        /// @brief Number of boundary iterations in grid/net orthogonalization within itatp (25)
        int boundary_iterations;

        /// @brief Number of inner iterations in grid/net orthogonalization within itbnd (25)
        int inner_iterations;

        /// @brief Factor from 0 to 1. between grid smoothing and grid orthogonality (0.975)
        double orthogonalization_to_smoothing_factor;

        /// @brief Minimum ATPF on the boundary (1.0)
        double orthogonalization_to_smoothing_factor_at_boundary;

        /// @brief Factor between smoother 1d0 and area-homogenizer 0d0 (1.0)
        double areal_to_angle_smoothing_factor;
    };
} // namespace meshkernelapi
