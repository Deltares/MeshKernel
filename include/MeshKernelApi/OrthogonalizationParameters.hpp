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
    /// @brief A struct used to describe the orthogonalization parameters in a C-compatible manner
    struct OrthogonalizationParameters
    {
        /// @brief Number of outer iterations in orthogonalization. Increase this parameter for complex grids (2)
        int OuterIterations;

        /// @brief Number of boundary iterations in grid/net orthogonalization within itatp (25)
        int BoundaryIterations;

        /// @brief Number of inner iterations in grid/net orthogonalization within itbnd (25)
        int InnerIterations;

        /// @brief Factor from 0 to 1. between grid smoothing and grid orthogonality (0.975)
        double OrthogonalizationToSmoothingFactor;

        /// @brief Minimum ATPF on the boundary (1.0)
        double OrthogonalizationToSmoothingFactorBoundary;

        /// @brief Factor to weight between circumcentres 1.0 and masscentre 0.0 (1.0)
        double CircumCenterOrMassCenter;

        /// @brief Factor between smoother 1d0 and area-homogenizer 0d0 (1.0)
        double SmoothAngleOrSmoothArea;

        /// @brief Mesh2D-adaptation method; 0: Winslow, 1: arc-length, 2: harmonic map (1)
        int AdaptMethod;

        /// @brief Mesh2D-refinement factor; between 0d0 and 1d0 (0.0)
        double AdaptBeta;

        /// @brief Number of smoothing iterations of `solution` u in adaptation (0)
        int AdaptNiterU;

        /// @brief Number of smoothing iterations of monitor matrix G in adaptation (4)
        int AdaptNiterG;

        /// @brief Curvi-linear-like 0d0 or pure 1d0 orthogonalisation (0.5)
        double OrthoPure;
    };
} // namespace meshkernelapi
