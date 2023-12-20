//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2023.
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

#include "MeshKernel/Utilities/LinearAlgebra.hpp"

#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Point.hpp"

#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"

namespace meshkernel
{

    /// @brief Compute the smoothness of a grid.
    class CurvilinearGridSmoothness
    {
    public:
        /// @brief Compute the smoothness of the grid in the direction specified.
        ///
        /// @param grid Grid for which the smoothness is to be calculated.
        /// @param direction Direction of the smoothness required.
        /// @param smoothness Matrix of smoothness values, will be resized to the size of the grid.
        static void Compute(const CurvilinearGrid& grid, const CurvilinearDirection direction, lin_alg::Matrix<double>& smoothness);

    private:
        /// @brief Compute the smoothness for the node.
        static double ComputeNodeSmoothness(const Point& p0, const Point& p1, const Point& p2);
    };

} // namespace meshkernel
