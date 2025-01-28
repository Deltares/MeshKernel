//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2025.
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

#include <span>
#include <vector>

#include "MeshKernel/Mesh2D.hpp"

namespace meshkernel
{
    /// @brief Compute the smoothness value for the edges
    class MeshSmoothness
    {
    public:
        /// @brief Compute the smoothness values returning values in a vector
        static std::vector<double> Compute(const Mesh2D& mesh);

        /// @brief Compute the smoothness values overwriting the values in an array
        static void Compute(const Mesh2D& mesh, std::span<double> smoothness);

    private:
        static constexpr double m_minimumCellArea = 1e-12; ///< Minimum cell area
    };

} // namespace meshkernel
