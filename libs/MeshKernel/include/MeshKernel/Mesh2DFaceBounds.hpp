//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2026.
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

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Mesh.hpp"
#include "MeshKernel/Point.hpp"

namespace meshkernel::algo
{

    /// @brief Compute the nodes for each face in the grid, will be saved as facebounds
    class Mesh2DFaceBounds
    {
    public:
        /// @brief Compute the nodes for each face in the grid
        static std::vector<Point> Compute(const Mesh& mesh);

    private:
        /// @brief Compute the face bounds for a single face
        static void ComputeBoundsForFace(const Mesh& mesh, UInt faceId, std::span<Point> faceBounds);
    };

} // namespace meshkernel::algo
