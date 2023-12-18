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
#include <array>
#include <utility>

#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/Point.hpp"

namespace meshkernel
{
    /// @brief Construct a global grid in spherical coordinates, as a base for later mesh refinements
    class Mesh2DGenerateGlobalGrid
    {
    public:
        /// @brief Compute the grid
        ///
        /// @param [in] numX The Number of x points
        /// @param [in] numY The Number of y points
        /// @param [in] polygon The polygon where to generate the mesh
        static std::unique_ptr<Mesh2D> Compute(const UInt numX, const UInt numY, const Polygons& polygon);

    private:
        /// @brief Compute the latitude increment given the current latitude and the longitude discretization
        static double getDeltaLatitude(const double currentLatitude, const double longitudeDiscretization);

        /// @brief Gets the node index from a give position
        static UInt getNodeIndexFromPosition(const Mesh& mesh, const Point& position);

        /// @brief Add a face to an existing mesh towards a specific direction
        static void addFace(Mesh& mesh, const std::array<Point, 8>& points, const double latitudeDirection, const UInt numNodes);
    };
} // namespace meshkernel
