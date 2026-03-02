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

#include <array>
#include <span>
#include <vector>

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Mesh.hpp"
#include "MeshKernel/Point.hpp"

namespace meshkernel::algo
{

    /// @brief Compute the polygon surrounding an edge with width to the circumcentres of the connected elements
    class NetlinkContourPolygons
    {
    public:
        /// @brief Compute the netlink contour polygons
        std::vector<Point> Compute(const Mesh& mesh) const;

    private:
        /// @brief Compute the netlink contour polygon for a single edge
        ///
        /// @note circumcentre1 must be a valid circumcentre point.
        void ComputePolygonForEdge(const Point& start, const Point& end, const Point& curcumcentre1, const Point& circumcentre2,
                                   const Projection projection,
                                   std::span<Point> polygon) const;
    };

} // namespace meshkernel::algo
