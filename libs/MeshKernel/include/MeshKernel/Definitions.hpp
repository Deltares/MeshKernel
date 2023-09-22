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

#include <map>
#include <string>

#include <cstdint>

namespace meshkernel
{
    /// @brief Integer type used when indexing mesh graph entities.
    using UInt = std::uint32_t;

    /// @enum Location
    /// @brief Mesh locations enumeration
    enum class MeshLocation
    {
        Faces = 0,  ///< Faces
        Nodes = 1,  ///< Nodes
        Edges = 2,  ///< Edges
        Unknown = 3 ///< Unknown
    };

    /// @brief Maps Location enumeration to a string
    inline static std::map<MeshLocation, std::string> const LocationToString = {
        {MeshLocation::Faces, "Faces"},
        {MeshLocation::Nodes, "Nodes"},
        {MeshLocation::Edges, "Edges"},
        {MeshLocation::Unknown, "Unknown"}};

    /// @brief Indicator for traversal direction of the points specifying a polygon
    // PolygonTraversalDirection? too long
    // PolygonOrientation
    enum class TraversalDirection
    {
        Clockwise,    ///< Points define a clockwise traversal of the polygon
        AntiClockwise ///< Points define a anti-clockwise (counter-clockwise) traversal of the polygon
    };

} // namespace meshkernel
