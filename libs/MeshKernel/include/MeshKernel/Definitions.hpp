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

#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace meshkernel
{
    /// @brief Integer type used when indexing mesh graph entities.
    using UInt = std::uint32_t;

    /// @brief Enumerator describing the supported projections
    enum class Projection
    {
        cartesian = 0,        // jsferic  = 0
        spherical = 1,        // jsferic  = 1
        sphericalAccurate = 2 // jasfer3D = 1
    };

    /// @brief Gets the valid projectionbs as vector of integers
    const std::vector<int>& GetValidProjections();

    /// @brief Convert an integer value to the Projection enumeration type
    ///
    /// If the integer projection value does not correspond to an enumeration
    /// value then a ConstraintError will be thrown
    Projection GetProjectionValue(int projection);

    /// @brief Get the string representation of the Projection enumeration values.
    const std::string& ProjectionToString(Projection projection);

    /// @brief Indicator for traversal direction of the points specifying a polygon
    // PolygonTraversalDirection? too long
    // PolygonOrientation
    enum class TraversalDirection
    {
        Clockwise,    ///< Points define a clockwise traversal of the polygon
        AntiClockwise ///< Points define a anti-clockwise (counter-clockwise) traversal of the polygon
    };

    /// @enum Location
    /// @brief Mesh locations enumeration
    enum class Location
    {
        Faces = 0,  ///< Faces
        Nodes = 1,  ///< Nodes
        Edges = 2,  ///< Edges
        Unknown = 3 ///< Unknown
    };

    /// @brief Maps Location enumeration to a string
    inline static std::map<Location, std::string> const LocationToString = {
        {Location::Faces, "Faces"},
        {Location::Nodes, "Nodes"},
        {Location::Edges, "Edges"},
        {Location::Unknown, "Unknown"}};

    /// @brief Direction to use in curvilinear grid algorithms
    enum class CurvilinearDirection
    {
        M, ///< M-direction
        N  ///< N-direction
    };

    /// @enum InterpolationValuesTypes
    /// @brief The possible types of the values to be interpolated in the gridded sample 
    enum class InterpolationValuesTypes
    {
        shortType = 0, ///< short type
        intType = 1,   ///< int type
        floatType = 2, ///< float type
        doubleType = 3 ///< double type
    };

    /// @brief Convert an integer value to the CurvilinearDirection enumeration type
    ///
    /// If the integer direction value does not correspond to an enumeration
    /// value then a ConstraintError will be thrown
    CurvilinearDirection GetCurvilinearDirectionValue(int direction);

    /// @brief Get the string representation of the CurvilinearDirection enumeration values.
    const std::string& CurvilinearDirectionToString(CurvilinearDirection direction);

} // namespace meshkernel
