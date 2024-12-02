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

#include <concepts>
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

    /// @brief Gets the valid deletion options as vector of integers
    const std::vector<int>& GetValidDeletionOptions();

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

    /// @enum InterpolationDataTypes
    /// @brief The possible types of the values to be interpolated in the gridded sample
    enum class InterpolationDataTypes
    {
        Short = 0,  ///< short type
        Float = 1,  ///< float type
        Int = 2,    ///< int type
        Double = 3, ///< double type
    };

    /// @brief Convert an integer value to the CurvilinearDirection enumeration type
    ///
    /// If the integer direction value does not correspond to an enumeration
    /// value then a ConstraintError will be thrown
    CurvilinearDirection GetCurvilinearDirectionValue(int direction);

    /// @brief Get the string representation of the CurvilinearDirection enumeration values.
    const std::string& CurvilinearDirectionToString(CurvilinearDirection direction);

    /// @brief The concept specifies that the array type must have an access operator returning the array element type or can be converted to one
    template <typename ArrayType, typename ResultType>
    concept ArrayConstAccessConcept = requires(const ArrayType& array, const size_t i) {
        { array[i] } -> std::convertible_to<ResultType>;
    };

    /// @brief The concept specifies that the array type must have an access operator returning a reference to the array element type
    template <typename ArrayType, typename ResultType>
    concept ArrayNonConstAccessConcept = requires(ArrayType& array, const size_t i) {
        { array[i] } -> std::same_as<ResultType&>;
    };

    /// @brief A concept that specifies that the array must have a size function return the number of elements in the array
    template <typename ArrayType>
    concept ArraySizeConcept = requires(const ArrayType& array) {
        { array.size() } -> std::same_as<size_t>;
    };

    /// @brief A concept that specifies that the array must have a begin and end function.
    ///
    /// Would like to also specify the return type here, but span needs some c++23 functionality here.
    /// Then change all iterator usage to cbegin and cend returning a const_iterator
    /// std::same_as<typename ArrayType::const_iterator>
    template <typename ArrayType>
    concept ArrayConstIteratorsConcept = requires(const ArrayType& array) {
        { array.begin() };
        { array.end() };
    };

    /// @brief A concept that specifies that the array must have a begin and end function.
    ///
    /// Would like to also specify the return type here, but span needs some c++23 functionality here.
    /// Then change all iterator usage to cbegin and cend returning a const_iterator
    /// std::same_as<typename ArrayType::const_iterator>
    template <typename ArrayType>
    concept ArrayNonConstIteratorsConcept = requires(ArrayType& array) {
        { array.begin() } -> std::same_as<typename ArrayType::iterator>;
        { array.end() } -> std::same_as<typename ArrayType::iterator>;
    };

    /// @brief A concept that specifies all the functionality required to be usable as a constant array of doubles.
    template <typename ArrayType>
    concept ValidConstDoubleArray = ArrayConstAccessConcept<ArrayType, double> &&
                                    ArrayConstIteratorsConcept<ArrayType> &&
                                    ArraySizeConcept<ArrayType>;

    /// @brief A concept that specifies all the functionality required to be usable as a constant array of doubles.
    template <typename ArrayType>
    concept ValidNonConstDoubleArray = ArrayNonConstAccessConcept<ArrayType, double> &&
                                       ArrayNonConstIteratorsConcept<ArrayType> &&
                                       ArraySizeConcept<ArrayType>;

} // namespace meshkernel
