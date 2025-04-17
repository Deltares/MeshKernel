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

#include "MeshKernel/Point.hpp"

namespace meshkernel
{

    /// @brief A struct describing the three coordinates in a cartesian projection.
    struct Cartesian3DPoint
    {
        double x; ///< X-coordinate
        double y; ///< Y-coordinate
        double z; ///< Z-coordinate
    };

    /// @brief Converts geographic coordinates (latitude and longitude) to 3D spherical Cartesian coordinates.
    /// @param[in] point The input point specified in latitude and longitude (in degrees).
    /// @return A tuple containing the corresponding (x, y, z) coordinates on the sphere.
    std::tuple<double, double, double> ComputeSphericalCoordinatesFromLatitudeAndLongitude(const Point& point);

    /// @brief Defines vector product for cartesian 3D-space
    /// @param[in] a The first cartesian 3D point
    /// @param[in] b The second cartesian 3D point
    /// @return The vector product
    [[nodiscard]] Cartesian3DPoint VectorProduct(const Cartesian3DPoint& a, const Cartesian3DPoint& b);

    /// @brief Defines inner product in cartesian 3D-space
    /// @param[in] a The first cartesian 3D point
    /// @param[in] b The second cartesian 3D point
    /// @return The resulting inner product
    [[nodiscard]] double InnerProduct(const Cartesian3DPoint& a, const Cartesian3DPoint& b);

    /// @brief Multiply Cartesian point p by a scalar value
    Cartesian3DPoint operator*(const Cartesian3DPoint& p, const double value);

    /// @brief Multiply Cartesian point p by a scalar value
    Cartesian3DPoint operator*(const double value, const Cartesian3DPoint& p);

    /// @brief Divide Cartesian point p by a scalar value
    Cartesian3DPoint operator/(const Cartesian3DPoint& p, const double value);

    /// @brief Add Cartesian point p2 to p1.
    Cartesian3DPoint operator+(const Cartesian3DPoint& p1, const Cartesian3DPoint& p2);

    /// @brief Subtract Cartesian point p2 from p1.
    Cartesian3DPoint operator-(const Cartesian3DPoint& p1, const Cartesian3DPoint& p2);

    /// @brief Transforms 2D point in spherical coordinates to 3D cartesian coordinates.
    /// @param[in] sphericalPoint The current spherical point (2 coordinates).
    /// @returns The converted cartesian 3d point.
    [[nodiscard]] Cartesian3DPoint SphericalToCartesian3D(const Point& sphericalPoint);

    /// @brief Transforms 3D cartesian coordinates to 2D point in spherical coordinates
    /// @param[in] cartesianPoint The 3d cartesian point
    /// @returns The spherical coordinate
    [[nodiscard]] Point Cartesian3DToSpherical(const Cartesian3DPoint& cartesianPoint);

} // end namespace meshkernel

inline meshkernel::Cartesian3DPoint meshkernel::operator+(const Cartesian3DPoint& p1, const Cartesian3DPoint& p2)
{
    return {p1.x + p2.x, p1.y + p2.y, p1.z + p2.z};
}

inline meshkernel::Cartesian3DPoint meshkernel::operator-(const Cartesian3DPoint& p1, const Cartesian3DPoint& p2)
{
    return {p1.x - p2.x, p1.y - p2.y, p1.z - p2.z};
}

inline meshkernel::Cartesian3DPoint meshkernel::operator/(const Cartesian3DPoint& p, const double value)
{
    return {p.x / value, p.y / value, p.z / value};
}

inline meshkernel::Cartesian3DPoint meshkernel::operator*(const double value, const Cartesian3DPoint& p)
{
    return {value * p.x, value * p.y, value * p.z};
}

inline meshkernel::Cartesian3DPoint meshkernel::operator*(const Cartesian3DPoint& p, const double value)
{
    return value * p;
}

inline double meshkernel::InnerProduct(const Cartesian3DPoint& a, const Cartesian3DPoint& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
