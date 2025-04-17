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

#include <cmath>

#include "MeshKernel/Cartesian3DPoint.hpp"
#include "MeshKernel/Constants.hpp"

std::tuple<double, double, double> meshkernel::ComputeSphericalCoordinatesFromLatitudeAndLongitude(const Point& point)
{
    const double theta = (90.0 - point.y) * constants::conversion::degToRad;
    const double phi = point.x * constants::conversion::degToRad;
    const double r = constants::geometric::earth_radius;

    const double x = r * std::sin(theta) * std::cos(phi);
    const double y = r * std::sin(theta) * std::sin(phi);
    const double z = r * std::cos(theta);
    return {x, y, z};
}

meshkernel::Cartesian3DPoint meshkernel::VectorProduct(const Cartesian3DPoint& a, const Cartesian3DPoint& b)
{
    return Cartesian3DPoint{
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x};
}

meshkernel::Point meshkernel::Cartesian3DToSpherical(const Cartesian3DPoint& cartesianPoint)
{
    const double r = constants::geometric::earth_radius;
    const double theta = std::acos(cartesianPoint.z / r);
    const double phi = std::atan2(cartesianPoint.y, cartesianPoint.x);

    const double latitude = 90.0 - theta * constants::conversion::radToDeg;
    const double longitude = phi * constants::conversion::radToDeg;

    return {longitude, latitude};
}

meshkernel::Cartesian3DPoint meshkernel::SphericalToCartesian3D(const Point& sphericalPoint)
{
    const auto [x, y, z] = meshkernel::ComputeSphericalCoordinatesFromLatitudeAndLongitude(sphericalPoint);
    Cartesian3DPoint result;
    result.x = x;
    result.y = y;
    result.z = z;
    return result;
}
