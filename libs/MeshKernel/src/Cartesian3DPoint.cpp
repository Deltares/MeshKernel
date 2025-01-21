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

meshkernel::Cartesian3DPoint meshkernel::VectorProduct(const Cartesian3DPoint& a, const Cartesian3DPoint& b)
{
    return Cartesian3DPoint{
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x};
}

meshkernel::Cartesian3DPoint meshkernel::SphericalToCartesian3D(const Point& sphericalPoint)
{
    Cartesian3DPoint result;
    result.z = constants::geometric::earth_radius * std::sin(sphericalPoint.y * constants::conversion::degToRad);
    const double rr = constants::geometric::earth_radius * std::cos(sphericalPoint.y * constants::conversion::degToRad);
    result.x = rr * std::cos(sphericalPoint.x * constants::conversion::degToRad);
    result.y = rr * std::sin(sphericalPoint.x * constants::conversion::degToRad);
    return result;
}

meshkernel::Point meshkernel::Cartesian3DToSpherical(const Cartesian3DPoint& cartesianPoint, double referenceLongitude)
{
    Point sphericalPoint;
    const double angle = std::atan2(cartesianPoint.y, cartesianPoint.x) * constants::conversion::radToDeg;
    sphericalPoint.y = std::atan2(cartesianPoint.z, sqrt(cartesianPoint.x * cartesianPoint.x + cartesianPoint.y * cartesianPoint.y)) * constants::conversion::radToDeg;
    sphericalPoint.x = angle + static_cast<double>(std::lround((referenceLongitude - angle) / 360.0)) * 360.0;
    return sphericalPoint;
}
