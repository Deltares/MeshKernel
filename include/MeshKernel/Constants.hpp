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

#include <cmath>
#include <limits>
#include <math.h>

namespace meshkernel
{
    namespace constants
    {
        // missing values
        namespace missing
        {
            constexpr double innerOuterSeparator = -998.0;                    ///< Double value used to separate the inner part of a polygon from its outer part
            constexpr double doubleValue = -999.0;                            ///< Double value used as missing value
            constexpr int intValue = -999;                                    ///< Integer value used as missing value
            constexpr size_t sizetValue = std::numeric_limits<size_t>::max(); ///< std::size_t missing value used for invalid indices
        }                                                                     // namespace missing
    }                                                                         // namespace constants

    // often used values
    static double const squareRootOfThree = std::sqrt(3.0); ///< The result of sqrt(3)
    constexpr double oneThird = 1.0 / 3.0;                  ///< The result of 1 / 3

    // geometric constants
    constexpr double degrad_hp = M_PI / 180.0;                   ///< Conversion factor from degrees to radians(pi / 180)
    constexpr double raddeg_hp = 180.0 / M_PI;                   ///< Conversion factor from radians to degrees(180 / pi)
    constexpr double earth_radius = 6378137.0;                   ///< Earth radius(m)
    constexpr double one_over_earth_radius = 1.0 / earth_radius; ///< One over earth_radius(m-1);
    constexpr double absLatitudeAtPoles = 0.0001;                ///< Pole tolerance in degrees

    // mesh constants
    constexpr double minimumDeltaCoordinate = 1e-14;                                       ///< Minimum delta coordinate
    constexpr std::size_t maximumNumberOfEdgesPerNode = 12;                                ///< Maximum number of edges per node
    constexpr std::size_t maximumNumberOfEdgesPerFace = 6;                                 ///< Maximum number of edges per face
    constexpr std::size_t maximumNumberOfNodesPerFace = 8;                                 ///< Maximum number of nodes per face
    constexpr std::size_t maximumNumberOfConnectedNodes = maximumNumberOfEdgesPerNode * 4; ///< Maximum number of connected nodes
    constexpr double minimumCellArea = 1e-12;                                              ///< Minimum cell area
    constexpr double weightCircumCenter = 1.0;                                             ///< Weight circum center
    constexpr std::size_t numNodesQuads = 4;                                               ///< Number of nodes in a quadrilateral
    constexpr std::size_t numNodesInTriangle = 3;                                          ///< Number of nodes in a triangle

    // orthogonalization
    constexpr double minimumEdgeLength = 1e-4;                   ///< Minimum edge length
    constexpr double curvilinearToOrthogonalRatio = 0.5;         ///< Ratio determining curvilinear-like(0.0) to pure(1.0) orthogonalization
    constexpr double orthogonalizationToSmoothingFactor = 0.975; ///< Factor between grid smoothing and grid ortho resp (0.<=ATPF<=1.)

} // namespace meshkernel
