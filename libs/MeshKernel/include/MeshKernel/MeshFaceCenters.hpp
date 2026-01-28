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
    /// @brief Compute the circum-center point of each of the faces
    std::vector<Point> ComputeFaceCircumcenters(const Mesh& mesh, CircumCentreMethod circumcentreMethod = constants::geometric::defaultCircumCentreMethod,
                                                double circumCentreWeight = constants::geometric::circumCentreWeight);

    /// @brief Compute the circum-center point of each of the faces overwriting the values in an array
    void ComputeFaceCircumcenters(const Mesh& mesh, std::span<Point> faceCenters,
                                  CircumCentreMethod circumcentreMethod = constants::geometric::defaultCircumCentreMethod,
                                  double circumCentreWeight = constants::geometric::circumCentreWeight);

    /// @brief Compute the circumcenter of a polygon
    Point ComputeFaceCircumenter(std::vector<Point>& polygon,
                                 const std::vector<UInt>& edgesNumFaces,
                                 const Projection projection,
                                 double circumCentreWeight = constants::geometric::circumCentreWeight,
                                 CircumCentreMethod circumcentreMethod = constants::geometric::defaultCircumCentreMethod);

    /// @brief Compute the mass centre element.
    Point ComputeMassCentre(const std::vector<Point>& polygon);

} // namespace meshkernel::algo
