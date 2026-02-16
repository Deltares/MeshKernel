//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2024.
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

#include "MeshKernelApi/ApiCache/SplineIntersectionCache.hpp"

int meshkernelapi::SplineIntersectionCache::NumberOfIntersections() const
{
    return static_cast<int>(m_splineIndices.size());
}

void meshkernelapi::SplineIntersectionCache::Set(const std::vector<int>& splineIndices,
                                                 const std::vector<double>& intersectionAngles,
                                                 const std::vector<double>& intersectionCoordX,
                                                 const std::vector<double>& intersectionCoordY)
{
    m_splineIndices = splineIndices;
    m_intersectionAngles = intersectionAngles;
    m_intersectionCoordinateX = intersectionCoordX;
    m_intersectionCoordinateY = intersectionCoordY;
}

void meshkernelapi::SplineIntersectionCache::Copy(SplineIntersections& intersections) const
{
    if (NumberOfIntersections() > 0)
    {
        std::memcpy(intersections.spline_index, m_splineIndices.data(), sizeof(int) * NumberOfIntersections());
        std::memcpy(intersections.intersection_angle, m_intersectionAngles.data(), sizeof(double) * NumberOfIntersections());
        std::memcpy(intersections.intersection_x, m_intersectionCoordinateX.data(), sizeof(double) * NumberOfIntersections());
        std::memcpy(intersections.intersection_y, m_intersectionCoordinateY.data(), sizeof(double) * NumberOfIntersections());
    }
}
