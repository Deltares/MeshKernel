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

#include <cstring>
#include <utility>

#include "MeshKernelApi/CachedPointValues.hpp"

meshkernelapi::CachedPointValues::CachedPointValues(const std::vector<meshkernel::Point>& coordinates)
    : m_coordsX(coordinates.size()),
      m_coordsY(coordinates.size())
{
    for (size_t i = 0; i < coordinates.size(); ++i)
    {
        m_coordsX[i] = coordinates[i].x;
        m_coordsY[i] = coordinates[i].y;
    }
}

void meshkernelapi::CachedPointValues::Copy(const GeometryList& geometry) const
{
    size_t valueCount = sizeof(double) * m_coordsX.size();

    std::memcpy(geometry.coordinates_x, m_coordsX.data(), valueCount);
    std::memcpy(geometry.coordinates_y, m_coordsY.data(), valueCount);
}

void meshkernelapi::CachedPointValues::Reset(std::vector<double>&& xValues,
                                             std::vector<double>&& yValues)
{
    m_coordsX = std::move(xValues);
    m_coordsY = std::move(yValues);
}
