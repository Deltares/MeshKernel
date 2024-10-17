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

#include <utility>

#include "MeshKernelApi/FacePolygonPropertyCache.hpp"
#include "MeshKernelApi/Utils.hpp"

meshkernelapi::FacePolygonPropertyCache::FacePolygonPropertyCache(const int propertyValue,
                                                                  const double minValue,
                                                                  const double maxValue,
                                                                  const meshkernel::Mesh2D& mesh,
                                                                  const int validSize,
                                                                  const std::vector<bool>& filterMask)
    : m_propertyValue(propertyValue),
      m_minimumValue(minValue),
      m_maximumValue(maxValue)
{
    std::vector<double> xCoordinates(validSize);
    std::vector<double> yCoordinates(validSize);

    GeometryList facePolygons;

    facePolygons.coordinates_x = xCoordinates.data();
    facePolygons.coordinates_y = yCoordinates.data();

    FillFacePolygons(mesh, filterMask, facePolygons);

    CachedPointValues::Reset(std::move(xCoordinates), std::move(yCoordinates));
}

bool meshkernelapi::FacePolygonPropertyCache::ValidOptions(const int propertyValue,
                                                           const double minValue,
                                                           const double maxValue) const
{
    return propertyValue == m_propertyValue &&
           minValue == m_minimumValue &&
           maxValue == m_maximumValue;
}
