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

#include <algorithm>
#include <cstring>
#include <utility>

#include "MeshKernel/Exceptions.hpp"

#include "MeshKernelApi/ApiCache/NodeInPolygonCache.hpp"

meshkernelapi::NodeInPolygonCache::NodeInPolygonCache(const std::vector<int>& nodeMask,
                                                      const std::vector<meshkernel::Point>& polygonPoints,
                                                      const int inside)
    : m_polygonPoints(polygonPoints), m_inside(inside)
{
    std::vector<int> nodeIndices;

    nodeIndices.reserve(m_polygonPoints.size());

    for (size_t i = 0; i < nodeMask.size(); ++i)
    {
        if (nodeMask[i] > 0)
        {
            nodeIndices.push_back(static_cast<int>(i));
        }
    }

    Reset(std::move(nodeIndices));
}

bool meshkernelapi::NodeInPolygonCache::ValidOptions(const std::vector<meshkernel::Point>& polygonPoints, const int inside) const
{
    return inside == m_inside &&
           polygonPoints.size() == m_polygonPoints.size() &&
           std::equal(polygonPoints.begin(), polygonPoints.end(), m_polygonPoints.begin());
}
