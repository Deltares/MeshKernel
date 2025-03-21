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

#include "MeshKernelApi/ApiCache/MeshBoundariesAsPolygonCache.hpp"

meshkernelapi::MeshBoundariesAsPolygonCache::MeshBoundariesAsPolygonCache(const std::vector<meshkernel::Point>& selectionPolygonPoints,
                                                                          const std::vector<meshkernel::Point>& boundaryPoints)
    : CachedPointValues(boundaryPoints), m_selectionPolygonPoints(selectionPolygonPoints)
{
}

bool meshkernelapi::MeshBoundariesAsPolygonCache::ValidOptions(const std::vector<meshkernel::Point>& selectionPolygonPoints) const
{

    if (selectionPolygonPoints.size() != m_selectionPolygonPoints.size())
    {
        return false;
    }

    constexpr double tolerance = 1e-9;
    for (auto i = 0u; i < selectionPolygonPoints.size(); ++i)
    {
        if (!IsEqual(selectionPolygonPoints[i], m_selectionPolygonPoints[i], tolerance))
        {
            return false;
        }
    }
    return true;
}
