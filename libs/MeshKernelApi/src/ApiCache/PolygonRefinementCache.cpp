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

#include "MeshKernelApi/ApiCache/PolygonRefinementCache.hpp"

meshkernelapi::PolygonRefinementCache::PolygonRefinementCache(const std::vector<meshkernel::Point>& polyPoints,
                                                              const int firstIndex,
                                                              const int secondIndex,
                                                              const double edgeLength,
                                                              const std::vector<meshkernel::Point>& refinedPoints)
    : CachedPointValues(refinedPoints),
      m_polygonPoints(polyPoints),
      m_firstNodeIndex(firstIndex),
      m_secondNodeIndex(secondIndex),
      m_targetEdgeLength(edgeLength) {}

bool meshkernelapi::PolygonRefinementCache::ValidOptions(const std::vector<meshkernel::Point>& polyPoints,
                                                         const int firstIndex,
                                                         const int secondIndex,
                                                         const double edgeLength) const
{
    return firstIndex == m_firstNodeIndex &&
           secondIndex == m_secondNodeIndex &&
           edgeLength == m_targetEdgeLength &&
           polyPoints.size() == m_polygonPoints.size() &&
           std::equal(polyPoints.begin(), polyPoints.end(), m_polygonPoints.begin());
}
