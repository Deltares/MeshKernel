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

#pragma once

#include <cstring>
#include <vector>

#include "MeshKernel/Point.hpp"

#include "MeshKernelApi/GeometryList.hpp"

namespace meshkernelapi
{

    class NodeInPolygonCache
    {
    public:
        NodeInPolygonCache(const std::vector<int>& nodeMask,
                           const std::vector<meshkernel::Point>& polygonPoints,
                           const int inside);

        bool ValidOptions(const std::vector<meshkernel::Point>& polygonPoints, const int inside) const;

        int Size() const;

        void Copy(int* selectedNodes) const;

    private:
        /// &brief Points making up the polygon
        std::vector<meshkernel::Point> m_polygonPoints;

        /// &brief Indicates if the points are inside or outside of the polygon
        int m_inside = -1;

        /// &brief Indices of nodes in the polygon
        std::vector<int> m_nodeIndices;
    };

} // namespace meshkernelapi

inline int meshkernelapi::NodeInPolygonCache::Size() const
{
    return static_cast<int>(m_nodeIndices.size());
}
