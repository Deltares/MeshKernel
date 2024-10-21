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

#include "MeshKernelApi/CachedIntegerValues.hpp"
#include "MeshKernelApi/GeometryList.hpp"

namespace meshkernelapi
{

    /// @brief Cache node indices contained in a polygon
    class NodeInPolygonCache : public CachedIntegerValues
    {
    public:
        /// @brief Constructor
        NodeInPolygonCache(const std::vector<int>& nodeMask,
                           const std::vector<meshkernel::Point>& polygonPoints,
                           const int inside);

        /// @brief Determine if current options match those used to construct the object
        bool ValidOptions(const std::vector<meshkernel::Point>& polygonPoints, const int inside) const;

    private:
        /// @brief Points making up the polygon
        std::vector<meshkernel::Point> m_polygonPoints;

        /// @brief Indicates if the points are inside or outside of the polygon
        int m_inside = -1;
    };

} // namespace meshkernelapi
