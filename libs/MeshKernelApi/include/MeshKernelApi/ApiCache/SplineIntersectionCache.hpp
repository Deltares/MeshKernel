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

#include <vector>

#include "MeshKernelApi/SplineIntersections.hpp"

namespace meshkernelapi
{

    /// @brief Cache spline intersection data
    class SplineIntersectionCache
    {
    public:
        /// @brief Get the number of spline intersections found
        int NumberOfIntersections() const;

        /// @brief Set the spline intersection data
        void Set(const std::vector<int>& splineIndices,
                 const std::vector<double>& intersectionAngles,
                 const std::vector<double>& intersectionCoordX,
                 const std::vector<double>& intersectionCoordY);

        /// @brief Copy the spline intersection data
        void Copy(SplineIntersections& intersections) const;

    private:
        std::vector<int> m_splineIndices;              ///< The indices of the spline intersected
        std::vector<double> m_intersectionAngles;      ///< The angles of the intersections
        std::vector<double> m_intersectionCoordinateX; ///< The x-coordinate of the intersection point
        std::vector<double> m_intersectionCoordinateY; ///< The y-coordinate of the intersection point
    };

} // namespace meshkernelapi
