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

#include <memory>
#include <vector>

#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridMeshExpansionCalculator.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridSnapping.hpp"
#include "MeshKernel/Point.hpp"
#include "MeshKernel/Splines.hpp"

namespace meshkernel
{

    /// @brief Smoothly snap the grid to a spline.
    class CurvilinearGridSnapGridToSpline : public CurvilinearGridSnapping
    {
    public:
        /// @brief constructor
        /// @param [in] grid   The input curvilinear grid
        /// @param [in] spline The spline to which the grid is to be snapped.
        /// @param [in] points The points used to control the snapping and expansion.
        CurvilinearGridSnapGridToSpline(CurvilinearGrid& grid,
                                        const Splines& spline,
                                        const std::vector<Point>& points);

    private:
        /// @brief Allocate a mesh expansion calculator
        std::unique_ptr<CurvilinearGridMeshExpansionCalculator>
        AllocateCurvilinearGridMeshExpansionCalculator(const CurvilinearGrid& originalGrid,
                                                       const CurvilinearGrid& snappedGrid) const override;

        /// @brief Find the nearest point to the spline
        Point FindNearestPoint(const Point& currentPoint) const override;

        /// @brief The spline
        const Splines& m_spline;
    };

} // namespace meshkernel
