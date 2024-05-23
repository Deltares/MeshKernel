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
#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/LandBoundary.hpp"
#include "MeshKernel/Point.hpp"

namespace meshkernel
{

    /// @brief Smoothly snap the grid to a land boundary.
    class CurvilinearGridSnapGridToLandBoundary : public CurvilinearGridSnapping
    {
    public:
        /// @brief constructor
        /// @param [in] grid         The input curvilinear grid
        /// @param [in] landBoundary The land boundary to which the grid is to be snapped.
        /// @param [in] points       The points used to control the snapping and expansion.
        CurvilinearGridSnapGridToLandBoundary(CurvilinearGrid& grid,
                                              const LandBoundary& landBounday,
                                              const std::vector<Point>& points);

    private:
        /// @brief Allocate the grid expansion calculator
        std::unique_ptr<CurvilinearGridMeshExpansionCalculator>
        AllocateCurvilinearGridMeshExpansionCalculator(const CurvilinearGrid& originalGrid,
                                                       const CurvilinearGrid& snappedGrid) const override;

        /// @brief Find the nearest point to the land boundary
        Point FindNearestPoint(const Point& currentPoint) const override;

        /// @brief Coordinate system
        const Projection m_projection;

        /// @brief The land boundary to which the grid is to be snapped.
        const LandBoundary& m_landBoundary;
    };

} // namespace meshkernel
