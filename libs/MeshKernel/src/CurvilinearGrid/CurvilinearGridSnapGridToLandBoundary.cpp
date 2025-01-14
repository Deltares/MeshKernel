//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2025.
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

#include "MeshKernel/CurvilinearGrid/CurvilinearGridSnapGridToLandBoundary.hpp"

meshkernel::CurvilinearGridSnapGridToLandBoundary::CurvilinearGridSnapGridToLandBoundary(CurvilinearGrid& grid,
                                                                                         const LandBoundary& landBoundary,
                                                                                         const std::vector<Point>& points) : CurvilinearGridSnapping(grid, points),
                                                                                                                             m_projection(grid.projection()),
                                                                                                                             m_landBoundary(landBoundary)
{
}

std::unique_ptr<meshkernel::CurvilinearGridMeshExpansionCalculator> meshkernel::CurvilinearGridSnapGridToLandBoundary::AllocateCurvilinearGridMeshExpansionCalculator(const CurvilinearGrid& originalGrid,
                                                                                                                                                                      const CurvilinearGrid& snappedGrid) const
{
    return std::make_unique<DefaultRegionExpasionCalculator>(originalGrid, snappedGrid, m_landBoundary);
}

meshkernel::Point meshkernel::CurvilinearGridSnapGridToLandBoundary::FindNearestPoint(const Point& currentPoint) const
{
    return m_landBoundary.FindNearestPoint(currentPoint, m_projection);
}
