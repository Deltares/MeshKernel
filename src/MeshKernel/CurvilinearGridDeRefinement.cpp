//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2021.
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

#include <utility>

#include <MeshKernel/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGridDeRefinement.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Splines.hpp>

meshkernel::CurvilinearGridDeRefinement::CurvilinearGridDeRefinement(std::shared_ptr<CurvilinearGrid> grid, const Point& firstPoint, const Point& secondPoint)
    : m_grid(std::move(grid)),
      m_firstPoint(firstPoint),
      m_secondPoint(secondPoint)
{
}

meshkernel::CurvilinearGrid meshkernel::CurvilinearGridDeRefinement::Compute()
{
    // Get the m and n indices from the point coordinates
    auto const firstNode = m_grid->GetNodeIndices(m_firstPoint);
    auto const secondNode = m_grid->GetNodeIndices(m_secondPoint);

    /// The points must lie on the same gridline
    if (secondNode.m - firstNode.m != 0 && secondNode.n - firstNode.n != 0)
    {
        throw std::invalid_argument("CurvilinearGridDeRefinement::Compute: The selected curvilinear grid nodes are not on the same grid-line");
    }

    /// estimate the dimension of the refined grid
    const auto numMToDeRefine = secondNode.m > firstNode.m ? secondNode.m - firstNode.m : 1;
    const auto numNToDeRefine = secondNode.n > firstNode.n ? secondNode.n - firstNode.n : 1;

    // the de-refined grid
    std::vector<std::vector<Point>> deRefinedGrid;
    deRefinedGrid.reserve(m_grid->m_numM);

    size_t mIndexOriginalGrid = 0;
    while (mIndexOriginalGrid < m_grid->m_numM)
    {
        size_t localMDeRefinement = 1;
        if (mIndexOriginalGrid >= firstNode.m && mIndexOriginalGrid < secondNode.m)
        {
            localMDeRefinement = numMToDeRefine;
        }
        deRefinedGrid.emplace_back(std::vector<Point>());
        deRefinedGrid.back().reserve(m_grid->m_numN);

        size_t nIndexOriginalGrid = 0;
        while (nIndexOriginalGrid < m_grid->m_numN)
        {
            size_t localNDeRefinement = 1;
            if (nIndexOriginalGrid >= firstNode.n && nIndexOriginalGrid < secondNode.n)
            {
                localNDeRefinement = numNToDeRefine;
            }
            deRefinedGrid.back().emplace_back(m_grid->m_gridNodes[mIndexOriginalGrid][nIndexOriginalGrid]);
            nIndexOriginalGrid += localNDeRefinement;
        }
        mIndexOriginalGrid += localMDeRefinement;
    }
    // substitute original grid with the derefined one
    return CurvilinearGrid(std::move(deRefinedGrid), m_grid->m_projection);
}
