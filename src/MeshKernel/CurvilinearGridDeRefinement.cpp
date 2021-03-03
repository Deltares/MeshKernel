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
    auto [mFirstNode, nFirstNode] = m_grid->GetNodeIndices(m_firstPoint);
    auto [mSecondNode, nSecondNode] = m_grid->GetNodeIndices(m_secondPoint);

    /// The points must lie on the same gridline
    if (mSecondNode - mFirstNode != 0 && nSecondNode - nFirstNode != 0)
    {
        throw std::invalid_argument("CurvilinearGridDeRefinement::Compute: The selected curvilinear grid nodes are not on the same grid-line");
    }

    /// estimate the dimension of the refined grid
    const auto numMToDeRefine = mSecondNode - mFirstNode == 0 ? 1 : mSecondNode - mFirstNode;
    const auto numNToDeRefine = nSecondNode - nFirstNode == 0 ? 1 : nSecondNode - nFirstNode;

    // the de-refined grid
    std::vector<std::vector<Point>> deRefinedGrid;
    deRefinedGrid.reserve(m_grid->m_numM);

    size_t localMDeRefinement;
    for (auto m = 0, mIndexOriginalGrid = 0; m < m_grid->m_numM - 1, mIndexOriginalGrid < m_grid->m_numM; ++m, mIndexOriginalGrid += localMDeRefinement)
    {
        localMDeRefinement = 1;
        if (mIndexOriginalGrid >= mFirstNode && mIndexOriginalGrid < mSecondNode)
        {
            localMDeRefinement = numMToDeRefine;
        }

        deRefinedGrid.emplace_back(std::vector<Point>());
        deRefinedGrid.back().reserve(m_grid->m_numN);
        size_t localNDeRefinement;
        for (auto n = 0, nIndexOriginalGrid = 0; n < m_grid->m_numN - 1, nIndexOriginalGrid < m_grid->m_numN; ++n, nIndexOriginalGrid += localNDeRefinement)
        {

            localNDeRefinement = 1;
            if (nIndexOriginalGrid >= nFirstNode && nIndexOriginalGrid < nSecondNode)
            {
                localNDeRefinement = numNToDeRefine;
            }
            deRefinedGrid.back().emplace_back(m_grid->m_gridNodes[mIndexOriginalGrid][nIndexOriginalGrid]);
        }
    }
    // substitute original grid with the derefined one
    return CurvilinearGrid(deRefinedGrid, m_grid->m_projection);
}