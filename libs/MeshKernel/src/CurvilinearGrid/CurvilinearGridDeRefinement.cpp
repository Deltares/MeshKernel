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

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridDeRefinement.hpp>
#include <MeshKernel/CurvilinearGrid/UndoActions/CurvilinearGridRefinementUndoAction.hpp>
#include <MeshKernel/Entities.hpp>

using meshkernel::CurvilinearGrid;
using meshkernel::CurvilinearGridDeRefinement;

CurvilinearGridDeRefinement::CurvilinearGridDeRefinement(CurvilinearGrid& grid, int derefinementFactor) : CurvilinearGridAlgorithm(grid), m_derefinementFactor(derefinementFactor)
{
}

meshkernel::UndoActionPtr CurvilinearGridDeRefinement::Compute()
{
    if (!m_lowerLeft.IsValid() || !m_upperRight.IsValid())
    {
        throw std::invalid_argument("CurvilinearGridDeRefinement::Compute: lower left and upper right corners defining the curvilinear grid block are not set");
    }

    // estimate the dimension of the refined grid
    const auto numMToDeRefine = m_upperRight.m_m > m_lowerLeft.m_m ? m_upperRight.m_m - m_lowerLeft.m_m : 1;
    const auto numNToDeRefine = m_upperRight.m_n > m_lowerLeft.m_n ? m_upperRight.m_n - m_lowerLeft.m_n : 1;

    const int mDeRefineFactor = numMToDeRefine > 1 ? m_derefinementFactor : 1;
    const int nDeRefineFactor = numNToDeRefine > 1 ? m_derefinementFactor : 1;

    std::unique_ptr<CurvilinearGridRefinementUndoAction> undoAction = CurvilinearGridRefinementUndoAction::Create(m_grid);

    // the de-refined grid
    std::vector<std::vector<Point>> deRefinedGrid;
    deRefinedGrid.reserve(m_grid.NumN());

    UInt nCurrentIndex = 0;
    UInt nLastUsedIndex = 0;
    while (nCurrentIndex < m_grid.NumN())
    {
        int localNDeRefinement = 1;
        if (nCurrentIndex >= m_lowerLeft.m_n && nCurrentIndex < m_upperRight.m_n)
        {
            localNDeRefinement = nDeRefineFactor;
            if (nCurrentIndex + localNDeRefinement > m_upperRight.m_n)
            {
                localNDeRefinement = static_cast<int>(m_upperRight.m_n) -
                                     static_cast<int>(nCurrentIndex);
            }
        }
        deRefinedGrid.emplace_back(std::vector<Point>());
        deRefinedGrid.back().reserve(m_grid.NumM());
        nLastUsedIndex = nCurrentIndex;

        UInt mCurrentIndex = 0;
        UInt mLastUsedIndex = 0;
        while (mCurrentIndex < m_grid.NumM())
        {
            int localMDeRefinement = 1;
            if (mCurrentIndex >= m_lowerLeft.m_m && mCurrentIndex < m_upperRight.m_m)
            {
                localMDeRefinement = mDeRefineFactor;
                if (mCurrentIndex + localMDeRefinement > m_upperRight.m_m)
                {
                    localMDeRefinement = static_cast<int>(m_upperRight.m_m) -
                                         static_cast<int>(mCurrentIndex);
                }
            }
            const auto currentNode = m_grid.GetNode(nCurrentIndex, mCurrentIndex);
            mLastUsedIndex = mCurrentIndex;
            deRefinedGrid.back().emplace_back(currentNode);
            mCurrentIndex += localMDeRefinement;
        }

        // Emplace the last column to avoid reducing the curvilinear grid size
        if (mLastUsedIndex < m_grid.NumM() - 1)
        {
            deRefinedGrid.back().emplace_back(m_grid.GetNode(nCurrentIndex, m_grid.NumM() - 1));
        }
        nCurrentIndex += localNDeRefinement;
    }

    // Emplace the last row to avoid reducing the grid size
    if (nLastUsedIndex < m_grid.NumN() - 1)
    {
        deRefinedGrid.emplace_back(std::vector<Point>());
        for (auto m = 0u; m < m_grid.NumM(); ++m)
        {
            deRefinedGrid.back().emplace_back(m_grid.GetNode(m_grid.NumN() - 1, m));
        }
    }

    // Substitute original grid with the derefined one
    m_grid.SetGridNodes(lin_alg::STLVectorOfVectorsToMatrix(deRefinedGrid));

    return undoAction;
}
