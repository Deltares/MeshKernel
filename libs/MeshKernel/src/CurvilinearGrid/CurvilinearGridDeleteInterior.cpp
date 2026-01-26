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

#include "MeshKernel/CurvilinearGrid/CurvilinearGridDeleteInterior.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridLine.hpp"
#include "MeshKernel/CurvilinearGrid/UndoActions/CurvilinearGridBlockUndoAction.hpp"

using meshkernel::CurvilinearGrid;
using meshkernel::CurvilinearGridLine;
using meshkernel::CurvilinearGridNodeIndices;

meshkernel::CurvilinearGridDeleteInterior::CurvilinearGridDeleteInterior(CurvilinearGrid& grid)
    : CurvilinearGridAlgorithm(grid)
{
}

meshkernel::UndoActionPtr meshkernel::CurvilinearGridDeleteInterior::Compute()
{
    const UInt lowerLimitI = m_lowerLeft.m_n;
    const UInt upperLimitI = m_upperRight.m_n;

    const UInt lowerLimitJ = m_lowerLeft.m_m;
    const UInt upperLimitJ = m_upperRight.m_m;

    std::unique_ptr<CurvilinearGridBlockUndoAction> undoAction = CurvilinearGridBlockUndoAction::Create(m_grid, m_lowerLeft, m_upperRight);

    bool awayFromUpperNBoundary = m_upperRight.m_n < m_grid.NumN() - 2;
    bool awayFromUpperMBoundary = m_upperRight.m_m < m_grid.NumM() - 2;

    bool awayFromLowerNBoundary = m_lowerLeft.m_n > 1;
    bool awayFromLowerMBoundary = m_lowerLeft.m_m > 1;

    bool onUpperNBoundary = m_upperRight.m_n == m_grid.NumN() - 1;
    bool onUpperMBoundary = m_upperRight.m_m == m_grid.NumM() - 1;

    bool onLowerNBoundary = m_lowerLeft.m_n == 0;
    bool onLowerMBoundary = m_lowerLeft.m_m == 0;

    for (UInt n = lowerLimitI + 1; n < upperLimitI; ++n)
    {
        for (UInt m = lowerLimitJ + 1; m < upperLimitJ; ++m)
        {
            m_grid.GetNode(n, m).SetInvalid();

            if (awayFromUpperNBoundary && n == upperLimitI - 1 && !m_grid.GetNode(n + 2, m).IsValid())
            {
                m_grid.GetNode(n + 1, m).SetInvalid();
            }

            if (awayFromLowerNBoundary && n == lowerLimitI + 1 && !m_grid.GetNode(n - 2, m).IsValid())
            {
                m_grid.GetNode(n - 1, m).SetInvalid();
            }

            if (onUpperNBoundary && n == upperLimitI - 1)
            {
                m_grid.GetNode(n + 1, m).SetInvalid();
            }

            if (onLowerNBoundary && n == 1)
            {
                m_grid.GetNode(n - 1, m).SetInvalid();
            }

            if (awayFromUpperMBoundary && m == upperLimitJ - 1 && !m_grid.GetNode(n, m + 2).IsValid())
            {
                m_grid.GetNode(n, m + 1).SetInvalid();
            }

            if (awayFromLowerMBoundary && m == lowerLimitJ + 1 && !m_grid.GetNode(n, m - 2).IsValid())
            {
                m_grid.GetNode(n, m - 1).SetInvalid();
            }

            if (onUpperMBoundary && m == upperLimitJ - 1)
            {
                m_grid.GetNode(n, upperLimitJ).SetInvalid();
            }

            if (onLowerMBoundary && m == 1)
            {
                m_grid.GetNode(n, m - 1).SetInvalid();
            }
        }
    }

    return undoAction;
}
