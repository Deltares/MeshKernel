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
#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLine.hpp>

using meshkernel::CurvilinearGrid;
using meshkernel::CurvilinearGridLine;
using meshkernel::CurvilinearGridNodeIndices;

meshkernel::CurvilinearGridDeleteInterior::CurvilinearGridDeleteInterior(CurvilinearGrid& grid)
    : CurvilinearGridAlgorithm(grid)
{
}

void meshkernel::CurvilinearGridDeleteInterior::Compute()
{

    if (m_lowerLeft.m_m >= m_grid.m_numM || m_lowerLeft.m_n >= m_grid.m_numN)
    {
        throw ConstraintError("Invalid index: first index {{{}, {}}} not in mesh limits {{{}, {}}}", m_lowerLeft.m_m, m_lowerLeft.m_n, m_grid.m_numM, m_grid.m_numN);
    }

    if (m_upperRight.m_m >= m_grid.m_numM || m_upperRight.m_n >= m_grid.m_numN)
    {
        throw ConstraintError("Invalid index: second index {{{}, {}}} not in mesh limits {{{}, {}}}", m_upperRight.m_m, m_upperRight.m_n, m_grid.m_numM, m_grid.m_numN);
    }

    UInt lowerLimitI = m_lowerLeft.m_n;
    UInt upperLimitI = m_upperRight.m_n;

    UInt lowerLimitJ = m_lowerLeft.m_m;
    UInt upperLimitJ = m_upperRight.m_m;

    for (UInt n = lowerLimitI + 1; n < upperLimitI; ++n)
    {
        for (UInt m = lowerLimitJ + 1; m < upperLimitJ; ++m)
        {
            m_grid.m_gridNodes(n, m).SetInvalid();
        }
    }
}
