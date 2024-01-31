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

#include "MeshKernel/CurvilinearGrid/CurvilinearGridDeleteExterior.hpp"

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLine.hpp>

using meshkernel::CurvilinearGrid;
using meshkernel::CurvilinearGridLine;
using meshkernel::CurvilinearGridNodeIndices;

meshkernel::CurvilinearGridDeleteExterior::CurvilinearGridDeleteExterior(CurvilinearGrid& grid)
    : CurvilinearGridAlgorithm(grid)
{
}

void meshkernel::CurvilinearGridDeleteExterior::Compute()
{
    const UInt lowerLimitI = m_lowerLeft.m_n;
    const UInt upperLimitI = m_upperRight.m_n;

    const UInt lowerLimitJ = m_lowerLeft.m_m;
    const UInt upperLimitJ = m_upperRight.m_m;

    // Split into 4 regions, setting the nodes in each region to invalid
    //
    // First region: all nodes "south" the designated box
    for (UInt n = 0; n < m_grid.NumN(); ++n)
    {
        for (UInt m = 0; m < lowerLimitJ; ++m)
        {
            m_grid.m_gridNodes(n, m).SetInvalid();
        }
    }

    // Second region: all nodes "directly west of" the designated box
    for (UInt n = 0; n < lowerLimitI; ++n)
    {
        for (UInt m = lowerLimitJ; m <= upperLimitJ; ++m)
        {
            m_grid.m_gridNodes(n, m).SetInvalid();
        }
    }

    // Third region: all nodes "directly east of" the designated box
    for (UInt n = upperLimitI + 1; n < m_grid.NumN(); ++n)
    {
        for (UInt m = lowerLimitJ; m <= upperLimitJ; ++m)
        {
            m_grid.m_gridNodes(n, m).SetInvalid();
        }
    }

    // Fourth region: all nodes "north" the designated box
    for (UInt n = 0; n < m_grid.NumN(); ++n)
    {
        for (UInt m = upperLimitJ + 1; m < m_grid.NumM(); ++m)
        {
            m_grid.m_gridNodes(n, m).SetInvalid();
        }
    }
}
