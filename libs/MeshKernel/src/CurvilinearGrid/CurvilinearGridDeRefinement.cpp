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
#include <MeshKernel/Entities.hpp>

using meshkernel::CurvilinearGrid;
using meshkernel::CurvilinearGridDeRefinement;

CurvilinearGridDeRefinement::CurvilinearGridDeRefinement(CurvilinearGrid& grid) : CurvilinearGridAlgorithm(grid)
{
}

void CurvilinearGridDeRefinement::Compute()
{
    if (!m_lowerLeft.IsValid() || !m_upperRight.IsValid())
    {
        throw std::invalid_argument("CurvilinearGridDeRefinement::Compute: lower left and upper right corners defining the curvilinear grid block are not set");
    }

    /// estimate the dimension of the refined grid
    const auto numMToDeRefine = m_upperRight.m_m > m_lowerLeft.m_m ? m_upperRight.m_m - m_lowerLeft.m_m : 1;
    const auto numNToDeRefine = m_upperRight.m_n > m_lowerLeft.m_n ? m_upperRight.m_n - m_lowerLeft.m_n : 1;

    // the de-refined grid
    std::vector<std::vector<Point>> deRefinedGrid;
    deRefinedGrid.reserve(m_grid.NumM());

    UInt mIndexOriginalGrid = 0;
    while (mIndexOriginalGrid < m_grid.NumM())
    {
        UInt localMDeRefinement = 1;
        if (mIndexOriginalGrid >= m_lowerLeft.m_m && mIndexOriginalGrid < m_upperRight.m_m)
        {
            localMDeRefinement = numMToDeRefine;
        }
        deRefinedGrid.emplace_back(std::vector<Point>());
        deRefinedGrid.back().reserve(m_grid.NumN());

        UInt nIndexOriginalGrid = 0;
        while (nIndexOriginalGrid < m_grid.NumN())
        {
            UInt localNDeRefinement = 1;
            if (nIndexOriginalGrid >= m_lowerLeft.m_n && nIndexOriginalGrid < m_upperRight.m_n)
            {
                localNDeRefinement = numNToDeRefine;
            }
            deRefinedGrid.back().emplace_back(m_grid.m_gridNodes(mIndexOriginalGrid, nIndexOriginalGrid));
            nIndexOriginalGrid += localNDeRefinement;
        }
        mIndexOriginalGrid += localMDeRefinement;
    }

    // substitute original grid with the derefined one
    m_grid.SetGridNodes(lin_alg::STLVectorOfVectorsToMatrix(deRefinedGrid));
}
