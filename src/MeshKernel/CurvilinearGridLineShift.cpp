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

#include <MeshKernel/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGridLineShift.hpp>
#include <MeshKernel/Entities.hpp>

using meshkernel::CurvilinearGrid;
using meshkernel::CurvilinearGridLineShift;

CurvilinearGridLineShift::CurvilinearGridLineShift(std::shared_ptr<CurvilinearGrid> grid) : CurvilinearGridAlgorithm(grid)

{
    // Allocate cache for storing grid nodes values
    m_gridNodesCache.resize(m_grid->m_gridNodes.size(), std::vector<Point>(m_grid->m_gridNodes[0].size()));
    // Compute the grid node types
    m_grid->ComputeGridNodeTypes();
}

std::shared_ptr<CurvilinearGrid> CurvilinearGridLineShift::Compute()
{
    return nullptr;
}

void CurvilinearGridLineShift::MoveNode(Point const& fromPoint, Point const& toPoint)
{
}

void CurvilinearGridLineShift::SnapGrid()
{
}