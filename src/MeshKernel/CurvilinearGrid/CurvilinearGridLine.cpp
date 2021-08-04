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
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLine.hpp>

using meshkernel::CurvilinearGrid;
using meshkernel::CurvilinearGridLine;
using meshkernel::CurvilinearGridNodeIndices;

CurvilinearGridLine::CurvilinearGridLine(CurvilinearGridNodeIndices const& startNode, CurvilinearGridNodeIndices const& endNode)
    : m_startNode(startNode),
      m_endNode(endNode)
{
    if (m_startNode == m_endNode)
    {
        throw std::invalid_argument("CurvilinearGridLine::CurvilinearGridLine Cannot construct a grid line with coinciding nodes.");
    }

    m_gridLineType = m_startNode.m_m == m_endNode.m_m ? GridLineDirection::NDirection : GridLineDirection::MDirection;
    m_startCoordinate = IsNGridLine() ? m_startNode.m_n : m_startNode.m_m;
    m_endCoordinate = IsNGridLine() ? m_endNode.m_n : m_endNode.m_m;
    m_constantCoordinate = IsNGridLine() ? m_startNode.m_m : m_startNode.m_n;
}

bool CurvilinearGridLine::IsNodeOnLine(CurvilinearGridNodeIndices const& node) const
{
    for (auto i = m_startCoordinate; i <= m_endCoordinate; ++i)
    {
        if (IsMGridLine() && node.m_m == i && node.m_n == m_constantCoordinate)
        {
            return true;
        }
        if (IsNGridLine() && node.m_n == i && node.m_m == m_constantCoordinate)
        {
            return true;
        }
    }
    return false;
}

CurvilinearGridNodeIndices CurvilinearGridLine::GetNodeIndexFromCoordinate(size_t const& coordinate) const
{
    auto const mCoordinate = IsMGridLine() ? coordinate : m_constantCoordinate;
    auto const nCoordinate = IsMGridLine() ? m_constantCoordinate : coordinate;

    return {mCoordinate, nCoordinate};
}