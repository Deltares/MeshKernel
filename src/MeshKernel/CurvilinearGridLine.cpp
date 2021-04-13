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

#include <MeshKernel/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGridLine.hpp>

using meshkernel::CurvilinearGrid;
using meshkernel::CurvilinearGridLine;

CurvilinearGridLine::CurvilinearGridLine(CurvilinearGrid::NodeIndices const& startNode, CurvilinearGrid::NodeIndices const& endNode)
    : m_startNode(startNode),
      m_endNode(endNode)
{
    if (m_startNode == m_endNode)
    {
        throw std::invalid_argument("CurvilinearGridLine::CurvilinearGridLine Cannot construct a grid line with coinciding nodes.");
    }

    m_gridLineType = m_startNode.m_m == m_endNode.m_m ? GridLineType::NGridLine : GridLineType::MGridLine;
    m_startCoordinate = m_gridLineType == GridLineType::NGridLine ? m_startNode.m_n : m_startNode.m_m;
    m_endCoordinate = m_gridLineType == GridLineType::NGridLine ? m_endNode.m_n : m_endNode.m_m;
    m_constantCoordinate = m_gridLineType == GridLineType::NGridLine ? m_startNode.m_m : m_startNode.m_n;
}

bool CurvilinearGridLine::IsNodeOnLine(CurvilinearGrid::NodeIndices const& node) const
{
    for (auto i = m_startCoordinate; i < m_endCoordinate; ++i)
    {
        if (m_gridLineType == GridLineType::MGridLine && node.m_m == i && node.m_n == m_constantCoordinate)
        {
            return true;
        }
        if (m_gridLineType == GridLineType::NGridLine && node.m_n == i && node.m_m == m_constantCoordinate)
        {
            return true;
        }
    }
    return false;
}

CurvilinearGrid::NodeIndices CurvilinearGridLine::GetNodeindexFromCoordinate(size_t const& coordinate) const
{
    auto const mCoordinate = m_gridLineType == GridLineType::MGridLine ? coordinate : m_constantCoordinate;
    auto const nCoordinate = m_gridLineType == GridLineType::MGridLine ? m_constantCoordinate : coordinate;

    return {mCoordinate, nCoordinate};
}
