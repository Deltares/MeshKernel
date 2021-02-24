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

#include <vector>

#include <MeshKernel/CurvilinearGrid.hpp>

meshkernel::CurvilinearGrid::CurvilinearGrid(size_t m, size_t n)
{
    m_nodes.resize(m + 1, std::vector<Point>(n + 1, {doubleMissingValue, doubleMissingValue}));
    m_numM = m;
    m_numN = n;
}

meshkernel::CurvilinearGrid::CurvilinearGrid(const std::vector<std::vector<Point>>& grid, Projection projection) : m_projection(projection)
{
    m_numM = grid.size();
    m_numN = grid[0].size();
    m_nodes = grid;
}

void meshkernel::CurvilinearGrid::BuildTree()
{
    m_flattenNodes.resize(m_nodes.size() * m_nodes[0].size(), {doubleMissingValue, doubleMissingValue, sizetMissingValue, sizetMissingValue});

    size_t index = 0;
    for (auto n = 0; n < m_nodes[0].size(); ++n)
    {
        for (auto m = 0; m < m_nodes.size(); ++m)
        {
            m_flattenNodes[index].x = m_nodes[m][n].x;
            m_flattenNodes[index].y = m_nodes[m][n].y;
            m_flattenNodes[index].m = m;
            m_flattenNodes[index].n = n;
            index++;
        }
    }

    m_nodesRTree.BuildTree(m_flattenNodes);
}

std::tuple<int, int> meshkernel::CurvilinearGrid::GetNodeIndices(Point point)
{
    m_nodesRTree.NearestNeighbors(point);
    if (m_nodesRTree.GetQueryResultSize() <= 0)
    {
        return {sizetMissingValue, sizetMissingValue};
    }
    const auto nodeIndex = m_nodesRTree.GetQueryResult(0);
    return {m_flattenNodes[nodeIndex].m, m_flattenNodes[nodeIndex].n};
}
