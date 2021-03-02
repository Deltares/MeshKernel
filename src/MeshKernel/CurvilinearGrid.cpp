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

meshkernel::CurvilinearGrid::CurvilinearGrid(size_t m, size_t n) : m_numM(m), m_numN(n)
{
    m_gridNodes.resize(m + 1, std::vector<Point>(n + 1, {doubleMissingValue, doubleMissingValue}));
}

meshkernel::CurvilinearGrid::CurvilinearGrid(std::vector<std::vector<Point>> grid, Projection projection) : m_projection(projection),
                                                                                                            m_gridNodes(std::move(grid))
{
    m_numM = m_gridNodes.size();
    m_numN = m_gridNodes[0].size();
}

std::tuple<std::vector<meshkernel::Point>, std::vector<meshkernel::Edge>> meshkernel::CurvilinearGrid::ConvertCurvilinearToNodesAndEdges()
{
    if (m_gridNodes.empty())
    {
        throw std::invalid_argument("Mesh2D::Mesh2D: The curvilinear grid is empty.");
    }

    std::vector<Point> nodes(m_gridNodes.size() * m_gridNodes[0].size());
    std::vector<Edge> edges(m_gridNodes.size() * (m_gridNodes[0].size() - 1) + (m_gridNodes.size() - 1) * m_gridNodes[0].size());
    std::vector<std::vector<size_t>> indices(m_gridNodes.size(), std::vector<size_t>(m_gridNodes[0].size(), sizetMissingValue));

    size_t ind = 0;
    for (auto m = 0; m < m_gridNodes.size(); m++)
    {
        for (auto n = 0; n < m_gridNodes[0].size(); n++)
        {
            if (m_gridNodes[m][n].IsValid())
            {
                nodes[ind] = m_gridNodes[m][n];
                indices[m][n] = ind;
                ind++;
            }
        }
    }
    nodes.resize(ind);

    ind = 0;
    for (auto m = 0; m < m_gridNodes.size() - 1; m++)
    {
        for (auto n = 0; n < m_gridNodes[0].size(); n++)
        {
            if (indices[m][n] != sizetMissingValue && indices[m + 1][n] != sizetMissingValue)
            {
                edges[ind].first = indices[m][n];
                edges[ind].second = indices[m + 1][n];
                ind++;
            }
        }
    }

    for (auto m = 0; m < m_gridNodes.size(); m++)
    {
        for (auto n = 0; n < m_gridNodes[0].size() - 1; n++)
        {
            if (indices[m][n] != sizetMissingValue && indices[m][n + 1] != sizetMissingValue)
            {
                edges[ind].first = indices[m][n];
                edges[ind].second = indices[m][n + 1];
                ind++;
            }
        }
    }
    edges.resize(ind);

    return {nodes, edges};
}

void meshkernel::CurvilinearGrid::SetFlatCopies()
{
    const auto [nodes, edges] = ConvertCurvilinearToNodesAndEdges();
    m_nodes = nodes;
    m_edges = edges;
    Mesh::SetFlatCopies();
}

void meshkernel::CurvilinearGrid::BuildTree()
{
    m_flattenNodes.resize(m_gridNodes.size() * m_gridNodes[0].size(), {doubleMissingValue, doubleMissingValue, sizetMissingValue, sizetMissingValue});

    size_t index = 0;
    for (auto n = 0; n < m_gridNodes[0].size(); ++n)
    {
        for (auto m = 0; m < m_gridNodes.size(); ++m)
        {
            m_flattenNodes[index].x = m_gridNodes[m][n].x;
            m_flattenNodes[index].y = m_gridNodes[m][n].y;
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
