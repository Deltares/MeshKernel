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

meshkernel::CurvilinearGrid::CurvilinearGrid(std::vector<std::vector<Point>> grid, Projection projection) : m_projection(projection),
                                                                                                            m_gridNodes(std::move(grid))
{
    m_numM = m_gridNodes.size();
    m_numN = m_gridNodes[0].size();
}

std::tuple<std::vector<meshkernel::Point>, std::vector<meshkernel::Edge>, std::vector<std::pair<size_t, size_t>>> meshkernel::CurvilinearGrid::ConvertCurvilinearToNodesAndEdges()
{
    if (m_gridNodes.empty())
    {
        throw std::invalid_argument("Mesh2D::Mesh2D: The curvilinear grid is empty.");
    }

    std::vector<Point> nodes(m_gridNodes.size() * m_gridNodes[0].size());
    std::vector<Edge> edges(m_gridNodes.size() * (m_gridNodes[0].size() - 1) + (m_gridNodes.size() - 1) * m_gridNodes[0].size());
    std::vector<std::vector<size_t>> nodeIndices(m_gridNodes.size(), std::vector<size_t>(m_gridNodes[0].size(), sizetMissingValue));
    std::vector<std::pair<size_t, size_t>> gridIndices(nodes.size(), std::pair<size_t, size_t>(sizetMissingValue, sizetMissingValue));

    size_t ind = 0;
    for (auto m = 0; m < m_gridNodes.size(); m++)
    {
        for (auto n = 0; n < m_gridNodes[0].size(); n++)
        {
            if (m_gridNodes[m][n].IsValid())
            {
                nodes[ind] = m_gridNodes[m][n];
                nodeIndices[m][n] = ind;
                gridIndices[ind] = {m, n};
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
            if (nodeIndices[m][n] != sizetMissingValue && nodeIndices[m + 1][n] != sizetMissingValue)
            {
                edges[ind].first = nodeIndices[m][n];
                edges[ind].second = nodeIndices[m + 1][n];
                ind++;
            }
        }
    }

    for (auto m = 0; m < m_gridNodes.size(); m++)
    {
        for (auto n = 0; n < m_gridNodes[0].size() - 1; n++)
        {
            if (nodeIndices[m][n] != sizetMissingValue && nodeIndices[m][n + 1] != sizetMissingValue)
            {
                edges[ind].first = nodeIndices[m][n];
                edges[ind].second = nodeIndices[m][n + 1];
                ind++;
            }
        }
    }
    edges.resize(ind);

    return {nodes, edges, gridIndices};
}

void meshkernel::CurvilinearGrid::SetFlatCopies()
{
    const auto [nodes, edges, gridIndices] = ConvertCurvilinearToNodesAndEdges();
    m_nodes = nodes;
    m_edges = edges;
    m_gridIndices = gridIndices;
    Mesh::SetFlatCopies();
}

void meshkernel::CurvilinearGrid::BuildTree()
{
    SetFlatCopies();
    m_nodesRTree.BuildTree(m_nodes);
}

std::tuple<int, int> meshkernel::CurvilinearGrid::GetNodeIndices(Point point)
{
    SearchNearestNeighbors(point, MeshLocations::Nodes);
    if (GetNumNearestNeighbors(MeshLocations::Nodes) == 0)
    {
        return {sizetMissingValue, sizetMissingValue};
    }

    const auto nodeIndex = GetNearestNeighborIndex(0, MeshLocations::Nodes);
    return {m_gridIndices[nodeIndex].first, m_gridIndices[nodeIndex].second};
}
