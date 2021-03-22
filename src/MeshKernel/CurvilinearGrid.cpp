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

#include <stdexcept>
#include <vector>

#include <MeshKernel/CurvilinearGrid.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Splines.hpp>

meshkernel::CurvilinearGrid::CurvilinearGrid(std::vector<std::vector<Point>>&& grid, Projection projection) : m_gridNodes(std::move(grid))
{
    if (m_gridNodes.empty())
    {
        throw std::invalid_argument("CurvilinearGrid::CurvilinearGrid: m_gridNodes is empty. ");
    }
    if (m_gridNodes[0].empty())
    {
        throw std::invalid_argument("CurvilinearGrid::CurvilinearGrid: m_gridNodes[0] is empty");
    }
    if (m_gridNodes.size() < 2)
    {
        throw std::invalid_argument("CurvilinearGrid::CurvilinearGrid: m_gridNodes.size() < 2");
    }
    if (m_gridNodes[0].size() < 2)
    {
        throw std::invalid_argument("CurvilinearGrid::CurvilinearGrid: m_gridNodes[0].size() < 2");
    }

    m_projection = projection;
    m_numM = m_gridNodes.size();
    m_numN = m_gridNodes[0].size();

    SetFlatCopies();
}

void meshkernel::CurvilinearGrid::SetFlatCopies()
{
    const auto [nodes, edges, gridIndices] = ConvertCurvilinearToNodesAndEdges();
    m_nodes = nodes;
    m_edges = edges;
    m_gridIndices = gridIndices;
}

std::tuple<std::vector<meshkernel::Point>, std::vector<meshkernel::Edge>, std::vector<meshkernel::CurvilinearGrid::NodeIndices>> meshkernel::CurvilinearGrid::ConvertCurvilinearToNodesAndEdges()
{
    if (m_gridNodes.empty())
    {
        throw std::invalid_argument("CurvilinearGrid::ConvertCurvilinearToNodesAndEdges: m_gridNodes is empty ");
    }
    if (m_gridNodes[0].empty())
    {
        throw std::invalid_argument("CurvilinearGrid::ConvertCurvilinearToNodesAndEdges: m_gridNodes[0] is empty");
    }
    if (m_gridNodes.size() < 2)
    {
        throw std::invalid_argument("CurvilinearGrid::ConvertCurvilinearToNodesAndEdges: m_gridNodes.size() < 2");
    }
    if (m_gridNodes[0].size() < 2)
    {
        throw std::invalid_argument("CurvilinearGrid::ConvertCurvilinearToNodesAndEdges: m_gridNodes[0].size() < 2");
    }

    std::vector<Point> nodes(m_gridNodes.size() * m_gridNodes[0].size());
    std::vector<Edge> edges(m_gridNodes.size() * (m_gridNodes[0].size() - 1) + (m_gridNodes.size() - 1) * m_gridNodes[0].size());
    std::vector<std::vector<size_t>> nodeIndices(m_gridNodes.size(), std::vector<size_t>(m_gridNodes[0].size(), sizetMissingValue));
    std::vector<meshkernel::CurvilinearGrid::NodeIndices> gridIndices(nodes.size(), meshkernel::CurvilinearGrid::NodeIndices{sizetMissingValue, sizetMissingValue});

    size_t ind = 0;
    for (size_t m = 0; m < m_gridNodes.size(); m++)
    {
        for (size_t n = 0; n < m_gridNodes[0].size(); n++)
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

meshkernel::CurvilinearGrid::NodeIndices meshkernel::CurvilinearGrid::GetNodeIndices(Point point)
{
    SearchNearestNeighbors(point, MeshLocations::Nodes);
    if (GetNumNearestNeighbors(MeshLocations::Nodes) == 0)
    {
        return {sizetMissingValue, sizetMissingValue};
    }

    const auto nodeIndex = GetNearestNeighborIndex(0, MeshLocations::Nodes);
    return m_gridIndices[nodeIndex];
}

bool meshkernel::CurvilinearGrid::IsValidFace(size_t m, size_t n) const
{
    return m_gridNodes[m][n].IsValid() &&
           m_gridNodes[m + 1][n].IsValid() &&
           m_gridNodes[m][n + 1].IsValid() &&
           m_gridNodes[m + 1][n + 1].IsValid();
};

void meshkernel::CurvilinearGrid::ComputeGridFacesMask()
{
    // Flag valid faces
    m_gridFacesMask.resize(m_numM - 1, std::vector<bool>(m_numN - 1, false));
    for (auto m = 0; m < m_numM - 1; ++m)
    {
        for (auto n = 0; n < m_numN - 1; ++n)
        {
            // Only if all grid nodes of the face are valid, the face is valid
            if (!IsValidFace(m, n))
            {
                continue;
            }
            m_gridFacesMask[m][n] = true;
        }
    }
}

void meshkernel::CurvilinearGrid::RemoveInvalidNodes(bool invalidNodesToRemove)
{

    if (!invalidNodesToRemove)
    {
        return;
    }

    // Compute the face mask
    ComputeGridFacesMask();

    invalidNodesToRemove = false;
    // Flag nodes not connected to valid faces
    for (auto m = 1; m < m_numM - 1; ++m)
    {
        for (auto n = 1; n < m_numN - 1; ++n)
        {
            if (m_gridNodes[m][n].IsValid() &&
                !m_gridFacesMask[m][n] &&
                !m_gridFacesMask[m - 1][n] &&
                !m_gridFacesMask[m - 1][n - 1] &&
                !m_gridFacesMask[m][n - 1])
            {
                m_gridNodes[m][n] = {doubleMissingValue, doubleMissingValue};
                invalidNodesToRemove = true;
            }
        }
    }

    for (auto m = 1; m < m_numM - 1; ++m)
    {
        if (m_gridNodes[m][0].IsValid() &&
            !m_gridFacesMask[m - 1][0] &&
            !m_gridFacesMask[m][0])
        {
            m_gridNodes[m][0] = {doubleMissingValue, doubleMissingValue};
        }
    }

    for (auto n = 1; n < m_numN - 1; ++n)
    {
        if (m_gridNodes[0][n].IsValid() &&
            !m_gridFacesMask[0][n - 1] &&
            !m_gridFacesMask[0][n])
        {
            m_gridNodes[0][n] = {doubleMissingValue, doubleMissingValue};
        }
    }

    if (m_gridNodes[0][0].IsValid() && !m_gridFacesMask[0][0])
    {
        m_gridNodes[0][0] = {doubleMissingValue, doubleMissingValue};
    }

    RemoveInvalidNodes(invalidNodesToRemove);
}

void meshkernel::CurvilinearGrid::ComputeGridNodeTypes()
{
    RemoveInvalidNodes(true);
    m_gridNodesMask.resize(m_numM, std::vector<NodeType>(m_numN, NodeType::Invalid));

    // Flag faces based on boundaries
    for (size_t m = 0; m < m_numM; ++m)
    {
        for (size_t n = 0; n < m_numN; ++n)
        {

            if (!m_gridNodes[m][n].IsValid())
            {
                continue;
            }

            // Left side
            if (m == 0 && n == 0)
            {
                m_gridNodesMask[m][n] = NodeType::BottomLeft;
                continue;
            }
            if (m == 0 && n == m_numN - 1)
            {
                m_gridNodesMask[m][n] = NodeType::UpperLeft;
                continue;
            }
            if (m == 0 && !m_gridNodes[m][n - 1].IsValid())
            {
                m_gridNodesMask[m][n] = NodeType::BottomLeft;
                continue;
            }
            if (m == 0 && !m_gridNodes[m][n + 1].IsValid())
            {
                m_gridNodesMask[m][n] = NodeType::UpperLeft;
                continue;
            }
            if (m == 0)
            {
                m_gridNodesMask[m][n] = NodeType::Left;
                continue;
            }
            // Right side
            if (m == m_numM - 1 && n == 0)
            {
                m_gridNodesMask[m][n] = NodeType::BottomRight;
                continue;
            }
            if (m == m_numM - 1 && n == m_numN - 1)
            {
                m_gridNodesMask[m][n] = NodeType::UpperRight;
                continue;
            }
            if (m == m_numM - 1 && !m_gridNodes[m][n - 1].IsValid())
            {
                m_gridNodesMask[m][n] = NodeType::BottomRight;
                continue;
            }
            if (m == m_numM - 1 && !m_gridNodes[m][n + 1].IsValid())
            {
                m_gridNodesMask[m][n] = NodeType::UpperRight;
                continue;
            }
            if (m == m_numM - 1)
            {
                m_gridNodesMask[m][n] = NodeType::Right;
                continue;
            }
            // Bottom side
            if (n == 0 && !m_gridNodes[m - 1][n].IsValid())
            {
                m_gridNodesMask[m][n] = NodeType::BottomLeft;
                continue;
            }
            if (n == 0 && !m_gridNodes[m + 1][n].IsValid())
            {
                m_gridNodesMask[m][n] = NodeType::BottomRight;
                continue;
            }
            if (n == 0)
            {
                m_gridNodesMask[m][n] = NodeType::Bottom;
                continue;
            }
            // Upper side
            if (n == m_numN - 1 && !m_gridNodes[m - 1][n].IsValid())
            {
                m_gridNodesMask[m][n] = NodeType::UpperLeft;
                continue;
            }
            if (n == m_numN - 1 && !m_gridNodes[m + 1][n].IsValid())
            {
                m_gridNodesMask[m][n] = NodeType::UpperRight;
                continue;
            }
            if (n == m_numN - 1)
            {
                m_gridNodesMask[m][n] = NodeType::Up;
                continue;
            }

            const auto isTopLeftFaceValid = m_gridFacesMask[m - 1][n];
            const auto isTopRightFaceValid = m_gridFacesMask[m][n];
            const auto isBottomLeftFaceValid = m_gridFacesMask[m - 1][n - 1];
            const auto isBottomRightFaceValid = m_gridFacesMask[m][n - 1];

            if (isTopRightFaceValid &&
                isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                m_gridNodesMask[m][n] = NodeType::InternalValid;
                continue;
            }
            if (!isTopRightFaceValid &&
                isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                m_gridNodesMask[m][n] = NodeType::BottomLeft;
                continue;
            }
            if (isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                m_gridNodesMask[m][n] = NodeType::BottomRight;
                continue;
            }
            if (isTopRightFaceValid &&
                isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                m_gridNodesMask[m][n] = NodeType::UpperRight;
                continue;
            }
            if (isTopRightFaceValid &&
                isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                m_gridNodesMask[m][n] = NodeType::UpperLeft;
                continue;
            }

            if (isTopRightFaceValid &&
                isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                m_gridNodesMask[m][n] = NodeType::Bottom;
                continue;
            }
            if (isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                m_gridNodesMask[m][n] = NodeType::Left;
                continue;
            }

            if (!isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                m_gridNodesMask[m][n] = NodeType::Up;
                continue;
            }

            if (!isTopRightFaceValid &&
                isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                m_gridNodesMask[m][n] = NodeType::Right;
                continue;
            }

            if (isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                m_gridNodesMask[m][n] = NodeType::BottomLeft;
                continue;
            }

            if (!isTopRightFaceValid &&
                isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                m_gridNodesMask[m][n] = NodeType::BottomRight;
                continue;
            }

            if (!isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                isBottomLeftFaceValid &&
                !isBottomRightFaceValid)
            {
                m_gridNodesMask[m][n] = NodeType::UpperRight;
                continue;
            }

            if (!isTopRightFaceValid &&
                !isTopLeftFaceValid &&
                !isBottomLeftFaceValid &&
                isBottomRightFaceValid)
            {
                m_gridNodesMask[m][n] = NodeType::UpperLeft;
            }
        }
    }
}

void meshkernel::CurvilinearGrid::InsertFace(Point firstPoint, Point secondPoint)
{
    // Compute the grid node types
    ComputeGridNodeTypes();

    // Get the indices of the two closest nodes to `point`
    auto [mFirstNode, nFirstNode] = GetNodeIndices(firstPoint);
    auto [mSecondNode, nSecondNode] = GetNodeIndices(secondPoint);
}

bool meshkernel::CurvilinearGrid::AreNeighbors(meshkernel::CurvilinearGrid::NodeIndices firstNode,
                                               meshkernel::CurvilinearGrid::NodeIndices secondNode) const
{
    if (firstNode.m == secondNode.m)
    {
        if (AbsoluteDifference(firstNode.n, secondNode.n) == 1)
        {
            return true;
        }
    }
    if (firstNode.n == secondNode.n)
    {
        if (AbsoluteDifference(firstNode.m, secondNode.m) == 1)
        {
            return true;
        }
    }
    return false;
}
