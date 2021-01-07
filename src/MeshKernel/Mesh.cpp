#pragma once
#include "MeshKernel/Mesh.hpp"
#include <MeshKernel/Entities.hpp>
#include <vector>

meshkernel::Mesh::Mesh(const std::vector<Edge>& edges,
                       const std::vector<Point>& nodes,
                       Projection projection) : m_edges(edges), m_nodes(nodes), m_projection(projection){};

void meshkernel::Mesh::NodeAdministration()
{
    // assume no duplicated links
    for (auto e = 0; e < GetNumEdges(); e++)
    {
        const auto firstNode = m_edges[e].first;
        const auto secondNode = m_edges[e].second;

        if (firstNode == sizetMissingValue || secondNode == sizetMissingValue)
        {
            continue;
        }

        if (m_nodesNumEdges[firstNode] >= maximumNumberOfEdgesPerNode || m_nodesNumEdges[secondNode] >= maximumNumberOfEdgesPerNode)
        {
            continue;
        }

        // Search for previously connected edges
        auto alreadyAddedEdge = false;
        for (auto i = 0; i < m_nodesNumEdges[firstNode]; ++i)
        {
            const auto currentEdge = m_edges[m_nodesEdges[firstNode][i]];
            if (currentEdge.first == secondNode || currentEdge.second == secondNode)
            {
                alreadyAddedEdge = true;
                break;
            }
        }
        if (!alreadyAddedEdge)
        {
            m_nodesEdges[firstNode][m_nodesNumEdges[firstNode]] = e;
            m_nodesNumEdges[firstNode]++;
        }

        // Search for previously connected edges
        alreadyAddedEdge = false;
        for (auto i = 0; i < m_nodesNumEdges[secondNode]; ++i)
        {
            const auto currentEdge = m_edges[m_nodesEdges[secondNode][i]];
            if (currentEdge.first == firstNode || currentEdge.second == firstNode)
            {
                alreadyAddedEdge = true;
                break;
            }
        }
        if (!alreadyAddedEdge)
        {
            m_nodesEdges[secondNode][m_nodesNumEdges[secondNode]] = e;
            m_nodesNumEdges[secondNode]++;
        }
    }

    // resize
    for (auto n = 0; n < GetNumNodes(); n++)
    {
        m_nodesEdges[n].resize(m_nodesNumEdges[n]);
    }
}

void meshkernel::Mesh::DeleteInvalidNodesAndEdges()
{

    // Mask nodes connected to valid edges
    std::vector<bool> connectedNodes(m_nodes.size(), false);
    size_t numInvalidEdges = 0;

    for (const auto& edge : m_edges)
    {
        auto const firstNode = edge.first;
        auto const secondNode = edge.second;

        if (firstNode == sizetMissingValue || secondNode == sizetMissingValue)
        {
            numInvalidEdges++;
            continue;
        }

        connectedNodes[firstNode] = true;
        connectedNodes[secondNode] = true;
    }

    // Count all invalid nodes (note: there might be nodes that are not connected to an edge)
    size_t numInvalidNodes = 0;
    for (auto n = 0; n < m_nodes.size(); ++n)
    {
        // invalidate nodes that are not connected
        if (!connectedNodes[n])
        {
            m_nodes[n] = {doubleMissingValue, doubleMissingValue};
        }

        if (!m_nodes[n].IsValid())
        {
            numInvalidNodes++;
        }
    }

    // If nothing to invalidate return
    if (numInvalidEdges == 0 && numInvalidNodes == 0)
    {
        m_numNodes = m_nodes.size();
        m_numEdges = m_edges.size();
        return;
    }

    // Flag invalid nodes
    std::vector<size_t> validNodesIndices(m_nodes.size());
    std::fill(validNodesIndices.begin(), validNodesIndices.end(), sizetMissingValue);
    size_t validIndex = 0;
    for (auto n = 0; n < m_nodes.size(); ++n)
    {
        if (m_nodes[n].IsValid())
        {
            validNodesIndices[n] = validIndex;
            validIndex++;
        }
    }

    // Flag invalid edges
    for (auto& edge : m_edges)
    {
        auto const firstNode = edge.first;
        auto const secondNode = edge.second;

        if (firstNode != sizetMissingValue && secondNode != sizetMissingValue && validNodesIndices[firstNode] != sizetMissingValue && validNodesIndices[secondNode] != sizetMissingValue)
        {
            edge.first = validNodesIndices[firstNode];
            edge.second = validNodesIndices[secondNode];
            continue;
        }

        edge.first = sizetMissingValue;
        edge.second = sizetMissingValue;
    }

    // Remove invalid nodes, without reducing capacity
    const auto endNodeVector = std::remove_if(m_nodes.begin(), m_nodes.end(), [](const Point& n) { return !n.IsValid(); });
    m_nodes.erase(endNodeVector, m_nodes.end());
    m_numNodes = m_nodes.size();

    // Remove invalid edges, without reducing capacity
    const auto endEdgeVector = std::remove_if(m_edges.begin(), m_edges.end(), [](const Edge& e) { return e.first == sizetMissingValue || e.second == sizetMissingValue; });
    m_edges.erase(endEdgeVector, m_edges.end());
    m_numEdges = m_edges.size();
}