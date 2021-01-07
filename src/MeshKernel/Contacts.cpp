#pragma once

#include <memory>
#include <vector>

#include "MeshKernel/Exceptions.hpp"
#include <MeshKernel/Contacts.hpp>
#include <MeshKernel/Operations.hpp>

meshkernel::Contacts::Contacts(std::shared_ptr<Mesh1D> mesh1d,
                               std::shared_ptr<Mesh2D> mesh) : m_mesh1d(mesh1d), m_mesh2d(mesh)
{
    // assert mesh1d and mesh have the same projection!
    if (m_mesh1d->m_projection != m_mesh2d->m_projection)
    {
        throw AlgorithmError("meshkernel::Contacts::Contacts: m_mesh1d and m_mesh2d projections are different");
    }
}

void meshkernel::Contacts::ComputeSingleConnections(const Polygons& polygons)
{
    m_mesh2d->FindFaces();
    m_mesh1d->FindFaces();

    const auto pointFaceIndices = m_mesh2d->PointFaceIndices(m_mesh1d->m_nodes);
    m_mesh1dIndices.reserve(m_mesh1d->m_nodes.size());
    m_mesh2dIndices.reserve(m_mesh1d->m_nodes.size());
    const double distanceFactor = 5.0;
    for (size_t n = 0; n < m_mesh1d->m_nodes.size(); ++n)
    {
        // Do not connect nodes at boundary of the network
        if (!m_mesh1d->m_nodeMask[n] || m_mesh1d->IsNodeOnBoundary(n))
        {
            continue;
        }

        // Inquire if a node is inside a the polygons
        if (!polygons.IsPointInPolygons(m_mesh1d->m_nodes[n]))
        {
            continue;
        }

        if (pointFaceIndices[n] != sizetMissingValue)
        {
            m_mesh1dIndices.emplace_back(n);
            m_mesh2dIndices.emplace_back(pointFaceIndices[n]);
            continue;
        }

        const auto left1dEdge = m_mesh1d->m_nodesEdges[n][0];
        const auto right1dEdge = m_mesh1d->m_nodesEdges[n][1];

        const auto otherLeft1dNode = m_mesh1d->m_edges[left1dEdge].first == n ? m_mesh1d->m_edges[left1dEdge].second : m_mesh1d->m_edges[left1dEdge].first;
        const auto otherRight1dNode = m_mesh1d->m_edges[right1dEdge].first == n ? m_mesh1d->m_edges[right1dEdge].second : m_mesh1d->m_edges[right1dEdge].first;

        const auto normalVector = NormalVectorOutside(m_mesh1d->m_nodes[otherLeft1dNode], m_mesh1d->m_nodes[otherRight1dNode], m_mesh1d->m_projection);
        const auto edgeLength = ComputeDistance(m_mesh1d->m_nodes[otherLeft1dNode], m_mesh1d->m_nodes[otherRight1dNode], m_mesh1d->m_projection);

        const auto rightProjectedNode = m_mesh1d->m_nodes[n] + normalVector * edgeLength * distanceFactor;
        size_t rightIntersectedFace;
        size_t rightIntersectedEdge;
        const auto isConnectionIntersectingAFace = m_mesh2d->IsSegmentCrossingAFace(m_mesh1d->m_nodes[n], rightProjectedNode, rightIntersectedFace, rightIntersectedEdge);
        if (isConnectionIntersectingAFace && !IsConnectionIntersectingMesh1d(n, rightIntersectedFace) && !IsContactIntersectingContact(n, rightIntersectedFace))
        {
            m_mesh1dIndices.emplace_back(n);
            m_mesh2dIndices.emplace_back(rightIntersectedFace);
        }

        const auto leftProjectedNode = m_mesh1d->m_nodes[n] - normalVector * edgeLength * distanceFactor;
        size_t leftIntersectedFace;
        size_t leftIntersectedEdge;
        const auto isLeftProjectedNodeIntersected = m_mesh2d->IsSegmentCrossingAFace(m_mesh1d->m_nodes[n], leftProjectedNode, leftIntersectedFace, leftIntersectedEdge);
        if (isLeftProjectedNodeIntersected && !IsConnectionIntersectingMesh1d(n, rightIntersectedFace) && !IsContactIntersectingContact(n, rightIntersectedFace))
        {
            m_mesh1dIndices.emplace_back(n);
            m_mesh2dIndices.emplace_back(leftIntersectedFace);
        }
    }
};

bool meshkernel::Contacts::IsConnectionIntersectingMesh1d(size_t node, size_t face) const
{
    for (size_t i = 0; i < m_mesh1dIndices.size(); ++i)
    {
        Point intersectionPoint;
        double crossProduct;
        double ratioFirstSegment;
        double ratioSecondSegment;
        const auto areSegmentCrossing = AreSegmentsCrossing(m_mesh1d->m_nodes[node],
                                                            m_mesh2d->m_facesCircumcenters[face],
                                                            m_mesh1d->m_nodes[m_mesh1dIndices[i]],
                                                            m_mesh2d->m_facesCircumcenters[m_mesh2dIndices[i]],
                                                            false,
                                                            m_mesh1d->m_projection,
                                                            intersectionPoint,
                                                            crossProduct,
                                                            ratioFirstSegment,
                                                            ratioSecondSegment);

        if (areSegmentCrossing)
        {
            return true;
        }
    }
    return false;
}

bool meshkernel::Contacts::IsContactIntersectingContact(size_t node, size_t face) const
{

    for (size_t e = 0; e < m_mesh1d->GetNumEdges(); ++e)
    {

        Point intersectionPoint;
        double crossProduct;
        double ratioFirstSegment;
        double ratioSecondSegment;
        const auto areSegmentCrossing = AreSegmentsCrossing(m_mesh1d->m_nodes[node],
                                                            m_mesh2d->m_facesCircumcenters[face],
                                                            m_mesh1d->m_nodes[m_mesh1d->m_edges[e].first],
                                                            m_mesh1d->m_nodes[m_mesh1d->m_edges[e].second],
                                                            false,
                                                            m_mesh1d->m_projection,
                                                            intersectionPoint,
                                                            crossProduct,
                                                            ratioFirstSegment,
                                                            ratioSecondSegment);

        if (areSegmentCrossing)
        {
            return true;
        }
    }
    return false;
}

void meshkernel::Contacts::ComputeMultipleConnections(){};

void meshkernel::Contacts::ComputeConnectionsWithPolygons(const Polygons& polygons){};

void meshkernel::Contacts::ComputeConnectionsWithPoints(const std::vector<Point>& points){};

void meshkernel::Contacts::ComputeBoundaryConnections(){};
