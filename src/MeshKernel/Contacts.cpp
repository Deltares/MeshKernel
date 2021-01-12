#pragma once

#include <iostream>
#include <memory>
#include <vector>

#include <MeshKernel/Contacts.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Polygons.hpp>

meshkernel::Contacts::Contacts(std::shared_ptr<Mesh1D> mesh1d,
                               std::shared_ptr<Mesh2D> mesh2d,
                               const std::vector<bool>& oneDNodeMask) : m_mesh1d(mesh1d), m_mesh2d(mesh2d), m_oneDNodeMask(oneDNodeMask)
{
    // assert mesh1d and mesh have the same projection!
    if (m_mesh1d->m_projection != m_mesh2d->m_projection)
    {
        throw AlgorithmError("meshkernel::Contacts::Contacts: m_mesh1d and m_mesh2d projections are different");
    }
}

void meshkernel::Contacts::ComputeSingleConnections(const Polygons& polygons)
{
    m_mesh2d->Administrate(Mesh2D::AdministrationOptions::AdministrateMeshEdgesAndFaces);
    m_mesh1d->AdministrateNodesEdges();

    const auto node1dFaceIndices = m_mesh2d->PointFaceIndices(m_mesh1d->m_nodes);
    m_mesh1dIndices.reserve(m_mesh1d->m_nodes.size());
    m_mesh2dIndices.reserve(m_mesh1d->m_nodes.size());
    const double distanceFactor = 5.0;
    for (size_t n = 0; n < m_mesh1d->m_nodes.size(); ++n)
    {
        // connect only nodes included in the polygons
        if (!polygons.IsPointInPolygons(m_mesh1d->m_nodes[n]))
        {
            continue;
        }

        // do not connect nodes at boundary of the 1d mesh
        if (m_mesh1d->IsNodeOnBoundary(n))
        {
            continue;
        }

        // if m_oneDNodeMask is not empty, connect only if the mask value for the current node is true
        if (!m_oneDNodeMask.empty() && !m_oneDNodeMask[n])
        {
            continue;
        }

        // if a node is inside a face, connect the 1d node with the face including the node. No more work to do
        if (node1dFaceIndices[n] != sizetMissingValue)
        {
            m_mesh1dIndices.emplace_back(n);
            m_mesh2dIndices.emplace_back(node1dFaceIndices[n]);
            continue;
        }

        // Connect faces crossing the right projected segment
        Connect1dNodesWithCrossingFace(n, 5.0);
        // Connect faces crossing the left projected segment
        Connect1dNodesWithCrossingFace(n, -5.0);
    }
};

void meshkernel::Contacts::Connect1dNodesWithCrossingFace(size_t node, double distanceFactor)
{

    const auto left1dEdge = m_mesh1d->m_nodesEdges[node][0];
    const auto right1dEdge = m_mesh1d->m_nodesEdges[node][1];

    const auto otherLeft1dNode = m_mesh1d->m_edges[left1dEdge].first == node ? m_mesh1d->m_edges[left1dEdge].second : m_mesh1d->m_edges[left1dEdge].first;
    const auto otherRight1dNode = m_mesh1d->m_edges[right1dEdge].first == node ? m_mesh1d->m_edges[right1dEdge].second : m_mesh1d->m_edges[right1dEdge].first;

    const auto normalVector = NormalVectorOutside(m_mesh1d->m_nodes[otherLeft1dNode], m_mesh1d->m_nodes[otherRight1dNode], m_mesh1d->m_projection);
    const auto edgeLength = ComputeDistance(m_mesh1d->m_nodes[otherLeft1dNode], m_mesh1d->m_nodes[otherRight1dNode], m_mesh1d->m_projection);

    const auto projectedNode = m_mesh1d->m_nodes[node] + normalVector * edgeLength * distanceFactor;
    size_t intersectedFace;
    size_t intersectedEdge;

    const auto isConnectionIntersectingAFace = m_mesh2d->IsSegmentCrossingAFace(m_mesh1d->m_nodes[node], projectedNode, intersectedFace, intersectedEdge);
    if (isConnectionIntersectingAFace &&
        !IsConnectionIntersectingMesh1d(node, intersectedFace) &&
        !IsContactIntersectingContact(node, intersectedFace))
    {
        m_mesh1dIndices.emplace_back(node);
        m_mesh2dIndices.emplace_back(intersectedFace);
    }
}

bool meshkernel::Contacts::IsConnectionIntersectingMesh1d(size_t node, size_t face) const
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

        if (areSegmentCrossing &&
            ratioFirstSegment > 0.0 && ratioFirstSegment < 1.0 &&
            ratioSecondSegment > 0.0 && ratioSecondSegment < 1.0)
        {
            return true;
        }
    }
    return false;
}

bool meshkernel::Contacts::IsContactIntersectingContact(size_t node, size_t face) const
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
        if (areSegmentCrossing &&
            ratioFirstSegment > 0.0 && ratioFirstSegment < 1.0 &&
            ratioSecondSegment > 0.0 && ratioSecondSegment < 1.0)
        {
            return true;
        }
    }

    return false;
}

void meshkernel::Contacts::ComputeMultipleConnections()
{
    // Perform mesh2d administration
    m_mesh2d->Administrate(Mesh2D::AdministrationOptions::AdministrateMeshEdgesAndFaces);

    // Perform mesh1d administration
    m_mesh1d->AdministrateNodesEdges();
    m_mesh1d->ComputeEdgesLengths();

    // Compute the indices of the faces including the 1d nodes
    const auto node1dFaceIndices = m_mesh2d->PointFaceIndices(m_mesh1d->m_nodes);

    // Build mesh2d face circumcenters r-tree
    RTree faceCircumcentersRTree;
    faceCircumcentersRTree.BuildTree(m_mesh2d->m_facesCircumcenters);
    std::vector<bool> isFaceAlreadyConnected(m_mesh2d->GetNumFaces(), false);

    for (auto i = 0; i < m_mesh1d->GetNumEdges(); ++i)
    {

        const auto firstNode1dMeshEdge = m_mesh1d->m_edges[i].first;
        const auto secondNode1dMeshEdge = m_mesh1d->m_edges[i].second;

        // loop over all edges connected to the first node, to determine the longest edge
        auto maxEdgeLenght = std::numeric_limits<double>::lowest();
        for (auto j = 0; j < m_mesh1d->m_nodesNumEdges[firstNode1dMeshEdge]; ++j)
        {
            const auto edge = m_mesh1d->m_nodesEdges[firstNode1dMeshEdge][j];
            maxEdgeLenght = std::max(maxEdgeLenght, m_mesh1d->m_edgeLengths[edge]);
        }

        // Computes the nearest 2d face indices
        faceCircumcentersRTree.NearestNeighborsOnSquaredDistance(m_mesh1d->m_nodes[firstNode1dMeshEdge], 1.1 * maxEdgeLenght * maxEdgeLenght);

        // For each face determine if it is crossing the current 1d edge
        for (auto k = 0; k < faceCircumcentersRTree.GetQueryResultSize(); ++k)
        {
            const auto face = faceCircumcentersRTree.GetQueryResult(k);

            // The face is already connected to a 1d node, nothing to do
            if (isFaceAlreadyConnected[face])
            {
                continue;
            }

            // Determine which of the mesh2d edges is crossing the current 1d edge
            for (auto l = 0; l < m_mesh2d->m_numFacesNodes[face]; ++l)
            {
                Point intersectionPoint;
                double crossProduct;
                double ratioFirstSegment;
                double ratioSecondSegment;

                const auto edge = m_mesh2d->m_facesEdges[face][l];
                const auto firstNode2dMeshEdge = m_mesh2d->m_edges[edge].first;
                const auto secondNode2dMeshEdge = m_mesh2d->m_edges[edge].second;

                const auto areSegmentCrossing = AreSegmentsCrossing(m_mesh1d->m_nodes[firstNode1dMeshEdge],
                                                                    m_mesh1d->m_nodes[secondNode1dMeshEdge],
                                                                    m_mesh2d->m_nodes[firstNode2dMeshEdge],
                                                                    m_mesh2d->m_nodes[secondNode2dMeshEdge],
                                                                    false,
                                                                    m_mesh1d->m_projection,
                                                                    intersectionPoint,
                                                                    crossProduct,
                                                                    ratioFirstSegment,
                                                                    ratioSecondSegment);
                // Nothing is crossing, thus continue
                if (!areSegmentCrossing)
                {
                    continue;
                }

                //compute the distance between the face circumcenter and the crossed 1d edge nodes.
                const auto leftDistance = ComputeDistance(m_mesh1d->m_nodes[firstNode1dMeshEdge], m_mesh2d->m_facesCircumcenters[face], m_mesh1d->m_projection);
                const auto rightDistance = ComputeDistance(m_mesh1d->m_nodes[secondNode1dMeshEdge], m_mesh2d->m_facesCircumcenters[face], m_mesh1d->m_projection);
                const auto nodeToConnect = leftDistance <= rightDistance ? firstNode1dMeshEdge : secondNode1dMeshEdge;

                // if m_oneDNodeMask is not empty, connect only if the mask value for the current node is true
                if (!m_oneDNodeMask.empty() && !m_oneDNodeMask[nodeToConnect])
                {
                    continue;
                }

                // the 1d mesh node to be connected needs to be included in the 2d mesh
                if (node1dFaceIndices[nodeToConnect] == sizetMissingValue)
                {
                    continue;
                }

                m_mesh1dIndices.emplace_back(nodeToConnect);
                m_mesh2dIndices.emplace_back(face);
                isFaceAlreadyConnected[face] = true;
                break;
            }
        }
    }
};

void meshkernel::Contacts::ComputeConnectionsWithPolygons(const Polygons& polygons){
    // Complete implementation
};

void meshkernel::Contacts::ComputeConnectionsWithPoints(const std::vector<Point>& points){
    // Complete implementation
};

void meshkernel::Contacts::ComputeBoundaryConnections(){
    // Complete implementation
};
