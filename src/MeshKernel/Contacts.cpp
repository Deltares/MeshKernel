#include <memory>
#include <stdexcept>
#include <vector>

#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Contacts.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Polygons.hpp>

using meshkernel::Contacts;

Contacts::Contacts(std::shared_ptr<Mesh1D> mesh1d,
                   std::shared_ptr<Mesh2D> mesh2d) : m_mesh1d(mesh1d), m_mesh2d(mesh2d)
{
    // assert mesh1d and mesh have the same projection
    if (m_mesh1d->m_projection != m_mesh2d->m_projection)
    {
        throw AlgorithmError("meshkernel::Contacts::Contacts: m_mesh1d and m_mesh2d projections are different");
    }
}

void Contacts::ComputeSingleContacts(const std::vector<bool>& oneDNodeMask,
                                     const Polygons& polygons)
{
    // assert oneDNodeMask and m_mesh1d have the same number of nodes
    if (oneDNodeMask.size() != m_mesh1d->m_nodes.size())
    {
        throw std::invalid_argument("meshkernel::Contacts::ComputeSingleContacts: oneDNodeMask and m_mesh1d do not have the same number of nodes");
    }

    m_mesh2d->Administrate(Mesh2D::AdministrationOption::AdministrateMeshEdgesAndFaces);
    m_mesh1d->AdministrateNodesEdges();

    const auto node1dFaceIndices = m_mesh2d->PointFaceIndices(m_mesh1d->m_nodes);
    m_mesh1dIndices.reserve(m_mesh1d->m_nodes.size());
    m_mesh2dIndices.reserve(m_mesh1d->m_nodes.size());

    const auto nodePolygonIndices = polygons.PolygonIndices(m_mesh1d->m_nodes);

    for (size_t n = 0; n < m_mesh1d->m_nodes.size(); ++n)
    {
        // connect only nodes included in the polygons
        if (nodePolygonIndices[n] == sizetMissingValue)
        {
            continue;
        }

        // do not connect nodes at boundary of the 1d mesh
        if (m_mesh1d->IsNodeOnBoundary(n))
        {
            continue;
        }

        // if oneDNodeMask is not empty, connect only if the mask value for the current node is true
        if (!oneDNodeMask.empty() && !oneDNodeMask[n])
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

        // connect faces crossing the right projected segment
        Connect1dNodesWithCrossingFaces(n, 5.0);
        // connect faces crossing the left projected segment
        Connect1dNodesWithCrossingFaces(n, -5.0);
    }
};

void Contacts::Connect1dNodesWithCrossingFaces(size_t node,
                                               double distanceFactor)
{
    const auto left1dEdge = m_mesh1d->m_nodesEdges[node][0];
    const auto right1dEdge = m_mesh1d->m_nodesEdges[node][1];

    const auto otherLeft1dNode = m_mesh1d->m_edges[left1dEdge].first == node ? m_mesh1d->m_edges[left1dEdge].second : m_mesh1d->m_edges[left1dEdge].first;
    const auto otherRight1dNode = m_mesh1d->m_edges[right1dEdge].first == node ? m_mesh1d->m_edges[right1dEdge].second : m_mesh1d->m_edges[right1dEdge].first;

    const auto normalVector = NormalVectorOutside(m_mesh1d->m_nodes[otherLeft1dNode], m_mesh1d->m_nodes[otherRight1dNode], m_mesh1d->m_projection);
    const auto edgeLength = ComputeDistance(m_mesh1d->m_nodes[otherLeft1dNode], m_mesh1d->m_nodes[otherRight1dNode], m_mesh1d->m_projection);

    const auto projectedNode = m_mesh1d->m_nodes[node] + normalVector * edgeLength * distanceFactor;

    const auto [intersectedFace, intersectedEdge] = m_mesh2d->IsSegmentCrossingABoundaryEdge(m_mesh1d->m_nodes[node], projectedNode);
    if (intersectedFace != sizetMissingValue &&
        intersectedEdge != sizetMissingValue &&
        !IsContactIntersectingMesh1d(node, intersectedFace) &&
        !IsContactIntersectingContact(node, intersectedFace))
    {
        m_mesh1dIndices.emplace_back(node);
        m_mesh2dIndices.emplace_back(intersectedFace);
    }
}

bool Contacts::IsContactIntersectingMesh1d(size_t node,
                                           size_t face) const
{
    for (size_t e = 0; e < m_mesh1d->GetNumEdges(); ++e)
    {

        Point intersectionPoint;
        double crossProduct;
        double ratioFirstSegment;
        double ratioSecondSegment;
        if (AreSegmentsCrossing(m_mesh1d->m_nodes[node],
                                m_mesh2d->m_facesCircumcenters[face],
                                m_mesh1d->m_nodes[m_mesh1d->m_edges[e].first],
                                m_mesh1d->m_nodes[m_mesh1d->m_edges[e].second],
                                false,
                                m_mesh1d->m_projection,
                                intersectionPoint,
                                crossProduct,
                                ratioFirstSegment,
                                ratioSecondSegment) &&
            ratioFirstSegment > 0.0 && ratioFirstSegment < 1.0 &&
            ratioSecondSegment > 0.0 && ratioSecondSegment < 1.0)
        {
            return true;
        }
    }
    return false;
}

bool Contacts::IsContactIntersectingContact(size_t node, size_t face) const
{
    for (size_t i = 0; i < m_mesh1dIndices.size(); ++i)
    {
        Point intersectionPoint;
        double crossProduct;
        double ratioFirstSegment;
        double ratioSecondSegment;

        if (AreSegmentsCrossing(m_mesh1d->m_nodes[node],
                                m_mesh2d->m_facesCircumcenters[face],
                                m_mesh1d->m_nodes[m_mesh1dIndices[i]],
                                m_mesh2d->m_facesCircumcenters[m_mesh2dIndices[i]],
                                false,
                                m_mesh1d->m_projection,
                                intersectionPoint,
                                crossProduct,
                                ratioFirstSegment,
                                ratioSecondSegment) &&
            ratioFirstSegment > 0.0 && ratioFirstSegment < 1.0 &&
            ratioSecondSegment > 0.0 && ratioSecondSegment < 1.0)
        {
            return true;
        }
    }

    return false;
}

void Contacts::ComputeMultipleContacts(const std::vector<bool>& oneDNodeMask)
{
    // assert oneDNodeMask and m_mesh1d have the same number of nodes
    if (oneDNodeMask.size() != m_mesh1d->m_nodes.size())
    {
        throw std::invalid_argument("meshkernel::Contacts::ComputeSingleContacts: oneDNodeMask and m_mesh1d do not have the same number of nodes");
    }

    // perform mesh2d administration
    m_mesh2d->Administrate(Mesh2D::AdministrationOption::AdministrateMeshEdgesAndFaces);

    // perform mesh1d administration
    m_mesh1d->AdministrateNodesEdges();

    // compute the indices of the faces including the 1d nodes
    const auto node1dFaceIndices = m_mesh2d->PointFaceIndices(m_mesh1d->m_nodes);

    // build mesh2d face circumcenters r-tree
    std::vector<bool> isFaceAlreadyConnected(m_mesh2d->GetNumFaces(), false);

    // loop over 1d mesh edges
    for (auto e = 0; e < m_mesh1d->GetNumEdges(); ++e)
    {
        // get the mesh1d edge nodes
        const auto firstNode1dMeshEdge = m_mesh1d->m_edges[e].first;
        const auto secondNode1dMeshEdge = m_mesh1d->m_edges[e].second;

        // computes the maximum edge length
        const auto maxEdgeLength = m_mesh1d->ComputeMaxLengthSurroundingEdges(firstNode1dMeshEdge);

        // compute the nearest 2d face indices
        m_mesh2d->SearchPointsWithinSquaredRadius(m_mesh1d->m_nodes[firstNode1dMeshEdge], 1.1 * maxEdgeLength * maxEdgeLength, MeshLocations::Faces);

        // for each face determine if it is crossing the current 1d edge
        const auto numNeighbours = m_mesh2d->GetNumNearestNeighbors(MeshLocations::Faces);
        for (auto f = 0; f < numNeighbours; ++f)
        {
            const auto face = m_mesh2d->GetNearestNeighborIndex(f, MeshLocations::Faces);

            // the face is already connected to a 1d node, nothing to do
            if (isFaceAlreadyConnected[face])
            {
                continue;
            }

            // determine which of the mesh2d edges is crossing the current 1d edge
            for (auto ee = 0; ee < m_mesh2d->m_numFacesNodes[face]; ++ee)
            {
                Point intersectionPoint;
                double crossProduct;
                double ratioFirstSegment;
                double ratioSecondSegment;

                const auto edge = m_mesh2d->m_facesEdges[face][ee];
                const auto firstNode2dMeshEdge = m_mesh2d->m_edges[edge].first;
                const auto secondNode2dMeshEdge = m_mesh2d->m_edges[edge].second;

                // nothing is crossing, continue
                if (!AreSegmentsCrossing(m_mesh1d->m_nodes[firstNode1dMeshEdge],
                                         m_mesh1d->m_nodes[secondNode1dMeshEdge],
                                         m_mesh2d->m_nodes[firstNode2dMeshEdge],
                                         m_mesh2d->m_nodes[secondNode2dMeshEdge],
                                         false,
                                         m_mesh1d->m_projection,
                                         intersectionPoint,
                                         crossProduct,
                                         ratioFirstSegment,
                                         ratioSecondSegment))
                {
                    continue;
                }

                // compute the distance between the face circumcenter and the crossed 1d edge nodes.
                const auto leftDistance = ComputeDistance(m_mesh1d->m_nodes[firstNode1dMeshEdge], m_mesh2d->m_facesCircumcenters[face], m_mesh1d->m_projection);
                const auto rightDistance = ComputeDistance(m_mesh1d->m_nodes[secondNode1dMeshEdge], m_mesh2d->m_facesCircumcenters[face], m_mesh1d->m_projection);
                const auto nodeToConnect = leftDistance <= rightDistance ? firstNode1dMeshEdge : secondNode1dMeshEdge;

                // if oneDNodeMask is not empty, connect only if the mask value for the current node is true
                if (!oneDNodeMask.empty() && !oneDNodeMask[nodeToConnect])
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

void Contacts::ComputeContactsWithPolygons(const std::vector<bool>& oneDNodeMask,
                                           const Polygons& polygons)
{
    // assert oneDNodeMask and m_mesh1d have the same number of nodes
    if (oneDNodeMask.size() != m_mesh1d->m_nodes.size())
    {
        throw std::invalid_argument("meshkernel::Contacts::ComputeSingleContacts: oneDNodeMask and m_mesh1d do not have the same number of nodes");
    }

    // no valid polygons provided
    if (polygons.IsEmpty())
    {
        return;
    }

    // perform mesh2d administration
    m_mesh2d->Administrate(Mesh2D::AdministrationOption::AdministrateMeshEdgesAndFaces);

    // perform mesh1d administration
    m_mesh1d->AdministrateNodesEdges();

    // for each mesh2d face, store polygon index
    std::vector<size_t> polygonIndices(m_mesh2d->GetNumFaces(), sizetMissingValue);
    for (auto faceIndex = 0; faceIndex < m_mesh2d->GetNumFaces(); ++faceIndex)
    {
        polygonIndices[faceIndex] = polygons.PolygonIndex(m_mesh2d->m_facesMassCenters[faceIndex]);
    }

    // for each polygon, find closest 1d node to any 2d mass center within the polygon
    std::vector<double> minimalDistance(polygons.GetNumPolygons(), doubleMissingValue);
    std::vector<size_t> closest1dNodeIndices(polygons.GetNumPolygons(), sizetMissingValue);
    std::vector<size_t> closest2dNodeIndices(polygons.GetNumPolygons(), sizetMissingValue);
    for (auto faceIndex = 0; faceIndex < m_mesh2d->GetNumFaces(); ++faceIndex)
    {
        const auto polygonIndex = polygonIndices[faceIndex];
        // if face is not within a polygon, continue
        if (polygonIndex == sizetMissingValue)
        {
            continue;
        }
        const auto faceMassCenter = m_mesh2d->m_facesMassCenters[faceIndex];

        const auto close1DNodeIndex = m_mesh1d->FindNodeCloseToAPoint(faceMassCenter, oneDNodeMask);

        const auto close1DNode = m_mesh1d->m_nodes[close1DNodeIndex];
        const auto squaredDistance = ComputeSquaredDistance(faceMassCenter, close1DNode, m_mesh2d->m_projection);
        // if it is the first found node of this polygon or
        // there is already a distance stored, but ours is smaller
        // -> store
        if (IsEqual(minimalDistance[polygonIndex], doubleMissingValue) || squaredDistance < minimalDistance[polygonIndex])
        {
            closest1dNodeIndices[polygonIndex] = close1DNodeIndex;
            closest2dNodeIndices[polygonIndex] = faceIndex;
            minimalDistance[polygonIndex] = squaredDistance;
        }
    }

    // connect 1D nodes to closest 2d node in a polygon
    for (auto polygonIndex = 0; polygonIndex < polygons.GetNumPolygons(); ++polygonIndex)
    {
        m_mesh1dIndices.emplace_back(closest1dNodeIndices[polygonIndex]);
        m_mesh2dIndices.emplace_back(closest2dNodeIndices[polygonIndex]);
    }
};

void Contacts::ComputeContactsWithPoints(const std::vector<bool>& oneDNodeMask,
                                         const std::vector<Point>& points)
{
    // assert oneDNodeMask and m_mesh1d have the same number of nodes
    if (oneDNodeMask.size() != m_mesh1d->m_nodes.size())
    {
        throw std::invalid_argument("meshkernel::Contacts::ComputeSingleContacts: oneDNodeMask and m_mesh1d do not have the same number of nodes");
    }

    // perform mesh2d administration
    m_mesh2d->Administrate(Mesh2D::AdministrationOption::AdministrateMeshEdgesAndFaces);

    // perform mesh1d administration (m_nodesRTree will also be build if necessary)
    m_mesh1d->AdministrateNodesEdges();

    // find the face indices containing the 1d points
    const auto pointsFaceIndices = m_mesh2d->PointFaceIndices(points);

    // for each 1d node in the 2d mesh, find the closest 1d node.
    for (auto i = 0; i < points.size(); ++i)
    {
        // point not in the mesh
        if (pointsFaceIndices[i] == sizetMissingValue)
        {
            continue;
        }

        // get the closest 1d node
        m_mesh1d->SearchNearestNeighbor(points[i], MeshLocations::Nodes);

        // if nothing found continue
        if (m_mesh1d->GetNumNearestNeighbors(MeshLocations::Nodes) == 0)
        {
            continue;
        }

        // form the 1d-2d contact
        m_mesh1dIndices.emplace_back(m_mesh1d->GetNearestNeighborIndex(0, MeshLocations::Nodes));
        m_mesh2dIndices.emplace_back(pointsFaceIndices[i]);
    }
};

void Contacts::ComputeBoundaryContacts(const std::vector<bool>& oneDNodeMask,
                                       const Polygons& polygons,
                                       double searchRadius)
{
    // assert oneDNodeMask and m_mesh1d have the same number of nodes
    if (oneDNodeMask.size() != m_mesh1d->m_nodes.size())
    {
        throw std::invalid_argument("meshkernel::Contacts::ComputeSingleContacts: oneDNodeMask and m_mesh1d do not have the same number of nodes");
    }

    // perform mesh2d administration
    m_mesh2d->Administrate(Mesh2D::AdministrationOption::AdministrateMeshEdgesAndFaces);

    // perform mesh1d administration
    m_mesh1d->AdministrateNodesEdges();

    // build mesh2d face circumcenters r-tree
    RTree faceCircumcentersRTree;
    faceCircumcentersRTree.BuildTree(m_mesh2d->m_facesCircumcenters);

    // get the indices
    const auto facePolygonIndices = polygons.PolygonIndices(m_mesh2d->m_facesCircumcenters);

    // Loop over 1d edges
    std::vector<bool> isValidFace(m_mesh2d->GetNumFaces(), true);
    std::vector<size_t> faceTo1DNode(m_mesh2d->GetNumFaces(), sizetMissingValue);
    for (auto n = 0; n < m_mesh1d->GetNumNodes(); ++n)
    {

        // account for the 1d node mask if present
        if (!oneDNodeMask.empty() && !oneDNodeMask[n])
        {
            continue;
        }

        auto localSearchRadius = searchRadius;
        if (IsEqual(searchRadius, doubleMissingValue))
        {
            localSearchRadius = m_mesh1d->ComputeMaxLengthSurroundingEdges(n);
        }

        // compute the nearest 2d face indices
        faceCircumcentersRTree.PointsWithinSearchRadius(m_mesh1d->m_nodes[n], localSearchRadius * localSearchRadius);

        //// for each face determine if it is crossing the current 1d edge
        //auto numResults = faceCircumcentersRTree.GetQueryResultSize();
        //std::cout << "The number of results is " << numResults << std::endl;

        for (auto f = 0; f < faceCircumcentersRTree.GetQueryResultSize(); ++f)
        {
            const auto face = faceCircumcentersRTree.GetQueryResult(f);

            // the face is already marked as invalid, nothing to do
            if (!isValidFace[face])
            {
                continue;
            }

            // not a boundary face
            if (!m_mesh2d->IsFaceOnBoundary(face))
            {
                isValidFace[face] = false;
                continue;
            }

            // the face is not inside a polygon
            if (facePolygonIndices[face] == sizetMissingValue)
            {
                isValidFace[face] = false;
                continue;
            }

            // a candidate contact does not exist
            if (faceTo1DNode[face] == sizetMissingValue)
            {
                faceTo1DNode[face] = n;
                continue;
            }

            const auto currentSquaredDistance = ComputeSquaredDistance(m_mesh1d->m_nodes[n],
                                                                       m_mesh2d->m_facesCircumcenters[face],
                                                                       m_mesh2d->m_projection);
            const auto previousSquaredDistance = ComputeSquaredDistance(m_mesh1d->m_nodes[faceTo1DNode[face]],
                                                                        m_mesh2d->m_facesCircumcenters[face],
                                                                        m_mesh2d->m_projection);

            if (currentSquaredDistance < previousSquaredDistance)
            {
                faceTo1DNode[face] = n;
            }
        }
    }

    for (auto f = 0; f < m_mesh2d->GetNumFaces(); ++f)
    {
        if (faceTo1DNode[f] != sizetMissingValue && isValidFace[f])
        {

            m_mesh1dIndices.emplace_back(faceTo1DNode[f]);
            m_mesh2dIndices.emplace_back(f);
        }
    }
};
