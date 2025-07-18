//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2025.
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

#include "MeshKernel/Contacts.hpp"

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/MeshEdgeLength.hpp"
#include "MeshKernel/MeshFaceCenters.hpp"
#include "MeshKernel/Operations.hpp"
#include "MeshKernel/Polygons.hpp"
#include "MeshKernel/Utilities/RTreeFactory.hpp"

using meshkernel::Contacts;

Contacts::Contacts(Mesh1D& mesh1d, Mesh2D& mesh2d)
    : m_mesh1d(mesh1d),
      m_mesh2d(mesh2d)
{
    // assert mesh1d and mesh have the same projection
    if (m_mesh1d.m_projection != m_mesh2d.m_projection)
    {
        throw AlgorithmError("meshkernel::Contacts::Contacts: m_mesh1d and m_mesh2d projections are different");
    }

    m_facesCircumcenters = algo::ComputeFaceCircumcenters(mesh2d);
}

void Contacts::ComputeSingleContacts(const std::vector<bool>& oneDNodeMask,
                                     const Polygons& polygons,
                                     double projectionFactor)
{

    // assert oneDNodeMask and m_mesh1d have the same number of nodes
    if (oneDNodeMask.size() != m_mesh1d.GetNumNodes())
    {
        throw AlgorithmError("oneDNodeMask and m_mesh1d do not have the same number of nodes ({} and {}, respectively)",
                             oneDNodeMask.size(),
                             m_mesh1d.GetNumNodes());
    }

    m_mesh1d.AdministrateNodesEdges();

    Validate();

    const auto node1dFaceIndices = m_mesh2d.PointFaceIndices(m_mesh1d.Nodes());
    m_mesh1dIndices.reserve(m_mesh1d.GetNumNodes());
    m_mesh2dIndices.reserve(m_mesh1d.GetNumNodes());

    const auto nodePolygonIndices = polygons.PointsInPolygons(m_mesh1d.Nodes());

    for (UInt n = 0; n < m_mesh1d.GetNumNodes(); ++n)
    {
        // connect only nodes included in the polygons
        if (!nodePolygonIndices[n])
        {
            continue;
        }

        // if oneDNodeMask is not empty, connect only if the mask value for the current node is true
        if (!oneDNodeMask.empty() && !oneDNodeMask[n])
        {
            continue;
        }

        // if a node is inside a face, connect the 1d node with the face including the node. No more work to do
        if (node1dFaceIndices[n] != constants::missing::uintValue)
        {
            m_mesh1dIndices.emplace_back(n);
            m_mesh2dIndices.emplace_back(node1dFaceIndices[n]);
            continue;
        }

        // connect faces crossing the right projected segment
        Connect1dNodesWithCrossingFaces(n, projectionFactor);
        // connect faces crossing the left projected segment
        Connect1dNodesWithCrossingFaces(n, -projectionFactor);
    }
}

void Contacts::Connect1dNodesWithCrossingFaces(UInt node,
                                               double projectionFactor)
{
    const auto projectedNode = m_mesh1d.ComputeProjectedNode(node, projectionFactor);

    const auto [intersectedFace, intersectedEdge] = m_mesh2d.IsSegmentCrossingABoundaryEdge(m_mesh1d.Node(node), projectedNode);
    if (intersectedFace != constants::missing::uintValue &&
        intersectedEdge != constants::missing::uintValue &&
        !IsContactIntersectingMesh1d(node, intersectedFace) &&
        !IsContactIntersectingContact(node, intersectedFace))
    {
        m_mesh1dIndices.emplace_back(node);
        m_mesh2dIndices.emplace_back(intersectedFace);
    }
}

bool Contacts::IsContactIntersectingMesh1d(UInt node,
                                           UInt face) const
{
    for (UInt e = 0; e < m_mesh1d.GetNumEdges(); ++e)
    {

        const auto [areSegmentCrossing,
                    intersectionPoint,
                    crossProduct,
                    ratioFirstSegment,
                    ratioSecondSegment] = AreSegmentsCrossing(m_mesh1d.Node(node),
                                                              m_facesCircumcenters[face],
                                                              m_mesh1d.Node(m_mesh1d.GetEdge(e).first),
                                                              m_mesh1d.Node(m_mesh1d.GetEdge(e).second),
                                                              false,
                                                              m_mesh1d.m_projection);

        if (areSegmentCrossing &&
            ratioFirstSegment > 0.0 && ratioFirstSegment < 1.0 &&
            ratioSecondSegment > 0.0 && ratioSecondSegment < 1.0)
        {
            return true;
        }
    }
    return false;
}

bool Contacts::IsContactIntersectingContact(UInt node, UInt face) const
{
    for (UInt i = 0; i < m_mesh1dIndices.size(); ++i)
    {
        const auto [areSegmentCrossing,
                    intersectionPoint,
                    crossProduct,
                    ratioFirstSegment,
                    ratioSecondSegment] = AreSegmentsCrossing(m_mesh1d.Node(node),
                                                              m_facesCircumcenters[face],
                                                              m_mesh1d.Node(m_mesh1dIndices[i]),
                                                              m_facesCircumcenters[m_mesh2dIndices[i]],
                                                              false,
                                                              m_mesh1d.m_projection);
        if (areSegmentCrossing &&
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
    if (oneDNodeMask.size() != m_mesh1d.GetNumNodes())
    {
        throw AlgorithmError("oneDNodeMask and m_mesh1d do not have the same number of nodes ({} and {}, respectively)",
                             oneDNodeMask.size(),
                             m_mesh1d.GetNumNodes());
    }

    // perform mesh1d administration
    m_mesh1d.AdministrateNodesEdges();
    std::vector<double> mesh1dEdgeLengths = algo::ComputeMeshEdgeLength(m_mesh1d);

    Validate();

    // compute the indices of the faces including the 1d nodes
    const auto node1dFaceIndices = m_mesh2d.PointFaceIndices(m_mesh1d.Nodes());

    // build mesh2d face circumcenters r-tree
    std::vector<bool> isFaceAlreadyConnected(m_mesh2d.GetNumFaces(), false);
    m_mesh2d.BuildTree(Location::Faces);
    auto& rtree = m_mesh2d.GetRTree(Location::Faces);

    // loop over 1d mesh edges
    for (UInt e = 0; e < m_mesh1d.GetNumEdges(); ++e)
    {
        // get the mesh1d edge nodes
        const auto firstNode1dMeshEdge = m_mesh1d.GetEdge(e).first;
        const auto secondNode1dMeshEdge = m_mesh1d.GetEdge(e).second;

        // computes the maximum edge length
        const auto maxEdgeLength = algo::MaxLengthSurroundingEdges(m_mesh1d, firstNode1dMeshEdge, mesh1dEdgeLengths);

        // compute the nearest 2d face indices
        rtree.SearchPoints(m_mesh1d.Node(firstNode1dMeshEdge), 1.1 * maxEdgeLength * maxEdgeLength);

        // for each face determine if it is crossing the current 1d edge
        const auto numNeighbours = rtree.GetQueryResultSize();
        for (UInt f = 0; f < numNeighbours; ++f)
        {
            const auto face = rtree.GetQueryResult(f);

            // the face is already connected to a 1d node, nothing to do
            if (isFaceAlreadyConnected[face])
            {
                continue;
            }

            // determine which of the mesh2d edges is crossing the current 1d edge
            for (UInt ee = 0; ee < m_mesh2d.m_numFacesNodes[face]; ++ee)
            {
                const auto edge = m_mesh2d.m_facesEdges[face][ee];
                const auto firstNode2dMeshEdge = m_mesh2d.GetEdge(edge).first;
                const auto secondNode2dMeshEdge = m_mesh2d.GetEdge(edge).second;

                const auto [areCrossing,
                            intersectionPoint,
                            crossProduct,
                            ratioFirstSegment,
                            ratioSecondSegment] =
                    AreSegmentsCrossing(m_mesh1d.Node(firstNode1dMeshEdge),
                                        m_mesh1d.Node(secondNode1dMeshEdge),
                                        m_mesh2d.Node(firstNode2dMeshEdge),
                                        m_mesh2d.Node(secondNode2dMeshEdge),
                                        false,
                                        m_mesh1d.m_projection);

                // nothing is crossing, continue
                if (!areCrossing)
                {
                    continue;
                }

                // compute the distance between the face circumcenter and the crossed 1d edge nodes.
                const auto leftDistance = ComputeDistance(m_mesh1d.Node(firstNode1dMeshEdge), m_facesCircumcenters[face], m_mesh1d.m_projection);
                const auto rightDistance = ComputeDistance(m_mesh1d.Node(secondNode1dMeshEdge), m_facesCircumcenters[face], m_mesh1d.m_projection);
                const auto nodeToConnect = leftDistance <= rightDistance ? firstNode1dMeshEdge : secondNode1dMeshEdge;

                // if oneDNodeMask is not empty, connect only if the mask value for the current node is true
                if (!oneDNodeMask.empty() && !oneDNodeMask[nodeToConnect])
                {
                    continue;
                }

                // the 1d mesh node to be connected needs to be included in the 2d mesh
                if (node1dFaceIndices[nodeToConnect] == constants::missing::uintValue)
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
}

void Contacts::ComputeContactsWithPolygons(const std::vector<bool>& oneDNodeMask,
                                           const Polygons& polygons)
{
    // assert oneDNodeMask and m_mesh1d have the same number of nodes
    if (oneDNodeMask.size() != m_mesh1d.GetNumNodes())
    {
        throw AlgorithmError("oneDNodeMask and m_mesh1d do not have the same number of nodes ({} and {}, respectively)",
                             oneDNodeMask.size(),
                             m_mesh1d.GetNumNodes());
    }

    // no valid polygons provided
    if (polygons.IsEmpty())
    {
        return;
    }

    // perform mesh1d administration
    m_mesh1d.AdministrateNodesEdges();
    Validate();

    // for each mesh2d face, store polygon index
    std::vector<UInt> facePolygonIndex(m_mesh2d.GetNumFaces(), constants::missing::uintValue);
    std::vector<bool> faceInPolygon(m_mesh2d.GetNumFaces(), false);
    for (UInt faceIndex = 0; faceIndex < m_mesh2d.GetNumFaces(); ++faceIndex)
    {
        auto [isInPolygon, polygonIndex] = polygons.IsPointInPolygons(m_mesh2d.m_facesMassCenters[faceIndex]);
        faceInPolygon[faceIndex] = isInPolygon;
        facePolygonIndex[faceIndex] = polygonIndex;
    }

    // for each polygon, find closest 1d node to any 2d mass center within the polygon
    std::vector<double> minimalDistance(polygons.GetNumPolygons(), constants::missing::doubleValue);
    std::vector<UInt> closest1dNodeIndices(polygons.GetNumPolygons(), constants::missing::uintValue);
    std::vector<UInt> closest2dNodeIndices(polygons.GetNumPolygons(), constants::missing::uintValue);
    m_mesh1d.BuildTree(Location::Nodes);

    for (UInt faceIndex = 0; faceIndex < m_mesh2d.GetNumFaces(); ++faceIndex)
    {
        // if face is not within a polygon, continue
        if (!faceInPolygon[faceIndex])
        {
            continue;
        }
        const auto polygonIndex = facePolygonIndex[faceIndex];
        const auto faceMassCenter = m_mesh2d.m_facesMassCenters[faceIndex];
        const auto close1DNodeIndex = m_mesh1d.FindLocationIndex(faceMassCenter, Location::Nodes, oneDNodeMask);

        const auto close1DNode = m_mesh1d.Node(close1DNodeIndex);
        const auto squaredDistance = ComputeSquaredDistance(faceMassCenter, close1DNode, m_mesh2d.m_projection);
        // if it is the first found node of this polygon or
        // there is already a distance stored, but ours is smaller
        // -> store
        if (IsEqual(minimalDistance[polygonIndex], constants::missing::doubleValue) || squaredDistance < minimalDistance[polygonIndex])
        {
            closest1dNodeIndices[polygonIndex] = close1DNodeIndex;
            closest2dNodeIndices[polygonIndex] = faceIndex;
            minimalDistance[polygonIndex] = squaredDistance;
        }
    }

    // connect 1D nodes to closest 2d node in a polygon
    for (UInt polygonIndex = 0; polygonIndex < polygons.GetNumPolygons(); ++polygonIndex)
    {
        m_mesh1dIndices.emplace_back(closest1dNodeIndices[polygonIndex]);
        m_mesh2dIndices.emplace_back(closest2dNodeIndices[polygonIndex]);
    }
}

void Contacts::ComputeContactsWithPoints(const std::vector<bool>& oneDNodeMask,
                                         const std::vector<Point>& points)
{
    // assert oneDNodeMask and m_mesh1d have the same number of nodes
    if (oneDNodeMask.size() != m_mesh1d.GetNumNodes())
    {
        throw AlgorithmError("oneDNodeMask and m_mesh1d do not have the same number of nodes ({} and {}, respectively)",
                             oneDNodeMask.size(),
                             m_mesh1d.GetNumNodes());
    }

    if (m_mesh1d.GetNumNodes() == 0)
    {
        return;
    }

    // perform mesh1d administration (m_nodesRTree will also be build if necessary)
    m_mesh1d.AdministrateNodesEdges();

    Validate();

    m_mesh1d.BuildTree(Location::Nodes);
    auto& rtree = m_mesh1d.GetRTree(Location::Nodes);

    // find the face indices containing the 1d points
    const auto pointsFaceIndices = m_mesh2d.PointFaceIndices(points);

    // for each 1d node in the 2d mesh, find the closest 1d node.
    for (UInt i = 0; i < points.size(); ++i)
    {

        // point not in the mesh
        if (pointsFaceIndices[i] == constants::missing::uintValue)
        {
            continue;
        }

        // get the closest 1d node
        rtree.SearchNearestPoint(points[i]);

        if (rtree.GetQueryResultSize() == 0)
        {
            continue;
        }

        UInt closest1dNodeindex = rtree.GetQueryResult(0);

        // form the 1d-2d contact
        // Account for 1d node mask
        if (closest1dNodeindex != constants::missing::uintValue && oneDNodeMask[closest1dNodeindex])
        {
            m_mesh1dIndices.emplace_back(closest1dNodeindex);
            m_mesh2dIndices.emplace_back(pointsFaceIndices[i]);
        }
    }
}

void Contacts::ComputeBoundaryContacts(const std::vector<bool>& oneDNodeMask,
                                       const Polygons& polygons,
                                       double searchRadius)
{
    // assert oneDNodeMask and m_mesh1d have the same number of nodes
    if (oneDNodeMask.size() != m_mesh1d.GetNumNodes())
    {
        throw AlgorithmError("oneDNodeMask and m_mesh1d do not have the same number of nodes ({} and {}, respectively)",
                             oneDNodeMask.size(),
                             m_mesh1d.GetNumNodes());
    }

    // perform mesh1d administration
    m_mesh1d.AdministrateNodesEdges();

    std::vector<double> mesh1dEdgeLengths = algo::ComputeMeshEdgeLength(m_mesh1d);

    Validate();

    // build mesh2d face circumcenters r-tree
    const auto faceCircumcentersRTree = RTreeFactory::Create(m_mesh2d.m_projection);
    faceCircumcentersRTree->BuildTree(m_facesCircumcenters);

    // get the indices
    const auto facePolygonIndices = polygons.PointsInPolygons(m_facesCircumcenters);

    bool computeLocalSearchRadius = true;
    double localSearchRadius = 0.0;
    if (!IsEqual(searchRadius, constants::missing::doubleValue))
    {
        computeLocalSearchRadius = false;
        localSearchRadius = searchRadius;
    }

    // Loop over 1d edges
    std::vector<bool> isValidFace(m_mesh2d.GetNumFaces(), true);
    std::vector<UInt> faceTo1DNode(m_mesh2d.GetNumFaces(), constants::missing::uintValue);
    for (UInt n = 0; n < m_mesh1d.GetNumNodes(); ++n)
    {
        // Account for 1d node mask if present
        if (!oneDNodeMask.empty() && !oneDNodeMask[n])
        {
            continue;
        }

        if (computeLocalSearchRadius)
        {
            localSearchRadius = algo::MaxLengthSurroundingEdges(m_mesh1d, n, mesh1dEdgeLengths);
        }

        // compute the nearest 2d face indices
        faceCircumcentersRTree->SearchPoints(m_mesh1d.Node(n), localSearchRadius * localSearchRadius);

        for (UInt f = 0; f < faceCircumcentersRTree->GetQueryResultSize(); ++f)
        {
            const auto face = faceCircumcentersRTree->GetQueryResult(f);

            // the face is already marked as invalid, nothing to do
            if (!isValidFace[face])
            {
                continue;
            }

            // not a boundary face
            if (!m_mesh2d.IsFaceOnBoundary(face))
            {
                isValidFace[face] = false;
                continue;
            }

            // the face is not inside a polygon
            if (!facePolygonIndices[face])
            {
                isValidFace[face] = false;
                continue;
            }

            // a candidate contact does not exist
            if (faceTo1DNode[face] == constants::missing::uintValue)
            {
                faceTo1DNode[face] = n;
                continue;
            }

            const auto currentSquaredDistance = ComputeSquaredDistance(m_mesh1d.Node(n),
                                                                       m_facesCircumcenters[face],
                                                                       m_mesh2d.m_projection);
            const auto previousSquaredDistance = ComputeSquaredDistance(m_mesh1d.Node(faceTo1DNode[face]),
                                                                        m_facesCircumcenters[face],
                                                                        m_mesh2d.m_projection);

            if (currentSquaredDistance < previousSquaredDistance)
            {
                faceTo1DNode[face] = n;
            }
        }
    }

    for (UInt f = 0; f < m_mesh2d.GetNumFaces(); ++f)
    {
        if (faceTo1DNode[f] != constants::missing::uintValue && isValidFace[f])
        {

            m_mesh1dIndices.emplace_back(faceTo1DNode[f]);
            m_mesh2dIndices.emplace_back(f);
        }
    }
}

void Contacts::SetIndices(const std::vector<meshkernel::UInt>& mesh1dIndices,
                          const std::vector<meshkernel::UInt>& mesh2dIndices)
{

    if (mesh1dIndices.size() != mesh2dIndices.size())
    {
        throw AlgorithmError("The size of the 1d mesh indices ({}) and that of the 2d mesh indices ({}) are not equal",
                             mesh1dIndices.size(),
                             mesh2dIndices.size());
    }

    m_mesh1dIndices = mesh1dIndices;
    m_mesh2dIndices = mesh2dIndices;
}

void Contacts::Validate() const
{

    if (m_mesh1d.GetNumNodes() == 0)
    {
        throw AlgorithmError("m_mesh1d has no nodes");
    }

    if (m_mesh2d.GetNumNodes() == 0)
    {
        throw AlgorithmError("m_mesh2d has no nodes");
    }
}
