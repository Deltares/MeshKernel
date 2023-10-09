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

#include <queue>

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/Operations.hpp"
#include "MeshKernel/Polygon.hpp"
#include "MeshKernel/Polygons.hpp"
#include "MeshKernel/TriangulationWrapper.hpp"
#include "MeshKernel/Utilities/RTree.hpp"

using meshkernel::Mesh2D;

Mesh2D::Mesh2D(Projection projection) : Mesh(projection) {}

Mesh2D::Mesh2D(const std::vector<Edge>& edges,
               const std::vector<Point>& nodes,
               Projection projection) : Mesh(edges, nodes, projection)
{
    Administrate();
}

Mesh2D::Mesh2D(const std::vector<Edge>& edges,
               const std::vector<Point>& nodes,
               const std::vector<std::vector<UInt>>& faceNodes,
               const std::vector<UInt>& numFaceNodes,
               Projection projection) : Mesh(edges, nodes, projection)
{

    AdministrateNodesEdges();

    ResizeAndInitializeFaceVectors();

    // The face nodes and the num face nodes
    m_facesNodes = faceNodes;
    m_numFacesNodes = numFaceNodes;

    std::vector<UInt> local_edges;
    std::vector<Point> local_nodes;
    std::vector<UInt> local_node_indices;
    for (UInt f = 0; f < m_facesNodes.size(); ++f)
    {
        local_edges.clear();
        local_nodes.clear();
        local_node_indices.clear();

        auto node = m_facesNodes[f][0];
        local_nodes.emplace_back(m_nodes[node]);
        local_node_indices.emplace_back(node);
        UInt index = 0;
        while (index < m_nodesEdges[node].size())
        {
            const auto edge = m_nodesEdges[node][index];
            const auto other_node = OtherNodeOfEdge(m_edges[edge], node);
            const auto edge_iter = std::find(local_edges.begin(), local_edges.end(), edge);
            if (std::find(m_facesNodes[f].begin(), m_facesNodes[f].end(), other_node) != m_facesNodes[f].end() &&
                edge_iter == local_edges.end() &&
                std::find(local_node_indices.begin(), local_node_indices.end(), other_node) == local_node_indices.end())

            {

                local_edges.emplace_back(edge);
                local_nodes.emplace_back(m_nodes[other_node]);
                local_node_indices.emplace_back(other_node);
                node = other_node;
                index = 0;
                continue;
            }
            if (other_node == m_facesNodes[f][0] && edge_iter == local_edges.end())
            {
                // loop closed
                local_edges.emplace_back(edge);
                break;
            }

            index++;
        }

        m_facesEdges.emplace_back(local_edges);
        for (const auto& e : local_edges)
        {
            if (m_edgesNumFaces[e] > 2)
            {
                throw AlgorithmError("AdministrateFromFaceNodes: m_edgesNumFaces > 2.");
            }
            m_edgesFaces[e][m_edgesNumFaces[e]] = f;
            m_edgesNumFaces[e] += 1;
        }

        local_nodes.emplace_back(local_nodes.front());

        auto [face_area, center_of_mass, direction] = Polygon::FaceAreaAndCenterOfMass(local_nodes, m_projection);

        m_faceArea.emplace_back(face_area);
        m_facesMassCenters.emplace_back(center_of_mass);
    }

    ComputeCircumcentersMassCentersAndFaceAreas();

    ClassifyNodes();
}

void Mesh2D::Administrate()
{
    AdministrateNodesEdges();

    // face administration
    ResizeAndInitializeFaceVectors();

    // find faces
    FindFaces();

    // find mesh circumcenters
    ComputeCircumcentersMassCentersAndFaceAreas();

    // classify node types
    ClassifyNodes();
}

Mesh2D::Mesh2D(const std::vector<Point>& inputNodes, const Polygons& polygons, Projection projection)
{
    m_projection = projection;
    // compute triangulation
    TriangulationWrapper triangulationWrapper;
    const auto numberOfTriangles = static_cast<UInt>(inputNodes.size()) * 6 + 10;
    triangulationWrapper.Compute(inputNodes,
                                 TriangulationWrapper::TriangulationOptions::TriangulatePointsAndGenerateFaces,
                                 0.0,
                                 numberOfTriangles);

    triangulationWrapper.BuildTriangulation();

    // For each triangle check
    // 1. Validity of its internal angles
    // 2. Is inside the polygon
    // If so we mark the edges and we add them m_edges
    std::vector<bool> edgeNodesFlag(triangulationWrapper.GetNumEdges(), false);
    for (auto i = 0; i < triangulationWrapper.GetNumFaces(); ++i)
    {
        const auto goodTriangle = HasTriangleNoAcuteAngles(triangulationWrapper.GetFaceNodes(i), inputNodes);

        if (!goodTriangle)
        {
            continue;
        }
        const Point approximateCenter = (inputNodes[triangulationWrapper.GetFaceNode(i, 0)] + inputNodes[triangulationWrapper.GetFaceNode(i, 1)] + inputNodes[triangulationWrapper.GetFaceNode(i, 2)]) * constants::numeric::oneThird;

        const auto isTriangleInPolygon = polygons.IsPointInPolygon(approximateCenter, 0);
        if (!isTriangleInPolygon)
        {
            continue;
        }

        // mark all edges of this triangle as good ones
        for (UInt j = 0; j < m_numNodesInTriangle; ++j)
        {
            edgeNodesFlag[triangulationWrapper.GetFaceEdge(i, j)] = true;
        }
    }

    // now add all points and all valid edges
    m_nodes = inputNodes;
    UInt validEdgesCount = 0;
    for (auto i = 0; i < triangulationWrapper.GetNumEdges(); ++i)
    {
        if (!edgeNodesFlag[i])
            continue;
        validEdgesCount++;
    }

    std::vector<Edge> edges(validEdgesCount);
    validEdgesCount = 0;
    for (auto i = 0; i < triangulationWrapper.GetNumEdges(); ++i)
    {
        if (!edgeNodesFlag[i])
            continue;

        edges[validEdgesCount].first = triangulationWrapper.GetEdgeNode(i, 0);
        edges[validEdgesCount].second = triangulationWrapper.GetEdgeNode(i, 1);
        validEdgesCount++;
    }

    m_nodesRTreeRequiresUpdate = true;
    m_edgesRTreeRequiresUpdate = true;

    *this = Mesh2D(edges, inputNodes, projection);
}

bool Mesh2D::HasTriangleNoAcuteAngles(const std::vector<UInt>& faceNodes, const std::vector<Point>& nodes) const
{
    // Used for triangular grids
    constexpr double triangleMinimumAngle = 5.0;
    constexpr double triangleMaximumAngle = 150.0;

    double phiMin = 1e3;
    double phiMax = 0.0;

    static std::array<std::array<UInt, 3>, 3> nodePermutations{{{2, 0, 1}, {0, 1, 2}, {1, 2, 0}}};

    for (UInt i = 0; i < faceNodes.size(); ++i)
    {
        Point x0 = nodes[faceNodes[nodePermutations[i][0]]];
        Point x1 = nodes[faceNodes[nodePermutations[i][1]]];
        Point x2 = nodes[faceNodes[nodePermutations[i][2]]];

        const auto cosphi = NormalizedInnerProductTwoSegments(x1, x0, x1, x2, m_projection);
        const auto phi = std::acos(std::min(std::max(cosphi, -1.0), 1.0)) * constants::conversion::radToDeg;
        phiMin = std::min(phiMin, phi);
        phiMax = std::max(phiMax, phi);
        if (phi < triangleMinimumAngle || phi > triangleMaximumAngle)
        {
            return false;
        }
    }
    return true;
}

void Mesh2D::DeleteDegeneratedTriangles()
{
    Administrate();

    // assume the max amount of degenerated triangles is 10% of the actual faces
    std::vector<UInt> degeneratedTriangles;
    degeneratedTriangles.reserve(static_cast<UInt>(static_cast<double>(GetNumFaces()) * 0.1));
    for (UInt f = 0; f < GetNumFaces(); ++f)
    {
        const auto numFaceNodes = m_numFacesNodes[f];
        if (numFaceNodes != m_numNodesInTriangle)
        {
            continue;
        }
        auto firstNode = m_facesNodes[f][0];
        auto secondNode = m_facesNodes[f][1];
        auto thirdNode = m_facesNodes[f][2];

        // account for periodic spherical coordinate
        if ((m_projection == Projection::spherical || m_projection == Projection::sphericalAccurate) && IsPointOnPole(m_nodes[firstNode]))
        {
            const auto saveFirstNode = firstNode;
            firstNode = secondNode;
            secondNode = thirdNode;
            thirdNode = saveFirstNode;
        }

        // compute coordinate differences, to check for collinearity
        const auto dx2 = GetDx(m_nodes[firstNode], m_nodes[secondNode], m_projection);
        const auto dy2 = GetDy(m_nodes[firstNode], m_nodes[secondNode], m_projection);
        const auto dx3 = GetDx(m_nodes[firstNode], m_nodes[thirdNode], m_projection);
        const auto dy3 = GetDy(m_nodes[firstNode], m_nodes[thirdNode], m_projection);

        const auto den = dy2 * dx3 - dy3 * dx2;

        if (IsEqual(den, 0.0))
        {
            // Flag edges to remove
            for (UInt e = 0; e < m_numNodesInTriangle; ++e)
            {
                const auto edge = m_facesEdges[f][e];
                m_edges[edge] = {constants::missing::uintValue, constants::missing::uintValue};
            }
            // save degenerated face index
            degeneratedTriangles.emplace_back(f);
        }
    }

    // collapse secondNode and thirdNode into firstNode, change coordinate of the firstNode to triangle center of mass
    for (auto const& face : degeneratedTriangles)
    {
        const auto firstNode = m_facesNodes[face][0];
        const auto secondNode = m_facesNodes[face][1];
        const auto thirdNode = m_facesNodes[face][2];

        m_nodes[thirdNode] = m_facesMassCenters[face];
        MergeTwoNodes(secondNode, firstNode);
        MergeTwoNodes(thirdNode, firstNode);
    }

    Administrate();
}

void Mesh2D::FindFacesRecursive(UInt startNode,
                                UInt node,
                                UInt previousEdge,
                                UInt numClosingEdges,
                                std::vector<UInt>& edges,
                                std::vector<UInt>& nodes,
                                std::vector<UInt>& sortedEdgesFaces,
                                std::vector<UInt>& sortedNodes,
                                std::vector<Point>& nodalValues)
{
    // The selected edge does not exist.
    if (nodes.size() >= numClosingEdges)
        return;

    if (m_edges[previousEdge].first == constants::missing::uintValue || m_edges[previousEdge].second == constants::missing::uintValue)
        throw std::invalid_argument("Mesh2D::FindFacesRecursive: The selected edge is invalid. This should not happen since all invalid edges should have been cleaned up.");

    // Check if the faces are already found
    if (m_edgesNumFaces[previousEdge] >= 2)
        return;

    edges.push_back(previousEdge);
    nodes.push_back(node);

    const auto otherNode = OtherNodeOfEdge(m_edges[previousEdge], node);

    // enclosure found
    if (otherNode == startNode && nodes.size() == numClosingEdges)
    {
        // no duplicated nodes allowed
        sortedNodes.clear();
        sortedNodes.reserve(nodes.size());
        std::copy(nodes.begin(), nodes.end(), std::back_inserter(sortedNodes));
        std::sort(sortedNodes.begin(), sortedNodes.end());
        for (UInt n = 0; n < sortedNodes.size() - 1; n++)
        {
            if (sortedNodes[n + 1] == sortedNodes[n])
            {
                return;
            }
        }

        // we need to add a face when at least one edge has no faces
        auto oneEdgeHasNoFace = false;
        for (const auto& edge : edges)
        {
            if (m_edgesNumFaces[edge] == 0)
            {
                oneEdgeHasNoFace = true;
                break;
            }
        }

        // check if least one edge has no face
        if (!oneEdgeHasNoFace)
        {
            sortedEdgesFaces.clear();
            sortedEdgesFaces.reserve(edges.size());
            // is an internal face only if all edges have a different face
            for (UInt ee = 0; ee < edges.size(); ee++)
            {
                sortedEdgesFaces.push_back(m_edgesFaces[edges[ee]][0]);
            }

            std::sort(sortedEdgesFaces.begin(), sortedEdgesFaces.end());

            for (UInt n = 0; n < sortedEdgesFaces.size() - 1; n++)
            {
                if (sortedEdgesFaces[n + 1] == sortedEdgesFaces[n])
                    return;
            }
        }

        // the order of the edges in a new face must be counterclockwise
        // in order to evaluate the clockwise order, the signed face area is computed
        nodalValues.clear();
        for (const auto& n : nodes)
        {
            nodalValues.emplace_back(m_nodes[n]);
        }
        nodalValues.emplace_back(nodalValues.front());

        auto const [area, center_of_mass, direction] = Polygon::FaceAreaAndCenterOfMass(nodalValues, m_projection);

        if (direction == TraversalDirection::Clockwise)
        {
            return;
        }

        // increase m_edgesNumFaces
        for (const auto& edge : edges)
        {
            // Increment the number of shared faces for the edge.
            ++m_edgesNumFaces[edge];
            const auto numFace = m_edgesNumFaces[edge];
            m_edgesFaces[edge][numFace - 1] = GetNumFaces();
        }

        // store the result
        m_facesNodes.emplace_back(nodes);
        m_facesEdges.emplace_back(edges);
        m_faceArea.emplace_back(area);
        m_facesMassCenters.emplace_back(center_of_mass);
        m_numFacesNodes.emplace_back(static_cast<UInt>(nodes.size()));

        return;
    }

    UInt edgeIndexOtherNode = 0;
    for (UInt e = 0; e < m_nodesNumEdges[otherNode]; e++)
    {
        if (m_nodesEdges[otherNode][e] == previousEdge)
        {
            edgeIndexOtherNode = e;
            break;
        }
    }

    if (edgeIndexOtherNode == 0)
    {
        edgeIndexOtherNode = m_nodesNumEdges[otherNode] - 1;
    }
    else if (edgeIndexOtherNode > m_nodesNumEdges[otherNode])
    {
        edgeIndexOtherNode = edgeIndexOtherNode - m_nodesNumEdges[otherNode] - 1;
    }
    else
    {
        edgeIndexOtherNode = edgeIndexOtherNode - 1;
    }

    auto const edge = m_nodesEdges[otherNode][edgeIndexOtherNode];
    FindFacesRecursive(startNode, otherNode, edge, numClosingEdges, edges, nodes, sortedEdgesFaces, sortedNodes, nodalValues);
}

void Mesh2D::FindFaces()
{
    std::vector<UInt> sortedEdgesFaces(m_maximumNumberOfEdgesPerFace);
    std::vector<UInt> sortedNodes(m_maximumNumberOfEdgesPerFace);
    std::vector<Point> nodalValues(m_maximumNumberOfEdgesPerFace);
    std::vector<UInt> edges(m_maximumNumberOfEdgesPerFace);
    std::vector<UInt> nodes(m_maximumNumberOfEdgesPerFace);

    for (UInt numEdgesPerFace = 3; numEdgesPerFace <= m_maximumNumberOfEdgesPerFace; numEdgesPerFace++)
    {
        for (UInt n = 0; n < GetNumNodes(); n++)
        {
            if (!m_nodes[n].IsValid())
            {
                continue;
            }

            for (UInt e = 0; e < m_nodesNumEdges[n]; e++)
            {
                nodes.clear();
                edges.clear();
                FindFacesRecursive(n, n, m_nodesEdges[n][e], numEdgesPerFace, edges, nodes, sortedEdgesFaces, sortedNodes, nodalValues);
            }
        }
    }
}

void Mesh2D::ComputeCircumcentersMassCentersAndFaceAreas(bool computeMassCenters)
{

    auto const numFaces = static_cast<int>(GetNumFaces());
    m_facesCircumcenters.resize(numFaces);
    m_faceArea.resize(numFaces);
    m_facesMassCenters.resize(numFaces);

    std::vector<UInt> numEdgeFacesCache;
    numEdgeFacesCache.reserve(m_maximumNumberOfEdgesPerFace);
    std::vector<Point> polygonNodesCache;
#pragma omp parallel for private(numEdgeFacesCache, polygonNodesCache)
    for (int f = 0; f < numFaces; f++)
    {
        // need to account for spherical coordinates. Build a polygon around a face
        ComputeFaceClosedPolygon(f, polygonNodesCache);

        if (computeMassCenters)
        {
            const auto [area, centerOfMass, direction] = Polygon::FaceAreaAndCenterOfMass(polygonNodesCache, m_projection);
            m_faceArea[f] = area;
            m_facesMassCenters[f] = centerOfMass;
        }

        UInt numberOfInteriorEdges = 0;
        const auto numberOfFaceNodes = static_cast<int>(GetNumFaceEdges(f));
        for (int n = 0; n < numberOfFaceNodes; n++)
        {
            if (!IsEdgeOnBoundary(m_facesEdges[f][n]))
            {
                numberOfInteriorEdges += 1;
            }
        }
        if (numberOfInteriorEdges == 0)
        {
            m_facesCircumcenters[f] = m_facesMassCenters[f];
            continue;
        }
        numEdgeFacesCache.clear();
        for (int n = 0; n < numberOfFaceNodes; n++)
        {
            numEdgeFacesCache.emplace_back(m_edgesNumFaces[m_facesEdges[f][n]]);
        }

        m_facesCircumcenters[f] = ComputeFaceCircumenter(polygonNodesCache, numEdgeFacesCache);
    }
}

void Mesh2D::ClassifyNodes()
{
    m_nodesTypes.resize(GetNumNodes(), 0);
    std::ranges::fill(m_nodesTypes, 0);

    for (UInt e = 0; e < GetNumEdges(); e++)
    {
        const auto firstNode = m_edges[e].first;
        const auto secondNode = m_edges[e].second;

        if (firstNode == constants::missing::uintValue || secondNode == constants::missing::uintValue)
        {
            continue;
        }

        if (m_nodesTypes[firstNode] == -1 || m_nodesTypes[secondNode] == -1)
        {
            continue;
        }

        if (m_edgesNumFaces[e] == 0)
        {
            m_nodesTypes[firstNode] = -1;
            m_nodesTypes[secondNode] = -1;
        }
        if (IsEdgeOnBoundary(e))
        {
            m_nodesTypes[firstNode] += 1;
            m_nodesTypes[secondNode] += 1;
        }
    }

    for (UInt n = 0; n < GetNumNodes(); n++)
    {
        if (m_nodesTypes[n] == 1 || m_nodesTypes[n] == 2)
        {
            if (m_nodesNumEdges[n] == 2)
            {
                // corner point
                m_nodesTypes[n] = 3;
            }
            else
            {
                UInt firstNode = constants::missing::uintValue;
                UInt secondNode = constants::missing::uintValue;
                for (UInt i = 0; i < m_nodesNumEdges[n]; ++i)
                {
                    const auto edgeIndex = m_nodesEdges[n][i];
                    if (!IsEdgeOnBoundary(edgeIndex))
                    {
                        continue;
                    }
                    if (firstNode == 0)
                    {
                        firstNode = OtherNodeOfEdge(m_edges[edgeIndex], n);
                    }
                    else
                    {
                        secondNode = OtherNodeOfEdge(m_edges[edgeIndex], n);
                        break;
                    }
                }

                // point at the border
                m_nodesTypes[n] = 2;
                if (firstNode != constants::missing::uintValue && secondNode != constants::missing::uintValue)
                {
                    const double cosPhi = NormalizedInnerProductTwoSegments(m_nodes[n], m_nodes[firstNode], m_nodes[n], m_nodes[secondNode], m_projection);

                    // threshold for corner points
                    const double cornerCosine = 0.25;
                    if (cosPhi > -cornerCosine)
                    {
                        // void angle
                        m_nodesTypes[n] = 3;
                    }
                }
            }
        }
        else if (m_nodesTypes[n] > 2)
        {
            // corner point
            m_nodesTypes[n] = 3;
        }
        else if (m_nodesTypes[n] != -1)
        {
            // internal node
            m_nodesTypes[n] = 1;
        }
        if (m_nodesNumEdges[n] < 2)
        {
            // hanging node
            m_nodesTypes[n] = -1;
        }
    }
}

void Mesh2D::ComputeFaceClosedPolygonWithLocalMappings(UInt faceIndex,
                                                       std::vector<Point>& polygonNodesCache,
                                                       std::vector<UInt>& localNodeIndicesCache,
                                                       std::vector<UInt>& globalEdgeIndicesCache) const
{
    const auto numFaceNodes = GetNumFaceEdges(faceIndex);
    polygonNodesCache.reserve(numFaceNodes + 1);
    polygonNodesCache.clear();
    localNodeIndicesCache.reserve(numFaceNodes + 1);
    localNodeIndicesCache.clear();
    globalEdgeIndicesCache.reserve(numFaceNodes + 1);
    globalEdgeIndicesCache.clear();

    for (UInt n = 0; n < numFaceNodes; n++)
    {
        polygonNodesCache.emplace_back(m_nodes[m_facesNodes[faceIndex][n]]);
        localNodeIndicesCache.emplace_back(n);
        globalEdgeIndicesCache.emplace_back(m_facesEdges[faceIndex][n]);
    }
    polygonNodesCache.emplace_back(polygonNodesCache.front());
    localNodeIndicesCache.emplace_back(0);
    globalEdgeIndicesCache.emplace_back(globalEdgeIndicesCache.front());
}

void Mesh2D::ComputeFaceClosedPolygon(UInt faceIndex, std::vector<Point>& polygonNodesCache) const
{
    const auto numFaceNodes = GetNumFaceEdges(faceIndex);
    polygonNodesCache.clear();
    polygonNodesCache.reserve(numFaceNodes + 1);

    for (UInt n = 0; n < numFaceNodes; n++)
    {
        polygonNodesCache.push_back(m_nodes[m_facesNodes[faceIndex][n]]);
    }

    polygonNodesCache.push_back(polygonNodesCache.front());
}

void Mesh2D::OffsetSphericalCoordinates(double minx, double maxx)
{
    if (m_projection == Projection::spherical && maxx - minx > 180.0)
    {
        for (UInt n = 0; n < GetNumNodes(); ++n)
        {
            if (m_nodes[n].x - 360.0 >= minx)
            {
                m_nodes[n].x -= 360.0;
            }

            if (m_nodes[n].x < minx)
            {
                m_nodes[n].x += 360.0;
            }
        }
    }
}

meshkernel::Point Mesh2D::ComputeFaceCircumenter(std::vector<Point>& polygon,
                                                 const std::vector<UInt>& edgesNumFaces) const
{
    const UInt maximumNumberCircumcenterIterations = 100;
    const double eps = m_projection == Projection::cartesian ? 1e-3 : 9e-10; // 111km = 0-e digit.

    std::array<Point, m_maximumNumberOfNodesPerFace> middlePoints;
    std::array<Point, m_maximumNumberOfNodesPerFace> normals;
    UInt pointCount = 0;

    const auto numNodes = static_cast<UInt>(polygon.size()) - 1;

    Point centerOfMass{0.0, 0.0};
    for (UInt n = 0; n < numNodes; n++)
    {
        centerOfMass.x += polygon[n].x;
        centerOfMass.y += polygon[n].y;
    }

    centerOfMass /= static_cast<double>(numNodes);

    auto result = centerOfMass;
    if (numNodes == m_numNodesInTriangle)
    {
        result = CircumcenterOfTriangle(polygon[0], polygon[1], polygon[2], m_projection);
    }
    else if (!edgesNumFaces.empty())
    {
        UInt numValidEdges = 0;
        for (UInt n = 0; n < numNodes; ++n)
        {
            if (edgesNumFaces[n] == 2)
            {
                numValidEdges++;
            }
        }

        if (numValidEdges > 1)
        {
            for (UInt n = 0; n < numNodes; n++)
            {
                if (edgesNumFaces[n] != 2)
                {
                    continue;
                }
                const auto nextNode = NextCircularForwardIndex(n, numNodes);

                middlePoints[pointCount] = ((polygon[n] + polygon[nextNode]) * 0.5);
                normals[pointCount] = NormalVector(polygon[n], polygon[nextNode], middlePoints[pointCount], m_projection);
                ++pointCount;
            }

            Point estimatedCircumCenter = centerOfMass;
            for (UInt iter = 0; iter < maximumNumberCircumcenterIterations; ++iter)
            {
                const Point previousCircumCenter = estimatedCircumCenter;
                for (UInt n = 0; n < pointCount; n++)
                {
                    const auto dx = GetDx(middlePoints[n], estimatedCircumCenter, m_projection);
                    const auto dy = GetDy(middlePoints[n], estimatedCircumCenter, m_projection);
                    const auto increment = -0.1 * DotProduct(dx, normals[n].x, dy, normals[n].y);
                    AddIncrementToPoint(normals[n], increment, centerOfMass, m_projection, estimatedCircumCenter);
                }
                if (iter > 0 &&
                    abs(estimatedCircumCenter.x - previousCircumCenter.x) < eps &&
                    abs(estimatedCircumCenter.y - previousCircumCenter.y) < eps)
                {
                    result = estimatedCircumCenter;
                    break;
                }
            }
        }
    }

    for (UInt n = 0; n < numNodes; n++)
    {
        polygon[n].x = m_weightCircumCenter * polygon[n].x + (1.0 - m_weightCircumCenter) * centerOfMass.x;
        polygon[n].y = m_weightCircumCenter * polygon[n].y + (1.0 - m_weightCircumCenter) * centerOfMass.y;
    }

    // The circumcenter is included in the face, then return the calculated circumcentre
    if (IsPointInPolygonNodes(result, polygon, m_projection))
    {
        return result;
    }

    // If the circumcenter is not included in the face,
    // the circumcenter will be placed at the intersection between an edge and the segment connecting the mass center with the circumcentre.
    for (UInt n = 0; n < numNodes; n++)
    {
        const auto nextNode = NextCircularForwardIndex(n, numNodes);

        const auto [areLineCrossing,
                    intersection,
                    crossProduct,
                    firstRatio,
                    secondRatio] = AreSegmentsCrossing(centerOfMass, result, polygon[n], polygon[nextNode], false, m_projection);

        if (areLineCrossing)
        {
            result = intersection;
            break;
        }
    }

    return result;
}

std::vector<meshkernel::Point> Mesh2D::GetObtuseTrianglesCenters()
{
    Administrate();
    std::vector<Point> result;
    result.reserve(GetNumFaces());
    for (UInt f = 0; f < GetNumFaces(); ++f)
    {
        // a triangle
        if (m_numFacesNodes[f] == 3)
        {
            const auto firstNode = m_facesNodes[f][0];
            const auto secondNode = m_facesNodes[f][1];
            const auto thirdNode = m_facesNodes[f][2];
            // compute squared edge lengths
            const auto firstEdgeSquaredLength = ComputeSquaredDistance(m_nodes[secondNode], m_nodes[firstNode], m_projection);
            const auto secondEdgeSquaredLength = ComputeSquaredDistance(m_nodes[thirdNode], m_nodes[firstNode], m_projection);
            const auto thirdEdgeSquaredLength = ComputeSquaredDistance(m_nodes[thirdNode], m_nodes[secondNode], m_projection);

            if (firstEdgeSquaredLength > secondEdgeSquaredLength + thirdEdgeSquaredLength ||
                secondEdgeSquaredLength > firstEdgeSquaredLength + thirdEdgeSquaredLength ||
                thirdEdgeSquaredLength > secondEdgeSquaredLength + firstEdgeSquaredLength)
            {
                result.emplace_back(m_facesMassCenters[f]);
            }
        }
    }
    return result;
}

std::vector<meshkernel::UInt> Mesh2D::GetEdgesCrossingSmallFlowEdges(double smallFlowEdgesThreshold)
{
    Administrate();
    std::vector<UInt> result;
    result.reserve(GetNumEdges());
    for (UInt e = 0; e < GetNumEdges(); ++e)
    {
        const auto firstFace = m_edgesFaces[e][0];
        const auto secondFace = m_edgesFaces[e][1];

        if (firstFace != constants::missing::uintValue && secondFace != constants::missing::uintValue)
        {
            const auto flowEdgeLength = ComputeDistance(m_facesCircumcenters[firstFace], m_facesCircumcenters[secondFace], m_projection);
            const double cutOffDistance = smallFlowEdgesThreshold * 0.5 * (std::sqrt(m_faceArea[firstFace]) + std::sqrt(m_faceArea[secondFace]));

            if (flowEdgeLength < cutOffDistance)
            {
                result.emplace_back(e);
            }
        }
    }
    return result;
}

std::vector<meshkernel::Point> Mesh2D::GetFlowEdgesCenters(const std::vector<UInt>& edges) const
{
    std::vector<Point> result;
    result.reserve(GetNumEdges());
    for (const auto& edge : edges)
    {
        const auto firstFace = m_edgesFaces[edge][0];
        const auto secondFace = m_edgesFaces[edge][1];
        result.emplace_back((m_facesCircumcenters[firstFace] + m_facesCircumcenters[secondFace]) * 0.5);
    }

    return result;
}

void Mesh2D::DeleteSmallFlowEdges(double smallFlowEdgesThreshold)
{
    DeleteDegeneratedTriangles();

    auto edges = GetEdgesCrossingSmallFlowEdges(smallFlowEdgesThreshold);
    if (!edges.empty())
    {
        // invalidate the edges
        for (const auto& e : edges)
        {
            m_edges[e] = {constants::missing::uintValue, constants::missing::uintValue};
        }
        Administrate();
    }
}

void Mesh2D::DeleteSmallTrianglesAtBoundaries(double minFractionalAreaTriangles)
{
    // On the second part, the small triangles at the boundaries are checked
    std::vector<std::vector<UInt>> smallTrianglesNodes;
    for (UInt face = 0; face < GetNumFaces(); ++face)
    {
        if (m_numFacesNodes[face] != m_numNodesInTriangle || m_faceArea[face] <= 0.0 || !IsFaceOnBoundary(face))
        {
            continue;
        }

        // compute the average area of neighboring faces
        double averageOtherFacesArea = 0.0;
        UInt numNonBoundaryFaces = 0;
        for (UInt e = 0; e < m_numNodesInTriangle; ++e)
        {
            // the edge must not be at the boundary, otherwise there is no "other" face
            const auto edge = m_facesEdges[face][e];
            if (IsEdgeOnBoundary(edge))
            {
                continue;
            }
            const auto otherFace = NextFace(face, edge);
            if (m_numFacesNodes[otherFace] > m_numNodesInTriangle)
            {
                averageOtherFacesArea += m_faceArea[otherFace];
                numNonBoundaryFaces++;
            }
        }

        if (numNonBoundaryFaces == 0 || m_faceArea[face] / (averageOtherFacesArea / static_cast<double>(numNonBoundaryFaces)) > minFractionalAreaTriangles)
        {
            // no valid boundary faces, the area of the current triangle is larger enough compared to the neighbors
            continue;
        }

        double minCosPhiSmallTriangle = 1.0;
        UInt nodeToPreserve = constants::missing::uintValue;
        UInt firstNodeToMerge = constants::missing::uintValue;
        UInt secondNodeToMerge = constants::missing::uintValue;
        UInt thirdEdgeSmallTriangle = constants::missing::uintValue;
        for (UInt e = 0; e < m_numNodesInTriangle; ++e)
        {
            const auto previousEdge = NextCircularBackwardIndex(e, m_numNodesInTriangle);
            const auto nextEdge = NextCircularForwardIndex(e, m_numNodesInTriangle);

            const auto k0 = m_facesNodes[face][previousEdge];
            const auto k1 = m_facesNodes[face][e];
            const auto k2 = m_facesNodes[face][nextEdge];

            // compute the angles between the edges
            const auto cosphi = std::abs(NormalizedInnerProductTwoSegments(m_nodes[k0], m_nodes[k1], m_nodes[k1], m_nodes[k2], m_projection));

            if (cosphi < minCosPhiSmallTriangle)
            {
                minCosPhiSmallTriangle = cosphi;
                nodeToPreserve = k1;
                firstNodeToMerge = k0;
                secondNodeToMerge = k2;
                thirdEdgeSmallTriangle = m_facesEdges[face][nextEdge];
            }
        }

        if (thirdEdgeSmallTriangle != constants::missing::uintValue && IsEdgeOnBoundary(thirdEdgeSmallTriangle))
        {
            smallTrianglesNodes.emplace_back(std::initializer_list<UInt>{nodeToPreserve, firstNodeToMerge, secondNodeToMerge});
        }
    }

    bool nodesMerged = false;
    for (const auto& triangleNodes : smallTrianglesNodes)
    {
        const auto nodeToPreserve = triangleNodes[0];
        const auto firstNodeToMerge = triangleNodes[1];
        const auto secondNodeToMerge = triangleNodes[2];

        // only
        UInt numInternalEdges = 0;
        for (UInt e = 0; e < m_nodesNumEdges[firstNodeToMerge]; ++e)
        {
            if (!IsEdgeOnBoundary(m_nodesEdges[firstNodeToMerge][e]))
            {
                numInternalEdges++;
            }
        }

        if (numInternalEdges == 1)
        {
            MergeTwoNodes(firstNodeToMerge, nodeToPreserve);
            nodesMerged = true;
        }

        // corner point of a triangle
        numInternalEdges = 0;
        for (UInt e = 0; e < m_nodesNumEdges[secondNodeToMerge]; ++e)
        {
            if (!IsEdgeOnBoundary(m_nodesEdges[secondNodeToMerge][e]))
            {
                numInternalEdges++;
            }
        }

        if (numInternalEdges == 1)
        {
            MergeTwoNodes(secondNodeToMerge, nodeToPreserve);
            nodesMerged = true;
        }
    }

    if (nodesMerged)
    {
        Administrate();
    }
}

void Mesh2D::ComputeNodeNeighbours()
{
    m_maxNumNeighbours = *(std::max_element(m_nodesNumEdges.begin(), m_nodesNumEdges.end()));
    m_maxNumNeighbours += 1;

    ResizeAndFill2DVector(m_nodesNodes, GetNumNodes(), m_maxNumNeighbours, true, constants::missing::uintValue);
    // for each node, determine the neighboring nodes
    for (UInt n = 0; n < GetNumNodes(); n++)
    {
        for (UInt nn = 0; nn < m_nodesNumEdges[n]; nn++)
        {
            const auto edge = m_edges[m_nodesEdges[n][nn]];
            m_nodesNodes[n][nn] = OtherNodeOfEdge(edge, n);
        }
    }
}

std::vector<double> Mesh2D::GetOrthogonality()
{
    std::vector<double> result(GetNumEdges());
    const auto numEdges = GetNumEdges();
    for (UInt e = 0; e < numEdges; e++)
    {
        auto val = constants::missing::doubleValue;
        const auto firstNode = m_edges[e].first;
        const auto secondNode = m_edges[e].second;
        const auto firstFaceIndex = m_edgesFaces[e][0];
        const auto secondFaceIndex = m_edgesFaces[e][1];

        if (firstNode != constants::missing::uintValue &&
            secondNode != constants::missing::uintValue &&
            firstFaceIndex != constants::missing::uintValue &&
            secondFaceIndex != constants::missing::uintValue && !IsEdgeOnBoundary(e))
        {
            val = NormalizedInnerProductTwoSegments(m_nodes[firstNode],
                                                    m_nodes[secondNode],
                                                    m_facesCircumcenters[firstFaceIndex],
                                                    m_facesCircumcenters[secondFaceIndex],
                                                    m_projection);

            if (val != constants::missing::doubleValue)
            {
                val = std::abs(val);
            }
        }
        result[e] = val;
    }
    return result;
}

std::vector<double> Mesh2D::GetSmoothness()
{
    std::vector<double> result(GetNumEdges());
    const auto numEdges = GetNumEdges();
    for (UInt e = 0; e < numEdges; e++)
    {
        auto val = constants::missing::doubleValue;
        const auto firstNode = m_edges[e].first;
        const auto secondNode = m_edges[e].second;
        const auto firstFaceIndex = m_edgesFaces[e][0];
        const auto secondFaceIndex = m_edgesFaces[e][1];

        if (firstNode != constants::missing::uintValue &&
            secondNode != constants::missing::uintValue &&
            firstFaceIndex != constants::missing::uintValue &&
            secondFaceIndex != constants::missing::uintValue && !IsEdgeOnBoundary(e))
        {
            const auto leftFaceArea = m_faceArea[firstFaceIndex];
            const auto rightFaceArea = m_faceArea[secondFaceIndex];

            if (leftFaceArea < m_minimumCellArea || rightFaceArea < m_minimumCellArea)
            {
                val = rightFaceArea / leftFaceArea;
            }
            if (val < 1.0)
            {
                val = 1.0 / val;
            }
        }
        result[e] = val;
    }
    return result;
}

void Mesh2D::ComputeAspectRatios(std::vector<double>& aspectRatios)
{
    std::vector<std::vector<double>> averageEdgesLength(GetNumEdges(), std::vector<double>(2, constants::missing::doubleValue));
    std::vector<double> averageFlowEdgesLength(GetNumEdges(), constants::missing::doubleValue);
    std::vector<bool> curvilinearGridIndicator(GetNumNodes(), true);
    std::vector<double> edgesLength(GetNumEdges(), 0.0);
    aspectRatios.resize(GetNumEdges(), 0.0);

    for (UInt e = 0; e < GetNumEdges(); e++)
    {
        const auto first = m_edges[e].first;
        const auto second = m_edges[e].second;

        if (first == second)
            continue;
        const double edgeLength = ComputeDistance(m_nodes[first], m_nodes[second], m_projection);
        edgesLength[e] = edgeLength;

        Point leftCenter;
        Point rightCenter;
        if (m_edgesNumFaces[e] > 0)
        {
            leftCenter = m_facesCircumcenters[m_edgesFaces[e][0]];
        }
        else
        {
            leftCenter = m_nodes[first];
        }

        // find right cell center, if it exists
        if (m_edgesNumFaces[e] == 2)
        {
            rightCenter = m_facesCircumcenters[m_edgesFaces[e][1]];
        }
        else
        {
            // otherwise, make ghost node by imposing boundary condition
            double dinry = InnerProductTwoSegments(m_nodes[first], m_nodes[second], m_nodes[first], leftCenter, m_projection);
            dinry = dinry / std::max(edgeLength * edgeLength, m_minimumEdgeLength);

            const double x0_bc = (1.0 - dinry) * m_nodes[first].x + dinry * m_nodes[second].x;
            const double y0_bc = (1.0 - dinry) * m_nodes[first].y + dinry * m_nodes[second].y;
            rightCenter.x = 2.0 * x0_bc - leftCenter.x;
            rightCenter.y = 2.0 * y0_bc - leftCenter.y;
        }

        averageFlowEdgesLength[e] = ComputeDistance(leftCenter, rightCenter, m_projection);
    }

    // Compute normal length
    for (UInt f = 0; f < GetNumFaces(); f++)
    {
        const auto numberOfFaceNodes = GetNumFaceEdges(f);
        if (numberOfFaceNodes < m_numNodesInTriangle)
            continue;

        for (UInt n = 0; n < numberOfFaceNodes; n++)
        {
            if (numberOfFaceNodes != m_numNodesQuads)
                curvilinearGridIndicator[m_facesNodes[f][n]] = false;
            const auto edgeIndex = m_facesEdges[f][n];

            if (m_edgesNumFaces[edgeIndex] == 0)
            {
                continue;
            }

            double edgeLength = edgesLength[edgeIndex];
            if (edgeLength != 0.0)
            {
                aspectRatios[edgeIndex] = averageFlowEdgesLength[edgeIndex] / edgeLength;
            }

            // quads
            if (numberOfFaceNodes == m_numNodesQuads)
            {
                UInt kkp2 = n + 2;
                if (kkp2 >= numberOfFaceNodes)
                {
                    kkp2 = kkp2 - numberOfFaceNodes;
                }

                const auto klinkp2 = m_facesEdges[f][kkp2];
                edgeLength = 0.5 * (edgesLength[edgeIndex] + edgesLength[klinkp2]);
            }

            if (IsEqual(averageEdgesLength[edgeIndex][0], constants::missing::doubleValue))
            {
                averageEdgesLength[edgeIndex][0] = edgeLength;
            }
            else
            {
                averageEdgesLength[edgeIndex][1] = edgeLength;
            }
        }
    }

    if (m_curvilinearToOrthogonalRatio == 1.0)
        return;

    for (UInt e = 0; e < GetNumEdges(); e++)
    {
        const auto first = m_edges[e].first;
        const auto second = m_edges[e].second;

        if (first == constants::missing::uintValue || second == constants::missing::uintValue)
            continue;
        if (m_edgesNumFaces[e] == 0)
            continue;
        // Consider only quads
        if (!curvilinearGridIndicator[first] || !curvilinearGridIndicator[second])
            continue;

        if (IsEdgeOnBoundary(e))
        {
            if (averageEdgesLength[e][0] > 0.0 &&
                IsEqual(averageEdgesLength[e][0], constants::missing::doubleValue))
            {
                aspectRatios[e] = averageFlowEdgesLength[e] / averageEdgesLength[e][0];
            }
        }
        else
        {
            if (averageEdgesLength[e][0] > 0.0 &&
                averageEdgesLength[e][1] > 0.0 &&
                IsEqual(averageEdgesLength[e][0], constants::missing::doubleValue) &&
                IsEqual(averageEdgesLength[e][1], constants::missing::doubleValue))
            {
                aspectRatios[e] = m_curvilinearToOrthogonalRatio * aspectRatios[e] +
                                  (1.0 - m_curvilinearToOrthogonalRatio) * averageFlowEdgesLength[e] / (0.5 * (averageEdgesLength[e][0] + averageEdgesLength[e][1]));
            }
        }
    }
}

void Mesh2D::TriangulateFaces()
{
    for (UInt i = 0; i < GetNumFaces(); ++i)
    {
        const auto NumEdges = GetNumFaceEdges(i);

        if (NumEdges < 4)
        {
            continue;
        }

        const auto indexFirstNode = m_facesNodes[i][0];
        for (UInt j = 2; j < NumEdges - 1; j++)
        {
            const auto nodeIndex = m_facesNodes[i][j];
            ConnectNodes(indexFirstNode, nodeIndex);
        }
    }

    m_edgesRTreeRequiresUpdate = true;
}

void Mesh2D::MakeDualFace(UInt node, double enlargementFactor, std::vector<Point>& dualFace)
{
    const auto sortedFacesIndices = SortedFacesAroundNode(node);
    const auto numEdges = m_nodesNumEdges[node];
    dualFace.reserve(m_maximumNumberOfEdgesPerNode);
    dualFace.clear();

    for (UInt e = 0; e < numEdges; ++e)
    {
        const auto edgeIndex = m_nodesEdges[node][e];
        auto edgeCenter = m_edgesCenters[edgeIndex];

        if (m_projection == Projection::spherical)
        {
            const auto firstNodeIndex = m_edges[edgeIndex].first;
            const auto secondNodeIndex = m_edges[edgeIndex].second;

            if (firstNodeIndex != constants::missing::uintValue && secondNodeIndex != constants::missing::uintValue)
            {
                const auto diff = m_nodes[firstNodeIndex].x - m_nodes[secondNodeIndex].x;

                if (diff > 180.0)
                {
                    edgeCenter.x = edgeCenter.x - 180.0;
                }
                if (diff < -180.0)
                {
                    edgeCenter.x = edgeCenter.x + 180.0;
                }
            }
        }
        dualFace.emplace_back(edgeCenter);

        const auto faceIndex = sortedFacesIndices[e];
        if (faceIndex != constants::missing::uintValue)
        {
            dualFace.emplace_back(m_facesMassCenters[faceIndex]);
        }
        else
        {
            dualFace.emplace_back(m_nodes[node]);
        }
    }
    dualFace.emplace_back(dualFace[0]);

    // now we can compute the mass center of the dual face
    auto [area, centerOfMass, direction] = Polygon::FaceAreaAndCenterOfMass(dualFace, m_projection);

    if (m_projection == Projection::spherical)
    {
        if (centerOfMass.x - m_nodes[node].x > 180.0)
        {
            centerOfMass.x -= 360.0;
        }
        if (centerOfMass.x - m_nodes[node].x < -180.0)
        {
            centerOfMass.x += 360.0;
        }
    }

    for (auto& v : dualFace)
    {
        v = centerOfMass + (v - centerOfMass) * enlargementFactor;
    }
}

std::vector<meshkernel::UInt> Mesh2D::SortedFacesAroundNode(UInt node) const
{

    const auto numEdges = m_nodesNumEdges[node];
    std::vector<UInt> result;
    for (UInt e = 0; e < numEdges; ++e)
    {
        const auto firstEdge = m_nodesEdges[node][e];

        // no faces for this edge
        if (m_edgesNumFaces[firstEdge] == 0)
        {
            continue;
        }

        auto const ee = NextCircularForwardIndex(e, numEdges);
        const auto secondEdge = m_nodesEdges[node][ee];
        const auto firstFace = m_edgesFaces[firstEdge][0];

        UInt secondFace = constants::missing::uintValue;
        if (m_edgesNumFaces[firstEdge] > 1)
        {
            secondFace = m_edgesFaces[firstEdge][1];
        }

        // check if the first face contains the first edge
        UInt firstEdgeIndexInFirstFace = 0;
        for (UInt n = 0; n < m_numFacesNodes[firstFace]; ++n)
        {
            if (m_facesEdges[firstFace][n] == firstEdge)
            {
                firstEdgeIndexInFirstFace = n;
                break;
            }
        }

        // check if previous edge in firstFace is secondEdge (so at least two edges share the same edge)
        auto const secondEdgeindexInFirstFace = NextCircularBackwardIndex(firstEdgeIndexInFirstFace, m_numFacesNodes[firstFace]);

        if (m_facesEdges[firstFace][secondEdgeindexInFirstFace] == secondEdge)
        {
            result.emplace_back(firstFace);
        }
        else
        {
            result.emplace_back(secondFace);
        }
    }

    return result;
}

std::vector<meshkernel::Point> Mesh2D::MeshBoundaryToPolygon(const std::vector<Point>& polygonNodes)
{

    Polygon polygon(polygonNodes, m_projection);

    // Find faces
    Administrate();
    std::vector<bool> isVisited(GetNumEdges(), false);
    std::vector<Point> meshBoundaryPolygon;
    meshBoundaryPolygon.reserve(GetNumNodes());

    for (UInt e = 0; e < GetNumEdges(); e++)
    {
        if (isVisited[e] || !IsEdgeOnBoundary(e))
        {
            continue;
        }

        const auto firstNodeIndex = m_edges[e].first;
        const auto secondNodeIndex = m_edges[e].second;
        const auto firstNode = m_nodes[firstNodeIndex];
        const auto secondNode = m_nodes[secondNodeIndex];

        bool firstNodeInPolygon = polygon.Contains(m_nodes[firstNodeIndex]);
        bool secondNodeInPolygon = polygon.Contains(m_nodes[secondNodeIndex]);

        if (!firstNodeInPolygon && !secondNodeInPolygon)
        {
            continue;
        }

        // Start a new polyline
        if (!meshBoundaryPolygon.empty())
        {
            meshBoundaryPolygon.emplace_back(constants::missing::doubleValue, constants::missing::doubleValue);
        }

        // Put the current edge on the mesh boundary, mark it as visited
        const auto startPolygonEdges = static_cast<UInt>(meshBoundaryPolygon.size());
        meshBoundaryPolygon.emplace_back(firstNode);
        meshBoundaryPolygon.emplace_back(secondNode);
        isVisited[e] = true;

        // walk the current mesh boundary
        auto currentNode = secondNodeIndex;
        WalkBoundaryFromNode(polygon, isVisited, currentNode, meshBoundaryPolygon);

        const auto numNodesFirstTail = static_cast<UInt>(meshBoundaryPolygon.size());

        // if the boundary polygon is not closed
        if (currentNode != firstNodeIndex)
        {
            // Now grow a polyline starting at the other side of the original link L, i.e., the second tail
            currentNode = firstNodeIndex;
            WalkBoundaryFromNode(polygon, isVisited, currentNode, meshBoundaryPolygon);
        }

        // There is a nonempty second tail, so reverse the first tail, so that they connect.
        if (meshBoundaryPolygon.size() > numNodesFirstTail)
        {
            const auto start = startPolygonEdges + static_cast<UInt>(std::ceil((numNodesFirstTail - startPolygonEdges + static_cast<UInt>(1)) * 0.5));
            for (auto n = start; n < numNodesFirstTail; n++)
            {
                const auto backupPoint = meshBoundaryPolygon[n];
                const auto replaceIndex = numNodesFirstTail - n + firstNodeIndex;
                meshBoundaryPolygon[n] = meshBoundaryPolygon[replaceIndex];
                meshBoundaryPolygon[replaceIndex] = backupPoint;
            }
        }
    }
    return meshBoundaryPolygon;
}

void Mesh2D::WalkBoundaryFromNode(const Polygon& polygon,
                                  std::vector<bool>& isVisited,
                                  UInt& currentNode,
                                  std::vector<Point>& meshBoundaryPolygon) const
{
    UInt e = 0;
    bool currentNodeInPolygon = false;
    while (e < m_nodesNumEdges[currentNode])
    {
        if (!currentNodeInPolygon)
        {
            currentNodeInPolygon = polygon.Contains(m_nodes[currentNode]);
        }

        if (!currentNodeInPolygon)
        {
            break;
        }

        const auto currentEdge = m_nodesEdges[currentNode][e];
        if (isVisited[currentEdge] || !IsEdgeOnBoundary(currentEdge))
        {
            e++;
            continue;
        }

        currentNode = OtherNodeOfEdge(m_edges[currentEdge], currentNode);
        e = 0;
        currentNodeInPolygon = false;

        meshBoundaryPolygon.emplace_back(m_nodes[currentNode]);
        isVisited[currentEdge] = true;
    }
}

std::vector<meshkernel::UInt> Mesh2D::GetHangingEdges() const
{
    std::vector<UInt> result;
    for (UInt e = 0; e < GetNumEdges(); e++)
    {
        const auto firstNode = m_edges[e].first;
        const auto secondNode = m_edges[e].second;

        if (firstNode != constants::missing::uintValue && secondNode != constants::missing::uintValue)
        {
            // if one of the nodes has no other attached edges, the current edge is an hanging edge
            if (m_nodesNumEdges[firstNode] > 1 && m_nodesNumEdges[secondNode] > 1)
            {
                continue;
            }
            result.emplace_back(e);
        }
    }

    return result;
}

void Mesh2D::DeleteMesh(const Polygons& polygon, int deletionOption, bool invertDeletion)
{
    // Find crossed faces
    const auto [edgeIntersections, faceIntersections] = GetPolygonIntersections(polygon);

    // Find faces with all nodes inside the polygon
    std::vector<bool> isNodeInsidePolygon(GetNumNodes(), false);
    for (UInt n = 0; n < GetNumNodes(); ++n)
    {
        auto [isInPolygon, polygonIndex] = polygon.IsPointInPolygons(m_nodes[n]);
        if (isInPolygon)
        {
            isNodeInsidePolygon[n] = true;
        }
    }

    std::vector<bool> isFaceCompletlyIncludedInPolygon(GetNumFaces(), true);
    for (UInt f = 0; f < GetNumFaces(); ++f)
    {
        for (UInt n = 0; n < GetNumFaceEdges(f); ++n)
        {
            const auto nodeIndex = m_facesNodes[f][n];
            if (!isNodeInsidePolygon[nodeIndex])
            {
                isFaceCompletlyIncludedInPolygon[f] = false;
                break;
            }
        }
    }

    std::function<bool(UInt)> excludedFace;
    if (deletionOption == InsideNotIntersected && !invertDeletion)
    {
        excludedFace = [&](UInt f)
        { return !isFaceCompletlyIncludedInPolygon[f] || faceIntersections[f].faceIndex != constants::missing::uintValue; };
    }
    else if (deletionOption == InsideNotIntersected && invertDeletion)
    {
        excludedFace = [&](UInt f)
        { return isFaceCompletlyIncludedInPolygon[f] && faceIntersections[f].faceIndex == constants::missing::uintValue; };
    }
    else if (deletionOption == InsideAndIntersected && !invertDeletion)
    {
        excludedFace = [&](UInt f)
        { return !isFaceCompletlyIncludedInPolygon[f] && faceIntersections[f].faceIndex == constants::missing::uintValue; };
    }
    else if (deletionOption == InsideAndIntersected && invertDeletion)
    {
        excludedFace = [&](UInt f)
        { return isFaceCompletlyIncludedInPolygon[f] || faceIntersections[f].faceIndex != constants::missing::uintValue; };
    }

    // Mark edges for deletion
    for (UInt e = 0; e < GetNumEdges(); ++e)
    {
        const auto numEdgeFaces = GetNumEdgesFaces(e);
        bool deleteEdge = true;

        if (numEdgeFaces == 1 && excludedFace(m_edgesFaces[e][0]))
        {
            deleteEdge = false;
        }
        if (numEdgeFaces == 2 && (excludedFace(m_edgesFaces[e][0]) || excludedFace(m_edgesFaces[e][1])))
        {
            deleteEdge = false;
        }

        if (deleteEdge)
        {
            m_edges[e].first = constants::missing::uintValue;
            m_edges[e].second = constants::missing::uintValue;
        }
    }

    m_nodesRTreeRequiresUpdate = true;
    m_edgesRTreeRequiresUpdate = true;

    AdministrateNodesEdges();
}

void Mesh2D::DeleteHangingEdges()
{
    const auto hangingEdges = GetHangingEdges();
    for (const auto& hangingEdge : hangingEdges)
    {
        DeleteEdge(hangingEdge);
    }
}

std::vector<meshkernel::UInt> Mesh2D::PointFaceIndices(const std::vector<Point>& points)
{
    const auto numPoints = static_cast<UInt>(points.size());
    std::vector<UInt> result;
    result.resize(numPoints, constants::missing::uintValue);
    std::vector<Point> polygonNodesCache;
    BuildTree(Location::Edges);

    for (UInt i = 0; i < numPoints; ++i)
    {
        const auto edgeIndex = FindEdgeCloseToAPoint(points[i]);

        if (edgeIndex == constants::missing::uintValue)
        {
            result[i] = constants::missing::uintValue;
            continue;
        }

        for (UInt e = 0; e < m_edgesNumFaces[edgeIndex]; ++e)
        {
            const auto faceIndex = m_edgesFaces[edgeIndex][e];
            ComputeFaceClosedPolygon(faceIndex, polygonNodesCache);
            const auto isPointInFace = IsPointInPolygonNodes(points[i], polygonNodesCache, m_projection);
            if (isPointInFace)
            {
                result[i] = faceIndex;
                break;
            }
        }
    }
    return result;
}

std::tuple<meshkernel::UInt, meshkernel::UInt> Mesh2D::IsSegmentCrossingABoundaryEdge(const Point& firstPoint,
                                                                                      const Point& secondPoint) const
{
    double intersectionRatio = std::numeric_limits<double>::max();
    UInt intersectedFace = constants::missing::uintValue;
    UInt intersectedEdge = constants::missing::uintValue;
    for (UInt e = 0; e < GetNumEdges(); ++e)
    {
        if (!IsEdgeOnBoundary(e))
        {
            continue;
        }

        const auto [areSegmentCrossing,
                    intersectionPoint,
                    crossProduct,
                    ratioFirstSegment,
                    ratioSecondSegment] = AreSegmentsCrossing(firstPoint,
                                                              secondPoint,
                                                              m_nodes[m_edges[e].first],
                                                              m_nodes[m_edges[e].second],
                                                              false,
                                                              m_projection);

        if (areSegmentCrossing && ratioFirstSegment < intersectionRatio)
        {
            intersectionRatio = ratioFirstSegment;
            intersectedFace = m_edgesFaces[e][0];
            intersectedEdge = e;
        }
    }

    return {intersectedFace, intersectedEdge};
}

std::tuple<meshkernel::UInt, meshkernel::UInt> Mesh2D::GetIntersectionSeed(const std::vector<Point>& polyLine,
                                                                           const std::vector<bool>& vistedEdges) const
{
    UInt crossedSegmentIndex = constants::missing::uintValue;
    UInt crossedEdgeIndex = constants::missing::uintValue;

    // Find starting edge and segment
    for (UInt segmentIndex = 0; segmentIndex < polyLine.size() - 1; ++segmentIndex)
    {
        for (UInt edgeIndex = 0; edgeIndex < GetNumEdges(); ++edgeIndex)
        {
            // edge already crossed, nothing to do
            if (vistedEdges[edgeIndex])
            {
                continue;
            }

            const auto [isEdgeCrossed,
                        intersectionPoint,
                        crossProductValue,
                        adimensionalPolylineSegmentDistance,
                        adimensionalEdgeDistance] = AreSegmentsCrossing(polyLine[segmentIndex],
                                                                        polyLine[segmentIndex + 1],
                                                                        m_nodes[m_edges[edgeIndex].first],
                                                                        m_nodes[m_edges[edgeIndex].second],
                                                                        false,
                                                                        m_projection);
            if (!isEdgeCrossed)
            {
                continue;
            }

            crossedSegmentIndex = segmentIndex;
            crossedEdgeIndex = edgeIndex;
            break;
        }

        if (crossedSegmentIndex != constants::missing::uintValue)
        {
            break;
        }
    }

    return {crossedEdgeIndex, crossedSegmentIndex};
}

void Mesh2D::GetPolylineIntersection(const std::vector<Point>& polyLine,
                                     std::vector<EdgeMeshPolylineIntersection>& edgesIntersectionsCache,
                                     std::vector<FaceMeshPolylineIntersection>& facesIntersectionsCache,
                                     std::vector<EdgeMeshPolylineIntersection>& edgesIntersectionsResult,
                                     std::vector<FaceMeshPolylineIntersection>& faceIntersectionsResult) const
{
    // 1. Find the intersection of any segment of the polyline with the mesh return if nothing is found
    std::ranges::fill(edgesIntersectionsCache, EdgeMeshPolylineIntersection());
    std::ranges::fill(facesIntersectionsCache, FaceMeshPolylineIntersection());

    const auto polylineSize = static_cast<UInt>(polyLine.size());

    std::vector<double> cumulativeLength(polyLine.size(), 0.0);
    for (UInt i = 1; i < polylineSize; ++i)
    {
        cumulativeLength[i] = cumulativeLength[i - 1] + ComputeDistance(polyLine[i], polyLine[i - 1], m_projection);
    }

    std::queue<std::array<UInt, 3>> crossingEdges;
    std::vector<bool> vistedEdges(GetNumEdges(), false);
    std::vector<bool> vistedFace(GetNumEdges(), false);

    // keep traversing the polyline as long crossed edges are found
    while (true)
    {
        // find a crossed edge on a non-visited segment
        const auto intersectionSeed = GetIntersectionSeed(polyLine,
                                                          vistedEdges);

        const auto crossedEdgeIndex = std::get<0>(intersectionSeed);

        // no valid seed found in the entire mesh, we are done
        if (crossedEdgeIndex == constants::missing::uintValue)
        {
            break;
        }

        const auto crossingSegmentIndex = std::get<1>(intersectionSeed);
        const auto crossingNextSegmentIndex = crossingSegmentIndex + 1;
        crossingEdges.push({crossedEdgeIndex, crossingSegmentIndex, crossingNextSegmentIndex});

        // use breadth search along the current polyline
        while (!crossingEdges.empty())
        {
            auto [currentCrossingEdge, segmentIndex, nextSegmentIndex] = crossingEdges.front();
            crossingEdges.pop();

            for (UInt f = 0; f < m_edgesFaces[currentCrossingEdge].size(); ++f)
            {
                const auto currentFaceIndex = m_edgesFaces[currentCrossingEdge][f];

                if (currentFaceIndex == constants::missing::uintValue)
                {
                    continue;
                }

                for (UInt e = 0; e < m_facesEdges[currentFaceIndex].size(); ++e)
                {
                    const auto edgeIndex = m_facesEdges[currentFaceIndex][e];
                    if (vistedEdges[edgeIndex] && vistedFace[currentFaceIndex])
                    {
                        continue;
                    }

                    UInt firstIndex = segmentIndex;
                    UInt secondIndex = nextSegmentIndex;

                    auto intersection = AreSegmentsCrossing(polyLine[firstIndex],
                                                            polyLine[secondIndex],
                                                            m_nodes[m_edges[edgeIndex].first],
                                                            m_nodes[m_edges[edgeIndex].second],
                                                            false,
                                                            m_projection);

                    auto intersectionFound = std::get<0>(intersection);
                    if (!intersectionFound)
                    {
                        UInt numForwardSteps = 0;
                        while (!intersectionFound && firstIndex < polylineSize - 2 && numForwardSteps < m_maxSteps)
                        {
                            firstIndex = secondIndex;
                            secondIndex = firstIndex + 1;
                            intersection = AreSegmentsCrossing(polyLine[firstIndex],
                                                               polyLine[secondIndex],
                                                               m_nodes[m_edges[edgeIndex].first],
                                                               m_nodes[m_edges[edgeIndex].second],
                                                               false,
                                                               m_projection);

                            intersectionFound = std::get<0>(intersection);
                            numForwardSteps++;
                        }
                    }

                    if (!intersectionFound)
                    {
                        firstIndex = segmentIndex;
                        secondIndex = nextSegmentIndex;
                        UInt numBackwardSteps = 0;
                        while (!intersectionFound && firstIndex >= 1 && numBackwardSteps < m_maxSteps)
                        {
                            secondIndex = firstIndex;
                            firstIndex = firstIndex - 1;
                            intersection = AreSegmentsCrossing(polyLine[firstIndex],
                                                               polyLine[secondIndex],
                                                               m_nodes[m_edges[edgeIndex].first],
                                                               m_nodes[m_edges[edgeIndex].second],
                                                               false,
                                                               m_projection);

                            intersectionFound = std::get<0>(intersection);
                            numBackwardSteps++;
                        }
                    }

                    // none of the polyline intersect the current edge
                    if (!intersectionFound)
                    {
                        continue;
                    }

                    const double crossProductValue = std::get<2>(intersection);
                    const double adimensionalPolylineSegmentDistance = std::get<3>(intersection);
                    const double adimensionalEdgeDistance = std::get<4>(intersection);

                    EdgeMeshPolylineIntersection::updateIntersections(
                        firstIndex,
                        edgeIndex,
                        m_edges[edgeIndex].first,
                        m_edges[edgeIndex].second,
                        cumulativeLength,
                        crossProductValue,
                        adimensionalEdgeDistance,
                        adimensionalPolylineSegmentDistance,
                        edgesIntersectionsCache);

                    FaceMeshPolylineIntersection::updateIntersections(currentFaceIndex, edgeIndex, facesIntersectionsCache);

                    if (edgeIndex != currentCrossingEdge)
                    {
                        crossingEdges.push({edgeIndex, firstIndex, secondIndex});
                    }
                    vistedEdges[edgeIndex] = true;
                }
                vistedFace[currentFaceIndex] = true;
            }
        }
    }

    // compute polylineDistance, sort the edges for each face
    for (auto& facesIntersection : facesIntersectionsCache)
    {
        if (facesIntersection.edgeIndexses.empty())
        {
            continue;
        }

        if (facesIntersection.edgeIndexses.size() == 1)
        {
            const auto edgeIndex = facesIntersection.edgeIndexses[0];
            facesIntersection.polylineDistance = edgesIntersectionsCache[edgeIndex].polylineDistance;
        }

        if (facesIntersection.edgeIndexses.size() == 2)
        {
            const auto firstEdgeIndex = facesIntersection.edgeIndexses[0];
            const auto secondEdgeIndex = facesIntersection.edgeIndexses[1];

            // swap the edge indexes if needed
            if (edgesIntersectionsCache[firstEdgeIndex].adimensionalPolylineSegmentDistance > edgesIntersectionsCache[secondEdgeIndex].adimensionalPolylineSegmentDistance)
            {
                std::swap(facesIntersection.edgeIndexses[0], facesIntersection.edgeIndexses[1]);
            }

            // compute the polylineDistance for the face
            facesIntersection.polylineDistance = 0.5 * (edgesIntersectionsCache[firstEdgeIndex].polylineDistance + edgesIntersectionsCache[secondEdgeIndex].polylineDistance);
        }

        // push back the face intersection edge nodes
        for (UInt e = 0; e < facesIntersection.edgeIndexses.size(); ++e)
        {
            const auto edgeIndex = facesIntersection.edgeIndexses[e];
            facesIntersection.edgeNodes.emplace_back(edgesIntersectionsCache[edgeIndex].edgeFirstNode);
            facesIntersection.edgeNodes.emplace_back(edgesIntersectionsCache[edgeIndex].edgeSecondNode);
        }
    }

    // edge intersections are unique
    for (UInt e = 0; e < GetNumEdges(); ++e)
    {
        if (edgesIntersectionsResult[e].polylineDistance < 0)
        {
            edgesIntersectionsResult[e] = edgesIntersectionsCache[e];
        }
    }

    // face intersections are not unique and a face could have been intersected already
    for (UInt f = 0; f < GetNumFaces(); ++f)
    {
        if (!faceIntersectionsResult[f].edgeNodes.empty() &&
            !facesIntersectionsCache[f].edgeNodes.empty())
        {
            faceIntersectionsResult[f].edgeIndexses.insert(faceIntersectionsResult[f].edgeIndexses.end(), facesIntersectionsCache[f].edgeIndexses.begin(), facesIntersectionsCache[f].edgeIndexses.end());
            faceIntersectionsResult[f].edgeNodes.insert(faceIntersectionsResult[f].edgeNodes.end(), facesIntersectionsCache[f].edgeNodes.begin(), facesIntersectionsCache[f].edgeNodes.end());
            faceIntersectionsResult[f].polylineDistance = 0.5 * (faceIntersectionsResult[f].polylineDistance + facesIntersectionsCache[f].polylineDistance);
        }
        else if (!facesIntersectionsCache[f].edgeNodes.empty())
        {
            faceIntersectionsResult[f] = facesIntersectionsCache[f];
        }
    }
}

std::tuple<std::vector<meshkernel::EdgeMeshPolylineIntersection>, std::vector<meshkernel::FaceMeshPolylineIntersection>>
Mesh2D::GetPolygonIntersections(const Polygons& polygon)
{

    // Intersection results
    std::vector<EdgeMeshPolylineIntersection> edgesIntersectionsResult(GetNumEdges());
    std::vector<FaceMeshPolylineIntersection> faceIntersectionsResult(GetNumFaces());

    // No polygon, nothing is crossed
    if (polygon.IsEmpty())
    {
        return {edgesIntersectionsResult, faceIntersectionsResult};
    }

    // Declare local caches
    std::vector<EdgeMeshPolylineIntersection> edgesIntersectionsCache(GetNumEdges());
    std::vector<FaceMeshPolylineIntersection> facesIntersectionsCache(GetNumFaces());

    // Make sure face information is available
    Administrate();

    // Multiple polygons here
    for (auto outer = 0u; outer < polygon.GetNumPolygons(); ++outer)
    {
        std::vector<std::vector<Point>> allPolylines;

        allPolylines.emplace_back(polygon.Enclosure(outer).Outer().Nodes());

        for (auto inner = 0u; inner < polygon.Enclosure(outer).NumberOfInner(); ++inner)
        {
            allPolylines.emplace_back(polygon.Enclosure(outer).Inner(inner).Nodes());
        }

        for (const auto& polyLine : allPolylines)
        {
            GetPolylineIntersection(polyLine,
                                    edgesIntersectionsCache,
                                    facesIntersectionsCache,
                                    edgesIntersectionsResult,
                                    faceIntersectionsResult);
        }
    }

    return {edgesIntersectionsResult, faceIntersectionsResult};
}

std::vector<int> Mesh2D::MaskEdgesOfFacesInPolygon(const Polygons& polygons, bool invertSelection, bool includeIntersected) const
{
    // mark all nodes in polygon with 1
    std::vector<int> nodeMask(GetNumNodes(), 0);
    for (UInt n = 0; n < GetNumNodes(); ++n)
    {
        const auto [isInPolygon, polygonIndex] = polygons.IsPointInPolygons(m_nodes[n]);
        if (isInPolygon)
        {
            nodeMask[n] = 1;
        }
    }

    // mark all edges with both start end end nodes included with 1
    std::vector<int> edgeMask(m_edges.size(), 0);
    for (UInt e = 0; e < GetNumEdges(); ++e)
    {
        const auto firstNodeIndex = m_edges[e].first;
        const auto secondNodeIndex = m_edges[e].second;

        int isEdgeIncluded;
        if (includeIntersected)
        {
            isEdgeIncluded = ((firstNodeIndex != constants::missing::uintValue && nodeMask[firstNodeIndex] == 1) ||
                              (secondNodeIndex != constants::missing::uintValue && nodeMask[secondNodeIndex] == 1))
                                 ? 1
                                 : 0;
        }
        else
        {
            isEdgeIncluded = (firstNodeIndex != constants::missing::uintValue && nodeMask[firstNodeIndex] == 1 &&
                              secondNodeIndex != constants::missing::uintValue && nodeMask[secondNodeIndex] == 1)
                                 ? 1
                                 : 0;
        }

        edgeMask[e] = isEdgeIncluded;
    }

    // if one edge of the face is not included do not include all the edges of that face
    auto secondEdgeMask = edgeMask;
    if (!includeIntersected)
    {
        for (UInt f = 0; f < GetNumFaces(); ++f)
        {
            bool isOneEdgeNotIncluded = false;
            for (UInt n = 0; n < GetNumFaceEdges(f); ++n)
            {
                const auto edgeIndex = m_facesEdges[f][n];
                if (edgeIndex != constants::missing::uintValue && edgeMask[edgeIndex] == 0)
                {
                    isOneEdgeNotIncluded = true;
                    break;
                }
            }

            if (isOneEdgeNotIncluded)
            {
                for (UInt n = 0; n < GetNumFaceEdges(f); ++n)
                {
                    const auto edgeIndex = m_facesEdges[f][n];
                    if (edgeIndex != constants::missing::uintValue)
                    {
                        secondEdgeMask[edgeIndex] = 0;
                    }
                }
            }
        }
    }

    // if the selection is inverted, do not delete the edges of included faces
    if (invertSelection)
    {
        for (UInt e = 0; e < GetNumEdges(); ++e)
        {
            if (secondEdgeMask[e] == 0)
            {
                secondEdgeMask[e] = 1;
            }

            if (edgeMask[e] == 1)
            {
                secondEdgeMask[e] = 0;
            }
        }
    }

    return secondEdgeMask;
}

std::vector<int> Mesh2D::NodeMaskFromEdgeMask(std::vector<int> const& edgeMask) const
{
    if (edgeMask.size() != GetNumEdges())
    {
        throw std::invalid_argument("Mesh2D::NodeMaskFromEdgeMask:The dimension of the edge mask do not fit the mesh.");
    }

    // fill node mask to 0
    std::vector<int> nodeMask(GetNumNodes(), 0);

    // compute node mask from edge mask
    for (UInt e = 0; e < GetNumEdges(); ++e)
    {
        if (edgeMask[e] != 1)
            continue;

        const auto firstNodeIndex = m_edges[e].first;
        const auto secondNodeIndex = m_edges[e].second;

        if (firstNodeIndex != constants::missing::uintValue)
        {
            nodeMask[firstNodeIndex] = 1;
        }
        if (secondNodeIndex != constants::missing::uintValue)
        {
            nodeMask[secondNodeIndex] = 1;
        }
    }
    return nodeMask;
}

std::vector<int> Mesh2D::NodeMaskFromPolygon(const Polygons& polygon, bool inside) const
{
    std::vector<int> nodeMask(GetNumNodes(), 0);
    const auto nodePolygonIndices = polygon.PointsInPolygons(m_nodes);

    for (UInt i = 0; i < nodeMask.size(); ++i)
    {
        auto isInPolygon = nodePolygonIndices[i];
        if (!inside)
        {
            isInPolygon = !isInPolygon;
        }
        nodeMask[i] = 0;
        if (isInPolygon)
        {
            nodeMask[i] = 1;
        }
    }
    return nodeMask;
}

meshkernel::UInt Mesh2D::FindOppositeEdge(const UInt faceId, const UInt edgeId) const
{
    if (m_numFacesNodes[faceId] != 4)
    {
        throw NotImplementedError("FindOppositeEdge only works for quadrilateral elements, request is for element with {} edges",
                                  m_numFacesNodes[faceId]);
    }

    UInt position = constants::missing::uintValue;

    // Find the corresponding position of edge
    for (UInt i = 0; i < m_numFacesNodes[faceId]; ++i)
    {
        if (m_facesEdges[faceId][i] == edgeId)
        {
            position = i;
            break;
        }
    }

    UInt opposite;

    switch (position)
    {
    case 0:
        opposite = 2;
        break;
    case 1:
        opposite = 3;
        break;
    case 2:
        opposite = 0;
        break;
    case 3:
        opposite = 1;
        break;
    default:
        opposite = constants::missing::uintValue;
    }

    if (opposite != constants::missing::uintValue)
    {
        return m_facesEdges[faceId][opposite];
    }

    return constants::missing::uintValue;
}

meshkernel::UInt Mesh2D::NextFace(const UInt faceId, const UInt edgeId) const
{
    if (faceId != constants::missing::uintValue)
    {
        if (m_edgesFaces[edgeId][0] == faceId)
        {
            return m_edgesFaces[edgeId][1];
        }

        if (m_edgesFaces[edgeId][1] == faceId)
        {
            return m_edgesFaces[edgeId][0];
        }
    }

    return constants::missing::uintValue;
}
