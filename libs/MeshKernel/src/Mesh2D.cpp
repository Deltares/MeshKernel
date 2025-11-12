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

#include <numeric>

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/Mesh2DIntersections.hpp"
#include "MeshKernel/MeshFaceCenters.hpp"
#include "MeshKernel/MeshOrthogonality.hpp"
#include "MeshKernel/Operations.hpp"
#include "MeshKernel/Polygon.hpp"
#include "MeshKernel/Polygons.hpp"
#include "MeshKernel/TriangulationWrapper.hpp"
#include "MeshKernel/UndoActions/CompoundUndoAction.hpp"
#include "MeshKernel/Utilities/RTreeBase.hpp"

using meshkernel::Mesh2D;

Mesh2D::Mesh2D()
    : Mesh()
{
}

Mesh2D::Mesh2D(Projection projection)
    : Mesh(projection)
{
}

Mesh2D::Mesh2D(const std::vector<Edge>& edges,
               const std::vector<Point>& nodes,
               Projection projection)
    : Mesh(edges, nodes, projection)
{
    DoAdministration();

    if (InvalidateEdgesWithNoFace() > 0)
    {
        DeleteInvalidNodesAndEdges();
        DoAdministration();
    }
}

Mesh2D::Mesh2D(const std::vector<Edge>& edges,
               const std::vector<Point>& nodes,
               const std::vector<std::vector<UInt>>& faceNodes,
               const std::vector<std::uint8_t>& numFaceNodes,
               Projection projection)
    : Mesh(edges, nodes, projection)
{
    DoAdministrationGivenFaceNodesMapping(faceNodes, numFaceNodes);
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

        const Point approximateCenter = (inputNodes[triangulationWrapper.GetFaceNode(i, 0)] +
                                         inputNodes[triangulationWrapper.GetFaceNode(i, 1)] +
                                         inputNodes[triangulationWrapper.GetFaceNode(i, 2)]) *
                                        constants::numeric::oneThird;

        const auto isTriangleInPolygon = polygons.IsPointInPolygon(approximateCenter, 0);
        if (!isTriangleInPolygon)
        {
            continue;
        }

        // mark all edges of this triangle as good ones
        for (UInt j = 0; j < constants::geometric::numNodesInTriangle; ++j)
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
        {
            continue;
        }
        validEdgesCount++;
    }

    std::vector<Edge> edges(validEdgesCount);
    validEdgesCount = 0;
    for (auto i = 0; i < triangulationWrapper.GetNumEdges(); ++i)
    {
        if (!edgeNodesFlag[i])
        {
            continue;
        }

        edges[validEdgesCount].first = triangulationWrapper.GetEdgeNode(i, 0);
        edges[validEdgesCount].second = triangulationWrapper.GetEdgeNode(i, 1);
        validEdgesCount++;
    }

    SetNodesRTreeRequiresUpdate(true);
    SetEdgesRTreeRequiresUpdate(true);

    m_edges = edges;
    m_nodes = inputNodes;
    m_projection = projection;

    DeleteInvalidNodesAndEdges();
    DoAdministration();
}

meshkernel::UInt Mesh2D::InvalidateEdgesWithNoFace()
{
    UInt numberInvalidated = 0;

    for (UInt e = 0; e < m_edgesFaces.size(); ++e)
    {
        if (m_edgesNumFaces[e] == 0)
        {
            m_edges[e].first = constants::missing::uintValue;
            m_edges[e].second = constants::missing::uintValue;
            ++numberInvalidated;
        }
    }

    return numberInvalidated;
}

void Mesh2D::DoAdministration(CompoundUndoAction* undoAction)
{
    if (!AdministrationRequired())
    {
        return;
    }

    AdministrateNodesEdges(undoAction);

    // face administration
    ResizeAndInitializeFaceVectors();

    // find faces
    FindFaces();

    // Compute face mass centers
    ComputeFaceAreaAndMassCenters();

    DeleteMeshHoles(undoAction);

    // classify node types
    ClassifyNodes();

    SetAdministrationRequired(false);
}

void Mesh2D::DoAdministrationGivenFaceNodesMapping(const std::vector<std::vector<UInt>>& faceNodes,
                                                   const std::vector<std::uint8_t>& numFaceNodes)
{
    AdministrateNodesEdges();

    // face administration
    ResizeAndInitializeFaceVectors();

    // find faces
    FindFacesGivenFaceNodesMapping(faceNodes, numFaceNodes);

    // Compute face mass centers
    ComputeFaceAreaAndMassCenters(true);

    m_invalidCellPolygons = ComputeInnerBoundaryPolygons();

    // classify node types
    ClassifyNodes();
}

void Mesh2D::Administrate(CompoundUndoAction* undoAction)
{
    DoAdministration(undoAction);
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

std::unique_ptr<meshkernel::UndoAction> Mesh2D::DeleteDegeneratedTriangles()
{
    std::unique_ptr<CompoundUndoAction> undoAction = CompoundUndoAction::Create();

    Administrate(undoAction.get());

    // assume the max amount of degenerated triangles is 10% of the actual faces
    std::vector<UInt> degeneratedTriangles;
    degeneratedTriangles.reserve(static_cast<UInt>(static_cast<double>(GetNumFaces()) * 0.1));
    for (UInt f = 0; f < GetNumFaces(); ++f)
    {
        const auto numFaceNodes = m_numFacesNodes[f];
        if (numFaceNodes != constants::geometric::numNodesInTriangle)
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
            for (UInt e = 0; e < constants::geometric::numNodesInTriangle; ++e)
            {
                const auto edge = m_facesEdges[f][e];
                undoAction->Add(ResetEdge(edge, {constants::missing::uintValue, constants::missing::uintValue}));
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

        undoAction->Add(ResetNode(thirdNode, m_facesMassCenters[face]));
        undoAction->Add(MergeTwoNodes(secondNode, firstNode));
        undoAction->Add(MergeTwoNodes(thirdNode, firstNode));
    }

    Administrate(undoAction.get());
    return undoAction;
}

bool Mesh2D::HasDuplicateNodes(const UInt numClosingEdges, const std::vector<UInt>& nodes, std::vector<UInt>& sortedNodes) const
{

    if (numClosingEdges == constants::geometric::numNodesInTriangle)
    {
        if (nodes[0] == nodes[1] || nodes[0] == nodes[2] || nodes[1] == nodes[2])
        {
            return true;
        }
    }
    else if (numClosingEdges == constants::geometric::numNodesInQuadrilateral)
    {
        if (nodes[0] == nodes[1] || nodes[0] == nodes[2] || nodes[0] == nodes[3] ||
            nodes[1] == nodes[2] || nodes[1] == nodes[3] ||
            nodes[2] == nodes[3])
        {
            return true;
        }
    }
    else
    {
        sortedNodes.clear();
        sortedNodes.reserve(nodes.size());
        std::copy(nodes.begin(), nodes.end(), std::back_inserter(sortedNodes));
        std::sort(sortedNodes.begin(), sortedNodes.end());
        for (UInt n = 0; n < sortedNodes.size() - 1; n++)
        {
            if (sortedNodes[n + 1] == sortedNodes[n])
            {
                return true;
            }
        }
    }

    return false;
}

bool Mesh2D::HasDuplicateEdgeFaces(const UInt numClosingEdges, const std::vector<UInt>& edges, std::vector<UInt>& sortedEdgesFaces) const
{

    // The number of edges is the same as the number of nodes for both triangles and quadrilateral
    if (numClosingEdges == constants::geometric::numNodesInTriangle)
    {
        if (m_edgesFaces[edges[0]][0] == m_edgesFaces[edges[1]][0] ||
            m_edgesFaces[edges[0]][0] == m_edgesFaces[edges[2]][0] ||
            m_edgesFaces[edges[1]][0] == m_edgesFaces[edges[2]][0])
        {
            return true;
        }
    }
    else if (numClosingEdges == constants::geometric::numNodesInQuadrilateral)
    {
        if (m_edgesFaces[edges[0]][0] == m_edgesFaces[edges[1]][0] ||
            m_edgesFaces[edges[0]][0] == m_edgesFaces[edges[2]][0] ||
            m_edgesFaces[edges[0]][0] == m_edgesFaces[edges[3]][0] ||

            m_edgesFaces[edges[1]][0] == m_edgesFaces[edges[2]][0] ||
            m_edgesFaces[edges[1]][0] == m_edgesFaces[edges[3]][0] ||

            m_edgesFaces[edges[2]][0] == m_edgesFaces[edges[3]][0])
        {
            return true;
        }
    }
    else
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
            {
                return true;
            }
        }
    }

    return false;
}

void Mesh2D::ResizeAndInitializeFaceVectors()
{
    // face administration
    m_edgesNumFaces.resize(m_edges.size());
    std::ranges::fill(m_edgesNumFaces, 0);

    m_edgesFaces.resize(m_edges.size());
    std::ranges::fill(m_edgesFaces, std::array{constants::missing::uintValue,
                                               constants::missing::uintValue});

    m_facesMassCenters.clear();
    m_faceArea.clear();
    m_facesNodes.clear();
    m_facesEdges.clear();
    m_numFacesNodes.clear();

    m_facesMassCenters.reserve(GetNumNodes());
    m_faceArea.reserve(GetNumNodes());
    m_facesNodes.reserve(GetNumNodes());
    m_facesEdges.reserve(GetNumNodes());
    m_numFacesNodes.reserve(GetNumNodes());
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
    {
        return;
    }

    if (m_edges[previousEdge].first == constants::missing::uintValue || m_edges[previousEdge].second == constants::missing::uintValue)
    {
        throw std::invalid_argument("Mesh2D::FindFacesRecursive: The selected edge is invalid. This should not happen since all invalid edges should have been cleaned up.");
    }

    // Check if the faces are already found
    if (m_edgesNumFaces[previousEdge] >= 2)
    {
        return;
    }

    edges.push_back(previousEdge);
    nodes.push_back(node);

    const auto otherNode = OtherNodeOfEdge(m_edges[previousEdge], node);

    // enclosure found
    if (otherNode == startNode && nodes.size() == numClosingEdges)
    {
        // no duplicated nodes allowed
        if (HasDuplicateNodes(numClosingEdges, nodes, sortedNodes))
        {
            return;
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

        // check if least one edge has no face and there are no duplicate edge-faces.
        if (!oneEdgeHasNoFace && HasDuplicateEdgeFaces(numClosingEdges, edges, sortedEdgesFaces))
        {
            return;
        }

        // the order of the edges in a new face must be counterclockwise
        // in order to evaluate the clockwise order, the signed face area is computed

        // The nodes array does not represent a closed polygon.
        auto const [area, center_of_mass, direction] = Polygon::FaceAreaAndCenterOfMass(m_nodes, nodes, m_projection, /* isClosed = */ false);

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
    std::vector<UInt> sortedEdgesFaces(constants::geometric::maximumNumberOfEdgesPerFace);
    std::vector<UInt> sortedNodes(constants::geometric::maximumNumberOfEdgesPerFace);
    std::vector<Point> nodalValues(constants::geometric::maximumNumberOfEdgesPerFace);
    std::vector<UInt> edges(constants::geometric::maximumNumberOfEdgesPerFace);
    std::vector<UInt> nodes(constants::geometric::maximumNumberOfEdgesPerFace);

    for (UInt numEdgesPerFace = constants::geometric::numNodesInTriangle; numEdgesPerFace <= constants::geometric::maximumNumberOfEdgesPerFace; ++numEdgesPerFace)
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

void Mesh2D::FindFacesGivenFaceNodesMapping(const std::vector<std::vector<UInt>>& faceNodes,
                                            const std::vector<std::uint8_t>& numFaceNodes)
{
    m_facesNodes = faceNodes;
    m_numFacesNodes = numFaceNodes;

    std::vector<UInt> local_edges;
    std::vector<Point> local_nodes;
    std::vector<UInt> local_node_indices;
    for (UInt f = 0; f < m_facesNodes.size(); ++f)
    {

        if (numFaceNodes[f] == 0)
        {
            continue;
        }

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
                throw AlgorithmError("FindFacesGivenMappings: m_edgesNumFaces > 2.");
            }
            m_edgesFaces[e][m_edgesNumFaces[e]] = f;
            m_edgesNumFaces[e] += 1;
        }

        local_nodes.emplace_back(local_nodes.front());

        auto [face_area, center_of_mass, direction] = Polygon::FaceAreaAndCenterOfMass(local_nodes, m_projection);

        m_faceArea.emplace_back(face_area);
        m_facesMassCenters.emplace_back(center_of_mass);
    }
}

void Mesh2D::ComputeFaceAreaAndMassCenters(bool computeMassCenters)
{

    if (!computeMassCenters)
    {
        return;
    }

    auto const numFaces = static_cast<int>(GetNumFaces());
    m_facesMassCenters.resize(numFaces);

    std::vector<Point> polygonNodesCache;
    // #pragma omp parallel for private(polygonNodesCache)
    for (int f = 0; f < numFaces; f++)
    {
        // need to account for spherical coordinates. Build a polygon around a face
        ComputeFaceClosedPolygon(f, polygonNodesCache);

        const auto [area, centerOfMass, direction] = Polygon::FaceAreaAndCenterOfMass(polygonNodesCache, m_projection);
        m_faceArea[f] = area;
        m_facesMassCenters[f] = centerOfMass;
    }
}

void Mesh2D::InitialiseBoundaryNodeClassification(std::vector<int>& intNodeType) const
{

    for (UInt e = 0; e < GetNumEdges(); e++)
    {
        const auto firstNode = m_edges[e].first;
        const auto secondNode = m_edges[e].second;

        if (firstNode == constants::missing::uintValue || secondNode == constants::missing::uintValue) [[unlikely]]
        {
            continue;
        }

        if (intNodeType[firstNode] == -1 || intNodeType[secondNode] == -1)
        {
            continue;
        }

        if (m_edgesNumFaces[e] == 0)
        {
            intNodeType[firstNode] = -1;
            intNodeType[secondNode] = -1;
        }
        if (IsEdgeOnBoundary(e))
        {
            intNodeType[firstNode] += 1;
            intNodeType[secondNode] += 1;
        }
    }
}

meshkernel::MeshNodeType Mesh2D::ClassifyNode(const UInt nodeId) const
{
    using enum MeshNodeType;

    MeshNodeType nodeType = Unspecified;

    if (m_nodesNumEdges[nodeId] == 2)
    {
        nodeType = Corner;
    }
    else
    {
        UInt firstNode = constants::missing::uintValue;
        UInt secondNode = constants::missing::uintValue;
        for (UInt i = 0; i < m_nodesNumEdges[nodeId]; ++i)
        {
            const auto edgeIndex = m_nodesEdges[nodeId][i];

            if (!IsEdgeOnBoundary(edgeIndex))
            {
                continue;
            }

            if (firstNode == constants::missing::uintValue)
            {
                firstNode = OtherNodeOfEdge(m_edges[edgeIndex], nodeId);
            }
            else
            {
                secondNode = OtherNodeOfEdge(m_edges[edgeIndex], nodeId);
                break;
            }
        }

        // point at the border
        nodeType = Boundary;

        if (firstNode != constants::missing::uintValue && secondNode != constants::missing::uintValue)
        {
            const double cosPhi = NormalizedInnerProductTwoSegments(m_nodes[nodeId], m_nodes[firstNode], m_nodes[nodeId], m_nodes[secondNode], m_projection);

            // threshold for corner points
            const double cornerCosine = 0.25;
            if (cosPhi > -cornerCosine)
            {
                // void angle
                nodeType = Corner;
            }
        }
    }

    return nodeType;
}

void Mesh2D::ClassifyNodes()
{
    using enum MeshNodeType;

    std::vector<int> intNodeType(GetNumNodes(), 0);
    m_nodesTypes.resize(GetNumNodes(), Unspecified);
    std::ranges::fill(m_nodesTypes, Unspecified);

    InitialiseBoundaryNodeClassification(intNodeType);

    for (UInt n = 0; n < GetNumNodes(); n++)
    {

        if (!Node(n).IsValid()) [[unlikely]]
        {
            continue;
        }

        if (intNodeType[n] == 1 || intNodeType[n] == 2)
        {
            m_nodesTypes[n] = ClassifyNode(n);
        }
        else if (intNodeType[n] > 2)
        {
            // corner point
            m_nodesTypes[n] = Corner;
        }
        else if (intNodeType[n] != -1)
        {
            // internal node
            m_nodesTypes[n] = Internal;
        }
        else
        {
            m_nodesTypes[n] = static_cast<MeshNodeType>(intNodeType[n]);
        }

        if (m_nodesNumEdges[n] < 2)
        {
            // hanging node
            m_nodesTypes[n] = Hanging;
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

std::unique_ptr<meshkernel::SphericalCoordinatesOffsetAction> Mesh2D::OffsetSphericalCoordinates(double minx, double maxx)
{
    // The nodes change in value, but not in any connectivity
    // So it is unnecessary to redo administration
    std::unique_ptr<SphericalCoordinatesOffsetAction> undoAction;

    if (m_projection == Projection::spherical && maxx - minx > 180.0)
    {
        undoAction = SphericalCoordinatesOffsetAction::Create(*this, minx, maxx);

        for (UInt n = 0; n < GetNumNodes(); ++n)
        {
            if (m_nodes[n].x - 360.0 >= minx)
            {
                undoAction->AddDecrease(n);
            }

            if (m_nodes[n].x < minx)
            {
                undoAction->AddIncrease(n);
            }
        }

        CommitAction(*undoAction);
    }

    return undoAction;
}

void Mesh2D::CommitAction(const SphericalCoordinatesOffsetAction& undoAction)
{
    undoAction.ApplyOffset(m_nodes);
}

void Mesh2D::CommitAction(PointArrayUndo& undoAction)
{
    undoAction.Swap(m_invalidCellPolygons);
    SetAdministrationRequired(true);
    SetNodesRTreeRequiresUpdate(true);
    SetEdgesRTreeRequiresUpdate(true);
    SetFacesRTreeRequiresUpdate(true);
}

void Mesh2D::RestoreAction(const SphericalCoordinatesOffsetAction& undoAction)
{
    undoAction.UndoOffset(m_nodes);
}

void Mesh2D::RestoreAction(PointArrayUndo& undoAction)
{
    undoAction.Swap(m_invalidCellPolygons);
    SetAdministrationRequired(true);
    SetNodesRTreeRequiresUpdate(true);
    SetEdgesRTreeRequiresUpdate(true);
    SetFacesRTreeRequiresUpdate(true);
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
    std::vector<Point> faceCircumcenters = algo::ComputeFaceCircumcenters(*this);

    for (UInt e = 0; e < GetNumEdges(); ++e)
    {
        const auto firstFace = m_edgesFaces[e][0];
        const auto secondFace = m_edgesFaces[e][1];

        if (firstFace != constants::missing::uintValue && secondFace != constants::missing::uintValue)
        {
            const auto flowEdgeLength = ComputeDistance(faceCircumcenters[firstFace], faceCircumcenters[secondFace], m_projection);
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
    std::vector<Point> faceCircumcenters = algo::ComputeFaceCircumcenters(*this);

    for (const auto& edge : edges)
    {
        const auto firstFace = m_edgesFaces[edge][0];
        const auto secondFace = m_edgesFaces[edge][1];
        result.emplace_back((faceCircumcenters[firstFace] + faceCircumcenters[secondFace]) * 0.5);
    }

    return result;
}

std::unique_ptr<meshkernel::UndoAction> Mesh2D::DeleteSmallFlowEdges(double smallFlowEdgesThreshold)
{
    std::unique_ptr<CompoundUndoAction> undoAction = CompoundUndoAction::Create();

    undoAction->Add(DeleteDegeneratedTriangles());

    if (auto edges = GetEdgesCrossingSmallFlowEdges(smallFlowEdgesThreshold); !edges.empty())
    {
        // invalidate the edges
        for (const auto& e : edges)
        {
            undoAction->Add(DeleteEdge(e));
        }

        Administrate(undoAction.get());
    }

    return undoAction;
}

void Mesh2D::ComputeAverageAreOfNeighbouringFaces(const UInt faceId, UInt& numNonBoundaryFaces, double& averageOtherFacesArea) const
{
    numNonBoundaryFaces = 0;
    averageOtherFacesArea = 0.0;

    if (m_numFacesNodes[faceId] != constants::geometric::numNodesInTriangle)
    {
        return;
    }

    for (UInt e = 0; e < constants::geometric::numNodesInTriangle; ++e)
    {
        // the edge must not be at the boundary, otherwise there is no "other" face
        const auto edge = m_facesEdges[faceId][e];
        if (IsEdgeOnBoundary(edge))
        {
            continue;
        }
        const auto otherFace = NextFace(faceId, edge);
        if (m_numFacesNodes[otherFace] > constants::geometric::numNodesInTriangle)
        {
            averageOtherFacesArea += m_faceArea[otherFace];
            numNonBoundaryFaces++;
        }
    }

    if (numNonBoundaryFaces != 0)
    {
        averageOtherFacesArea /= static_cast<double>(numNonBoundaryFaces);
    }
}

void Mesh2D::FindSmallestCornerAngle(const UInt faceId,
                                     double& minCosPhiSmallTriangle,
                                     UInt& nodeToPreserve,
                                     UInt& firstNodeToMerge,
                                     UInt& secondNodeToMerge,
                                     UInt& thirdEdgeSmallTriangle) const
{
    minCosPhiSmallTriangle = 1.0;
    nodeToPreserve = constants::missing::uintValue;
    firstNodeToMerge = constants::missing::uintValue;
    secondNodeToMerge = constants::missing::uintValue;
    thirdEdgeSmallTriangle = constants::missing::uintValue;

    for (UInt e = 0; e < constants::geometric::numNodesInTriangle; ++e)
    {
        const auto previousEdge = NextCircularBackwardIndex(e, constants::geometric::numNodesInTriangle);
        const auto nextEdge = NextCircularForwardIndex(e, constants::geometric::numNodesInTriangle);

        const auto k0 = m_facesNodes[faceId][previousEdge];
        const auto k1 = m_facesNodes[faceId][e];
        const auto k2 = m_facesNodes[faceId][nextEdge];

        // compute the angles between the edges
        const auto cosphi = std::abs(NormalizedInnerProductTwoSegments(m_nodes[k0], m_nodes[k1], m_nodes[k1], m_nodes[k2], m_projection));

        if (cosphi < minCosPhiSmallTriangle)
        {
            minCosPhiSmallTriangle = cosphi;
            nodeToPreserve = k1;
            firstNodeToMerge = k0;
            secondNodeToMerge = k2;
            thirdEdgeSmallTriangle = m_facesEdges[faceId][nextEdge];
        }
    }
}

void Mesh2D::DeleteSmallTriangle(const UInt nodeToPreserve,
                                 const UInt firstNodeToMerge,
                                 const UInt secondNodeToMerge,
                                 bool& nodesMerged,
                                 CompoundUndoAction& undoAction)
{
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
        undoAction.Add(MergeTwoNodes(firstNodeToMerge, nodeToPreserve));
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
        undoAction.Add(MergeTwoNodes(secondNodeToMerge, nodeToPreserve));
        nodesMerged = true;
    }
}

std::unique_ptr<meshkernel::UndoAction> Mesh2D::DeleteSmallTrianglesAtBoundaries(double minFractionalAreaTriangles)
{
    std::unique_ptr<CompoundUndoAction> undoAction = CompoundUndoAction::Create();

    // On the second part, the small triangles at the boundaries are checked
    std::vector<std::array<UInt, 3>> smallTrianglesNodes;
    for (UInt face = 0; face < GetNumFaces(); ++face)
    {
        if (m_numFacesNodes[face] != constants::geometric::numNodesInTriangle || m_faceArea[face] <= 0.0 || !IsFaceOnBoundary(face))
        {
            continue;
        }

        // compute the average area of neighboring faces
        double averageOtherFacesArea = 0.0;
        UInt numNonBoundaryFaces = 0;
        ComputeAverageAreOfNeighbouringFaces(face, numNonBoundaryFaces, averageOtherFacesArea);

        if (numNonBoundaryFaces == 0 || m_faceArea[face] / averageOtherFacesArea > minFractionalAreaTriangles)
        {
            // no valid boundary faces, the area of the current triangle is larger enough compared to the neighbors
            continue;
        }

        double minCosPhiSmallTriangle = 1.0;
        UInt nodeToPreserve = constants::missing::uintValue;
        UInt firstNodeToMerge = constants::missing::uintValue;
        UInt secondNodeToMerge = constants::missing::uintValue;
        UInt thirdEdgeSmallTriangle = constants::missing::uintValue;

        FindSmallestCornerAngle(face, minCosPhiSmallTriangle, nodeToPreserve, firstNodeToMerge, secondNodeToMerge, thirdEdgeSmallTriangle);

        if (thirdEdgeSmallTriangle != constants::missing::uintValue && IsEdgeOnBoundary(thirdEdgeSmallTriangle))
        {
            smallTrianglesNodes.emplace_back(std::array<UInt, 3>{nodeToPreserve, firstNodeToMerge, secondNodeToMerge});
        }
    }

    bool nodesMerged = false;
    for (const auto& triangleNodes : smallTrianglesNodes)
    {
        const auto nodeToPreserve = triangleNodes[0];
        const auto firstNodeToMerge = triangleNodes[1];
        const auto secondNodeToMerge = triangleNodes[2];

        DeleteSmallTriangle(nodeToPreserve, firstNodeToMerge, secondNodeToMerge, nodesMerged, *undoAction);
    }

    if (nodesMerged)
    {
        Administrate(undoAction.get());
    }

    return undoAction;
}

void Mesh2D::ComputeNodeNeighbours(std::vector<std::vector<UInt>>& nodesNodes, UInt& maxNumNeighbours) const
{
    maxNumNeighbours = *(std::max_element(m_nodesNumEdges.begin(), m_nodesNumEdges.end()));
    maxNumNeighbours += 1;

    ResizeAndFill2DVector(nodesNodes, GetNumNodes(), maxNumNeighbours, true, constants::missing::uintValue);

    // for each node, determine the neighboring nodes
    for (UInt n = 0; n < GetNumNodes(); n++)
    {

        for (UInt nn = 0; nn < m_nodesNumEdges[n]; nn++)
        {
            const auto edge = m_edges[m_nodesEdges[n][nn]];
            nodesNodes[n][nn] = OtherNodeOfEdge(edge, n);
        }
    }
}

void Mesh2D::ComputeAverageFlowEdgesLength(std::vector<double>& edgesLength,
                                           std::vector<double>& averageFlowEdgesLength) const
{
    std::vector<Point> faceCircumcenters = algo::ComputeFaceCircumcenters(*this);

    for (UInt e = 0; e < GetNumEdges(); e++)
    {
        const auto first = m_edges[e].first;
        const auto second = m_edges[e].second;

        if (first == second)
        {
            continue;
        }

        const double edgeLength = ComputeDistance(m_nodes[first], m_nodes[second], m_projection);
        edgesLength[e] = edgeLength;

        Point leftCenter;
        Point rightCenter;

        if (m_edgesNumFaces[e] > 0)
        {
            leftCenter = faceCircumcenters[m_edgesFaces[e][0]];
        }
        else
        {
            leftCenter = m_nodes[first];
        }

        // find right cell center, if it exists
        if (m_edgesNumFaces[e] == 2)
        {
            rightCenter = faceCircumcenters[m_edgesFaces[e][1]];
        }
        else
        {
            // otherwise, make ghost node by imposing boundary condition
            double dinry = InnerProductTwoSegments(m_nodes[first], m_nodes[second], m_nodes[first], leftCenter, m_projection);
            dinry = dinry / std::max(edgeLength * edgeLength, m_minimumEdgeLength);

            const Point bc = (1.0 - dinry) * m_nodes[first] + dinry * m_nodes[second];
            rightCenter = 2.0 * bc - leftCenter;
        }

        averageFlowEdgesLength[e] = ComputeDistance(leftCenter, rightCenter, m_projection);
    }
}

void Mesh2D::ComputeAverageEdgeLength(const std::vector<double>& edgesLength,
                                      const std::vector<double>& averageFlowEdgesLength,
                                      std::vector<bool>& curvilinearGridIndicator,
                                      std::vector<std::array<double, 2>>& averageEdgesLength,
                                      std::vector<double>& aspectRatios) const
{
    // Compute normal length
    for (UInt f = 0; f < GetNumFaces(); f++)
    {
        const auto numberOfFaceNodes = GetNumFaceEdges(f);
        if (numberOfFaceNodes < constants::geometric::numNodesInTriangle)
            continue;

        for (UInt n = 0; n < numberOfFaceNodes; n++)
        {
            if (numberOfFaceNodes != constants::geometric::numNodesInQuadrilateral)
            {
                curvilinearGridIndicator[m_facesNodes[f][n]] = false;
            }
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
            if (numberOfFaceNodes == constants::geometric::numNodesInQuadrilateral)
            {
                UInt kkp2 = n + 2;
                if (kkp2 >= numberOfFaceNodes)
                {
                    kkp2 = kkp2 - numberOfFaceNodes;
                }

                const auto klinkp2 = m_facesEdges[f][kkp2];
                edgeLength = 0.5 * (edgesLength[edgeIndex] + edgesLength[klinkp2]);
            }

            if (averageEdgesLength[edgeIndex][0] == constants::missing::doubleValue)
            {
                averageEdgesLength[edgeIndex][0] = edgeLength;
            }
            else
            {
                averageEdgesLength[edgeIndex][1] = edgeLength;
            }
        }
    }
}

void Mesh2D::ComputeAspectRatios(std::vector<double>& aspectRatios) const
{
    std::vector<std::array<double, 2>> averageEdgesLength(GetNumEdges(), std::array<double, 2>{constants::missing::doubleValue, constants::missing::doubleValue});
    std::vector<double> averageFlowEdgesLength(GetNumEdges(), constants::missing::doubleValue);
    std::vector<double> edgesLength(GetNumEdges(), 0.0);
    std::vector<bool> curvilinearGridIndicator(GetNumNodes(), true);
    aspectRatios.resize(GetNumEdges(), 0.0);

    ComputeAverageFlowEdgesLength(edgesLength, averageFlowEdgesLength);
    ComputeAverageEdgeLength(edgesLength, averageFlowEdgesLength,
                             curvilinearGridIndicator, averageEdgesLength, aspectRatios);

    if (m_curvilinearToOrthogonalRatio == 1.0)
    {
        return;
    }

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
            if (averageEdgesLength[e][0] != 0.0 && averageEdgesLength[e][0] != constants::missing::doubleValue)
            {
                aspectRatios[e] = averageFlowEdgesLength[e] / averageEdgesLength[e][0];
            }
        }
        else
        {
            if (averageEdgesLength[e][0] != 0.0 &&
                averageEdgesLength[e][1] != 0.0 &&
                averageEdgesLength[e][0] != constants::missing::doubleValue &&
                averageEdgesLength[e][1] != constants::missing::doubleValue)
            {
                aspectRatios[e] = m_curvilinearToOrthogonalRatio * aspectRatios[e] +
                                  (1.0 - m_curvilinearToOrthogonalRatio) * averageFlowEdgesLength[e] / (0.5 * (averageEdgesLength[e][0] + averageEdgesLength[e][1]));
            }
        }
    }
}

std::unique_ptr<meshkernel::UndoAction> Mesh2D::TriangulateFaces()
{
    // Create empty polygon
    Polygons polygon;
    return TriangulateFaces(polygon);
}

std::unique_ptr<meshkernel::UndoAction> Mesh2D::TriangulateFaces(const Polygons& polygon)
{
    std::unique_ptr<meshkernel::CompoundUndoAction> triangulationAction = CompoundUndoAction::Create();
    std::vector<Boolean> nodeInsidePolygon(IsLocationInPolygon(polygon, Location::Nodes));

    for (UInt i = 0; i < GetNumFaces(); ++i)
    {
        const UInt NumEdges = GetNumFaceEdges(i);

        if (NumEdges < 4)
        {
            continue;
        }

        bool elementIsOutsidePolygon = false;

        // Determine if any node is outside the polygon
        for (UInt j = 0; j < m_facesNodes[i].size(); ++j)
        {
            if (!nodeInsidePolygon[m_facesNodes[i][j]])
            {
                elementIsOutsidePolygon = true;
            }
        }

        // If any node of the element lies outside the polygon then do not triangulate.
        if (elementIsOutsidePolygon)
        {
            continue;
        }

        const auto indexFirstNode = m_facesNodes[i][0];
        for (UInt j = 2; j < NumEdges - 1; j++)
        {
            const auto nodeIndex = m_facesNodes[i][j];

            auto [edgeId, edgeConnectionAction] = ConnectNodes(indexFirstNode, nodeIndex);
            triangulationAction->Add(std::move(edgeConnectionAction));
        }
    }

    SetEdgesRTreeRequiresUpdate(true);

    return triangulationAction;
}

void Mesh2D::MakeDualFace(const std::span<const Point> edgeCentres,
                          UInt node,
                          double enlargementFactor,
                          std::vector<Point>& dualFace) const
{
    if (edgeCentres.size() != GetNumEdges())
    {
        throw ConstraintError("edgeCentre array is not the correct size: {} != {}",
                              edgeCentres.size(), GetNumEdges());
    }

    const auto sortedFacesIndices = SortedFacesAroundNode(node);
    const auto numEdges = m_nodesNumEdges[node];

    dualFace.reserve(constants::geometric::maximumNumberOfEdgesPerNode);
    dualFace.clear();

    if (sortedFacesIndices.empty())
    {
        return;
    }

    UInt sortedFacesCount = 0;

    for (UInt e = 0; e < numEdges; ++e)
    {
        const auto edgeIndex = m_nodesEdges[node][e];

        if (m_edgesNumFaces[edgeIndex] == 0)
        {
            continue;
        }

        auto edgeCenter = edgeCentres[edgeIndex];

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

        const auto faceIndex = sortedFacesIndices[sortedFacesCount];
        ++sortedFacesCount;

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

std::vector<meshkernel::Point> Mesh2D::ComputeBoundaryPolygons(const std::vector<Point>& polygonNodes)
{
    const Polygon polygon(polygonNodes, m_projection);

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

std::vector<meshkernel::Point> Mesh2D::ComputeInnerBoundaryPolygons() const
{
    if (GetNumFaces() == 0)
    {
        return std::vector<meshkernel::Point>();
    }

    std::vector<Point> illegalCells;
    illegalCells.reserve(GetNumNodes());
    std::vector<Point> meshBoundaryPolygon;
    meshBoundaryPolygon.reserve(GetNumNodes());
    std::vector<Point> subSequence;
    subSequence.reserve(GetNumNodes());

    std::vector<bool> edgeIsVisited(GetNumEdges(), false);
    std::vector<bool> nodeIsVisited(GetNumNodes(), false);

    std::vector<UInt> nodeIds;
    nodeIds.reserve(GetNumNodes());

    for (UInt e = 0; e < GetNumEdges(); e++)
    {
        if (edgeIsVisited[e] || !IsEdgeOnBoundary(e))
        {
            continue;
        }

        const auto firstNodeIndex = m_edges[e].first;
        const auto secondNodeIndex = m_edges[e].second;
        const auto firstNode = m_nodes[firstNodeIndex];
        const auto secondNode = m_nodes[secondNodeIndex];

        // Start a new polyline
        if (!subSequence.empty())
        {
            subSequence.emplace_back(constants::missing::doubleValue, constants::missing::doubleValue);
            nodeIds.emplace_back(constants::missing::uintValue);
        }

        // Put the current edge on the mesh boundary, mark it as visited
        const auto startPolygonEdges = static_cast<UInt>(subSequence.size());
        subSequence.emplace_back(firstNode);
        subSequence.emplace_back(secondNode);
        nodeIds.emplace_back(firstNodeIndex);
        nodeIds.emplace_back(secondNodeIndex);
        edgeIsVisited[e] = true;
        nodeIsVisited[firstNodeIndex] = true;
        nodeIsVisited[secondNodeIndex] = true;

        // walk the current mesh boundary
        auto currentNode = secondNodeIndex;
        WalkMultiBoundaryFromNode(edgeIsVisited, nodeIsVisited, currentNode, subSequence, nodeIds, meshBoundaryPolygon, illegalCells);

        const auto numNodesFirstTail = static_cast<UInt>(subSequence.size());

        // if the boundary polygon is not closed
        if (currentNode != firstNodeIndex)
        {
            // Now grow a polyline starting at the other side of the original link L, i.e., the second tail
            currentNode = firstNodeIndex;
            WalkMultiBoundaryFromNode(edgeIsVisited, nodeIsVisited, currentNode, subSequence, nodeIds, meshBoundaryPolygon, illegalCells);
        }

        // There is a nonempty second tail, so reverse the first tail, so that they connect.
        if (subSequence.size() > numNodesFirstTail)
        {
            const auto start = startPolygonEdges + static_cast<UInt>(std::ceil((numNodesFirstTail - startPolygonEdges + static_cast<UInt>(1)) * 0.5));

            for (auto n = start; n < numNodesFirstTail; n++)
            {
                const auto backupPoint = subSequence[n];
                const auto replaceIndex = numNodesFirstTail - n + firstNodeIndex;
                subSequence[n] = subSequence[replaceIndex];
                subSequence[replaceIndex] = backupPoint;

                const UInt backupPointIndex = nodeIds[n];
                nodeIds[n] = nodeIds[replaceIndex];
                nodeIds[replaceIndex] = backupPointIndex;
            }
        }
    }

    OrientatePolygonsAntiClockwise(illegalCells);

    return illegalCells;
}

void Mesh2D::OrientatePolygonsAntiClockwise(std::vector<Point>& polygonNodes) const
{
    UInt polygonStart = 0;
    UInt polygonLength = 0;
    UInt index = 0;

    while (index < polygonNodes.size())
    {
        polygonStart = index;
        polygonLength = 0;

        for (UInt i = polygonStart; i < polygonNodes.size(); ++i)
        {
            ++index;

            if (!polygonNodes[i].IsValid())
            {
                polygonLength = i - polygonStart;
                break;
            }

            if (index == polygonNodes.size())
            {
                ++polygonLength;
            }
        }

        if (polygonLength > 0)
        {
            const Point inValidPoint = {constants::missing::doubleValue, constants::missing::doubleValue};
            Point zeroPoint{0.0, 0.0};

            Point midPoint = std::accumulate(polygonNodes.begin() + polygonStart, polygonNodes.begin() + polygonStart + polygonLength - 1, zeroPoint) / static_cast<double>(polygonLength - 1);

            if (!IsPointInPolygonNodes(midPoint, polygonNodes, m_projection, inValidPoint, polygonStart, polygonStart + polygonLength))
            {
                // reverse order of polygon nodes
                if (polygonLength - 1 == 3)
                {
                    // Only the second and third points need be swapped to reverse the points in a triangle polygon
                    std::swap(polygonNodes[polygonStart + 1], polygonNodes[polygonStart + 2]);
                }
                else if (polygonLength - 1 == 4)
                {
                    // Only the second and fourth points need be swapped to reverse the points in a quadrilateral polygon
                    std::swap(polygonNodes[polygonStart + 1], polygonNodes[polygonStart + 3]);
                }
                else
                {
                    std::reverse(polygonNodes.begin() + polygonStart, polygonNodes.begin() + polygonStart + polygonLength);
                }
            }
        }
    }
}

std::vector<meshkernel::Point> Mesh2D::RemoveOuterDomainBoundaryPolygon(const std::vector<Point>& polygonNodes) const
{
    // Remove outer boundary.
    // It is assumed that the boundary polygon with the most points is from the outer domain boundary

    UInt outerPolygonLength = 0;
    UInt outerPolygonStart = 0;

    UInt outerPolygonLengthIntermediate = 0;
    UInt outerPolygonStartIntermediate = 0;

    UInt index = 0;

    // Find start index of outer boundary and the length of the polygon
    while (index < polygonNodes.size())
    {
        outerPolygonStartIntermediate = index;

        for (UInt i = outerPolygonStartIntermediate; i < polygonNodes.size(); ++i)
        {
            ++index;

            if (!polygonNodes[i].IsValid())
            {
                outerPolygonLengthIntermediate = i - outerPolygonStartIntermediate;
                break;
            }
        }

        if (outerPolygonLengthIntermediate > outerPolygonLength)
        {
            outerPolygonLength = outerPolygonLengthIntermediate;
            outerPolygonStart = outerPolygonStartIntermediate;
        }
    }

    if (outerPolygonLength > 0 && outerPolygonLength < polygonNodes.size())
    {
        ++outerPolygonLength;
    }

    std::vector<Point> innerBoundaryNodes;

    // Gather all node except those from the outer domain boundary polygon (this is assumed to be the polygon with the largest number of nodes)
    if (outerPolygonStart > 0)
    {
        innerBoundaryNodes.insert(innerBoundaryNodes.begin(), polygonNodes.begin(), polygonNodes.begin() + outerPolygonStart - 1);
    }

    if (outerPolygonLength > 0 && outerPolygonStart + outerPolygonLength < polygonNodes.size())
    {
        innerBoundaryNodes.insert(innerBoundaryNodes.end(), polygonNodes.begin() + outerPolygonStart + outerPolygonLength, polygonNodes.end());
    }

    return innerBoundaryNodes;
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

void Mesh2D::WalkMultiBoundaryFromNode(std::vector<bool>& edgeIsVisited,
                                       std::vector<bool>& nodeIsVisited,
                                       UInt& currentNode,
                                       std::vector<Point>& subSequence,
                                       std::vector<UInt>& nodeIds,
                                       std::vector<Point>& meshBoundaryPolygon,
                                       std::vector<Point>& illegalCells) const
{
    UInt e = 0;

    while (e < m_nodesNumEdges[currentNode])
    {
        const auto currentEdge = m_nodesEdges[currentNode][e];

        if (edgeIsVisited[currentEdge] || !IsEdgeOnBoundary(currentEdge))
        {
            e++;
            continue;
        }

        UInt nextNode = OtherNodeOfEdge(m_edges[currentEdge], currentNode);
        e = 0;

        if (nodeIsVisited[nextNode])
        {
            UInt lastIndex = constants::missing::uintValue;

            // Find index of last time node was added
            for (size_t ii = nodeIds.size(); ii >= 1; --ii)
            {
                UInt i = static_cast<UInt>(ii) - 1;

                if (nodeIds[i] == constants::missing::uintValue)
                {
                    break;
                }

                if (nodeIds[i] == nextNode)
                {
                    lastIndex = i;
                    break;
                }
            }

            if (lastIndex != constants::missing::uintValue)
            {
                size_t start = meshBoundaryPolygon.size();

                if (!meshBoundaryPolygon.empty())
                {
                    meshBoundaryPolygon.emplace_back(constants::missing::doubleValue, constants::missing::doubleValue);
                    ++start;
                }

                for (size_t i = lastIndex; i < nodeIds.size(); ++i)
                {
                    meshBoundaryPolygon.emplace_back(subSequence[i]);
                }

                meshBoundaryPolygon.emplace_back(subSequence[lastIndex]);

                std::span<const Point> currentPolygon(meshBoundaryPolygon.data() + start, meshBoundaryPolygon.data() + meshBoundaryPolygon.size());
                // Since the edge lies on a boundary, there will be only 1 attached element.
                // This element will be in the 0th position
                UInt connectedFace = m_edgesFaces[currentEdge][0];
                bool isInPolygon = IsPointInPolygonNodes(m_facesMassCenters[connectedFace], currentPolygon, m_projection);

                if (!isInPolygon)
                {

                    if (!illegalCells.empty())
                    {
                        // If illegal cells array is not empty then add the polygon separator
                        illegalCells.emplace_back(constants::missing::doubleValue, constants::missing::doubleValue);
                    }

                    illegalCells.insert (illegalCells.end (), currentPolygon.begin (), currentPolygon.end ());
                }

                subSequence.resize(lastIndex);
                nodeIds.resize(lastIndex);
            }
        }

        currentNode = nextNode;
        subSequence.emplace_back(m_nodes[currentNode]);
        edgeIsVisited[currentEdge] = true;
        nodeIsVisited[currentNode] = true;
        nodeIds.emplace_back(currentNode);
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

std::vector<bool> Mesh2D::FilterBasedOnMetric(Location location,
                                              Property property,
                                              double minValue,
                                              double maxValue) const
{
    if (location != Location::Faces)
    {
        throw ConstraintError("Unsupported location. Only location faces is supported");
    }
    if (property != Property::Orthogonality)
    {
        throw ConstraintError("Unsupported metric. Only orthogonality metric is supported");
    }

    const auto numFaces = GetNumFaces();

    // Pre-allocate memory for result vector and set all elements to false
    std::vector<bool> result(numFaces, false);

    // Retrieve orthogonality values
    MeshOrthogonality meshOrthogonality;
    const std::vector<double> metricValues(meshOrthogonality.Compute(*this));

    // Loop through faces and compute how many edges have the metric within the range
    for (UInt f = 0; f < numFaces; ++f)
    {
        const UInt numFaceEdges = GetNumFaceEdges(f);
        UInt numEdgesFiltered = 0;
        for (UInt e = 0; e < numFaceEdges; ++e)
        {
            const auto edge = m_facesEdges[f][e];
            const double metricValue = metricValues[edge];
            if (metricValue < minValue || metricValue > maxValue)
            {
                break;
            }
            numEdgesFiltered = numEdgesFiltered + 1;
        }

        // If all edges have the metric within the range, the face is masked
        if (numEdgesFiltered == numFaceEdges)
        {
            result[f] = true;
        }
    }

    return result;
}

void Mesh2D::FindNodesToDelete(const Polygons& polygon,
                               const bool invertDeletion,
                               std::vector<bool>& isNodeInsidePolygon,
                               std::vector<bool>& deleteNode) const
{
    for (UInt n = 0; n < GetNumNodes(); ++n)
    {
        if (m_nodes[n].IsValid())
        {
            auto [isInPolygon, polygonIndex] = polygon.IsPointInPolygons(m_nodes[n]);

            if (isInPolygon)
            {
                isNodeInsidePolygon[n] = true;
                deleteNode[n] = !invertDeletion;
            }
        }
        else
        {
            isNodeInsidePolygon[n] = false;
            deleteNode[n] = false;
        }
    }
}

std::vector<bool> Mesh2D::FindFacesEntirelyInsidePolygon(const std::vector<bool>& isNodeInsidePolygon) const
{

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

    return isFaceCompletlyIncludedInPolygon;
}

void Mesh2D::DeletedMeshNodesAndEdges(const std::function<bool(UInt)>& excludedFace,
                                      std::vector<bool>& deleteNode,
                                      CompoundUndoAction& deleteMeshAction)
{
    std::vector<std::uint8_t> nodeEdgeCount(m_nodesNumEdges);
    std::vector<bool> includeFacePolygon(GetNumFaces(), false);

    for (UInt e = 0; e < GetNumEdges(); ++e)
    {
        const auto numEdgeFaces = GetNumEdgesFaces(e);

        if (numEdgeFaces == 1 && excludedFace(m_edgesFaces[e][0]))
        {
            deleteNode[m_edges[e].first] = false;
            deleteNode[m_edges[e].second] = false;
            continue;
        }

        if (numEdgeFaces == 2 && (excludedFace(m_edgesFaces[e][0]) || excludedFace(m_edgesFaces[e][1])))
        {
            deleteNode[m_edges[e].first] = false;
            deleteNode[m_edges[e].second] = false;
            continue;
        }

        if (m_edges[e].first != constants::missing::uintValue && nodeEdgeCount[m_edges[e].first] > 0)
        {
            --nodeEdgeCount[m_edges[e].first];
        }

        if (m_edges[e].second != constants::missing::uintValue && nodeEdgeCount[m_edges[e].second] > 0)
        {
            --nodeEdgeCount[m_edges[e].second];
        }

        if (m_edgesFaces[e][0] != constants::missing::uintValue && !excludedFace(m_edgesFaces[e][0]))
        {
            includeFacePolygon[m_edgesFaces[e][0]] = true;
        }

        if (m_edgesFaces[e][1] != constants::missing::uintValue && !excludedFace(m_edgesFaces[e][1]))
        {
            includeFacePolygon[m_edgesFaces[e][1]] = true;
        }

        deleteMeshAction.Add(DeleteEdge(e));
    }

    std::vector<UInt> nodesToDelete;
    nodesToDelete.reserve(GetNumNodes());

    for (UInt i = 0; i < m_nodes.size(); ++i)
    {
        if ((deleteNode[i] || nodeEdgeCount[i] == 0) && m_nodes[i].IsValid())
        {

            for (UInt edgeId = 0; edgeId < m_nodesNumEdges[i]; ++edgeId)
            {

                if (m_edgesFaces[edgeId][0] != constants::missing::uintValue && !excludedFace(m_edgesFaces[edgeId][0]))
                {
                    includeFacePolygon[m_edgesFaces[edgeId][0]] = true;
                }

                if (m_edgesFaces[edgeId][1] != constants::missing::uintValue && !excludedFace(m_edgesFaces[edgeId][1]))
                {
                    includeFacePolygon[m_edgesFaces[edgeId][1]] = true;
                }
            }

            nodesToDelete.push_back(i);
        }
    }

    deleteMeshAction.Add(PointArrayUndo::Create(*this, m_invalidCellPolygons));

    // Must append cell polygons before deleting node, otherwise the node value will be invalid.
    for (UInt faceId = 0; faceId < includeFacePolygon.size(); ++faceId)
    {
        if (includeFacePolygon[faceId])
        {
            AppendCellPolygon(faceId);
        }
    }

    for (const UInt nodeId : nodesToDelete)
    {
        deleteMeshAction.Add(DeleteNode(nodeId));
    }
}

void Mesh2D::AppendCellPolygon(const UInt faceId)
{
    if (faceId == constants::missing::uintValue)
    {
        throw ConstraintError("Face index is invalid");
    }

    size_t pointIndex = m_invalidCellPolygons.size();

    if (m_invalidCellPolygons.empty())
    {
        // +1 for the closing node
        m_invalidCellPolygons.resize(m_numFacesNodes[faceId] + 1);
    }
    else
    {
        // +2: 1 for the closing node and another for the separation value
        m_invalidCellPolygons.resize(m_invalidCellPolygons.size() + m_numFacesNodes[faceId] + 2);
        // Add separator
        m_invalidCellPolygons[pointIndex] = Point{constants::missing::doubleValue, constants::missing::doubleValue};
        ++pointIndex;
    }

    for (UInt i = 0; i < m_numFacesNodes[faceId]; ++i)
    {
        m_invalidCellPolygons[pointIndex] = m_nodes[m_facesNodes[faceId][i]];
        ++pointIndex;
    }

    // Close the polygon
    m_invalidCellPolygons[pointIndex] = m_nodes[m_facesNodes[faceId][0]];
}

void Mesh2D::DeleteMeshHoles(CompoundUndoAction* undoAction)
{
    if (!m_invalidCellPolygons.empty())
    {
        Polygons invalidElementPolygon(m_invalidCellPolygons, m_projection);
        auto deleteFacesUndo = DeleteMeshFacesInPolygon(invalidElementPolygon, false);

        if (undoAction != nullptr)
        {
            undoAction->Add(std::move(deleteFacesUndo));
        }
    }
}

std::unique_ptr<meshkernel::UndoAction> Mesh2D::DeleteMeshFacesInPolygon(const Polygons& polygon, const bool appendDeletedFaces)
{

    // A mapping between the old face-id and the new face-id.
    // Any deleted elements will be the invalid uint value.
    std::vector<UInt> faceIndices(GetNumFaces(), constants::missing::uintValue);
    // Indicate if the cell has a valid cell centre and is inside the polygons or regions to be deleted.
    std::vector<Boolean> validAndInside(GetNumFaces(), false);
    UInt faceIndex = 0;

    // Compute the expensive part (is point in polygon for all cells) in parallel
#pragma omp parallel for
    for (int f = 0; f < static_cast<int>(GetNumFaces()); ++f)
    {
        validAndInside[f] = m_facesMassCenters[f].IsValid() && polygon.IsPointInAnyPolygon(m_facesMassCenters[f]);
    }

    for (UInt f = 0u; f < GetNumFaces(); ++f)
    {
        if (!validAndInside[f])
        {
            faceIndices[f] = faceIndex;
            ++faceIndex;
        }
    }

    return UpdateFaceInformation(faceIndices, appendDeletedFaces);
}

std::unique_ptr<meshkernel::UndoAction> Mesh2D::UpdateFaceInformation(const std::vector<UInt>& faceIndices, const bool appendDeletedFaces)
{
    std::vector<UInt> facesToDelete;
    facesToDelete.reserve(GetNumFaces());

    // Collect face-ids of faces that are marked for deletion
    for (UInt i = 0; i < faceIndices.size(); ++i)
    {
        if (faceIndices[i] == constants::missing::uintValue)
        {
            facesToDelete.push_back(i);
        }
    }

    if (facesToDelete.empty())
    {
        // No elements to be deleted
        return nullptr;
    }

    std::unique_ptr<meshkernel::CompoundUndoAction> deleteMeshAction = CompoundUndoAction::Create();
    deleteMeshAction->Add(PointArrayUndo::Create(*this, m_invalidCellPolygons));

    // The edges of one of more faces will be deleted, so indicate that an administrate will be required.
    SetAdministrationRequired(true);

    if (appendDeletedFaces)
    {
        for (UInt faceId : facesToDelete)
        {
            AppendCellPolygon(faceId);
        }
    }

    // Remove deleted faces from edge-face connectivity
    for (UInt faceId : facesToDelete)
    {
        for (UInt e = 0u; e < m_facesEdges[faceId].size(); ++e)
        {
            UInt edge = m_facesEdges[faceId][e];
            bool removedFace = false;

            if (m_edgesFaces[edge][0] == faceId)
            {
                m_edgesFaces[edge][0] = constants::missing::uintValue;
                --m_edgesNumFaces[edge];
                // Ensure the first connected face is on the 0th side
                // If the 1st side is already the invalid value then this operation does not change anything.
                std::swap(m_edgesFaces[edge][0], m_edgesFaces[edge][1]);
                removedFace = true;
            }
            else if (m_edgesFaces[edge][1] == faceId)
            {
                m_edgesFaces[edge][1] = constants::missing::uintValue;
                --m_edgesNumFaces[edge];
                removedFace = true;
            }

            if (removedFace && m_edgesNumFaces[edge] == 0)
            {
                deleteMeshAction->Add(DeleteEdge(edge));
            }
        }
    }

    // Renumber existing edge-face ids to match new face-ids
    for (size_t edge = 0; edge < m_edgesFaces.size(); ++edge)
    {
        if (m_edgesFaces[edge][0] != constants::missing::uintValue)
        {
            m_edgesFaces[edge][0] = faceIndices[m_edgesFaces[edge][0]];
        }

        if (m_edgesFaces[edge][1] != constants::missing::uintValue)
        {
            m_edgesFaces[edge][1] = faceIndices[m_edgesFaces[edge][1]];
        }
    }

    // Shift face connectivity in arrays where deleted faces have been removed from arrays
    // Loop must be iterated in reverse order, from highest face id value to lowest.
    for (UInt faceId : facesToDelete | std::views::reverse)
    {
        m_facesNodes.erase(m_facesNodes.begin() + faceId);
        m_numFacesNodes.erase(m_numFacesNodes.begin() + faceId);
        m_facesEdges.erase(m_facesEdges.begin() + faceId);
        m_facesMassCenters.erase(m_facesMassCenters.begin() + faceId);
        m_faceArea.erase(m_faceArea.begin() + faceId);
    }

    return deleteMeshAction;
}

std::unique_ptr<meshkernel::UndoAction> Mesh2D::DeleteMesh(const Polygons& polygon, DeleteMeshOptions deletionOption, bool invertDeletion)
{

    if (deletionOption == FacesWithIncludedCircumcenters)
    {
        return DeleteMeshFaces(polygon, invertDeletion);
    }

    std::unique_ptr<CompoundUndoAction> deleteMeshAction = CompoundUndoAction::Create();

    // Find crossed faces
    Mesh2DIntersections mesh2DIntersections(*this);
    mesh2DIntersections.Compute(polygon);
    const auto& faceIntersections = mesh2DIntersections.FaceIntersections();

    // Find faces with all nodes inside the polygon
    std::vector<bool> isNodeInsidePolygon(GetNumNodes(), false);
    std::vector<bool> deleteNode(GetNumNodes(), invertDeletion);

    FindNodesToDelete(polygon, invertDeletion, isNodeInsidePolygon, deleteNode);

    std::vector<bool> isFaceCompletlyIncludedInPolygon(FindFacesEntirelyInsidePolygon(isNodeInsidePolygon));

    std::function<bool(UInt)> excludedFace;

    if (deletionOption == InsideNotIntersected && !invertDeletion)
    {
        excludedFace = [&isFaceCompletlyIncludedInPolygon, &faceIntersections](UInt f)
        { return !isFaceCompletlyIncludedInPolygon[f] || faceIntersections[f].faceIndex != constants::missing::uintValue; };
    }
    else if (deletionOption == InsideNotIntersected && invertDeletion)
    {
        excludedFace = [&isFaceCompletlyIncludedInPolygon, &faceIntersections](UInt f)
        { return isFaceCompletlyIncludedInPolygon[f] && faceIntersections[f].faceIndex == constants::missing::uintValue; };
    }
    else if (deletionOption == InsideAndIntersected && !invertDeletion)
    {
        excludedFace = [&isFaceCompletlyIncludedInPolygon, &faceIntersections](UInt f)
        { return !isFaceCompletlyIncludedInPolygon[f] && faceIntersections[f].faceIndex == constants::missing::uintValue; };
    }
    else if (deletionOption == InsideAndIntersected && invertDeletion)
    {
        excludedFace = [&isFaceCompletlyIncludedInPolygon, &faceIntersections](UInt f)
        { return isFaceCompletlyIncludedInPolygon[f] || faceIntersections[f].faceIndex != constants::missing::uintValue; };
    }

    // Delete nodes and edges marked for deletion
    DeletedMeshNodesAndEdges(excludedFace, deleteNode, *deleteMeshAction);

    SetNodesRTreeRequiresUpdate(true);
    SetEdgesRTreeRequiresUpdate(true);

    Administrate(deleteMeshAction.get());
    return deleteMeshAction;
}

std::unique_ptr<meshkernel::UndoAction> Mesh2D::DeleteMeshFaces(const Polygons& polygon, bool invertDeletion)
{
    std::unique_ptr<meshkernel::CompoundUndoAction> deleteMeshAction = CompoundUndoAction::Create();

    Administrate(deleteMeshAction.get());
    std::vector<Point> faceCircumcenters = algo::ComputeFaceCircumcenters(*this);
    std::vector<bool> includeFace(GetNumFaces(), false);

    for (UInt e = 0u; e < GetNumEdges(); ++e)
    {
        bool allFaceCircumcentersInPolygon = true;

        for (UInt f = 0u; f < GetNumEdgesFaces(e); ++f)
        {
            const auto faceIndex = m_edgesFaces[e][f];
            if (faceIndex == constants::missing::uintValue)
            {
                continue;
            }

            auto [isInPolygon, polygonIndex] = polygon.IsPointInPolygons(faceCircumcenters[faceIndex]);

            if (invertDeletion)
            {
                isInPolygon = !isInPolygon;
            }
            if (!isInPolygon)
            {
                allFaceCircumcentersInPolygon = false;
                break;
            }
        }

        // 2D edge without surrounding faces.
        if (GetNumEdgesFaces(e) == 0)
        {
            const auto firstNodeIndex = m_edges[e].first;
            const auto secondNodeIndex = m_edges[e].second;

            if (firstNodeIndex == constants::missing::uintValue || secondNodeIndex == constants::missing::uintValue)
            {
                continue;
            }

            const auto edgeCenter = (m_nodes[firstNodeIndex] + m_nodes[secondNodeIndex]) / 2.0;

            auto [isInPolygon, polygonIndex] = polygon.IsPointInPolygons(edgeCenter);
            allFaceCircumcentersInPolygon = isInPolygon;

            if (invertDeletion)
            {
                allFaceCircumcentersInPolygon = !allFaceCircumcentersInPolygon;
            }
        }

        if (allFaceCircumcentersInPolygon)
        {
            if (m_edgesFaces[e][0] != constants::missing::uintValue)
            {
                includeFace[m_edgesFaces[e][0]] = true;
            }

            if (m_edgesFaces[e][1] != constants::missing::uintValue)
            {
                includeFace[m_edgesFaces[e][1]] = true;
            }

            deleteMeshAction->Add(DeleteEdge(e));
        }
    }

    deleteMeshAction->Add(PointArrayUndo::Create(*this, m_invalidCellPolygons));

    for (UInt faceId = 0; faceId < includeFace.size(); ++faceId)
    {
        if (includeFace[faceId])
        {
            AppendCellPolygon(faceId);
        }
    }

    Administrate(deleteMeshAction.get());

    return deleteMeshAction;
}

std::unique_ptr<meshkernel::UndoAction> Mesh2D::DeleteHangingEdges()
{
    std::unique_ptr<meshkernel::CompoundUndoAction> deleteAction = CompoundUndoAction::Create();

    for (const auto& hangingEdge : GetHangingEdges())
    {
        deleteAction->Add(DeleteEdge(hangingEdge));
    }

    return deleteAction;
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
        const auto edgeIndex = FindLocationIndex(points[i], Location::Edges);

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

std::vector<int> Mesh2D::ComputeNodeMask(const Polygons& polygons) const
{
    std::vector<int> nodeMask(GetNumNodes(), 0);

    for (UInt n = 0; n < GetNumNodes(); ++n)
    {
        const auto [isInPolygon, polygonIndex] = polygons.IsPointInPolygons(m_nodes[n]);

        if (isInPolygon)
        {
            nodeMask[n] = 1;
        }
    }

    return nodeMask;
}

std::vector<int> Mesh2D::ComputeEdgeMask(const std::vector<int>& nodeMask,
                                         bool includeIntersected) const
{
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

    return edgeMask;
}

void Mesh2D::RemoveIntersected(const std::vector<int>& edgeMask,
                               std::vector<int>& secondEdgeMask) const
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

void Mesh2D::InvertSelection(const std::vector<int>& edgeMask,
                             std::vector<int>& secondEdgeMask) const
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

std::vector<int> Mesh2D::MaskEdgesOfFacesInPolygon(const Polygons& polygons,
                                                   bool invertSelection,
                                                   bool includeIntersected) const
{
    // mark all nodes in polygon with 1
    std::vector<int> nodeMask(ComputeNodeMask(polygons));

    // mark all edges with both start end end nodes included with 1
    std::vector<int> edgeMask(ComputeEdgeMask(nodeMask, includeIntersected));

    // if one edge of the face is not included do not include all the edges of that face
    auto secondEdgeMask = edgeMask;

    if (!includeIntersected)
    {
        RemoveIntersected(edgeMask, secondEdgeMask);
    }

    // if the selection is inverted, do not delete the edges of included faces
    if (invertSelection)
    {
        InvertSelection(edgeMask, secondEdgeMask);
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

std::unique_ptr<Mesh2D> Mesh2D::Merge(const Mesh2D& mesh1, const Mesh2D& mesh2)
{

    if (mesh1.m_projection != mesh2.m_projection)
    {
        throw MeshKernelError("The two meshes cannot be merged: they have different projections");
    }

    if ((mesh2.GetNumNodes() == 0 || mesh2.GetNumEdges() == 0) && (mesh1.GetNumNodes() == 0 || mesh1.GetNumEdges() == 0))
    {
        throw MeshKernelError("The two meshes cannot be merged: both meshes are empty");
    }

    if ((mesh1.GetNumNodes() == 0 || mesh1.GetNumEdges() == 0) && (mesh2.GetNumNodes() > 0 && mesh2.GetNumEdges() > 0))
    {
        return std::make_unique<Mesh2D>(mesh2.m_edges,
                                        mesh2.m_nodes,
                                        mesh2.m_projection);
    }

    if ((mesh2.GetNumNodes() == 0 || mesh2.GetNumEdges() == 0) && (mesh1.GetNumNodes() > 0 && mesh1.GetNumEdges() > 0))
    {
        return std::make_unique<Mesh2D>(mesh1.m_edges,
                                        mesh1.m_nodes,
                                        mesh1.m_projection);
    }

    // Initialise with mesh1,
    Mesh2D mergedMesh(mesh1.m_edges, mesh1.m_nodes, mesh1.m_projection);

    const auto mesh1NodeOffset = static_cast<UInt>(mergedMesh.m_nodes.size());
    const auto mesh1EdgeOffset = static_cast<UInt>(mergedMesh.m_edges.size());
    const auto mesh1FaceOffset = static_cast<UInt>(mergedMesh.m_numFacesNodes.size());

    // Merge node arrays
    mergedMesh.m_nodes.insert(mergedMesh.m_nodes.end(), mesh2.m_nodes.begin(), mesh2.m_nodes.end());

    // Merge edge arrays
    mergedMesh.m_edges.insert(mergedMesh.m_edges.end(), mesh2.m_edges.begin(), mesh2.m_edges.end());

    // Update edge-node indices
    for (UInt i = 0; i < mesh2.m_edges.size(); ++i)
    {
        IncrementValidValue(mergedMesh.m_edges[i + mesh1EdgeOffset].first, mesh1NodeOffset);
        IncrementValidValue(mergedMesh.m_edges[i + mesh1EdgeOffset].second, mesh1NodeOffset);
    }

    //--------------------------------

    // Merge node-edge arrays
    mergedMesh.m_nodesEdges.insert(mergedMesh.m_nodesEdges.end(), mesh2.m_nodesEdges.begin(), mesh2.m_nodesEdges.end());

    for (UInt i = 0; i < mesh2.m_nodesEdges.size(); ++i)
    {
        for (UInt j = 0; j < mesh2.m_nodesEdges[i].size(); ++j)
        {
            IncrementValidValue(mergedMesh.m_nodesEdges[i + mesh1NodeOffset][j], mesh1EdgeOffset);
        }
    }

    //--------------------------------

    // Merge face-node arrays
    mergedMesh.m_facesNodes.insert(mergedMesh.m_facesNodes.end(), mesh2.m_facesNodes.begin(), mesh2.m_facesNodes.end());

    for (UInt i = 0; i < mesh2.m_facesNodes.size(); ++i)
    {
        for (UInt j = 0; j < mesh2.m_facesNodes[i].size(); ++j)
        {
            IncrementValidValue(mergedMesh.m_facesNodes[i + mesh1FaceOffset][j], mesh1NodeOffset);
        }
    }

    //--------------------------------

    // Merge edge-face arrays
    mergedMesh.m_edgesFaces.insert(mergedMesh.m_edgesFaces.end(), mesh2.m_edgesFaces.begin(), mesh2.m_edgesFaces.end());

    for (UInt i = 0; i < mesh2.m_edgesFaces.size(); ++i)
    {
        IncrementValidValue(mergedMesh.m_edgesFaces[i + mesh1EdgeOffset][0], mesh1FaceOffset);
        IncrementValidValue(mergedMesh.m_edgesFaces[i + mesh1EdgeOffset][1], mesh1FaceOffset);
    }

    //--------------------------------

    // Merge face-edge arrays
    mergedMesh.m_facesEdges.insert(mergedMesh.m_facesEdges.end(), mesh2.m_facesEdges.begin(), mesh2.m_facesEdges.end());

    for (UInt i = 0; i < mesh2.m_facesEdges.size(); ++i)
    {
        for (UInt j = 0; j < mesh2.m_facesEdges[i].size(); ++j)
        {
            IncrementValidValue(mergedMesh.m_facesEdges[i + mesh1FaceOffset][j], mesh1EdgeOffset);
        }
    }

    //--------------------------------

    // Now merge remaining arrays

    mergedMesh.m_nodesNumEdges.insert(mergedMesh.m_nodesNumEdges.end(), mesh2.m_nodesNumEdges.begin(), mesh2.m_nodesNumEdges.end());
    mergedMesh.m_nodesTypes.insert(mergedMesh.m_nodesTypes.end(), mesh2.m_nodesTypes.begin(), mesh2.m_nodesTypes.end());

    mergedMesh.m_edgesNumFaces.insert(mergedMesh.m_edgesNumFaces.end(), mesh2.m_edgesNumFaces.begin(), mesh2.m_edgesNumFaces.end());

    mergedMesh.m_numFacesNodes.insert(mergedMesh.m_numFacesNodes.end(), mesh2.m_numFacesNodes.begin(), mesh2.m_numFacesNodes.end());
    mergedMesh.m_facesMassCenters.insert(mergedMesh.m_facesMassCenters.end(), mesh2.m_facesMassCenters.begin(), mesh2.m_facesMassCenters.end());
    mergedMesh.m_faceArea.insert(mergedMesh.m_faceArea.end(), mesh2.m_faceArea.begin(), mesh2.m_faceArea.end());

    // Indicate that the mesh state has changed and the r-trees will need to be re-computed when required.
    mergedMesh.SetNodesRTreeRequiresUpdate(true);
    mergedMesh.SetEdgesRTreeRequiresUpdate(true);
    mergedMesh.SetFacesRTreeRequiresUpdate(true);

    mergedMesh.SetAdministrationRequired(true);

    return std::make_unique<Mesh2D>(mergedMesh.m_edges,
                                    mergedMesh.m_nodes,
                                    mergedMesh.m_projection);
}

std::unique_ptr<meshkernel::Mesh2D> Mesh2D::Merge(const std::span<const Point>& mesh1Nodes,
                                                  const std::span<const Edge>& mesh1Edges,
                                                  const std::span<const Point>& mesh2Nodes,
                                                  const std::span<const Edge>& mesh2Edges,
                                                  const Projection projection)
{
    std::vector<Point> mergedNodes(mesh1Nodes.size() + mesh2Nodes.size());
    std::vector<Edge> mergedEdges(mesh1Edges.size() + mesh2Edges.size());

    if (!mesh1Nodes.empty())
    {
        // Merge node array from mesh1 nodes
        std::ranges::copy(mesh1Nodes, mergedNodes.begin());

        // Merge edge array from mesh1 edges
        std::ranges::copy(mesh1Edges, mergedEdges.begin());
    }

    if (!mesh2Nodes.empty())
    {
        // Merge node array from mesh2 nodes
        std::ranges::copy(mesh2Nodes, mergedNodes.begin() + mesh1Nodes.size());

        // Merge edge array from mesh2 edges
        std::ranges::copy(mesh2Edges, mergedEdges.begin() + mesh1Edges.size());

        if (!mesh1Nodes.empty())
        {
            const UInt nodeOffset = static_cast<UInt>(mesh1Nodes.size());

            // Re-assign the node ids for the second mesh data set
            for (size_t i = mesh1Edges.size(); i < mergedEdges.size(); ++i)
            {
                IncrementValidValue(mergedEdges[i].first, nodeOffset);
                IncrementValidValue(mergedEdges[i].second, nodeOffset);
            }
        }
    }

    return std::make_unique<Mesh2D>(mergedEdges, mergedNodes, projection);
}

meshkernel::BoundingBox Mesh2D::GetBoundingBox() const
{
    Point lowerLeft(std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
    Point upperRight(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest());

    const auto numNodes = GetNumNodes();
    for (UInt n = 0; n < numNodes; ++n)
    {
        [[unlikely]] if (!m_nodes[n].IsValid())
        {
            continue;
        }

        lowerLeft.x = std::min(lowerLeft.x, m_nodes[n].x);
        lowerLeft.y = std::min(lowerLeft.y, m_nodes[n].y);
        upperRight.x = std::max(upperRight.x, m_nodes[n].x);
        upperRight.y = std::max(upperRight.y, m_nodes[n].y);
    }

    return {lowerLeft, upperRight};
}

std::vector<meshkernel::BoundingBox> Mesh2D::GetEdgesBoundingBoxes() const
{
    std::vector<BoundingBox> result;
    result.reserve(GetNumEdges());
    for (const auto& e : m_edges)
    {
        if (e.first == constants::missing::uintValue || e.second == constants::missing::uintValue)
        {
            result.emplace_back(CreateNonOverlappingBoundingBox());
            continue;
        }

        result.emplace_back(BoundingBox::CreateBoundingBox(m_nodes[e.first], m_nodes[e.second]));
    }

    return result;
}

void Mesh2D::FindFacesConnectedToNode(UInt nodeIndex, std::vector<UInt>& sharedFaces) const
{
    sharedFaces.clear();

    UInt newFaceIndex = constants::missing::uintValue;
    for (UInt e = 0; e < m_nodesNumEdges[nodeIndex]; e++)
    {
        const auto firstEdge = m_nodesEdges[nodeIndex][e];

        UInt secondEdgeIndex = e + 1;
        if (secondEdgeIndex >= m_nodesNumEdges[nodeIndex])
        {
            secondEdgeIndex = 0;
        }

        const auto secondEdge = m_nodesEdges[nodeIndex][secondEdgeIndex];
        if (m_edgesNumFaces[firstEdge] == 0 || m_edgesNumFaces[secondEdge] == 0)
        {
            continue;
        }

        // find the face shared by the two edges
        const auto firstFace = std::max(std::min<UInt>(m_edgesNumFaces[firstEdge], 2U), 1U) - 1U;
        const auto secondFace = std::max(std::min<UInt>(m_edgesNumFaces[secondEdge], 2U), 1U) - 1U;

        if (m_edgesFaces[firstEdge][0] != newFaceIndex &&
            (m_edgesFaces[firstEdge][0] == m_edgesFaces[secondEdge][0] ||
             m_edgesFaces[firstEdge][0] == m_edgesFaces[secondEdge][secondFace]))
        {
            newFaceIndex = m_edgesFaces[firstEdge][0];
        }
        else if (m_edgesFaces[firstEdge][firstFace] != newFaceIndex &&
                 (m_edgesFaces[firstEdge][firstFace] == m_edgesFaces[secondEdge][0] ||
                  m_edgesFaces[firstEdge][firstFace] == m_edgesFaces[secondEdge][secondFace]))
        {
            newFaceIndex = m_edgesFaces[firstEdge][firstFace];
        }
        else
        {
            newFaceIndex = constants::missing::uintValue;
        }

        // corner face (already found in the first iteration)
        if (m_nodesNumEdges[nodeIndex] == 2 &&
            e == 1 &&
            m_nodesTypes[nodeIndex] == MeshNodeType::Corner &&
            !sharedFaces.empty() &&
            sharedFaces[0] == newFaceIndex)
        {
            newFaceIndex = constants::missing::uintValue;
        }
        sharedFaces.emplace_back(newFaceIndex);
    }
}

void Mesh2D::GetConnectingNodes(UInt nodeIndex, std::vector<UInt>& connectedNodes) const
{
    connectedNodes.clear();
    connectedNodes.emplace_back(nodeIndex);

    // edge connected nodes
    for (UInt e = 0; e < m_nodesNumEdges[nodeIndex]; e++)
    {
        const auto edgeIndex = m_nodesEdges[nodeIndex][e];
        const auto node = OtherNodeOfEdge(m_edges[edgeIndex], nodeIndex);
        connectedNodes.emplace_back(node);
    }
}

void Mesh2D::FindNodesSharedByFaces(UInt nodeIndex, const std::vector<UInt>& sharedFaces, std::vector<UInt>& connectedNodes, std::vector<std::vector<UInt>>& faceNodeMapping) const
{

    // for each face store the positions of the its nodes in the connectedNodes (compressed array)
    if (faceNodeMapping.size() < sharedFaces.size())
    {
        ResizeAndFill2DVector(faceNodeMapping, static_cast<UInt>(sharedFaces.size()), constants::geometric::maximumNumberOfNodesPerFace);
    }

    // Find all nodes shared by faces
    for (UInt f = 0; f < sharedFaces.size(); f++)
    {
        const auto faceIndex = sharedFaces[f];

        if (faceIndex == constants::missing::uintValue)
        {
            continue;
        }

        // find the stencil node position in the current face
        UInt faceNodeIndex = GetLocalFaceNodeIndex(faceIndex, nodeIndex);
        const auto numFaceNodes = GetNumFaceEdges(faceIndex);

        if (faceNodeIndex == constants::missing::uintValue)
        {
            continue;
        }

        for (UInt n = 0; n < numFaceNodes; n++)
        {

            if (faceNodeIndex >= numFaceNodes)
            {
                faceNodeIndex -= numFaceNodes;
            }

            const auto node = m_facesNodes[faceIndex][faceNodeIndex];

            bool isNewNode = true;

            // Find if node of face is already in connectedNodes list
            for (UInt i = 0; i < connectedNodes.size(); i++)
            {
                if (node == connectedNodes[i])
                {
                    isNewNode = false;
                    faceNodeMapping[f][faceNodeIndex] = static_cast<UInt>(i);
                    break;
                }
            }

            // If node is not already contained in connectedNodes list, then add it to the list.
            if (isNewNode)
            {
                connectedNodes.emplace_back(node);
                faceNodeMapping[f][faceNodeIndex] = static_cast<UInt>(connectedNodes.size() - 1);
            }

            // update node index
            faceNodeIndex += 1;
        }
    }
}

meshkernel::UInt Mesh2D::IsStartOrEnd(const UInt edgeId, const UInt nodeId) const
{
    UInt isStartEnd = constants::missing::uintValue;

    if (m_edges[edgeId].first == nodeId)
    {
        isStartEnd = 0;
    }
    else if (m_edges[edgeId].second == nodeId)
    {
        isStartEnd = 1;
    }

    return isStartEnd;
}

meshkernel::UInt Mesh2D::IsLeftOrRight(const UInt elementId, const UInt edgeId) const
{
    UInt edgeIndex = constants::missing::uintValue;
    UInt nextEdgeIndex = constants::missing::uintValue;
    UInt endNodeIndex = m_edges[edgeId].second;

    for (UInt i = 0; i < m_facesEdges[elementId].size(); ++i)
    {
        UInt faceEdgeId = m_facesEdges[elementId][i];

        if (faceEdgeId == edgeId)
        {
            edgeIndex = i;
        }
        else if (m_edges[faceEdgeId].first == endNodeIndex || m_edges[faceEdgeId].second == endNodeIndex)
        {
            nextEdgeIndex = i;
        }
    }

    if (edgeIndex == constants::missing::uintValue || nextEdgeIndex == constants::missing::uintValue)
    {
        // EdgeId was not found
        return constants::missing::uintValue;
    }

    UInt isLeftRight = constants::missing::uintValue;

    if (nextEdgeIndex == edgeIndex + 1 || nextEdgeIndex + m_numFacesNodes[elementId] == edgeIndex + 1)
    {
        isLeftRight = 0;
    }
    else if (edgeIndex == nextEdgeIndex + 1 || edgeIndex + m_numFacesNodes[elementId] == nextEdgeIndex + 1)
    {
        isLeftRight = 1;
    }

    return isLeftRight;
}

meshkernel::UInt Mesh2D::FindCommonFace(const UInt edge1, const UInt edge2) const
{
    for (UInt i = 0; i < m_edgesNumFaces[edge1]; ++i)
    {
        for (UInt j = 0; j < m_edgesNumFaces[edge2]; ++j)
        {
            if (m_edgesFaces[edge1][i] == m_edgesFaces[edge2][j])
            {
                return m_edgesFaces[edge1][i];
            }
        }
    }

    return constants::missing::uintValue;
}
