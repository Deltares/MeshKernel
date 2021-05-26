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

#include <algorithm>
#include <array>
#include <cmath>
#include <initializer_list>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/RTree.hpp>
#include <MeshKernel/TriangulationWrapper.hpp>

using meshkernel::Mesh2D;

Mesh2D::Mesh2D(const std::vector<Edge>& edges,
               const std::vector<Point>& nodes,
               Projection projection,
               AdministrationOption administration) : Mesh(edges, nodes, projection)
{

    Administrate(administration);

    //no polygon involved, so node mask is 1 everywhere
    m_nodeMask.resize(m_nodes.size());
    std::fill(m_nodeMask.begin(), m_nodeMask.end(), 1);
};

void Mesh2D::Administrate(AdministrationOption administrationOption)
{
    AdministrateNodesEdges();

    if (administrationOption == AdministrationOption::AdministrateMeshEdges)
    {
        return;
    }

    // face administration
    m_edgesNumFaces.resize(m_edges.size());
    std::fill(m_edgesNumFaces.begin(), m_edgesNumFaces.end(), 0);

    m_edgesFaces.resize(m_edges.size());
    std::fill(m_edgesFaces.begin(), m_edgesFaces.end(), std::vector<size_t>(2, sizetMissingValue));

    m_facesMassCenters.clear();
    m_faceArea.clear();
    m_facesNodes.clear();
    m_facesEdges.clear();
    m_facesCircumcenters.clear();

    m_facesMassCenters.reserve(GetNumNodes());
    m_faceArea.reserve(GetNumNodes());
    m_facesNodes.reserve(GetNumNodes());
    m_facesEdges.reserve(GetNumNodes());
    m_facesCircumcenters.reserve(GetNumNodes());

    // find faces
    FindFaces();

    // find mesh circumcenters
    ComputeFaceCircumcentersMassCentersAndAreas();

    // classify node types
    ClassifyNodes();
}

Mesh2D::Mesh2D(const std::vector<Point>& inputNodes, const Polygons& polygons, Projection projection)
{
    m_projection = projection;
    // compute triangulation
    TriangulationWrapper triangulationWrapper;
    const auto numberOfTriangles = inputNodes.size() * 6 + 10;
    triangulationWrapper.Compute(inputNodes,
                                 TriangulationWrapper::TriangulationOptions::TriangulatePointsAndGenerateFaces,
                                 0.0,
                                 numberOfTriangles);

    // For each triangle check
    // 1. Validity of its internal angles
    // 2. Is inside the polygon
    // If so we mark the edges and we add them m_edges
    std::vector<bool> edgeNodesFlag(triangulationWrapper.m_numEdges, false);
    for (auto i = 0; i < triangulationWrapper.m_numFaces; ++i)
    {
        const auto goodTriangle = CheckTriangle(triangulationWrapper.m_faceNodes[i], inputNodes);

        if (!goodTriangle)
        {
            continue;
        }
        const Point approximateCenter = (inputNodes[triangulationWrapper.m_faceNodes[i][0]] + inputNodes[triangulationWrapper.m_faceNodes[i][1]] + inputNodes[triangulationWrapper.m_faceNodes[i][2]]) * oneThird;

        const auto isTriangleInPolygon = polygons.IsPointInPolygon(approximateCenter, 0);
        if (!isTriangleInPolygon)
        {
            continue;
        }

        // mark all edges of this triangle as good ones
        for (auto j = 0; j < numNodesInTriangle; ++j)
        {
            edgeNodesFlag[triangulationWrapper.m_faceEdges[i][j]] = true;
        }
    }

    // now add all points and all valid edges
    m_nodes = inputNodes;
    size_t validEdgesCount = 0;
    for (auto i = 0; i < triangulationWrapper.m_numEdges; ++i)
    {
        if (!edgeNodesFlag[i])
            continue;
        validEdgesCount++;
    }

    std::vector<Edge> edges(validEdgesCount);
    validEdgesCount = 0;
    for (auto i = 0; i < triangulationWrapper.m_numEdges; ++i)
    {
        if (!edgeNodesFlag[i])
            continue;

        edges[validEdgesCount].first = triangulationWrapper.m_edgeNodes[i][0];
        edges[validEdgesCount].second = triangulationWrapper.m_edgeNodes[i][1];
        validEdgesCount++;
    }

    m_nodesRTreeRequiresUpdate = true;
    m_edgesRTreeRequiresUpdate = true;

    *this = Mesh2D(edges, inputNodes, projection, AdministrationOption::AdministrateMeshEdges);
}

bool Mesh2D::CheckTriangle(const std::vector<size_t>& faceNodes, const std::vector<Point>& nodes) const
{
    // Used for triangular grids
    constexpr double triangleMinimumAngle = 5.0;
    constexpr double triangleMaximumAngle = 150.0;

    double phiMin = 1e3;
    double phiMax = 0.0;

    static std::array<std::array<size_t, 3>, 3> nodePermutations{{{2, 0, 1}, {0, 1, 2}, {1, 2, 0}}};

    for (auto i = 0; i < faceNodes.size(); ++i)
    {
        Point x0 = nodes[faceNodes[nodePermutations[i][0]]];
        Point x1 = nodes[faceNodes[nodePermutations[i][1]]];
        Point x2 = nodes[faceNodes[nodePermutations[i][2]]];

        const auto cosphi = NormalizedInnerProductTwoSegments(x1, x0, x1, x2, m_projection);
        const auto phi = std::acos(std::min(std::max(cosphi, -1.0), 1.0)) * raddeg_hp;
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
    Administrate(AdministrationOption::AdministrateMeshEdgesAndFaces);

    // assume the max amount of degenerated triangles is 10% of the actual faces
    std::vector<size_t> degeneratedTriangles;
    degeneratedTriangles.reserve(static_cast<size_t>(static_cast<double>(GetNumFaces()) * 0.1));
    for (auto f = 0; f < GetNumFaces(); ++f)
    {
        const auto numFaceNodes = m_numFacesNodes[f];
        if (numFaceNodes != numNodesInTriangle)
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
            for (auto e = 0; e < numNodesInTriangle; ++e)
            {
                const auto edge = m_facesEdges[f][e];
                m_edges[edge] = {sizetMissingValue, sizetMissingValue};
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

    Administrate(AdministrationOption::AdministrateMeshEdgesAndFaces);
}

void Mesh2D::FindFacesRecursive(size_t startingNode,
                                size_t node,
                                size_t index,
                                size_t previousEdge,
                                size_t numClosingEdges,
                                std::vector<size_t>& edges,
                                std::vector<size_t>& nodes,
                                std::vector<size_t>& sortedEdgesFaces,
                                std::vector<size_t>& sortedNodes,
                                std::vector<Point>& nodalValues)
{
    // The selected edge does not exist.
    if (index >= numClosingEdges)
        return;

    if (m_edges[previousEdge].first == sizetMissingValue || m_edges[previousEdge].second == sizetMissingValue)
        throw std::invalid_argument("Mesh2D::FindFacesRecursive: The selected edge is invalid. This should not happen since all invalid edges should have been cleaned up.");

    // Check if the faces are already found
    if (m_edgesNumFaces[previousEdge] >= 2)
        return;

    edges[index] = previousEdge;
    nodes[index] = node;
    const auto otherNode = OtherNodeOfEdge(m_edges[previousEdge], node);

    // enclosure found
    if (otherNode == startingNode && index == edges.size() - 1)
    {
        // no duplicated nodes allowed
        sortedNodes = nodes;
        std::sort(sortedNodes.begin(), sortedNodes.end());
        for (auto n = 0; n < sortedNodes.size() - 1; n++)
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
            // is an internal face only if all edges have a different face
            for (auto ee = 0; ee < edges.size(); ee++)
            {
                sortedEdgesFaces[ee] = m_edgesFaces[edges[ee]][0];
            }
            std::sort(sortedEdgesFaces.begin(), sortedEdgesFaces.end());
            for (auto n = 0; n < sortedEdgesFaces.size() - 1; n++)
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

        double area;
        Point centerOfMass;
        bool isCounterClockWise;
        FaceAreaAndCenterOfMass(nodalValues, m_projection, area, centerOfMass, isCounterClockWise);
        if (!isCounterClockWise)
        {
            return;
        }

        // increase m_edgesNumFaces
        for (const auto& edge : edges)
        {
            m_edgesNumFaces[edge] += 1;
            const auto numFace = m_edgesNumFaces[edge];
            m_edgesFaces[edge][numFace - 1] = GetNumFaces();
        }

        // store the result
        m_facesNodes.emplace_back(nodes);
        m_facesEdges.emplace_back(edges);
        m_faceArea.emplace_back(area);
        m_facesMassCenters.emplace_back(centerOfMass);

        return;
    }

    size_t edgeIndexOtherNode = 0;
    for (auto e = 0; e < m_nodesNumEdges[otherNode]; e++)
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

    const auto edge = m_nodesEdges[otherNode][edgeIndexOtherNode];
    FindFacesRecursive(startingNode, otherNode, index + 1, edge, numClosingEdges, edges, nodes, sortedEdgesFaces, sortedNodes, nodalValues);
}

void Mesh2D::FindFaces()
{
    for (auto numEdgesPerFace = 3; numEdgesPerFace <= maximumNumberOfEdgesPerFace; numEdgesPerFace++)
    {
        std::vector<size_t> edges(numEdgesPerFace);
        std::vector<size_t> nodes(numEdgesPerFace);
        std::vector<size_t> sortedEdgesFaces(numEdgesPerFace);
        std::vector<size_t> sortedNodes(numEdgesPerFace);
        std::vector<Point> nodalValues;
        nodalValues.reserve(maximumNumberOfEdgesPerFace);
        for (auto n = 0; n < GetNumNodes(); n++)
        {
            if (!m_nodes[n].IsValid())
            {
                continue;
            }

            for (auto e = 0; e < m_nodesNumEdges[n]; e++)
            {
                FindFacesRecursive(n, n, 0, m_nodesEdges[n][e], numEdgesPerFace, edges, nodes, sortedEdgesFaces, sortedNodes, nodalValues);
            }
        }
    }

    m_numFacesNodes.resize(GetNumFaces());
    for (auto f = 0; f < GetNumFaces(); ++f)
    {
        m_numFacesNodes[f] = m_facesNodes[f].size();
    }
}

void Mesh2D::ComputeFaceCircumcentersMassCentersAndAreas(bool computeMassCenters)
{
    m_facesCircumcenters.resize(GetNumFaces());
    m_faceArea.resize(GetNumFaces());
    m_facesMassCenters.resize(GetNumFaces());

    std::vector<size_t> numEdgeFacesCache;
    numEdgeFacesCache.reserve(maximumNumberOfEdgesPerFace);
    for (auto f = 0; f < GetNumFaces(); f++)
    {
        //need to account for spherical coordinates. Build a polygon around a face
        ComputeFaceClosedPolygon(f, m_polygonNodesCache);

        if (computeMassCenters)
        {
            double area;
            Point centerOfMass;
            bool isCounterClockWise;
            FaceAreaAndCenterOfMass(m_polygonNodesCache, m_projection, area, centerOfMass, isCounterClockWise);
            m_faceArea[f] = area;
            m_facesMassCenters[f] = centerOfMass;
        }

        size_t numberOfInteriorEdges = 0;
        const auto numberOfFaceNodes = GetNumFaceEdges(f);
        for (auto n = 0; n < numberOfFaceNodes; n++)
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
        for (auto n = 0; n < numberOfFaceNodes; n++)
        {
            numEdgeFacesCache.emplace_back(m_edgesNumFaces[m_facesEdges[f][n]]);
        }

        m_facesCircumcenters[f] = ComputeFaceCircumenter(m_polygonNodesCache,
                                                         numEdgeFacesCache);
    }
}

void Mesh2D::ClassifyNodes()
{
    m_nodesTypes.resize(GetNumNodes(), 0);
    std::fill(m_nodesTypes.begin(), m_nodesTypes.end(), 0);

    for (auto e = 0; e < GetNumEdges(); e++)
    {
        const auto firstNode = m_edges[e].first;
        const auto secondNode = m_edges[e].second;

        if (firstNode == sizetMissingValue || secondNode == sizetMissingValue)
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

    for (auto n = 0; n < GetNumNodes(); n++)
    {
        if (m_nodesTypes[n] == 1 || m_nodesTypes[n] == 2)
        {
            if (m_nodesNumEdges[n] == 2)
            {
                //corner point
                m_nodesTypes[n] = 3;
            }
            else
            {
                size_t firstNode = sizetMissingValue;
                size_t secondNode = sizetMissingValue;
                for (auto i = 0; i < m_nodesNumEdges[n]; i++)
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
                if (firstNode != sizetMissingValue && secondNode != sizetMissingValue)
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
            //internal node
            m_nodesTypes[n] = 1;
        }
        if (m_nodesNumEdges[n] < 2)
        {
            //hanging node
            m_nodesTypes[n] = -1;
        }
    }
}

void Mesh2D::ComputeFaceClosedPolygonWithLocalMappings(size_t faceIndex,
                                                       std::vector<Point>& polygonNodesCache,
                                                       std::vector<size_t>& localNodeIndicesCache,
                                                       std::vector<size_t>& globalEdgeIndicesCache) const
{
    const auto numFaceNodes = GetNumFaceEdges(faceIndex);
    polygonNodesCache.reserve(numFaceNodes + 1);
    polygonNodesCache.clear();
    localNodeIndicesCache.reserve(numFaceNodes + 1);
    localNodeIndicesCache.clear();
    globalEdgeIndicesCache.reserve(numFaceNodes + 1);
    globalEdgeIndicesCache.clear();

    for (auto n = 0; n < numFaceNodes; n++)
    {
        polygonNodesCache.emplace_back(m_nodes[m_facesNodes[faceIndex][n]]);
        localNodeIndicesCache.emplace_back(n);
        globalEdgeIndicesCache.emplace_back(m_facesEdges[faceIndex][n]);
    }
    polygonNodesCache.emplace_back(polygonNodesCache.front());
    localNodeIndicesCache.emplace_back(0);
    globalEdgeIndicesCache.emplace_back(globalEdgeIndicesCache.front());
}

void Mesh2D::ComputeFaceClosedPolygon(size_t faceIndex, std::vector<Point>& polygonNodesCache) const
{
    const auto numFaceNodes = GetNumFaceEdges(faceIndex);
    polygonNodesCache.clear();
    polygonNodesCache.reserve(numFaceNodes);
    for (auto n = 0; n < numFaceNodes; n++)
    {
        polygonNodesCache.push_back(m_nodes[m_facesNodes[faceIndex][n]]);
    }
    polygonNodesCache.push_back(polygonNodesCache.front());
}

void Mesh2D::MaskNodesInPolygons(const Polygons& polygon, bool inside)
{
    std::fill(m_nodeMask.begin(), m_nodeMask.end(), 0);
    const auto nodePolygonIndices = polygon.PolygonIndices(m_nodes);

    for (auto i = 0; i < GetNumNodes(); ++i)
    {
        auto isInPolygon = nodePolygonIndices[i] == sizetMissingValue ? false : true;
        if (!inside)
        {
            isInPolygon = !isInPolygon;
        }
        m_nodeMask[i] = 0;
        if (isInPolygon)
        {
            m_nodeMask[i] = 1;
        }
    }
}

bool Mesh2D::IsFullFaceNotInPolygon(size_t faceIndex) const
{
    for (auto n = 0; n < GetNumFaceEdges(faceIndex); n++)
    {
        if (m_nodeMask[m_facesNodes[faceIndex][n]] != 1)
        {
            return true;
        }
    }
    return false;
}

void Mesh2D::OffsetSphericalCoordinates(double minx, double maxx)
{
    if (m_projection == Projection::spherical && maxx - minx > 180.0)
    {
        for (auto n = 0; n < GetNumNodes(); ++n)
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

void Mesh2D::MaskFaceEdgesInPolygon(const Polygons& polygons, bool invertSelection, bool includeIntersected)
{
    Administrate(AdministrationOption::AdministrateMeshEdgesAndFaces);

    // mark all nodes in polygon with 1
    std::fill(m_nodeMask.begin(), m_nodeMask.end(), 0);
    for (auto n = 0; n < GetNumNodes(); ++n)
    {
        const auto isInPolygon = polygons.IsPointInPolygon(m_nodes[n], 0);
        if (isInPolygon)
        {
            m_nodeMask[n] = 1;
        }
    }

    // mark all edges with both start end end nodes included with 1
    std::vector<int> edgeMask(m_edges.size(), 0);
    for (auto e = 0; e < GetNumEdges(); ++e)
    {
        const auto firstNodeIndex = m_edges[e].first;
        const auto secondNodeIndex = m_edges[e].second;

        int isEdgeIncluded;
        if (includeIntersected)
        {
            isEdgeIncluded = (firstNodeIndex != sizetMissingValue && m_nodeMask[firstNodeIndex] == 1 ||
                              secondNodeIndex != sizetMissingValue && m_nodeMask[secondNodeIndex] == 1)
                                 ? 1
                                 : 0;
        }
        else
        {
            isEdgeIncluded = (firstNodeIndex != sizetMissingValue && m_nodeMask[firstNodeIndex] == 1 &&
                              secondNodeIndex != sizetMissingValue && m_nodeMask[secondNodeIndex] == 1)
                                 ? 1
                                 : 0;
        }

        edgeMask[e] = isEdgeIncluded;
    }

    // if one edge of the face is not included do not include all the edges of that face
    auto secondEdgeMask = edgeMask;
    if (!includeIntersected)
    {
        for (auto f = 0; f < GetNumFaces(); ++f)
        {
            bool isOneEdgeNotIncluded = false;
            for (auto n = 0; n < GetNumFaceEdges(f); ++n)
            {
                const auto edgeIndex = m_facesEdges[f][n];
                if (edgeIndex != sizetMissingValue && edgeMask[edgeIndex] == 0)
                {
                    isOneEdgeNotIncluded = true;
                    break;
                }
            }

            if (isOneEdgeNotIncluded)
            {
                for (auto n = 0; n < GetNumFaceEdges(f); ++n)
                {
                    const auto edgeIndex = m_facesEdges[f][n];
                    if (edgeIndex != sizetMissingValue)
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
        for (auto e = 0; e < GetNumEdges(); ++e)
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

    m_edgeMask = std::move(secondEdgeMask);
}

void Mesh2D::ComputeNodeMaskFromEdgeMask()
{
    if (m_edgeMask.size() != GetNumEdges() || m_nodeMask.size() != GetNumNodes())
    {
        throw std::invalid_argument("Mesh2D::ComputeNodeMaskFromEdgeMask:The dimension of the masks do not fit the mesh.");
    }

    // fill node mask to 0
    std::fill(m_nodeMask.begin(), m_nodeMask.end(), 0);

    // compute node mask from edge mask
    for (auto e = 0; e < GetNumEdges(); ++e)
    {
        if (m_edgeMask[e] != 1)
            continue;

        const auto firstNodeIndex = m_edges[e].first;
        const auto secondNodeIndex = m_edges[e].second;

        if (firstNodeIndex != sizetMissingValue)
        {
            m_nodeMask[firstNodeIndex] = 1;
        }
        if (secondNodeIndex != sizetMissingValue)
        {
            m_nodeMask[secondNodeIndex] = 1;
        }
    }
}

meshkernel::Point Mesh2D::ComputeFaceCircumenter(std::vector<Point>& polygon,
                                                 const std::vector<size_t>& edgesNumFaces) const
{
    const size_t maximumNumberCircumcenterIterations = 100;
    const double eps = m_projection == Projection::cartesian ? 1e-3 : 9e-10; //111km = 0-e digit.
    std::vector<Point> middlePoints;
    middlePoints.reserve(maximumNumberOfNodesPerFace);
    std::vector<Point> normals;
    normals.reserve(maximumNumberOfNodesPerFace);
    const auto numNodes = polygon.size() - 1;

    Point centerOfMass{0.0, 0.0};
    for (auto n = 0; n < numNodes; n++)
    {
        centerOfMass.x += polygon[n].x;
        centerOfMass.y += polygon[n].y;
    }
    centerOfMass = centerOfMass / static_cast<double>(numNodes);

    auto result = centerOfMass;
    if (numNodes == numNodesInTriangle)
    {
        result = CircumcenterOfTriangle(polygon[0], polygon[1], polygon[2], m_projection);
    }
    else if (!edgesNumFaces.empty())
    {
        size_t numValidEdges = 0;
        for (auto n = 0; n < numNodes; ++n)
        {
            if (edgesNumFaces[n] == 2)
            {
                numValidEdges++;
            }
        }

        if (numValidEdges > 1)
        {
            for (auto n = 0; n < numNodes; n++)
            {
                if (edgesNumFaces[n] != 2)
                {
                    continue;
                }
                const auto nextNode = NextCircularForwardIndex(n, numNodes);
                middlePoints.emplace_back((polygon[n] + polygon[nextNode]) * 0.5);
                normals.emplace_back(NormalVector(polygon[n], polygon[nextNode], middlePoints.back(), m_projection));
            }

            Point estimatedCircumCenter = centerOfMass;
            for (auto iter = 0; iter < maximumNumberCircumcenterIterations; ++iter)
            {
                const Point previousCircumCenter = estimatedCircumCenter;
                for (auto n = 0; n < middlePoints.size(); n++)
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

    for (auto n = 0; n < numNodes; n++)
    {
        polygon[n].x = weightCircumCenter * polygon[n].x + (1.0 - weightCircumCenter) * centerOfMass.x;
        polygon[n].y = weightCircumCenter * polygon[n].y + (1.0 - weightCircumCenter) * centerOfMass.y;
    }

    if (IsPointInPolygonNodes(result, polygon, m_projection))
    {
        return result;
    }

    for (auto n = 0; n < numNodes; n++)
    {
        const auto nextNode = NextCircularForwardIndex(n, numNodes);
        Point intersection;
        double crossProduct;
        double firstRatio;
        double secondRatio;
        const auto areLineCrossing = AreSegmentsCrossing(centerOfMass, result, polygon[n], polygon[nextNode], false, m_projection, intersection, crossProduct, firstRatio, secondRatio);
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
    Administrate(AdministrationOption::AdministrateMeshEdgesAndFaces);
    std::vector<Point> result;
    result.reserve(GetNumFaces());
    for (auto f = 0; f < GetNumFaces(); ++f)
    {
        // a triangle
        if (m_numFacesNodes[f] == 3)
        {
            const auto firstNode = m_facesNodes[f][0];
            const auto secondNode = m_facesNodes[f][1];
            const auto thirdNode = m_facesNodes[f][2];
            //compute squared edge lengths
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

std::vector<size_t> Mesh2D::GetEdgesCrossingSmallFlowEdges(double smallFlowEdgesThreshold)
{
    Administrate(AdministrationOption::AdministrateMeshEdgesAndFaces);
    std::vector<size_t> result;
    result.reserve(GetNumEdges());
    for (auto e = 0; e < GetNumEdges(); ++e)
    {
        const auto firstFace = m_edgesFaces[e][0];
        const auto secondFace = m_edgesFaces[e][1];

        if (firstFace != sizetMissingValue && secondFace != sizetMissingValue)
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

std::vector<meshkernel::Point> Mesh2D::GetFlowEdgesCenters(const std::vector<size_t>& edges) const
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
            m_edges[e] = {sizetMissingValue, sizetMissingValue};
        }
        Administrate(AdministrationOption::AdministrateMeshEdgesAndFaces);
    }
}

void Mesh2D::DeleteSmallTrianglesAtBoundaries(double minFractionalAreaTriangles)
{
    // On the second part, the small triangles at the boundaries are checked
    const double minCosPhi = 0.2;
    std::vector<std::vector<size_t>> smallTrianglesNodes;
    for (auto face = 0; face < GetNumFaces(); ++face)
    {
        if (m_numFacesNodes[face] != numNodesInTriangle || m_faceArea[face] <= 0.0 || !IsFaceOnBoundary(face))
        {
            continue;
        }

        // compute the average area of neighboring faces
        double averageOtherFacesArea = 0.0;
        size_t numNonBoundaryFaces = 0;
        for (auto e = 0; e < numNodesInTriangle; ++e)
        {
            // the edge must not be at the boundary, otherwise there is no "other" face
            const auto edge = m_facesEdges[face][e];
            if (IsEdgeOnBoundary(edge))
            {
                continue;
            }
            const auto otherFace = face == m_edgesFaces[edge][0] ? m_edgesFaces[edge][1] : m_edgesFaces[edge][0];
            if (m_numFacesNodes[otherFace] > numNodesInTriangle)
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
        size_t nodeToPreserve = sizetMissingValue;
        size_t firstNodeToMerge;
        size_t secondNodeToMerge;
        size_t thirdEdgeSmallTriangle = sizetMissingValue;
        for (auto e = 0; e < numNodesInTriangle; ++e)
        {
            const auto previousEdge = NextCircularBackwardIndex(e, numNodesInTriangle);
            const auto nextEdge = NextCircularForwardIndex(e, numNodesInTriangle);

            const auto k0 = m_facesNodes[face][previousEdge];
            const auto k1 = m_facesNodes[face][e];
            const auto k2 = m_facesNodes[face][nextEdge];

            // compute the angles between the edges
            const auto cosphi = std::abs(NormalizedInnerProductTwoSegments(m_nodes[k0], m_nodes[k1], m_nodes[k1], m_nodes[k2], m_projection));

            if (cosphi < minCosPhiSmallTriangle)
            {
                minCosPhiSmallTriangle = cosphi;
                firstNodeToMerge = k0;
                nodeToPreserve = k1;
                secondNodeToMerge = k2;
                thirdEdgeSmallTriangle = m_facesEdges[face][nextEdge];
            }
        }

        if (minCosPhiSmallTriangle < minCosPhi && thirdEdgeSmallTriangle != sizetMissingValue && IsEdgeOnBoundary(thirdEdgeSmallTriangle))
        {
            smallTrianglesNodes.emplace_back(std::initializer_list<size_t>{nodeToPreserve, firstNodeToMerge, secondNodeToMerge});
        }
    }

    bool nodesMerged = false;
    for (const auto& triangleNodes : smallTrianglesNodes)
    {
        const auto nodeToPreserve = triangleNodes[0];
        const auto firstNodeToMerge = triangleNodes[1];
        const auto secondNodeToMerge = triangleNodes[2];

        // only
        size_t numInternalEdges = 0;
        for (auto e = 0; e < m_nodesNumEdges[firstNodeToMerge]; ++e)
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
        for (auto e = 0; e < m_nodesNumEdges[secondNodeToMerge]; ++e)
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
        Administrate(AdministrationOption::AdministrateMeshEdgesAndFaces);
    }
}

void Mesh2D::ComputeNodeNeighbours()
{
    m_maxNumNeighbours = *(std::max_element(m_nodesNumEdges.begin(), m_nodesNumEdges.end()));
    m_maxNumNeighbours += 1;

    m_nodesNodes.resize(GetNumNodes(), std::vector<size_t>(m_maxNumNeighbours, sizetMissingValue));
    //for each node, determine the neighbouring nodes
    for (auto n = 0; n < GetNumNodes(); n++)
    {
        for (auto nn = 0; nn < m_nodesNumEdges[n]; nn++)
        {
            const auto edge = m_edges[m_nodesEdges[n][nn]];
            m_nodesNodes[n][nn] = OtherNodeOfEdge(edge, n);
        }
    }
}

std::vector<double> Mesh2D::GetOrthogonality()
{
    std::vector<double> result;
    result.reserve(GetNumEdges());
    for (auto e = 0; e < GetNumEdges(); e++)
    {
        auto val = doubleMissingValue;
        const auto firstNode = m_edges[e].first;
        const auto secondNode = m_edges[e].second;

        if (firstNode != sizetMissingValue && secondNode != sizetMissingValue && !IsEdgeOnBoundary(e))
        {
            val = NormalizedInnerProductTwoSegments(m_nodes[firstNode],
                                                    m_nodes[secondNode],
                                                    m_facesCircumcenters[m_edgesFaces[e][0]],
                                                    m_facesCircumcenters[m_edgesFaces[e][1]],
                                                    m_projection);
            if (val != doubleMissingValue)
            {
                val = std::abs(val);
            }
        }
        result.emplace_back(val);
    }
    return result;
}

std::vector<double> Mesh2D::GetSmoothness()
{
    std::vector<double> result;
    result.reserve(GetNumEdges());
    for (auto e = 0; e < GetNumEdges(); e++)
    {
        auto val = doubleMissingValue;
        const auto firstNode = m_edges[e].first;
        const auto secondNode = m_edges[e].second;

        if (firstNode != sizetMissingValue && secondNode != sizetMissingValue && !IsEdgeOnBoundary(e))
        {
            const auto leftFace = m_edgesFaces[e][0];
            const auto rightFace = m_edgesFaces[e][1];
            const auto leftFaceArea = m_faceArea[leftFace];
            const auto rightFaceArea = m_faceArea[rightFace];

            if (leftFaceArea < minimumCellArea || rightFaceArea < minimumCellArea)
            {
                val = rightFaceArea / leftFaceArea;
            }
            if (val < 1.0)
            {
                val = 1.0 / val;
            }
        }
        result.emplace_back(val);
    }
    return result;
}

void Mesh2D::ComputeAspectRatios(std::vector<double>& aspectRatios)
{
    std::vector<std::vector<double>> averageEdgesLength(GetNumEdges(), std::vector<double>(2, doubleMissingValue));
    std::vector<double> averageFlowEdgesLength(GetNumEdges(), doubleMissingValue);
    std::vector<bool> curvilinearGridIndicator(GetNumNodes(), true);
    std::vector<double> edgesLength(GetNumEdges(), 0.0);
    aspectRatios.resize(GetNumEdges(), 0.0);

    for (auto e = 0; e < GetNumEdges(); e++)
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

        //find right cell center, if it exists
        if (m_edgesNumFaces[e] == 2)
        {
            rightCenter = m_facesCircumcenters[m_edgesFaces[e][1]];
        }
        else
        {
            //otherwise, make ghost node by imposing boundary condition
            double dinry = InnerProductTwoSegments(m_nodes[first], m_nodes[second], m_nodes[first], leftCenter, m_projection);
            dinry = dinry / std::max(edgeLength * edgeLength, minimumEdgeLength);

            const double x0_bc = (1.0 - dinry) * m_nodes[first].x + dinry * m_nodes[second].x;
            const double y0_bc = (1.0 - dinry) * m_nodes[first].y + dinry * m_nodes[second].y;
            rightCenter.x = 2.0 * x0_bc - leftCenter.x;
            rightCenter.y = 2.0 * y0_bc - leftCenter.y;
        }

        averageFlowEdgesLength[e] = ComputeDistance(leftCenter, rightCenter, m_projection);
    }

    // Compute normal length
    for (auto f = 0; f < GetNumFaces(); f++)
    {
        const auto numberOfFaceNodes = GetNumFaceEdges(f);
        if (numberOfFaceNodes < numNodesInTriangle)
            continue;

        for (auto n = 0; n < numberOfFaceNodes; n++)
        {
            if (numberOfFaceNodes != numNodesQuads)
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

            //quads
            if (numberOfFaceNodes == numNodesQuads)
            {
                size_t kkp2 = n + 2;
                if (kkp2 >= numberOfFaceNodes)
                {
                    kkp2 = kkp2 - numberOfFaceNodes;
                }

                const auto klinkp2 = m_facesEdges[f][kkp2];
                edgeLength = 0.5 * (edgesLength[edgeIndex] + edgesLength[klinkp2]);
            }

            if (IsEqual(averageEdgesLength[edgeIndex][0], doubleMissingValue))
            {
                averageEdgesLength[edgeIndex][0] = edgeLength;
            }
            else
            {
                averageEdgesLength[edgeIndex][1] = edgeLength;
            }
        }
    }

    if (curvilinearToOrthogonalRatio == 1.0)
        return;

    for (auto e = 0; e < GetNumEdges(); e++)
    {
        const auto first = m_edges[e].first;
        const auto second = m_edges[e].second;

        if (first == sizetMissingValue || second == sizetMissingValue)
            continue;
        if (m_edgesNumFaces[e] == 0)
            continue;
        // Consider only quads
        if (!curvilinearGridIndicator[first] || !curvilinearGridIndicator[second])
            continue;

        if (IsEdgeOnBoundary(e))
        {
            if (averageEdgesLength[e][0] > 0.0 &&
                IsEqual(averageEdgesLength[e][0], doubleMissingValue))
            {
                aspectRatios[e] = averageFlowEdgesLength[e] / averageEdgesLength[e][0];
            }
        }
        else
        {
            if (averageEdgesLength[e][0] > 0.0 &&
                averageEdgesLength[e][1] > 0.0 &&
                IsEqual(averageEdgesLength[e][0], doubleMissingValue) &&
                IsEqual(averageEdgesLength[e][1], doubleMissingValue))
            {
                aspectRatios[e] = curvilinearToOrthogonalRatio * aspectRatios[e] +
                                  (1.0 - curvilinearToOrthogonalRatio) * averageFlowEdgesLength[e] / (0.5 * (averageEdgesLength[e][0] + averageEdgesLength[e][1]));
            }
        }
    }
}

void Mesh2D::TriangulateFaces()
{
    for (auto i = 0; i < GetNumFaces(); i++)
    {
        const auto NumEdges = GetNumFaceEdges(i);

        if (NumEdges < 4)
        {
            continue;
        }

        const auto indexFirstNode = m_facesNodes[i][0];
        for (auto j = 2; j < NumEdges - 1; j++)
        {
            const auto nodeIndex = m_facesNodes[i][j];
            ConnectNodes(indexFirstNode, nodeIndex);
        }
    }

    m_edgesRTreeRequiresUpdate = true;
}

void Mesh2D::MakeDualFace(size_t node, double enlargementFactor, std::vector<Point>& dualFace)
{
    const auto sortedFacesIndices = SortedFacesAroundNode(node);
    const auto numEdges = m_nodesNumEdges[node];
    dualFace.reserve(maximumNumberOfEdgesPerNode);
    dualFace.clear();

    for (auto e = 0; e < numEdges; ++e)
    {
        const auto edgeIndex = m_nodesEdges[node][e];
        auto edgeCenter = m_edgesCenters[edgeIndex];

        if (m_projection == Projection::spherical)
        {
            const auto firstNodeIndex = m_edges[edgeIndex].first;
            const auto secondNodeIndex = m_edges[edgeIndex].second;

            if (firstNodeIndex != sizetMissingValue && secondNodeIndex != sizetMissingValue)
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
        if (faceIndex != sizetMissingValue)
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
    double area;
    Point centerOfMass;
    bool isCounterClockWise;
    FaceAreaAndCenterOfMass(dualFace, m_projection, area, centerOfMass, isCounterClockWise);

    if (m_projection == Projection::spherical)
    {
        if (centerOfMass.x - m_nodes[node].x > 180.0)
        {
            centerOfMass.x -= 360.0;
        }
        if (centerOfMass.x - m_nodes[node].x < 180.0)
        {
            centerOfMass.x += 360.0;
        }
    }

    for (auto& v : dualFace)
    {
        v = centerOfMass + (v - centerOfMass) * enlargementFactor;
    }
}

std::vector<size_t> Mesh2D::SortedFacesAroundNode(size_t node) const
{

    const auto numEdges = m_nodesNumEdges[node];
    std::vector<size_t> result;
    for (auto e = 0; e < numEdges; ++e)
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

        size_t secondFace = sizetMissingValue;
        if (m_edgesNumFaces[firstEdge] > 1)
        {
            secondFace = m_edgesFaces[firstEdge][1];
        }

        // check if the first face contains the first edge
        size_t firstEdgeIndexInFirstFace = 0;
        for (auto n = 0; n < m_numFacesNodes[firstFace]; ++n)
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

std::vector<meshkernel::Point> Mesh2D::MeshBoundaryToPolygon(const std::vector<Point>& polygon)
{

    // Find faces
    Administrate(AdministrationOption::AdministrateMeshEdgesAndFaces);
    std::vector<bool> isVisited(GetNumEdges(), false);
    std::vector<Point> meshBoundaryPolygon;
    meshBoundaryPolygon.reserve(GetNumNodes());

    for (auto e = 0; e < GetNumEdges(); e++)
    {
        if (isVisited[e] || !IsEdgeOnBoundary(e))
        {
            continue;
        }

        const auto firstNodeIndex = m_edges[e].first;
        const auto secondNodeIndex = m_edges[e].second;
        const auto firstNode = m_nodes[firstNodeIndex];
        const auto secondNode = m_nodes[secondNodeIndex];

        const auto firstNodeInPolygon = IsPointInPolygonNodes(m_nodes[firstNodeIndex], polygon, m_projection);
        const auto secondNodeInPolygon = IsPointInPolygonNodes(m_nodes[secondNodeIndex], polygon, m_projection);

        if (!firstNodeInPolygon && !secondNodeInPolygon)
        {
            continue;
        }

        //Start a new polyline
        if (!meshBoundaryPolygon.empty())
        {
            meshBoundaryPolygon.emplace_back(doubleMissingValue, doubleMissingValue);
        }

        // Put the current edge on the mesh boundary, mark it as visited
        const auto startPolygonEdges = meshBoundaryPolygon.size();
        meshBoundaryPolygon.emplace_back(firstNode);
        meshBoundaryPolygon.emplace_back(secondNode);
        isVisited[e] = true;

        // walk the current mesh boundary
        auto currentNode = secondNodeIndex;
        WalkBoundaryFromNode(polygon, isVisited, currentNode, meshBoundaryPolygon);

        const auto numNodesFirstTail = meshBoundaryPolygon.size();

        // if the boundary polygon is not closed
        if (currentNode != firstNodeIndex)
        {
            //Now grow a polyline starting at the other side of the original link L, i.e., the second tail
            currentNode = firstNodeIndex;
            WalkBoundaryFromNode(polygon, isVisited, currentNode, meshBoundaryPolygon);
        }

        // There is a nonempty second tail, so reverse the first tail, so that they connect.
        if (meshBoundaryPolygon.size() > numNodesFirstTail)
        {
            const auto start = startPolygonEdges + static_cast<size_t>(std::ceil((numNodesFirstTail - startPolygonEdges + static_cast<size_t>(1)) * 0.5));
            for (auto n = start; n < numNodesFirstTail; n++)
            {
                const auto backupPoint = meshBoundaryPolygon[n];
                const auto replaceIndex = numNodesFirstTail - n + firstNodeIndex;
                meshBoundaryPolygon[n] = meshBoundaryPolygon[replaceIndex];
                meshBoundaryPolygon[replaceIndex] = backupPoint;
            }
        }

        //Start a new polyline
        meshBoundaryPolygon.emplace_back(doubleMissingValue, doubleMissingValue);
    }
    return meshBoundaryPolygon;
}

void Mesh2D::WalkBoundaryFromNode(const std::vector<Point>& polygonNodes,
                                  std::vector<bool>& isVisited,
                                  size_t& currentNode,
                                  std::vector<Point>& meshBoundaryPolygon) const
{
    size_t e = 0;
    bool currentNodeInPolygon = false;
    while (e < m_nodesNumEdges[currentNode])
    {
        if (!currentNodeInPolygon)
        {
            currentNodeInPolygon = IsPointInPolygonNodes(m_nodes[currentNode], polygonNodes, m_projection);
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

std::vector<size_t> Mesh2D::GetHangingEdges() const
{
    std::vector<size_t> result;
    for (auto e = 0; e < GetNumEdges(); e++)
    {
        const auto firstNode = m_edges[e].first;
        const auto secondNode = m_edges[e].second;

        if (firstNode != sizetMissingValue && secondNode != sizetMissingValue)
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

void Mesh2D::DeleteMesh(const Polygons& polygons, int deletionOption, bool invertDeletion)
{
    if (deletionOption == AllNodesInside)
    {
        for (auto n = 0; n < GetNumNodes(); ++n)
        {
            auto isInPolygon = polygons.IsPointInPolygon(m_nodes[n], 0);
            if (invertDeletion)
            {
                isInPolygon = !isInPolygon;
            }
            if (isInPolygon)
            {
                m_nodes[n] = {doubleMissingValue, doubleMissingValue};
            }
        }
    }

    if (deletionOption == FacesWithIncludedCircumcenters)
    {
        Administrate(AdministrationOption::AdministrateMeshEdgesAndFaces);

        for (auto e = 0; e < GetNumEdges(); ++e)
        {
            bool allFaceCircumcentersInPolygon = true;

            for (auto f = 0; f < GetNumEdgesFaces(e); ++f)
            {
                const auto faceIndex = m_edgesFaces[e][f];
                if (faceIndex == sizetMissingValue)
                {
                    continue;
                }

                const auto faceCircumcenter = m_facesCircumcenters[faceIndex];
                auto isInPolygon = polygons.IsPointInPolygon(faceCircumcenter, 0);
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

                if (firstNodeIndex == sizetMissingValue || secondNodeIndex == sizetMissingValue)
                {
                    continue;
                }

                const auto edgeCenter = (m_nodes[firstNodeIndex] + m_nodes[secondNodeIndex]) / 2.0;

                allFaceCircumcentersInPolygon = polygons.IsPointInPolygon(edgeCenter, 0);
                if (invertDeletion)
                {
                    allFaceCircumcentersInPolygon = !allFaceCircumcentersInPolygon;
                }
            }

            if (allFaceCircumcentersInPolygon)
            {
                m_edges[e].first = sizetMissingValue;
                m_edges[e].second = sizetMissingValue;
            }
        }
    }

    if (deletionOption == FacesCompletelyIncluded)
    {
        MaskFaceEdgesInPolygon(polygons, invertDeletion, false);

        // mark the edges for deletion
        for (auto e = 0; e < GetNumEdges(); ++e)
        {
            if (m_edgeMask[e] == 1)
            {
                m_edges[e].first = sizetMissingValue;
                m_edges[e].second = sizetMissingValue;
            }
        }
    }

    m_nodesRTreeRequiresUpdate = true;
    m_edgesRTreeRequiresUpdate = true;

    Administrate(AdministrationOption::AdministrateMeshEdges);
}

void Mesh2D::DeleteHangingEdges()
{
    const auto hangingEdges = GetHangingEdges();
    for (const auto& hangingEdge : hangingEdges)
    {
        DeleteEdge(hangingEdge);
    }
}

std::vector<size_t> Mesh2D::PointFaceIndices(const std::vector<Point>& points)
{
    std::vector<size_t> result;
    result.resize(points.size(), sizetMissingValue);
    for (auto i = 0; i < points.size(); ++i)
    {
        const auto edgeIndex = FindEdgeCloseToAPoint(points[i]);

        if (edgeIndex == sizetMissingValue)
        {
            result[i] = sizetMissingValue;
            continue;
        }

        for (auto e = 0; e < m_edgesNumFaces[edgeIndex]; ++e)
        {
            const auto faceIndex = m_edgesFaces[edgeIndex][e];
            ComputeFaceClosedPolygon(faceIndex, m_polygonNodesCache);
            const auto isPointInFace = IsPointInPolygonNodes(points[i], m_polygonNodesCache, m_projection);
            if (isPointInFace)
            {
                result[i] = faceIndex;
                break;
            }
        }
    }
    return result;
}

std::tuple<size_t, size_t> Mesh2D::IsSegmentCrossingABoundaryEdge(const Point& firstPoint,
                                                                  const Point& secondPoint) const
{
    double intersectionRatio = std::numeric_limits<double>::max();
    size_t intersectedFace = sizetMissingValue;
    size_t intersectedEdge = sizetMissingValue;
    for (auto e = 0; e < GetNumEdges(); ++e)
    {
        if (!IsEdgeOnBoundary(e))
        {
            continue;
        }

        Point intersectionPoint;
        double crossProduct;
        double ratioFirstSegment;
        double ratioSecondSegment;
        const auto areSegmentCrossing = AreSegmentsCrossing(firstPoint,
                                                            secondPoint,
                                                            m_nodes[m_edges[e].first],
                                                            m_nodes[m_edges[e].second],
                                                            false,
                                                            m_projection,
                                                            intersectionPoint,
                                                            crossProduct,
                                                            ratioFirstSegment,
                                                            ratioSecondSegment);

        if (areSegmentCrossing && ratioFirstSegment < intersectionRatio)
        {
            intersectionRatio = ratioFirstSegment;
            intersectedFace = m_edgesFaces[e][0];
            intersectedEdge = e;
        }
    }

    return {intersectedFace, intersectedEdge};
}
