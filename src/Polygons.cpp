//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2020.
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

#pragma once

#include <utility>
#include <vector>
#include <stdexcept>
#include <MeshKernel/Mesh.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Operations.cpp>
#include <MeshKernel/TriangulationWrapper.hpp>
#include <MeshKernel/Exceptions.hpp>

namespace meshkernel
{
    Polygons::Polygons()
    {
    }

    Polygons::Polygons(const std::vector<Point>& polygon, Projections projection) : m_nodes(polygon), m_projection(projection)
    {
        // find the polygons in the current list of points
        m_indices = FindIndexes(polygon, 0, polygon.size(), doubleMissingValue);
    }

    std::vector<Point> Polygons::MeshBoundaryToPolygon(Mesh& mesh) const
    {

        // Find faces
        mesh.Administrate(Mesh::AdministrationOptions::AdministrateMeshEdgesAndFaces);
        std::vector<bool> isVisited(mesh.GetNumEdges(), false);
        std::vector<Point> meshBoundaryPolygon;
        meshBoundaryPolygon.reserve(mesh.GetNumNodes());

        for (int e = 0; e < mesh.GetNumEdges(); e++)
        {
            if (isVisited[e] || !mesh.IsEdgeOnBoundary(e))
            {
                continue;
            }

            const auto firstNodeIndex = mesh.m_edges[e].first;
            const auto secondNodeIndex = mesh.m_edges[e].second;
            const auto firstNode = mesh.m_nodes[firstNodeIndex];
            const auto secondNode = mesh.m_nodes[secondNodeIndex];

            const auto firstNodeInPolygon = IsPointInPolygonNodes(mesh.m_nodes[firstNodeIndex], m_nodes, 0, GetNumNodes() - 1, mesh.m_projection);
            const auto secondNodeInPolygon = IsPointInPolygonNodes(mesh.m_nodes[secondNodeIndex], m_nodes, 0, GetNumNodes() - 1, mesh.m_projection);

            if (!firstNodeInPolygon && !secondNodeInPolygon)
            {
                continue;
            }

            //Start a new polyline
            if (!meshBoundaryPolygon.empty())
            {
                meshBoundaryPolygon.emplace_back(doubleMissingValue,doubleMissingValue);
            }

            // Put the current edge on the mesh boundary, mark it as visited
            const auto startPolygonEdges = meshBoundaryPolygon.size();
            meshBoundaryPolygon.emplace_back(firstNode);
            meshBoundaryPolygon.emplace_back(secondNode);
            isVisited[e] = true;

            // walk the current mesh boundary
            auto currentNode = secondNodeIndex;
            WalkBoundaryFromNode(mesh, isVisited, currentNode, meshBoundaryPolygon);

            const auto numNodesFirstTail = meshBoundaryPolygon.size();

            // if the boundary polygon is not closed
            if (currentNode != firstNodeIndex)
            {
                //Now grow a polyline starting at the other side of the original link L, i.e., the second tail
                currentNode = firstNodeIndex;
                WalkBoundaryFromNode(mesh, isVisited, currentNode, meshBoundaryPolygon);
            }

            // There is a nonempty second tail, so reverse the first tail, so that they connect.
            if (meshBoundaryPolygon.size() > numNodesFirstTail)
            {
                const auto start = startPolygonEdges + static_cast<size_t>(std::ceil((numNodesFirstTail - startPolygonEdges + 1) * 0.5));
                for (auto n = start; n < numNodesFirstTail; n++)
                {
                    const auto backupPoint = meshBoundaryPolygon[n];
                    const int replaceIndex = numNodesFirstTail - n + static_cast<size_t>(firstNodeIndex);
                    meshBoundaryPolygon[n] = meshBoundaryPolygon[replaceIndex];
                    meshBoundaryPolygon[replaceIndex] = backupPoint;
                }
            }

            //Start a new polyline
            meshBoundaryPolygon.emplace_back(doubleMissingValue, doubleMissingValue);
        }
        return meshBoundaryPolygon;
    }

    void Polygons::WalkBoundaryFromNode(const Mesh& mesh,
                                        std::vector<bool>& isVisited,
                                        int& currentNode,
                                        std::vector<Point>& meshBoundaryPolygon) const
    {
        int e = 0;
        bool currentNodeInPolygon = false;
        while (e < mesh.m_nodesNumEdges[currentNode])
        {
            if (!currentNodeInPolygon)
            {
                currentNodeInPolygon = IsPointInPolygonNodes(mesh.m_nodes[currentNode], m_nodes, 0, GetNumNodes() - 1, m_projection);
            }

            if (!currentNodeInPolygon)
            {
                break;
            }

            const auto currentEdge = mesh.m_nodesEdges[currentNode][e];
            if (isVisited[currentEdge] || !mesh.IsEdgeOnBoundary(currentEdge))
            {
                e++;
                continue;
            }

            const auto firstNode = mesh.m_edges[currentEdge].first;
            const auto secondNode = mesh.m_edges[currentEdge].second;

            currentNode = secondNode == currentNode ? firstNode : secondNode;
            e = 0;
            currentNodeInPolygon = false;

            meshBoundaryPolygon.emplace_back(mesh.m_nodes[currentNode]);
            isVisited[currentEdge] = true;
        }
    }

    std::vector<std::vector<Point>> Polygons::ComputePointsInPolygons() const
    {

        std::vector<std::vector<Point>> generatedPoints;
        generatedPoints.reserve(m_indices.size());
        std::vector<Point> localPolygon(GetNumNodes());

        for (int i = 0; i < m_indices.size(); ++i)
        {
            localPolygon.clear();
            for (auto j = m_indices[i][0]; j <= m_indices[i][1]; ++j)
            {
                localPolygon.emplace_back(m_nodes[j]);
            }

            // not a closed polygon
            const auto numLocalPoints = localPolygon.size();
            if (localPolygon[numLocalPoints - 1] != localPolygon[0])
            {
                continue;
            }

            double localPolygonArea = 0.0;
            Point centerOfMass;
            bool isCounterClockWise;
            FaceAreaAndCenterOfMass(localPolygon, numLocalPoints - 1, m_projection, localPolygonArea, centerOfMass, isCounterClockWise);

            const auto perimeter = PerimeterClosedPolygon(localPolygon);

            double maximumEdgeLength;
            MaximumEdgeLength(localPolygon, numLocalPoints, maximumEdgeLength);

            // average triangle size
            const double averageEdgeLength = perimeter / static_cast<double>(numLocalPoints - 1);
            const double averageTriangleArea = 0.25 * squareRootOfThree * averageEdgeLength * averageEdgeLength;

            // estimated number of triangles
            const int SafetySize = 11;
            const auto numberOfTriangles = int(SafetySize * localPolygonArea / averageTriangleArea);
            if (numberOfTriangles <= 0)
            {
                throw AlgorithmError("Polygons::ComputePointsInPolygons: The number of triangles is <= 0.");
            }

            TriangulationWrapper triangulationWrapper;

            const auto numPolygonNodes = static_cast<int>(localPolygon.size() - 1); // open polygon

            triangulationWrapper.Compute(localPolygon,
                                         numPolygonNodes,
                                         TriangulationWrapper::TriangulationOptions::GeneratePoints,
                                         averageTriangleArea,
                                         numberOfTriangles);

            generatedPoints.emplace_back(triangulationWrapper.m_nodes);
        }

        return generatedPoints;
    }

    void Polygons::RefinePolygonPart(int startIndex,
                                     int endIndex,
                                     double refinementDistance,
                                     std::vector<Point>& refinedPolygon)
    {
        if (m_indices.empty())
        {
            throw std::invalid_argument("Polygons::RefinePolygonPart: No nodes in polygon.");
        }

        if (startIndex == 0 && endIndex == 0)
        {
            startIndex = m_indices[0][0];
            endIndex = m_indices[0][1];
        }

        if (endIndex <= startIndex)
        {
            throw std::invalid_argument("Polygons::RefinePolygonPart: The end index is smaller than the start index.");
        }

        bool areIndicesValid = false;
        int polygonIndex;
        for (int i = 0; i < m_indices.size(); ++i)
        {
            if (startIndex >= m_indices[i][0] && endIndex <= m_indices[i][1])
            {
                areIndicesValid = true;
                polygonIndex = i;
                break;
            }
        }

        if (!areIndicesValid)
        {
            throw std::invalid_argument("Polygons::RefinePolygonPart: The indices are not valid.");
        }

        const auto edgeLengths = PolygonEdgeLengths(m_nodes);
        std::vector<double> nodeLengthCoordinate(edgeLengths.size());
        nodeLengthCoordinate[0] = 0.0;
        for (int i = 1; i < edgeLengths.size(); ++i)
        {
            nodeLengthCoordinate[i] = nodeLengthCoordinate[i - 1] + edgeLengths[i - 1];
        }

        auto numNodesRefinedPart = int(std::ceil((nodeLengthCoordinate[endIndex] - nodeLengthCoordinate[startIndex]) / refinementDistance) + (endIndex - startIndex));
        int numNodesNotRefinedPart = startIndex - m_indices[polygonIndex][0] + m_indices[polygonIndex][1] - endIndex;
        int totalNumNodes = numNodesRefinedPart + numNodesNotRefinedPart;
        refinedPolygon.resize(totalNumNodes);

        // before refinement
        int refinedNodeIndex = 0;
        for (int i = m_indices[polygonIndex][0]; i <= startIndex; ++i)
        {
            refinedPolygon[refinedNodeIndex] = m_nodes[i];
            refinedNodeIndex++;
        }

        // refined part
        int nodeIndex = startIndex;
        int nextNodeIndex = nodeIndex + 1;
        Point p0 = m_nodes[nodeIndex];
        Point p1 = m_nodes[nextNodeIndex];
        double pointLengthCoordinate = nodeLengthCoordinate[startIndex];
        bool snappedToLastPoint = false;
        while (nodeIndex < endIndex)
        {
            // initial point already accounted for
            pointLengthCoordinate += refinementDistance;
            if (pointLengthCoordinate > nodeLengthCoordinate[nextNodeIndex])
            {
                // if not snapped to the original last polygon point, snap it
                if (!snappedToLastPoint)
                {
                    refinedPolygon[refinedNodeIndex] = m_nodes[nextNodeIndex];
                    refinedNodeIndex++;
                }

                // find the next point
                bool nextNodeFound = false;
                for (int i = nextNodeIndex + 1; i <= endIndex; ++i)
                {
                    if (nodeLengthCoordinate[i] > pointLengthCoordinate)
                    {
                        nextNodeFound = true;
                        nodeIndex = i - 1;
                        nextNodeIndex = i;
                        break;
                    }
                }
                if (nextNodeIndex > endIndex || !nextNodeFound)
                {
                    break;
                }

                p0 = m_nodes[nodeIndex];
                p1 = m_nodes[nextNodeIndex];
                pointLengthCoordinate = nodeLengthCoordinate[nodeIndex] + refinementDistance;
                snappedToLastPoint = false;
            }
            double distanceFromLastNode = pointLengthCoordinate - nodeLengthCoordinate[nodeIndex];
            double factor = distanceFromLastNode / edgeLengths[nodeIndex];
            Point p;
            if (std::abs(factor - 1.0) <= std::numeric_limits<double>::epsilon())
            {
                snappedToLastPoint = true;
                p = p1;
            }
            else
            {
                p = p0 + (p1 - p0) * distanceFromLastNode / edgeLengths[nodeIndex];
            }
            refinedPolygon[refinedNodeIndex] = p;
            refinedNodeIndex++;
        }

        // after refinement
        for (int i = endIndex + 1; i <= m_indices[polygonIndex][1]; ++i)
        {
            refinedPolygon[refinedNodeIndex] = m_nodes[i];
            refinedNodeIndex++;
        }
        refinedPolygon.resize(refinedNodeIndex);
    }

    void Polygons::MaximumEdgeLength(const std::vector<Point>& localPolygon, size_t numPoints, double& maximumEdgeLength) const
    {

        if (localPolygon[0].x != localPolygon[numPoints - 1].x)
        {
            throw std::invalid_argument("Polygons::MaximumEdgeLength: The first and last point of the polygon is not the same.");
        }

        maximumEdgeLength = std::numeric_limits<double>::lowest();
        for (int p = 0; p < numPoints - 1; ++p)
        {
            double edgeLength = ComputeDistance(m_nodes[p], m_nodes[p + 1], m_projection);
            maximumEdgeLength = std::max(maximumEdgeLength, edgeLength);
        }
    }

    void Polygons::OffsetCopy(double distance, bool innerAndOuter, Polygons& newPolygon)
    {
        int sizenewPolygon = GetNumNodes();
        if (innerAndOuter)
        {
            sizenewPolygon += GetNumNodes() + 1;
        }

        std::vector<Point> normalVectors(sizenewPolygon);
        double dxNormalPreviousEdge = 0.0;
        double dyNormalPreviuosEdge = 0.0;
        double dxNormal = 0.0;
        double dyNormal = 0.0;
        for (int n = 0; n < GetNumNodes(); n++)
        {
            if (n < GetNumNodes() - 1)
            {
                auto dx = GetDx(m_nodes[n], m_nodes[n + 1], m_projection);
                auto dy = GetDy(m_nodes[n], m_nodes[n + 1], m_projection);
                const auto nodeDistance = std::sqrt(dx * dx + dy * dy);
                dxNormal = -dy / nodeDistance;
                dyNormal = dx / nodeDistance;
            }
            else
            {
                dxNormal = dxNormalPreviousEdge;
                dyNormal = dyNormalPreviuosEdge;
            }

            if (n == 0)
            {
                dxNormalPreviousEdge = dxNormal;
                dyNormalPreviuosEdge = dyNormal;
            }

            double factor = 1.0 / (1.0 + dxNormalPreviousEdge * dxNormal + dyNormalPreviuosEdge * dyNormal);
            normalVectors[n].x = factor * (dxNormalPreviousEdge + dxNormal);
            normalVectors[n].y = factor * (dyNormalPreviuosEdge + dyNormal);

            dxNormalPreviousEdge = dxNormal;
            dyNormalPreviuosEdge = dyNormal;
        }

        // negative sign introduced because normal vector pointing inward
        distance = -distance;
        if (m_projection == Projections::spherical)
        {
            distance = distance / (earth_radius * degrad_hp);
        }

        std::vector<Point> newPolygonPoints(sizenewPolygon, {doubleMissingValue, doubleMissingValue});
        for (auto i = 0; i < GetNumNodes(); i++)
        {
            auto dx = normalVectors[i].x * distance;
            auto dy = normalVectors[i].y * distance;
            if (m_projection == Projections::spherical)
            {
                dx = dx / std::cos((m_nodes[i].y + 0.5 * dy) * degrad_hp);
            }
            newPolygonPoints[i].x = m_nodes[i].x + dx;
            newPolygonPoints[i].y = m_nodes[i].y + dy;

            if (innerAndOuter)
            {
                newPolygonPoints[i + GetNumNodes() + 1].x = m_nodes[i].x - dx;
                newPolygonPoints[i + GetNumNodes() + 1].y = m_nodes[i].y - dy;
            }
        }

        // set the new polygon
        newPolygon = {newPolygonPoints, m_projection};
    }

    bool Polygons::IsPointInPolygon(const Point& point, int polygonIndex) const
    {
        // empty polygon means everything is included
        if (m_indices.empty())
        {
            return true;
        }

        if (polygonIndex >= m_indices.size())
        {
            throw std::invalid_argument("Polygons::IsPointInPolygon: Invalid polygon index.");
        }

        bool inPolygon = IsPointInPolygonNodes(point, m_nodes, m_indices[polygonIndex][0], m_indices[polygonIndex][1], m_projection);
        return inPolygon;
    }

    bool Polygons::IsPointInPolygons(const Point& point) const
    {
        // empty polygon means everything is included
        if (m_indices.empty())
        {
            return true;
        }

        bool inPolygon = false;

        for (const auto& indexes : m_indices)
        {
            // Calculate the bounding box
            double XMin = std::numeric_limits<double>::max();
            double XMax = std::numeric_limits<double>::lowest();
            double YMin = std::numeric_limits<double>::max();
            double YMax = std::numeric_limits<double>::lowest();

            for (int n = indexes[0]; n <= indexes[1]; n++)
            {
                XMin = std::min(XMin, m_nodes[n].x);
                XMax = std::max(XMax, m_nodes[n].x);
                YMin = std::min(YMin, m_nodes[n].y);
                YMax = std::max(YMax, m_nodes[n].y);
            }

            if ((point.x >= XMin && point.x <= XMax) && (point.y >= YMin && point.y <= YMax))
            {
                inPolygon = IsPointInPolygonNodes(point, m_nodes, indexes[0], indexes[1], m_projection);
            }

            if (inPolygon)
            {
                return true;
            }
        }

        return inPolygon;
    }

    bool Polygons::IsEmpty() const
    {
        return m_indices.empty();
    }

    double Polygons::PerimeterClosedPolygon(const std::vector<Point>& polygonNodes) const
    {
        if (polygonNodes.front() != polygonNodes.back())
        {
            throw std::invalid_argument("Polygons::PerimeterClosedPolygon: The first and last point of the polygon is not the same.");
        }

        const auto edgeLengths = PolygonEdgeLengths(polygonNodes);
        return std::accumulate(edgeLengths.begin(), edgeLengths.end(), 0.0);
    }

    std::vector<double> Polygons::PolygonEdgeLengths(const std::vector<Point>& polygonNodes) const
    {
        std::vector<double> edgeLengths;
        edgeLengths.reserve(polygonNodes.size());

        for (auto p = 0; p < polygonNodes.size(); ++p)
        {
            const auto firstNode = p;
            auto secondNode = p + 1;
            if (secondNode == polygonNodes.size())
            {
                secondNode = 0;
            }
            edgeLengths.emplace_back(ComputeDistance(polygonNodes[firstNode], polygonNodes[secondNode], m_projection));
        }
        return edgeLengths;
    }

}; // namespace meshkernel
