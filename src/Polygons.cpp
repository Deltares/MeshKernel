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

    Polygons::Polygons(const std::vector<Point>& polygon, Projections projection) : m_projection(projection)
    {
        ResizeVectorIfNeededWithMinimumSize(m_numAllocatedNodes, m_nodes, m_allocationSize);

        // find the polygons in the current list of points
        const auto indexes = FindIndexes(polygon, 0, polygon.size(), doubleMissingValue);

        // resize if necessary
        int numNodes = GetNumNodes();
        int currentNodePosition = numNodes;
        if (currentNodePosition != 0)
        {
            ResizeVectorIfNeeded(numNodes + 1 + int(polygon.size()), m_nodes);
            m_nodes[currentNodePosition] = {doubleMissingValue, doubleMissingValue};
        }
        else
        {
            ResizeVectorIfNeeded(int(polygon.size()), m_nodes);
        }

        m_numAllocatedNodes = int(m_nodes.size());
        m_numNodes = int(m_nodes.size());

        auto numPolygons = int(m_indices.size());
        ResizeVectorIfNeeded(numPolygons + int(indexes.size()), m_indices, std::vector<int>(2, 0));
        for (int p = 0; p < indexes.size(); p++)
        {
            int indexInIndexes = 0;
            if (currentNodePosition > 0)
            {
                m_nodes[currentNodePosition] = {doubleMissingValue, doubleMissingValue};
                currentNodePosition++;
            }

            for (auto n = indexes[p][0]; n <= indexes[p][1]; n++)
            {
                m_nodes[currentNodePosition] = polygon[indexes[p][0] + indexInIndexes];
                indexInIndexes++;
                currentNodePosition++;
            }
            m_indices[numPolygons + p][0] = int(indexes[p][0]) + numNodes;
            m_indices[numPolygons + p][1] = int(indexes[p][1]) + numNodes;
        }
    }

    void Polygons::MeshBoundaryToPolygon(Mesh& mesh,
                                         std::vector<Point>& meshBoundaryPolygon,
                                         int& numNodesBoundaryPolygons) const
    {
        numNodesBoundaryPolygons = 0;

        // Find faces
        mesh.Administrate(Mesh::AdministrationOptions::AdministrateMeshEdgesAndFaces);
        std::vector<bool> isVisited(mesh.GetNumEdges(), false);
        meshBoundaryPolygon.resize(mesh.GetNumNodes(), {doubleMissingValue, doubleMissingValue});

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

            bool firstNodeInPolygon = IsPointInPolygonNodes(mesh.m_nodes[firstNodeIndex], m_nodes, 0, GetNumNodes() - 1, mesh.m_projection);
            bool secondNodeInPolygon = IsPointInPolygonNodes(mesh.m_nodes[secondNodeIndex], m_nodes, 0, GetNumNodes() - 1, mesh.m_projection);

            if (!firstNodeInPolygon && !secondNodeInPolygon)
            {
                continue;
            }

            ResizeVectorIfNeeded(numNodesBoundaryPolygons + 3, meshBoundaryPolygon, {doubleMissingValue, doubleMissingValue});

            //Start a new polyline
            if (numNodesBoundaryPolygons > 0)
            {
                numNodesBoundaryPolygons++;
            }

            const int startPolygonEdges = numNodesBoundaryPolygons;

            meshBoundaryPolygon[numNodesBoundaryPolygons] = firstNode;
            numNodesBoundaryPolygons++;
            meshBoundaryPolygon[numNodesBoundaryPolygons] = secondNode;

            isVisited[e] = true;

            // walk the current mesh boundary
            auto currentNode = secondNodeIndex;
            WalkBoundaryFromNode(mesh, isVisited, numNodesBoundaryPolygons, currentNode, meshBoundaryPolygon);

            const auto numNodesFirstTail = numNodesBoundaryPolygons;

            // if the boundary polygon is not closed
            if (currentNode != firstNodeIndex)
            {
                //Now grow a polyline starting at the other side of the original link L, i.e., the second tail
                currentNode = firstNodeIndex;
                WalkBoundaryFromNode(mesh, isVisited, numNodesBoundaryPolygons, currentNode, meshBoundaryPolygon);
            }

            // There is a nonempty second tail, so reverse the first tail, so that they connect.
            if (numNodesBoundaryPolygons > numNodesFirstTail)
            {
                const int start = startPolygonEdges + int(std::ceil((numNodesFirstTail - startPolygonEdges + 1) / 2.0));
                Point backupPoint;
                for (int n = start; n < numNodesFirstTail; n++)
                {
                    backupPoint = meshBoundaryPolygon[n];
                    const int replaceIndex = numNodesFirstTail - n + firstNodeIndex;
                    meshBoundaryPolygon[n] = meshBoundaryPolygon[replaceIndex];
                    meshBoundaryPolygon[replaceIndex] = backupPoint;
                }
            }

            //Start a new polyline
            numNodesBoundaryPolygons++;
        }
    }

    void Polygons::WalkBoundaryFromNode(const Mesh& mesh,
                                        std::vector<bool>& isVisited,
                                        int& nodeIndex,
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

            nodeIndex++;
            ResizeVectorIfNeeded(nodeIndex + 1, meshBoundaryPolygon, {doubleMissingValue, doubleMissingValue});

            meshBoundaryPolygon[nodeIndex] = mesh.m_nodes[currentNode];

            isVisited[currentEdge] = true;
        }
    }

    void Polygons::CreatePointsInPolygons(std::vector<std::vector<Point>>& generatedPoints)
    {

        std::vector<Point> localPolygon(GetNumNodes());

        const int SafetySize = 11;
        bool isOnePolygonClosed = false;
        generatedPoints.resize(m_indices.size());
        for (int i = 0; i < m_indices.size(); ++i)
        {
            localPolygon.clear();
            for (int j = m_indices[i][0]; j <= m_indices[i][1]; ++j)
            {
                localPolygon.emplace_back(m_nodes[j]);
            }

            // not a closed polygon
            const auto numLocalPoints = localPolygon.size();
            if (localPolygon[numLocalPoints - 1] != localPolygon[0])
            {
                continue;
            }

            isOnePolygonClosed = true;
            double localPolygonArea = 0.0;
            Point centerOfMass;
            bool isCounterClockWise;
            FaceAreaAndCenterOfMass(localPolygon, numLocalPoints - 1, m_projection, localPolygonArea, centerOfMass, isCounterClockWise);

            double perimeter;
            PerimeterClosedPolygon(localPolygon, numLocalPoints, perimeter);

            double maximumEdgeLength;
            MaximumEdgeLength(localPolygon, numLocalPoints, maximumEdgeLength);

            // average triangle size
            const double averageEdgeLength = perimeter / static_cast<double>(numLocalPoints - 1);
            const double averageTriangleArea = 0.25 * squareRootOfThree * averageEdgeLength * averageEdgeLength;

            // estimated number of triangles
            const auto numberOfTriangles = int(SafetySize * localPolygonArea / averageTriangleArea);
            if (numberOfTriangles <= 0)
            {
                throw AlgorithmError("Polygons::CreatePointsInPolygons: The number of triangles is <= 0.");
            }

            TriangulationWrapper triangulationWrapper;

            const auto numPolygonNodes = static_cast<int>(localPolygon.size() - 1); // open polygon

            triangulationWrapper.Compute(localPolygon,
                                         numPolygonNodes,
                                         TriangulationWrapper::TriangulationOptions::GeneratePoints,
                                         averageTriangleArea,
                                         numberOfTriangles);

            generatedPoints[i] = std::move(triangulationWrapper.m_nodes);
        }
        if (!isOnePolygonClosed)
        {
            throw AlgorithmError("Polygons::CreatePointsInPolygons: There is no closed polygon.");
        }
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

        std::vector<double> edgeLengths;
        PolygonEdgeLengths(m_nodes, edgeLengths);
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

    void Polygons::PerimeterClosedPolygon(const std::vector<Point>& localPolygon, size_t numPoints, double& perimeter) const
    {
        if (localPolygon[0] != localPolygon[numPoints - 1])
        {
            throw std::invalid_argument("Polygons::PerimeterClosedPolygon: The first and last point of the polygon is not the same.");
        }

        perimeter = 0.0;
        std::vector<double> edgeLengths;
        PolygonEdgeLengths(localPolygon, edgeLengths);
        perimeter = std::accumulate(edgeLengths.begin(), edgeLengths.end(), 0.0);
    }

    void Polygons::PolygonEdgeLengths(const std::vector<Point>& localPolygon, std::vector<double>& edgeLengths) const
    {
        edgeLengths.resize(localPolygon.size());
        for (int p = 0; p < localPolygon.size(); ++p)
        {
            int firstNode = p;
            int secondNode = p + 1;
            if (secondNode == localPolygon.size())
            {
                secondNode = 0;
            }
            double edgeLength = ComputeDistance(localPolygon[firstNode], localPolygon[secondNode], m_projection);
            edgeLengths[p] = edgeLength;
        }
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

}; // namespace meshkernel
