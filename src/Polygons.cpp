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
#include "Mesh.hpp"
#include "Polygons.hpp"
#include "Constants.cpp"
#include "Operations.cpp"

namespace MeshKernel
{
    Polygons::Polygons() : m_numNodes(0), m_numAllocatedNodes(0)
    {
        ResizeVectorIfNeededWithMinimumSize(m_numAllocatedNodes, m_nodes, m_allocationSize);
        m_numAllocatedNodes = m_nodes.size();
    }

    bool Polygons::Set(const std::vector<Point>& polygon, Projections projection)
    {
        // always set the projection first
        m_projection = projection;
        
        // in no points to add, return
        if (polygon.empty())
        {
            return true;
        }

        // find the polygons in the current list of points 
        std::vector<std::vector<size_t>> indexes(polygon.size(), std::vector<size_t>(2,0));
        int pos = FindIndexes(polygon, 0, polygon.size(), doubleMissingValue, indexes);
        indexes.resize(pos);

        // resize if necessary
        int numNodes = GetNumNodes();
        int currentNodePosition = numNodes;
        if (currentNodePosition != 0)
        {
            ResizeVectorIfNeeded(numNodes + 1 + polygon.size(), m_nodes);
            m_nodes[currentNodePosition] = { doubleMissingValue,doubleMissingValue };
        }
        else
        {
            ResizeVectorIfNeeded(polygon.size(), m_nodes);
        }

        m_numAllocatedNodes = m_nodes.size();
        m_numNodes = m_nodes.size();

        int numPolygons = m_indexses.size();
        ResizeVectorIfNeeded(numPolygons + indexes.size(), m_indexses, std::vector<int>(2, 0));
        for (int p = 0; p < indexes.size(); p++)
        {
            int indexInIndexes = 0;
            if (currentNodePosition > 0) 
            {
                m_nodes[currentNodePosition] = { doubleMissingValue,doubleMissingValue };
                currentNodePosition++;
            }
            for (int n = indexes[p][0]; n <= indexes[p][1]; n++)
            {
                m_nodes[currentNodePosition] = polygon[indexes[p][0] + indexInIndexes];
                indexInIndexes++;
                currentNodePosition++;
            }
            m_indexses[numPolygons + p][0] = indexes[p][0] + numNodes;
            m_indexses[numPolygons + p][1] = indexes[p][1] + numNodes;
        }

        return true;
    }

    /// copynetboundstopol
    bool Polygons::MeshBoundaryToPolygon(Mesh& mesh,
        int counterClockWise,
        std::vector<Point>& meshBoundaryPolygon,
        int& numNodesBoundaryPolygons)
    {
        // Find faces
        mesh.Administrate(Mesh::AdministrationOptions::AdministrateMeshEdgesAndFaces);

        std::vector<bool> isVisited(mesh.GetNumEdges(), false);
        numNodesBoundaryPolygons = 0;

        meshBoundaryPolygon.resize(mesh.GetNumNodes() , { doubleMissingValue ,doubleMissingValue });
        int meshBoundaryPolygonSize = meshBoundaryPolygon.size();

        for (int e = 0; e < mesh.GetNumEdges(); e++)
        {
            if (isVisited[e] || mesh.m_edgesNumFaces[e] != 1)
            {
                continue;
            }

            const auto firstNodeIndex = mesh.m_edges[e].first;
            const auto secondNodeIndex = mesh.m_edges[e].second;
            const auto firstNode = mesh.m_nodes[firstNodeIndex];
            const auto secondNode = mesh.m_nodes[secondNodeIndex];

            bool firstNodeInPolygon = IsPointInPolygonNodes(mesh.m_nodes[firstNodeIndex], m_nodes, 0,GetNumNodes()-1);
            bool secondNodeInPolygon = IsPointInPolygonNodes(mesh.m_nodes[secondNodeIndex], m_nodes, 0, GetNumNodes()-1);

            if (!firstNodeInPolygon && !secondNodeInPolygon)
            {
                continue;
            }

            ResizeVectorIfNeeded(numNodesBoundaryPolygons + 3, meshBoundaryPolygon,{doubleMissingValue, doubleMissingValue});

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
                const int start = startPolygonEdges + std::ceil((numNodesFirstTail - startPolygonEdges + 1) / 2.0);
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
       

        return true;
    }


    bool Polygons::WalkBoundaryFromNode(const Mesh& mesh,
        std::vector<bool>& isVisited,
        int& nodeIndex,
        int& currentNode,
        std::vector<Point>& meshBoundaryPolygon)
    {
        int e = 0;
        bool currentNodeInPolygon = false;
        while (e < mesh.m_nodesNumEdges[currentNode])
        {
            if (!currentNodeInPolygon)
            {
                currentNodeInPolygon = IsPointInPolygonNodes(mesh.m_nodes[currentNode], m_nodes, 0, GetNumNodes()-1);
            }

            if (!currentNodeInPolygon)
            {
                break;
            }

            const auto currentEdge = mesh.m_nodesEdges[currentNode][e];
            if (isVisited[currentEdge] || mesh.m_edgesNumFaces[currentEdge] != 1)
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
            ResizeVectorIfNeeded(nodeIndex + 1, meshBoundaryPolygon,{doubleMissingValue,doubleMissingValue});

            meshBoundaryPolygon[nodeIndex] = mesh.m_nodes[currentNode];

            isVisited[currentEdge] = true;
        }
        return true;
    }

    /// triangulate..
    bool Polygons::CreatePointsInPolygons(std::vector<std::vector<Point>>& generatedPoints)
    {        
        generatedPoints.resize(m_indexses.size());
        std::vector<Point> localPolygon(GetNumNodes());
        std::vector<double> xLocalPolygon(GetNumNodes());
        std::vector<double> yLocalPolygon(GetNumNodes());
        std::vector<int> faceNodes;
        std::vector<int> edgeNodes;
        std::vector<int> faceEdges;
        std::vector<double> xPoint;
        std::vector<double> yPoint;
        const int safetySize = 11;
        bool isOnePolygonClosed = false;
        for (int i = 0; i < m_indexses.size(); ++i)
        {
            int numLocalPoints = 0;
            for (int j = m_indexses[i][0]; j <= m_indexses[i][1]; ++j)
            {
                localPolygon[numLocalPoints] = m_nodes[j];
                xLocalPolygon[numLocalPoints] = m_nodes[j].x;
                yLocalPolygon[numLocalPoints] = m_nodes[j].y;
                numLocalPoints++;
            }

            // not a closed polygon
            if(localPolygon[numLocalPoints - 1] != localPolygon[0])
            {
                continue;
            }

            isOnePolygonClosed = true;
            double localPolygonArea = 0.0;
            Point centerOfMass;
            bool success = FaceAreaAndCenterOfMass(localPolygon, numLocalPoints - 1, m_projection, localPolygonArea, centerOfMass);
            if(!success)
            {
                return false;
            }

            double perimeter;
            success = PerimeterClosedPolygon(localPolygon, numLocalPoints, perimeter);
            if (!success)
            {
                return false;
            }

            double maximumEdgeLength;
            success = MaximumEdgeLength(localPolygon, numLocalPoints, maximumEdgeLength);
            if (!success)
            {
                return false;
            }

            // average edge size 
            double averageEdgeLength = perimeter / (numLocalPoints-1);

            // average triangle size
            double averageTriangleArea = 0.25 * squareRootOfThree * averageEdgeLength * averageEdgeLength;

            int numberOfTriangles =  safetySize * localPolygonArea / averageTriangleArea;

            if(numberOfTriangles<=0)
            {
                return false;
            }

            int numtri = -1;
            int jatri = 2;
            int numedge = 0;
            int numPoints = 0;
            int numLocalPointsOpenPolygon = numLocalPoints - 1;
            // if the number of estimated triangles is not sufficent, tricall must be repeated
            while(numtri < 0)
            {
                numtri = numberOfTriangles;
                faceNodes.resize(int(numberOfTriangles) * 3);
                edgeNodes.resize(int(numberOfTriangles) * 2);
                faceEdges.resize(int(numberOfTriangles) * 3);
                xPoint.resize(int(numberOfTriangles) * 3);
                yPoint.resize(int(numberOfTriangles) * 3);
                Triangulation(&jatri,
                    &xLocalPolygon[0],
                    &yLocalPolygon[0],
                    &numLocalPointsOpenPolygon,
                    &faceNodes[0],
                    &numtri,
                    &edgeNodes[0],
                    &numedge,
                    &faceEdges[0],
                    &xPoint[0],
                    &yPoint[0],
                    &numPoints,
                    &averageTriangleArea);
                if (numberOfTriangles)
                {
                    numberOfTriangles = -numtri;
                }
            }
            generatedPoints[i].resize(numPoints);
            for (int j = 0; j < numPoints; ++j)
            {
                generatedPoints[i][j] = { xPoint[j], yPoint[j] };
            }
            //TODO: CHECK POINTS ARE INSIDE THE POLYGON?
         }

        return isOnePolygonClosed ? true : false;
    }

    bool Polygons::RefinePolygonPart(int startIndex,
                                     int endIndex,
                                     double refinementDistance,
                                     std::vector<Point>& refinedPolygon)
    {
        if (m_indexses.empty())
        {
            return false;
        }

        if (startIndex == 0 && endIndex == 0)
        {
            startIndex = m_indexses[0][0];
            endIndex = m_indexses[0][1];
        }


        if (endIndex <= startIndex)
        {
            return false;
        }

        bool areIndexesValid = false;
        int polygonIndex;
        for (int i = 0; i < m_indexses.size(); ++i)
        {
            if (startIndex >= m_indexses[i][0] && endIndex <= m_indexses[i][1])
            {
                areIndexesValid = true;
                polygonIndex = i;
                break;
            }
        }

        if (!areIndexesValid)
        {
            return false;
        }

        std::vector<double> edgeLengths;
        PolygonEdgeLengths(m_nodes, edgeLengths);
        std::vector<double> nodeLengthCoordinate(edgeLengths.size());
        nodeLengthCoordinate[0] = 0.0;
        for (int i = 1; i < edgeLengths.size(); ++i)
        {
            nodeLengthCoordinate[i] = nodeLengthCoordinate[i-1] + edgeLengths[i-1];
        }

        int numNodesRefinedPart = std::ceil((nodeLengthCoordinate[endIndex] - nodeLengthCoordinate[startIndex]) / refinementDistance) + (int(endIndex) - int(startIndex));
        int numNodesNotRefinedPart = startIndex - m_indexses[polygonIndex][0] + m_indexses[polygonIndex][1] - endIndex;
        int totalNumNodes = numNodesRefinedPart + numNodesNotRefinedPart;
        refinedPolygon.resize(totalNumNodes);

        // before refinement
        int refinedNodeIndex = 0;
        for (int i = m_indexses[polygonIndex][0]; i <= startIndex; ++i)
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
        for (int i = endIndex + 1; i <= m_indexses[polygonIndex][1]; ++i)
        {
            refinedPolygon[refinedNodeIndex] = m_nodes[i];
            refinedNodeIndex++;
        }
        refinedPolygon.resize(refinedNodeIndex);

        return true;
    }

    bool Polygons::PerimeterClosedPolygon(const std::vector<Point>& localPolygon, const int numPoints, double& perimeter)
    {

        if(numPoints < 0 || localPolygon[0] != localPolygon[numPoints - 1])
        {
            return false;
        }

        perimeter = 0.0;
        std::vector<double> edgeLengths;
        PolygonEdgeLengths(localPolygon, edgeLengths);
        perimeter = std::accumulate(edgeLengths.begin(), edgeLengths.end(), 0.0);
        return true;
    }

    bool Polygons::PolygonEdgeLengths(const std::vector<Point>& localPolygon, std::vector<double>& edgeLengths)
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
            double edgeLength = Distance(localPolygon[firstNode], localPolygon[secondNode], m_projection);
            edgeLengths[p] = edgeLength;
        }

        return true;
    }

    bool Polygons::MaximumEdgeLength(const std::vector<Point>& localPolygon, const int numPoints, double& maximumEdgeLength)
    {

        if (numPoints < 0 || localPolygon[0].x != localPolygon[numPoints - 1].x)
        {
            return false;
        }

        maximumEdgeLength = std::numeric_limits<double>::min();
        for (int p = 0; p < numPoints - 1; ++p)
        {
            double edgeLength = Distance(m_nodes[p], m_nodes[p + 1], m_projection);
            maximumEdgeLength = std::max(maximumEdgeLength, edgeLength);
        }

        return true;
    }

    //copypol: look how the layer thickness is determined when the input distance is not given, but th coordinate of another point
    bool Polygons::OffsetCopy(int nodeIndex, double distance, bool innerAndOuter, Polygons& newPolygon)
    {        
        int sizenewPolygon = GetNumNodes();
        if (innerAndOuter) 
        {
            sizenewPolygon += GetNumNodes() + 1;
        }
        
        std::vector<Point> normalVectors(sizenewPolygon);
        double dxNormalPreviousEdge;
        double dyNormalPreviuosEdge;
        for (int n = 0; n < GetNumNodes(); n++)
        {
            double dxNormal;
            double dyNormal;
            if (n < GetNumNodes() - 1)
            {
                auto dx = GetDx(m_nodes[n], m_nodes[n + 1], m_projection);
                auto dy = GetDy(m_nodes[n], m_nodes[n + 1], m_projection);
                auto distance = std::sqrt(dx*dx + dy*dy);
                dxNormal = -dy / distance;
                dyNormal = dx / distance;
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

            double factor = 1.0 / (1.0 + dxNormalPreviousEdge*dxNormal + dyNormalPreviuosEdge*dyNormal);
            normalVectors[n].x = factor *(dxNormalPreviousEdge + dxNormal);
            normalVectors[n].y = factor *(dyNormalPreviuosEdge + dyNormal);

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
            if (m_projection== Projections::spherical)
            {
                dx = dx / std::cos((m_nodes[i].y + 0.5 *dy)*degrad_hp);
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
        newPolygon.Set(newPolygonPoints, m_projection);

        return true;
    }

    bool Polygons::IsPointInPolygon(const Point& point, int polygonIndex) const
    {
        if (polygonIndex >= m_indexses.size())
        {
            return true;
        }

        bool inPolygon = IsPointInPolygonNodes(point, m_nodes, m_indexses[polygonIndex][0], m_indexses[polygonIndex][1]);

        return inPolygon;
    }

    bool Polygons::IsPointInPolygons(const Point& point) const
    {
        
        // empty polygon, always included
        if (m_indexses.empty())
        {
            return true;
        }

        bool inPolygon = false;
        for (int p = 0; p < m_indexses.size(); p++)
        {
            // Calculate the bounding box
            double XMin = std::numeric_limits<double>::max();
            double XMax = std::numeric_limits<double>::min();
            double YMin = std::numeric_limits<double>::max();
            double YMax = std::numeric_limits<double>::min();

            for (int n = m_indexses[p][0]; n <= m_indexses[p][1]; n++)
            {
                XMin = std::min(XMin, m_nodes[n].x);
                XMax = std::max(XMax, m_nodes[n].x);
                YMin = std::min(YMin, m_nodes[n].y);
                YMax = std::max(YMax, m_nodes[n].y);
            }

            
            if ((point.x >= XMin && point.x <= XMax) && (point.y >= YMin && point.y <= YMax))
            {
                inPolygon = IsPointInPolygonNodes(point, m_nodes, m_indexses[p][0], m_indexses[p][1]);
            }

            if (inPolygon) 
            {
                return true;
            }  
        }

        return inPolygon;
    }




};