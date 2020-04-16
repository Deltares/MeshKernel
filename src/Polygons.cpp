#pragma once

#include <utility>
#include <vector>
#include <iostream>
#include "Mesh.hpp"
#include "Polygons.hpp"
#include "Constants.cpp"
#include "Operations.cpp"

namespace GridGeom
{
    Polygons::Polygons() : m_numNodes(0), m_numAllocatedNodes(0)
    {
        AllocateVector(m_numAllocatedNodes, m_nodes, m_allocationSize);
        m_numAllocatedNodes = m_nodes.size();
    }

    bool Polygons::Set(const std::vector<Point>& polygon, Projections projection)
    {
        m_projection = projection;
        // resize if necessary
        int startNewNodes = m_numNodes;
        if (startNewNodes != 0)
        {
            ResizeVector(m_numNodes + 1 + polygon.size(), m_nodes);
            m_nodes[startNewNodes] = { doubleMissingValue,doubleMissingValue };
            startNewNodes += 1;
        }
        else
        {
            ResizeVector(polygon.size(), m_nodes);
        }
        m_numAllocatedNodes = m_nodes.size();
        m_numNodes = m_nodes.size();

        for (int n = startNewNodes, nn = 0; nn < startNewNodes + polygon.size(); ++n, ++nn)
        {
            m_nodes[n] = polygon[nn];
        }
        return true;
    }

    /// copynetboundstopol
    bool Polygons::MeshBoundaryToPolygon(const Mesh& mesh,
        int counterClockWise,
        int setMeshState,
        std::vector<Point>& meshBoundaryPolygon,
        int& numNodesBoundaryPolygons)
    {
        std::vector<bool> isVisited(mesh.GetNumEdges());
        std::vector<int> boundaryPolygonStarts(mesh.GetNumEdges());
        int numBoundaryPolygons = 0;
        numNodesBoundaryPolygons = 0;

        meshBoundaryPolygon.resize(mesh.GetNumNodes() , { doubleMissingValue ,doubleMissingValue });
        int meshBoundaryPolygonSize = meshBoundaryPolygon.size();

        for (int e = 0; e < mesh.GetNumEdges(); e++)
        {
            if (isVisited[e] || mesh.m_edgesNumFaces[e] != 1)
            {
                continue;
            }

            const int first = mesh.m_edges[e].first;
            const int second = mesh.m_edges[e].second;
            const auto firstPoint = mesh.m_nodes[first];
            const auto secondPoint = mesh.m_nodes[second];

            bool inHullFirst = IsPointInPolygon(mesh.m_nodes[first], m_nodes, m_numNodes);
            bool inHullSecond = IsPointInPolygon(mesh.m_nodes[second], m_nodes, m_numNodes);

            if (!inHullFirst && !inHullSecond)
            {
                continue;
            }

            ResizeVector(numNodesBoundaryPolygons + 3, meshBoundaryPolygon);

            //Start a new polyline
            if (numNodesBoundaryPolygons > 0)
            {
                numNodesBoundaryPolygons++;
            }

            boundaryPolygonStarts[numBoundaryPolygons] = numNodesBoundaryPolygons;
            numBoundaryPolygons++;

            const int startPolygonEdges = numNodesBoundaryPolygons;
            const int nodeStart = first;

            meshBoundaryPolygon[numNodesBoundaryPolygons] = firstPoint;
            numNodesBoundaryPolygons++;
            meshBoundaryPolygon[numNodesBoundaryPolygons] = secondPoint;

            isVisited[e] = true;

            // walk the current mesh boundary
            int currentNode = second;
            WalkBoundary(mesh, isVisited, numNodesBoundaryPolygons, currentNode, meshBoundaryPolygonSize, meshBoundaryPolygon);

            const int numNodesFirstTail = numNodesBoundaryPolygons;

            // if the boundary polygon is not closed
            if (currentNode != nodeStart)
            {
                //Now grow a polyline starting at the other side of the original link L, i.e., the second tail
                currentNode = nodeStart;
                WalkBoundary(mesh, isVisited, numNodesBoundaryPolygons, currentNode, meshBoundaryPolygonSize, meshBoundaryPolygon);
            }

            // There is a nonempty second tail, so reverse the first tail, so that they connect.
            if (numNodesBoundaryPolygons > numNodesFirstTail)
            {
                const int start = startPolygonEdges + std::ceil((numNodesFirstTail - startPolygonEdges + 1) / 2.0);
                Point backupPoint;
                for (int n = start; n < numNodesFirstTail; n++)
                {
                    backupPoint = meshBoundaryPolygon[n];
                    const int replaceIndex = numNodesFirstTail - n + nodeStart;
                    meshBoundaryPolygon[n] = meshBoundaryPolygon[replaceIndex];
                    meshBoundaryPolygon[replaceIndex] = backupPoint;
                }
            }
        }

        boundaryPolygonStarts[numBoundaryPolygons] = numNodesBoundaryPolygons + 1;

        return true;
    }


    bool Polygons::WalkBoundary(const Mesh& mesh,
        std::vector<bool>& isVisited,
        int& nodeIndex,
        int& currentNode,
        int meshBoundaryPolygonSize,
        std::vector<Point>& meshBoundaryPolygon)
    {
        int ee = 0;
        while (ee < mesh.m_nodesNumEdges[currentNode])
        {
            bool inHull = IsPointInPolygon(mesh.m_nodes[currentNode], m_nodes, m_numNodes);

            if (!inHull)
            {
                break;
            }

            const int currentEdge = mesh.m_nodesEdges[currentNode][ee];

            if (isVisited[currentEdge] == 1 || mesh.m_edgesNumFaces[currentEdge] != 1)
            {
                ee++;
                continue;
            }

            const int firstNode = mesh.m_edges[currentEdge].first;
            const int secondNode = mesh.m_edges[currentEdge].second;

            if (secondNode == currentNode)
            {
                currentNode = firstNode;
                ee = 0;
            }
            else
            {
                currentNode = secondNode;
                ee = 0;
            }

            nodeIndex++;

            ResizeVector(nodeIndex + 1, meshBoundaryPolygon);

            meshBoundaryPolygon[nodeIndex] = mesh.m_nodes[currentNode];

            isVisited[currentEdge] = true;
        }
        return true;
    }

    /// triangulate..
    bool Polygons::CreatePointsInPolygons(std::vector<std::vector<GridGeom::Point>>& generatedPoints)
    {        
        std::vector<std::vector<int>> indexes(m_numNodes, std::vector<int>(2));
        int pos = FindIndexes(m_nodes, 0, m_numNodes, doubleMissingValue, indexes);
        indexes.resize(pos);

        generatedPoints.resize(pos);
        std::vector<Point> localPolygon(m_numNodes);
        std::vector<double> xLocalPolygon(m_numNodes);
        std::vector<double> yLocalPolygon(m_numNodes);
        std::vector<int> faceNodes;
        std::vector<int> edgeNodes;
        std::vector<int> faceEdges;
        std::vector<double> xPoint;
        std::vector<double> yPoint;
        const int safetySize = 11;
        for (int i = 0; i < indexes.size(); ++i)
        {
            int numLocalPoints = 0;
            for (int j = indexes[i][0]; j <= indexes[i][1]; ++j)
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

            double localPolygonArea = 0.0;
            Point centerOfMass;
            bool success = faceAreaAndCenterOfMass(localPolygon, numLocalPoints - 1, localPolygonArea, centerOfMass, m_projection);
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
                faceNodes.resize(numberOfTriangles * 3);
                edgeNodes.resize(numberOfTriangles * 2);
                faceEdges.resize(numberOfTriangles * 3);
                xPoint.resize(numberOfTriangles * 3);
                yPoint.resize(numberOfTriangles * 3);
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

        return true;
    }

    bool Polygons::RefinePart(int startIndex, int endIndex, double refinementDistance, std::vector<Point>& refinedPolygon)
    {
        std::vector<std::vector<int>> indexes(m_numNodes, std::vector<int>(2));
        int pos = FindIndexes(m_nodes, 0, m_numNodes, doubleMissingValue, indexes);
        indexes.resize(pos);

        if(startIndex==0 && endIndex==0)
        {
            startIndex = indexes[0][0];
            endIndex = indexes[0][1];
        }

        if (endIndex <= startIndex)
        {
            return false;
        }

        bool polygonFound = false;
        int polygonIndex;
        for (int i = 0; i < indexes.size(); ++i)
        {
            if (startIndex >= indexes[i][0] && endIndex <= indexes[i][1])
            {
                polygonFound = true;
                polygonIndex = i;
                break;
            }
        }

        if (!polygonFound)
        {
            return false;
        }

        std::vector<double> edgeLengths;
        EdgeLengths(m_nodes, edgeLengths);
        std::vector<double> nodeLengthCoordinate(edgeLengths.size());
        nodeLengthCoordinate[0] = 0.0;
        for (int i = 1; i < edgeLengths.size(); ++i)
        {
            nodeLengthCoordinate[i] = nodeLengthCoordinate[i - 1] + edgeLengths[i-1];
        }

        int numNodesRefinedPart = std::ceil((nodeLengthCoordinate[endIndex] - nodeLengthCoordinate[startIndex]) / refinementDistance) + (endIndex - startIndex);
        int numNodesNotRefinedPart = startIndex - indexes[polygonIndex][0] + indexes[polygonIndex][1] - endIndex;
        int totalNumNodes = numNodesRefinedPart + numNodesNotRefinedPart;
        refinedPolygon.resize(totalNumNodes);

        // before refinement
        int refinedNodeIndex = 0;
        for (int i = indexes[polygonIndex][0]; i <= startIndex; ++i)
        {
            refinedPolygon[refinedNodeIndex] = m_nodes[i];
            refinedNodeIndex++;
        }

        // refined part
        int nodeIndex = startIndex;
        int nextNodeIndex = nodeIndex + 1;
        Point p0 = m_nodes[nodeIndex];
        Point p1 = m_nodes[nodeIndex + 1];
        double pointLengthCoordinate = nodeLengthCoordinate[startIndex];
        while (nodeIndex < endIndex)
        {
            pointLengthCoordinate += refinementDistance;
            if (pointLengthCoordinate > nodeLengthCoordinate[nextNodeIndex])
            {
                // find next point 
                bool nextNodeFound = false;
                for (int i = nextNodeIndex + 1; i <= endIndex; ++i)
                {
                    if(nodeLengthCoordinate[i]>pointLengthCoordinate)
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
            }
            double distanceFromLastNode = pointLengthCoordinate - nodeLengthCoordinate[nodeIndex];
            Point p = p0 + (p1 - p0) *  distanceFromLastNode / edgeLengths[nodeIndex];
            refinedPolygon[refinedNodeIndex] = p;
            refinedNodeIndex++;
        }

        // after refinement
        for (int i = endIndex + 1; i <= indexes[polygonIndex][1]; ++i)
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
        EdgeLengths(localPolygon, edgeLengths);
        perimeter = std::accumulate(edgeLengths.begin(), edgeLengths.end(), 0.0);
        return true;
    }

    bool Polygons::EdgeLengths(const std::vector<Point>& localPolygon, std::vector<double>& edgeLengths)
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
        std::vector<std::vector<int>> indexes(m_numNodes, std::vector<int>(2));
        int pos = FindIndexes(m_nodes, 0, m_numNodes, doubleMissingValue, indexes);
        indexes.resize(pos);

        
        int sizenewPolygon = m_numNodes;
        if (innerAndOuter) 
        {
            sizenewPolygon += m_numNodes + 1;
        }
        
        std::vector<Point> normalVectors(sizenewPolygon);
        double dxNormalPreviusEdge;
        double dyNormalPreviusEdge;
        double dxNormalNodeIndex;
        double dyNormalNodeIndex;
        double dxNormalPreviusEdgeNodeIndex;
        double dyNormalPreviusEdgeNodeIndex;
        Point normalVectorNodeIndex;
        for (int n = 0; n < m_numNodes; n++)
        {
            double dxNormal;
            double dyNormal;
            if (n < m_numNodes - 1)
            {
                auto dx = GetDx(m_nodes[n], m_nodes[n + 1], m_projection);
                auto dy = GetDy(m_nodes[n], m_nodes[n + 1], m_projection);
                auto distance = std::sqrt(dx*dx + dy*dy);
                dxNormal = -dy / distance;
                dyNormal = dx / distance;
            }
            else
            {
                dxNormal = dxNormalPreviusEdge;
                dyNormal = dyNormalPreviusEdge;
            }

            if (n == 0)
            {
                dxNormalPreviusEdge = dxNormal;
                dyNormalPreviusEdge = dyNormal;
            }

            double factor = 1.0 / (1.0 + dxNormalPreviusEdge*dxNormal + dyNormalPreviusEdge*dyNormal);
            normalVectors[n].x = factor *(dxNormalPreviusEdge + dxNormal);
            normalVectors[n].y = factor *(dyNormalPreviusEdge + dyNormal);

            if (n == nodeIndex) 
            {
                dxNormalNodeIndex = dxNormal;
                dyNormalNodeIndex = dyNormal;
                dxNormalPreviusEdgeNodeIndex = dxNormalPreviusEdge;
                dyNormalPreviusEdgeNodeIndex = dyNormalPreviusEdge;
                normalVectorNodeIndex = normalVectors[n];
            }

            dxNormalPreviusEdge = dxNormal;
            dyNormalPreviusEdge = dyNormal;
        }

        // negative sign introduced because normal vector pointing inward
        distance = -distance;
        if (m_projection == Projections::spherical) 
        {
            distance = distance / (earth_radius * degrad_hp);
        }

        std::vector<Point> newPolygonPoints(sizenewPolygon, {doubleMissingValue, doubleMissingValue});
        for (int i = 0; i < m_nodes.size(); i++)
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
                newPolygonPoints[i + m_numNodes + 1].x = m_nodes[i].x - dx;
                newPolygonPoints[i + m_numNodes + 1].y = m_nodes[i].y - dy;
            }
        }

        // set the new polygon
        newPolygon.Set(newPolygonPoints, m_projection);

        return true;
    }



};