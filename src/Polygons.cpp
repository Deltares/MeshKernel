#pragma once

#include <utility>
#include <vector>
#include "Mesh.hpp"
#include "Polygons.hpp"
#include "Constants.cpp"
#include "Operations.cpp"

namespace GridGeom
{
    Polygons::Polygons() : m_numNodes(0), m_numAllocatedNodes(0)
    {
        ResizeVectorIfNeededWithMinimumSize(m_numAllocatedNodes, m_nodes, m_allocationSize);
        m_numAllocatedNodes = m_nodes.size();
    }

    bool Polygons::Set(const std::vector<Point>& polygon, Projections projection)
    {
        m_projection = projection;
        // resize if necessary
        int numNodes = GetNumNodes();
        int startNewNodes = numNodes;
        if (numNodes != 0)
        {
            ResizeVectorIfNeeded(numNodes + 1 + polygon.size(), m_nodes);
            m_nodes[numNodes] = { doubleMissingValue,doubleMissingValue };
            startNewNodes += 1;
        }
        else
        {
            ResizeVectorIfNeeded(polygon.size(), m_nodes);
        }
        m_numAllocatedNodes = m_nodes.size();
        m_numNodes = m_nodes.size();

        for (int n = startNewNodes, nn = 0; n < startNewNodes + polygon.size(); ++n, ++nn)
        {
            m_nodes[n] = polygon[nn];
        }
        return true;
    }

    /// copynetboundstopol
    bool Polygons::MeshBoundaryToPolygon(const Mesh& mesh,
        int counterClockWise,
        std::vector<Point>& meshBoundaryPolygon,
        int& numNodesBoundaryPolygons)
    {
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

            bool firstNodeInPolygon = IsPointInPolygon(mesh.m_nodes[firstNodeIndex], m_nodes, GetNumNodes());
            bool secondNodeInPolygon = IsPointInPolygon(mesh.m_nodes[secondNodeIndex], m_nodes, GetNumNodes());

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
                currentNodeInPolygon = IsPointInPolygon(mesh.m_nodes[currentNode], m_nodes, GetNumNodes());
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
    bool Polygons::CreatePointsInPolygons(std::vector<std::vector<GridGeom::Point>>& generatedPoints)
    {        
        std::vector<std::vector<int>> indexes(GetNumNodes(), std::vector<int>(2));
        int pos = FindIndexes(m_nodes, 0, GetNumNodes(), doubleMissingValue, indexes);
        indexes.resize(pos);

        generatedPoints.resize(pos);
        std::vector<Point> localPolygon(GetNumNodes());
        std::vector<double> xLocalPolygon(GetNumNodes());
        std::vector<double> yLocalPolygon(GetNumNodes());
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
        std::vector<std::vector<int>> indexes(GetNumNodes(), std::vector<int>(2));
        int pos = FindIndexes(m_nodes, 0, GetNumNodes(), doubleMissingValue, indexes);
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
        std::vector<std::vector<int>> indexes(GetNumNodes(), std::vector<int>(2));
        int pos = FindIndexes(m_nodes, 0, GetNumNodes(), doubleMissingValue, indexes);
        indexes.resize(pos);

        
        int sizenewPolygon = GetNumNodes();
        if (innerAndOuter) 
        {
            sizenewPolygon += GetNumNodes() + 1;
        }
        
        std::vector<Point> normalVectors(sizenewPolygon);
        double dxNormalPreviusEdge;
        double dyNormalPreviusEdge;
        double dxNormalNodeIndex;
        double dyNormalNodeIndex;
        double dxNormalPreviusEdgeNodeIndex;
        double dyNormalPreviusEdgeNodeIndex;
        Point normalVectorNodeIndex;
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



};