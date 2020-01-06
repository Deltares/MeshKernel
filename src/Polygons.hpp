#pragma once

#include <utility>
#include <vector>
#include "Mesh.hpp"
#include "Constants.cpp"
#include "Operations.cpp"

namespace GridGeom
{
    class Polygons
    {

    public:

        Polygons() : m_numNodes(0), m_numAllocatedNodes(0)
        {
            AllocateVector(m_numAllocatedNodes, m_nodes);
        }

        bool Set(const std::vector<Point>& polygon)
        {
            // resize if necessary
            ResizeVector(m_numNodes + 1 + polygon.size(), m_nodes);

            m_nodes[m_numNodes] = { doubleMissingValue,doubleMissingValue };
            for (int n = m_numNodes + 1, nn = 0; nn < polygon.size(); nn++)
            {
                m_nodes[n] = polygon[nn];
            }
            return true;
        }

        /// copynetboundstopol
        bool MeshBoundaryToPolygon(const Mesh& mesh,
            int counterClockWise,
            int setMeshState,
            std::vector<Point>& meshBoundaryPolygon,
            int& numNodesBoundaryPolygons)
        {
            std::vector<bool> isVisited(mesh.m_edges.size());
            std::vector<int> boundaryPolygonStarts(mesh.m_edges.size());
            int numBoundaryPolygons = 0;
            numNodesBoundaryPolygons = 0;

            meshBoundaryPolygon.resize(mesh.m_nodes.size(), { doubleMissingValue ,doubleMissingValue });
            int meshBoundaryPolygonSize = meshBoundaryPolygon.size();

            for (int e = 0; e < mesh.m_edges.size(); e++)
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

        std::vector<Point> m_nodes;             // Polygon nodes
        int m_numNodes;                         // NPL
        int m_numAllocatedNodes;                // MAXPOL

    private:

        bool WalkBoundary(const Mesh& mesh,
            std::vector<bool>& isVisited,
            int& numNodesBoundaryPolygon,
            int& currentNode,
            int meshBoundaryPolygonSize,
            std::vector<Point>& meshBoundaryPolygon)
        {
            int ee = 0;
            while(ee < mesh.m_nodesNumEdges[currentNode])
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

                const int firstCurrentEdge = mesh.m_edges[currentEdge].first;
                const int secondCurrentEdge = mesh.m_edges[currentEdge].second;

                if (secondCurrentEdge == currentNode)
                {
                    currentNode = firstCurrentEdge;
                    ee = 0;
                }
                else
                {
                    currentNode = secondCurrentEdge;
                    ee = 0;
                }

                numNodesBoundaryPolygon++;

                ResizeVector(numNodesBoundaryPolygon, meshBoundaryPolygon);

                meshBoundaryPolygon[numNodesBoundaryPolygon] = mesh.m_nodes[currentNode];

                isVisited[currentEdge] = true;
            }
            return true;
        }
    };

}