#pragma once

#include <utility>
#include <vector>
#include <algorithm>
#include "Mesh.hpp"
#include "Constants.cpp"
#include "Operations.cpp"

namespace GridGeom
{
    class Polygons
    {

    public:

        Polygons() : m_numPolygonNodes(0), m_allocatedSize(0)
        {
            increasePolygon(0);
        }

        bool SetPolygon(const std::vector<Point>& polygon)
        {
            // resize if necessary
            resizeVector(m_numPolygonNodes + 1 + polygon.size(), m_polygonNodes, { doubleMissingValue ,doubleMissingValue });

            m_polygonNodes[m_numPolygonNodes] = { doubleMissingValue,doubleMissingValue };
            for (int n = m_numPolygonNodes + 1, nn = 0; nn < polygon.size(); nn++)
            {
                m_polygonNodes[n] = polygon[nn];
            }
            return true;
        }


        // copynetboundstopol
        bool MeshBoundaryToPolygon(Mesh& mesh,
            int counterClockWise,
            int setMeshState,
            std::vector<Point>& meshBoundaryPolygon)
        {
            std::vector<bool> isVisited(mesh.m_edges.size());
            std::vector<int> startingPolygonEdges(mesh.m_edges.size());

            meshBoundaryPolygon.resize(mesh.m_nodes.size(), { doubleMissingValue ,doubleMissingValue });
            int meshBoundaryPolygonSize = meshBoundaryPolygon.size();
            int numNodesBoundaryPolygon = 0;
            int numPolygonEdges = 0;

            for (int e = 0; e < mesh.m_edges.size(); e++)
            {
                if (isVisited[e])
                {
                    continue;
                }

                if (mesh.m_edgesNumFaces[e] != 1)
                {
                    continue;
                }

                int first = mesh.m_edges[e].first;
                int second = mesh.m_edges[e].second;
                auto firstPoint = mesh.m_nodes[first];
                auto secondPoint = mesh.m_nodes[second];

                bool inHullFirst = pointInPolygon(mesh.m_nodes[first], m_polygonNodes, m_numPolygonNodes);
                bool inHullSecond = pointInPolygon(mesh.m_nodes[second], m_polygonNodes, m_numPolygonNodes);

                if (!inHullFirst && !inHullSecond)
                {
                    continue;
                }

                resizeVector(numNodesBoundaryPolygon + 3, meshBoundaryPolygon, { doubleMissingValue, doubleMissingValue });

                //Start a new polyline
                if (numNodesBoundaryPolygon > 0)
                {
                    numNodesBoundaryPolygon++;
                }

                startingPolygonEdges[numPolygonEdges] = numNodesBoundaryPolygon;

                const int startPolygonEdges = numNodesBoundaryPolygon;
                const int nodeStart = first;


                meshBoundaryPolygon[numNodesBoundaryPolygon] = firstPoint;
                numNodesBoundaryPolygon++;
                meshBoundaryPolygon[numNodesBoundaryPolygon] = secondPoint;

                isVisited[e] = true;
                numPolygonEdges++;

                int currentNode = second;

                walkBoundary(mesh, isVisited, numNodesBoundaryPolygon, currentNode, meshBoundaryPolygonSize, meshBoundaryPolygon);

                int numNodesFirstTail = numNodesBoundaryPolygon;

                // polyline is closed itself
                if (currentNode == nodeStart)
                {

                }
                else
                {
                    //Now grow a polyline starting at the other side of the original link L, i.e., the second tail
                    currentNode = nodeStart;
                    walkBoundary(mesh, isVisited, numNodesBoundaryPolygon, currentNode, meshBoundaryPolygonSize, meshBoundaryPolygon);
                }

                // There is a nonempty second tail, so reverse the first tail, so that they connect.
                if (numNodesBoundaryPolygon > numNodesFirstTail)
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

            startingPolygonEdges[numPolygonEdges + 1] = numNodesBoundaryPolygon + 2;

            return true;
        }

    private:

        bool walkBoundary(const Mesh& mesh,
            std::vector<bool>& isVisited,
            int& numNodesBoundaryPolygon,
            int& currentNode,
            int meshBoundaryPolygonSize,
            std::vector<Point>& meshBoundaryPolygon)
        {
            int ee = 0;
            while(ee < mesh.m_nodesNumEdges[currentNode])
            {
                bool inHull = pointInPolygon(mesh.m_nodes[currentNode], m_polygonNodes, m_numPolygonNodes);

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

                resizeVector(numNodesBoundaryPolygon, meshBoundaryPolygon, { doubleMissingValue , doubleMissingValue });

                meshBoundaryPolygon[numNodesBoundaryPolygon] = mesh.m_nodes[currentNode];

                isVisited[currentEdge] = true;
            }
            return true;
        }

        template<typename T>
        bool resizeVector(int newSize, std::vector<T>& vectorToResize, T fillValue = T())
        {
            int currentSize = vectorToResize.size();
            if (newSize > currentSize)
            {
                currentSize = std::max(newSize, int(currentSize * 1.2));
                vectorToResize.resize(currentSize, fillValue);
            }
            return true;
        }

        bool increasePolygon(int newSize)
        {
            int currentSize = m_polygonNodes.size();
            if (newSize < currentSize)
            {
                return true;
            }

            m_allocatedSize = std::max(10000, int(5 * newSize));
            m_polygonNodes.resize(m_allocatedSize, { doubleMissingValue,doubleMissingValue });
            return true;
        }


        std::vector<Point> m_polygonNodes;  // Polygon nodes
        int m_numPolygonNodes;              // NPL
        int m_allocatedSize;                // MAXPOL

    };

}