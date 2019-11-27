#pragma once

#include <vector>
#include "Mesh.hpp"
#include "Constants.cpp"
#include "Operations.cpp"
#include "Polygons.hpp"
#include "Polygons.hpp"

namespace GridGeom
{
    class LandBoundaries
    {

    public:

        LandBoundaries(const Projections& projection) :
            m_numAllocatedNodes(0)
        {
            AllocateVector(m_numAllocatedNodes, m_nodes, { doubleMissingValue,doubleMissingValue });
            m_numNode = 0;
        }

        /// admin_landboundary_segments
        bool Set(const std::vector<Point>& landBoundary)
        {
            int newSize = m_numNode + landBoundary.size();

            bool successful = AllocateVector(newSize, m_nodes, { doubleMissingValue,doubleMissingValue });

            m_numAllocatedNodes = newSize;

            for (int i = 0; i < landBoundary.size(); i++)
            {
                m_nodes[m_numNode] = landBoundary[i];
                m_numNode++;
            }

            return true;
        };


        bool FindNearestMeshBoundary(const Mesh& mesh, const Polygons& polygon, bool meshBoundOnly)
        {
            bool successful = false;
            m_nodeMask.resize(mesh.m_nodes.size(), intMissingValue);
            m_faceMask.resize(mesh.m_numFaces, intMissingValue);
            m_edgeMask.resize(mesh.m_edges.size(), intMissingValue);
            polygonCache.resize(maximumNumberOfNodesPerFace + 1);
            m_landMask = true;
            m_nodeLandBoundarySegments.resize(mesh.m_nodes.size(), -1);
            m_nodesMinDistances.resize(mesh.m_nodes.size(), doubleMissingValue);

            for (int n = 0; n < m_numSegmentIndexses; n++)
            {
                int numPaths;
                int numRejectedPaths;
                
                successful = MakePath(mesh, polygon, n, meshBoundOnly, numPaths, numRejectedPaths);
                if (!successful)
                    return false;

                if (meshBoundOnly && numRejectedPaths > 0) 
                {
                    successful = MakePath(mesh, polygon, n, false, numPaths, numRejectedPaths);
                    if (!successful)
                        return false;
                }
            }
            //m_nodeLandBoundarySegments

            //if (meshBoundOnly) 
            //{
            //    for (int e = 0; e < mesh.m_edges.size(); e++)
            //    {
            //        successful = connect_boundary_paths(L, nodemask, 1, numnodes, nodelist);
            //        if (!successful)
            //            return false;
            //    }            
            //}

            return successful;
        };

        // admin_landboundary_segments
        bool Administrate(Mesh& mesh, Polygons& polygons)
        {

            std::vector<int> landBoundaryMask(m_numNode - 1, 0);
            //mask the landboundary that is inside the selecting polygon
            for (int n = 0; n < m_numNode - 1; n++)
            {
                if (m_nodes[n].x != doubleMissingValue && m_nodes[n + 1].x != doubleMissingValue)
                {
                    bool firstPointInPolygon = IsPointInPolygon(m_nodes[n], polygons.m_nodes, polygons.m_numNodes);
                    bool secondPointInPolygon = IsPointInPolygon(m_nodes[n + 1], polygons.m_nodes, polygons.m_numNodes);

                    if (firstPointInPolygon || secondPointInPolygon)
                    {
                        landBoundaryMask[n] = -1;
                    }
                }
            }

            // network boundary to polygon
            int counterClockWise = 0;
            int setMeshState = 0;
            std::vector<Point> meshBoundaryPolygon;
            int numNodesBoundaryPolygons;
            const bool successful = polygons.MeshBoundaryToPolygon(mesh, counterClockWise, setMeshState, meshBoundaryPolygon, numNodesBoundaryPolygons);
            if (!successful)
            {
                return false;
            }

            // Mask nodes close enough to land boundary segments 
            for (int n = 0; n < m_numNode - 1; n++)
            {
                if (landBoundaryMask[n] != 0)
                {

                    Point firstPoint = m_nodes[n];
                    Point secondPoint = m_nodes[n + 1];
                    const double landBoundaryLength = Distance(firstPoint, secondPoint, mesh.m_projection);

                    bool landBoundaryIsClose = false;
                    for (int nn = 0; nn < numNodesBoundaryPolygons - 1; nn++)
                    {
                        Point firstMeshBoundaryNode = meshBoundaryPolygon[nn];
                        Point secondMeshBoundaryNode = meshBoundaryPolygon[nn + 1];

                        if (firstMeshBoundaryNode.x == doubleMissingValue || secondMeshBoundaryNode.x == doubleMissingValue)
                        {
                            continue;
                        }

                        const double edgeLength = Distance(firstMeshBoundaryNode, secondMeshBoundaryNode, mesh.m_projection);
                        const double minDistance = m_closeToLandBoundary * edgeLength;

                        Point normalPoint;
                        double rlout;
                        const double distanceFirstMeshNode = DistanceFromLine(firstMeshBoundaryNode, firstPoint, secondPoint, normalPoint, rlout, mesh.m_projection);
                        const double distanceSecondMeshNode = DistanceFromLine(secondMeshBoundaryNode, firstPoint, secondPoint, normalPoint, rlout, mesh.m_projection);

                        if (distanceFirstMeshNode <= minDistance || distanceFirstMeshNode <= distanceSecondMeshNode)
                        {
                            landBoundaryIsClose = true;
                            break;
                        }
                    }

                    if (landBoundaryIsClose)
                    {
                        landBoundaryMask[n] = 1;
                    }
                }
            }

            // In m_segmentIndices store the starting and ending indexses of the land boundaries having
            // 1. the same landBoundaryMask
            // 2. are inside polygons
            m_segmentIndices.resize(m_numNode, std::vector<int>(2, -1));
            int begin = 0;
            int end = m_numNode + 1;
            while (begin < end)
            {
                int found = begin;
                for (int n = begin; n < end; n++)
                {
                    if (m_nodes[n].x == doubleMissingValue)
                    {
                        found = n;
                        break;
                    }
                }

                int startInd = begin;
                int endInd = found - 1;
                if (endInd > startInd && landBoundaryMask[startInd] != 0)
                {
                    m_segmentIndices[m_numSegmentIndexses] = { startInd, endInd };
                    m_numSegmentIndexses++;
                }
                begin = found + 1;
            }

            // Generate two segments for closed land boundaries
            int numSegmentIndexsesBeforeSplitting = m_numSegmentIndexses;
            for (int i = 0; i < numSegmentIndexsesBeforeSplitting; i++)
            {
                int start = m_segmentIndices[i][0];
                int end = m_segmentIndices[i][1];
                if (end > start)
                {
                    int split = int(start + (end - start) / 2);
                    m_segmentIndices[i][1] = split;
                    m_segmentIndices[m_numSegmentIndexses][0] = split;
                    m_segmentIndices[m_numSegmentIndexses][1] = end;
                    m_numSegmentIndexses++;
                }
            }

            return true;
        };


        /// Assigns to each node a land boundary a segment index (m_nodeLandBoundarySegments) 
        bool MakePath(const Mesh& mesh, const Polygons& polygons, 
            int landBoundarySegment, bool meshBoundOnly, int& numNodesInPath, int& numRejectedNodesInPath)
        {
            int startLandBoundaryIndex = m_segmentIndices[landBoundarySegment][0];
            int endLandBoundaryIndex = m_segmentIndices[landBoundarySegment][1];

            // these values are updated and used in FindStartEndNodes
            m_LeftIndex = endLandBoundaryIndex - 1;
            m_RightIndex = startLandBoundaryIndex;
            m_leftEdgeRatio = 1.0;
            m_rightEdgeRatio = 0.0;

            if (startLandBoundaryIndex < 0 || startLandBoundaryIndex >= m_numNode || startLandBoundaryIndex >= endLandBoundaryIndex)
                return false;

            bool successful = ComputeMask(mesh, polygons, landBoundarySegment, meshBoundOnly, startLandBoundaryIndex, endLandBoundaryIndex);
            if (!successful)
                return false;

            int startMeshNode;
            int endMeshNode;
            successful = FindStartEndMeshNodes(mesh, polygons, startLandBoundaryIndex, endLandBoundaryIndex, startMeshNode, endMeshNode);
            if (!successful)
                return false;

            if (startMeshNode < 0 || endMeshNode < 0 || startMeshNode == endMeshNode)
                return false;

            std::vector<int> connectedNodes;
            successful = ShortestPath(mesh, polygons, landBoundarySegment, startLandBoundaryIndex, endLandBoundaryIndex, startMeshNode, meshBoundOnly, connectedNodes);
            if (!successful)
                return false;

            int lastSegment = m_nodeLandBoundarySegments[endMeshNode];
            int lastNode = -1;
            int currentNode = endMeshNode;
            double minDinstanceFromLandBoundary;
            double distanceFromLandBoundary;
            Point nodeOnLandBoundary;
            int currentNodeLandBoundaryNodeIndex;
            double currentNodeEdgeRatio;
            bool stopPathSearch; //the path has been temporarily stopped(.true.) or not (.false.)
            int numConnectedNodes = 0;
            numRejectedNodesInPath = 0;
            numNodesInPath = 0;

            while (true)
            {
                stopPathSearch = true;

                if (m_nodeLandBoundarySegments[currentNode] >= 0)
                {
                    // Multiple boundary segments: take the nearest
                    int previousLandBoundarySegment = m_nodeLandBoundarySegments[currentNode];
                    int previousStartMeshNode = m_segmentIndices[previousLandBoundarySegment][0];
                    int previousEndMeshNode = m_segmentIndices[previousLandBoundarySegment][1];
                    int previousEndLandBoundaryIndex;
                    double previousMinDistance;
                  
                    successful = NearestPointOnLandBoundary(mesh, mesh.m_nodes[currentNode], previousStartMeshNode, previousEndMeshNode,
                        previousMinDistance, nodeOnLandBoundary, currentNodeLandBoundaryNodeIndex, currentNodeEdgeRatio);
                    if (!successful)
                        return false;

                    successful = NearestPointOnLandBoundary(mesh, mesh.m_nodes[currentNode], startLandBoundaryIndex, endLandBoundaryIndex,
                        distanceFromLandBoundary, nodeOnLandBoundary, currentNodeLandBoundaryNodeIndex, currentNodeEdgeRatio);
                    if (!successful)
                        return false;

                    double minDinstanceFromLandBoundaryCurrentNode = m_nodesMinDistances[currentNode];

                    if (distanceFromLandBoundary <= previousMinDistance && distanceFromLandBoundary < m_minDistanceFromLandFactor * minDinstanceFromLandBoundaryCurrentNode)
                    {
                        stopPathSearch = false;
                    }
                }
                else
                {
                    if (m_nodesMinDistances[currentNode] == doubleMissingValue)
                    {
                        successful = NearestPointOnLandBoundary(mesh, mesh.m_nodes[currentNode], 0, m_numNode,
                            minDinstanceFromLandBoundary, nodeOnLandBoundary, currentNodeLandBoundaryNodeIndex, currentNodeEdgeRatio);
                        if (!successful)
                            return false;
                        m_nodesMinDistances[currentNode] = minDinstanceFromLandBoundary;
                    }
                    else
                    {
                        minDinstanceFromLandBoundary = m_nodesMinDistances[currentNode];
                    }

                    successful = NearestPointOnLandBoundary(mesh, mesh.m_nodes[currentNode], startLandBoundaryIndex, endLandBoundaryIndex,
                        distanceFromLandBoundary, nodeOnLandBoundary, currentNodeLandBoundaryNodeIndex, currentNodeEdgeRatio);
                    if (!successful)
                        return false;

                    if (distanceFromLandBoundary < m_minDistanceFromLandFactor * minDinstanceFromLandBoundary)
                    {
                        if (meshBoundOnly == false || mesh.m_nodesTypes[currentNode] == 2 || mesh.m_nodesTypes[currentNode] == 3)
                        {
                            stopPathSearch = false;
                        }
                    }
                }

                if (stopPathSearch)
                {

                    if (numConnectedNodes == 1 && lastSegment != -1)
                    {
                        m_nodeLandBoundarySegments[lastNode] = lastSegment;
                    }
                    numConnectedNodes = 0;
                    numRejectedNodesInPath += 1;
                }
                else
                {
                    lastSegment = m_nodeLandBoundarySegments[currentNode];
                    lastNode = currentNode;

                    numNodesInPath += 1;
                    numConnectedNodes += 1;

                    m_nodeLandBoundarySegments[lastNode] = landBoundarySegment;
                }

                if (stopPathSearch && numConnectedNodes == 1)
                {
                    m_nodeLandBoundarySegments[lastSegment];
                }

                if (currentNode == startMeshNode)
                {
                    // Exit the while loop
                    break;
                }

                int nextEdgeIndex = connectedNodes[currentNode];
                if (nextEdgeIndex < 0 || nextEdgeIndex >= mesh.m_edges.size())
                {
                    break;
                }

                currentNode = mesh.m_edges[nextEdgeIndex].first + mesh.m_edges[nextEdgeIndex].second - currentNode;

                if (currentNode < 0 || currentNode >= mesh.m_nodes.size())
                {
                    break;
                }
            }

            if (numConnectedNodes == 1)
            {
                m_nodeLandBoundarySegments[lastNode] = lastSegment;
            }

            return true;
        }

        /// masks the nodes that are considered in the shortest path algorithm
        bool ComputeMask(const Mesh& mesh, const Polygons& polygon,
            int segmentIndex, bool meshBoundOnly, int& startLandBoundaryIndex, int& endLandBoundaryIndex)
        {
            std::fill(m_nodeMask.begin(), m_nodeMask.end(), doubleMissingValue);

            // check if any of the land boundary node is inside a mesh face
            bool nodeInFace = false;
            for (int i = startLandBoundaryIndex; i < endLandBoundaryIndex; i++)
            {
                for (int f = 0; f < mesh.m_numFaces; f++)
                {
                    auto numFaceNodes = mesh.m_facesNodes[f].size();

                    if (numFaceNodes == 0)
                        continue;

                    for (int n = 0; n < numFaceNodes; n++)
                    {
                        polygonCache[n] = mesh.m_nodes[mesh.m_facesNodes[f][n]];
                    }
                    polygonCache[numFaceNodes] = polygonCache[0];

                    nodeInFace = IsPointInPolygon(m_nodes[i], polygonCache, numFaceNodes);
                    if (nodeInFace)
                        break;
                }

                if (nodeInFace)
                    break;
            }

            // try to find a boundary cell that is crossed by the land boundary segment
            int crossedFaceIndex = -1;
            if (!nodeInFace)
            {
                for (int e = 0; e < mesh.m_edges.size(); e++)
                {
                    if (mesh.m_edgesNumFaces[e] != 1)
                        continue;

                    bool isCellCrossed = IsFaceCrossedByLandBoundaries(mesh, mesh.m_edgesFaces[e][0], startLandBoundaryIndex, endLandBoundaryIndex);

                    if (isCellCrossed)
                    {
                        crossedFaceIndex = mesh.m_edgesFaces[e][0];
                        break;
                    }
                }
            }

            if (m_landMask)
            {
                std::fill(m_faceMask.begin(), m_faceMask.end(), intMissingValue);
                std::fill(m_edgeMask.begin(), m_edgeMask.end(), intMissingValue);
                //m_faceMask assumes crossedFace has already been done.
                if (crossedFaceIndex >= 0 && crossedFaceIndex < mesh.m_numFaces)
                {
                    m_faceMask[crossedFaceIndex] = 1;
                }

                m_numFacesMasked = 0;
                m_maskDepth = 0;
                std::vector<int> landBoundaryFaces{ crossedFaceIndex };
                MaskFaces(mesh, meshBoundOnly, landBoundaryFaces, startLandBoundaryIndex, endLandBoundaryIndex);

                // Mask all nodes of the masked faces
                for (int f = 0; f < mesh.m_numFaces; f++)
                {
                    if (m_faceMask[f] == 1)
                    {
                        for (int n = 0; n < mesh.m_facesNodes[f].size(); n++)
                        {
                            m_nodeMask[mesh.m_facesNodes[f][n]] = segmentIndex;
                        }
                    }
                }
            }
            else
            {
                for (auto& e : m_nodeMask)
                    e = segmentIndex;
            }

            for (int n = 0; n < mesh.m_nodes.size(); n++)
            {
                if (m_nodeMask[n] > 0)
                {
                    bool inPolygon = IsPointInPolygon(mesh.m_nodes[n], polygon.m_nodes, polygon.m_numNodes);
                    if (!inPolygon)
                    {
                        m_nodeMask[n] = 0;
                    }
                }
            }

            return true;
        }

        /// maskcells, mask the cells that are intersected by the land boundary
        bool MaskFaces(const Mesh& mesh, const bool& meshBoundOnly, 
            std::vector<int>& landBoundaryFaces, int& startNodeLandBoundary, int &endNodeLandBoundary)
        {
            int numNextFaces = 0;
            std::vector<int> nextFaces(landBoundaryFaces.size(), intMissingValue);
            for (int f = 0; f < landBoundaryFaces.size(); f++)
            {
                int face = landBoundaryFaces[f];

                // no face was crossed by the land boundary: mask boundary faces only
                // these are the faces that are close (up to a certain tolerance) by a land boundary
                if (face < 0)
                {
                    int endNode = 0;
                    for (int e = 0; e < mesh.m_edges.size(); e++)
                    {
                        if (mesh.m_edgesNumFaces[e] != 1)
                            continue;

                        int otherFace = mesh.m_edgesFaces[e][0];

                        if (mesh.m_edges[e].first < 0 || mesh.m_edges[e].second < 0)
                            continue;

                        bool isClose = false;
                        int landBoundaryNode = 0;
                        for (int ee = 0; ee < mesh.m_facesEdges[otherFace].size(); ee++)
                        {
                            int edge = mesh.m_facesEdges[otherFace][ee];
                            isClose = IsEdgeCloseToLandBoundaries(mesh, edge, startNodeLandBoundary, endNodeLandBoundary, meshBoundOnly, landBoundaryNode);

                            if (isClose)
                            {
                                break;
                            }
                        }

                        if (isClose)
                        {
                            m_faceMask[otherFace] = 1;
                        }
                    }
                }
                else
                {
                    // face is crossed
                    if (mesh.m_facesEdges.size() < 3)
                        continue;

                    int isFaceFound = 0;
                    for (int e = 0; e <  mesh.m_facesEdges.size(); e++)
                    {
                        // is a boundary edge, continue
                        int currentEdge = mesh.m_facesEdges[face][e];
                        if (mesh.m_edgesNumFaces[currentEdge] <= 1)
                            continue;

                        int otherFace = mesh.m_edgesFaces[e][0] + mesh.m_edgesFaces[e][1] - face;

                        // already masked
                        if (m_faceMask[otherFace] != intMissingValue)
                            continue;

                        for (int ee = 0; ee < mesh.m_facesEdges[otherFace].size(); ee++)
                        {
                            int edgeOtherFace = mesh.m_facesEdges[otherFace][ee];

                            if (m_edgeMask[edgeOtherFace] == 1)
                            {
                                // previously visited crossed edge
                                isFaceFound = 1;
                            }
                            else if (m_edgeMask[edgeOtherFace] == 0)
                            {
                                // previously visited uncrossed edge - check next (otherFace) edge
                                continue;
                            }
                            else
                            {
                                // visited edge
                                m_edgeMask[edgeOtherFace] = 0;
                                int landBoundaryNode = 0;
                                bool isClose = IsEdgeCloseToLandBoundaries(mesh, edgeOtherFace, landBoundaryNode, startNodeLandBoundary, meshBoundOnly, endNodeLandBoundary);

                                if (isClose)
                                {
                                    m_edgeMask[edgeOtherFace] = 1;
                                    isFaceFound = 1;
                                }
                            }
                        }

                        // the rest of the loop
                        m_faceMask[otherFace] = isFaceFound;

                        if (isFaceFound == 1)
                        {
                            m_numFacesMasked += 1;
                            numNextFaces += 1;
                            if (numNextFaces > nextFaces.size())
                            {
                                nextFaces.resize(std::max(int(numNextFaces *1.2), 10));
                            }
                            nextFaces[numNextFaces] = otherFace;
                        }
                    }
                }
            }

            landBoundaryFaces.resize(0);

            if (numNextFaces > 0)
            {
                m_maskDepth += 1;
                MaskFaces(mesh, meshBoundOnly, nextFaces, startNodeLandBoundary, endNodeLandBoundary);
                m_maskDepth += 1;
            }

            return true;
        }

        /// linkcrossedbyland, check if a mesh edge is close to a land boundary segment
        bool IsEdgeCloseToLandBoundaries(const Mesh& mesh,
            int edgeIndex, int startNodeLandBoundary, int endNodeLandBoundary, bool meshBoundOnly, 
            int& landBoundaryNode)
        {
            bool isClose = false;

            const int startNode = std::max(std::min(landBoundaryNode, endNodeLandBoundary - 1), startNodeLandBoundary);

            if (mesh.m_edges[edgeIndex].first < 0 || mesh.m_edges[edgeIndex].second < 0)
                return false;

            const auto firstMeshNode = mesh.m_nodes[mesh.m_edges[edgeIndex].first];
            const auto secondMeshNode = mesh.m_nodes[mesh.m_edges[edgeIndex].second];

            const double meshEdgeLength = Distance(firstMeshNode, secondMeshNode, mesh.m_projection);

            const double distanceFactor = meshBoundOnly ? m_closeToLandBoundary : m_closeWholeMesh;
            const double closeDistance = meshEdgeLength * distanceFactor;

            double ratioFirstMeshNode;
            double ratioSecondMeshNode;
            Point normalPoint;

            int currentNode = startNode;
            int searchIter = 0;
            int stepNode = 0;
            while (searchIter < 3)
            {
                const double landBoundaryLength = Distance(m_nodes[currentNode], m_nodes[currentNode + 1], mesh.m_projection);        

                if (landBoundaryLength > 0 )
                {
                    ratioFirstMeshNode = 0.0;
                    ratioSecondMeshNode = 1.0;

                    const double distanceFromLandBoundaryFirstMeshNode = DistanceFromLine(firstMeshNode, m_nodes[currentNode], m_nodes[currentNode + 1], normalPoint, ratioFirstMeshNode, mesh.m_projection);

                    if (distanceFromLandBoundaryFirstMeshNode < closeDistance)
                    {
                        isClose = true;
                        landBoundaryNode = currentNode;
                        // the projection of firstMeshNode is within the segment currentNode / currentNode + 1
                        if (ratioFirstMeshNode >= 0.0 && ratioFirstMeshNode <= 1.0)
                        {
                            break;
                        }
                    }
                    else
                    {
                        // check the second point
                        double distanceFromLandBoundarySecondMeshNode = DistanceFromLine(secondMeshNode, m_nodes[currentNode], m_nodes[currentNode + 1], normalPoint, ratioSecondMeshNode, mesh.m_projection);

                        if (distanceFromLandBoundarySecondMeshNode < closeDistance)
                        {
                            isClose = true;
                            landBoundaryNode = currentNode;
                            // the projection of secondMeshNode is within the segment currentNode / currentNode + 1
                            if (ratioSecondMeshNode >= 0.0 && ratioSecondMeshNode <= 1.0)
                            {
                                break;
                            }
                        }
                    }
                }

                // search the next land boundary edge if projections are not within is within the segment currentNode / currentNode + 1
                searchIter = 0;
                while ((searchIter == 0 || currentNode < startNodeLandBoundary || currentNode > endNodeLandBoundary - 1) && searchIter < 3)
                {
                    searchIter += 1;
                    if (stepNode < 0)
                    {
                        stepNode = -stepNode + 1;
                    }
                    else
                    {
                        stepNode = -stepNode - 1;
                    }
                    currentNode = currentNode + stepNode;
                }
            }

            // the edge is close to the land boundary segment
            if (isClose)
            {
                // minimum
                double minimumRatio = std::min(ratioFirstMeshNode, ratioSecondMeshNode);
                if (landBoundaryNode < m_LeftIndex)
                {
                    m_LeftIndex = landBoundaryNode;
                    m_leftEdgeRatio = std::min(std::max(minimumRatio, 0.0), 1.0);
                }
                else if (landBoundaryNode == m_LeftIndex)
                {
                    m_leftEdgeRatio = std::min(std::max(minimumRatio, 0.0), m_leftEdgeRatio);
                }

                // maximum
                double maximumRatio = std::max(ratioFirstMeshNode, ratioSecondMeshNode);
                if (landBoundaryNode > m_RightIndex)
                {
                    m_RightIndex = landBoundaryNode;
                    m_rightEdgeRatio = std::min(std::max(maximumRatio, 0.0), 1.0);
                }
                else if (landBoundaryNode == m_RightIndex)
                {
                    m_rightEdgeRatio = std::max(std::min(maximumRatio, 1.0), m_rightEdgeRatio);
                }
            }

            return isClose;
        }

        /// get_kstartend2, find start and end nodes for a land boundary segment. These are nodes that are
        /// on a edge that is closest to the start and end node of the boundary segment respectively
        bool FindStartEndMeshNodes(const Mesh& mesh, const Polygons& polygons, int startLandBoundaryIndex, int endLandBoundaryIndex, int& startLandBoundaryNode, int& endLandBoundaryNode)
        {
            // compute the start and end point of the land boundary respectively
            int nextLeftIndex = std::min(m_LeftIndex + 1, endLandBoundaryIndex);
            Point startPoint =
            {
                m_nodes[m_LeftIndex].x + m_leftEdgeRatio * (m_nodes[nextLeftIndex].x - m_nodes[m_LeftIndex].x),
                m_nodes[m_LeftIndex].y + m_leftEdgeRatio * (m_nodes[nextLeftIndex].y - m_nodes[m_LeftIndex].y)
            };

            int nextRightIndex = std::min(m_RightIndex + 1, endLandBoundaryIndex);
            Point endPoint =
            {
                m_nodes[m_RightIndex].x + m_rightEdgeRatio * (m_nodes[nextRightIndex].x - m_nodes[m_RightIndex].x),
                m_nodes[m_RightIndex].y + m_rightEdgeRatio * (m_nodes[nextRightIndex].y - m_nodes[m_RightIndex].y)
            };

            bool isStartPointInsideAPolygon = IsPointInPolygons(startPoint, polygons.m_nodes, polygons.m_numNodes);
            bool isEndPointInsideAPolygon = IsPointInPolygons(endPoint, polygons.m_nodes, polygons.m_numNodes);

            if (!isStartPointInsideAPolygon)
            {
                startPoint.x = m_nodes[m_LeftIndex + 1].x;
                startPoint.y = m_nodes[m_LeftIndex + 1].y;
            }
            if (!isEndPointInsideAPolygon)
            {
                endPoint.x = m_nodes[m_RightIndex].x;
                endPoint.y = m_nodes[m_RightIndex].y;
            }

            // Get the edges that are closest the land boundary
            Point normalPoint;
            double ratio;

            double minDistStart = std::numeric_limits<double>::max();
            double minDistEnd = std::numeric_limits<double>::max();
            int startEdge = -1;
            int endEdge = -1;

            for (int e = 0; e < mesh.m_edges.size(); e++)
            {
                if (mesh.m_edges[e].first < 0 || mesh.m_edges[e].second < 0)
                    continue;

                if (m_nodeMask[mesh.m_edges[e].first] < 0 || m_nodeMask[mesh.m_edges[e].second] < 0)
                    continue;

                const double distanceFromFirstMeshNode = DistanceFromLine(startPoint, mesh.m_nodes[mesh.m_edges[e].first], mesh.m_nodes[mesh.m_edges[e].second], normalPoint, ratio, mesh.m_projection);
                const double distanceFromSecondMeshNode = DistanceFromLine(endPoint, mesh.m_nodes[mesh.m_edges[e].first], mesh.m_nodes[mesh.m_edges[e].second], normalPoint, ratio, mesh.m_projection);

                if (distanceFromFirstMeshNode < minDistStart)
                {
                    startEdge = e;
                    minDistStart = distanceFromFirstMeshNode;
                }
                if (distanceFromSecondMeshNode < minDistEnd)
                {
                    endEdge = e;
                    minDistEnd = distanceFromSecondMeshNode;
                }
            }

            if (startEdge == -1 || endEdge == -1)
            {
                return false;
            }

            // Find start and end node on the found edges
            int firstMeshNodeIndex = mesh.m_edges[startEdge].first;
            int secondMeshNodeIndex = mesh.m_edges[startEdge].second;
            double firstDinstance = Distance(mesh.m_nodes[firstMeshNodeIndex], startPoint, mesh.m_projection);
            double secondDinstance = Distance(mesh.m_nodes[secondMeshNodeIndex], startPoint, mesh.m_projection);

            if (firstDinstance <= secondDinstance)
            {
                startLandBoundaryNode = firstMeshNodeIndex;
            }
            else
            {
                startLandBoundaryNode = secondMeshNodeIndex;
            }

            firstMeshNodeIndex = mesh.m_edges[endEdge].first;
            secondMeshNodeIndex = mesh.m_edges[endEdge].second;
            firstDinstance = Distance(mesh.m_nodes[firstMeshNodeIndex], endPoint, mesh.m_projection);
            secondDinstance = Distance(mesh.m_nodes[secondMeshNodeIndex], endPoint, mesh.m_projection);

            if (firstDinstance <= secondDinstance)
            {
                endLandBoundaryNode = firstMeshNodeIndex;
            }
            else
            {
                endLandBoundaryNode = secondMeshNodeIndex;
            }

            return true;
        }



        // Shortest_path
        // connect mesh nodes starting from startMeshNode, using Dijkstra's shortest path algorithm
        // the distance of each edge is the edge length multiplied by the distance from the land boundary
        bool ShortestPath(const Mesh& mesh, const Polygons& polygons, int landBoundarySegment,
            int startLandBoundaryIndex, int endLandBoundaryIndex, int startMeshNode, bool meshBoundOnly, std::vector<int>& connectedNodes)
        {
            connectedNodes.resize(mesh.m_nodes.size(), -1);
            // infinite distance for all nodes
            std::vector<double> nodeDistances(mesh.m_nodes.size(), std::numeric_limits<double>::max());
            std::vector<bool> isVisited(mesh.m_nodes.size(), false);

            int currentNodeIndex = startMeshNode;
            nodeDistances[startMeshNode] = 0.0;
            while (true)
            {
                isVisited[currentNodeIndex] = true;
                Point currentNode = mesh.m_nodes[currentNodeIndex];
                Point currentNodeOnLandBoundary;
                int currentNodeLandBoundaryNodeIndex;
                double currentNodeEdgeRatio;
                double currentNodeDistance;
                bool successful = NearestPointOnLandBoundary(mesh, currentNode, startLandBoundaryIndex, endLandBoundaryIndex,
                    currentNodeDistance, currentNodeOnLandBoundary, currentNodeLandBoundaryNodeIndex, currentNodeEdgeRatio);

                if (!successful || currentNodeLandBoundaryNodeIndex < 0)
                    return false;

                for (int e = 0; e < mesh.m_nodesEdges[currentNodeIndex].size(); e++)
                {
                    int edgeIndex = mesh.m_nodesEdges[currentNodeIndex][e];

                    if (mesh.m_edges[edgeIndex].first < 0 || mesh.m_edges[edgeIndex].second < 0)
                    {
                        continue;
                    }

                    int neighbouringNodeIndex = mesh.m_edges[edgeIndex].first + mesh.m_edges[edgeIndex].second - currentNodeIndex;

                    if (isVisited[neighbouringNodeIndex])
                    {
                        continue;
                    }

                    Point neighbouringNode = mesh.m_nodes[neighbouringNodeIndex];

                    Point neighbouringNodeOnLandBoundary;
                    int neighbouringNodeLandBoundaryNodeIndex;
                    double neighbouringNodeEdgeRatio;
                    double neighbouringNodeDistance;
                    successful = NearestPointOnLandBoundary(mesh, neighbouringNode, startLandBoundaryIndex, endLandBoundaryIndex,
                        neighbouringNodeDistance, neighbouringNodeOnLandBoundary, neighbouringNodeLandBoundaryNodeIndex, neighbouringNodeEdgeRatio);
                    if (!successful)
                        return false;

                    Point middlePoint
                    {
                        (currentNode.x + neighbouringNode.x) * 0.5,
                        (currentNode.y + neighbouringNode.y) * 0.5
                    };

                    Point middlePointOnLandBoundary;
                    int middlePointLandBoundaryNodeIndex;
                    double middlePointEdgeRatio;
                    double middlePointDistance;
                    successful = NearestPointOnLandBoundary(mesh, middlePoint, startLandBoundaryIndex, endLandBoundaryIndex,
                        middlePointDistance, middlePointOnLandBoundary, middlePointLandBoundaryNodeIndex, middlePointEdgeRatio);
                    if (!successful)
                        return false;


                    double firstNodeDistance = Distance(currentNodeOnLandBoundary, middlePointOnLandBoundary, mesh.m_projection);
                    double secondNodeDistance = Distance(neighbouringNodeOnLandBoundary, middlePointOnLandBoundary, mesh.m_projection);
                    double maximumDistance = std::max(currentNodeDistance, neighbouringNodeDistance);

                    if (currentNodeLandBoundaryNodeIndex < neighbouringNodeLandBoundaryNodeIndex)
                    {
                        for (int n = currentNodeLandBoundaryNodeIndex + 1; n < neighbouringNodeLandBoundaryNodeIndex; n++)
                        {
                            double ratio;
                            Point projectedPoint;
                            middlePointDistance = DistanceFromLine(m_nodes[n], currentNode, neighbouringNode, middlePointOnLandBoundary, ratio, mesh.m_projection);
                            if (middlePointDistance > maximumDistance)
                                maximumDistance = middlePointDistance;
                        }
                    }
                    else if (currentNodeLandBoundaryNodeIndex > neighbouringNodeLandBoundaryNodeIndex)
                    {
                        for (int n = neighbouringNodeLandBoundaryNodeIndex + 1; n < currentNodeLandBoundaryNodeIndex; n++)
                        {
                            double ratio;
                            Point projectedPoint;
                            middlePointDistance = DistanceFromLine(m_nodes[n], currentNode, neighbouringNode, middlePointOnLandBoundary, ratio, mesh.m_projection);
                            if (middlePointDistance > maximumDistance)
                                maximumDistance = middlePointDistance;
                        }
                    }

                    // In case of netboundaries only: set penalty when edge is not on the boundary
                    if (meshBoundOnly && mesh.m_edgesNumFaces[edgeIndex] != 1)
                    {
                        maximumDistance = 1e6 * maximumDistance;
                    }

                    double edgeLength = Distance(currentNode, neighbouringNode, mesh.m_projection);
                    double correctedDistance = nodeDistances[currentNodeIndex] + edgeLength * maximumDistance;

                    if (correctedDistance < nodeDistances[neighbouringNodeIndex])
                    {
                        nodeDistances[neighbouringNodeIndex] = correctedDistance;
                        connectedNodes[neighbouringNodeIndex] = edgeIndex;
                    }
                }

                // linear search with masking
                currentNodeIndex = 0;
                double minValue = std::numeric_limits<double>::max();
                for (int n = 0; n < mesh.m_nodes.size(); n++)
                {
                    if (m_nodeMask[n] == landBoundarySegment && !isVisited[n] && nodeDistances[n] < minValue)
                    {
                        currentNodeIndex = n;
                        minValue = nodeDistances[n];

                    }
                }

                if (currentNodeIndex < 0 ||
                    currentNodeIndex >= mesh.m_nodes.size() ||
                    nodeDistances[currentNodeIndex] == std::numeric_limits<double>::max() ||
                    isVisited[currentNodeIndex])
                {
                    break;
                }
            }
            return true;
        }

        bool NearestPointOnLandBoundary(const Mesh& mesh, const Point& node, int startLandBoundaryIndex, int endLandBoundaryIndex,
            double& minimumDistance, Point& pointOnLandBoundary, int& nearestLandBoundaryNodeIndex, double& edgeRatio)
        {

            minimumDistance = std::numeric_limits<double>::max();
            nearestLandBoundaryNodeIndex = -1;
            edgeRatio = -1.0;
            pointOnLandBoundary = node;
            for (int n = startLandBoundaryIndex; n < endLandBoundaryIndex; n++)
            {
                if (m_nodes[n].x == doubleMissingValue || m_nodes[n + 1].x == doubleMissingValue)
                    continue;

                Point normalPoint;
                double ratio;
                const double distanceFromLandBoundary = DistanceFromLine(node, m_nodes[n], m_nodes[n + 1], normalPoint, ratio, mesh.m_projection);

                if (distanceFromLandBoundary > 0.0 && distanceFromLandBoundary < minimumDistance)
                {
                    minimumDistance = distanceFromLandBoundary;
                    pointOnLandBoundary = normalPoint;
                    nearestLandBoundaryNodeIndex = n;
                    edgeRatio = ratio;
                }
            }
            return true;
        }


        // cellcrossedbyland
        // TODO: it could be moved to generic operations
        bool IsFaceCrossedByLandBoundaries(const Mesh& mesh, int face, int startLandBoundaryIndex, int endLandBoundaryIndex)
        {
            bool isCrossed = false;
            bool areCrossing = false;

            for (int e = 0; e < mesh.m_facesEdges[face].size(); e++)
            {
                auto edge = mesh.m_facesEdges[face][e];

                for (int i = startLandBoundaryIndex; i < endLandBoundaryIndex; i++)
                {
                    auto firstMeshNode = mesh.m_nodes[mesh.m_edges[edge].first];
                    auto secondMeshNode = mesh.m_nodes[mesh.m_edges[edge].second];

                    auto firstNode = m_nodes[i];
                    auto secondNode = m_nodes[i + 1];
                    bool adimensional = false;
                    Point intersection;
                    double crossProduct;

                    bool areCrossing = AreLinesCrossing(firstMeshNode, secondMeshNode, firstNode, secondNode, adimensional, intersection, crossProduct, mesh.m_projection);

                    if (areCrossing)
                    {
                        return true;
                    }
                }
            }
            return false;
        }


        std::vector<Point> m_nodes;                     // XLAN, YLAN, ZLAN
        int m_numAllocatedNodes;                        // MAXLAN
        int m_numNodesLoc;                              // MXLAN_loc

        // number of land boundary segments
        int m_numSegments;                               // Nlanseg
        int m_numNode;                                   // actual MXLAN

        int m_numSegmentIndexses = 0;
        std::vector<std::vector<int>> m_segmentIndices;  // lanseg_startend
        std::vector<std::vector<double>> m_nodesLand;    // !node to land boundary segment mapping
        std::vector<int> m_nodeMask;                     // masking the net nodes
        std::vector<int> m_faceMask;                     // masking faces
        std::vector<int> m_edgeMask;                     // masking edges
        std::vector<Point> polygonCache;                 // array of points (e.g. points of a face)

        std::vector<int> m_nodeLandBoundarySegments;
        std::vector<double> m_nodesMinDistances;

        bool m_landMask = true;
        int m_numFacesMasked = 0;
        int m_maskDepth = 0;

        int m_LeftIndex;                                // outer land boundary segments in projection
        int m_RightIndex;
        double m_leftEdgeRatio;                         // fractional location of the projected outer nodes(min and max) on the land boundary segment
        double m_rightEdgeRatio;

        // Parameters
        const double m_closeToLandBoundary = 5.0;      // close - to - landboundary tolerance, measured in number of meshwidths
        const double m_closeWholeMesh = 1.0;           // close - to - landboundary tolerance, measured in number of meshwidths
        const double m_minDistanceFromLandFactor = 2.0;
    };
}
