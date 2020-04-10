#pragma once

#include <vector>
#include <iostream>
#include "Mesh.hpp"
#include "Entities.hpp"
#include "Constants.cpp"
#include "Operations.cpp"
#include "Polygons.hpp"
#include "LandBoundaries.hpp"

namespace GridGeom
{
    LandBoundaries::LandBoundaries() :
        m_numAllocatedNodes(0)
    {
        AllocateVector(m_numAllocatedNodes, m_nodes, m_allocationSize, { doubleMissingValue,doubleMissingValue });
        m_numAllocatedNodes = m_nodes.size();
        m_numNode = 0;
    }

    /// admin_landboundary_segments
    bool LandBoundaries::Set(const std::vector<Point>& landBoundary)
    {

        bool successful = AllocateVector(m_numNode + landBoundary.size(), m_nodes, m_allocationSize,{ doubleMissingValue,doubleMissingValue });
        m_numAllocatedNodes = m_nodes.size();

        for (int i = 0; i < landBoundary.size(); i++)
        {
            m_nodes[m_numNode] = landBoundary[i];
            m_numNode++;
        }

        return true;
    };

    /// admin_landboundary_segments
    /// The land boundary will be split into segments that are within the polygon, and either close or not to the mesh boundary
    /// TODO: ? Why splitting in two segments is required?
    bool LandBoundaries::Administrate(const Mesh& mesh, Polygons& polygons)
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
                    const double minDistance = m_closeToLandBoundaryFactor * edgeLength;

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

        // In m_segmentIndices: start and ending indexses of segments inside polygon
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
                m_segmentIndices[m_numSegments] = { startInd, endInd };
                m_numSegments++;
            }
            begin = found + 1;
        }

        // Generate two segments for closed land boundaries
        int numSegmentIndexsesBeforeSplitting = m_numSegments;
        for (int i = 0; i < numSegmentIndexsesBeforeSplitting; i++)
        {
            int start = m_segmentIndices[i][0];
            int end = m_segmentIndices[i][1];
            if (end > start)
            {
                int split = int(start + (end - start) / 2);
                m_segmentIndices[i][1] = split;
                m_segmentIndices[m_numSegments][0] = split;
                m_segmentIndices[m_numSegments][1] = end;
                m_numSegments++;
            }
        }

        return true;
    };


    /// find_nearest_meshline
    /// Find the mesh boundary line closest to the land boundary
    bool LandBoundaries::FindNearestMeshBoundary(const Mesh& mesh, const Polygons& polygon, int snapping)
    {
        bool successful = false;
        bool meshBoundOnly = false;

        if (snapping == 2 || snapping == 3)
        {
            meshBoundOnly = true;
        }

        m_nodeMask.resize(mesh.m_nodes.size(), intMissingValue);
        m_faceMask.resize(mesh.m_numFaces, intMissingValue);
        m_edgeMask.resize(mesh.m_edges.size(), intMissingValue);
        polygonCache.resize(maximumNumberOfNodesPerFace + 1);
        m_meshNodesLandBoundarySegments.resize(mesh.m_nodes.size(), -1);
        m_nodesMinDistances.resize(mesh.m_nodes.size(), doubleMissingValue);

        //landBoundarySegment
        for (int landBoundarySegment = 0; landBoundarySegment < m_numSegments; landBoundarySegment++)
        {
            int numPaths;
            int numRejectedPaths;
            successful = MakePath(mesh, polygon, landBoundarySegment, meshBoundOnly, numPaths, numRejectedPaths);

            if (!successful)
                return false;

            if (numRejectedPaths > 0 && snapping == 3)
            {
                successful = MakePath(mesh, polygon, landBoundarySegment, false, numPaths, numRejectedPaths);
                if (!successful)
                    return false;
            }
        }

        // connect the mesh nodes
        if (meshBoundOnly)
        {
            std::vector<int> connectedNodes;
            int numConnectedNodes = 0;
            for (int e = 0; e < mesh.m_edges.size(); e++)
            {
                if (mesh.m_edgesNumFaces[e] != 1)
                    continue;
                successful = AssignSegmentsToAllMeshNodes(mesh, e, true, connectedNodes, numConnectedNodes);
                if (!successful)
                    return false;
            }
        }

        return successful;
    };


    /// connect_boundary_paths, build an additional boundary for not assigned nodes  
    bool LandBoundaries::AssignSegmentsToAllMeshNodes(const Mesh& mesh, int edgeIndex, bool initialize, std::vector<int>& nodes, int numNodes)
    {
        bool successful = false;
        std::vector<int> nodesLoc;
        int numNodesLoc;

        if (initialize)
        {
            if (mesh.m_edgesNumFaces[edgeIndex] != 1 || mesh.m_edges[edgeIndex].first < 0 || mesh.m_edges[edgeIndex].second < 0)
                return true;

            int firstMeshNode = mesh.m_edges[edgeIndex].first;
            int secondMeshNode = mesh.m_edges[edgeIndex].second;

            if (m_meshNodesLandBoundarySegments[firstMeshNode] >= 0 &&
                m_meshNodesLandBoundarySegments[secondMeshNode] < 0 &&
                m_nodeMask[firstMeshNode] > 0 &&
                m_nodeMask[secondMeshNode] > 0)
            {
                nodesLoc.resize(3);
                nodesLoc[0] = firstMeshNode;
                nodesLoc[1] = secondMeshNode;
                numNodesLoc = 2;
            }
            else if (m_meshNodesLandBoundarySegments[firstMeshNode] < 0 &&
                m_meshNodesLandBoundarySegments[secondMeshNode] >= 0 &&
                m_nodeMask[firstMeshNode] > 0 &&
                m_nodeMask[secondMeshNode] > 0)
            {
                nodesLoc.resize(3);
                nodesLoc[0] = secondMeshNode;
                nodesLoc[1] = firstMeshNode;
                numNodesLoc = 2;
            }
            else
            {
                //not a valid edge
                return true;
            }
        }
        else
        {
            nodesLoc.resize(numNodes + 1);
            numNodesLoc = numNodes;
            std::copy(nodes.begin(), nodes.end(), nodesLoc.begin());
        }

        int maxNodes = *std::max_element(nodesLoc.begin(), nodesLoc.end() - 1);
        if (numNodesLoc > maxNodes)
            return true;

        int lastVisitedNode = nodesLoc[numNodesLoc - 1];

        for (int e = 0; e < mesh.m_nodesNumEdges[lastVisitedNode]; e++)
        {
            int edge = mesh.m_nodesEdges[lastVisitedNode][e];

            if (mesh.m_edgesNumFaces[edge] != 1)
                continue;

            int otherNode = mesh.m_edges[edge].first + mesh.m_edges[edge].second - lastVisitedNode;

            // path stopped
            if (m_nodeMask[otherNode] < 0)
                break;

            bool otherNodeAlreadyVisited = false;
            for (int n = numNodesLoc - 1; n >= 0; n--)
            {
                if (otherNode == nodesLoc[n])
                {
                    otherNodeAlreadyVisited = true;
                    break;
                }
            }

            if (otherNodeAlreadyVisited)
                continue;

            nodesLoc[numNodesLoc] = otherNode;

            if (m_meshNodesLandBoundarySegments[otherNode] >= 0)
            {
                // Now check landboundary for otherNode
                for (int n = 1; n < numNodesLoc; n++)
                {
                    int meshNode = nodesLoc[n];
                    double minimumDistance;
                    Point pointOnLandBoundary;
                    int nearestLandBoundaryNodeIndex;
                    double edgeRatio;

                    successful = NearestLandBoundaryNode(mesh.m_projection, mesh.m_nodes[meshNode], 0, m_numNode,
                        minimumDistance, pointOnLandBoundary, nearestLandBoundaryNodeIndex, edgeRatio);

                    // find the segment index of the found point
                    int landboundarySegmentIndex = -1;
                    for (int s = 0; s < m_numSegments; s++)
                    {
                        if (nearestLandBoundaryNodeIndex >= m_segmentIndices[s][0] && nearestLandBoundaryNodeIndex < m_segmentIndices[s][1])
                        {
                            landboundarySegmentIndex = s;
                            break;
                        }
                    }

                    //TODO: maybe this is not completly correct
                    if (landboundarySegmentIndex == -1)
                        return false;

                    if (meshNode == 0) 
                    {
                        std::cout << "D";
                    
                    }

                    if ((nearestLandBoundaryNodeIndex == m_segmentIndices[landboundarySegmentIndex][0] && edgeRatio < 0.0) ||
                        (nearestLandBoundaryNodeIndex == m_segmentIndices[landboundarySegmentIndex][1] - 1 && edgeRatio > 1.0))
                    {
                        if (m_addLandboundaries)
                        {
                            AddLandBoundary(mesh, nodesLoc, numNodesLoc, lastVisitedNode);
                            m_meshNodesLandBoundarySegments[meshNode] = m_numSegments - 1; //last added ;and boundary
                        }
                    }
                    else
                    {
                        m_meshNodesLandBoundarySegments[meshNode] = landboundarySegmentIndex;
                    }
                }
            }
            else
            {
                AssignSegmentsToAllMeshNodes(mesh, edge, false, nodesLoc, numNodesLoc + 1);
            }
        }
        return true;
    }

    /// add_land, add new land boundary segment that connects two others
    bool LandBoundaries::AddLandBoundary(const Mesh& mesh, const std::vector<int>& nodesLoc, int numNodesLoc, int nodeIndex)
    {
        bool successful = false;

        int startSegmentIndex = m_meshNodesLandBoundarySegments[nodesLoc[0]];
        int endSegmentIndex = m_meshNodesLandBoundarySegments[nodesLoc[numNodesLoc]];

        if (startSegmentIndex < 0 || startSegmentIndex >= m_numSegments ||
            endSegmentIndex < 0 || endSegmentIndex >= m_numSegments)
        {
            return false;
        }

        // find start/end
        int startNode = m_segmentIndices[startSegmentIndex][0];
        int endNode = m_segmentIndices[startSegmentIndex][1];

        Point newNodeLeft;
        if (Distance(mesh.m_nodes[nodeIndex], m_nodes[startNode], mesh.m_projection) <= Distance(mesh.m_nodes[nodeIndex], m_nodes[endNode], mesh.m_projection))
        {
            newNodeLeft = m_nodes[startNode];
        }
        else
        {
            newNodeLeft = m_nodes[endNode];
        }

        Point newNodeRight;
        if (endSegmentIndex == startSegmentIndex)
        {
            newNodeRight = m_nodes[startNode] + m_nodes[endNode] - newNodeLeft;
        }
        else
        {
            // find start/end
            startNode = m_segmentIndices[endSegmentIndex][0];
            endNode = m_segmentIndices[endSegmentIndex][1];
            if (Distance(mesh.m_nodes[nodeIndex], m_nodes[startNode], mesh.m_projection) <= Distance(mesh.m_nodes[nodeIndex], m_nodes[endNode], mesh.m_projection))
            {
                newNodeRight = m_nodes[startNode];
            }
            else
            {
                newNodeRight = m_nodes[endNode];
            }
        }

        // Update  nodes
        int minSize = m_numNode + 3;
        AllocateVector(minSize, m_nodes, m_allocationSize, { doubleMissingValue,doubleMissingValue });
        m_numNode += 3;
        m_nodes[m_numNode - 2] = newNodeLeft;
        m_nodes[m_numNode - 1] = newNodeRight;

        // Update segment indexses
        minSize = m_numSegments + 1;
        AllocateVector(minSize, m_segmentIndices, m_allocationSize, { -1, -1 });
        m_segmentIndices[m_numSegments][0] = m_numNode - 2;
        m_segmentIndices[m_numSegments][1] = m_numNode - 1;
        m_numSegments++;

        return successful;
    }

    /// make_path
    /// Assigns to each mesh node a land boundary a segment index (m_nodeLandBoundarySegments) 
    bool LandBoundaries::MakePath(const Mesh& mesh,
        const Polygons& polygons,
        int landBoundarySegment,
        bool meshBoundOnly,
        int& numNodesInPath,
        int& numRejectedNodesInPath)
    {
        int startLandBoundaryIndex = m_segmentIndices[landBoundarySegment][0];
        int endLandBoundaryIndex = m_segmentIndices[landBoundarySegment][1];

        if (startLandBoundaryIndex < 0 || startLandBoundaryIndex >= m_numNode || startLandBoundaryIndex >= endLandBoundaryIndex)
            return false;

        // fractional location of the projected outer nodes(min and max) on the land boundary segment
        double leftEdgeRatio = 1.0;
        double rightEdgeRatio = 0.0;
        int leftIndex = endLandBoundaryIndex - 1;
        int rightIndex = startLandBoundaryIndex;

        bool successful = ComputeMask(mesh,
            polygons,
            landBoundarySegment,
            meshBoundOnly,
            leftIndex,
            rightIndex,
            leftEdgeRatio,
            rightEdgeRatio,
            startLandBoundaryIndex,
            endLandBoundaryIndex);

        if (!successful)
            return false;

        int startMeshNode;
        int endMeshNode;
        successful = FindStartEndMeshNodes(mesh,
            polygons,
            startLandBoundaryIndex,
            endLandBoundaryIndex,
            leftIndex,
            rightIndex,
            leftEdgeRatio,
            rightEdgeRatio,
            startMeshNode,
            endMeshNode);

        if (!successful)
            return false;

        if (startMeshNode < 0 || endMeshNode < 0 || startMeshNode == endMeshNode)
            return false;

        std::vector<int> connectedNodes;
        successful = ShortestPath(mesh, polygons, landBoundarySegment, startLandBoundaryIndex, endLandBoundaryIndex, startMeshNode, meshBoundOnly, connectedNodes);
        if (!successful)
            return false;

        int lastSegment = m_meshNodesLandBoundarySegments[endMeshNode];
        int lastNode = -1;
        int currentNode = endMeshNode;
        double minDinstanceFromLandBoundary;
        double distanceFromLandBoundary;
        Point nodeOnLandBoundary;
        int currentNodeLandBoundaryNodeIndex;
        double currentNodeEdgeRatio;
        bool stopPathSearch;    //the path has been temporarily stopped(.true.) or not (.false.)
        int numConnectedNodes = 0;
        numRejectedNodesInPath = 0;
        numNodesInPath = 0;

        while (true)
        {
            stopPathSearch = true;

            if (m_meshNodesLandBoundarySegments[currentNode] >= 0)
            {
                // Multiple boundary segments: take the nearest
                int previousLandBoundarySegment = m_meshNodesLandBoundarySegments[currentNode];
                int previousStartMeshNode = m_segmentIndices[previousLandBoundarySegment][0];
                int previousEndMeshNode = m_segmentIndices[previousLandBoundarySegment][1];
                int previousEndLandBoundaryIndex;
                double previousMinDistance;

                successful = NearestLandBoundaryNode(mesh.m_projection, mesh.m_nodes[currentNode], previousStartMeshNode, previousEndMeshNode,
                    previousMinDistance, nodeOnLandBoundary, currentNodeLandBoundaryNodeIndex, currentNodeEdgeRatio);
                if (!successful)
                    return false;

                successful = NearestLandBoundaryNode(mesh.m_projection, mesh.m_nodes[currentNode], startLandBoundaryIndex, endLandBoundaryIndex,
                    distanceFromLandBoundary, nodeOnLandBoundary, currentNodeLandBoundaryNodeIndex, currentNodeEdgeRatio);
                if (!successful)
                    return false;

                const double minDinstanceFromLandBoundaryCurrentNode = m_nodesMinDistances[currentNode];

                if (distanceFromLandBoundary <= previousMinDistance &&
                    distanceFromLandBoundary < m_minDistanceFromLandFactor * minDinstanceFromLandBoundaryCurrentNode)
                {
                    stopPathSearch = false;
                }
            }
            else
            {
                if (m_nodesMinDistances[currentNode] == doubleMissingValue)
                {
                    successful = NearestLandBoundaryNode(mesh.m_projection, mesh.m_nodes[currentNode], 0, m_numNode,
                        minDinstanceFromLandBoundary, nodeOnLandBoundary, currentNodeLandBoundaryNodeIndex, currentNodeEdgeRatio);
                    if (!successful)
                        return false;
                    m_nodesMinDistances[currentNode] = minDinstanceFromLandBoundary;
                }
                else
                {
                    minDinstanceFromLandBoundary = m_nodesMinDistances[currentNode];
                }

                successful = NearestLandBoundaryNode(mesh.m_projection, mesh.m_nodes[currentNode], startLandBoundaryIndex, endLandBoundaryIndex,
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
                    m_meshNodesLandBoundarySegments[lastNode] = lastSegment;
                }
                numConnectedNodes = 0;
                numRejectedNodesInPath += 1;
            }
            else
            {
                lastSegment = m_meshNodesLandBoundarySegments[currentNode];
                lastNode = currentNode;

                numNodesInPath += 1;
                numConnectedNodes += 1;

                m_meshNodesLandBoundarySegments[lastNode] = landBoundarySegment;
            }

            if (stopPathSearch && numConnectedNodes == 1)
            {
                m_meshNodesLandBoundarySegments[lastSegment];
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
            m_meshNodesLandBoundarySegments[lastNode] = lastSegment;
        }

        return true;
    }

    /// masknodes
    /// mask the mesh nodes to be considered in the shortest path algorithm for the current segmentIndex
    /// is setting leftIndex, rightIndex, leftEdgeRatio, rightEdgeRatio 
    bool LandBoundaries::ComputeMask(const Mesh& mesh,
        const Polygons& polygon,
        int segmentIndex,
        bool meshBoundOnly,
        int& leftIndex,
        int& rightIndex,
        double& leftEdgeRatio,
        double& rightEdgeRatio,
        int& startLandBoundaryIndex,
        int& endLandBoundaryIndex)
    {
        std::fill(m_nodeMask.begin(), m_nodeMask.end(), doubleMissingValue);

        // check if any of the land boundary node is inside a mesh face
        bool nodeInFace = false;
        int crossedFaceIndex = -1;
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
                {
                    crossedFaceIndex = f;
                    break;
                }
            }

            if (nodeInFace)
                break;
        }

        // try to find a boundary cell that is crossed by the land boundary segment
        if (!nodeInFace)
        {
            crossedFaceIndex = -1;
            for (int e = 0; e < mesh.m_edges.size(); e++)
            {
                if (mesh.m_edgesNumFaces[e] != 1)
                    continue;

                bool isCellCrossed = IsFaceCrossedByLandBoundaries(mesh,
                    mesh.m_edgesFaces[e][0],
                    startLandBoundaryIndex,
                    endLandBoundaryIndex);

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

            MaskFaces(mesh,
                meshBoundOnly,
                landBoundaryFaces,
                startLandBoundaryIndex,
                endLandBoundaryIndex,
                leftIndex,
                rightIndex,
                leftEdgeRatio,
                rightEdgeRatio);

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

    /// maskcells
    /// mask the faces that are intersected by the land boundary
    bool LandBoundaries::MaskFaces(const Mesh& mesh,
        const bool& meshBoundOnly,
        std::vector<int>& landBoundaryFaces,
        int startNodeLandBoundaryIndex,
        int endNodeLandBoundaryindex,
        int& leftIndex,
        int& rightIndex,
        double& leftEdgeRatio,
        double& rightEdgeRatio)
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
                    if (mesh.m_edgesNumFaces[e] != 1 || mesh.m_edges[e].first < 0 || mesh.m_edges[e].second < 0)
                        continue;

                    bool isClose = false;
                    int landBoundaryNode = 0;
                    int otherFace = mesh.m_edgesFaces[e][0];
                    for (int ee = 0; ee < mesh.m_facesEdges[otherFace].size(); ee++)
                    {
                        int edge = mesh.m_facesEdges[otherFace][ee];
                        isClose = IsMeshEdgeCloseToLandBoundaries(mesh,
                            edge,
                            startNodeLandBoundaryIndex,
                            endNodeLandBoundaryindex,
                            meshBoundOnly,
                            leftIndex,
                            rightIndex,
                            leftEdgeRatio,
                            rightEdgeRatio,
                            landBoundaryNode);

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
                for (int e = 0; e < mesh.m_facesEdges[face].size(); e++)
                {
                    // is a boundary edge, continue
                    int currentEdge = mesh.m_facesEdges[face][e];
                    if (mesh.m_edgesNumFaces[currentEdge] <= 1)
                        continue;

                    int otherFace = mesh.m_edgesFaces[currentEdge][0] + mesh.m_edgesFaces[currentEdge][1] - face;

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
                            bool isClose = IsMeshEdgeCloseToLandBoundaries(mesh,
                                edgeOtherFace,
                                startNodeLandBoundaryIndex,
                                endNodeLandBoundaryindex,
                                meshBoundOnly,
                                leftIndex,
                                rightIndex,
                                leftEdgeRatio,
                                rightEdgeRatio,
                                landBoundaryNode);

                            if (isClose)
                            {
                                m_edgeMask[edgeOtherFace] = 1;
                                isFaceFound = 1;
                            }
                        }
                    }

                    m_faceMask[otherFace] = isFaceFound;

                    if (isFaceFound == 1)
                    {
                        m_numFacesMasked += 1;
                        numNextFaces += 1;
                        if (numNextFaces >= nextFaces.size())
                        {
                            nextFaces.resize(std::max(int(numNextFaces *1.2), 10));
                        }
                        nextFaces[numNextFaces - 1] = otherFace;
                    }
                }
            }
        }

        landBoundaryFaces.resize(0);

        if (numNextFaces > 0)
        {
            m_maskDepth += 1;

            MaskFaces(mesh,
                meshBoundOnly,
                nextFaces,
                startNodeLandBoundaryIndex,
                endNodeLandBoundaryindex,
                leftIndex,
                rightIndex,
                leftEdgeRatio,
                rightEdgeRatio);

            m_maskDepth -= 1;
        }

        return true;
    }

    /// linkcrossedbyland
    /// check if a mesh edge is close to a land boundary segment
    bool LandBoundaries::IsMeshEdgeCloseToLandBoundaries(const Mesh& mesh,
        int edgeIndex,
        int startNodeLandBoundaryIndex,
        int endNodeLandBoundaryIndex,
        bool meshBoundOnly,
        int& leftIndex,
        int& rightIndex,
        double& leftEdgeRatio,
        double& rightEdgeRatio,
        int& landBoundaryNode)
    {
        bool isClose = false;

        const int startNode = std::max(std::min(landBoundaryNode, endNodeLandBoundaryIndex - 1), startNodeLandBoundaryIndex);

        if (mesh.m_edges[edgeIndex].first < 0 || mesh.m_edges[edgeIndex].second < 0)
            return false;

        const auto firstMeshNode = mesh.m_nodes[mesh.m_edges[edgeIndex].first];
        const auto secondMeshNode = mesh.m_nodes[mesh.m_edges[edgeIndex].second];

        const double meshEdgeLength = Distance(firstMeshNode, secondMeshNode, mesh.m_projection);

        const double distanceFactor = meshBoundOnly ? m_closeToLandBoundaryFactor : m_closeWholeMeshFactor;
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

            if (landBoundaryLength > 0)
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
            while ((searchIter == 0 || currentNode < startNodeLandBoundaryIndex || currentNode > endNodeLandBoundaryIndex - 1) && searchIter < 3)
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
            if (landBoundaryNode < leftIndex)
            {
                leftIndex = landBoundaryNode;
                leftEdgeRatio = std::min(std::max(minimumRatio, 0.0), 1.0);
            }
            else if (landBoundaryNode == leftIndex)
            {
                leftEdgeRatio = std::min(std::max(minimumRatio, 0.0), leftEdgeRatio);
            }

            // maximum
            double maximumRatio = std::max(ratioFirstMeshNode, ratioSecondMeshNode);
            if (landBoundaryNode > rightIndex)
            {
                rightIndex = landBoundaryNode;
                rightEdgeRatio = std::min(std::max(maximumRatio, 0.0), 1.0);
            }
            else if (landBoundaryNode == rightIndex)
            {
                rightEdgeRatio = std::max(std::min(maximumRatio, 1.0), rightEdgeRatio);
            }
        }

        return isClose;
    }

    /// get_kstartend2
    /// Finds the start and end mesh node. These are the nodes that are
    /// on a edge close to the land boundary segment
    bool LandBoundaries::FindStartEndMeshNodes(const Mesh& mesh,
        const Polygons& polygons,
        int startLandBoundaryIndex,
        int endLandBoundaryIndex,
        int leftIndex,
        int rightIndex,
        double leftEdgeRatio,
        double rightEdgeRatio,
        int& startMeshNode,
        int& endMeshNode)
    {
        // compute the start and end point of the land boundary respectively
        int nextLeftIndex = std::min(leftIndex + 1, endLandBoundaryIndex);
        Point startPoint =
        {
            m_nodes[leftIndex].x + leftEdgeRatio * (m_nodes[nextLeftIndex].x - m_nodes[leftIndex].x),
            m_nodes[leftIndex].y + leftEdgeRatio * (m_nodes[nextLeftIndex].y - m_nodes[leftIndex].y)
        };

        int nextRightIndex = std::min(rightIndex + 1, endLandBoundaryIndex);
        Point endPoint =
        {
            m_nodes[rightIndex].x + rightEdgeRatio * (m_nodes[nextRightIndex].x - m_nodes[rightIndex].x),
            m_nodes[rightIndex].y + rightEdgeRatio * (m_nodes[nextRightIndex].y - m_nodes[rightIndex].y)
        };

        bool isStartPointInsideAPolygon = IsPointInPolygons(startPoint, polygons.m_nodes, polygons.m_numNodes);
        bool isEndPointInsideAPolygon = IsPointInPolygons(endPoint, polygons.m_nodes, polygons.m_numNodes);

        if (!isStartPointInsideAPolygon)
        {
            startPoint.x = m_nodes[leftIndex + 1].x;
            startPoint.y = m_nodes[leftIndex + 1].y;
        }
        if (!isEndPointInsideAPolygon)
        {
            endPoint.x = m_nodes[rightIndex].x;
            endPoint.y = m_nodes[rightIndex].y;
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
            startMeshNode = firstMeshNodeIndex;
        }
        else
        {
            startMeshNode = secondMeshNodeIndex;
        }

        firstMeshNodeIndex = mesh.m_edges[endEdge].first;
        secondMeshNodeIndex = mesh.m_edges[endEdge].second;
        firstDinstance = Distance(mesh.m_nodes[firstMeshNodeIndex], endPoint, mesh.m_projection);
        secondDinstance = Distance(mesh.m_nodes[secondMeshNodeIndex], endPoint, mesh.m_projection);

        if (firstDinstance <= secondDinstance)
        {
            endMeshNode = firstMeshNodeIndex;
        }
        else
        {
            endMeshNode = secondMeshNodeIndex;
        }

        return true;
    }



    /// Shortest_path
    /// connect mesh nodes starting from startMeshNode, using Dijkstra's shortest path algorithm
    /// the distance of each edge is the edge length multiplied by the distance from the land boundary
    bool LandBoundaries::ShortestPath(const Mesh& mesh, const Polygons& polygons, int landBoundarySegment,
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
            bool successful = NearestLandBoundaryNode(mesh.m_projection, currentNode, startLandBoundaryIndex, endLandBoundaryIndex,
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
                successful = NearestLandBoundaryNode(mesh.m_projection, neighbouringNode, startLandBoundaryIndex, endLandBoundaryIndex,
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
                successful = NearestLandBoundaryNode(mesh.m_projection, middlePoint, startLandBoundaryIndex, endLandBoundaryIndex,
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

    /// toland, compute the nearest point on the land boundary
    bool LandBoundaries::NearestLandBoundaryNode(const Projections& projection,
        const Point& node,
        int startLandBoundaryIndex,
        int endLandBoundaryIndex,
        double& minimumDistance,
        Point& pointOnLandBoundary,
        int& nearestLandBoundaryNodeIndex,
        double& edgeRatio)
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
            const double distanceFromLandBoundary = DistanceFromLine(node, m_nodes[n], m_nodes[n + 1], normalPoint, ratio, projection);

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


    /// cellcrossedbyland
    /// TODO: it could be moved to generic operations
    bool LandBoundaries::IsFaceCrossedByLandBoundaries(const Mesh& mesh, int face, int startLandBoundaryIndex, int endLandBoundaryIndex)
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
                double firstRatio;
                double secondRatio;
                bool areCrossing = AreLinesCrossing(firstMeshNode, secondMeshNode, firstNode, secondNode, adimensional, intersection, crossProduct, firstRatio, secondRatio, mesh.m_projection);

                if (areCrossing)
                {
                    return true;
                }
            }
        }
        return false;
    }


    /// snap_to_landboundary
    /// snap netnodes to land boundary segment
    bool LandBoundaries::SnapMeshToLandBoundaries(Mesh& mesh)
    {

        for (int n = 0; n < mesh.m_nodes.size(); n++)
        {
            if (mesh.m_nodesTypes[n] == 1 || mesh.m_nodesTypes[n] == 2 || mesh.m_nodesTypes[n] == 3)
            {
                int meshNodeToLandBoundarySegment = m_meshNodesLandBoundarySegments[n];
                if (meshNodeToLandBoundarySegment < 0)
                {
                    continue;
                }

                double minimumDistance;
                Point pointOnLandBoundary;
                int nearestLandBoundaryNodeIndex;
                double edgeRatio;

                bool successful = NearestLandBoundaryNode(mesh.m_projection,
                    mesh.m_nodes[n],
                    m_segmentIndices[meshNodeToLandBoundarySegment][0],
                    m_segmentIndices[meshNodeToLandBoundarySegment][1],
                    minimumDistance,
                    pointOnLandBoundary,
                    nearestLandBoundaryNodeIndex,
                    edgeRatio);

                if (!successful)
                {
                    return false;
                }

                mesh.m_nodes[n] = pointOnLandBoundary;
            }
        }
        return true;
    }
};
