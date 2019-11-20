#pragma once

#include <vector>
#include "Mesh.hpp"
#include "Constants.cpp"
#include "Operations.cpp"
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

        // admin_landboundary_segments
        bool Set(const std::vector<Point>& landBoundary)
        {
            int newSize = m_numNode + landBoundary.size();

            bool state = AllocateVector(newSize, m_nodes, { doubleMissingValue,doubleMissingValue });

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
            m_meshBoundOnly = meshBoundOnly;
            if(m_meshBoundOnly)
            {
                m_relativeCloseDistance = 5;
            }
            else
            {
                m_relativeCloseDistance = 1;
            }
            
            m_nodeMask.resize(mesh.m_nodes.size(), intMissingValue);
            m_faceMask.resize(mesh.m_numFaces, intMissingValue);
            m_edgeMask.resize(mesh.m_edges.size(), intMissingValue);

            polygonCache.resize(maximumNumberOfNodesPerFace + 1);
            m_landMask = true;
            for (int n = 0; n < m_numSegmentIndexses; n++)
            {
                int numPaths;
                int numRejectedPaths;
                bool successful = MakePath(mesh, polygon, n, numPaths, numRejectedPaths);

            }

            //

            return true;
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
            const bool success = polygons.MeshBoundaryToPolygon(mesh, counterClockWise, setMeshState, meshBoundaryPolygon, numNodesBoundaryPolygons);
            if (!success)
            {
                return false;
            }

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

            // in m_segmentIndices store the starting and ending indexses of the land boundaries having
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


            int numSegmentIndexsesBeforeSplitting = m_numSegmentIndexses;
            // split the line segments into two to accommodate closed segments
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


        //make_path
        bool MakePath(const Mesh& mesh, const Polygons& polygons, int landBoundarySegment, int& numPaths, int& numRejectedPaths)
        {

            int startLandBoundaryIndex = m_segmentIndices[landBoundarySegment][0];
            int endLandBoundaryIndex = m_segmentIndices[landBoundarySegment][1];

            m_LeftIndex = endLandBoundaryIndex - 1;
            m_RightIndex = startLandBoundaryIndex;

            m_leftEdgeRatio = 1.0;
            m_rightEdgeRatio = 0.0;

            if (startLandBoundaryIndex < 0 || startLandBoundaryIndex >= m_numNode || startLandBoundaryIndex >= endLandBoundaryIndex)
                return true;

            bool success = ComputeMask(mesh, polygons, landBoundarySegment, startLandBoundaryIndex, endLandBoundaryIndex);

            if (!success)
                return false;

            int startMeshNode;
            int endMeshNode;
            success = FindStartEndNodes(mesh, polygons, startLandBoundaryIndex, endLandBoundaryIndex, startMeshNode, endMeshNode);
            if (!success)
                return false;

            if (startMeshNode < 0 || endMeshNode < 0 || startMeshNode == endMeshNode)
                return false;

            std::vector<int> connectedEdges;
            success = ShortestPath(mesh, polygons,
                startLandBoundaryIndex, endLandBoundaryIndex, startMeshNode, connectedEdges);
            if (!success)
                return false;



            return true;
        }


        // Shortest_path
        // connect all edges closest to the land boundary using Dijkstra's shortest path algorithm
        // the distance of each edge is the actual edge length multiplied by the distance from the land boundary
        bool ShortestPath(const Mesh& mesh, const Polygons& polygons, 
            int startLandBoundaryIndex, int endLandBoundaryIndex, int startMeshNode, std::vector<int>& connectedEdges)
        {
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
                bool success = NearestPointOnLandBoundary(mesh, currentNode, startLandBoundaryIndex, endLandBoundaryIndex,
                    currentNodeDistance, currentNodeOnLandBoundary, currentNodeLandBoundaryNodeIndex, currentNodeEdgeRatio);
                
                if (!success || currentNodeLandBoundaryNodeIndex < 0)
                    return false;

                for (int e = 0; e < mesh.m_nodesEdges[currentNodeIndex].size(); e++) 
                {
                    int edgeIndex = mesh.m_nodesEdges[currentNodeIndex][e];
                    auto edge = mesh.m_edges[edgeIndex];

                    if (edge.first < 0 || edge.second < 0) 
                    {
                        continue;
                    }

                    int neighbouringNodeIndex = edge.first + edge.second - currentNodeIndex;
                    Point neighbouringNode = mesh.m_nodes[neighbouringNodeIndex];

                    Point neighbouringNodeOnLandBoundary;
                    int neighbouringNodeLandBoundaryNodeIndex;
                    double neighbouringNodeEdgeRatio;
                    double neighbouringNodeDistance;
                    success = NearestPointOnLandBoundary(mesh, neighbouringNode, startLandBoundaryIndex, endLandBoundaryIndex,
                        neighbouringNodeDistance, neighbouringNodeOnLandBoundary, neighbouringNodeLandBoundaryNodeIndex, neighbouringNodeEdgeRatio);
                    if (!success) 
                        return false;

                    Point middlePoint
                    {
                        (currentNode.x + neighbouringNode.x) / 2.0,
                        (currentNode.y + neighbouringNode.y) / 2.0
                    };

                    Point middlePointOnLandBoundary;
                    int middlePointLandBoundaryNodeIndex;
                    double middlePointEdgeRatio;
                    double middlePointDistance;
                    success = NearestPointOnLandBoundary(mesh, middlePoint, startLandBoundaryIndex, endLandBoundaryIndex,
                        middlePointDistance, middlePointOnLandBoundary, middlePointLandBoundaryNodeIndex, middlePointEdgeRatio);
                    if (!success) 
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
                    if (m_meshBoundOnly && mesh.m_edgesNumFaces[edgeIndex] != 1)
                    {
                        maximumDistance = 1e6 * maximumDistance;
                    }

                    double edgeLength = Distance(currentNode, neighbouringNode, mesh.m_projection);
                    double correctedDistance = nodeDistances[currentNodeIndex] + edgeLength * maximumDistance;

                    if (correctedDistance < nodeDistances[neighbouringNodeIndex]) 
                    {
                        nodeDistances[neighbouringNodeIndex] = correctedDistance;
                        connectedEdges[neighbouringNodeIndex] = edgeIndex;
                    }
                }

                // for the next node find the node index with the minimum distance 
                currentNodeIndex = std::min(nodeDistances.begin(), nodeDistances.end()) - nodeDistances.begin();
                if (currentNodeIndex < 0 ||
                    currentNodeIndex >= mesh.m_nodes.size() ||
                    nodeDistances[currentNodeIndex] == std::numeric_limits<double>::max() ||
                    m_nodeMask[currentNodeIndex] < 0) 
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
                auto firstNode = m_nodes[n];
                auto secondNode = m_nodes[n + 1];

                if (firstNode.x == doubleMissingValue || firstNode.x == doubleMissingValue)
                    continue;
                
                Point normalPoint;
                double ratio;
                const double distanceFromLandBoundary = DistanceFromLine(node, firstNode, secondNode, normalPoint, ratio, mesh.m_projection);

                if (distanceFromLandBoundary > 0.0 && distanceFromLandBoundary<minimumDistance)
                {
                    minimumDistance = distanceFromLandBoundary;
                    pointOnLandBoundary = normalPoint;
                    nearestLandBoundaryNodeIndex = n;
                    edgeRatio = ratio;
                }
            }
            return true;
        }



        //get_kstartend2
        bool FindStartEndNodes(const Mesh& mesh, const Polygons& polygons, int startLandBoundaryIndex, int endLandBoundaryIndex, int& startLandBoundaryNode, int& endLandBoundaryNode)
        {
            //compute the start and end point of the land boundary respectively
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

            // Get the links that are closest the land boundary start and end respectively
            Point normalPoint;
            double ratio;

            double minDistStart = std::numeric_limits<double>::max();
            double minDistEnd = std::numeric_limits<double>::max();
            int startEdge = -1;
            int endEdge = -1;

            for (int e = 0; e < mesh.m_edges.size(); e++)
            {

                int firstMeshNodeIndex = mesh.m_edges[e].first;
                int secondMeshNodeIndex = mesh.m_edges[e].second;

                if (firstMeshNodeIndex < 0 || secondMeshNodeIndex < 0)
                    continue;

                if (m_nodeMask[firstMeshNodeIndex] < 0 || m_nodeMask[secondMeshNodeIndex] < 0) 
                    continue;

                const double distanceFromFirstMeshNode = DistanceFromLine(startPoint, mesh.m_nodes[firstMeshNodeIndex], mesh.m_nodes[secondMeshNodeIndex], normalPoint, ratio, mesh.m_projection);
                const double distanceFromSecondMeshNode = DistanceFromLine(endPoint, mesh.m_nodes[firstMeshNodeIndex], mesh.m_nodes[secondMeshNodeIndex], normalPoint, ratio, mesh.m_projection);

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


        //masknodes
        bool ComputeMask(const Mesh& mesh, const Polygons& polygon,int segmentIndex, int& startLandBoundaryIndex, int& endLandBoundaryIndex)
        {
            std::fill(m_nodeMask.begin(), m_nodeMask.end(), doubleMissingValue);

            // first node of land boundary segment in mesh
            bool nodeInFace = false;
            for (int i = startLandBoundaryIndex; i < endLandBoundaryIndex; i++)
            {
                auto node = m_nodes[i];

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

                    nodeInFace = IsPointInPolygon(node, polygonCache, numFaceNodes);
                    if (nodeInFace)
                        break;
                }

                if (nodeInFace)
                    break;
            }

            // try to find a boundary cell that is crossed by the land boundary segment
            int crossedFace = -1;
            if (!nodeInFace)
            {
                for (int e = 0; e < mesh.m_edges.size(); e++)
                {
                    if (mesh.m_edgesNumFaces[e] != 1)
                        continue;

                    int face = mesh.m_edgesFaces[e][0];

                    bool isCellCrossed = IsFaceCrossedByLandBoundaries(mesh, face, startLandBoundaryIndex, endLandBoundaryIndex);

                    if (isCellCrossed)
                    {
                        crossedFace = face;
                        break;
                    }
                }
            }

            if (m_landMask)
            {
                std::fill(m_faceMask.begin(), m_faceMask.end(), doubleMissingValue);
                std::fill(m_edgeMask.begin(), m_edgeMask.end(), doubleMissingValue);
                //m_faceMask assumes crossedFace has already been done.
                if (crossedFace >= 0 && crossedFace < mesh.m_numFaces)
                {
                    m_faceMask[crossedFace] = 1;
                }

                m_numFacesMasked = 0;
                m_maskDepth = 0;
                std::vector<int> landBoundaryFaces{ crossedFace };
                MaskFaces(mesh, landBoundaryFaces, startLandBoundaryIndex, endLandBoundaryIndex);

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
                if(m_nodeMask[n]> 0)
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

        //linkcrossedbyland
        bool IsEdgeCloseToLandBoundaries(const Mesh& mesh, int edgeIndex, int& landBoundaryNode, int startNodeLandBoundary, int endNodeLandBoundary)
        {
            bool isClose = false;

            const int startNode = std::max(std::min(landBoundaryNode, endNodeLandBoundary - 1), landBoundaryNode);

            const int firstMeshNodeIndex = mesh.m_edges[edgeIndex].first;
            const int secondMeshNodeIndex = mesh.m_edges[edgeIndex].second;

            if (firstMeshNodeIndex < 0 || secondMeshNodeIndex < 0)
                return isClose;

            const auto firstMeshNode = mesh.m_nodes[firstMeshNodeIndex];
            const auto secondMeshNode = mesh.m_nodes[secondMeshNodeIndex];

            const double meshEdgeLength = Distance(firstMeshNode, secondMeshNode, mesh.m_projection);
            const double closeDistance = meshEdgeLength * m_relativeCloseDistance;

            //TODO: this is a slightly different implementation 
            double ratioFirstMeshNode;
            double ratioSecondMeshNode;
            for (int n = startNode; n < endNodeLandBoundary; n++)
            {

                auto firstNode = m_nodes[n];
                auto secondNode = m_nodes[n + 1];

                if (firstNode.x == doubleMissingValue || firstNode.x == doubleMissingValue)
                    continue;

                const double landBoundaryLength = Distance(firstNode, secondNode, mesh.m_projection);

                if (landBoundaryLength <= std::numeric_limits<double>::min())
                    continue;

                Point normalPoint;
                ratioFirstMeshNode = 0.0;
                ratioSecondMeshNode = 1.0;

                const double distanceFromLandBoundaryFirstMeshNode = DistanceFromLine(firstMeshNode, firstNode, secondNode, normalPoint, ratioFirstMeshNode, mesh.m_projection);

                if (distanceFromLandBoundaryFirstMeshNode < closeDistance)
                {
                    isClose = true;
                    landBoundaryNode = n;
                    // the projection is within the land boundary segment
                    if (ratioFirstMeshNode >= 0.0 && ratioFirstMeshNode <= 1.0)
                    {
                        break;
                    }
                }
                else
                {
                    double distanceFromLandBoundarySecondMeshNode = DistanceFromLine(secondMeshNode, firstNode, secondNode, normalPoint, ratioSecondMeshNode, mesh.m_projection);

                    if (distanceFromLandBoundarySecondMeshNode < closeDistance)
                    {
                        isClose = true;
                        landBoundaryNode = n;
                        // the projection is within the land boundary segment
                        if (ratioSecondMeshNode >= 0.0 && ratioSecondMeshNode <= 1.0)
                        {
                            break;
                        }
                    }
                }
            }


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



        bool MaskFaces(const Mesh& mesh, std::vector<int>& landBoundaryFaces, int& startNodeLandBoundary, int &endNodeLandBoundary)
        {
            int numNextFaces = 0;
            std::vector<int> nextFaces(landBoundaryFaces.size(), intMissingValue);
            for (int f = 0; f < landBoundaryFaces.size(); f++)
            {
                // no start face specified: mask boundary faces only
                // these are the faces that are close (up to a certain tolerance) by a land boundary
                int face = landBoundaryFaces[f];
                if (face < 0)
                {
                    int endNode = 0;
                    for (int e = 0; e < mesh.m_edges.size(); e++)
                    {
                        if (mesh.m_edgesNumFaces[e] != 1)
                            continue;

                        int otherFace = mesh.m_edgesFaces[e][0];

                        int firstMeshNodeIndex = mesh.m_edges[e].first;
                        int secondMeshNodeIndex = mesh.m_edges[e].second;

                        if (firstMeshNodeIndex < 0 || secondMeshNodeIndex < 0)
                            continue;

                        bool isClose = false;
                        for (int ee = 0; ee < mesh.m_facesEdges[otherFace].size(); ee++)
                        {
                            int landBoundaryNode = 0;
                            int edge = mesh.m_facesEdges[otherFace][ee];

                            isClose = IsEdgeCloseToLandBoundaries(mesh, edge, landBoundaryNode, startNodeLandBoundary, endNodeLandBoundary);

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
                    const int numEdges = mesh.m_facesEdges.size();

                    if (numEdges < 3)
                        continue;

                    bool isClose = false;
                    int isFaceFound = 0;
                    for (int e = 0; e < numEdges; e++)
                    {
                        int edge = mesh.m_facesEdges[face][e];
                        // If boundary face : no further checking in that direction.
                        if (mesh.m_edgesNumFaces[edge] <= 1)
                            continue;

                        int otherFace = mesh.m_edgesFaces[e][0] + mesh.m_edgesFaces[e][1] - face;
                        if (m_faceMask[otherFace] != intMissingValue)
                            continue;


                        for (int ee = 0; ee < mesh.m_facesEdges[otherFace].size(); ee++)
                        {
                            int edge = mesh.m_facesEdges[otherFace][ee];

                            if (m_edgeMask[edge] == 1)
                            {
                                // previously visited crossed link
                                isClose = true;
                                isFaceFound = 1;
                            }
                            else if (m_edgeMask[edge] == 0)
                            {
                                // previously visited uncrossed link - check next (kothercell) link
                                continue;
                            }
                            else
                            {
                                int endNode = 0;
                                // linkmask is IMISS, i.e.the link is unvisited, unvisited links
                                m_edgeMask[edge] = 0;

                                int landBoundaryNode = 0;
                                //TODO: to check
                                isClose = IsEdgeCloseToLandBoundaries(mesh, edge, landBoundaryNode, startNodeLandBoundary, endNodeLandBoundary);

                                if (isClose)
                                {
                                    m_edgeMask[edge] = 1;
                                    isFaceFound = 1;
                                }
                            }
                        }

                        // the rest of the loop
                        m_faceMask[otherFace] = isFaceFound;

                        if(isFaceFound==1)
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

            if(numNextFaces> 0 )
            {
                m_maskDepth += 1;
                MaskFaces(mesh, nextFaces, startNodeLandBoundary, endNodeLandBoundary);
                m_maskDepth += 1;
            }

            return true;
        }

        std::vector<Point> m_nodes;                   // XLAN, YLAN, ZLAN
        std::vector<int>   m_nclan;                   // NCLAN

        // number of land boundary segments
        int m_numSegments;                            // Nlanseg
        int m_numNode;                                // actual MXLAN
        int m_numAllocatedNodes;                      // MAXLAN
        int m_numNodesLoc;                            // MXLAN_loc

        std::vector<std::vector<int>> m_segmentIndices;  // lanseg_startend
        int m_numSegmentIndexses = 0;
        std::vector<std::vector<double>> m_nodesLand;    // !node to land boundary segment mapping
        std::vector<int> m_nodeMask;                  // masking the net nodes
        std::vector<int> m_faceMask;                  // masking faces
        std::vector<int> m_edgeMask;                  // masking edges
        std::vector<Point> polygonCache;                 // array of points (e.g. points of a face)
        bool m_landMask = true;
        int m_numFacesMasked = 0;
        int m_maskDepth = 0;

        int m_LeftIndex;                                  // outer land boundary segments in projection
        int m_RightIndex;
        double m_leftEdgeRatio;                              // fractional location of the projected outer nodes(min and max) on the land boundary segment
        double m_rightEdgeRatio;

        double m_closeToLandBoundary = 5.0;           // close - to - landboundary tolerance, measured in number of meshwidths

        double m_relativeCloseDistance = 1.0;         // DCLOSE, close - to - landboundary tolerance, measured in number of meshwidths

        bool m_Ladd_land = true;                      // add land boundary between land boundary segments that are close to each other

        bool m_meshBoundOnly = false;
    };

}
