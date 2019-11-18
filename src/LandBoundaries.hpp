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


        bool FindNearestMeshBoundary(const Mesh& mesh, const Polygons& polygon, int meshLineOption)
        {
            m_nodeMask.resize(mesh.m_nodes.size(), intMissingValue);
            m_faceMask.resize(mesh.m_numFaces, intMissingValue);
            m_edgeMask.resize(mesh.m_edges.size(), intMissingValue);

            polygonCache.resize(maximumNumberOfNodesPerFace);
            m_landMask = true;
            for (int n = 0; n < m_numSegmentIndexses; n++)
            {
                int numPaths;
                int numRejectedPaths;
                bool successful = MakePath(mesh, polygon, n, numPaths, numRejectedPaths);



            }

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
                    bool firstPointInPolygon = pointInPolygon(m_nodes[n], polygons.m_nodes, polygons.m_numNodes);
                    bool secondPointInPolygon = pointInPolygon(m_nodes[n + 1], polygons.m_nodes, polygons.m_numNodes);

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
                    const double landBoundaryLength = distance(firstPoint, secondPoint, mesh.m_projection);

                    bool landBoundaryIsClose = false;
                    for (int nn = 0; nn < numNodesBoundaryPolygons - 1; nn++)
                    {
                        Point firstMeshBoundaryNode = meshBoundaryPolygon[nn];
                        Point secondMeshBoundaryNode = meshBoundaryPolygon[nn + 1];

                        if (firstMeshBoundaryNode.x == doubleMissingValue || secondMeshBoundaryNode.x == doubleMissingValue)
                        {
                            continue;
                        }

                        const double edgeLength = distance(firstMeshBoundaryNode, secondMeshBoundaryNode, mesh.m_projection);
                        const double minDistance = m_closeToLandBoundary * edgeLength;

                        Point normalPoint;
                        double rlout;
                        const double distanceFirstMeshNode = distanceFromLine(firstMeshBoundaryNode, firstPoint, secondPoint, normalPoint, rlout, mesh.m_projection);
                        const double distanceSecondMeshNode = distanceFromLine(secondMeshBoundaryNode, firstPoint, secondPoint, normalPoint, rlout, mesh.m_projection);

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


            int numSegmentIndexses = m_numSegmentIndexses;
            // split the line segments into two to accommodate closed segments
            for (int i = 0; i < numSegmentIndexses; i++)
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
        bool MakePath(const Mesh& mesh, const Polygons& polygon, int landBoundarySegment, int& numPaths, int& numRejectedPaths)
        {

            int startLandBoundaryIndex = m_segmentIndices[landBoundarySegment][0];
            int endLandBoundaryIndex = m_segmentIndices[landBoundarySegment][1];

            m_outerLeft = endLandBoundaryIndex - 1;
            m_outerRight = startLandBoundaryIndex;

            m_leftRatio = 1.0;
            m_rightRatio = 0.0;

            if (startLandBoundaryIndex < 0 || startLandBoundaryIndex >= m_numNode || startLandBoundaryIndex >= endLandBoundaryIndex)
                return true;

            bool state = ComputeMask(mesh, polygon, landBoundarySegment, startLandBoundaryIndex, endLandBoundaryIndex);

            //will set jleft, jright, rLleft and rLright
            //!write(6, *) numseg, jstart, jend, jleft, jright, rLleft, rLright

            return true;
        }

        //masknodes
        bool ComputeMask(const Mesh& mesh, const Polygons& polygon,int segmentIndex, int& jstart, int& jend)
        {
            std::fill(m_nodeMask.begin(), m_nodeMask.end(), doubleMissingValue);
            //find the start cell for node masking
            int kstart = 0;
            int in = -1;
            // first node of land boundary segment in mesh
            int jstart1 = jstart;
            bool nodeInFace = false;
            for (int l = jstart; l < jend; l++)
            {
                auto node = m_nodes[l];

                for (int f = 0; f < mesh.m_numFaces; f++)
                {
                    auto numFaceNodes = mesh.m_facesNodes[f].size();

                    if (numFaceNodes == 0)
                        continue;

                    for (int n = 0; n < numFaceNodes; n++)
                    {
                        polygonCache[n] = mesh.m_nodes[mesh.m_facesNodes[f][n]];
                    }

                    nodeInFace = pointInPolygon(node, polygonCache, numFaceNodes);
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

                    bool isCellCrossed = isCellCrossedByLandBoundaries(mesh, face, jstart, jend);

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
                maskFaces(mesh, landBoundaryFaces, jstart, jend);

                // Mask all nodes of the masked faces
                for (int f = 0; f < mesh.m_numFaces; f++)
                {
                    if (m_faceMask[f] == 1)
                    {
                        for (int n = 0; n < mesh.m_facesNodes[f].size(); n++)
                        {
                            m_faceMask[mesh.m_facesNodes[f][n]] = segmentIndex;
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
                    bool inPolygon = pointInPolygon(mesh.m_nodes[n], polygon.m_nodes, polygon.m_numNodes);
                    if (!inPolygon)
                    {
                        m_nodeMask[n] = 0;
                    }
                }
            }

            return true;
        }

        // cellcrossedbyland
        bool isCellCrossedByLandBoundaries(const Mesh& mesh, int face, int jstart, int jend)
        {
            bool isCrossed = false;
            bool areCrossing = false;

            for (int e = 0; e < mesh.m_facesEdges[face].size(); e++)
            {
                auto edge = mesh.m_facesEdges[face][e];

                for (int s = jstart; s < jend - 1; s++)
                {
                    auto firstMeshNode = mesh.m_nodes[mesh.m_edges[edge].first];
                    auto secondMeshNode = mesh.m_nodes[mesh.m_edges[edge].second];

                    auto firstNode = m_nodes[s];
                    auto secondNode = m_nodes[s + 1];
                    bool adimensional = false;
                    Point intersection;
                    double crossProduct;

                    bool areCrossing = linesCrossing(firstMeshNode, secondMeshNode, firstNode, secondNode, adimensional, intersection, crossProduct, mesh.m_projection);

                    if (areCrossing)
                    {
                        return true;
                    }
                }
            }
            return false;
        }

        //linkcrossedbyland
        bool isEdgeCloseToLandBoundaries(const Mesh& mesh, int edgeIndex, int endNode, int& landBoundaryNode, int& startNodeLandBoundary, int& endNodeLandBoundary)
        {
            bool isClose = false;

            int startNode = std::max(std::min(landBoundaryNode, endNode), landBoundaryNode);

            int firstMeshNodeIndex = mesh.m_edges[edgeIndex].first;
            int secondMeshNodeIndex = mesh.m_edges[edgeIndex].second;

            if (firstMeshNodeIndex < 0 || secondMeshNodeIndex < 0)
                return isClose;

            auto firstMeshNode = mesh.m_nodes[firstMeshNodeIndex];
            auto secondMeshNode = mesh.m_nodes[secondMeshNodeIndex];

            double edgeLength = distance(firstMeshNode, secondMeshNode, mesh.m_projection);
            double closeDistance = edgeLength * m_relativeCloseDistance;

            //TODO: this is a slightly different implementation 
            double ratioFirstMeshNode;
            double ratioSecondMeshNode;
            for (int n = startNode; n < endNode - 1; n++)
            {

                auto firstNode = m_nodes[n];
                auto secondNode = m_nodes[n + 1];

                if (firstNode.x == doubleMissingValue || firstNode.x == doubleMissingValue)
                    continue;

                double edgeLength = distance(firstNode, secondNode, mesh.m_projection);

                if (edgeLength <= std::numeric_limits<double>::min())
                    continue;

                Point normalPoint;
                ratioFirstMeshNode = 0.0;
                ratioSecondMeshNode = 1.0;

                const double distanceFromLandBoundaryFirstMeshNode = distanceFromLine(firstMeshNode, firstNode, secondNode, normalPoint, ratioFirstMeshNode, mesh.m_projection);

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
                    double distanceFromLandBoundarySecondMeshNode = distanceFromLine(secondMeshNode, firstNode, secondNode, normalPoint, ratioSecondMeshNode, mesh.m_projection);

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
                if (landBoundaryNode < startNodeLandBoundary)
                {
                    startNodeLandBoundary = landBoundaryNode;
                    m_leftRatio = std::min(std::max(minimumRatio, 0.0), 1.0);
                }
                else if (landBoundaryNode == startNodeLandBoundary)
                {
                    m_leftRatio = std::min(std::max(minimumRatio, 0.0), m_leftRatio);
                }

                // maximum
                double maximumRatio = std::max(ratioFirstMeshNode, ratioSecondMeshNode);
                if (landBoundaryNode > endNodeLandBoundary)
                {
                    endNodeLandBoundary = landBoundaryNode;
                    m_rightRatio = std::min(std::max(maximumRatio, 0.0), 1.0);
                }
                else if (landBoundaryNode == endNodeLandBoundary)
                {
                    m_rightRatio = std::max(std::min(minimumRatio, 1.0), m_rightRatio);
                }
            }
            return isClose;
        }



        bool maskFaces(const Mesh& mesh, std::vector<int>& landBoundaryFaces, int& startNodeLandBoundary, int &endNodeLandBoundary)
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
                            //TODO: to check
                            bool isClose =
                                isEdgeCloseToLandBoundaries(mesh, edge, endNode, landBoundaryNode, startNodeLandBoundary, endNodeLandBoundary);

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
                                isClose = isEdgeCloseToLandBoundaries(mesh, edge, endNode, landBoundaryNode, startNodeLandBoundary, endNodeLandBoundary);

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
                maskFaces(mesh, nextFaces, startNodeLandBoundary, endNodeLandBoundary);
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

        int m_outerLeft;                                  // outer land boundary segments in projection
        int m_outerRight;
        double m_leftRatio;                              // fractional location of the projected outer nodes(min and max) on the land boundary segment
        double m_rightRatio;

        double m_DCLOSE_bound = 5.0;                  // close - to - landboundary tolerance, netbound only, measured in number of meshwidths
        double m_DCLOSE_whole = 1.0;                  // close - to - landboundary tolerance, whole network, measured in number of meshwidths
        double m_closeToLandBoundary = 5.0;           // close - to - landboundary tolerance, measured in number of meshwidths

        double m_relativeCloseDistance = 1.0;         //close - to - landboundary tolerance, measured in number of meshwidths

        bool m_Ladd_land = true;                      // add land boundary between land boundary segments that are close to each other
    };

}
