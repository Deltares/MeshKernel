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
        
        LandBoundaries(const Projections& projection): m_numAllocatedNodes(0)
        {
            AllocateVector(m_numAllocatedNodes, m_nodes);
        }

        // admin_landboundary_segments
        bool Set(const std::vector<Point>& landBoundary)
        {
            int newSize = m_numAllocatedNodes + landBoundary.size();

            bool state = AllocateVector(newSize, m_nodes);

            m_numAllocatedNodes = newSize;

            for (int i= 0; i<landBoundary.size(); i++)
            {
                m_nodes[m_numNode] = landBoundary[i];
                m_numNode++;
            }

            return true;
        };


        bool FindNearestMeshBoundary(int meshLineOption)
        {

            return true;
        };

        // admin_landboundary_segments
        bool Administrate(Mesh& mesh, Polygons& polygons)
        {
            
            std::vector<int> landBoundaryMask(m_numAllocatedNodes - 1, 0);
            //mask the landboundary that is inside the selecting polygon
            for (int n = 0; n < m_numAllocatedNodes - 1; n++)
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
            const bool success = polygons.MeshBoundaryToPolygon(mesh, counterClockWise, setMeshState, meshBoundaryPolygon);
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
                    for (int nn = 0; nn < meshBoundaryPolygon.size() - 1; nn++)
                    {
                        Point firstMeshBoundaryNode = meshBoundaryPolygon[nn];
                        Point secondMeshBoundaryNode = meshBoundaryPolygon[nn + 1];

                        if(firstMeshBoundaryNode.x == doubleMissingValue || secondMeshBoundaryNode.x == doubleMissingValue)
                        {
                            continue;
                        }

                        const double edgeLength = distance(firstMeshBoundaryNode, secondMeshBoundaryNode, mesh.m_projection);
                        const double minDistance = m_closeToLandBoundary * edgeLength;

                        Point normalPoint;
                        double rlout;
                        const double distanceFirstMeshNode = distanceFromLine(firstMeshBoundaryNode, firstPoint, secondPoint, normalPoint, rlout, mesh.m_projection);
                        const double distanceSecondMeshNode = distanceFromLine(secondMeshBoundaryNode, firstPoint, secondPoint, normalPoint, rlout, mesh.m_projection);

                        if(distanceFirstMeshNode<= minDistance || distanceFirstMeshNode <= distanceSecondMeshNode)
                        {
                            landBoundaryIsClose = true;
                            break;
                        }
                    }

                    if(landBoundaryIsClose)
                    {
                        landBoundaryMask[n] = 1;
                    }
                }
            }

            // in m_segmentIndices store the starting and ending indexses of the land boundaries having
            // 1. the same landBoundaryMask
            // 2. are inside polygons
            // TODO: complete
            m_segmentIndices.resize(m_nodes.size(), std::vector<int>(2));
            auto begin = m_nodes.begin();
            auto end = m_nodes.end();
            while ( begin - end > 0)
            {
                auto found = std::find(begin, end,
                    [](auto const& node) { return node.x == doubleMissingValue; });

                int startInd = begin - m_nodes.begin();
                int endInd = found - m_nodes.begin();
                if( endInd > startInd && landBoundaryMask[startInd] != 0 )
                {
                    m_segmentIndices[m_numSegmentIndexses] = { startInd, endInd };
                    m_numSegmentIndexses++;
                }
                begin = found++;
            }

            return true;
        };

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

        int m_jleft;                                  // outer land boundary segments in projection
        int m_jright;

        double m_rLleft;                              // fractional location of the projected outer nodes(min and max) on the land boundary segment
        double m_rLright;        
        double m_DCLOSE_bound = 5.0;                  // close - to - landboundary tolerance, netbound only, measured in number of meshwidths
        double m_DCLOSE_whole = 1.0;                  // close - to - landboundary tolerance, whole network, measured in number of meshwidths
        double m_closeToLandBoundary = 1.0;                        // close - to - landboundary tolerance, measured in number of meshwidths
        bool m_Ladd_land = true;                      // add land boundary between land boundary segments that are close to each other
    };

}
