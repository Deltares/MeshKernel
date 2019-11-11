#pragma once

#include <utility>
#include <vector>
#include <algorithm>
#include "Mesh.hpp"
#include "Constants.cpp"
#include "Operations.cpp"

namespace GridGeom
{
    class LandBoundaries
    {
        // admin_landboundary_segments
        bool Set(const std::vector<double>& nodex, const std::vector<double>& nodey, const std::vector<double>& nodez)
        {
            PreAllocateLandBoundaries(nodex.size());

            for (int i=0; i<nodex.size(); i++)
            {
                m_nodex[i] = nodex[i];
                m_nodey[i] = nodey[i];
                m_nodez[i] = nodez[i];
            }

            return true;
        };


        bool FindNearestMeshBoundary(int meshLineOption)
        {

            return true;
        };

        // admin_landboundary_segments
        bool Administrate(std::vector<Point>& polygon)
        {
            m_indices.resize(2, std::vector<int>(1, 0));
            std::vector<int> inPolygon(m_numReserved-1, 0); //lanmask

            for (int n = 0; n < m_numReserved - 1; n++)
            {
                if (m_nodes[n].x != doubleMissingValue && m_nodes[n + 1].x != doubleMissingValue)
                {
                    bool firstPointInPolygon = pointInPolygon(m_nodes[n], polygon, polygon.size());
                    bool secondPointInPolygon = pointInPolygon(m_nodes[n + 1], polygon, polygon.size());

                    if (firstPointInPolygon || secondPointInPolygon)
                    {
                        inPolygon[n] = -1;
                    }
                }
            }

            // network boundary to polygon


            return true;
        };


        bool PreAllocateLandBoundaries(int size)
        {
            if (size < m_numReserved)
            {
                return true;
            }

            m_numReserved = std::max(50000, static_cast<int>(1.2)*size);

            m_nodex.resize(m_numReserved, doubleMissingValue);
            m_nodey.resize(m_numReserved, doubleMissingValue);
            m_nodez.resize(m_numReserved, doubleMissingValue);
            m_nclan.resize(m_numReserved, doubleMissingValue);
            return true;
        }

        std::vector<Point> m_nodes;                   // XLAN, YLAN, ZLAN
        std::vector<int>   m_nclan;                  // NCLAN

        // number of land boundary segments
        int m_numSegments;                            // Nlanseg
        int m_numNode;                                // actual MXLAN
        int m_numReserved;                            // MAXLAN
        int m_numNodesLoc;                            // MXLAN_loc

        std::vector<std::vector<int>> m_indices;   // lanseg_startend
        std::vector<std::vector<double>> m_nodesLand; // !node to land boundary segment mapping

        int m_jleft;                                  // outer land boundary segments in projection
        int m_jright;

        double m_rLleft;                              // fractional location of the projected outer nodes(min and max) on the land boundary segment
        double m_rLright;        
        double m_DCLOSE_bound = 5.0;                  // close - to - landboundary tolerance, netbound only, measured in number of meshwidths
        double m_DCLOSE_whole = 1.0;                  // close - to - landboundary tolerance, whole network, measured in number of meshwidths
        double m_DCLOSE = 1.0;                        // close - to - landboundary tolerance, measured in number of meshwidths
        bool m_Ladd_land = true;                      // add land boundary between land boundary segments that are close to each other

    };

}