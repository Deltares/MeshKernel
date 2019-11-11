#pragma once

#include <utility>
#include <vector>
#include "Mesh.hpp"

namespace GridGeom
{
    class LandBoundaries
    {
        //admin_landboundary_segments
        bool Set(const Mesh& mesh);
        bool PreAllocateLandBoundaries(int size);

        std::vector<double> m_nodex;
        std::vector<double> m_nodey;
        std::vector<double> m_nodez;

        // number of land boundary segments
        int m_numSegments;                            // Nlanseg
        int m_mxlan;                                  // actual MXLAN
        int m_maxlan;                                 // MAXLAN
        int m_maxlanloc;                              // MXLAN_loc

        std::vector<std::vector<double>> m_indices;   // lanseg_startend
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