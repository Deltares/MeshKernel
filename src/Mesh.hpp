#pragma once

#include <vector>
#include "Entities.hpp"

namespace GridGeom 
{
    class Mesh
    {

    public:

        Mesh(){}

        bool Set(const std::vector<Edge>& edges, const std::vector<Point>& nodes, Projections projection);
        
        bool SetFlatCopies();
        
        bool DeleteFlatCopies();

        void FacesAreasAndMassCenters();

        void FaceCircumcenters(const double& weightCircumCenter);

        void FindFaces();

        std::vector<Edge>  m_edges;                                 // KN
        std::vector<Point> m_nodes;                                 // KN
        std::vector<std::vector<std::size_t>> m_nodesEdges;         // NOD
        std::vector<std::size_t> m_nodesNumEdges;                   // NMK

        //edges
        std::vector<std::size_t> m_edgesNumFaces;                   // LNN
        std::vector<std::vector<int>> m_edgesFaces;                 // LNE

        // faces
        std::vector<std::vector<std::size_t>> m_facesNodes;         // netcell%Nod, the nodes composing the faces, in ccw order
        std::vector<std::vector<std::size_t>> m_facesEdges;         // netcell%lin
        std::vector<Point>   m_facesCircumcenters;                  // xz  the face circumcenter
        std::vector<Point>   m_facesMassCenters;                    // xzw the faces canters of mass

        std::size_t m_numFaces;                                     // NUMP
        std::vector<double> m_faceArea;                             // Face area
        
        //Used for internal state
        std::vector<double> m_nodex;
        std::vector<double> m_nodey;
        std::vector<double> m_nodez;
        std::vector<int>    m_edgeNodes;

        Projections m_projection;

    private:

        // Set node admin
        void NodeAdministration();

        // Sort_links_ccw
        void SortEdgesInCounterClockWiseOrder();

        // find cells
        void FindFaces(const int& numEdges);

        double m_dcenterinside = 1.0;

    };
}