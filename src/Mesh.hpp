#pragma once

#include <vector>
#include "Entities.hpp"
#include "IOperations.hpp"

namespace GridGeom 
{
    class Mesh
    {
    public:
        Mesh()
        {
        }

        explicit Mesh(IOperations* operations)
            : m_operations(operations)
        {
        }

        bool setMesh(const std::vector<Edge>& edges, const std::vector<Point>& nodes);
        bool setState();
        bool deleteState();

        double m_dcenterinside = 1.0;

        std::vector<Edge>  m_edges;                                 // KN
        std::vector<Point> m_nodes;                                 // KN
        std::vector<std::vector<size_t>> m_nodesEdges;              // NOD
        std::vector<size_t> m_nodesNumEdges;                        // NMK

        //edges
        std::vector<size_t> m_edgesNumFaces;                        // LNN
        std::vector<std::vector<int>> m_edgesFaces;                 // LNE

        // faces
        std::vector<std::vector<size_t>> m_facesNodes;              // netcell%Nod, the nodes composing the faces, in ccw order
        std::vector<std::vector<size_t>> m_facesEdges;              // netcell%lin
        std::vector<Point>   m_facesCircumcenters;                  // xz  the face circumcenter
        std::vector<Point>   m_facesMassCenters;                    // xzw the faces canters of mass

        size_t m_numFaces;                                          // NUMP
        std::vector<double> m_faceArea;                             // Face area

        IOperations* m_operations;                                  // Run-time selection of the operations to perform
        
        //Used for internal state
        std::vector<double> m_nodex;
        std::vector<double> m_nodey;
        std::vector<double> m_nodez;
        std::vector<int>    m_edgeNodes;

        void facesAreasAndMassCenters();

        void faceCircumcenters(const double& weightCircumCenter);

    private:

        // Set node admin
        void NodeAdministration();

        void SortEdgesInCounterClockWiseOrder();

        void findFaces(const int& numEdges);

    };
}