#ifndef MESH_HPP
#define MESH_HPP
#include <vector>
#include "Entities.hpp"

template<typename Point>
class Mesh
{
public:

    Mesh(const std::vector<GridGeom::Edge>& edges, const std::vector<Point>& nodes);

    const double m_dcenterinside = 1.0;

    std::vector<GridGeom::Edge> m_edges;                        //KN
    std::vector<Point> m_nodes;                                 //KN
    std::vector<std::vector<size_t>> m_nodesEdges;              //NOD
    std::vector<size_t> m_nodesNumEdges;                        //NMK

    //edges
    std::vector<size_t> m_edgesNumFaces;                        //LNN
    std::vector<std::vector<size_t>> m_edgesFaces;              //LNE

    // faces
    std::vector<std::vector<size_t>> m_facesNodes;              //netcell%Nod, the nodes composing the faces, in ccw order
    std::vector<std::vector<size_t>> m_facesEdges;              //netcell%lin
    std::vector<Point>   m_facesCircumcenters;                  //xw
    std::vector<Point>   m_facesMasscenters;                    //the faces canters of mass

    size_t m_numFaces;                                          //NUMP
    std::vector<double> m_faceArea;                             //Face area

private:

    // Set node admin
    void NodeAdministration();
    
    void SortEdgesInCounterClockWiseOrder();

    void findFaces(const int& numEdges);
    
    void faceCircumcenters(const double& weightCircumCenter);

    void facesAreasAndCentersOfMass();
};
#endif