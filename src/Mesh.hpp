#pragma once

#include <vector>
#include "MakeGridParametersNative.hpp"
#include "Entities.hpp"
#include "SpatialTrees.hpp"

namespace GridGeom 
{

    class CurvilinearGrid;
    class Polygons;
    class MakeGridParametersNative;
    class GeometryListNative;


    class Mesh
    {

    public:

        Mesh(){}

        //gridtonet
        Mesh(const CurvilinearGrid& curvilinearGrid, Projections projection);
        
        // triangulatesamplestonetwork
        Mesh(std::vector<Point>& nodes, const Polygons& polygons, Projections projection);

        
        bool Set(const std::vector<Edge>& edges, const std::vector<Point>& nodes, Projections projection);

        bool IsAdministrationDone() const;
        
        bool Administrate();
        
        bool SetFlatCopies();
        
        bool DeleteFlatCopies();

        void ComputeFaceCircumcentersMassCentersAreas();

        void FindFaces();

        int GetNumNodes() const { return m_nodes.size(); }

        int GetNumEdges() const { return m_edges.size(); }

        int GetNumFaces() const { return m_numFaces; }

        int GetNumFaceEdges(const int faceIndex) const { return m_numFacesNodes[faceIndex]; }

        int GetNumEdgesFaces(const int edgeIndex) const { return m_edgesNumFaces[edgeIndex]; }

        int GetRTreeSize() const { return m_rtree.Size(); }

        bool GetBoundingBox(Point& lowerLeft, Point& upperRight) const;

        bool OffsetSphericalCoordinates(double minx, double miny);

        ///MERGENODESINPOLYGON
        bool MergeNodesInPolygon(const Polygons& polygons);

        ///mergenodes
        inline bool MergeTwoNodes(int startNode, int endNode);

        ///makenet
        bool MakeMesh(const GridGeomApi::MakeGridParametersNative& makeGridParametersNative, const Polygons& polygons);

        ///connectdbn
        bool ConnectNodes(int startNode, int endNode, int& newEdgeIndex);

        ///setnewpoint
        bool InsertNode(const Point& newPoint, int& newNodeIndex);

        ///setnewpoint
        bool DeleteNode(int nodeIndex);

        ///DELELEM
        bool DeleteEdge(int startNode, int endNode);

        // find an edge
        bool FindEdge(int firstNodeIndex, int secondNodeIndex, int& edgeIndex) const;

        ///get_cellpolygon
        //need to account for spherical coordinates. Build a polygon around a face
        bool FacePolygon(int faceIndex, std::vector<Point>& polygonNodesCache, std::vector<int>& localNodeIndexsesCache, std::vector<int>& nodeIndexsesCache) const;

        bool FacePolygon(int faceIndex, std::vector<Point>& polygonNodesCache, int& numPolygonPoints) const;

        bool IsFullFaceNotInPolygon(int faceIndex) const;

        bool SelectNodesInPolygon(const Polygons& polygons, int inside);

        bool FindCommonNode(int firstEdgeIndex, int secondEdgeIndex, int& node) const;

        bool GetNodeIndex(Point point, double searchRadius, int& vertexIndex);

        bool BuildRTree();

        bool InsertMissingNodesInRTree();


        std::vector<Edge>  m_edges;                                 // KN
        std::vector<Point> m_nodes;                                 // KN
        std::vector<std::vector<int>> m_nodesEdges;                 // NOD
        std::vector<int> m_nodesNumEdges;                           // NMK
        std::vector<int> m_nodeMask;                                // KC

        //edges
        std::vector<int> m_edgesNumFaces;                           // LNN
        std::vector<std::vector<int>> m_edgesFaces;                 // LNE
        std::vector<double> m_edgeLengths;

        // faces
        std::vector<std::vector<int>> m_facesNodes;                 // netcell%Nod, the nodes composing the faces, in ccw order
        std::vector<int>              m_numFacesNodes;              // netcell%N
        std::vector<std::vector<int>> m_facesEdges;                 // netcell%lin
        std::vector<Point>            m_facesCircumcenters;         // xz  the face circumcenter
        std::vector<Point>            m_facesMassCenters;           // xzw the faces canters of mass

        std::vector<double> m_faceArea;                             // Face area   
        std::vector<int> m_nodesTypes;                              // Node types,  1=internal, 2=on ring, 3=corner point, 0/-1=other (e.g. 1d)

        // Used for internal state
        std::vector<double> m_nodex;
        std::vector<double> m_nodey;
        std::vector<double> m_nodez;
        std::vector<int>    m_edgeNodes;

        std::vector<Point> m_polygonNodesCache;                     // polygon cache

        // Used for triangular grids
        double m_triangleMinimumAngle = 5.0;                        // minimum angle of created triangles. If minimum angle > maximum angle, no check 
        double m_triangleMaximumAngle = 150.0;                      // maximum angle of created triangles

        // spatial tree to inquire node vertices
        SpatialTrees::RTree m_rtree;

        bool m_isAdministrationDone = false;


        Projections m_projection;

    private:

        // Set node admin
        void NodeAdministration();

        // Sort_links_ccw
        void SortEdgesInCounterClockWiseOrder();

        // find cells
        void FindFaces(int numEdges);

        // find cells recursive
        bool FindFacesRecursive(int startingNode, int node, int numEdges, int previousEdge, 
            std::vector<int>& edges, 
            std::vector<int>& nodes,
            std::vector<int>& sortedEdges,
            std::vector<int>& sortedNodes);

        /// @brief makenetnodescoding: computes node types
        bool ClassifyNodes();

        /// CHECKTRIANGLE
        bool CheckTriangle(const std::vector<int>& faceNodes, const std::vector<Point>& nodes);

        bool ComputeEdgeLengths();

        double m_dcenterinside = 1.0;

        int m_numFaces;                                             // NUMP
        int m_numNodes;                                             // valid nodes

    };
}
