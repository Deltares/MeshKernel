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

        enum DeleteMeshOptions
        {
            AllVerticesInside = 0,
            FacesWithIncludedCircumcenters = 1,
            FacesCompletelyIncluded = 2
        };

        enum class AdministrationOptions
        {
            AdministrateMeshEdges,
            AdministrateMeshEdgesAndFaces
        };

        enum class NodeTypes
        {
            internalNode,
            onRing,
            cornerNode,
            hangingNode,
            other              // e.g. 1d node
        };

        /// <summary>
        /// Default constructor
        /// </summary>
        /// <returns></returns>
        Mesh();
        
        /// <summary>
        /// Converting constructor, from curvilinear grid to mesh (gridtonet)
        /// </summary>
        /// <param name="curvilinearGrid"></param>
        /// <param name="projection"></param>
        /// <returns></returns>
        Mesh(const CurvilinearGrid& curvilinearGrid, Projections projection);
        
        /// <summary>
        /// Create triangular grid from nodes (triangulatesamplestonetwork)
        /// </summary>
        /// <param name="nodes">Input nodes</param>
        /// <param name="polygons">Selection polygon</param>
        /// <param name="projection">Projection to use</param>
        /// <returns></returns>
        Mesh(std::vector<Point>& nodes, const Polygons& polygons, Projections projection);

        /// <summary>
        /// Add meshes: result is an unique mesh composed of the additions 
        /// firstMesh += secondmesh results in the second mesh being added to the first
        /// </summary>
        /// <param name="rhs"></param>
        /// <returns></returns>
        Mesh& operator+=(Mesh const& rhs);

        /// <summary>
        /// Set the mesh starting from the edges and the nodes
        /// </summary>
        /// <param name="edges">The input edges</param>
        /// <param name="nodes">The input nodes</param>
        /// <param name="projection">Projection to use</param>
        /// <param name="administration">Type of administration to perform</param>
        /// <returns>If the operation succeeded</returns>
        bool Set(const std::vector<Edge>& edges, const std::vector<Point>& nodes, Projections projection, AdministrationOptions administration = AdministrationOptions::AdministrateMeshEdgesAndFaces);
        
        /// <summary>
        /// Set internal flat arrays copies of nodes and edges, so they can be communicated to the front-end
        /// </summary>
        /// <param name="administrationOption">Type of administration to perform</param>
        /// <returns>If the operation succeeded</returns>
        bool SetFlatCopies(AdministrationOptions administrationOption);

        /// <summary>
        /// Perform mesh administration
        /// </summary>
        /// <param name="administrationOption">Type of administration to perform</param>
        /// <returns>If the operation succeeded</returns>
        bool Administrate(AdministrationOptions administrationOption);

        /// <summary>
        /// Compute face circumcenters, centers of mass and face areas
        /// </summary>
        void ComputeFaceCircumcentersMassCentersAndAreas();

        /// <summary>
        /// Find faces: constructs node to faces mapping. (findcells)
        /// </summary>
        void FindFaces();

        /// <summary>
        /// Gets the corners of a box bounding the mesh
        /// </summary>
        /// <param name="lowerLeft">Lower left corner</param>
        /// <param name="upperRight">Upper right corner</param>
        /// <returns>If the operation succeeded</returns>
        bool GetBoundingBox(Point& lowerLeft, Point& upperRight) const;

        /// <summary>
        /// Offset the x coordinates if ptojection is spherical
        /// </summary>
        /// <param name="minx"></param>
        /// <param name="miny"></param>
        /// <returns>If the operation succeeded</returns>
        bool OffsetSphericalCoordinates(double minx, double miny);

        /// <summary>
        /// Merge mesh nodes in a polygon (MERGENODESINPOLYGON)
        /// </summary>
        /// <param name="polygons">Polygon where to perform the merging</param>
        /// <returns>If the operation succeeded</returns>
        bool MergeNodesInPolygon(const Polygons& polygons);

        /// <summary>
        /// Merges two mesh nodes
        /// </summary>
        /// <param name="startNode">The index of the first node to be merged</param>
        /// <param name="endNode">The second of the second node to be merged</param>
        /// <returns>If the operation succeeded</returns>
        bool MergeTwoNodes(int startNode, int endNode);

        /// <summary>
        /// Make a new mesh (makenet)
        /// </summary>
        /// <param name="makeGridParametersNative">The structure containing the make grid parameters </param>
        /// <param name="polygons">The polygon to account for</param>
        /// <returns>If the operation succeeded</returns>
        bool MakeMesh(const GridGeomApi::MakeGridParametersNative& makeGridParametersNative, const Polygons& polygons);

        /// <summary>
        /// Deletes a mesh in a polygon using several options (delnet)
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="polygons">The polygon where to perform the operation</param>
        /// <param name="deletionOption">The deletion option</param>
        /// <param name="invertDeletion">Inverts the deletion of selected features</param>
        /// <returns>If the operation succeeded</returns>
        bool DeleteMesh(const Polygons& polygons, int deletionOption, bool invertDeletion);

        /// <summary>
        /// Connect two existing nodes, forming a new wdge index (connectdbn)
        /// </summary>
        /// <param name="startNode"></param>
        /// <param name="endNode"></param>
        /// <param name="newEdgeIndex"></param>
        /// <returns>If the operation succeeded</returns>
        bool ConnectNodes(int startNode, int endNode, int& newEdgeIndex);

        /// <summary>
        /// Insert a new node in the mesh (setnewpoint)
        /// </summary>
        /// <param name="newPoint">The coordinate of the new point</param>
        /// <param name="newNodeIndex">The index of the new node</param>
        /// <param name="updateRTree">Update the nodes RTree </param>
        /// <returns>If the operation succeeded</returns>
        bool InsertNode(const Point& newPoint, int& newNodeIndex, bool updateRTree = false);

        /// <summary>
        /// Delete a node
        /// </summary>
        /// <param name="nodeIndex"></param>
        /// <param name="updateRTree"></param>
        /// <returns>If the operation succeeded</returns>
        bool DeleteNode(int nodeIndex, bool updateRTree = false);

        /// <summary>
        /// Find the edge sharing two nodes
        /// </summary>
        /// <param name="firstNodeIndex">The index of the first node to be merged</param>
        /// <param name="secondNodeIndex">The index of the second node to be merged</param>
        /// <param name="edgeIndex"></param>
        /// <returns>If the operation succeeded</returns>
        bool FindEdge(int firstNodeIndex, int secondNodeIndex, int& edgeIndex) const;

        /// <summary>
        /// Move a node to a new location
        /// </summary>
        /// <param name="newPoint">The new location</param>
        /// <param name="nodeindex">The index of the node to move</param>
        /// <returns>If the operation succeeded</returns>
        bool MoveNode(Point newPoint, int nodeindex);

        /// <summary>
        /// Get the index of a node close to a point
        /// </summary>
        /// <param name="point">The starting point from where to start the search </param>
        /// <param name="searchRadius">The search radius</param>
        /// <param name="nodeIndex">The node index (-1 if no node is found)</param>
        /// <returns>If the operation succeeded</returns>
        bool GetNodeIndex(Point point, double searchRadius, int& nodeIndex);

        /// <summary>
        /// Deletes an edge
        /// </summary>
        /// <param name="edgeIndex"></param>
        /// <returns>If the operation succeeded</returns>
        bool DeleteEdge(int edgeIndex);

        /// <summary>
        /// Finds an edge close to a point
        /// </summary>
        /// <param name="point">The starting point from where to start the search</param>
        /// <param name="searchRadius">The search radius</param>
        /// <returns>The edge index (-1 if no edges is found)</returns>
        int FindEdgeCloseToAPoint(Point point, double searchRadius);

        /// <summary>
        /// Masks the edges of the faces included in a polygon
        /// </summary>
        /// <param name="polygons">The selection polygon</param>
        /// <param name="invertSelection">Invert selection</param>
        /// <param name="includeIntersected">Included the edges intersected by the polygon</param>
        /// <returns>If the operation succeeded</returns>
        bool MaskFaceEdgesInPolygon(const Polygons& polygons, bool invertSelection, bool includeIntersected);

        /// <summary>
        /// From the masked edges compute the masked nodes
        /// </summary>
        /// <returns>If the operation succeeded</returns>
        bool ComputeNodeMaskFromEdgeMask();

        /// <summary>
        /// For a face, fills the local caches (get_cellpolygon)
        /// </summary>
        /// <param name="faceIndex">The face index</param>
        /// <param name="polygonNodesCache">The node cache array filled with the nodes values</param>
        /// <param name="localNodeIndexsesCache">The consecutive node index in polygonNodesCache (0, 1, 2,...)</param>
        /// <param name="edgeIndexsesCache">The edge cache array filled with edge indexses</param>
        /// <param name="numClosedPolygonNodes">The number of valid values in the array above</param>
        /// <returns>If the operation succeeded</returns>
        bool FaceClosedPolygon(int faceIndex, 
            std::vector<Point>& polygonNodesCache, 
            std::vector<int>& localNodeIndexsesCache,
            std::vector<int>& edgeIndexsesCache,
            int& numClosedPolygonNodes) const;

        /// <summary>
        /// For a face, fills the polygon nodes cache
        /// </summary>
        /// <param name="faceIndex">The face index</param>
        /// <param name="polygonNodesCache">The node cache array filled with the nodes values</param>
        /// <param name="numClosedPolygonNodes">The number of valid values in the array above</param>
        /// <returns>If the operation succeeded</returns>
        bool FaceClosedPolygon(int faceIndex, std::vector<Point>& polygonNodesCache, int& numClosedPolygonNodes) const;

        /// <summary>
        /// Determine if a face is fully contained in polygon or not, based on m_nodeMask
        /// </summary>
        /// <param name="faceIndex">The face index</param>
        /// <returns>If the operation succeeded</returns>
        bool IsFullFaceNotInPolygon(int faceIndex) const;

        /// <summary>
        /// Mask all nodes in a polygon
        /// </summary>
        /// <param name="polygons">The input polygon</param>
        /// <param name="inside">Inside/outside option</param>
        /// <returns>If the operation succeeded</returns>
        bool MaskNodesInPolygons(const Polygons& polygons, bool inside);

        /// <summary>
        /// Find the common node two edges share
        /// </summary>
        /// <param name="firstEdgeIndex">The index of the first edge</param>
        /// <param name="secondEdgeIndex">The index of the second edge</param>
        /// <param name="node">The shared node (-1 if no node is found)</param>
        /// <returns>If the operation succeeded</returns>
        bool FindCommonNode(int firstEdgeIndex, int secondEdgeIndex, int& node) const;

        /// <summary>
        /// Compute the lenghts of all edges in one go
        /// </summary>
        bool ComputeEdgeLengths();

        /// <summary>
        /// Get the number of valid nodes
        /// </summary>
        /// <returns>If the operation succeeded</returns>
        int GetNumNodes() const { return m_numNodes; }

        /// <summary>
        /// Get the number of valid edges
        /// </summary>
        /// <returns>If the operation succeeded</returns>
        int GetNumEdges() const { return m_numEdges; }

        /// <summary>
        /// Get the number of valid faces
        /// </summary>
        /// <returns>If the operation succeeded</returns>
        int GetNumFaces() const { return m_numFaces; }

        /// <summary>
        /// Get the number of edges for a face
        /// </summary>
        /// <param name="faceIndex"></param>
        /// <returns>If the operation succeeded</returns>
        int GetNumFaceEdges(const int faceIndex) const { return m_numFacesNodes[faceIndex]; }

        /// <summary>
        /// Get the number of faces an edges shares
        /// </summary>
        /// <param name="edgeIndex"></param>
        /// <returns>If the operation succeeded</returns>
        int GetNumEdgesFaces(const int edgeIndex) const { return m_edgesNumFaces[edgeIndex]; }

        /// <summary>
        ///  Circumcenter of a face (getcircumcenter)
        /// </summary>
        /// <param name="polygon">Caching array for face nodes</param>
        /// <param name="middlePoints">Caching array for the edges middle points</param>
        /// <param name="normals">Caching array for normals</param>
        /// <param name="numNodes">Number of nodes</param>
        /// <param name="edgesNumFaces">For meshes, the number of faces sharing the polygon edge</param>
        /// <param name="weightCircumCenter">circumcenter weight</param>
        /// <param name="result">The computed polygon result</param>
        /// <returns>If the operation succeeded</returns>
        bool ComputeFaceCircumenter(std::vector<Point>& polygon,
            std::vector<Point>& middlePoints,
            std::vector<Point>& normals,
            int numNodes,
            const std::vector<int>& edgesNumFaces,
            const double weightCircumCenter,
            Point& result);

        // nodes
        std::vector<Point>              m_nodes;                    // (xk, yk)
        std::vector<std::vector<int>>   m_nodesEdges;               // (nod)
        std::vector<int>                m_nodesNumEdges;            // (nmk)
        std::vector<int>                m_nodeMask;                 // (kc)

        // edges
        std::vector<Edge>               m_edges;                    // (kn)
        std::vector<int>                m_edgesNumFaces;            // (lnn)
        std::vector<std::vector<int>>   m_edgesFaces;               // (lne)
        std::vector<double>             m_edgeLengths;
        std::vector<int>                m_edgeMask;                 // (lc)

        // faces
        std::vector<std::vector<int>>   m_facesNodes;               // netcell%Nod, the nodes composing the faces, in ccw order
        std::vector<int>                m_numFacesNodes;            // netcell%N
        std::vector<std::vector<int>>   m_facesEdges;               // netcell%lin
        std::vector<Point>              m_facesCircumcenters;       // xz  the face circumcenter
        std::vector<Point>              m_facesMassCenters;         // xzw the faces canters of mass

        std::vector<double> m_faceArea;                             // Face area   
        std::vector<int> m_nodesTypes;                              // Node types

        // flat arrays for communication with the client
        std::vector<double>              m_nodex;
        std::vector<double>              m_nodey;
        std::vector<double>              m_nodez;
        std::vector<int>                 m_edgeNodes;
        std::vector<int>                 m_faceNodes;
        std::vector<double>              m_facesCircumcentersx;
        std::vector<double>              m_facesCircumcentersy;
        std::vector<double>              m_facesCircumcentersz;

        Projections                      m_projection;
        std::vector<Point>               m_polygonNodesCache;       // Cache to store polygon points
        SpatialTrees::RTree              m_nodesRTree;              // Spatial tree to inquire node vertices

    private:

        /// <summary>
        /// Node administration (setnodadmin)
        /// </summary>
        void NodeAdministration();

        /// <summary>
        /// Sort edges in conterclockwise orther (Sort_links_ccw)
        /// </summary>
        void SortEdgesInCounterClockWiseOrder();

        /// <summary>
        /// Find cells recursive, works with an arbitrary number of edges
        /// </summary>
        /// <param name="startingNode">Th starting node</param>
        /// <param name="node">The current node</param>
        /// <param name="numEdges">The number of edges visited so far</param>
        /// <param name="previousEdge">The previously visited edge</param>
        /// <param name="edges">The vector storing the current edges forming a face</param>
        /// <param name="nodes">The vector storing the current nodes forming a face</param>
        /// <param name="sortedEdges">A caching array used for sorting the edges and inquire if an edge has been already visited</param>
        /// <param name="sortedNodes">A caching array used for sorting the nodes and inquire if a node has been already visited</param>
        /// <returns>If the operation succeeded</returns>
        bool FindFacesRecursive(int startingNode, 
                                int node, 
                                int numEdges, 
                                int previousEdge, 
                                std::vector<int>& edges, 
                                std::vector<int>& nodes,
                                std::vector<int>& sortedEdges,
                                std::vector<int>& sortedNodes);

        /// <summary>
        /// Classifies the nodes (makenetnodescoding)
        /// </summary>
        /// <returns>If the operation succeeded</returns>
        bool ClassifyNodes();

        /// <summary>
        /// Checks if a triangle has an acute angle (checktriangle)
        /// </summary>
        /// <returns>If the operation succeeded</returns>
        bool CheckTriangle(const std::vector<int>& faceNodes, const std::vector<Point>& nodes);

        /// <summary>
        /// Removes all invalid nodes and edges
        /// </summary>
        /// <returns>If the operation succeeded</returns>
        bool RemoveInvalidNodesAndEdges();

        /// <summary>
        /// Refresh nodes RTree if needed
        /// </summary>
        /// <returns>If the operation succeeded</returns>
        bool RefreshNodesRTreeIfNeeded();

        int m_numFaces = 0;                                       // number of valid faces (nump)
        int m_numNodes = 0;                                       // Number of valid nodes in m_nodes
        int m_numEdges = 0;                                       // Number of valid edges in m_edges

    };
}
