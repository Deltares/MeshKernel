//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2020.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation version 3.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// contact: delft3d.support@deltares.nl
// Stichting Deltares
// P.O. Box 177
// 2600 MH Delft, The Netherlands
//
// All indications and logos of, and references to, "Delft3D" and "Deltares"
// are registered trademarks of Stichting Deltares, and remain the property of
// Stichting Deltares. All rights reserved.
//
//------------------------------------------------------------------------------

#pragma once

#include <vector>
#include "MakeGridParametersNative.hpp"
#include "Entities.hpp"
#include "SpatialTrees.hpp"

namespace meshkernel
{
    // Forward declarations
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
            other
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
        /// Add meshes: result is a mesh composed of the additions
        /// firstMesh += secondmesh results in the second mesh being added to the first
        /// </summary>
        /// <param name="rhs">The mesh to add</param>
        /// <returns>The resulting mesh</returns>
        Mesh& operator+=(Mesh const& rhs);

        /// <summary>
        /// Set the mesh starting from the edges and nodes
        /// </summary>
        /// <param name="edges">The input edges</param>
        /// <param name="nodes">The input nodes</param>
        /// <param name="projection">Projection to use</param>
        /// <param name="administration">Type of administration to perform</param>
        /// <returns>If the method succeeded</returns>
        bool Set(const std::vector<Edge>& edges, const std::vector<Point>& nodes, Projections projection, AdministrationOptions administration = AdministrationOptions::AdministrateMeshEdgesAndFaces);

        /// <summary>
        /// Set internal flat copies of nodes and edges, so the pointer to the first entry is communicated with the front-end
        /// </summary>
        /// <param name="administrationOption">Type of administration to perform</param>
        /// <returns>If the method succeeded</returns>
        bool SetFlatCopies(AdministrationOptions administrationOption);

        /// <summary>
        /// Perform mesh administration
        /// </summary>
        /// <param name="administrationOption">Type of administration to perform</param>
        /// <returns>If the method succeeded</returns>
        bool Administrate(AdministrationOptions administrationOption);

        /// <summary>
        /// Compute face circumcenters, centers of mass and face areas
        /// </summary>
        void ComputeFaceCircumcentersMassCentersAndAreas();

        /// <summary>
        /// Find faces: constructs the m_facesNodes mapping. (findcells)
        /// </summary>
        void FindFaces();

        /// <summary>
        /// Gets the corners of a box bounding the mesh
        /// </summary>
        /// <param name="lowerLeft">Lower left corner</param>
        /// <param name="upperRight">Upper right corner</param>
        /// <returns>If the method succeeded</returns>
        bool GetBoundingBox(Point& lowerLeft, Point& upperRight) const;

        /// <summary>
        /// Offset the x coordinates if m_projection is spherical
        /// </summary>
        /// <param name="minx"></param>
        /// <param name="miny"></param>
        /// <returns>If the method succeeded</returns>
        bool OffsetSphericalCoordinates(double minx, double miny);

        /// <summary>
        /// Merge close mesh nodes inside a polygon (MERGENODESINPOLYGON)
        /// </summary>
        /// <param name="polygons">Polygon where to perform the merging</param>
        /// <returns>If the method succeeded</returns>
        bool MergeNodesInPolygon(const Polygons& polygons);

        /// <summary>
        /// Merges two mesh nodes
        /// </summary>
        /// <param name="startNode">The index of the first node to be merged</param>
        /// <param name="endNode">The second of the second node to be merged</param>
        /// <returns>If the method succeeded</returns>
        bool MergeTwoNodes(int startNode, int endNode);

        /// <summary>
        /// Make a new rectangular mesh, composed of quads (makenet)
        /// </summary>
        /// <param name="makeGridParametersNative">The structure containing the make grid parameters </param>
        /// <param name="polygons">The polygon to account for</param>
        /// <returns>If the method succeeded</returns>
        bool MakeMesh(const meshkernelapi::MakeGridParametersNative& makeGridParametersNative, const Polygons& polygons);

        /// <summary>
        /// Deletes a mesh in a polygon, using several options (delnet)
        /// </summary>
        /// <param name="polygons">The polygon where to perform the operation</param>
        /// <param name="deletionOption">The deletion option</param>
        /// <param name="invertDeletion">Inverts the selected node to delete (instead of outside the polygon, inside the polygon) </param>
        /// <returns>If the method succeeded</returns>
        bool DeleteMesh(const Polygons& polygons, int deletionOption, bool invertDeletion);

        /// <summary>
        /// Connect two existing nodes, forming a new edge (connectdbn)
        /// </summary>
        /// <param name="startNode">The start node index</param>
        /// <param name="endNode">The end node index</param>
        /// <param name="newEdgeIndex">The index of the new edge</param>
        /// <returns>If the method succeeded</returns>
        bool ConnectNodes(int startNode, int endNode, int& newEdgeIndex);

        /// <summary>
        /// Insert a new node in the mesh (setnewpoint)
        /// </summary>
        /// <param name="newPoint">The coordinate of the new point</param>
        /// <param name="newNodeIndex">The index of the new node</param>
        /// <param name="updateRTree">Update m_nodesRTree</param>
        /// <returns>If the method succeeded</returns>
        bool InsertNode(const Point& newPoint, int& newNodeIndex, bool updateRTree = false);

        /// <summary>
        /// Delete a node
        /// </summary>
        /// <param name="nodeIndex">The index of the node to delete</param>
        /// <param name="updateRTree">Update m_nodesRTree</param>
        /// <returns>If the method succeeded</returns>
        bool DeleteNode(int nodeIndex);

        /// <summary>
        /// Find the edge sharing two nodes
        /// </summary>
        /// <param name="firstNodeIndex">The index of the first node</param>
        /// <param name="secondNodeIndex">The index of the second node</param>
        /// <param name="edgeIndex">The edge index</param>
        /// <returns>If the method succeeded</returns>
        bool FindEdge(int firstNodeIndex, int secondNodeIndex, int& edgeIndex) const;

        /// <summary>
        /// Move a node to a new location
        /// </summary>
        /// <param name="newPoint">The new location</param>
        /// <param name="nodeindex">The index of the node to move</param>
        /// <returns>If the method succeeded</returns>
        bool MoveNode(Point newPoint, int nodeindex);

        /// <summary>
        /// Get the index of a node close to a point
        /// </summary>
        /// <param name="point">The starting point from where to start the search </param>
        /// <param name="searchRadius">The search radius</param>
        /// <param name="nodeIndex">The node index (-1 if no node is found)</param>
        /// <returns>If the method succeeded</returns>
        bool GetNodeIndex(Point point, double searchRadius, int& nodeIndex);

        /// <summary>
        /// Deletes an edge
        /// </summary>
        /// <param name="edgeIndex">The edge index</param>
        /// <returns>If the method succeeded</returns>
        bool DeleteEdge(int edgeIndex);

        /// <summary>
        /// Finds the closest edge close to a point
        /// </summary>
        /// <param name="point">The starting point from where to start the search</param>
        /// <param>The edge index (-1 if no edges is found)</param>
        /// <returns>If the method succeeded</returns>
        bool FindEdgeCloseToAPoint(Point point, int& edgeIndex);

        /// <summary>
        /// Masks the edges of all faces included in a polygon
        /// </summary>
        /// <param name="polygons">The selection polygon</param>
        /// <param name="invertSelection">Invert selection</param>
        /// <param name="includeIntersected">Included the edges intersected by the polygon</param>
        /// <returns>If the method succeeded</returns>
        bool MaskFaceEdgesInPolygon(const Polygons& polygons, bool invertSelection, bool includeIntersected);

        /// <summary>
        /// From the masked edges compute the masked nodes
        /// </summary>
        /// <returns>If the method succeeded</returns>
        bool ComputeNodeMaskFromEdgeMask();

        /// <summary>
        /// For a face, fills the local caches (get_cellpolygon)
        /// </summary>
        /// <param name="faceIndex">The face index</param>
        /// <param name="polygonNodesCache">The node cache array filled with the nodes values</param>
        /// <param name="localNodeIndicesCache">The consecutive node index in polygonNodesCache (0, 1, 2,...)</param>
        /// <param name="edgeIndicesCache">The edge cache array filled with edge indices</param>
        /// <param name="numClosedPolygonNodes">The number of valid values in the array above</param>
        /// <returns>If the method succeeded</returns>
        bool FaceClosedPolygon(int faceIndex,
                               std::vector<Point>& polygonNodesCache,
                               std::vector<int>& localNodeIndicesCache,
                               std::vector<int>& edgeIndicesCache,
                               int& numClosedPolygonNodes) const;

        /// <summary>
        /// For a face, fills the polygonNodesCache with the face nodes
        /// </summary>
        /// <param name="faceIndex">The face index</param>
        /// <param name="polygonNodesCache">The cache array to be filled </param>
        /// <param name="numClosedPolygonNodes">The number of valid face nodes</param>
        /// <returns>If the method succeeded</returns>
        bool FaceClosedPolygon(int faceIndex,
                               std::vector<Point>& polygonNodesCache,
                               int& numClosedPolygonNodes) const;

        /// <summary>
        /// Determine if a face is fully contained in polygon or not, based on m_nodeMask
        /// </summary>
        /// <param name="faceIndex">The face index</param>
        /// <returns>If the method succeeded</returns>
        bool IsFullFaceNotInPolygon(int faceIndex) const;

        /// <summary>
        /// Mask all nodes in a polygon
        /// </summary>
        /// <param name="polygons">The input polygon</param>
        /// <param name="inside">Inside/outside option</param>
        /// <returns>If the method succeeded</returns>
        bool MaskNodesInPolygons(const Polygons& polygons, bool inside);

        /// <summary>
        /// Find the common node two edges share
        /// </summary>
        /// <param name="firstEdgeIndex">The index of the first edge</param>
        /// <param name="secondEdgeIndex">The index of the second edge</param>
        /// <param name="node">The shared node (-1 if no node is found)</param>
        /// <returns>If the method succeeded</returns>
        bool FindCommonNode(int firstEdgeIndex, int secondEdgeIndex, int& node) const;

        /// <summary>
        /// Compute the lengths of all edges in one go
        /// </summary>
        bool ComputeEdgeLengths();

        /// <summary>
        /// Computes the edges centers
        /// </summary>
        /// <returns>If the method succeeded</returns>
        bool ComputeEdgesCenters();

        /// <summary>
        /// Get the number of valid nodes
        /// </summary>
        /// <returns>The number of valid node</returns>
        int GetNumNodes() const { return m_numNodes; }

        /// <summary>
        /// Get the number of valid edges
        /// </summary>
        /// <returns>The number of valid edges</returns>
        int GetNumEdges() const { return m_numEdges; }

        /// <summary>
        /// Get the number of valid faces
        /// </summary>
        /// <returns>The number of valid faces</returns>
        int GetNumFaces() const { return m_numFaces; }

        /// <summary>
        /// Get the number of edges for a face
        /// </summary>
        /// <param name="faceIndex">The face index</param>
        /// <returns>The number of edges for a face</returns>
        int GetNumFaceEdges(const int faceIndex) const { return m_numFacesNodes[faceIndex]; }

        /// <summary>
        /// Get the number of faces an edges shares
        /// </summary>
        /// <param name="edgeIndex">The edge index</param>
        /// <returns>The number of faces an edges shares</returns>
        int GetNumEdgesFaces(const int edgeIndex) const { return m_edgesNumFaces[edgeIndex]; }

        /// <summary>
        ///  Circumcenter of a face (getcircumcenter)
        /// </summary>
        /// <param name="polygon">Cache storing the face nodes</param>
        /// <param name="middlePoints">Caching array for the edges middle points</param>
        /// <param name="normals">Caching array for normals</param>
        /// <param name="numNodes">Number of valid nodes in the cache</param>
        /// <param name="edgesNumFaces">For meshes, the number of faces sharing the edges</param>
        /// <param name="weightCircumCenter">Circumcenter weight</param>
        /// <param name="result">The computed circumcenter</param>
        /// <returns>If the method succeeded</returns>
        bool ComputeFaceCircumenter(std::vector<Point>& polygon,
                                    std::vector<Point>& middlePoints,
                                    std::vector<Point>& normals,
                                    int numNodes,
                                    const std::vector<int>& edgesNumFaces,
                                    double weightCircumCenter,
                                    Point& result) const;

        /// <summary>
        /// Computes m_nodesNodes, see class members
        /// </summary>
        /// <returns>If the method succeeded</returns>
        bool ComputeNodeNeighbours();

        /// <summary>
        /// Get the orthogonality values, the inner product of edges and segments connecting the face circumcenters
        /// </summary>
        /// <param name="orthogonality">The edge orthogonality, passed to the client</param>
        /// <returns>If the method succeeded</returns>
        bool GetOrthogonality(double* orthogonality);

        /// <summary>
        /// Gets the smoothness values, ratios of the face areas
        /// </summary>
        /// <param name="smoothness">The smoothness at the edges</param>
        /// <returns>If the method succeeded</returns>
        bool GetSmoothness(double* smoothness);

        /// <summary>
        /// Gets the aspect ratios, the ratio edges to segments connecting the face circumcenters lengths
        /// </summary>
        /// <param name="aspectRatio">The aspect ratios</param>
        /// <returns>If the method succeeded</returns>
        bool GetAspectRatios(std::vector<double>& aspectRatios);

        /// <summary>
        /// Classifies the nodes (makenetnodescoding)
        /// </summary>
        /// <returns>If the method succeeded</returns>
        bool ClassifyNodes();

        /// <summary>
        ///  Sort edges in conterclockwise orther (Sort_links_ccw)
        /// </summary>
        /// <param name="nodeIndex">The node index for which sorting should take place</param>
        void SortEdgesInCounterClockWiseOrder(int nodeIndex);

        /// <summary>
        /// Transform non-triangular faces in triangular faces
        /// </summary>
        /// <returns>If the method succeeded</returns>
        bool TriangulateFaces();

        // nodes
        std::vector<Point> m_nodes;                 // The mesh nodes (xk, yk)
        std::vector<std::vector<int>> m_nodesEdges; // For each node, the indices of connected edges (nod%lin)
        std::vector<int> m_nodesNumEdges;           // For each node, the number of connected edges (nmk)
        std::vector<int> m_nodeMask;                // The node mask (kc)
        std::vector<std::vector<int>> m_nodesNodes; // For each node, its neighbours
        std::vector<int> m_nodesTypes;              // The node types (nb)

        // edges
        std::vector<Edge> m_edges;                  // The edges, defined as first and second node(kn)
        std::vector<std::vector<int>> m_edgesFaces; // For each edge, the shared face index (lne)
        std::vector<int> m_edgesNumFaces;           // For each edge, the number of shared faces(lnn)
        std::vector<double> m_edgeLengths;          // The edge lenghts
        std::vector<int> m_edgeMask;                // The edge mask (lc)
        std::vector<Point> m_edgesCenters;          // The edges centers

        // faces
        std::vector<std::vector<int>> m_facesNodes; // The nodes composing the faces, in ccw order (netcell%Nod)
        std::vector<int> m_numFacesNodes;           // The number of nodes composing the face (netcell%N)
        std::vector<std::vector<int>> m_facesEdges; // The edge indices composing the face (netcell%lin)
        std::vector<Point> m_facesCircumcenters;    // The face circumcenters the face circumcenter (xz, yz)
        std::vector<Point> m_facesMassCenters;      // The faces centers of mass (xzw, yzw)
        std::vector<double> m_faceArea;             // The face area
        std::vector<Point> m_polygonNodesCache;     // Cache to store the face nodes

        // vectors for communicating with the client
        std::vector<double> m_nodex;               // The nodes x-coordinate
        std::vector<double> m_nodey;               // The nodes y-coordinate
        std::vector<double> m_nodez;               // The nodes z-coordinate
        std::vector<int> m_edgeNodes;              // For each edge, the nodes
        std::vector<int> m_faceNodes;              // For each face, the nodes
        std::vector<double> m_facesCircumcentersx; // The circumcenters x-coordinate
        std::vector<double> m_facesCircumcentersy; // The circumcenters y-coordinate
        std::vector<double> m_facesCircumcentersz; // The circumcenters z-coordinate

        Projections m_projection; // The projection used

        SpatialTrees::RTree m_nodesRTree; // Spatial R-Tree used to inquire node vertices
        SpatialTrees::RTree m_edgesRTree; // Spatial R-Tree used to inquire edges centers

        int m_maxNumNeighbours = 0;

    private:
        /// <summary>
        /// Node administration (setnodadmin)
        /// </summary>
        void NodeAdministration();

        /// <summary>
        /// Find cells recursive, works with an arbitrary number of edges
        /// </summary>
        /// <param name="startingNode">The starting node</param>
        /// <param name="node">The current node</param>
        /// <param name="numEdges">The number of edges visited so far</param>
        /// <param name="previousEdge">The previously visited edge</param>
        /// <param name="edges">The vector storing the current edges forming a face</param>
        /// <param name="nodes">The vector storing the current nodes forming a face</param>
        /// <param name="sortedEdges">The caching array used for sorting the edges, used to inquire if an edge has been already visited</param>
        /// <param name="sortedNodes">The caching array used for sorting the nodes, used to inquire if a node has been already visited</param>
        /// <returns>If the method succeeded</returns>
        bool FindFacesRecursive(int startingNode,
                                int node,
                                int numEdges,
                                int previousEdge,
                                std::vector<int>& edges,
                                std::vector<int>& nodes,
                                std::vector<int>& sortedEdges,
                                std::vector<int>& sortedNodes);

        /// <summary>
        /// Checks if a triangle has an acute angle (checktriangle)
        /// </summary>
        /// <returns>If the method succeeded</returns>
        bool CheckTriangle(const std::vector<int>& faceNodes, const std::vector<Point>& nodes) const;

        /// <summary>
        /// Removes all invalid nodes and edges
        /// </summary>
        /// <returns>If the method succeeded</returns>
        bool RemoveInvalidNodesAndEdges();

        int m_numFaces = 0;               // number of valid faces (nump)
        int m_numNodes = 0;               // Number of valid nodes in m_nodes
        int m_numEdges = 0;               // Number of valid edges in m_edges
        std::vector<double> m_edgeAngles; // internal cache for sorting the edges around nodes

        bool m_nodesRTreeRequiresUpdate = false; //m_nodesRTree requires an update
        bool m_edgesRTreeRequiresUpdate = false; //m_edgesRTree requires an update
    };
} // namespace meshkernel
