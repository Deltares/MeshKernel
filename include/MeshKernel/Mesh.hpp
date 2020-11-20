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
#include <MeshKernel/MakeGridParametersNative.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/SpatialTrees.hpp>

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

        /// @brief Set the mesh starting from the edges and nodes
        /// @param[in] edges">The input edges</param>
        /// @param[in] nodes The input nodes
        /// @param[in] projection Projection to use
        /// @param[in] administration Type of administration to perform
        void Set(const std::vector<Edge>& edges, const std::vector<Point>& nodes, Projections projection, AdministrationOptions administration = AdministrationOptions::AdministrateMeshEdgesAndFaces);

        /// @brief Set internal flat copies of nodes and edges, so the pointer to the first entry is communicated with the front-end
        /// @param administrationOption Type of administration to perform
        void SetFlatCopies(AdministrationOptions administrationOption);

        /// @brief Perform mesh administration
        /// @param administrationOption Type of administration to perform
        void Administrate(AdministrationOptions administrationOption);

        /// @brief Compute face circumcenters
        void ComputeFaceCircumcentersMassCentersAndAreas(bool computeMassCenters = false);

        /// <summary>
        /// Find faces: constructs the m_facesNodes mapping, face mass centers and areas (findcells)
        /// </summary>
        void FindFaces();

        /// @brief Gets the corners of a box bounding the mesh
        /// @param[out] lowerLeft Lower left corner
        /// @param[out] upperRight Upper right corner
        void GetBoundingBox(Point& lowerLeft, Point& upperRight) const;

        /// @brief Offset the x coordinates if m_projection is spherical
        /// @param[in] minx
        /// @param[in] miny
        void OffsetSphericalCoordinates(double minx, double miny);

        /// @brief Merge close mesh nodes inside a polygon (MERGENODESINPOLYGON)
        /// @param[in] polygons Polygon where to perform the merging
        void MergeNodesInPolygon(const Polygons& polygons);

        /// @brief Merges two mesh nodes
        /// @param[in] startNode The index of the first node to be merged
        /// @param[in] endNode The second of the second node to be merged
        void MergeTwoNodes(int startNode, int endNode);

        /// @brief Make a new rectangular mesh, composed of quads (makenet)
        /// @param[in] makeGridParametersNative The structure containing the make grid parameters
        /// @param[in] polygons The polygon to account for
        void MakeMesh(const meshkernelapi::MakeGridParametersNative& makeGridParametersNative, const Polygons& polygons);

        /// @brief Deletes a mesh in a polygon, using several options (delnet)
        /// @param[in] polygons The polygon where to perform the operation
        /// @param[in] deletionOption The deletion option
        /// @param[in] invertDeletion Inverts the selected node to delete (instead of outside the polygon, inside the polygon)
        void DeleteMesh(const Polygons& polygons, int deletionOption, bool invertDeletion);

        /// @brief Connect two existing nodes, forming a new edge (connectdbn)
        /// @param[in] startNode The start node index
        /// @param[in] endNode The end node index
        /// @param[out] newEdgeIndex The index of the new edge
        void ConnectNodes(int startNode, int endNode, int& newEdgeIndex);

        /// @brief Insert a new node in the mesh (setnewpoint)
        /// @param[in] newPoint The coordinate of the new point
        /// @param[out] newNodeIndex The index of the new node
        void InsertNode(const Point& newPoint, int& newNodeIndex);

        /// @brief Delete a node
        /// @param[in] nodeIndex The index of the node to delete
        void DeleteNode(int nodeIndex);

        /// @brief Find the edge sharing two nodes
        /// @param[in] firstNodeIndex The index of the first node
        /// @param[in] secondNodeIndex The index of the second node
        /// @param[out] edgeIndex The edge index
        void FindEdge(int firstNodeIndex, int secondNodeIndex, int& edgeIndex) const;

        /// @brief Move a node to a new location
        /// @param[in] newPoint The new location
        /// @param[in] nodeindex The index of the node to move
        void MoveNode(Point newPoint, int nodeindex);

        /// @brief Get the index of a node close to a point
        /// @param[in] point The starting point from where to start the search
        /// @param[in] searchRadius The search radius
        /// @returns The index of the closest node
        [[nodiscard]] int GetNodeIndex(Point point, double searchRadius);

        /// @brief Deletes an edge
        /// @param[in] edgeIndex The edge index
        void DeleteEdge(int edgeIndex);

        /// Finds the closest edge close to a point
        /// @param[in] point The starting point from where to start the search
        /// @returns The index of the closest edge
        [[nodiscard]] int FindEdgeCloseToAPoint(Point point);

        /// @brief Masks the edges of all faces included in a polygon
        /// @param polygons The selection polygon
        /// @param invertSelection Invert selection
        /// @param includeIntersected Included the edges intersected by the polygon
        void MaskFaceEdgesInPolygon(const Polygons& polygons, bool invertSelection, bool includeIntersected);

        /// @brief From the masked edges compute the masked nodes
        void ComputeNodeMaskFromEdgeMask();

        /// @brief For a face, fills the local caches (get_cellpolygon)
        /// @param[in] faceIndex The face index
        /// @param[out] polygonNodesCache The node cache array filled with the nodes values
        /// @param[out] localNodeIndicesCache The consecutive node index in polygonNodesCache (0, 1, 2,...)
        /// @param[out] edgeIndicesCache The edge cache array filled with edge indices
        /// @param[out] numClosedPolygonNodes The number of valid values in the array above
        void FaceClosedPolygon(int faceIndex,
                               std::vector<Point>& polygonNodesCache,
                               std::vector<int>& localNodeIndicesCache,
                               std::vector<int>& edgeIndicesCache,
                               int& numClosedPolygonNodes) const;

        /// @brief For a face, fills the polygonNodesCache with the face nodes
        /// @param[in] faceIndex The face index
        /// @param[out] polygonNodesCache The cache array to be filled
        /// @param[out] numClosedPolygonNodes The number of valid face nodes
        void FaceClosedPolygon(int faceIndex,
                               std::vector<Point>& polygonNodesCache,
                               int& numClosedPolygonNodes) const;

        /// @brief Determine if a face is fully contained in polygon or not, based on m_nodeMask
        /// @param[in] faceIndex The face index
        /// @returns If the face is fully contained in the polygon or not
        [[nodiscard]] bool IsFullFaceNotInPolygon(int faceIndex) const;

        /// @brief Mask all nodes in a polygon
        /// @param[in] polygons The input polygon
        /// @param[in] inside Inside/outside option
        void MaskNodesInPolygons(const Polygons& polygons, bool inside);

        /// @brief Find the common node two edges share
        /// This method uses return parameters since the success is evaluated in a hot loop
        /// @param[in] firstEdgeIndex The index of the first edge
        /// @param[in] secondEdgeIndex The index of the second edge
        /// @param[out] node The shared node (-1 if no node is found)
        /// @return If the node could be found
        [[nodiscard]] bool FindCommonNode(int firstEdgeIndex, int secondEdgeIndex, int& node) const;

        /// @brief Compute the lengths of all edges in one go
        void ComputeEdgeLengths();

        /// @brief Computes the edges centers
        void ComputeEdgesCenters();

        /// @brief Get the number of valid nodes
        /// @return The number of valid node
        [[nodiscard]] int GetNumNodes() const { return m_numNodes; }

        /// @brief Get the number of valid edges
        /// @return The number of valid edges
        [[nodiscard]] int GetNumEdges() const { return m_numEdges; }

        /// @brief Get the number of valid faces
        /// @return The number of valid faces
        [[nodiscard]] int GetNumFaces() const { return m_numFaces; }

        /// @brief Get the number of edges for a face
        /// @param[in] faceIndex The face index
        /// @return The number of valid faces
        [[nodiscard]] int GetNumFaceEdges(int faceIndex) const { return m_numFacesNodes[faceIndex]; }

        /// @brief Get the number of faces an edges shares
        /// @param[in] edgeIndex The edge index
        /// @return The number of faces an edges shares
        [[nodiscard]] int GetNumEdgesFaces(int edgeIndex) const { return m_edgesNumFaces[edgeIndex]; }

        /// @brief Inquire if an edge is on boundary
        /// @param the edge index
        /// @return if the edge os on boundary
        [[nodiscard]] bool IsEdgeOnBoundary(int edge) const { return m_edgesNumFaces[edge] == 1; }

        /// @brief Inquire if a face is on boundary
        /// @param the edge index
        /// @return if the edge os on boundary
        [[nodiscard]] bool IsFaceOnBoundary(int face) const;

        /// @brief Circumcenter of a face (getcircumcenter)
        /// @param[in,out] polygon Cache storing the face nodes
        /// @param[in,out] middlePoints Caching array for the edges middle points
        /// @param[in,out] normals Caching array for normals
        /// @param[in] numNodes Number of valid nodes in the cache
        /// @param[in] edgesNumFaces For meshes, the number of faces sharing the edges
        /// @param[in] weightCircumCenter Circumcenter weight
        /// @returns The computed circumcenter
        [[nodiscard]] Point ComputeFaceCircumenter(std::vector<Point>& polygon,
                                                   std::vector<Point>& middlePoints,
                                                   std::vector<Point>& normals,
                                                   int numNodes,
                                                   const std::vector<int>& edgesNumFaces,
                                                   double weightCircumCenter) const;

        /// @brief Gets the mass centers of obtuse triangles
        /// @returns The center of obtuse triangles
        [[nodiscard]] std::vector<Point> GetObtuseTrianglesCenters();

        /// @brief Gets the edges crossing the small flow edges
        /// @param[in] smallFlowEdgesThreshold The configurable threshold for detecting the small flow edges
        /// @returns The indexes of the edges crossing small flow edges
        [[nodiscard]] std::vector<int> GetEdgesCrossingSmallFlowEdges(double smallFlowEdgesThreshold);

        /// @brief Gets the flow edges centers from the crossing edges
        /// @param[in] edges The crossing edges indexes
        /// @returns The centers of the flow edges
        [[nodiscard]] std::vector<Point> GetFlowEdgesCenters(const std::vector<int>& edges) const;

        /// @brief Remove small flow edges (removesmallflowlinks, part 1)
        /// @param[in] smallFlowEdgesThreshold The configurable threshold for detecting the small flow edges
        void RemoveSmallFlowEdges(double smallFlowEdgesThreshold);

        /// @brief Remove small triangles at the boundaries (removesmallflowlinks, part 2)
        /// @param[in] minFractionalAreaTriangles Small triangles at the boundaries will be eliminated.
        /// This threshold is the ration of the face area to the average area of neighboring faces.
        void RemoveSmallTrianglesAtBoundaries(double minFractionalAreaTriangles);

        /// @brief Computes m_nodesNodes, see class members
        void ComputeNodeNeighbours();

        /// @brief Get the orthogonality values, the inner product of edges and segments connecting the face circumcenters
        /// @param[out] orthogonality The edge orthogonality, passed to the client
        void GetOrthogonality(double* orthogonality);

        /// @brief Gets the smoothness values, ratios of the face areas
        /// @param[out] smoothness The smoothness at the edges
        void GetSmoothness(double* smoothness);

        /// @brief Gets the aspect ratios, the ratio edges to segments connecting the face circumcenters lengths
        /// @param aspectRatio The aspect ratios
        void GetAspectRatios(std::vector<double>& aspectRatios);

        ///  @brief Classifies the nodes (makenetnodescoding)
        void ClassifyNodes();

        /// @brief Sort edges in conterclockwise orther (Sort_links_ccw)
        /// @param[in] nodeIndex The node index for which sorting should take place
        void SortEdgesInCounterClockWiseOrder(int nodeIndex);

        /// @brief Remove coinciding triangles
        void RemoveDegeneratedTriangles();

        /// @brief Transform non-triangular faces in triangular faces
        void TriangulateFaces();

        /// @brief Make a dual face around the node, enlarged by a factor
        /// @param nodeIndex
        /// @return
        bool MakeDualFace(int node, double enlargmentFactor, std::vector<Point>& dualFace);

        /// @brief Sorts the faces around a node, sorted in counter clock wise order
        /// @param[in] nodeIndex The node index
        /// @return The face indexses
        [[nodiscard]] std::vector<int> SortedFacesAroundNode(int node) const;

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
        std::vector<double> m_edgeLengths;          // The edge lengths
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
        /// @brief Node administration (setnodadmin)
        void NodeAdministration();

        /// @brief Find cells recursive, works with an arbitrary number of edges
        /// @param startingNode The starting node
        /// @param node The current node
        /// @param numEdges The number of edges visited so far
        /// @param previousEdge The previously visited edge
        /// @param edges The vector storing the current edges forming a face
        /// @param nodes The vector storing the current nodes forming a face
        /// @param sortedEdges The caching array used for sorting the edges, used to inquire if an edge has been already visited
        /// @param sortedNodes The caching array used for sorting the nodes, used to inquire if a node has been already visited
        void FindFacesRecursive(int startingNode,
                                int node,
                                int numEdges,
                                int previousEdge,
                                std::vector<int>& edges,
                                std::vector<int>& nodes,
                                std::vector<int>& sortedEdges,
                                std::vector<int>& sortedNodes,
                                std::vector<Point>& nodalValues);

        /// @brief Checks if a triangle has an acute angle (checktriangle)
        /// @param[in] faceNodes
        /// @param[in] nodes
        /// @returns If triangle is okay
        [[nodiscard]] bool CheckTriangle(const std::vector<int>& faceNodes, const std::vector<Point>& nodes) const;

        /// @brief Removes all invalid nodes and edges
        void RemoveInvalidNodesAndEdges();

        int m_numFaces = 0;               // number of valid faces (nump)
        int m_numNodes = 0;               // Number of valid nodes in m_nodes
        int m_numEdges = 0;               // Number of valid edges in m_edges
        std::vector<double> m_edgeAngles; // internal cache for sorting the edges around nodes

        bool m_nodesRTreeRequiresUpdate = false; //m_nodesRTree requires an update
        bool m_edgesRTreeRequiresUpdate = false; //m_edgesRTree requires an update
    };
} // namespace meshkernel
