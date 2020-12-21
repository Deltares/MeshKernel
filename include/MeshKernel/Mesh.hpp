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

#include <MeshKernel/Entities.hpp>
#include <MeshKernel/MakeMeshParameters.hpp>
#include <MeshKernel/SpatialTrees.hpp>

/// \namespace meshkernel
/// @brief Contains the logic of the C++ static library
namespace meshkernel
{
    // Forward declarations
    class CurvilinearGrid;
    class Polygons;
    class MakeMeshParameters;
    class GeometryList;

    /// @brief A class describing an unstructured mesh
    class Mesh
    {
    public:
        /// Enumerator describing the different options to delete a mesh
        enum DeleteMeshOptions
        {
            AllNodesInside = 0,
            FacesWithIncludedCircumcenters = 1,
            FacesCompletelyIncluded = 2
        };

        /// Enumerator describing the different options to administrate a mesh
        enum class AdministrationOptions
        {
            AdministrateMeshEdges,
            AdministrateMeshEdgesAndFaces
        };

        /// Enumerator describing the different node types
        enum class NodeTypes
        {
            internalNode,
            onRing,
            cornerNode,
            hangingNode,
            other
        };

        /// @brief Default constructor
        /// @returns
        Mesh() = default;

        /// @brief Converting constructor, from curvilinear grid to mesh (gridtonet)
        /// @param curvilinearGrid The curvilinear grid to create the mesh from
        /// @param projection The \ref Projection to use
        /// @returns
        Mesh(const CurvilinearGrid& curvilinearGrid, Projection projection);

        /// @brief Create triangular grid from nodes (triangulatesamplestonetwork)
        /// @param[in] nodes Input nodes
        /// @param[in] polygons Selection polygon
        /// @param[in] projection Projection to use
        Mesh(const std::vector<Point>& nodes, const Polygons& polygons, Projection projection);

        /// @brief Construct the mesh starting from the edges and nodes
        /// @param[in] edges The input edges
        /// @param[in] nodes The input nodes
        /// @param[in] projection Projection to use
        /// @param[in] administration Type of administration to perform
        Mesh(const std::vector<Edge>& edges, const std::vector<Point>& nodes, Projection projection, AdministrationOptions administration = AdministrationOptions::AdministrateMeshEdgesAndFaces);

        /// @brief Add meshes: result is a mesh composed of the additions
        /// firstMesh += secondmesh results in the second mesh being added to the first
        /// @param[in] rhs The mesh to add
        /// @returns The resulting mesh
        Mesh& operator+=(Mesh const& rhs);

        /// @brief Set internal flat copies of nodes and edges, so the pointer to the first entry is communicated with the front-end
        /// @param administrationOption Type of administration to perform
        void SetFlatCopies(AdministrationOptions administrationOption);

        /// @brief Perform mesh administration
        /// @param administrationOption Type of administration to perform
        void Administrate(AdministrationOptions administrationOption);

        /// @brief Compute face circumcenters
        void ComputeFaceCircumcentersMassCentersAndAreas(bool computeMassCenters = false);

        /// @brief Find faces: constructs the m_facesNodes mapping, face mass centers and areas (findcells)
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
        /// @param[in] MakeMeshParameters The structure containing the make grid parameters
        /// @param[in] polygons The polygon to account for
        void MakeMesh(const meshkernelapi::MakeMeshParameters& MakeMeshParameters, const Polygons& polygons);

        /// @brief Deletes a mesh in a polygon, using several options (delnet)
        /// @param[in] polygons The polygon where to perform the operation
        /// @param[in] deletionOption The deletion option
        /// @param[in] invertDeletion Inverts the selected node to delete (instead of outside the polygon, inside the polygon)
        void DeleteMesh(const Polygons& polygons, int deletionOption, bool invertDeletion);

        /// @brief Connect two existing nodes, forming a new edge (connectdbn)
        /// @param[in] startNode The start node index
        /// @param[in] endNode The end node index
        /// @return The index of the new edge
        size_t ConnectNodes(int startNode, int endNode);

        /// @brief Insert a new node in the mesh (setnewpoint)
        /// @param[in] newPoint The coordinate of the new point
        /// @return The index of the new node
        size_t InsertNode(const Point& newPoint);

        /// @brief Delete a node
        /// @param[in] nodeIndex The index of the node to delete
        void DeleteNode(int nodeIndex);

        /// @brief Find the edge sharing two nodes
        /// @param[in] firstNodeIndex The index of the first node
        /// @param[in] secondNodeIndex The index of the second node
        /// @return The edge index
        int FindEdge(int firstNodeIndex, int secondNodeIndex) const;

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

        /// @brief For a face create a closed polygon and fill local mapping caches (get_cellpolygon)
        /// @param[in] faceIndex The face index
        /// @param[out] polygonNodesCache The node cache array filled with the nodes values
        /// @param[out] localNodeIndicesCache The consecutive node index in polygonNodesCache (0, 1, 2,...)
        /// @param[out] globalEdgeIndicesCache The edge cache array filled with the global edge indices
        void ComputeFaceClosedPolygonWithLocalMappings(int faceIndex,
                                                       std::vector<Point>& polygonNodesCache,
                                                       std::vector<int>& localNodeIndicesCache,
                                                       std::vector<int>& globalEdgeIndicesCache) const;

        /// @brief For a face create a closed polygon
        /// @param[in] faceIndex The face index
        /// @param[in,out] polygonNodesCache The cache array to be filled with the nodes values
        void ComputeFaceClosedPolygon(int faceIndex, std::vector<Point>& polygonNodesCache) const;

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
        /// @return The shared node (-1 if no node is found)
        [[nodiscard]] int FindCommonNode(int firstEdgeIndex, int secondEdgeIndex) const;

        /// @brief Compute the lengths of all edges in one go
        void ComputeEdgesLengths();

        /// @brief Computes the edges centers  in one go
        void ComputeEdgesCenters();

        /// @brief Get the number of valid nodes
        /// @return The number of valid node
        [[nodiscard]] auto GetNumNodes() const { return m_numNodes; }

        /// @brief Get the number of valid edges
        /// @return The number of valid edges
        [[nodiscard]] auto GetNumEdges() const { return m_numEdges; }

        /// @brief Get the number of valid faces
        /// @return The number of valid faces
        [[nodiscard]] auto GetNumFaces() const { return m_numFaces; }

        /// @brief Get the number of edges for a face
        /// @param[in] faceIndex The face index
        /// @return The number of valid faces
        [[nodiscard]] auto GetNumFaceEdges(size_t faceIndex) const { return m_numFacesNodes[faceIndex]; }

        /// @brief Get the number of faces an edges shares
        /// @param[in] edgeIndex The edge index
        /// @return The number of faces an edges shares
        [[nodiscard]] auto GetNumEdgesFaces(size_t edgeIndex) const { return m_edgesNumFaces[edgeIndex]; }

        /// @brief Inquire if an edge is on boundary
        /// @param edge The edge index
        /// @return If the edge is on boundary
        [[nodiscard]] bool IsEdgeOnBoundary(size_t edge) const { return m_edgesNumFaces[edge] == 1; }

        /// @brief Inquire if a face is on boundary
        /// @param face The face index
        /// @return If the face is on boundary
        [[nodiscard]] bool IsFaceOnBoundary(size_t face) const;

        /// @brief For a closed polygon, compute the circumcenter of a face (getcircumcenter)
        /// @param[in,out] polygon Cache storing the face nodes
        /// @param[in] edgesNumFaces For meshes, the number of faces sharing the edges
        /// @returns The computed circumcenter
        [[nodiscard]] Point ComputeFaceCircumenter(std::vector<Point>& polygon,
                                                   const std::vector<size_t>& edgesNumFaces) const;

        /// @brief Gets the mass centers of obtuse triangles
        /// @returns The center of obtuse triangles
        [[nodiscard]] std::vector<Point> GetObtuseTrianglesCenters();

        /// @brief Gets the edges crossing the small flow edges
        /// @param[in] smallFlowEdgesThreshold The configurable threshold for detecting the small flow edges
        /// @returns The indices of the edges crossing small flow edges
        [[nodiscard]] std::vector<int> GetEdgesCrossingSmallFlowEdges(double smallFlowEdgesThreshold);

        /// @brief Gets the flow edges centers from the crossing edges
        /// @param[in] edges The crossing edges indices
        /// @returns The centers of the flow edges
        [[nodiscard]] std::vector<Point> GetFlowEdgesCenters(const std::vector<int>& edges) const;

        /// @brief Deletes small flow edges (removesmallflowlinks, part 1)
        /// @param[in] smallFlowEdgesThreshold The configurable threshold for detecting the small flow edges
        void DeleteSmallFlowEdges(double smallFlowEdgesThreshold);

        /// @brief Deletes small triangles at the boundaries (removesmallflowlinks, part 2)
        /// @param[in] minFractionalAreaTriangles Small triangles at the boundaries will be eliminated.
        /// This threshold is the ration of the face area to the average area of neighboring faces.
        void DeleteSmallTrianglesAtBoundaries(double minFractionalAreaTriangles);

        /// @brief Computes m_nodesNodes, see class members
        void ComputeNodeNeighbours();

        /// @brief Get the orthogonality values, the inner product of edges and segments connecting the face circumcenters
        /// @return The edge orthogonality
        [[nodiscard]] std::vector<double> GetOrthogonality();

        /// @brief Gets the smoothness values, ratios of the face areas
        /// @return The smoothness at the edges
        [[nodiscard]] std::vector<double> GetSmoothness();

        /// @brief Gets the aspect ratios (the ratios edges lengths to flow edges lengths)
        /// @param[in,out] aspectRatios The aspect ratios (passed as reference to avoid re-allocation)
        void ComputeAspectRatios(std::vector<double>& aspectRatios);

        ///  @brief Classifies the nodes (makenetnodescoding)
        void ClassifyNodes();

        /// @brief Sort edges in conterclockwise orther (Sort_links_ccw)
        /// @param[in] nodeIndex The node index for which sorting should take place
        void SortEdgesInCounterClockWiseOrder(int nodeIndex);

        /// @brief Deletes coinciding triangles
        void DeleteDegeneratedTriangles();

        /// @brief Transform non-triangular faces in triangular faces
        void TriangulateFaces();

        /// @brief Make a dual face around the node, enlarged by a factor
        /// @param[in] node The node index
        /// @param[in] enlargementFactor The factor by which the dual face is enlarged
        /// @param[out] dualFace The dual face to be calculated
        void MakeDualFace(int node, double enlargementFactor, std::vector<Point>& dualFace);

        /// @brief Sorts the faces around a node, sorted in counter clock wise order
        /// @param[in] node The node index
        /// @return The face indexses
        [[nodiscard]] std::vector<int> SortedFacesAroundNode(int node) const;

        /// @brief Convert all mesh boundaries to a vector of polygon nodes, including holes (copynetboundstopol)
        /// @return The resulting polygon mesh boundary
        [[nodiscard]] std::vector<Point> MeshBoundaryToPolygon(const std::vector<Point>& polygonNodes);

        /// @brief Constructs a polygon from the meshboundary, by walking through the mesh
        /// @param[in] polygonNodes The input mesh
        /// @param[in] isVisited the visited mesh nodes
        /// @param[in] currentNode the current node
        /// @param[out] meshBoundaryPolygon The resulting polygon points
        void WalkBoundaryFromNode(const std::vector<Point>& polygonNodes,
                                  std::vector<bool>& isVisited,
                                  int& currentNode,
                                  std::vector<Point>& meshBoundaryPolygon) const;

        /// @brief Gets the hanging edges
        /// @return A vector with the indices of the hanging edges
        std::vector<size_t> GetHangingEdges() const;

        /// @brief Deletes the hanging edges
        void DeleteHangingEdges();

        // nodes
        std::vector<Point> m_nodes;                    ///< The mesh nodes (xk, yk)
        std::vector<std::vector<size_t>> m_nodesEdges; ///< For each node, the indices of connected edges (nod%lin)
        std::vector<int> m_nodesNumEdges;              ///< For each node, the number of connected edges (nmk)
        std::vector<int> m_nodeMask;                   ///< The node mask (kc)
        std::vector<std::vector<int>> m_nodesNodes;    ///< For each node, its neighbours
        std::vector<int> m_nodesTypes;                 ///< The node types (nb)

        // edges
        std::vector<Edge> m_edges;                  ///< The edges, defined as first and second node(kn)
        std::vector<std::vector<int>> m_edgesFaces; ///< For each edge, the shared face index (lne)
        std::vector<int> m_edgesNumFaces;           ///< For each edge, the number of shared faces(lnn)
        std::vector<double> m_edgeLengths;          ///< The edge lengths
        std::vector<int> m_edgeMask;                ///< The edge mask (lc)
        std::vector<Point> m_edgesCenters;          ///< The edges centers

        // faces
        std::vector<std::vector<int>> m_facesNodes; ///< The nodes composing the faces, in ccw order (netcell%Nod)
        std::vector<int> m_numFacesNodes;           ///< The number of nodes composing the face (netcell%N)
        std::vector<std::vector<int>> m_facesEdges; ///< The edge indices composing the face (netcell%lin)
        std::vector<Point> m_facesCircumcenters;    ///< The face circumcenters the face circumcenter (xz, yz)
        std::vector<Point> m_facesMassCenters;      ///< The faces centers of mass (xzw, yzw)
        std::vector<double> m_faceArea;             ///< The face area
        std::vector<Point> m_polygonNodesCache;     ///< Cache to store the face nodes

        // vectors for communicating with the client
        std::vector<double> m_nodex;               ///< The nodes x-coordinate
        std::vector<double> m_nodey;               ///< The nodes y-coordinate
        std::vector<double> m_nodez;               ///< The nodes z-coordinate
        std::vector<int> m_edgeNodes;              ///< For each edge, the nodes
        std::vector<int> m_faceNodes;              ///< For each face, the nodes
        std::vector<double> m_facesCircumcentersx; ///< The circumcenters x-coordinate
        std::vector<double> m_facesCircumcentersy; ///< The circumcenters y-coordinate
        std::vector<double> m_facesCircumcentersz; ///< The circumcenters z-coordinate

        Projection m_projection; ///< The projection used

        SpatialTrees::RTree m_nodesRTree; ///< Spatial R-Tree used to inquire node nodes
        SpatialTrees::RTree m_edgesRTree; ///< Spatial R-Tree used to inquire edges centers

        int m_maxNumNeighbours = 0; ///< Maximum number of neighbors

    private:
        /// @brief Node administration (setnodadmin)
        void NodeAdministration();

        /// @brief Find cells recursive, works with an arbitrary number of edges
        /// @param startingNode The starting node
        /// @param node The current node
        /// @param numEdges The number of edges visited so far
        /// @param previousEdge The previously visited edge
        /// @param numClosingEdges The number of edges closing a face (3 for triangles, 4 for quads, etc)
        /// @param edges The vector storing the current edges forming a face
        /// @param nodes The vector storing the current nodes forming a face
        /// @param sortedEdges The caching array used for sorting the edges, used to inquire if an edge has been already visited
        /// @param sortedNodes The caching array used for sorting the nodes, used to inquire if a node has been already visited
        void FindFacesRecursive(int startingNode,
                                int node,
                                int numEdges,
                                size_t previousEdge,
                                size_t numClosingEdges,
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
        void DeleteInvalidNodesAndEdges();

        int m_numFaces = 0;               ///< Number of valid faces (nump)
        int m_numNodes = 0;               ///< Number of valid nodes in m_nodes
        int m_numEdges = 0;               ///< Number of valid edges in m_edges
        std::vector<double> m_edgeAngles; ///< Internal cache for sorting the edges around nodes

        bool m_nodesRTreeRequiresUpdate = true; ///< m_nodesRTree requires an update
        bool m_edgesRTreeRequiresUpdate = true; ///< m_edgesRTree requires an update
    };
} // namespace meshkernel
