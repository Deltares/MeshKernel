//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2021.
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

#include <MeshKernel/Entities.hpp>
#include <MeshKernel/RTree.hpp>

/// \namespace meshkernel
/// @brief Contains the logic of the C++ static library
namespace meshkernel
{
    class Polygons;

    /// @brief A class describing an unstructured mesh.
    /// This class contains the shared functionality between 1d or 2d meshes.
    ///
    /// MeshKernel can handle 2d meshes and 1d meshes.
    /// Algorithms require certain mappings to be available for both Mesh1D and Mesh2D,
    /// such as a mapping listing all edge indices connected to a particular node.
    /// The methods computing these mappings are shared between Mesh2D and Mesh1D, and implemented in the Mesh base class.
    /// The Mesh base class also contains other common data members,
    /// such as the node coordinate, the edges definitions, the face definitions and the mesh projection.
    /// The Mesh base class has the following responsibilities:
    ///
    /// -   Construct the mesh faces from the nodes and edges and other mesh
    ///     mappings required by all algorithms (Mesh::FindFaces).
    ///     Mesh::FindFaces is using recursion to find faces with up to 6 edges (meshkernel::Mesh::m_maximumNumberOfEdgesPerFace).
    ///
    /// -   Supporting mesh editing, namely:
    ///
    ///     -   Node merging
    ///
    ///     -   Node insertion
    ///
    ///     -   Moving a node
    ///
    ///     -   Inserting edges
    ///
    ///     -   Deleting edges
    ///
    ///     -   Merging nodes (merging two nodes at meshkernel::mergingDistance). This algorithm use an r-tree for inquiring
    /// adjacent nodes, see later.
    ///
    /// -   Converting a curvilinear grid to an unstructured mesh (converting
    ///     constructor).
    ///
    /// -   Holding the mesh projection (cartesian, spherical, or spherical
    ///     accurate).
    ///
    /// -   Making a quad mesh from a polygon or from parameters.
    ///
    /// -   Making a triangular mesh from a polygon. This algorithm introduces a dependency on
    /// the Richard Shewchuk Triangle.c library, added as an external component in extern/triangle folder.
    ///
    /// The public interface of the mesh class contains several algorithms,
    /// which modify the mesh class members when they are called.
    class Mesh
    {
    public:
        /// @brief Enumerator describing the different mesh types
        enum class Type
        {
            Mesh1D, ///< Mesh1D
            Mesh2D  ///< Mesh2D
        };

        /// @enum Location
        /// @brief Mesh locations enumeration
        enum class Location
        {
            Faces = 0,  ///< Faces
            Nodes = 1,  ///< Nodes
            Edges = 2,  ///< Edges
            Unknown = 3 ///< Unknown
        };

        /// edge-segment intersection
        struct EdgeMeshPolylineIntersection
        {
            int polylineSegmentIndex{constants::missing::intValue};          ///< The intersected segment index (a polyline can formed by several segments)
            double polylineDistance{constants::missing::doubleValue};        ///< The location of the intersection expressed as distance from the polyline start
            double adimensionalPolylineSegmentDistance{constants::missing::doubleValue}; ///< The location of the intersection expressed as an adimensional distance from the segment start
            size_t edgeIndex{constants::missing::sizetValue};                ///< The first node of the edge is on the left (the virtual node)
            size_t edgeFirstNode{constants::missing::sizetValue};            ///< The first node of the edge is on the left (the virtual node)
            size_t edgeSecondNode{constants::missing::sizetValue};           ///< The second node of the edge is on the right (the inner node)
            double edgeDistance{constants::missing::doubleValue};            ///< The location of the intersection expressed as an adimensional distance from the edge start
        };

        /// face-segment intersection
        struct FaceMeshPolylineIntersection
        {
            double polylineDistance{constants::missing::doubleValue};        ///< The location of the intersection expressed as an adimensional distance from the polyline start
            size_t faceIndex{constants::missing::sizetValue};                ///< The face index
            std::vector<size_t> edgeIndexses;                                ///< The indexes of crossed edges
            std::vector<size_t> edgeNodes;                                   ///< The indexes of the nodes defining the crossed edges
        };

        /// @brief Default constructor
        Mesh() = default;

        /// @brief  Constructs an empty mesh, sets only the projection
        /// @param[in] projection  The projection to use
        Mesh(Projection projection) : m_projection(projection) {}

        /// @brief Construct a mesh starting from the edges and nodes
        /// @param[in] edges The input edges
        /// @param[in] nodes The input nodes
        /// @param[in] projection  The projection to use
        Mesh(const std::vector<Edge>& edges,
             const std::vector<Point>& nodes,
             Projection projection);

        /// @brief Inquire if a node is on boundary
        /// @param[in] node The node index
        /// @return If the node is on boundary
        [[nodiscard]] bool IsNodeOnBoundary(size_t node) const { return m_nodesNumEdges[node] == 1; }

        /// @brief Get the number of valid nodes
        /// @return The number of valid node
        [[nodiscard]] auto GetNumNodes() const { return m_nodes.size(); }

        /// @brief Get the number of valid edges
        /// @return The number of valid edges
        [[nodiscard]] auto GetNumEdges() const { return m_edges.size(); }

        /// @brief Get the number of valid faces
        /// @return The number of valid faces
        [[nodiscard]] auto GetNumFaces() const { return m_facesNodes.size(); }

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
        /// @param[in] face The face index
        /// @return If the face is on boundary
        [[nodiscard]] bool IsFaceOnBoundary(size_t face) const;

        /// @brief Merges two mesh nodes
        /// @param[in] startNode The index of the first node to be merged
        /// @param[in] endNode The second of the second node to be merged
        void MergeTwoNodes(size_t startNode, size_t endNode);

        /// @brief Merge close mesh nodes inside a polygon (MERGENODESINPOLYGON)
        /// @param[in] polygons Polygon where to perform the merging
        /// @param[in] mergingDistance The distance below which two nodes will be merged
        void MergeNodesInPolygon(const Polygons& polygons, double mergingDistance);

        /// @brief Connect two existing nodes, forming a new edge (connectdbn)
        /// @param[in] startNode The start node index
        /// @param[in] endNode The end node index
        /// @return The index of the new edge
        size_t ConnectNodes(size_t startNode, size_t endNode);

        /// @brief Insert a new node in the mesh (setnewpoint)
        /// @param[in] newPoint The coordinate of the new point
        /// @return The index of the new node
        size_t InsertNode(const Point& newPoint);

        /// @brief Delete a node
        /// @param[in] node The index of the node to delete
        void DeleteNode(size_t node);

        /// @brief Find the edge sharing two nodes
        /// @param[in] firstNodeIndex The index of the first node
        /// @param[in] secondNodeIndex The index of the second node
        /// @return The edge index
        [[nodiscard]] size_t FindEdge(size_t firstNodeIndex, size_t secondNodeIndex) const;

        /// @brief Move a node to a new location
        /// @param[in] newPoint The new location
        /// @param[in] nodeindex The index of the node to move
        void MoveNode(Point newPoint, size_t nodeindex);

        /// @brief Get the index of a node close to a point
        /// @param[in] point The starting point from where to start the search
        /// @param[in] nodeMask The mask to apply to mesh nodes, if the mask value is false, the next closest node will be considered
        /// @returns The index of the closest node
        [[nodiscard]] size_t FindNodeCloseToAPoint(Point point, const std::vector<bool>& nodeMask);

        /// @brief Get the index of a node close to a point
        /// @param[in] point The starting point from where to start the search
        /// @param[in] searchRadius The search radius
        /// @returns The index of the closest node
        [[nodiscard]] size_t FindNodeCloseToAPoint(Point const& point, double searchRadius);

        /// @brief Deletes an edge
        /// @param[in] edge The edge index
        void DeleteEdge(size_t edge);

        /// Finds the closest edge close to a point
        /// @param[in] point The starting point from where to start the search
        /// @returns The index of the closest edge
        [[nodiscard]] size_t FindEdgeCloseToAPoint(Point point);

        /// @brief Find the common node two edges share
        /// This method uses return parameters since the success is evaluated in a hot loop
        /// @param[in] firstEdgeIndex The index of the first edge
        /// @param[in] secondEdgeIndex The index of the second edge
        /// @return The shared node (constants::missing::sizetValue if no node is found)
        [[nodiscard]] size_t FindCommonNode(size_t firstEdgeIndex, size_t secondEdgeIndex) const;

        /// @brief Compute the lengths of all edges in one go
        void ComputeEdgesLengths();

        /// @brief Computes the edges centers  in one go
        void ComputeEdgesCenters();

        /// @brief Node administration (setnodadmin)
        void NodeAdministration();

        /// @brief Removes all invalid nodes and edges
        void DeleteInvalidNodesAndEdges();

        /// @brief Perform node and edges administration
        void AdministrateNodesEdges();

        /// @brief Sort mesh edges around a node in counterclockwise order (Sort_links_ccw)
        /// @param[in] startNode The first node index where to perform edge sorting.
        /// @param[in] endNode   The last node index where to perform edge sorting.
        void SortEdgesInCounterClockWiseOrder(size_t startNode, size_t endNode);

        /// @brief Compute the max length of the edges connected to a node
        /// @param node The mesh node
        /// @return The max edge length
        double ComputeMaxLengthSurroundingEdges(size_t node);

        /// @brief Build the rtree for the corresponding location
        /// @param[in] meshLocation The mesh location for which the RTree is build
        void BuildTree(Location meshLocation);

        /// @brief Search all points sorted by proximity to another point.
        /// @param[in] point The reference point.
        /// @param[in] meshLocation The mesh location (e.g. nodes, edge centers or face circumcenters).
        void SearchNearestLocation(Point point, Location meshLocation);

        /// @brief Search the nearest point within a radius to another point.
        /// @param[in] point The reference point.
        /// @param[in] squaredRadius the squared value of the radius.
        /// @param[in] meshLocation The mesh location (e.g. nodes, edge centers or face circumcenters).
        void SearchNearestLocation(Point point, double squaredRadius, Location meshLocation);

        /// @brief Search the nearest points within a radius to another point.
        /// @param[in] point The reference point.
        /// @param[in] squaredRadius the squared value of the radius.
        /// @param[in] meshLocation The mesh location (e.g. nodes, edge centers or face circumcenters).
        void SearchLocations(Point point, double squaredRadius, Location meshLocation);

        /// @brief Gets the search results.
        /// To be used after \ref SearchLocations or \ref SearchNearestLocation.
        ///
        /// @param[in] meshLocation The mesh location (e.g. nodes, edge centers or face circumcenters).
        /// @return The number of found neighbors.
        size_t GetNumLocations(Location meshLocation) const;

        /// @brief Gets the index of the location, sorted by proximity. To be used after SearchNearestLocation or SearchNearestLocation.
        /// @param[in] index The closest neighbor index (index 0 corresponds to the closest).
        /// @param[in] meshLocation The mesh location (e.g. nodes, edge centers or face circumcenters).
        /// @return The index of the closest location.
        [[nodiscard]] size_t GetLocationsIndices(size_t index, Mesh::Location meshLocation);

        /// @brief Computes a vector with the mesh locations coordinates (nodes, edges or faces coordinates).
        ///
        /// @param[in] location The mesh location (e.g. nodes, edge centers or face circumcenters).
        /// @return The vector with the mesh locations.
        [[nodiscard]] std::vector<Point> ComputeLocations(Mesh::Location location) const;

        /// @brief Add meshes: result is a mesh composed of the additions
        /// firstMesh += secondmesh results in the second mesh being added to firstMesh
        /// @param[in] rhs The mesh to add
        /// @returns The resulting mesh
        Mesh& operator+=(Mesh const& rhs);

        // nodes
        std::vector<Point> m_nodes;                    ///< The mesh nodes (xk, yk)
        std::vector<std::vector<size_t>> m_nodesEdges; ///< For each node, the indices of connected edges (nod%lin)
        std::vector<size_t> m_nodesNumEdges;           ///< For each node, the number of connected edges (nmk)
        std::vector<std::vector<size_t>> m_nodesNodes; ///< For each node, its neighbors
        std::vector<int> m_nodesTypes;                 ///< The node types (nb)

        // edges
        std::vector<Edge> m_edges;                     ///< The edges, defined as first and second node(kn)
        std::vector<std::vector<size_t>> m_edgesFaces; ///< For each edge, the shared face index (lne)
        std::vector<size_t> m_edgesNumFaces;           ///< For each edge, the number of shared faces(lnn)
        std::vector<double> m_edgeLengths;             ///< The edge lengths
        std::vector<Point> m_edgesCenters;             ///< The edges centers

        // faces
        std::vector<std::vector<size_t>> m_facesNodes; ///< The nodes composing the faces, in ccw order (netcell%Nod)
        std::vector<size_t> m_numFacesNodes;           ///< The number of nodes composing the face (netcell%N)
        std::vector<std::vector<size_t>> m_facesEdges; ///< The edge indices composing the face (netcell%lin)
        std::vector<Point> m_facesCircumcenters;       ///< The face circumcenters the face circumcenter (xz, yz)
        std::vector<Point> m_facesMassCenters;         ///< The faces centers of mass (xzw, yzw)
        std::vector<double> m_faceArea;                ///< The face area

        Projection m_projection; ///< The projection used

        // counters
        bool m_nodesRTreeRequiresUpdate = true; ///< m_nodesRTree requires an update
        bool m_edgesRTreeRequiresUpdate = true; ///< m_edgesRTree requires an update
        RTree m_nodesRTree;                     ///< Spatial R-Tree used to inquire node nodes
        RTree m_edgesRTree;                     ///< Spatial R-Tree used to inquire edges centers
        RTree m_facesRTree;                     ///< Spatial R-Tree used to inquire face circumcenters

        // constants
        static constexpr size_t m_maximumNumberOfEdgesPerNode = 12;                                  ///< Maximum number of edges per node
        static constexpr size_t m_maximumNumberOfEdgesPerFace = 6;                                   ///< Maximum number of edges per face
        static constexpr size_t m_maximumNumberOfNodesPerFace = 8;                                   ///< Maximum number of nodes per face
        static constexpr size_t m_maximumNumberOfConnectedNodes = m_maximumNumberOfEdgesPerNode * 4; ///< Maximum number of connected nodes
        static constexpr size_t m_numNodesQuads = 4;                                                 ///< Number of nodes in a quadrilateral
        static constexpr size_t m_numNodesInTriangle = 3;                                            ///< Number of nodes in a triangle

    private:
        static double constexpr m_minimumDeltaCoordinate = 1e-14; ///< Minimum delta coordinate
    };
} // namespace meshkernel
