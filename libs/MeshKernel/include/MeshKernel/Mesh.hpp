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
#include <memory>

#include "MeshKernel/BoundingBox.hpp"
#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/UndoActions/AddEdgeAction.hpp"
#include "MeshKernel/UndoActions/AddNodeAction.hpp"
#include "MeshKernel/UndoActions/CompoundUndoAction.hpp"
#include "MeshKernel/UndoActions/DeleteEdgeAction.hpp"
#include "MeshKernel/UndoActions/DeleteNodeAction.hpp"
#include "MeshKernel/UndoActions/FullUnstructuredGridUndo.hpp"
#include "MeshKernel/UndoActions/MeshConversionAction.hpp"
#include "MeshKernel/UndoActions/NodeTranslationAction.hpp"
#include "MeshKernel/UndoActions/ResetEdgeAction.hpp"
#include "MeshKernel/UndoActions/ResetNodeAction.hpp"
#include "MeshKernel/UndoActions/UndoAction.hpp"
#include "Utilities/RTreeBase.hpp"

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
    ///
    class Mesh
    {
    public:
        /// @brief Enumerator describing the different mesh types
        enum class Type
        {
            Mesh1D, ///< Mesh1D
            Mesh2D  ///< Mesh2D
        };

        /// @brief Define virtual destructor
        virtual ~Mesh() = default;

        /// @brief Default constructor, setting a cartesian projection
        Mesh();

        /// @brief Delete assignment operator
        Mesh& operator=(const Mesh& mesh) = delete;

        /// @brief Delete move assignment operator
        Mesh& operator=(Mesh&& mesh) = delete;

        /// @brief Copy constructor taking only a mesh
        Mesh(const Mesh& mesh) = delete;

        /// @brief Move constructor taking only a mesh
        Mesh(Mesh&& mesh) = delete;

        /// @brief  Constructs an empty mesh, sets only the projection
        /// @param[in] projection  The projection to use
        explicit Mesh(Projection projection);

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
        [[nodiscard]] bool IsNodeOnBoundary(UInt node) const { return m_nodesNumEdges[node] == 1; }

        /// @brief Get the number of valid nodes
        /// @return The number of valid node
        [[nodiscard]] auto GetNumNodes() const { return static_cast<UInt>(m_nodes.size()); }

        /// @brief Get the number of valid edges
        /// @return The number of valid edges
        [[nodiscard]] auto GetNumEdges() const { return static_cast<UInt>(m_edges.size()); }

        /// @brief Get the number of valid faces
        /// @return The number of valid faces
        [[nodiscard]] auto GetNumFaces() const { return static_cast<UInt>(m_facesNodes.size()); }

        /// @brief Get the number of valid nodes
        /// @return The number of valid nodes
        [[nodiscard]] UInt GetNumValidNodes() const;

        /// @brief Get the number of valid edges
        /// @return The number of valid edges
        [[nodiscard]] UInt GetNumValidEdges() const;

        /// @brief Get the number of edges for a face
        /// @param[in] faceIndex The face index
        /// @return The number of valid faces
        [[nodiscard]] auto GetNumFaceEdges(UInt faceIndex) const { return m_numFacesNodes[faceIndex]; }

        /// @brief Get the number of faces an edges shares
        /// @param[in] edgeIndex The edge index
        /// @return The number of faces an edges shares
        [[nodiscard]] auto GetNumEdgesFaces(UInt edgeIndex) const { return m_edgesNumFaces[edgeIndex]; }

        /// @brief Get the local edge number for an element edge.
        // TODO add unit test and with all failing cases
        [[nodiscard]] UInt GetEdgeIndex(const UInt elementId, const UInt edgeId) const;

        /// @brief Get the local node number for an element node.
        // TODO add unit test and with all failing cases
        [[nodiscard]] UInt GetNodeIndex(const UInt elementId, const UInt nodeId) const;

        /// @brief Inquire if an edge is on boundary
        /// @param edge The edge index
        /// @return If the edge is on boundary
        [[nodiscard]] bool IsEdgeOnBoundary(UInt edge) const { return m_edgesNumFaces[edge] == 1; }

        /// @brief Inquire if a face is on boundary
        /// @param[in] face The face index
        /// @return If the face is on boundary
        [[nodiscard]] bool IsFaceOnBoundary(UInt face) const;

        /// @brief Get vector of all nodes
        // TODO Can this be removed?
        const std::vector<Point>& Nodes() const;

        /// @brief Get the node at the position
        const Point& Node(const UInt index) const;

        /// @brief Set all nodes to a new set of values.
        void SetNodes(const std::vector<Point>& newValues);

        /// @brief Set a node to a new value, bypassing the undo action.
        void SetNode(const UInt index, const Point& newValue);

        /// @brief Set the node to a new value, this value may be the in-valid value.
        [[nodiscard]] std::unique_ptr<ResetNodeAction> ResetNode(const UInt index, const Point& newValue);

        /// @brief Get constant reference to an edge
        const Edge& GetEdge(const UInt index) const;

        /// @brief Get a non-constant reference to an edge
        Edge& GetEdge(const UInt index);

        /// @brief Get all edges
        // TODO get rid of this function
        const std::vector<Edge>& Edges() const;

        /// @brief Set all edges to a new set of values.
        void SetEdges(const std::vector<Edge>& newValues);

        /// @brief Get the local index of the node belong to a face.
        ///
        /// If the node cannot be found the null value will be returned.
        UInt GetLocalFaceNodeIndex(const UInt faceIndex, const UInt nodeIndex) const;

        /// @brief Merges two mesh nodes
        /// @param[in] startNode The index of the first node to be merged
        /// @param[in] endNode The second of the second node to be merged
        [[nodiscard]] std::unique_ptr<UndoAction> MergeTwoNodes(UInt startNode, UInt endNode);

        /// @brief Merge close mesh nodes inside a polygon (MERGENODESINPOLYGON)
        /// @param[in] polygons Polygon where to perform the merging
        /// @param[in] mergingDistance The distance below which two nodes will be merged
        [[nodiscard]] std::unique_ptr<UndoAction> MergeNodesInPolygon(const Polygons& polygons, double mergingDistance);

        /// @brief Insert a new node in the mesh (setnewpoint)
        /// @param[in] newPoint The coordinate of the new point
        /// @return The index of the new node and the pointer to the undoAction
        [[nodiscard]] std::tuple<UInt, std::unique_ptr<AddNodeAction>> InsertNode(const Point& newPoint);

        /// @brief Connect two existing nodes, checking if the nodes are already connected.
        /// If the nodes are not connected a new edge is formed, otherwise UInt invalid value is returned. (connectdbn)
        /// @param[in] startNode The start node index
        /// @param[in] endNode The end node index
        /// @return The index of the new edge and the undoAction to connect two nodes
        [[nodiscard]] std::tuple<UInt, std::unique_ptr<AddEdgeAction>> ConnectNodes(UInt startNode, UInt endNode);

        /// @brief Change the nodes referenced by the edge.
        [[nodiscard]] std::unique_ptr<ResetEdgeAction> ResetEdge(UInt edgeId, const Edge& edge);

        /// @brief Deletes a node and removes any connected edges
        /// @param[in] node The node index
        /// @return The undoAction to delete the node and any connecting edges
        [[nodiscard]] std::unique_ptr<DeleteNodeAction> DeleteNode(UInt node);

        /// @brief Find the edge sharing two nodes
        /// @param[in] firstNodeIndex The index of the first node
        /// @param[in] secondNodeIndex The index of the second node
        /// @return The edge index
        [[nodiscard]] UInt FindEdge(UInt firstNodeIndex, UInt secondNodeIndex) const;

        /// @brief Find the edge using a linear search, without connectivity information (much slower than FindEdge)
        /// @param[in] firstNodeIndex The index of the first node
        /// @param[in] secondNodeIndex The index of the second node
        /// @return The edge index
        [[nodiscard]] UInt FindEdgeWithLinearSearch(UInt firstNodeIndex, UInt secondNodeIndex) const;

        /// @brief Move a node to a new location
        /// @param[in] newPoint The new location
        /// @param[in] nodeindex The index of the node to move
        [[nodiscard]] std::unique_ptr<UndoAction> MoveNode(Point newPoint, UInt nodeindex);

        /// @brief Get the index of a location (node/edge or face) close to a point
        /// @param[in] point The starting point from where to start the search
        /// @param[in] location The location
        /// @param[in] locationMask The mask to apply to each location
        /// @param[in] boundingBox The bounding box
        /// @returns The index of the closest node
        [[nodiscard]] UInt FindLocationIndex(Point point,
                                             Location location,
                                             const std::vector<bool>& locationMask = {},
                                             const BoundingBox& boundingBox = {});

        /// @brief Get the index of a node close to a point
        /// @param[in] point The starting point from where to start the search
        /// @param[in] searchRadius The search radius
        /// @returns The index of the closest node
        [[nodiscard]] UInt FindNodeCloseToAPoint(Point const& point, double searchRadius);

        /// @brief Deletes an edge
        /// @param[in] edge The edge index
        /// @return The undoAction to delete the edge
        [[nodiscard]] std::unique_ptr<DeleteEdgeAction> DeleteEdge(UInt edge);

        /// @brief Find the common node two edges share
        /// This method uses return parameters since the success is evaluated in a hot loop
        /// @param[in] firstEdgeIndex The index of the first edge
        /// @param[in] secondEdgeIndex The index of the second edge
        /// @return The shared node (constants::missing::sizetValue if no node is found)
        [[nodiscard]] UInt FindCommonNode(UInt firstEdgeIndex, UInt secondEdgeIndex) const;

        /// @brief Compute the lengths of all edges in one go
        void ComputeEdgesLengths();

        /// @brief Compute the minimum edge length of the edges included in the polygon.
        /// An edge is considered included if one of the two nodes is inside the polygon.
        /// @param[in] polygon The polygon for considering an edge included
        /// @return The minimum edge length
        [[nodiscard]] double ComputeMinEdgeLength(const Polygons& polygon) const;

        /// @brief Computes the edges centers  in one go
        void ComputeEdgesCenters();

        /// @brief Node administration (setnodadmin)
        /// @return An estimated indicator for a quadrilateral dominated mesh.
        bool NodeAdministration();

        /// @brief Removes all invalid nodes and edges
        void DeleteInvalidNodesAndEdges();

        /// @brief Perform complete administration
        virtual void Administrate(CompoundUndoAction* undoAction = nullptr);

        /// @brief Perform node and edges administration
        void AdministrateNodesEdges(CompoundUndoAction* undoAction = nullptr);

        /// @brief Sort mesh edges around a node in counterclockwise order (Sort_links_ccw)
        /// @param[in] startNode The first node index where to perform edge sorting.
        /// @param[in] endNode   The last node index where to perform edge sorting.
        void SortEdgesInCounterClockWiseOrder(UInt startNode, UInt endNode);

        /// @brief Compute the max length of the edges connected to a node
        /// @param node The mesh node
        /// @return The max edge length
        double ComputeMaxLengthSurroundingEdges(UInt node);

        /// @brief Build the rtree for the corresponding location, using only the locations inside the bounding box
        /// @param[in] location The mesh location for which the RTree is build
        /// @param[in] boundingBox The bounding box
        void BuildTree(Location location, const BoundingBox& boundingBox = {});

        /// @brief Computes a vector with the mesh locations coordinates (nodes, edges or faces coordinates).
        ///
        /// @param[in] location The mesh location (e.g. nodes, edge centers or face circumcenters).
        /// @return The vector with the mesh locations.
        [[nodiscard]] std::vector<Point> ComputeLocations(Location location) const;

        /// @brief Computes if a location is in polygon.
        /// @param[in] polygon The input polygon.
        /// @param[in] location The mesh location (e.g. nodes, edge centers or face circumcenters).
        /// @return A vector of booleans indicating if a location is in a polygon or not.
        [[nodiscard]] std::vector<bool> IsLocationInPolygon(const Polygons& polygon, Location location) const;

        /// @brief Add meshes: result is a mesh composed of the additions
        /// firstMesh += secondmesh results in the second mesh being added to firstMesh
        /// @param[in] rhs The mesh to add
        /// @returns The resulting mesh
        Mesh& operator+=(Mesh const& rhs);

        /// @brief Add meshes: result is a mesh composed of the additions
        /// firstMesh += secondmesh results in the second mesh being added to firstMesh
        /// @param[in] rhs The mesh to add
        /// @returns The undo action
        std::unique_ptr<UndoAction> Join(const Mesh& rhs);

        /// @brief Get the mapping/indexing from the node array mapped to valid nodes
        std::vector<UInt> GetValidNodeMapping() const;

        /// @brief Get the mapping/indexing from the edge array mapped to valid edges
        std::vector<UInt> GetValidEdgeMapping() const;

        /// @brief Indicate if the edge-id is a valid edge
        ///
        /// A valid edge satisfies four conditions:
        /// The start and end indices are not the null value, if neither is, then
        /// also the nodes indexed by the edge are valid nodes. If any of these conditions
        /// is false then the edge is in-valid.
        bool IsValidEdge(const UInt edgeId) const;

        /// @brief Apply the reset node action
        void CommitAction(const ResetNodeAction& undoAction);

        /// @brief Apply the add node action
        void CommitAction(const AddNodeAction& undoAction);

        /// @brief Apply the add edge action
        void CommitAction(const AddEdgeAction& undoAction);

        /// @brief Apply the reset edge action
        void CommitAction(const ResetEdgeAction& undoAction);

        /// @brief Apply the delete node action
        void CommitAction(const DeleteNodeAction& undoAction);

        /// @brief Apply the node translation action
        void CommitAction(NodeTranslationAction& undoAction);

        /// @brief Apply the node translation action
        void CommitAction(MeshConversionAction& undoAction);

        /// @brief Apply the delete edge action
        void CommitAction(const DeleteEdgeAction& undoAction);

        /// @brief Set the node and edge values.
        void CommitAction(FullUnstructuredGridUndo& undoAction);

        /// @brief Undo the reset node action
        ///
        /// Restore mesh to state before node was reset
        void RestoreAction(const ResetNodeAction& undoAction);

        /// @brief Undo the add node action
        ///
        /// Restore mesh to state before node was added
        void RestoreAction(const AddNodeAction& undoAction);

        /// @brief Undo the add edge action
        ///
        /// Restore mesh to state before edge was added
        void RestoreAction(const AddEdgeAction& undoAction);

        /// @brief Undo the reset edge action
        ///
        /// Restore mesh to state before edge was reset
        void RestoreAction(const ResetEdgeAction& undoAction);

        /// @brief Undo the delete node action
        ///
        /// Restore mesh to state before node was deleted
        void RestoreAction(const DeleteNodeAction& undoAction);

        /// @brief Undo the node translation action
        ///
        /// Restore mesh to state before node was translated
        void RestoreAction(NodeTranslationAction& undoAction);

        /// @brief Undo the node translation action
        ///
        /// Restore mesh to state before node was translated
        void RestoreAction(MeshConversionAction& undoAction);

        /// @brief Undo the delete edge action
        ///
        /// Restore mesh to state before edge was deleted
        void RestoreAction(const DeleteEdgeAction& undoAction);

        /// @brief Undo entire node and edge values
        ///
        /// Restore mesh to previous state.
        void RestoreAction(FullUnstructuredGridUndo& undoAction);

        /// @brief Get a reference to the RTree for a specific location
        RTreeBase& GetRTree(Location location) const { return *m_RTrees.at(location); }

        /// @brief Set the m_nodesRTreeRequiresUpdate flag
        /// @param[in] value The value of the flag
        void SetNodesRTreeRequiresUpdate(bool value) { m_nodesRTreeRequiresUpdate = value; }

        /// @brief Set the m_edgesRTreeRequiresUpdate flag
        /// @param[in] value The value of the flag
        void SetEdgesRTreeRequiresUpdate(bool value) { m_edgesRTreeRequiresUpdate = value; }

        /// @brief Set the m_facesRTreeRequiresUpdate flag
        /// @param[in] value The value of the flag
        void SetFacesRTreeRequiresUpdate(bool value) { m_facesRTreeRequiresUpdate = value; }

        // nodes
        std::vector<std::vector<UInt>> m_nodesEdges; ///< For each node, the indices of connected edges (nod%lin)
        std::vector<UInt> m_nodesNumEdges;           ///< For each node, the number of connected edges (nmk)
        std::vector<std::vector<UInt>> m_nodesNodes; ///< For each node, its neighbors
        std::vector<int> m_nodesTypes;               ///< The node types (nb)

        // edges
        std::vector<std::array<UInt, 2>> m_edgesFaces; ///< For each edge, the shared face index (lne)
        std::vector<UInt> m_edgesNumFaces;             ///< For each edge, the number of shared faces(lnn)
        std::vector<double> m_edgeLengths;             ///< The edge lengths
        std::vector<Point> m_edgesCenters;             ///< The edges centers

        // faces
        std::vector<std::vector<UInt>> m_facesNodes; ///< The nodes composing the faces, in ccw order (netcell%Nod)
        std::vector<UInt> m_numFacesNodes;           ///< The number of nodes composing the face (netcell%N)
        std::vector<std::vector<UInt>> m_facesEdges; ///< The edge indices composing the face (netcell%lin)
        std::vector<Point> m_facesCircumcenters;     ///< The face circumcenters the face circumcenter (xz, yz)
        std::vector<Point> m_facesMassCenters;       ///< The faces centers of mass (xzw, yzw)
        std::vector<double> m_faceArea;              ///< The face area

        Projection m_projection; ///< The projection used

        // constants
        static constexpr UInt m_maximumNumberOfEdgesPerNode = 16;                                  ///< Maximum number of edges per node
        static constexpr UInt m_maximumNumberOfEdgesPerFace = 6;                                   ///< Maximum number of edges per face
        static constexpr UInt m_maximumNumberOfNodesPerFace = 6;                                   ///< Maximum number of nodes per face
        static constexpr UInt m_maximumNumberOfConnectedNodes = m_maximumNumberOfEdgesPerNode * 4; ///< Maximum number of connected nodes

    protected:
        // Make private
        std::vector<Point> m_nodes; ///< The mesh nodes (xk, yk)
        std::vector<Edge> m_edges;  ///< The edges, defined as first and second node(kn)

    private:
        static double constexpr m_minimumDeltaCoordinate = 1e-14; ///< Minimum delta coordinate

        // RTrees
        bool m_nodesRTreeRequiresUpdate = true;                            ///< m_nodesRTree requires an update
        bool m_edgesRTreeRequiresUpdate = true;                            ///< m_edgesRTree requires an update
        bool m_facesRTreeRequiresUpdate = true;                            ///< m_facesRTree requires an update
        std::unordered_map<Location, std::unique_ptr<RTreeBase>> m_RTrees; ///< The RTrees to use
        BoundingBox m_boundingBoxCache;                                    ///< Caches the last bounding box used for selecting the locations

        /// @brief Set nodes and edges that are not connected to be invalid.
        void SetUnConnectedNodesAndEdgesToInvalid(CompoundUndoAction* undoAction);

        /// @brief Find all nodes that are connected to an edge.
        ///
        /// Also count the number of edges that have either invalid index values or
        /// reference invalid nodes
        void FindConnectedNodes(std::vector<bool>& connectedNodes,
                                UInt& numInvalidEdges) const;

        /// @brief Invalidate any not connected to any edge.
        void InvalidateUnConnectedNodes(const std::vector<bool>& connectedNodes,
                                        UInt& numInvalidNodes,
                                        CompoundUndoAction* undoAction = nullptr);
    };
} // namespace meshkernel

inline const std::vector<meshkernel::Point>& meshkernel::Mesh::Nodes() const
{
    return m_nodes;
}

inline const meshkernel::Point& meshkernel::Mesh::Node(const UInt index) const
{
    if (index >= GetNumNodes())
    {
        throw ConstraintError("The node index, {}, is not in range.", index);
    }

    return m_nodes[index];
}

inline void meshkernel::Mesh::SetNodes(const std::vector<Point>& newValues)
{
    m_nodes = newValues;
    m_nodesRTreeRequiresUpdate = true;
    m_edgesRTreeRequiresUpdate = true;
    m_facesRTreeRequiresUpdate = true;
}

inline const meshkernel::Edge& meshkernel::Mesh::GetEdge(const UInt index) const
{
    if (index >= GetNumEdges())
    {
        throw ConstraintError("The edge index, {}, is not in range.", index);
    }

    return m_edges[index];
}

inline meshkernel::Edge& meshkernel::Mesh::GetEdge(const UInt index)
{
    if (index >= GetNumEdges())
    {
        throw ConstraintError("The edge index, {}, is not in range.", index);
    }

    return m_edges[index];
}

inline const std::vector<meshkernel::Edge>& meshkernel::Mesh::Edges() const
{
    return m_edges;
}

inline void meshkernel::Mesh::SetEdges(const std::vector<Edge>& newValues)
{
    m_edges = newValues;
    m_nodesRTreeRequiresUpdate = true;
    m_edgesRTreeRequiresUpdate = true;
    m_facesRTreeRequiresUpdate = true;
}
