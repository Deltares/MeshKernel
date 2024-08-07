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

#include <set>
#include <vector>

#include <MeshKernel/BoundingBox.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLine.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridNodeIndices.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridUtilities.hpp>
#include <MeshKernel/CurvilinearGrid/UndoActions/AddGridLineUndoAction.hpp>
#include <MeshKernel/CurvilinearGrid/UndoActions/CurvilinearGridBlockUndoAction.hpp>
#include <MeshKernel/CurvilinearGrid/UndoActions/CurvilinearGridRefinementUndoAction.hpp>
#include <MeshKernel/CurvilinearGrid/UndoActions/ResetCurvilinearNodeAction.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/Mesh.hpp>
#include <MeshKernel/UndoActions/UndoAction.hpp>
#include <MeshKernel/Utilities/LinearAlgebra.hpp>

namespace meshkernel
{
    /// @brief A class representing a curvilinear grid
    class CurvilinearGrid final
    {

    public:
        /// @brief Typedef for edge indices
        using CurvilinearEdgeNodeIndices = std::array<CurvilinearGridNodeIndices, 2>;

        /// @brief Typedef for face indices
        using CurvilinearFaceNodeIndices = std::array<CurvilinearGridNodeIndices, 4>;

        /// @brief Typedef for defining a curvilinear edge
        using CurvilinearEdge = std::pair<CurvilinearGridNodeIndices, CurvilinearGridNodeIndices>;

        /// @brief An enum for boundary grid line types
        enum class BoundaryGridLineType
        {
            Bottom,
            Right,
            Top,
            Left
        };

        /// @brief Default constructor
        CurvilinearGrid() = default;

        /// @brief Copy constructor
        CurvilinearGrid(const CurvilinearGrid& grid);

        /// @brief Move constructor
        CurvilinearGrid(CurvilinearGrid&& grid) noexcept;

        /// @brief Constructor taking only a projection
        explicit CurvilinearGrid(Projection projection);

        /// @brief Lvalue constructor. Creates a new curvilinear grid from a given set of points
        /// @details The matrix row index corresponds to the CurvilinearGrid n index, the matrix column index corresponds to the CurvilinearGrid m index
        /// @param[in] grid       The input grid points
        /// @param[in] projection The projection to use
        CurvilinearGrid(lin_alg::Matrix<Point> const& grid, Projection projection);

        /// @brief Creates a new curvilinear grid from a given set of points
        /// @details The matrix row index corresponds to the CurvilinearGrid n index, the matrix column index corresponds to the CurvilinearGrid m index
        /// @param[in] grid       The input grid points
        /// @param[in] projection The projection to use
        CurvilinearGrid(lin_alg::Matrix<Point>&& grid, Projection projection);

        /// @brief Move assignment operator for CurvilinearGrid
        CurvilinearGrid& operator=(CurvilinearGrid&& copy) noexcept;

        /// @brief Ccopy assignment operator for CurvilinearGrid
        CurvilinearGrid& operator=(const CurvilinearGrid& copy);

        /// @brief Set the grid nodes of a curvilinear grid instance
        /// @details The matrix row index corresponds to the CurvilinearGrid n index, the matrix column index corresponds to the CurvilinearGrid m index
        /// @param[in] gridNodes The input grid points
        void SetGridNodes(const lin_alg::Matrix<Point>& gridNodes);

        /// @brief Set the grid nodes of a curvilinear grid instance
        /// @details The matrix row index corresponds to the CurvilinearGrid n index, the matrix column index corresponds to the CurvilinearGrid m index
        /// @param[in] gridNodes The input grid points
        void SetGridNodes(lin_alg::Matrix<Point>&& gridNodes);

        /// @brief Deletes a curvilinear grid inside a polygon
        /// @param[in] polygons The polygons
        /// @param[in] polygonIndex The index of the polygon to use for deletion
        void Delete(std::shared_ptr<Polygons> polygons, UInt polygonIndex);

        /// @brief Check if current curvilinear grid instance is valid
        /// @return True if valid, false otherwise
        [[nodiscard]] inline bool IsValid() const
        {
            return NumM() > 1 && NumN() > 1;
        }

        /// @brief Gets a reference to the grid node at the (m,n) location
        /// @param[in] n The n-dimension index
        /// @param[in] m The m-dimension index
        [[nodiscard]] inline Point& GetNode(const UInt n, const UInt m);

        /// @brief Gets a constant reference to the grid node at the (m,n) location
        /// @param[in] n The n-dimension index
        /// @param[in] m The m-dimension index
        [[nodiscard]] inline Point const& GetNode(const UInt n, const UInt m) const;

        /// @brief Gets a reference to the grid node at the location specified by the index.
        /// @note Exception will be raised for a non-valid index
        /// This is just a helper function, it calls GetNode with (index.m_m, index.m_n)
        [[nodiscard]] inline Point& GetNode(const CurvilinearGridNodeIndices& index);

        /// @brief Get a constant reference to the grid node at the location specified by the index.
        /// @note Exception will be raised for a non-valid index
        /// This is just a helper function, it calls GetNode with (index.m_m, index.m_n)
        [[nodiscard]] inline Point const& GetNode(const CurvilinearGridNodeIndices& index) const;

        /// @brief From a point gets the node indices of the closest edges
        /// @param[in] point The input point
        /// @return The curvilinear grid indices of the closest edge
        [[nodiscard]] std::tuple<CurvilinearGridNodeIndices, CurvilinearGridNodeIndices> GetEdgeNodeIndices(Point const& point);

        /// @brief Computes the grid nodes types and the faces masks
        void ComputeGridNodeTypes();

        /// @brief Determines the grid node type
        /// @param[in] n The n-dimension index
        /// @param[in] m The m-dimension index
        /// @return the node type
        NodeType GetNodeType(UInt n, UInt m) const
        {
            if (n >= m_gridNodesTypes.rows()) [[unlikely]]
            {
                throw ConstraintError("Invalid row index {} > {}", n, m_gridNodesTypes.rows());
            }

            if (m >= m_gridNodesTypes.cols()) [[unlikely]]
            {
                throw ConstraintError("Invalid column index {} > {}", m, m_gridNodesTypes.cols());
            }

            return m_gridNodesTypes(n + m_startOffset.m_n, m + m_startOffset.m_m);
        }

        /// @brief Determines the grid node type
        /// @param[in] n The n-dimension index
        /// @param[in] m The m-dimension index
        /// @return reference to the node type
        NodeType& GetNodeType(UInt n, UInt m)
        {
            if (n >= m_gridNodesTypes.rows()) [[unlikely]]
            {
                throw ConstraintError("Invalid row index {} > {}", n, m_gridNodesTypes.rows());
            }

            if (m >= m_gridNodesTypes.cols()) [[unlikely]]
            {
                throw ConstraintError("Invalid column index {} > {}", m, m_gridNodesTypes.cols());
            }

            return m_gridNodesTypes(n + m_startOffset.m_n, m + m_startOffset.m_m);
        }

        /// @brief Determines if all nodes of a face are valid.
        /// A face is valid if all its nodes are valid.
        /// @param[in] n The n-dimension index
        /// @param[in] m The m-dimension index
        /// @return True if the face is valid, false otherwise
        [[nodiscard]] bool AreFaceNodesValid(UInt n, UInt m) const;

        /// @brief Determines if the face mask is true (valid face) or false (invalid face)
        /// @param[in] n The n-dimension index
        /// @param[in] m The m-dimension index
        /// @return the face mask value (true/false)
        [[nodiscard]] bool IsFaceMaskValid(UInt n, UInt m) const
        {
            if (n >= m_gridFacesMask.rows()) [[unlikely]]
            {
                throw ConstraintError("Invalid row index {} > {}", n, m_gridFacesMask.rows());
            }

            if (m >= m_gridFacesMask.cols()) [[unlikely]]
            {
                throw ConstraintError("Invalid column index {} > {}", m, m_gridFacesMask.cols());
            }

            return m_gridFacesMask(n + m_startOffset.m_n, m + m_startOffset.m_m);
        }

        /// @brief Determines if the face mask is true (valid face) or false (invalid face)
        /// @param[in] n The n-dimension index
        /// @param[in] m The m-dimension index
        /// @return reference to the face mask value
        bool& IsFaceMaskValid(UInt n, UInt m)
        {
            if (n >= m_gridFacesMask.rows()) [[unlikely]]
            {
                throw ConstraintError("Invalid row index {} > {}", n, m_gridFacesMask.rows());
            }

            if (m >= m_gridFacesMask.cols()) [[unlikely]]
            {
                throw ConstraintError("Invalid column index {} > {}", m, m_gridFacesMask.cols());
            }

            return m_gridFacesMask(n + m_startOffset.m_n, m + m_startOffset.m_m);
        }

        /// @brief Inserts a new face. The new face will be inserted on top of the closest edge.
        /// @param[in] point  The point used for finding the closest edge.
        [[nodiscard]] UndoActionPtr InsertFace(Point const& point);

        /// @brief From two points expressed as CurvilinearGridNodeIndices, gets the two corner points defining a block in m and n coordinates
        /// @param[in] firstNode The node indices of the first node
        /// @param[in] secondNode The node indices of the second node
        /// @return The upper left and lower right of the box defined by the two points
        std::tuple<CurvilinearGridNodeIndices, CurvilinearGridNodeIndices> ComputeBlockFromCornerPoints(const CurvilinearGridNodeIndices& firstNode, const CurvilinearGridNodeIndices& secondNode) const;

        /// @brief From two points expressed in cartesian coordinates, get the two corner nodes defining a block in m and n coordinates
        /// @param[in] firstCornerPoint The first corner point
        /// @param[in] secondCornerPoint The second corner point
        /// @return The upper left and lower right nodes of the box
        [[nodiscard]] std::tuple<CurvilinearGridNodeIndices, CurvilinearGridNodeIndices> ComputeBlockFromCornerPoints(Point const& firstCornerPoint, Point const& secondCornerPoint);

        /// @brief Function for computing the smoothing factors at the current location given a line and a zone of influence (SMEERFUNCTIE)
        /// The smoothing factor is maximum at the line and 0 at the boundary of the smoothing zone.
        /// @param[in] currentPointIndices The indices of the current point
        /// @param[in] pointOnSmoothingLineIndices The indices of a point on the smoothing line
        /// @param[in] lowerLeftIndices The lower left indices of the smoothing area
        /// @param[in] upperRightIndices The upper right indices of the smoothing area
        /// @return A tuple containing the horizontal, the vertical and mixed smoothing factors
        [[nodiscard]] static std::tuple<double, double, double> ComputeDirectionalSmoothingFactors(CurvilinearGridNodeIndices const& currentPointIndices,
                                                                                                   CurvilinearGridNodeIndices const& pointOnSmoothingLineIndices,
                                                                                                   CurvilinearGridNodeIndices const& lowerLeftIndices,
                                                                                                   CurvilinearGridNodeIndices const& upperRightIndices);

        /// @brief Computes an average distance of the current node from the other nodes (DXB)
        /// @param[in] nodeIndex The current node index
        /// @param[in] direction The direction, either m or n
        /// @return The computed distance
        [[nodiscard]] double ComputeAverageNodalDistance(CurvilinearGridNodeIndices const& nodeIndex, CurvilinearGridLine::GridLineDirection direction);

        /// @brief Transform the displacement around a node to local or global (TOLOCL)
        /// @param[in] displacement The displacement to transform.
        /// @param[in] node The node position
        /// @param[in] toLocal A boolean to indicate whatever to transform the displacement to local grid (True) or to global grid (false)
        /// @return The new displacement
        [[nodiscard]] Point TransformDisplacement(Point const& displacement, CurvilinearGridNodeIndices const& node, bool isLocal) const;

        /// @brief Allocates a new grid line at the boundary of the curvilinear grid if needed.
        /// @param firstNode The first node of the boundary grid line.
        /// @param secondNode The second node of the boundary grid line.
        /// @return If a new grid line has been allocated
        std::tuple<bool, UndoActionPtr> AddGridLineAtBoundary(CurvilinearGridNodeIndices const& firstNode, CurvilinearGridNodeIndices const& secondNode);

        /// @brief Get the boundary grid line type: left, right, bottom or up
        /// @param[in] firstNode The first node of the grid line
        /// @param[in] secondNode The second node of the grid line
        /// @return The boundary grid line type
        [[nodiscard]] BoundaryGridLineType GetBoundaryGridLineType(CurvilinearGridNodeIndices const& firstNode, CurvilinearGridNodeIndices const& secondNode) const;

        /// @brief Delete a node at a specific location by setting it to an invalid point.
        /// @param[in] point The input point coordinate. The closest grid node will be deleted.
        [[nodiscard]] UndoActionPtr DeleteNode(Point const& point);

        /// @brief Moves a node from one position to another
        /// @param[in] nodeIndex The input index
        /// @param[in] toPoint The coordinates of the new position
        [[nodiscard]] UndoActionPtr MoveNode(const CurvilinearGridNodeIndices& fromPoint, Point const& toPoint);

        /// @brief Moves a node from one position to another
        /// @param[in] fromPoint The input position, the closest node will be used
        /// @param[in] toPoint The coordinates of the new position
        [[nodiscard]] UndoActionPtr MoveNode(Point const& fromPoint, Point const& toPoint);

        /// @brief Converting a curvilinear mesh to a set of nodes
        /// @details Nodes and grid indices from the matrix are serialized in row-major order (n runs fastest).
        /// @note Also invalid nodes are serialized to points
        /// @returns The nodes
        [[nodiscard]] std::vector<Point> ComputeNodes() const;

        /// @brief Computes the number of nodes
        [[nodiscard]] UInt GetNumNodes() const { return NumN() * NumM(); }

        /// @brief Computes the number of edges
        [[nodiscard]] UInt GetNumEdges() const { return NumM() * (NumN() - 1) +
                                                        (NumM() - 1) * NumN(); }

        /// @brief Computes the number of faces
        [[nodiscard]] UInt GetNumFaces() const { return (NumN() - 1) * (NumM() - 1); }

        /// @brief Computes the node from an index
        [[nodiscard]] const Point& Node(const UInt index) const
        {
            const auto nodeIndex = GetCurvilinearGridNodeIndices(index);
            return GetNode(nodeIndex);
        }

        /// @brief Get the node CurvilinearGridNodeIndices from a position
        /// @return The CurvilinearGridNodeIndices for the position
        CurvilinearGridNodeIndices GetCurvilinearGridNodeIndices(const UInt index) const
        {
            if (index >= NumM() * NumN())
            {
                throw AlgorithmError("Invalid index");
            }

            const UInt n = index / NumM();
            const UInt m = index % NumM();
            return {n, m};
        }

        /// @brief Get the node index from a CurvilinearGridNodeIndices
        /// @return The CurvilinearGridNodeIndices
        UInt GetNodeIndex(const CurvilinearGridNodeIndices& nodeIndex) const
        {
            if (!nodeIndex.IsValid())
            {
                throw AlgorithmError("CurvilinearGridNodeIndices ");
            }

            return nodeIndex.m_m + nodeIndex.m_n * NumM();
        }

        /// @brief Finds the index of the closest location
        /// @param[in] point the input point
        /// @param[in] location the location type
        /// @param[in] boundingBox The bounding box
        /// @returns The location index
        UInt FindLocationIndex(Point point,
                               Location location,
                               const BoundingBox& boundingBox = {});

        /// @brief Get the projection
        /// @return The curvilinear grid projection
        [[nodiscard]] const Projection& projection() const { return m_projection; }

        /// @brief Computes the edge centers in the correct order
        /// @returns the edge centers
        [[nodiscard]] std::vector<Point> ComputeEdgesCenters() const;

        /// @brief Compute the node indices in row-major order (n runs fastest).
        /// @returns The  mapping (m and n indices for each node)
        [[nodiscard]] std::vector<Point> ComputeFaceCenters() const;

        /// @brief Converting a curvilinear mesh to a set of edges
        /// @details Edges are serialized as follows: first all m-oriented edges ((m,n)-(m+1,n)) in row-major order, then all
        /// n-oriented edges ((m,n)-(m,n+1)), in row-major order.
        /// @returns The edges
        [[nodiscard]] std::vector<Edge> ComputeEdges() const;

        /// @brief Get the mesh bounding box.
        [[nodiscard]] BoundingBox GetBoundingBox() const;

        /// @brief The number of nodes M in the m dimension
        /// @return A number >= 2 for a valid curvilinear grid
        UInt NumM() const { return static_cast<UInt>(m_gridNodes.cols()) - m_startOffset.m_m - m_endOffset.m_m; }

        /// @brief The number of nodes N in the n dimension
        /// @return A number >= 2 for a valid curvilinear grid
        UInt NumN() const { return static_cast<UInt>(m_gridNodes.rows()) - m_startOffset.m_n - m_endOffset.m_n; }

        /// @brief Get the row and column start index offset
        CurvilinearGridNodeIndices StartOffset() const;

        /// @brief Get the row and column end index offset
        CurvilinearGridNodeIndices EndOffset() const;

        /// @brief Is the node matrix empty
        /// @return true iff the node matrix is empty
        bool IsEmpty() const { return lin_alg::MatrixIsEmpty(m_gridNodes); }

        /// @brief Get a copy of the nodes matrix
        /// @note:the m dimension is the first (or row) index, the n dimension is the second (or column) index
        /// @return a copy of the matrix
        lin_alg::Matrix<Point> GetNodes() const { return m_gridNodes; }

        /// @brief Get the array of nodes at an m-dimension index
        /// @param [in] m the m-dimension index
        /// @return a vector of N nodes
        std::vector<Point> GetNodeVectorAtM(UInt m) const { return lin_alg::MatrixRowToSTLVector(m_gridNodes, m + m_startOffset.m_m); }

        /// @brief Get the array of nodes at an n-dimension index
        /// @param [in] n the n-dimension index
        /// @return a vector of M nodes
        std::vector<Point> GetNodeVectorAtN(UInt n) const { return lin_alg::MatrixColToSTLVector(m_gridNodes, n + m_startOffset.m_n); }

        /// @brief Compute a sequence of one or more outer boundary polygons separated by geometry separators
        /// @param[in]  lowerLeft   The node index of the lower left corner
        /// @param[in]  upperRight  The node index of the upper right corner
        /// @returns The vector containing the boundary polygon points
        [[nodiscard]] std::vector<Point> ComputeBoundaryPolygons(const CurvilinearGridNodeIndices& lowerLeft, const CurvilinearGridNodeIndices& upperRight) const;

        /// @brief The number of nodes M in the m dimension
        /// @return A number >= 2 for a valid curvilinear grid
        UInt FullNumM() const { return static_cast<UInt>(m_gridNodes.cols()); }

        /// @brief The number of nodes N in the n dimension
        /// @return A number >= 2 for a valid curvilinear grid
        UInt FullNumN() const { return static_cast<UInt>(m_gridNodes.rows()); }

        /// @brief Restore grid to state before grid line was added
        void RestoreAction(const AddGridLineUndoAction& undoAction);

        /// @brief Restore grid to state after grid line was added
        void CommitAction(const AddGridLineUndoAction& undoAction);

        /// @brief Restore grid to state before grid block was modified
        ///
        /// The modification could be from e.g. orthogonalisation, delete interior, ...
        void RestoreAction(CurvilinearGridBlockUndoAction& undoAction);

        /// @brief Restore grid to state after grid block was modified
        void CommitAction(CurvilinearGridBlockUndoAction& undoAction);

        /// @brief Restore grid to state before refinement operation
        void RestoreAction(CurvilinearGridRefinementUndoAction& undoAction);

        /// @brief Restore grid to state after refinement operation
        void CommitAction(CurvilinearGridRefinementUndoAction& undoAction);

        /// @brief Restore grid to state before node was modified
        void RestoreAction(const ResetCurvilinearNodeAction& undoAction);

        /// @brief Restore grid to state after node was modified
        void CommitAction(const ResetCurvilinearNodeAction& undoAction);

    private:
        /// @brief Remove invalid nodes.
        /// This function is recursive
        /// @param[in] invalidNodesToRemove Whether there are still invalid nodes to remove
        void RemoveInvalidNodes(bool invalidNodesToRemove);

        /// @brief Adds an edge at the boundary forming a new face. Increase the grid if required (MODGR1)
        /// The position of the new edge depends on the type of \p firstNode or \p secondNode.
        /// For example, if one of the node types is 'Left' the new edge will be inserted on the left.
        /// The new node will be calculated by a first order approximation: x2 = x1 + (x1 - x0) = 2*x1 - x0
        /// @param[in] firstNode The indices of the first new node in the modified grid.
        /// @param[in] secondNode The indices of the second new node in the modified grid.
        [[nodiscard]] UndoActionPtr AddEdge(CurvilinearGridNodeIndices const& firstNode,
                                            CurvilinearGridNodeIndices const& secondNode);

        /// @brief Build the rtree for the corresponding location, using only the locations inside the bounding box
        /// @param[in] location The mesh location for which the RTree is build
        /// @param[in] boundingBox The bounding box
        void BuildTree(Location location, const BoundingBox& boundingBox = {});

        /// @brief Computes the valid grid faces
        void ComputeGridFacesMask();

        /// @brief Compute the node indices in row-major order (n runs fastest).
        /// @returns The  mapping (m and n indices for each node)
        [[nodiscard]] std::vector<CurvilinearGridNodeIndices> ComputeNodeIndices() const;

        /// @brief Compute the edge indices.
        /// @returns The  mapping (m and n indices for each node of the edge)
        [[nodiscard]] std::vector<CurvilinearEdgeNodeIndices> ComputeEdgeIndices() const;

        /// @brief Compute the boundary edges. While iterating over edges of valid faces, boundary edges are seen once, all internal edges are seen by two faces.
        /// @param[in]  lowerLeft   The node index of the lower left corner
        /// @param[in]  upperRight  The node index of the upper right corner
        /// @returns The set of boundary edges
        [[nodiscard]] std::set<CurvilinearEdge> ComputeBoundaryEdges(const CurvilinearGridNodeIndices& lowerLeft, const CurvilinearGridNodeIndices& upperRight) const;

        /// @brief Set the m_nodesRTreeRequiresUpdate flag
        /// @param[in] value The value of the flag
        void SetNodesRTreeRequiresUpdate(bool value) { m_nodesRTreeRequiresUpdate = value; }

        /// @brief Set the m_edgesRTreeRequiresUpdate flag
        /// @param[in] value The value of the flag
        void SetEdgesRTreeRequiresUpdate(bool value) { m_edgesRTreeRequiresUpdate = value; }

        /// @brief Set the m_facesRTreeRequiresUpdate flag
        /// @param[in] value The value of the flag
        void SetFacesRTreeRequiresUpdate(bool value) { m_facesRTreeRequiresUpdate = value; }

        Projection m_projection;                               ///< The curvilinear grid projection
        lin_alg::Matrix<Point> m_gridNodes;                    ///< Member variable storing the grid
        lin_alg::Matrix<bool> m_gridFacesMask;                 ///< The mask of the grid faces (true/false)
        lin_alg::Matrix<NodeType> m_gridNodesTypes;            ///< The grid node types
        std::vector<CurvilinearGridNodeIndices> m_gridIndices; ///< The original mapping of the flatten nodes in the curvilinear grid

        // RTrees
        bool m_nodesRTreeRequiresUpdate = true;                            ///< m_nodesRTree requires an update
        bool m_edgesRTreeRequiresUpdate = true;                            ///< m_edgesRTree requires an update
        bool m_facesRTreeRequiresUpdate = true;                            ///< m_facesRTree requires an update
        std::unordered_map<Location, std::unique_ptr<RTreeBase>> m_RTrees; ///< The RTrees to use
        BoundingBox m_boundingBoxCache;                                    ///< Caches the last bounding box used for selecting the locations

        std::vector<Edge> m_edges; ///< Member variable storing the edges

        /// @brief
        CurvilinearGridNodeIndices m_startOffset{0, 0}; ///< Row and column start index offset
        CurvilinearGridNodeIndices m_endOffset{0, 0};   ///< Row and column end index offset
    };
} // namespace meshkernel

inline meshkernel::CurvilinearGridNodeIndices meshkernel::CurvilinearGrid::StartOffset() const
{
    return m_startOffset;
}

inline meshkernel::CurvilinearGridNodeIndices meshkernel::CurvilinearGrid::EndOffset() const
{
    return m_endOffset;
}

inline meshkernel::Point& meshkernel::CurvilinearGrid::GetNode(const UInt n, const UInt m)
{

    if (n >= m_gridNodes.rows()) [[unlikely]]
    {
        throw ConstraintError("Invalid row index {} >= {}", n, m_gridNodes.rows());
    }

    if (m >= m_gridNodes.cols()) [[unlikely]]
    {
        throw ConstraintError("Invalid column index {} >= {}", m, m_gridNodes.cols());
    }

    m_nodesRTreeRequiresUpdate = true;
    m_edgesRTreeRequiresUpdate = true;
    m_facesRTreeRequiresUpdate = true;

    return m_gridNodes(n + m_startOffset.m_n, m + m_startOffset.m_m);
}

meshkernel::Point const& meshkernel::CurvilinearGrid::GetNode(const UInt n, const UInt m) const
{

    if (n >= m_gridNodes.rows()) [[unlikely]]
    {
        throw ConstraintError("Invalid row index {} >= {}", n, m_gridNodes.rows());
    }

    if (m >= m_gridNodes.cols()) [[unlikely]]
    {
        throw ConstraintError("Invalid column index {} >= {}", m, m_gridNodes.cols());
    }
    return m_gridNodes(n + m_startOffset.m_n, m + m_startOffset.m_m);
}

meshkernel::Point& meshkernel::CurvilinearGrid::GetNode(const CurvilinearGridNodeIndices& index)
{
    if (!index.IsValid()) [[unlikely]]
    {
        throw ConstraintError("Invalid node index");
    }
    m_nodesRTreeRequiresUpdate = true;
    m_edgesRTreeRequiresUpdate = true;
    m_facesRTreeRequiresUpdate = true;

    return GetNode(index.m_n, index.m_m);
}

meshkernel::Point const& meshkernel::CurvilinearGrid::GetNode(const CurvilinearGridNodeIndices& index) const
{
    if (!index.IsValid()) [[unlikely]]
    {
        throw ConstraintError("Invalid node index");
    }

    return GetNode(index.m_n, index.m_m);
}
