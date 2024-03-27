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

#include <vector>

#include <MeshKernel/BoundingBox.hpp>
#include <MeshKernel/CurvilinearGrid/AddGridLineUndoAction.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridBlockUndo.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLine.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridNodeIndices.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridRefinementUndoAction.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridUtilities.hpp>
#include <MeshKernel/CurvilinearGrid/ResetCurvilinearNodeAction.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/Mesh.hpp>
#include <MeshKernel/UndoActions/UndoAction.hpp>
#include <MeshKernel/Utilities/LinearAlgebra.hpp>

namespace meshkernel
{
    /// @brief A class representing a curvilinear grid
    class CurvilinearGrid final : public Mesh
    {

    public:
        // /// @brief An enum for curvilinear node types
        // enum class NodeType
        // {
        //     BottomLeft,    //(11)
        //     UpperLeft,     //(14)
        //     BottomRight,   //(12)
        //     UpperRight,    //(13)
        //     Left,          //(4)
        //     Right,         //(2)
        //     Bottom,        //(1)
        //     Up,            //(3)
        //     InternalValid, //(10)
        //     Invalid        //(0)
        // };

        /// @brief An enum for boundary grid line types
        enum class BoundaryGridLineType
        {
            Left,   ///< Bottom of domain
            Right,  ///<
            Bottom, ///<
            Up      ///< Right side of domain
        };

        /// @brief Default destructor
        ~CurvilinearGrid() override = default;

        /// @brief Default constructor
        CurvilinearGrid() = default;

        /// @brief Copy constructor taking only a curvilinear grid
        CurvilinearGrid(const CurvilinearGrid& grid);

        /// @brief Constructor taking only a projection
        explicit CurvilinearGrid(Projection projection);

        /// @brief Lvalue constructor. Creates a new curvilinear grid from a given set of points
        /// @details The matrix row index corresponds to the CurvilinearGrid m index, the matrix column index corresponds to the CurvilinearGrid n index
        /// @param[in] grid       The input grid points
        /// @param[in] projection The projection to use
        CurvilinearGrid(lin_alg::Matrix<Point> const& grid, Projection projection);

        /// @brief Set the grid nodes of a curvilinear grid instance
        /// @details The matrix row index corresponds to the CurvilinearGrid m index, the matrix column index corresponds to the CurvilinearGrid n index
        /// @param[in] gridNodes The input grid points
        void SetGridNodes(const lin_alg::Matrix<Point>& gridNodes);

        /// @brief Deletes a curvilinear grid inside a polygon
        /// @param[in] polygons The polygons
        /// @param[in] polygonIndex The index of the polygon to use for deletion
        void Delete(std::shared_ptr<Polygons> polygons, UInt polygonIndex);

        /// @brief Check if current curvilinear grid instance is valid
        /// @return True if valid, false otherwise
        [[nodiscard]] bool IsValid() const;

        /// @brief Converting a curvilinear mesh to a set of nodes, edges and returns the original mapping (gridtonet)
        /// @details Nodes and grid indices from the matrix are serialized in row-major order (n runs fastest).
        /// Edges are serialized as follows: first all m-oriented edges ((m,n)-(m+1,n)) in row-major order, then all
        /// n-oriented edges ((m,n)-(m,n+1)), in row-major order.
        /// @note Also invalid nodes are serialized to points, edges and grid indices
        /// @returns The nodes, the edges, and the original mapping (m and n indices for each node)
        [[nodiscard]] std::tuple<std::vector<Point>, std::vector<Edge>, std::vector<CurvilinearGridNodeIndices>> ConvertCurvilinearToNodesAndEdges() const;

        /// @brief Set internal flat copies of nodes and edges, so the pointer to the first entry is communicated with the front-end
        /// @details The Mesh nodes and edges arrays, and the grid node indices array are populated by the result of ConvertCurvilinearToNodesAndEdges.
        void SetFlatCopies();

        /// @brief Get the m and n indices of the node closest to the point
        /// @param[in] point       The input grid points
        [[nodiscard]] CurvilinearGridNodeIndices GetNodeIndices(Point point);

        /// @brief Gets a reference to the grid node at the (m,n) location
        /// @param[in] n The n-dimension index
        /// @param[in] m The m-dimension index
        [[nodiscard]] meshkernel::Point& GetNode(const UInt n, const UInt m) { return m_gridNodes(n + m_startOffset.m_n, m + m_startOffset.m_m); }

        /// @brief Gets a constant reference to the grid node at the (m,n) location
        /// @param[in] n The n-dimension index
        /// @param[in] m The m-dimension index
        [[nodiscard]] meshkernel::Point const& GetNode(const UInt n, const UInt m) const { return m_gridNodes(n + m_startOffset.m_n, m + m_startOffset.m_m); }

        /// @brief Gets a reference to the grid node at the location specified by the index.
        /// @note Exception will be raised for a non-valid index
        /// This is just a helper function, it calls GetNode with (index.m_m, index.m_n)
        [[nodiscard]] meshkernel::Point& GetNode(const CurvilinearGridNodeIndices& index)
        {
            if (!index.IsValid())
            {
                throw ConstraintError("Invalid node index");
            }
            return GetNode(index.m_n, index.m_m);
            // return m_gridNodes(index.m_n, index.m_m);
        }

        /// @brief Get a constant reference to the grid node at the location specified by the index.
        /// @note Exception will be raised for a non-valid index
        /// This is just a helper function, it calls GetNode with (index.m_m, index.m_n)
        [[nodiscard]] meshkernel::Point const& GetNode(const CurvilinearGridNodeIndices& index) const
        {
            if (!index.IsValid())
            {
                throw ConstraintError("Invalid node index");
            }
            return GetNode(index.m_n, index.m_m);
            // return m_gridNodes(index.m_n, index.m_m);
        }

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
        NodeType GetNodeType(UInt n, UInt m) const { return m_gridNodesTypes(n + m_startOffset.m_n, m + m_startOffset.m_m); }

        NodeType& GetNodeType(UInt n, UInt m) { return m_gridNodesTypes(n + m_startOffset.m_n, m + m_startOffset.m_m); }

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
        [[nodiscard]] bool IsFaceMaskValid(UInt n, UInt m) const { return m_gridFacesMask(n + m_startOffset.m_n, m + m_startOffset.m_m); }

        bool& IsFaceMaskValid(UInt n, UInt m) { return m_gridFacesMask(n + m_startOffset.m_n, m + m_startOffset.m_m); }

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

        void DeleteGridLineAtBoundary(CurvilinearGridNodeIndices const& firstNode, CurvilinearGridNodeIndices const& secondNode);

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
        [[nodiscard]] UndoActionPtr MoveNode(const CurvilinearGridNodeIndices& nodeIndex, Point const& toPoint);

        /// @brief Moves a node from one position to another
        /// @param[in] fromPoint The input position, the closest node will be used
        /// @param[in] toPoint The coordinates of the new position
        [[nodiscard]] UndoActionPtr MoveNode(Point const& fromPoint, Point const& toPoint);

        /// @brief Get the mesh bounding box.
        BoundingBox GetBoundingBox() const;

        /// @brief The number of nodes M in the m dimension
        /// @return A number >= 2 for a valid curvilinear grid
        UInt NumM() const { return static_cast<UInt>(m_gridNodes.cols()) - m_startOffset.m_m - m_endOffset.m_m; }

        /// @brief The number of nodes N in the n dimension
        /// @return A number >= 2 for a valid curvilinear grid
        UInt NumN() const { return static_cast<UInt>(m_gridNodes.rows()) - m_startOffset.m_n - m_endOffset.m_n; }

        /// @brief Is the node matrix empty
        /// @return true iff the node matrix is empty
        bool IsEmpty() const { return lin_alg::MatrixIsEmpty(m_gridNodes); }

        /// @brief Get a copy of the nodes matrix
        /// @note:the m dimension is the first (or row) index, the n dimension is the second (or column) index
        /// @return a copy of the matrix
        // TODO need to retuen a slice of the matrix m_gridNodes (m_startOffset.m_m .. m_endOffset.m_m, m_startOffset.m_n .. m_endOffset.m_n)
        lin_alg::Matrix<Point> GetNodes() const { return m_gridNodes; }

        /// @brief Get the array of nodes at an m-dimension index
        /// @param [in] m the m-dimension index
        /// @return a vector of N nodes
        // TODO Need to include the column offset too, m_startOffset.m_n
        std::vector<Point> GetNodeVectorAtM(UInt m) const { return lin_alg::MatrixRowToSTLVector(m_gridNodes, m + m_startOffset.m_m); }

        /// @brief Get the array of nodes at an n-dimension index
        /// @param [in] n the n-dimension index
        /// @return a vector of M nodes
        // TODO Need to include the column offset too, m_startOffset.m_m
        std::vector<Point> GetNodeVectorAtN(UInt n) const { return lin_alg::MatrixColToSTLVector(m_gridNodes, n + m_startOffset.m_n); }

        /// @brief The number of nodes M in the m dimension
        /// @return A number >= 2 for a valid curvilinear grid
        UInt FullNumM() const { return static_cast<UInt>(m_gridNodes.cols()); }

        /// @brief The number of nodes N in the n dimension
        /// @return A number >= 2 for a valid curvilinear grid
        UInt FullNumN() const { return static_cast<UInt>(m_gridNodes.rows()); }

        /// @brief Restore grid to state before grid line was added
        void RestoreAction(AddGridLineUndoAction& undoAction);

        /// @brief Restore grid to state after grid line was added
        void CommitAction(AddGridLineUndoAction& undoAction);

        /// @brief Restore grid to state before grid block was modified
        ///
        /// The modification could be from e.g. orthogonalisation, delete interior, ...
        void RestoreAction(CurvilinearGridBlockUndo& undoAction);

        /// @brief Restore grid to state after grid block was modified
        void CommitAction(CurvilinearGridBlockUndo& undoAction);

        /// @brief Restore grid to state before refinement operation
        void RestoreAction(CurvilinearGridRefinementUndoAction& undoAction);

        /// @brief Restore grid to state after refinement operation
        void CommitAction(CurvilinearGridRefinementUndoAction& undoAction);

        /// @brief Restore grid to state before node was modified
        void RestoreAction(const ResetCurvilinearNodeAction& undoAction);

        /// @brief Restore grid to state after node was modified
        void CommitAction(const ResetCurvilinearNodeAction& undoAction);

        // TODO make private;
        CurvilinearGridNodeIndices m_startOffset{0, 0};
        CurvilinearGridNodeIndices m_endOffset{0, 0};

        void print(std::ostream& out = std::cout);

        void printGraph(std::ostream& out = std::cout);

    private:
        /// @brief Remove invalid nodes.
        /// This function is recursive
        /// @param[in] invalidNodesToRemove Whether there are still invalid nodes to remove
        void RemoveInvalidNodes(bool invalidNodesToRemove);

        /// @brief Computes the valid grid faces
        void ComputeGridFacesMask();

        /// @brief Adds an edge at the boundary forming a new face. Increase the grid if required (MODGR1)
        /// The position of the new edge depends on the type of \p firstNode or \p secondNode.
        /// For example, if one of the node types is 'Left' the new edge will be inserted on the left.
        /// The new node will be calculated by a first order approximation: x2 = x1 + (x1 - x0) = 2*x1 - x0
        /// @param[in] firstNode The indices of the first new node in the modified grid.
        /// @param[in] secondNode The indices of the second new node in the modified grid.
        [[nodiscard]] UndoActionPtr AddEdge(CurvilinearGridNodeIndices const& firstNode,
                                            CurvilinearGridNodeIndices const& secondNode);

        lin_alg::Matrix<Point> m_gridNodes;                    ///< Member variable storing the grid
        lin_alg::Matrix<bool> m_gridFacesMask;                 ///< The mask of the grid faces (true/false)
        lin_alg::Matrix<NodeType> m_gridNodesTypes;            ///< The grid node types
        std::vector<CurvilinearGridNodeIndices> m_gridIndices; ///< The original mapping of the flatten nodes in the curvilinear grid
    };
} // namespace meshkernel
