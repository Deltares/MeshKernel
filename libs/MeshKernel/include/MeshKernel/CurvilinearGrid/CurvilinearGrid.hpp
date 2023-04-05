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

#include <MeshKernel/CurvilinearGrid/CurvilinearGridLine.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridNodeIndices.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh.hpp>

namespace meshkernel
{
    /// @brief A class representing a curvilinear grid
    class CurvilinearGrid : public Mesh
    {

    public:
        /// @brief An enum for curvilinear node types
        enum class NodeType
        {
            BottomLeft,    //(11)
            UpperLeft,     //(14)
            BottomRight,   //(12)
            UpperRight,    //(13)
            Left,          //(4)
            Right,         //(2)
            Bottom,        //(1)
            Up,            //(3)
            InternalValid, //(10)
            Invalid        //(0)
        };

        /// @brief An enum for boundary grid line types
        enum class BoundaryGridLineType
        {
            Left,
            Right,
            Bottom,
            Up
        };

        /// @brief Default constructor
        CurvilinearGrid() = default;

        /// @brief Constructor taking only a projection
        CurvilinearGrid(Projection projection);

        /// @brief Deletes a curvilinear grid inside a polygon
        /// @param[in] polygons The polygons
        /// @param[in] polygonIndex The index of the polygon to use for deletion
        void Delete(std::shared_ptr<Polygons> polygons, size_t polygonIndex);

        /// @brief Lvalue constructor. Creates a new curvilinear grid from a given set of points
        /// @param[in] grid       The input grid points
        /// @param[in] projection The projection to use
        CurvilinearGrid(std::vector<std::vector<Point>> const& grid, Projection projection);

        /// @brief Check if current curvilinear grid instance is valid
        /// @return True if valid, false otherwise
        [[nodiscard]] bool IsValid() const;

        /// @brief Converting a curvilinear mesh to a set of nodes, edges and returns the original mapping (gridtonet)
        /// @returns The nodes, the edges, and the original mapping (m and n indices for each node)
        [[nodiscard]] std::tuple<std::vector<Point>, std::vector<Edge>, std::vector<CurvilinearGridNodeIndices>> ConvertCurvilinearToNodesAndEdges() const;

        /// @brief Set internal flat copies of nodes and edges, so the pointer to the first entry is communicated with the front-end
        void SetFlatCopies();

        /// @brief Get the m and n indices of the node closest to the point
        /// @param[in] point       The input grid points
        [[nodiscard]] CurvilinearGridNodeIndices GetNodeIndices(Point point);

        /// @brief From a point gets the node indices of the closest edges
        /// @param[in] point The input point
        /// @return The curvilinear grid indices of the closest edge
        [[nodiscard]] std::tuple<CurvilinearGridNodeIndices, CurvilinearGridNodeIndices> GetEdgeNodeIndices(Point const& point);

        /// @brief Computes the grid nodes types and the faces masks
        void ComputeGridNodeTypes();

        /// @brief Determines if a face is valid.
        /// A face is valid if all its nodes are valid.
        /// @param[in] m The m coordinate
        /// @param[in] n The n coordinate
        /// @return True if the face is valid, false otherwise
        [[nodiscard]] bool IsValidFace(size_t m, size_t n) const;

        /// @brief Inserts a new face. The new face will be inserted on top of the closest edge.
        /// @param[in] point  The point used for finding the closest edge.
        void InsertFace(Point const& point);

        /// @brief From two points expressed as CurvilinearGridNodeIndices, gets the two corner points defining a block in m and n coordinates
        /// @param[in] firstNode The node indices of the first node
        /// @param[in] secondNode The node indices of the second node
        /// @return The upper left and lower right of the box defined by the two points
        [[nodiscard]] std::tuple<CurvilinearGridNodeIndices, CurvilinearGridNodeIndices> ComputeBlockFromCornerPoints(const CurvilinearGridNodeIndices& firstNode, const CurvilinearGridNodeIndices& secondNode) const;

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

        /// @brief Clones the curvilinear grid instance
        /// @return A pointer to a deep copy of current curvilinear grid instance
        [[nodiscard]] CurvilinearGrid CloneCurvilinearGrid() const;

        /// @brief Allocates a new grid line at the boundary of the curvilinear grid if needed.
        /// @param firstNode The first node of the boundary grid line.
        /// @param secondNode The second node of the boundary grid line.
        /// @return If a new grid line has been allocated
        bool AddGridLineAtBoundary(CurvilinearGridNodeIndices const& firstNode, CurvilinearGridNodeIndices const& secondNode);

        /// @brief Get the boundary grid line type: left, right, bottom or up
        /// @param[in] firstNode The first node of the grid line
        /// @param[in] secondNode The second node of the grid line
        /// @return The boundary grid line type
        [[nodiscard]] BoundaryGridLineType GetBoundaryGridLineType(CurvilinearGridNodeIndices const& firstNode, CurvilinearGridNodeIndices const& secondNode) const;

        /// @brief Delete a node at a specific location by setting it to an invalid point.
        /// @param[in] point The input point coordinate. The closest grid node will be deleted.
        void DeleteNode(Point const& point);

        /// @brief Moves a node from one position to another
        /// @param[in] fromPoint The input position, the closest node will be used
        /// @param[in] toPoint The coordinates of the new position
        void MoveNode(Point const& fromPoint, Point const& toPoint);

        size_t m_numM = 0;                                     ///< The number of m coordinates (vertical lines)
        size_t m_numN = 0;                                     ///< The number of n coordinates (horizontal lines)
        std::vector<std::vector<Point>> m_gridNodes;           ///< Member variable storing the grid
        std::vector<std::vector<bool>> m_gridFacesMask;        ///< The mask of the grid faces (true/false)
        std::vector<std::vector<NodeType>> m_gridNodesTypes;   ///< The grid node types
        std::vector<CurvilinearGridNodeIndices> m_gridIndices; ///< The original mapping of the flatten nodes in the curvilinear grid

    private:
        /// @brief Remove invalid nodes.
        /// This function is recursive
        /// @param[in] invalidNodesToRemove Whether there are still invalid nodes to remove
        void RemoveInvalidNodes(bool invalidNodesToRemove);

        /// @brief Computes the valid grid faces
        void ComputeGridFacesMask();

        /// @brief Compute spline derivatives along a gridline, also accounting for missing values
        /// @param[in] gridLine The input gridline
        /// @returns The spline derivatives
        [[nodiscard]] std::vector<Point> ComputeSplineDerivatesAlongGridLine(const std::vector<Point>& gridLine) const;

        /// @brief Adds an edge at the boundary forming a new face. Increase the grid if required (MODGR1)
        /// The position of the new edge depends on the type of \p firstNode or \p secondNode.
        /// For example, if one of the node types is 'Left' the new edge will be inserted on the left.
        /// The new node will be calculated by a first order approximation: x2 = x1 + (x1 - x0) = 2*x1 - x0
        /// @param[in] firstNode The indices of the first new node in the modified grid.
        /// @param[in] secondNode The indices of the second new node in the modified grid.
        void AddEdge(CurvilinearGridNodeIndices const& firstNode,
                     CurvilinearGridNodeIndices const& secondNode);
    };
} // namespace meshkernel
