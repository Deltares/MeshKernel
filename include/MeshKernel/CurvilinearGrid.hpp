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

#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh.hpp>
#include <MeshKernel/RTree.hpp>

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

        /// @brief A struct describing a node in the curvilinear grid in terms of node indices
        struct NodeIndices
        {
            /// @brief Default constructor sets the indices to invalid
            NodeIndices() : m_m(sizetMissingValue), m_n(sizetMissingValue){};

            /// @brief Constructor sets indices from values
            /// @param[in] m The m index
            /// @param[in] n The n index
            NodeIndices(size_t m, size_t n) : m_m(m), m_n(n){};

            /// @brief Determines if one of the indices  equals to \p missingValue
            [[nodiscard]] bool IsValid(const double missingValue = sizetMissingValue) const
            {
                const bool isInvalid = m_m == missingValue || m_n == missingValue;
                return !isInvalid;
            }

            /// @brief Overloads equality with another NodeIndices
            bool operator==(const NodeIndices& rhs) const
            {
                return m_m == rhs.m_m && m_n == rhs.m_n;
            }
            /// @brief Overloads negation with another NodeIndices
            bool operator!=(const NodeIndices& rhs) const
            {
                return !(*this == rhs);
            }

            /// @brief Inquires if another node is on the same grid line of the current node
            /// @param[in] rhs The node to inquire
            /// @return True if on the same grid line, false otherwise
            bool IsOnTheSameGridLine(const NodeIndices& rhs) const
            {
                if (m_m == rhs.m_m || m_n == rhs.m_n)
                {
                    return true;
                }
                return false;
            }

            size_t m_m; ///< Columns
            size_t m_n; ///< Rows
        };

        /// @brief Default constructor
        /// @returns
        CurvilinearGrid() = default;

        /// @brief Rvalue constructor. Creates a new curvilinear grid from a given set of points
        /// @param[in] grid       The input grid points
        /// @param[in] projection The projection to use
        explicit CurvilinearGrid(std::vector<std::vector<Point>>&& grid, Projection projection);

        /// @brief Lvalue constructor. Creates a new curvilinear grid from a given set of points
        /// @param[in] grid       The input grid points
        /// @param[in] projection The projection to use
        explicit CurvilinearGrid(std::vector<std::vector<Point>>& grid, Projection projection);

        /// @brief Check if current curvilinear grid instance is valid
        /// @return True if valid, false otherwise
        [[nodiscard]] bool IsValid() const;

        /// @brief Converting a curvilinear mesh to a set of nodes, edges and returns the original mapping (gridtonet)
        /// @returns The nodes, the edges, and the original mapping (m_m and m_n indices for each node)
        std::tuple<std::vector<Point>, std::vector<Edge>, std::vector<std::pair<size_t, size_t>>> ConvertCurvilinearToNodesAndEdges();

        /// @brief Set internal flat copies of nodes and edges, so the pointer to the first entry is communicated with the front-end
        void SetFlatCopies();

        /// @brief Get the m_m and m_n indices of the node closest to the point
        /// @param[in] point       The input grid points
        NodeIndices GetNodeIndices(Point point);

        /// @brief Computes the grid nodes types and the faces masks
        void ComputeGridNodeTypes();

        /// @brief If the face is valid. A face is valid if all its nodes are valid.
        /// @param[in] m The m coordinate
        /// @param[in] n The n coordinate
        /// @return True if the face is valid, false otherwise
        bool IsValidFace(size_t m, size_t n) const;

        /// @brief From two points expressed as NodeIndices, gets the two corner points defining a block in m_m and m_n coordinates
        /// @param[in] firstNode The node indices of the first node
        /// @param[in] secondNode The node indices of the second node
        /// @return The upper left and lower right of the box defined by the two points
        [[nodiscard]] std::tuple<NodeIndices, NodeIndices> ComputeBlockFromCornerPoints(const NodeIndices& firstNode, const NodeIndices& secondNode) const;

        /// @brief From two points expressed in cartesian coordinates, get the two corner nodes defining a block in m_m and m_n coordinates
        /// @param[in] firstCornerPoint The first corner point
        /// @param[in] secondCornerPoint The second corner point
        /// @return The upper left and lower right nodes of the box
        [[nodiscard]] std::tuple<NodeIndices, NodeIndices> ComputeBlockFromCornerPoints(Point const& firstCornerPoint, Point const& secondCornerPoint);

        /// @brief Function for computing the smoothing factors at the current location given a line and a zone of influence (SMEERFUNCTIE)
        /// The smoothing factor is maximum at the line and 0 at the boundary of the smoothing zone.
        /// @param[in] currentPointIndices The indices of the current point
        /// @param[in] pointOnSmoothingLineIndices The indices of a point on the smoothing line
        /// @param[in] lowerLeftIndices The lower left indices of the smoothing area
        /// @param[in] upperRightIndices The upper right indices of the smoothing area
        /// @return A tuple containing the horizontal, the vertical and mixed smoothing factors
        [[nodiscard]] static std::tuple<double, double, double> ComputeDirectionalSmoothingFactors(NodeIndices const& currentPointIndices,
                                                                                                   NodeIndices const& pointOnSmoothingLineIndices,
                                                                                                   NodeIndices const& lowerLeftIndices,
                                                                                                   NodeIndices const& upperRightIndices);

        size_t m_numM = 0;                                    ///< The number of m_m coordinates (vertical lines)
        size_t m_numN = 0;                                    ///< The number of m_n coordinates (horizontal lines)
        std::vector<std::vector<Point>> m_gridNodes;          ///< Member variable storing the grid
        std::vector<std::vector<bool>> m_gridFacesMask;       ///< The mask of the grid faces (true/false)
        std::vector<std::vector<NodeType>> m_gridNodesMask;   ///< The grid node types
        std::vector<std::pair<size_t, size_t>> m_gridIndices; ///< The original mapping of the flatten nodes in the curvilinear grid

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
        std::vector<Point> ComputeSplineDerivatesAlongGridLine(const std::vector<Point>& gridLine) const;
    };
} // namespace meshkernel
