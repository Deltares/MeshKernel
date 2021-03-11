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
        enum class NodeTypes
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

        /// @brief Default constructor
        /// @returns
        CurvilinearGrid() = default;

        /// @brief Creates a new curvilinear grid from a given set of points
        /// @param[in] grid       The input grid points
        /// @param[in] projection The projection to use
        CurvilinearGrid(std::vector<std::vector<Point>>&& grid, Projection projection);

        /// @brief Converting a curvilinear mesh to a set of nodes, edges and returns the original mapping (gridtonet)
        /// @returns The nodes, the edges, and the original mapping (m and n indices for each node)
        std::tuple<std::vector<Point>, std::vector<Edge>, std::vector<std::pair<size_t, size_t>>> ConvertCurvilinearToNodesAndEdges();

        /// @brief Set internal flat copies of nodes and edges, so the pointer to the first entry is communicated with the front-end
        void SetFlatCopies();

        /// @brief Get the m and n indices of the node closest to the point
        /// @param[in] point       The input grid points
        std::tuple<int, int> GetNodeIndices(Point point);

        /// @brief Computes the grid nodes types and the faces masks
        void ComputeGridNodeTypes();

        /// @brief If the face is valid. A face is valid if all its nodes are valid.
        /// @param[in] m the m coordinate
        /// @param[in] n the n coordinate
        /// @return true if theface is valid, false otherwise
        bool IsValidFace(size_t m, size_t n) const;

        size_t m_numM = 0;                                    ///< The number of m coordinates (vertical lines)
        size_t m_numN = 0;                                    ///< The number of n coordinates (horizontal lines)
        std::vector<std::vector<Point>> m_gridNodes;          ///< Member variable storing the grid
        std::vector<std::vector<bool>> m_gridFacesMask;       ///< The mask of the grid faces (true/false)
        std::vector<std::vector<NodeTypes>> m_gridNodesMask;  ///< The grid nodes types
        std::vector<std::pair<size_t, size_t>> m_gridIndices; ///< the original mapping of the flatten nodes in the curvilinear grid

    private:
        /// @brief Remove invalid nodes. Is a recursive function
        /// @param[in] invalidNodesToRemove If there are still invalid nodes to remove
        void RemoveInvalidNodes(bool invalidNodesToRemove);

        /// @brief Computes the valid grid faces
        void ComputeGridFacesMask();

        /// @brief Compute spline derivatives along a gridline, also accounting for missing values
        /// @param[in] gridLine The input gridline
        /// @returns The spline derivatives
        std::vector<Point> ComputeSplineDerivatesAlongGridLine(const std::vector<Point>& gridLine) const;
    };
} // namespace meshkernel
