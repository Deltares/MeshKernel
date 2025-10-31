//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2023.
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

#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"
#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/Point.hpp"
#include "Utilities/LinearAlgebra.hpp"

namespace meshkernel
{
    /// @brief Construct a curvilinear grid from an unstructured mesh
    class Mesh2DToCurvilinear
    {
    public:
        /// @brief Constructor
        /// @param[in] mesh The input unstructured mesh
        explicit Mesh2DToCurvilinear(Mesh2D& mesh);

        /// @brief Computes the curvilinear grid starting from a specific point
        /// @param[in] point The point from which start growing the curvilinear mesh. The point must be inside a quadrangular face
        /// @returns The curvilinear grid
        std::unique_ptr<CurvilinearGrid> Compute(const Point& point);

    private:
        /// @brief Matrix class supporting negative indices, mimicking behavior similar to Fortran.
        /// This class allows matrix access with negative indices, which is not possible in C++ by default.
        class MatrixWithNegativeIndices
        {
        public:
            /// @brief Retrieves the value at the given position (j, i) in the matrix.
            /// @param j Row index, can be negative.
            /// @param i Column index, can be negative.
            /// @return Value at the specified position.
            int getValue(int j, int i) const
            {
                const auto jIndex = j - m_minJ;
                const auto iIndex = i - m_minI;
                const auto value = m_matrix(jIndex, iIndex);
                return value;
            }
            /// @brief Sets the value at the given position (j, i) in the matrix.
            /// @param j Row index, can be negative.
            /// @param i Column index, can be negative.
            /// @param value Value to be set at the specified position.
            void setValue(int j, int i, int value)
            {
                const auto jIndex = j - m_minJ;
                const auto iIndex = i - m_minI;
                m_matrix(jIndex, iIndex) = value;
            }

            /// @brief Checks if the position (j, i) in the matrix is valid (not equal to missing value).
            /// @param j Row index, can be negative.
            /// @param i Column index, can be negative.
            /// @return True if the value at the position is not missing, otherwise false.
            bool IsValid(int j, int i) const
            {
                const auto jIndex = j - m_minJ;
                const auto iIndex = i - m_minI;
                return m_matrix(jIndex, iIndex) != constants::missing::intValue;
            }

            /// @brief Gets the number of rows in the matrix.
            /// @return Number of rows.
            [[nodiscard]] int rows() const
            {
                return static_cast<int>(m_matrix.rows());
            }

            /// @brief Gets the number of columns in the matrix.
            /// @return Number of columns.
            [[nodiscard]] int cols() const
            {
                return static_cast<int>(m_matrix.cols());
            }

            /// @brief Gets the minimum column index (i) allowed in the matrix.
            /// @return Minimum column index.
            [[nodiscard]] int minCol() const
            {
                return m_minI;
            }

            /// @brief Gets the maximum column index (i) allowed in the matrix.
            /// @return Maximum column index.
            [[nodiscard]] int maxCol() const
            {
                return m_maxI;
            }

            /// @brief Gets the minimum row index (j) allowed in the matrix.
            /// @return Minimum row index.
            [[nodiscard]] int minRow() const
            {
                return m_minJ;
            }

            /// @brief Gets the maximum row index (j) allowed in the matrix.
            /// @return Maximum row index.
            [[nodiscard]] int maxRow() const
            {
                return m_maxJ;
            }

            /// @brief Resizes the matrix to accommodate the new bounds.
            /// @param minJ New minimum row index.
            /// @param minI New minimum column index.
            /// @param maxJ New maximum row index.
            /// @param maxI New maximum column index.
            void resize(int minJ, int minI, int maxJ, int maxI)
            {
                // Determine the size change needed
                const int extraRowsTop = std::max(m_minJ - minJ, 0);
                const int extraRowsBottom = std::max(maxJ - m_maxJ, 0);
                const int extraColsLeft = std::max(m_minI - minI, 0);
                const int extraColsRight = std::max(maxI - m_maxI, 0);

                if (extraRowsTop == 0 && extraRowsBottom == 0 && extraColsLeft == 0 && extraColsRight == 0)
                {
                    return;
                }

                // Create new matrix with the new size and initialize to missing value
                const auto newRows = static_cast<int>(m_matrix.rows()) + extraRowsTop + extraRowsBottom;
                const auto newCols = static_cast<int>(m_matrix.cols()) + extraColsLeft + extraColsRight;
                lin_alg::Matrix<int> newMatrix(newRows, newCols);
                newMatrix.setConstant(constants::missing::intValue);

                // Copy the existing matrix into the new one
                newMatrix.block(extraRowsTop, extraColsLeft, m_matrix.rows(), m_matrix.cols()) = m_matrix;

                // Update matrix and bounds
                m_matrix.swap(newMatrix);
                m_minI = std::min(m_minI, minI);
                m_minJ = std::min(m_minJ, minJ);
                m_maxI = std::max(m_maxI, maxI);
                m_maxJ = std::max(m_maxJ, maxJ);
            }

        private:
            lin_alg::Matrix<int> m_matrix = lin_alg::Matrix<int>(1, 1); ///< Underlying matrix storage.
            int m_minI = 0;                                             ///< Minimum column index.
            int m_minJ = 0;                                             ///< Minimum row index.
            int m_maxI = 0;                                             ///< Maximum column index.
            int m_maxJ = 0;                                             ///< Maximum row index.
        };

        /// @brief Computes the local mapping of the nodes composing the face
        [[nodiscard]] Eigen::Matrix<UInt, 2, 2> ComputeLocalNodeMapping(UInt face) const;

        /// @brief Computes the node indices of the neighbouring faces
        [[nodiscard]] UInt ComputeNeighbouringFaceNodes(const UInt face,
                                                        const Eigen::Matrix<UInt, 2, 2>& localNodeMapping,
                                                        const UInt d);

        /// @brief Checks the valid node and the candidate node are connected and are part of a quadrangular face
        [[nodiscard]] bool CheckGridLine(const UInt validNode, const UInt candidateNode) const;

        /// @brief Computes the final curvilinear matrix
        bool IsConnectionValid(const UInt candidateNode, const int iCandidate, const int jCandidate);

        /// @brief Computes the final curvilinear matrix
        [[nodiscard]] lin_alg::Matrix<Point> ComputeCurvilinearMatrix();

        Mesh2D& m_mesh;                     ///< The mesh to convert
        std::vector<int> m_i;               ///< The i indices of each node on the curvilinear grid
        std::vector<int> m_j;               ///< The j indices of each node on the curvilinear grid
        std::vector<bool> m_visitedFace;    ///< Mask for faces already visited in th breadth-first search
        std::vector<bool> m_convertedFaces; ///< The faces in the mesh converted to curvilinear

        const std::array<std::array<int, 2>, 4> m_nodeFrom = {{{0, 0},
                                                               {0, 0},
                                                               {1, 0},
                                                               {1, 1}}}; ///< starting edge node indices for each direction in the local mapping

        const std::array<std::array<int, 2>, 4> m_nodeTo = {{{0, 1},
                                                             {1, 0},
                                                             {1, 1},
                                                             {0, 1}}}; ///< ending edge node indices for each direction in the local mapping

        const std::array<std::array<int, 2>, 4> m_directionsDeltas = {{{-1, 0},
                                                                       {0, -1},
                                                                       {1, 0},
                                                                       {0, 1}}}; ///< increments for the new nodes depending on the node direction

        
        MatrixWithNegativeIndices m_mapping; ///< Unstructured node indices in the curvilinear grid
    };

} // namespace meshkernel
