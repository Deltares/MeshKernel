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

using namespace meshkernel::constants;

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
        class CurviMask
        {
        public:
            int getValue(int i, int j) const
            {
                const auto jIndex = j - m_minJ;
                const auto iIndex = i - m_minI;
                const auto value = m_matrix(jIndex, iIndex);
                return value;
            }

            void setValue(int i, int j, int value)
            {
                const auto jIndex = j - m_minJ;
                const auto iIndex = i - m_minI;
                m_matrix(jIndex, iIndex) = value;
            }

            [[nodiscard]] int rows() const
            {
                return static_cast<int>(m_matrix.rows());
            }

            [[nodiscard]] int cols() const
            {
                return static_cast<int>(m_matrix.cols());
            }

            [[nodiscard]] int minI() const
            {
                return m_minI;
            }

            [[nodiscard]] int maxI() const
            {
                return m_maxI;
            }
            [[nodiscard]] int minJ() const
            {
                return m_minJ;
            }

            [[nodiscard]] int maxJ() const
            {
                return m_maxJ;
            }

            void resize(int minI, int minJ, int maxI, int maxJ)
            {
                const int numLeftColumnsToAdd = std::max(m_minI - minI, 0);
                const int numRightColumnsToAdd = std::max(maxI - m_maxI, 0);
                const int numBottomRowsToAdd = std::max(m_minJ - minJ, 0);
                const int numUpperRowsToAdd = std::max(maxJ - m_maxJ, 0);

                const int newRows = numBottomRowsToAdd + numUpperRowsToAdd;
                const int newCols = numLeftColumnsToAdd + numRightColumnsToAdd;

                if (newRows == 0 && newCols == 0)
                {
                    return;
                }

                // Create a new matrix with the new size and initialize it to uintValue
                lin_alg::Matrix<int> newMatrix(m_matrix.rows() + newRows, m_matrix.cols() + newCols);
                newMatrix.setConstant(missing::intValue);

                // Copy the existing matrix to the appropriate position in the new matrix
                newMatrix.block(numBottomRowsToAdd, numLeftColumnsToAdd, m_matrix.rows(), m_matrix.cols()) = m_matrix;

                // Update the matrix and the bounds
                m_matrix.swap(newMatrix);
                m_minI = std::min(m_minI, minI);
                m_minJ = std::min(m_minJ, minJ);
                m_maxI = std::max(m_maxI, maxI);
                m_maxJ = std::max(m_maxJ, maxJ);
            }

        private:
            lin_alg::Matrix<int> m_matrix = lin_alg::Matrix<int>(1, 1);
            int m_minI = 0;
            int m_minJ = 0;
            int m_maxI = 0;
            int m_maxJ = 0;
        };

        /// @brief Computes the local mapping of the nodes composing the face
        [[nodiscard]] Eigen::Matrix<UInt, 2, 2> ComputeLocalNodeMapping(UInt face) const;

        /// @brief Computes the node indices of the neighbouring faces
        [[nodiscard]] UInt ComputeNeighbouringFaceNodes(const UInt face,
                                                        const Eigen::Matrix<UInt, 2, 2>& localNodeMapping,
                                                        const UInt d,
                                                        const std::vector<bool>& visitedFace);

        [[nodiscard]] bool CheckGridLine(const UInt validNode, const UInt candidateNode) const;

        /// @brief Computes the final curvilinear matrix
        bool IsConnectionValid(const UInt candidateNode, const int iCandidate, const int jCandidate);

        /// @brief Computes the final curvilinear matrix
        [[nodiscard]] lin_alg::Matrix<Point> ComputeCurvilinearMatrix();

        Mesh2D& m_mesh;       ///< The mesh to convert
        std::vector<int> m_i; ///< The i indices of each node on the curvilinear grid
        std::vector<int> m_j; ///< The j indices of each node on the curvilinear grid

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

        const int n_maxNumRowsColumns = 1000000; ///< The maximum number of allowed rows or columns

        CurviMask m_mapping;
    };

} // namespace meshkernel
