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
        /// @brief Computes the local mapping of the nodes composing the face
        [[nodiscard]] Eigen::Matrix<UInt, 2, 2> ComputeLocalNodeMapping(UInt face) const;

        /// @brief Computes the node indices of the neighbouring faces
        [[nodiscard]] UInt ComputeNeighbouringFaceNodes(const UInt face,
                                                        const Eigen::Matrix<UInt, 2, 2>& localNodeMapping,
                                                        const UInt d,
                                                        std::vector<bool>& visitedFace);

        /// @brief Computes the final curvilinear matrix
        [[nodiscard]] lin_alg::Matrix<Point> ComputeCurvilinearMatrix();

        Mesh2D& m_mesh;       ///< The mesh to convert
        std::vector<int> m_i; ///< The i indices of each node on the curvilinear grid
        std::vector<int> m_j; ///< The j indices of each node on the curvilinear grid

        std::array<std::array<int, 2>, 4> m_nodeFrom = {{{0, 0},
                                                         {0, 0},
                                                         {1, 0},
                                                         {1, 1}}}; ///< The mesh where the edges should be found

        std::array<std::array<int, 2>, 4> m_nodeTo = {{{0, 1},
                                                       {1, 0},
                                                       {1, 1},
                                                       {0, 1}}}; ///< The mesh where the edges should be found

        std::array<std::array<int, 2>, 4> m_directionsDeltas = {{{-1, 0},
                                                                 {0, -1},
                                                                 {1, 0},
                                                                 {0, 1}}}; ///< The mesh where the edges should be found

        int n_maxNumRowsColumns = 1000000; ///< The mesh where the edges should be found
    };

} // namespace meshkernel
