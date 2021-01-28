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

namespace meshkernel
{
    class Mesh2D;

    /// <summary>
    /// Orthogonalizion (optimize the aspect ratios) and mesh smoothing (optimize internal face angles or area).
    /// </summary>
    class Orthogonalizer
    {

    public:
        /// @brief Ctor
        /// @returns
        explicit Orthogonalizer(std::shared_ptr<Mesh2D> mesh);

        /// @brief Computes the smoother weights and the right hans side
        void Compute();

        /// @brief Gets the weight for a certain node and connected node
        /// @brief node
        /// @brief connectedNode
        /// @returns The contribution of orthogonalizer to the left hand side, the linear system
        [[nodiscard]] double GetWeight(size_t node, size_t connectedNode)
        {
            return m_weights[node][connectedNode];
        }

        /// @brief Gets the weight for a certain node and connected node
        /// @brief node
        /// @brief connectedNode
        /// @returns The contribution of orthogonalizer to the right hand size
        [[nodiscard]] double GetRightHandSide(size_t node, size_t connectedNode)
        {
            return m_rhs[node][connectedNode];
        }

    private:
        /// @brief Computes the aspect ratio of each edge (orthonet_compute_aspect)
        /// @brief mesh
        /// @returns If the method succeeded
        bool AspectRatio(const Mesh2D& mesh);

        std::shared_ptr<Mesh2D> m_mesh;             ///< Pointer to mesh
        std::vector<double> m_aspectRatios;         ///< Aspect ratios
        std::vector<std::vector<double>> m_weights; ///< Weights
        std::vector<std::vector<double>> m_rhs;     ///< Right hand side
    };
} // namespace meshkernel
