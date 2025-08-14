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

#include <array>
#include <vector>

#include "MeshKernel/Definitions.hpp"

namespace meshkernel
{
    class Mesh2D;

    /// <summary>
    /// Orthogonalizion (optimize the aspect ratios) and mesh smoothing (optimize internal face angles or area).
    /// </summary>
    class Orthogonalizer
    {

    public:
        /// @brief Constructor
        explicit Orthogonalizer(const Mesh2D& mesh,
                                const std::vector<std::vector<UInt>>& nodesNodes,
                                const std::vector<MeshNodeType>& nodeTypes,
                                const std::vector<Point>& faceCircumCentres);

        /// @brief Computes the smoother weights and the right hand side
        void Compute();

        /// @brief Gets the weight for a certain node and connected node
        /// @brief node
        /// @brief connectedNode
        /// @returns The contribution of orthogonalizer to the left hand side, the linear system
        [[nodiscard]] double GetWeight(UInt node, UInt connectedNode) const
        {
            return m_weights[node][connectedNode];
        }

        /// @brief Gets the weight for a certain node and connected node
        /// @brief node
        /// @brief connectedNode
        /// @returns The contribution of orthogonalizer to the right hand size
        [[nodiscard]] double GetRightHandSide(UInt node, UInt connectedNode) const
        {
            return m_rhs[node][connectedNode];
        }

    private:
        /// @brief Computes the aspect ratio of each edge (orthonet_compute_aspect)
        /// @brief mesh
        /// @returns If the method succeeded
        bool AspectRatio(const Mesh2D& mesh);

        const Mesh2D& m_mesh;                               ///< A reference to mesh
        const std::vector<std::vector<UInt>>& m_nodesNodes; ///< Node-node connectivity
        const std::vector<MeshNodeType>& m_nodeType;        ///< Type of each node
        const std::vector<Point>& m_faceCircumCentres;      ///< The circumcentre of each face in the mesh
        std::vector<double> m_aspectRatios;                 ///< Aspect ratios
        std::vector<std::vector<double>> m_weights;         ///< Weights
        std::vector<std::array<double, 2>> m_rhs;           ///< Right hand side
    };
} // namespace meshkernel
