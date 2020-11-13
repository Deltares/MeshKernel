//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2020.
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
    class Mesh;

    /// <summary>
    /// Orthogonalizion (optimize the aspect ratios) and and mesh smoothing (optimize internal face angles or area).
    /// </summary>
    class Orthogonalizer
    {

    public:
        /// <summary>
        /// Ctor
        /// </summary>
        /// <returns></returns>
        explicit Orthogonalizer(std::shared_ptr<Mesh> mesh);

        /// @brief Computes the smoother weights and the right hans side
        void Compute();

        /// <summary>
        /// Gets the weight for a certain node and connected node
        /// </summary>
        /// <param name="node"></param>
        /// <param name="connectedNode"></param>
        /// <returns>The contribution of orthogonalizer to the left hand side, the linear system</returns>
        [[nodiscard]] inline double GetWeight(int node, int connectedNode)
        {
            return m_weights[node][connectedNode];
        }

        /// <summary>
        /// Gets the weight for a certain node and connected node
        /// </summary>
        /// <param name="node"></param>
        /// <param name="connectedNode"></param>
        /// <returns>The contribution of orthogonalizer to the right hand size</returns>
        [[nodiscard]] inline double GetRightHandSide(int node, int connectedNode)
        {
            return m_rhs[node][connectedNode];
        }

    private:
        /// <summary>
        /// Computes the aspect ratio of each edge (orthonet_compute_aspect)
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns>If the method succeeded</returns>
        bool AspectRatio(const Mesh& mesh);

        std::shared_ptr<Mesh> m_mesh;
        std::vector<double> m_aspectRatios;
        std::vector<std::vector<double>> m_weights;
        std::vector<std::vector<double>> m_rhs;
    };
} // namespace meshkernel
