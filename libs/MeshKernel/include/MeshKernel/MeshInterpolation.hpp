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

#include "MeshKernel/Constants.hpp"

namespace meshkernel
{
    /// @brief The interpolant types
    enum class Interpolants
    {
        AveragingInterpolation = 0,
        BilinearInterpolationOnGriddedSamples = 1,
    };

    /// @brief Interface for interpolation methods
    class MeshInterpolation
    {
    public:
        /// @brief Virtual destructor
        virtual ~MeshInterpolation() = default;

        /// @brief Compute
        virtual void Compute() = 0;

        /// @brief Gets the interpolation value at a specific node
        /// @param[in] node The node index
        /// @return The interpolated value
        [[nodiscard]] double GetNodeResult(Index node) const { return m_nodeResults[node]; }

        /// @brief Gets the interpolation value at a specific edge
        /// @param[in] edge The edge index
        /// @return The interpolated value
        [[nodiscard]] double GetEdgeResult(Index edge) const { return m_edgeResults[edge]; }

        /// @brief Gets the interpolation value at a specific face
        /// @param[in] face The face index
        /// @return  The interpolated value
        [[nodiscard]] double GetFaceResult(Index face) const { return m_faceResults[face]; }

        /// @brief Gets all interpolated values at nodes
        /// @return The interpolated values
        [[nodiscard]] const std::vector<double>& GetNodeResults() const { return m_nodeResults; }

        /// @brief Gets all interpolated values at edges
        /// @return The interpolated values
        [[nodiscard]] const std::vector<double>& GetEdgeResults() const { return m_edgeResults; }

        /// @brief Gets all interpolated values at faces
        /// @return The interpolated values
        [[nodiscard]] const std::vector<double>& GetFaceResults() const { return m_faceResults; }

    protected:
        std::vector<double> m_nodeResults; ///< The interpolation results at nodes
        std::vector<double> m_edgeResults; ///< The interpolation results at edges
        std::vector<double> m_faceResults; ///< The interpolation results at faces
    };
} // namespace meshkernel
