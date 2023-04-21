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

namespace meshkernel
{
    /// @brief The interpolant types
    enum class Interpolants
    {
        AveragingInterpolation = 0,
        BilinearInterpolationOnGriddedSamples = 1,
    };

    /// @brief Interface for interpolation methods
    class MeshInterpolationInterface
    {
    public:
        /// @brief Virtual destructor
        virtual ~MeshInterpolationInterface() = default;

        /// @brief Compute
        virtual void Compute() = 0;

        /// @brief Gets the interpolation value at a specific node
        /// @param[in] node The node index
        /// @return The interpolated value
        [[nodiscard]] virtual double GetNodeResult(size_t node) const = 0;

        /// @brief Gets the interpolation value at a specific edge
        /// @param[in] edge The edge index
        /// @return The interpolated value
        [[nodiscard]] virtual double GetEdgeResult(size_t edge) const = 0;

        /// @brief Gets the interpolation value at a specific face
        /// @param[in] face The face index
        /// @return  The interpolated value
        [[nodiscard]] virtual double GetFaceResult(size_t face) const = 0;

        /// @brief Gets all interpolated values at nodes
        /// @return The interpolated values
        [[nodiscard]] virtual const std::vector<double>& GetNodeResults() const = 0;

        /// @brief Gets all interpolated values at edges
        /// @return The interpolated values
        [[nodiscard]] virtual const std::vector<double>& GetEdgeResults() const = 0;

        /// @brief Gets all interpolated values at faces
        /// @return The interpolated values
        [[nodiscard]] virtual const std::vector<double>& GetFaceResults() const = 0;
    };
} // namespace meshkernel
