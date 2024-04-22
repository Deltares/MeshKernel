//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2024.
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

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Point.hpp"
#include "MeshKernel/UndoActions/BaseMeshUndoAction.hpp"

namespace meshkernel
{
    /// @brief Forward declaration of the unstructured mesh
    class Mesh2D;

    /// @brief Action to add a node to an unstructured mesh.
    class SphericalCoordinatesOffsetAction : public BaseMeshUndoAction<SphericalCoordinatesOffsetAction, Mesh2D>
    {
    public:
        /// @brief Allocate a AddNodeAction and return a unique_ptr to the newly create object.
        static std::unique_ptr<SphericalCoordinatesOffsetAction> Create(Mesh2D& mesh, const double minx, const double maxx);

        /// @brief Constructor
        SphericalCoordinatesOffsetAction(Mesh2D& mesh, const double minx, const double maxx);

        /// @brief Get the minimum x-value
        double MinX() const;

        /// @brief Get the maximum x-value
        double MaxX() const;

        /// @brief Add index of node whose value is to be decreased
        void AddDecrease(const UInt nodeId);

        /// @brief Add index of node whose value is to be increased
        void AddIncrease(const UInt nodeId);

        /// @brief Apply the offset action to the nodes
        void ApplyOffset(std::vector<Point>& nodes) const;

        /// @brief Revert the offset action to the nodes
        void UndoOffset(std::vector<Point>& nodes) const;

        /// \brief Compute the approximate amount of memory being used, in bytes.
        std::uint64_t MemorySize() const override;

    private:
        /// @brief Minimum x-value
        double m_xMin;

        /// @brief Maximum x-value
        double m_xMax;

        /// @brief List of identifiers for nodes greater than minimum + 360 deg;
        std::vector<UInt> m_offsetNodesDecrease;

        /// @brief List of identifiers of nodes less that the minimum-x
        std::vector<UInt> m_offsetNodesIncrease;
    };

} // namespace meshkernel
