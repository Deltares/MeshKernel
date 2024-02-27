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

#include "MeshKernel/BaseMeshUndoAction.hpp"
#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Point.hpp"

namespace meshkernel
{
    /// @brief Forward declaration of the unstructured mesh
    class Mesh2D;

    /// @brief Action to add a node to an unstructured mesh.
    class SphericalCoordinatesOffsetAction : public BaseMeshUndoAction<SphericalCoordinatesOffsetAction, Mesh2D>
    {
    public:
        /// @brief Typedef of const iterator.
        using const_iterator = std::vector<UInt>::const_iterator;

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

        /// @brief Return iterator pointing to the first element of the nodes
        const_iterator BeginDecrease() const;

        /// @brief Return iterator pointing to one-past-the-end element of the
        const_iterator EndDecrease() const;

        /// @brief Return iterator pointing to the first element of the nodes
        const_iterator BeginIncrease() const;

        /// @brief Return iterator pointing to one-past-the-end element of the
        const_iterator EndIncrease() const;

        /// \brief Compute the approximate amount of memory being used, in bytes.
        std::uint64_t MemorySize() const override;

        /// @brief Print the add node action to the stream
        void Print(std::ostream& out = std::cout) const override;

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

inline double meshkernel::SphericalCoordinatesOffsetAction::MinX() const
{
    return m_xMin;
}

inline double meshkernel::SphericalCoordinatesOffsetAction::MaxX() const
{
    return m_xMax;
}

inline meshkernel::SphericalCoordinatesOffsetAction::const_iterator meshkernel::SphericalCoordinatesOffsetAction::BeginDecrease() const
{
    return m_offsetNodesDecrease.begin();
}

inline meshkernel::SphericalCoordinatesOffsetAction::const_iterator meshkernel::SphericalCoordinatesOffsetAction::EndDecrease() const
{
    return m_offsetNodesDecrease.end();
}

inline meshkernel::SphericalCoordinatesOffsetAction::const_iterator meshkernel::SphericalCoordinatesOffsetAction::BeginIncrease() const
{
    return m_offsetNodesIncrease.begin();
}

inline meshkernel::SphericalCoordinatesOffsetAction::const_iterator meshkernel::SphericalCoordinatesOffsetAction::EndIncrease() const
{
    return m_offsetNodesIncrease.end();
}
