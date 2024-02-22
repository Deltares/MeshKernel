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
    class SphericalCoordinateOffsetAction : public BaseMeshUndoAction<AddNodeAction, Mesh2D>
    {
    public:
        /// @brief Allocate a AddNodeAction and return a unique_ptr to the newly create object.
        static std::unique_ptr<SphericalCoordinateOffsetAction> Create(Mesh2D& mesh, const double minx, const double maxx);

        /// @brief Constructor
        SphericalCoordinateOffsetAction(Mesh& mesh, const double minx, const double maxx);

        /// @brief Get the minimum x-value
        double MinX() const;

        /// @brief Get the maximum x-value
        double MaxX() const;

        /// @brief Return iterator pointing to the first element of the displacement vector
        const_iterator begin() const;

        /// @brief Return iterator pointing to one-past-the-end element of the displacement vector
        const_iterator end() const;

        /// \brief Compute the approximate amount of memory being used, in bytes.
        std::uint64_t MemorySize() const override;

        /// @brief Print the add node action to the stream
        void Print(std::ostream& out = std::cout) const override;

    private:
        /// @brief Minimum x-value
        double m_xMin;

        /// @brief Maximum x-value
        double m_xMax;

        /// @brief List of identifiers of transformed node.
        std::vector<UInt> m_offsetNodes;
    };

} // namespace meshkernel

inline double meshkernel::SphericalCoordinateOffsetAction::MinX() const
{
    return m_xMin;
}

inline double meshkernel::SphericalCoordinateOffsetAction::MaxX() const
{
    return m_xMax;
}

inline meshkernel::SphericalCoordinateOffsetAction::const_iterator meshkernel::SphericalCoordinateOffsetAction::begin() const
{
    return m_offsetNode.begin();
}

inline meshkernel::SphericalCoordinateOffsetAction::const_iterator meshkernel::SphericalCoordinateOffsetAction::end() const
{
    return m_offsetNode.end();
}
