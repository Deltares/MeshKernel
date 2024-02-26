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
#include "MeshKernel/Entities.hpp"

namespace meshkernel
{
    /// @brief Forward declaration of the unstructured mesh
    class Mesh2D;

    /// @brief Action to add an node to an unstructured mesh.
    class OrthogonalizationAndSmoothingAction : public BaseMeshUndoAction<OrthogonalizationAndSmoothingAction, Mesh2D>
    {
    public:
        /// @brief Allocate a ResetNodeAction and return a unique_ptr to the newly create object.
        static std::unique_ptr<OrthogonalizationAndSmoothingAction> Create(Mesh2D& mesh, const std::vector<Point>& nodes, const std::vector<Edge>& edges);

        /// @brief Constructor
        OrthogonalizationAndSmoothingAction(Mesh2D& mesh, const std::vector<Point>& nodes, const std::vector<Edge>& edges);

        const std::vector<Point>& OriginalNodes() const;

        const std::vector<Edge>& OriginalEdges() const;

        void SetTransformed

            /// @brief Get the number of bytes used by this object.
            std::uint64_t
            MemorySize() const override;

        /// @brief Print the reset node action to the stream
        void Print(std::ostream& out = std::cout) const override;

    private:
        std::vector<Point> m_originalNodes;
        std::vector<Point> m_originalEdges;
    };

} // namespace meshkernel

inline const std::vector<Point>& meshkernel::OrthogonalizationAndSmoothingAction::OriginalNodes() const
{
    return m_originalNodes;
}

inline const std::vector<Edge>& meshkernel::OrthogonalizationAndSmoothingAction::OriginalEdges() const
{
    return m_originalEdges;
}
