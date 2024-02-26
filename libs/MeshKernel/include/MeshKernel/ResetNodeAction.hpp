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
    class Mesh;

    /// @brief Action to add an node to an unstructured mesh.
    class ResetNodeAction : public BaseMeshUndoAction<ResetNodeAction, Mesh>
    {
    public:
        /// @brief Allocate a ResetNodeAction and return a unique_ptr to the newly create object.
        static std::unique_ptr<ResetNodeAction> Create(Mesh& mesh, const UInt id, const Point& initial, const Point& updated);

        /// @brief Constructor
        ResetNodeAction(Mesh& mesh, const UInt id, const Point& initial, const Point& updated);

        /// @brief Get the node identifier
        UInt NodeId() const;

        /// @brief Get the initial node
        const Point& InitialNode() const;

        /// @brief Get the initial node
        const Point& UpdatedNode() const;

        /// @brief Print the reset node action to the stream
        void Print(std::ostream& out = std::cout) const override;

    private:
        /// @brief The node identifier
        UInt m_nodeId;

        /// @brief The initial node
        Point m_initialNode;

        /// @brief The updated node
        Point m_updatedNode;
    };

} // namespace meshkernel
