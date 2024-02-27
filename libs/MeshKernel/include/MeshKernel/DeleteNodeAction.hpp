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
#include <vector>

#include "MeshKernel/BaseMeshUndoAction.hpp"
#include "MeshKernel/Constants.hpp"
#include "MeshKernel/DeleteEdgeAction.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Point.hpp"
#include "MeshKernel/UndoAction.hpp"

namespace meshkernel
{
    /// @brief Forward declaration of the unstructured mesh
    class Mesh;

    /// @brief Action to delete a node from an unstructured mesh.
    class DeleteNodeAction : public UndoAction
    {
    public:
        /// @brief Allocate a DeleteNodeAction and return a unique_ptr to the newly create object.
        static std::unique_ptr<DeleteNodeAction> Create(Mesh& mesh, const UInt id, const Point& node);

        /// @brief Constructor
        DeleteNodeAction(Mesh& mesh, const UInt id, const Point& node);

        void Add(std::unique_ptr<DeleteEdgeAction>&& action);

        /// @brief Get the node identifier
        UInt NodeId() const;

        /// @brief Get the node location
        const Point& Node() const;

        /// \brief Compute the approximate amount of memory being used, in bytes.
        std::uint64_t MemorySize() const override;

        /// @brief Print the delete node action to the stream
        void Print(std::ostream& out = std::cout) const override;

    private:
        /// @brief Commit the action of deleting a node and all connecting edges
        void DoCommit() override;

        /// @brief The action of restoring a node and all connecting edges
        void DoRestore() override;

        /// @brief The unstructured mesh from which the node will be deleted.
        Mesh& m_mesh;

        /// @brief The node identifier
        UInt m_nodeId;

        /// @brief The deleted node location
        Point m_node;

        /// @brief Sequence of delete edge actions.
        std::vector<std::unique_ptr<DeleteEdgeAction>> m_deletedEdges;
    };

} // namespace meshkernel
