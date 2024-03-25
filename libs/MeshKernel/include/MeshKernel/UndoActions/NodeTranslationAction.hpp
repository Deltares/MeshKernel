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

#include "MeshKernel/Point.hpp"
#include "MeshKernel/UndoActions/BaseMeshUndoAction.hpp"

namespace meshkernel
{
    /// @brief Forward declaration of the unstructured mesh
    class Mesh;

    /// @brief Undo action for all operations that translate nodes
    ///
    /// \note This does not keep track of any changes in edge information
    class NodeTranslationAction : public BaseMeshUndoAction<NodeTranslationAction, Mesh>
    {
    public:
        /// @brief Allocate a ResetNodeAction and return a unique_ptr to the newly create object.
        ///
        /// Assumes all nodes are to be translated.
        /// \note A NodeTranslationAction needs to be constructed with the untranslated node values.
        static std::unique_ptr<NodeTranslationAction> Create(Mesh& mesh);

        /// @brief Allocate a ResetNodeAction and return a unique_ptr to the newly create object.
        ///
        /// \note A NodeTranslationAction needs to be constructed with the untranslated node values.
        static std::unique_ptr<NodeTranslationAction> Create(Mesh& mesh, const std::vector<UInt>& nodeIndices);

        /// @brief Constructor
        explicit NodeTranslationAction(Mesh& mesh);

        /// @brief Constructor
        NodeTranslationAction(Mesh& mesh, const std::vector<UInt>& nodeIndices);

        /// @brief Swap the stored nodes with those given.
        void Swap(std::vector<Point>& nodes);

        /// @brief Get the number of bytes used by this object.
        std::uint64_t MemorySize() const override;

        /// @brief Print the reset node action to the stream
        void Print(std::ostream& out = std::cout) const override;

    protected:
        /// @brief Get the number of nodes
        UInt NumberOfNodes() const;

    private:
        /// @brief Set of nodes moved during the translation procedure.
        ///
        /// Value of the node depends on the state of the action:
        /// 1. When ActionState::Committed: m_nodes contain mesh node state before translation
        /// 2. When ActionState::Restored:  m_nodes contain mesh node state after translation
        std::vector<Point> m_nodes;

        /// @brief Indices of nodes moved during the translation procedure.
        ///
        /// This array may be empty if all nodes have been moved.
        std::vector<UInt> m_nodeIndices;
    };

} // namespace meshkernel
