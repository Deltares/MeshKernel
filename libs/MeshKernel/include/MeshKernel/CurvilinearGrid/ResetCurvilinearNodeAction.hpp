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

#include "MeshKernel/CurvilinearGrid/CurvilinearGridNodeIndices.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/UndoActions/BaseMeshUndoAction.hpp"

namespace meshkernel
{
    /// @brief Forward declaration of the curvilinear grid
    class CurvilinearGrid;

    /// @brief Action to add an node to an unstructured mesh.
    class ResetCurvilinearNodeAction : public BaseMeshUndoAction<ResetCurvilinearNodeAction, CurvilinearGrid>
    {
    public:
        /// @brief Return unique pointer to newly created ResetCurvilinearNodeAction object
        static std::unique_ptr<ResetCurvilinearNodeAction> Create(CurvilinearGrid& grid,
                                                                  const CurvilinearGridNodeIndices nodeId,
                                                                  const Point& initial,
                                                                  const Point& updated,
                                                                  const bool recalculateNodeTypes = false);

        /// @brief Constructor
        ResetCurvilinearNodeAction(CurvilinearGrid& grid,
                                   const CurvilinearGridNodeIndices nodeId,
                                   const Point& initial,
                                   const Point& updated,
                                   const bool recalculateNodeTypes = false);

        /// @brief Get the node identifier
        CurvilinearGridNodeIndices NodeId() const;

        /// @brief Get the initial node
        const Point& InitialNode() const;

        /// @brief Get the initial node
        const Point& UpdatedNode() const;

        /// @brief Get recalculateNodeTypes indicator
        bool RecalculateNodeTypes() const;

        /// @brief Print the reset node action to the stream
        void Print(std::ostream& out = std::cout) const override;

    private:
        /// @brief The node identifier
        CurvilinearGridNodeIndices m_nodeId;

        /// @brief The initial node
        Point m_initialNode;

        /// @brief The updated node
        Point m_updatedNode;

        /// @brief Indicate if the node types need to be recomputed after undo/redo
        bool m_recalculateNodeTypes;
    };

} // namespace meshkernel
