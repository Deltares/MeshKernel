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
#include <utility>
#include <vector>

#include "MeshKernel/CurvilinearGrid/CurvilinearGridBlock.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridNodeIndices.hpp"
#include "MeshKernel/UndoActions/BaseMeshUndoAction.hpp"

namespace meshkernel
{
    /// @brief Forward declaration of the curvilinear mesh
    class CurvilinearGrid;

    /// @brief Undo action for blocks of nodes in a curvilinear mesh.
    class CurvilinearGridBlockUndoAction : public BaseMeshUndoAction<CurvilinearGridBlockUndoAction, CurvilinearGrid>
    {
    public:
        /// @brief Return unique pointer to newly created CurvilinearGridBlockUndoAction object
        static std::unique_ptr<CurvilinearGridBlockUndoAction> Create(CurvilinearGrid& grid,
                                                                      const CurvilinearGridNodeIndices& startOffset,
                                                                      const CurvilinearGridNodeIndices& endOffset);

        /// @brief Return unique pointer to newly created CurvilinearGridBlockUndoAction object
        ///
        /// The entire grid is the block
        static std::unique_ptr<CurvilinearGridBlockUndoAction> Create(CurvilinearGrid& grid);

        /// @brief Constructor, node values are copied from the grid for the block specified
        CurvilinearGridBlockUndoAction(CurvilinearGrid& grid,
                                       const CurvilinearGridNodeIndices& startOffset,
                                       const CurvilinearGridNodeIndices& endOffset);

        /// @brief Swap the saved grid nodes with those from the mesh.
        void Swap(CurvilinearGrid& grid);

        /// \brief Compute the approximate amount of memory being used, in bytes.
        std::uint64_t MemorySize() const override;

    private:
        /// @brief The saved block of grid nodes.
        CurvilinearGridBlock m_block;
    };

} // namespace meshkernel
