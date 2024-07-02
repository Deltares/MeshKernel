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
#include <array>
#include <utility>
#include <vector>

#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/UndoActions/CompoundUndoAction.hpp"
#include "MeshKernel/UndoActions/UndoAction.hpp"

namespace meshkernel
{
    class SplitRowColumnOfMesh3 final
    {
    public:
        [[nodiscard]] std::unique_ptr<UndoAction> Compute(Mesh2D& mesh, const UInt edgeId) const;

    private:
        struct SplittingInfo
        {
            UInt elementId;
            UInt localStartEdgeId;
            UInt localEndEdgeId;
            UInt startEdgeId;
            UInt endEdgeId;
        };

        void SplitAlongRow(Mesh2D& mesh, const UInt edgeId, CompoundUndoAction& undoActions) const;

        void SplitAlongRow(Mesh2D& mesh, const std::vector<UInt>& elementIds, const std::vector<UInt>& edgeIds, CompoundUndoAction& undoActions, std::vector<UInt>& edgesToDelete) const;

        void SplitEdges(Mesh2D& mesh, std::vector<UInt>& elementIds, std::vector<UInt>& edgeIds, CompoundUndoAction& undoActions) const;

        void SplitEdge(Mesh2D& mesh, UInt elementId, UInt edgeId, UInt& previousNewNode, CompoundUndoAction& undoActions) const;

        void SplitEdge2(Mesh2D& mesh, UInt elementId, UInt edgeId, UInt& previousNewNode, CompoundUndoAction& undoActions) const;

        void SplitEdge3(Mesh2D& mesh, UInt elementId, UInt edgeId, UInt& previousNewNode, CompoundUndoAction& undoActions, std::vector<UInt>& edgesToDelete) const;

        void SplitEdge4(Mesh2D& mesh, UInt elementId, UInt edgeId, UInt& newNode, CompoundUndoAction& undoActions, std::vector<UInt>& edgesToDelete) const;

        bool IsValidEdge(const Mesh2D& mesh, const UInt edgeId) const;

        /// @brief Determine if it may be possible to split the edge
        ///
        /// If two element attached to edge then both must be quadrilaterals,
        /// otherwise the sole attached elememt must be quadrilateral
        bool MayBeSplit(const Mesh2D& mesh, const UInt edgeId) const;

        void CollectElementIdsToSplit(const Mesh2D& mesh, const UInt edgeId, std::vector<UInt>& elementIds, std::vector<UInt>& edgeIds) const;

        void GetNextElement(const Mesh2D& mesh, UInt& edgeId, UInt& elementId) const;
    };

} // namespace meshkernel
