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
#include <memory>
#include <utility>
#include <vector>

#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/UndoActions/CompoundUndoAction.hpp"
#include "MeshKernel/UndoActions/UndoAction.hpp"

namespace meshkernel
{
    /// @brief Split the row or column connected to the given edge.
    class SplitRowColumnOfMesh final
    {
    public:
        /// @brief Split the row or column connected to the given edge.
        ///
        /// The splitting will occur upto either the boundary (next elememnt is null value),
        /// when the next element is not a quadrilateral, or
        /// when the next element is the same as the first element, i.e. a loop has been detected
        [[nodiscard]] std::unique_ptr<UndoAction> Compute(Mesh2D& mesh, const UInt edgeId) const;

    private:
        /// @brief Split an edge in the middle, into two edges half the size returning the ID of the new node.
        ///
        /// The two half size edges are each connected to the new (mid point) node and the
        /// node on the respective end of the original edge.
        UInt SplitEdge(Mesh2D& mesh, const UInt edgeId, std::vector<UInt>& edgesToDelete, CompoundUndoAction& undoActions) const;

        /// @brief Split the element
        void SplitElement(Mesh2D& mesh,
                          const UInt elementId,
                          const UInt edgeId,
                          UInt& newNode,
                          CompoundUndoAction& undoActions,
                          std::vector<UInt>& edgesToDelete) const;

        /// @brief Split the first element of a loop of elements
        void SplitFirstLoopElement(Mesh2D& mesh,
                                   const UInt elementId,
                                   const UInt edgeId,
                                   UInt& firstNode,
                                   UInt& secondNode,
                                   CompoundUndoAction& undoActions,
                                   std::vector<UInt>& edgesToDelete) const;

        //// @brief Split the elements and edges along the row or column
        void SplitAlongRow(Mesh2D& mesh,
                           const std::vector<UInt>& elementIds,
                           const std::vector<UInt>& edgeIds,
                           CompoundUndoAction& undoActions,
                           std::vector<UInt>& edgesToDelete) const;

        /// @brief Get the element along the opposite edge
        UInt GetNextElement(const Mesh2D& mesh, const UInt elementId, const UInt edgeId) const;

        /// @brief Get the ID of the next edge
        UInt OppositeEdgeId(const Mesh2D& mesh, const UInt elementId, const UInt edgeId) const;

        /// @brief Determine is the edge is a valid edge or not.
        bool IsValidEdge(const Mesh2D& mesh, const UInt edgeId) const;

        /// @brief Determine is the element is a quadrilateral or not
        bool IsQuadrilateral(const Mesh2D& mesh, const UInt elementId) const;

        /// @brief Determine if it may be possible to split the edge
        ///
        /// If two elements are attached to edge then either must be a quadrilateral,
        /// otherwise the single attached elememt must be quadrilateral.
        bool MayBeSplit(const Mesh2D& mesh, const UInt edgeId) const;

        /// @brief Collect the ID's of the edges and elements to be split.
        void CollectElementsToSplit(const Mesh2D& mesh, const UInt edgeId, std::vector<UInt>& elementIds, std::vector<UInt>& edgeIds) const;

        /// @brief Get the next edge.
        void GetNextEdge(const Mesh2D& mesh, UInt& elementId, UInt& edgeId) const;
    };

} // namespace meshkernel
