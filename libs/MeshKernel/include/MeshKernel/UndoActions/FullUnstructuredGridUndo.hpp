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
#include <vector>

#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Point.hpp"
#include "MeshKernel/UndoActions/BaseMeshUndoAction.hpp"

namespace meshkernel
{
    /// @brief Forward declaration of the unstructured mesh
    class Mesh;

    /// @brief Action to save all the node and edge values from a mesh
    ///
    /// The undo will simply swap these values
    class FullUnstructuredGridUndo : public BaseMeshUndoAction<FullUnstructuredGridUndo, Mesh>
    {
    public:
        /// @brief Allocate a DeleteNodeAction and return a unique_ptr to the newly created object.
        static std::unique_ptr<FullUnstructuredGridUndo> Create(Mesh& mesh);

        /// @brief Constructor
        explicit FullUnstructuredGridUndo(Mesh& mesh);

        /// @brief Swap the nodes and edges with the saved values.
        void Swap(std::vector<Point>& nodes, std::vector<Edge>& edges);

        /// \brief Compute the approximate amount of memory being used, in bytes.
        std::uint64_t MemorySize() const override;

    private:
        /// @brief The saved node values.
        std::vector<Point> m_savedNodes;

        /// @brief The saved edge values.
        std::vector<Edge> m_savedEdges;
    };

} // namespace meshkernel
