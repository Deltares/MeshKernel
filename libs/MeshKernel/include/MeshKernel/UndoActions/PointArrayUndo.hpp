//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2025.
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

#include <cstdint>
#include <vector>

#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Point.hpp"
#include "MeshKernel/UndoActions/BaseMeshUndoAction.hpp"

namespace meshkernel
{
    /// @brief Forward declaration of the unstructured mesh
    class Mesh2D;

    /// @brief Action to save an array of points stored in the mesh
    ///
    /// The undo will simply swap these values
    class PointArrayUndo : public BaseMeshUndoAction<PointArrayUndo, Mesh2D>
    {
    public:
        /// @brief Allocate a DeleteNodeAction and return a unique_ptr to the newly created object.
        static std::unique_ptr<PointArrayUndo> Create(Mesh2D& mesh, const std::vector<Point>& nodes);

        /// @brief Constructor
        explicit PointArrayUndo(Mesh2D& mesh, const std::vector<Point>& nodes);

        /// @brief Swap the array of points with the saved values.
        void Swap(std::vector<Point>& nodes);

        /// \brief Compute the approximate amount of memory being used, in bytes.
        std::uint64_t MemorySize() const override;

    private:
        /// @brief The saved point values.
        std::vector<Point> m_savedNodes;
    };

} // namespace meshkernel
