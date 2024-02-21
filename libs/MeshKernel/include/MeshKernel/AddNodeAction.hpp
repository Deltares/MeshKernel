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
#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Point.hpp"

namespace meshkernel
{
    /// @brief Forward declaration of the unstructured mesh
    class Mesh;

    /// @brief Action to add a node to an unstructured mesh.
    class AddNodeAction : public BaseMeshUndoAction<AddNodeAction, Mesh>
    {
    public:
        /// @brief Allocate a AddNodeAction and return a unique_ptr to the newly create object.
        static std::unique_ptr<AddNodeAction> Create(Mesh& mesh, const UInt id, const Point& point);

        /// @brief Constructor
        AddNodeAction(Mesh& mesh, const UInt id, const Point& p);

        /// @brief Get the node identifier
        UInt NodeId() const;

        /// @brief Get the node location
        const Point& Node() const;

        /// @brief Print the add node action to the stream
        void Print(std::ostream& out = std::cout) const override;

    private:
        /// @brief The node identifier
        UInt m_nodeId;

        /// @brief The added node location
        Point m_node;
    };

} // namespace meshkernel
