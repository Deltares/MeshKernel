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
#include "MeshKernel/Vector.hpp"

namespace meshkernel
{
    /// @brief Forward declaration of the unstructured mesh
    class Mesh;

    /// @brief Action to move a node in an unstructured mesh.
    class MoveNodeAction : public BaseMeshUndoAction<MoveNodeAction, Mesh>
    {
    public:
        /// @brief Simple struct, combining the id of the node to be displaced and the displacement vector
        struct NodeDisplacement
        {
            /// @brief Displacement vector.
            Vector m_displacement;

            /// @brief Id of the node to be displaced.
            UInt m_nodeId;
        };

        /// @brief Typedef for array of node displacements
        using DisplacementVector = std::vector<NodeDisplacement>;

        /// @brief Const iterator over the displacement vector.
        using const_iterator = DisplacementVector::const_iterator;

        /// @brief Allocate a MoveNodeAction and return a unique_ptr to the newly create object.
        static std::unique_ptr<MoveNodeAction> Create(Mesh& mesh);

        /// @brief Constructor
        MoveNodeAction(Mesh& mesh);

        /// @brief Set the amount of displacement for node
        void AddDisplacement(const UInt nodeId, const double xDisplacement, const double yDisplacement);

        // TODO How best to apply and and revert displacements?
        // 1. iterate over internals, the mesh can then use use it as it wants
        // 2. implement 2 functions apply and revert taking node array and apply (or revert) displacements
        //    then all implementaiton details can be hidden

        /// @brief Return iterator pointing to the first element of the displacement vector
        const_iterator begin() const;

        /// @brief Return iterator pointing to one-past-the-end element of the displacement vector
        const_iterator end() const;

        /// @brief Print the move node action to the stream
        void Print(std::ostream& out = std::cout) const override;

    private:
        /// @brief Vector of node displacements.
        DisplacementVector m_displacements;
    };

} // namespace meshkernel

inline meshkernel::MoveNodeAction::const_iterator meshkernel::MoveNodeAction::begin() const
{
    return m_displacements.begin();
}

inline meshkernel::MoveNodeAction::const_iterator meshkernel::MoveNodeAction::end() const
{
    return m_displacements.end();
}
