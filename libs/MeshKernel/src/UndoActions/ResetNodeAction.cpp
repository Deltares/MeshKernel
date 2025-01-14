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

#include "MeshKernel/UndoActions/ResetNodeAction.hpp"
#include "MeshKernel/Mesh.hpp"

std::unique_ptr<meshkernel::ResetNodeAction> meshkernel::ResetNodeAction::Create(Mesh& mesh, const UInt id, const Point& initial, const Point& updated)
{
    return std::make_unique<ResetNodeAction>(mesh, id, initial, updated);
}

meshkernel::ResetNodeAction::ResetNodeAction(Mesh& mesh, const UInt id, const Point& initial, const Point& updated) : BaseMeshUndoAction<ResetNodeAction, Mesh>(mesh), m_nodeId(id), m_initialNode(initial), m_updatedNode(updated) {}

meshkernel::UInt meshkernel::ResetNodeAction::NodeId() const
{
    return m_nodeId;
}

const meshkernel::Point& meshkernel::ResetNodeAction::InitialNode() const
{
    return m_initialNode;
}

const meshkernel::Point& meshkernel::ResetNodeAction::UpdatedNode() const
{
    return m_updatedNode;
}
