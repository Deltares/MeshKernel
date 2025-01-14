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

#include "MeshKernel/UndoActions/ResetEdgeAction.hpp"
#include "MeshKernel/Mesh.hpp"

std::unique_ptr<meshkernel::ResetEdgeAction> meshkernel::ResetEdgeAction::Create(Mesh& mesh, const UInt id, const Edge& initial, const Edge& updated)
{
    return std::make_unique<ResetEdgeAction>(mesh, id, initial, updated);
}

meshkernel::ResetEdgeAction::ResetEdgeAction(Mesh& mesh, const UInt id, const Edge& initial, const Edge& updated) : BaseMeshUndoAction<ResetEdgeAction, Mesh>(mesh), m_edgeId(id), m_initialEdge(initial), m_updatedEdge(updated) {}

meshkernel::UInt meshkernel::ResetEdgeAction::EdgeId() const
{
    return m_edgeId;
}

const meshkernel::Edge& meshkernel::ResetEdgeAction::InitialEdge() const
{
    return m_initialEdge;
}

const meshkernel::Edge& meshkernel::ResetEdgeAction::UpdatedEdge() const
{
    return m_updatedEdge;
}
