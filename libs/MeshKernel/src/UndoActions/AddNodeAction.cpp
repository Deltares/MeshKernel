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

#include "MeshKernel/UndoActions/AddNodeAction.hpp"
#include "MeshKernel/Mesh.hpp"

std::unique_ptr<meshkernel::AddNodeAction> meshkernel::AddNodeAction::Create(Mesh& mesh, const UInt id, const Point& point)
{
    return std::make_unique<AddNodeAction>(mesh, id, point);
}

meshkernel::AddNodeAction::AddNodeAction(Mesh& mesh, const UInt id, const Point& p) : BaseMeshUndoAction<AddNodeAction, Mesh>(mesh), m_nodeId(id), m_node(p) {}

meshkernel::UInt meshkernel::AddNodeAction::NodeId() const
{
    return m_nodeId;
}

const meshkernel::Point& meshkernel::AddNodeAction::Node() const
{
    return m_node;
}
