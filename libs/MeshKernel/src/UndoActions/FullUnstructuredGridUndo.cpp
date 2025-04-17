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

#include "MeshKernel/UndoActions/FullUnstructuredGridUndo.hpp"
#include "MeshKernel/Mesh.hpp"

std::unique_ptr<meshkernel::FullUnstructuredGridUndo> meshkernel::FullUnstructuredGridUndo::Create(Mesh& mesh)
{
    return std::make_unique<FullUnstructuredGridUndo>(mesh);
}

meshkernel::FullUnstructuredGridUndo::FullUnstructuredGridUndo(Mesh& mesh) : BaseMeshUndoAction<FullUnstructuredGridUndo, Mesh>(mesh), m_savedNodes(mesh.Nodes()), m_savedEdges(mesh.Edges()) {}

void meshkernel::FullUnstructuredGridUndo::Swap(std::vector<Point>& nodes, std::vector<Edge>& edges)
{
    std::swap(nodes, m_savedNodes);
    std::swap(edges, m_savedEdges);
}

std::uint64_t meshkernel::FullUnstructuredGridUndo::MemorySize() const
{
    return sizeof(FullUnstructuredGridUndo) + sizeof(Point) * m_savedNodes.capacity() + sizeof(Edge) * m_savedEdges.capacity();
}
