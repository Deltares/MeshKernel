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

#include "MeshKernel/UndoActions/MeshConversionAction.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Mesh.hpp"

#include <algorithm>

std::unique_ptr<meshkernel::MeshConversionAction> meshkernel::MeshConversionAction::Create(Mesh& mesh)
{
    return std::make_unique<MeshConversionAction>(mesh);
}

meshkernel::MeshConversionAction::MeshConversionAction(Mesh& mesh)
    : BaseMeshUndoAction<MeshConversionAction, Mesh>(mesh), m_nodes(mesh.Nodes()), m_projection(mesh.m_projection) {}

void meshkernel::MeshConversionAction::Swap(std::vector<Point>& nodes, Projection& projection)
{
    if (nodes.size() < m_nodes.size())
    {
        throw ConstraintError("Number of nodes passed is less than nodes stored. {} < {}",
                              nodes.size(), m_nodes.size());
    }

    std::swap_ranges(m_nodes.begin(), m_nodes.end(), nodes.begin());
    std::swap(m_projection, projection);
}

std::uint64_t meshkernel::MeshConversionAction::MemorySize() const
{
    return sizeof(*this) + m_nodes.capacity() * sizeof(Point);
}
