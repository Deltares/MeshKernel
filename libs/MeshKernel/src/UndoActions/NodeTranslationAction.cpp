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

#include "MeshKernel/UndoActions/NodeTranslationAction.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Formatting.hpp"
#include "MeshKernel/Mesh.hpp"

#include <algorithm>
#include <ranges>

std::unique_ptr<meshkernel::NodeTranslationAction> meshkernel::NodeTranslationAction::Create(Mesh& mesh)
{
    return std::make_unique<NodeTranslationAction>(mesh);
}

std::unique_ptr<meshkernel::NodeTranslationAction> meshkernel::NodeTranslationAction::Create(Mesh& mesh, const std::vector<UInt>& nodeIndices)
{
    return std::make_unique<NodeTranslationAction>(mesh, nodeIndices);
}

meshkernel::NodeTranslationAction::NodeTranslationAction(Mesh& mesh) : NodeTranslationAction(mesh, std::vector<UInt>()) {}

meshkernel::NodeTranslationAction::NodeTranslationAction(Mesh& mesh, const std::vector<UInt>& nodeIndices)
    : BaseMeshUndoAction<NodeTranslationAction, Mesh>(mesh)
{

    if (nodeIndices.size() > mesh.GetNumNodes())
    {
        throw ConstraintError("Number of node indices is greater than the number of nodes in the mesh, {} > {}",
                              nodeIndices.size(), mesh.GetNumNodes());
    }

    if (nodeIndices.empty() || nodeIndices.size() * (sizeof(UInt) + sizeof(Point)) > mesh.GetNumNodes() * sizeof(Point))
    {
        // save all nodes and only the nodes. No need to save the node indices
        // Based on approximation of total memory required for the action.
        m_nodes = mesh.Nodes();
    }
    else
    {
        // save only nodes to be moved.
        m_nodes.resize(nodeIndices.size());
        m_nodeIndices = nodeIndices;

        // Copy nodes to be saved from mesh.
        std::ranges::transform(nodeIndices, m_nodes.begin(), [&m = mesh](UInt pos)
                               { return m.Node(pos); });
        // std::transform(nodeIndices.begin(), nodeIndices.end(), m_nodes.begin(), [&m = mesh](UInt pos)
        //                { return m.Node(pos); });
    }
}

meshkernel::UInt meshkernel::NodeTranslationAction::NumberOfNodes() const
{
    return static_cast<UInt>(m_nodes.size());
}

void meshkernel::NodeTranslationAction::Swap(std::vector<Point>& nodes)
{
    if (nodes.size() < m_nodes.size())
    {
        throw ConstraintError("Number of nodes passed is less than nodes stored. {} < {}",
                              nodes.size(), m_nodes.size());
    }

    if (m_nodeIndices.empty())
    {
        std::swap_ranges(m_nodes.begin(), m_nodes.end(), nodes.begin());
    }
    else
    {
        for (UInt i = 0; i < m_nodeIndices.size(); ++i)
        {
            std::swap(m_nodes[i], nodes[m_nodeIndices[i]]);
        }
    }
}

std::uint64_t meshkernel::NodeTranslationAction::MemorySize() const
{
    return sizeof(*this) + m_nodes.capacity() * sizeof(Point) + m_nodeIndices.capacity() * sizeof(UInt);
}
