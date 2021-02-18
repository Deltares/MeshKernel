//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2021.
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

#include "MeshKernel/Mesh1D.hpp"

#include <MeshKernel/Entities.hpp>
#include <vector>

meshkernel::Mesh1D::Mesh1D(const std::vector<Edge>& edges,
                           const std::vector<Point>& nodes,
                           Projection projection) : Mesh(edges, nodes, projection){};

void meshkernel::Mesh1D::SetFlatCopies()
{

    m_nodex.resize(GetNumNodes());
    m_nodey.resize(GetNumNodes());

    for (auto n = 0; n < GetNumNodes(); n++)
    {
        m_nodex[n] = m_nodes[n].x;
        m_nodey[n] = m_nodes[n].y;
    }

    size_t edgeIndex = 0;
    m_edgeNodes.resize(GetNumEdges() * 2);
    for (auto e = 0; e < GetNumEdges(); e++)
    {
        m_edgeNodes[edgeIndex] = static_cast<int>(m_edges[e].first);
        edgeIndex++;
        m_edgeNodes[edgeIndex] = static_cast<int>(m_edges[e].second);
        edgeIndex++;
    }

    // we always need to provide pointers to not empty memory
    if (m_nodex.empty())
    {
        m_nodex.resize(1);
    }
    if (m_nodey.empty())
    {
        m_nodey.resize(1);
    }
    if (m_edgeNodes.empty())
    {
        m_edgeNodes.resize(1);
    }
}
