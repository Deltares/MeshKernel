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

#pragma once
#include "MeshKernel/Mesh1D.hpp"
#include <MeshKernel/Exceptions.hpp>

#include <MeshKernel/Entities.hpp>
#include <vector>

meshkernel::Mesh1D::Mesh1D(const std::vector<Edge>& edges,
                           const std::vector<Point>& nodes,
                           Projection projection) : Mesh(edges, nodes, projection){};

size_t meshkernel::Mesh1D::Find1dNodeCloseToAPoint(Point point, std::vector<bool> oneDNodeMask)
{
    if (GetNumNodes() <= 0)
    {
        throw std::invalid_argument("Mesh1D::Find1dNodeCloseToAPoint: There are no valid nodes.");
    }

    // create rtree a first time
    if (m_nodesRTree.Empty())
    {
        m_nodesRTree.BuildTree(m_nodes);
        m_nodesRTreeRequiresUpdate = false;
    }

    m_nodesRTree.NearestNeighbors(point);
    const auto resultSize = m_nodesRTree.GetQueryResultSize();

    for (auto index = 0; index < resultSize; ++index)
    {
        const auto nodeIndex = m_nodesRTree.GetQueryResult(0);
        if (oneDNodeMask[nodeIndex])
        {
            return nodeIndex;
        }
    }

    throw AlgorithmError("Mesh1D::Find1dNodeCloseToAPoint: Could not find the node index close to a point.");
}
