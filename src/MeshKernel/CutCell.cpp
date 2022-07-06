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

#include "MeshKernel/CutCell.hpp"
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/RTree.hpp>

using namespace meshkernel;

CutCell::CutCell(std::shared_ptr<Mesh2D> mesh) : m_mesh(mesh)
{
};

std::vector<int> CutCell::ClassifyNodes(const std::vector<Point>& boundaryLines) const
{

    std::vector nodeMask(m_mesh->GetNumNodes(), NodeClasses::innerNodeFlag);

    // Mask all mesh nodes on the left of the boundary lines as inactive
    for (auto i = 0; i < boundaryLines.size() - 1; ++i)
    {
        for (auto n = 0; n < m_mesh->m_nodes.size(); ++n)
        {
            const bool isLeft = crossProduct(boundaryLines[i], boundaryLines[i + 1], boundaryLines[i], m_mesh->m_nodes[n], m_mesh->m_projection) > 0.0;
            if (isLeft)
            {
                nodeMask[n] = NodeClasses::inactiveFlag;  
            }
        }
    }

    // Mask all faces crossed by the boundary line
    std::vector faceMask(m_mesh->GetNumFaces(), false);
    for (auto i = 0; i < boundaryLines.size() - 1; ++i)
    {
        for (auto e = 0; e < m_mesh->GetNumEdges(); ++e)
        {

            Point intersectionPoint;
            double crossProductValue;
            double ratioFirstSegment;
            double ratioSecondSegment;
            const auto areSegmentsCrossing = AreSegmentsCrossing(boundaryLines[i],
                                                                 boundaryLines[i + 1],
                                                                 m_mesh->m_nodes[m_mesh->m_edges[e].first],
                                                                 m_mesh->m_nodes[m_mesh->m_edges[e].second],
                                                                 false,
                                                                 m_mesh->m_projection,
                                                                 intersectionPoint,
                                                                 crossProductValue,
                                                                 ratioFirstSegment,
                                                                 ratioSecondSegment);
            if (areSegmentsCrossing)
            {
                // mask all faces crossed by the boundary line
                for (auto f = 0; f < m_mesh->m_edgesNumFaces[e]; ++f)
                {
                    const auto face = m_mesh->m_edgesFaces[e][f];
                    faceMask[face] = true;
                }
            }
        }
    }

    // Loop over the edges and mark as virtual nodes the inactive nodes belonging to a masked face
    for (auto f = 0; f < m_mesh->GetNumFaces(); ++f)
    {
        if (!faceMask[f])
        {
            continue;
        }
        for (int n = 0; n < m_mesh->m_numFacesNodes[f]; ++n)
        {
            const auto node = m_mesh->m_facesNodes[f][n];
            if (nodeMask[node] == NodeClasses::inactiveFlag)
            {
                nodeMask[node] = NodeClasses::virtualNodeFlag;
            }
        }
    }

    return nodeMask;
}