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

#include "MeshKernel/TriangulationWrapper.hpp"
#include "MeshKernel/Mesh.hpp"
#include "MeshKernel/Operations.hpp"

using namespace meshkernel;

void TriangulationWrapper::BuildTriangulation()
{

    if (m_numFaces < 0)
    {
        m_numFaces = 0;
    }

    // Create nodes
    m_nodes.resize(m_numNodes);
    for (auto i = 0; i < m_numNodes; ++i)
    {
        m_nodes[i] = {m_xCoordFlat[i], m_yCoordFlat[i]};
    }

    // Create m_faceNodes
    ResizeAndFill2DVector(m_faceNodes, m_numFaces, 3, true, constants::missing::uintValue);
    ResizeAndFill2DVector(m_faceEdges, m_numFaces, 3, true, constants::missing::uintValue);
    UInt faceCounter = 0;
    for (int f = 0; f < m_numFaces; ++f)
    {
        for (int e = 0; e < 3; ++e)
        {
            m_faceNodes[f][e] = static_cast<UInt>(m_faceNodesFlat[faceCounter] - 1);
            m_faceEdges[f][e] = static_cast<UInt>(m_faceEdgesFlat[faceCounter] - 1);
            faceCounter++;
        }
    }

    // Create edges
    if (m_numEdges == 0)
    {
        return;
    }

    ResizeAndFill2DVector(m_edgeNodes, m_numEdges, 2, true, constants::missing::uintValue);
    UInt edgeCounter = 0;
    for (int e = 0; e < m_numEdges; ++e)
    {
        for (int n = 0; n < 2; ++n)
        {
            m_edgeNodes[e][n] = static_cast<UInt>(m_edgeNodesFlat[edgeCounter] - 1);
            edgeCounter++;
        }
    }

    ResizeAndFill2DVector(m_edgesFaces, m_numEdges, 2, true, constants::missing::uintValue);
    edgeCounter = 0;
    for (int f = 0; f < m_numFaces; ++f)
    {

        for (UInt n = 0; n < Mesh::m_numNodesInTriangle; ++n)
        {
            auto const edge = static_cast<UInt>(m_faceEdgesFlat[edgeCounter] - 1);
            edgeCounter++;
            // For each edge, the shared face index
            if (m_edgesFaces[edge][0] == constants::missing::uintValue)
            {
                m_edgesFaces[edge][0] = f;
            }
            else
            {
                m_edgesFaces[edge][1] = f;
            }
        }
    }
}
