//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2026.
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

#include "MeshKernel/Mesh2DFaceBounds.hpp"

#include <algorithm>

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/MeshFaceCenters.hpp"
#include "MeshKernel/Operations.hpp"

std::vector<meshkernel::Point> meshkernel::algo::Mesh2DFaceBounds::Compute(const Mesh& mesh)
{
    std::vector<Point> meshFaceBounds(constants::geometric::maximumNumberOfNodesPerFace * mesh.GetNumFaces());

#pragma omp parallel for
    for (int face = 0; face < static_cast<int>(mesh.GetNumFaces()); ++face)
    {
        const UInt pointCount = constants::geometric::maximumNumberOfNodesPerFace * static_cast<UInt>(face);

        std::span<Point> faceBounds(meshFaceBounds.data() + pointCount, meshFaceBounds.data() + pointCount + constants::geometric::maximumNumberOfNodesPerFace);
        ComputeBoundsForFace(mesh, static_cast<UInt>(face), faceBounds);
    }

    return meshFaceBounds;
}

void meshkernel::algo::Mesh2DFaceBounds::ComputeBoundsForFace(const Mesh& mesh, UInt faceId, std::span<Point> faceBounds)
{
    if (mesh.IsValidFace(faceId))
    {
        UInt numNodes = static_cast<UInt>(mesh.m_numFacesNodes[faceId]);

        for (UInt i = 0; i < numNodes; ++i)
        {
            faceBounds[i] = mesh.Node(mesh.m_facesNodes[faceId][i]);
        }

        for (UInt i = numNodes; i < constants::geometric::maximumNumberOfNodesPerFace; ++i)
        {
            faceBounds[i] = Point{constants::missing::doubleValue, constants::missing::doubleValue};
        }
    }
    else
    {
        std::ranges::fill(faceBounds, Point{constants::missing::doubleValue, constants::missing::doubleValue});
    }
}
