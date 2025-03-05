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

#include "MeshKernel/MeshFaceCenters.hpp"
#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Operations.hpp"

std::vector<meshkernel::Point> meshkernel::MeshFaceCenters::ComputeCircumcenters(const Mesh& mesh)
{
    std::vector<Point> faceCenters(mesh.GetNumFaces());
    ComputeCircumcenters(mesh, faceCenters);

    return faceCenters;
}

void meshkernel::MeshFaceCenters::ComputeCircumcenters(const Mesh& mesh, std::span<Point> faceCenters)
{
    if (faceCenters.size() != mesh.GetNumFaces())
    {
        throw ConstraintError("array for faceCenters values is not the correct size");
    }

    const auto numFaces = static_cast<int>(mesh.GetNumFaces());

    std::vector<UInt> numEdgeFacesCache;
    numEdgeFacesCache.reserve(constants::geometric::maximumNumberOfEdgesPerFace);
    std::vector<Point> polygonNodesCache;
    polygonNodesCache.reserve(constants::geometric::maximumNumberOfEdgesPerFace);

#pragma omp parallel for private(numEdgeFacesCache, polygonNodesCache)
    for (int f = 0; f < numFaces; f++)
    {

        UInt numberOfInteriorEdges = 0;
        const auto numberOfFaceNodes = mesh.GetNumFaceEdges(f);

        for (UInt n = 0; n < numberOfFaceNodes; ++n)
        {
            if (!mesh.IsEdgeOnBoundary(mesh.m_facesEdges[f][n]))
            {
                numberOfInteriorEdges += 1;
            }
        }

        if (numberOfInteriorEdges == 0)
        {
            faceCenters[f] = mesh.m_facesMassCenters[f];
        }
        else
        {
            // need to account for spherical coordinates. Build a polygon around a face
            mesh.ComputeFaceClosedPolygon(f, polygonNodesCache);
            numEdgeFacesCache.clear();

            for (UInt n = 0; n < numberOfFaceNodes; ++n)
            {
                numEdgeFacesCache.emplace_back(mesh.m_edgesNumFaces[mesh.m_facesEdges[f][n]]);
            }

            faceCenters[f] = ComputeFaceCircumenter(polygonNodesCache, numEdgeFacesCache, mesh.m_projection);
        }
    }
}
