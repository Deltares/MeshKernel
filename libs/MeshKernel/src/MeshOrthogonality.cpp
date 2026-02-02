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

#include "MeshKernel/MeshOrthogonality.hpp"
#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/MeshFaceCenters.hpp"
#include "MeshKernel/Operations.hpp"

std::vector<double> meshkernel::MeshOrthogonality::Compute(const Mesh2D& mesh)
{
    std::vector<double> orthogonality(mesh.GetNumEdges(), constants::missing::doubleValue);
    Compute(mesh, orthogonality);

    return orthogonality;
}

double meshkernel::MeshOrthogonality::ComputeValue(const Mesh2D& mesh, const std::vector<Point>& faceCircumcentres, const UInt edgeId)
{
    const auto [firstNode, secondNode] = mesh.GetEdge(edgeId);

    if (firstNode == constants::missing::uintValue ||
        secondNode == constants::missing::uintValue)
    {
        return constants::missing::doubleValue;
    }

    const auto firstFaceIndex = mesh.m_edgesFaces[edgeId][0];
    const auto secondFaceIndex = mesh.m_edgesFaces[edgeId][1];

    if (firstFaceIndex == constants::missing::uintValue ||
        secondFaceIndex == constants::missing::uintValue)
    {
        return constants::missing::doubleValue;
    }

    double val = constants::missing::doubleValue;

    if (!mesh.IsEdgeOnBoundary(edgeId))
    {
        val = NormalizedInnerProductTwoSegments(mesh.Node(firstNode),
                                                mesh.Node(secondNode),
                                                faceCircumcentres[firstFaceIndex],
                                                faceCircumcentres[secondFaceIndex],
                                                mesh.m_projection);

        if (val != constants::missing::doubleValue)
        {
            val = std::abs(val);
        }
    }

    return val;
}

void meshkernel::MeshOrthogonality::Compute(const Mesh2D& mesh, std::span<double> orthogonality)
{
    if (orthogonality.size() != mesh.GetNumEdges())
    {
        throw ConstraintError("array for orthogonality values is not the correct size");
    }

    std::vector<Point> faceCircumcentres = algo::ComputeFaceCircumcenters(mesh);

    const auto numEdges = mesh.GetNumEdges();

    for (UInt e = 0; e < numEdges; e++)
    {
        orthogonality[e] = ComputeValue(mesh, faceCircumcentres, e);
    }
}
