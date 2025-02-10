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

#include "MeshKernel/MeshEdgeCenters.hpp"
#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Operations.hpp"

std::vector<meshkernel::Point> meshkernel::MeshEdgeCenters::Compute(const Mesh& mesh)
{
    std::vector<Point> edgeCentres(mesh.GetNumEdges());
    Compute(mesh, edgeCentres);

    return edgeCentres;
}

void meshkernel::MeshEdgeCenters::Compute(const Mesh& mesh, std::span<Point> edgeCentres)
{
    if (edgeCentres.size() != mesh.GetNumEdges())
    {
        throw ConstraintError("array for edgeCentres values is not the correct size");
    }

    const auto numEdges = static_cast<int>(mesh.GetNumEdges());

#pragma omp parallel for
    for (int e = 0; e < numEdges; e++)
    {
        edgeCentres[e] = ComputeValue(mesh, static_cast<UInt>(e));
    }
}

meshkernel::Point meshkernel::MeshEdgeCenters::ComputeValue(const Mesh& mesh, const UInt edgeId)
{
    if (edgeId == constants::missing::uintValue)
    {
        throw ConstraintError("Invalid edgeId");
    }

    const auto [firstNode, secondNode] = mesh.GetEdge(edgeId);

    Point edgeCentre;

    if (firstNode == constants::missing::uintValue || secondNode == constants::missing::uintValue)
    {
        // Initialise with invalid data.
        edgeCentre = Point(constants::missing::doubleValue, constants::missing::doubleValue);
    }
    else
    {
        edgeCentre = (mesh.Node(firstNode) + mesh.Node(secondNode)) * 0.5;
    }

    return edgeCentre;
}
