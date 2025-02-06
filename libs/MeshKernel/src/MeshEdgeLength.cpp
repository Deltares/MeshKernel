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

#include "MeshKernel/MeshEdgeLength.hpp"
#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Operations.hpp"

std::vector<double> meshkernel::MeshEdgeLength::Compute(const Mesh& mesh)
{
    std::vector<double> length(mesh.GetNumEdges(), constants::missing::doubleValue);
    Compute(mesh, length);

    return length;
}

double meshkernel::MeshEdgeLength::ComputeValue(const Mesh& mesh, const UInt edgeId)
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
                                                mesh.m_facesCircumcenters[firstFaceIndex],
                                                mesh.m_facesCircumcenters[secondFaceIndex],
                                                mesh.m_projection);

        if (val != constants::missing::doubleValue)
        {
            val = std::abs(val);
        }
    }

    return val;
}

void meshkernel::MeshEdgeLength::Compute(const Mesh& mesh, std::span<double> length)
{
    if (length.size() != mesh.GetNumEdges())
    {
        throw ConstraintError("array for length values is not the correct size");
    }

    const auto numEdges = mesh.GetNumEdges();

#pragma omp parallel for
    for (UInt e = 0; e < numEdges; e++)
    {
        length[e] = ComputeValue(mesh, e);
    }
}

double meshkernel::MeshEdgeLength::MinEdgeLength(const Mesh& mesh, const Polygons& polygon, const std::span<const double> edgeLengths)
{
    auto const numEdges = mesh.GetNumEdges();

    const auto isNodeInPolygon = mesh.IsLocationInPolygon(polygon, Location::Nodes);
    auto result = std::numeric_limits<double>::max();
#pragma omp parallel for reduction(min : result)
    for (UInt e = 0; e < numEdges; e++)
    {
        const auto& [firstNode, secondNode] = mesh.GetEdge(e);

        if (isNodeInPolygon[firstNode] || isNodeInPolygon[secondNode])
        {
            result = std::min(result, edgeLengths[e]);
        }
    }

    return result;
}
