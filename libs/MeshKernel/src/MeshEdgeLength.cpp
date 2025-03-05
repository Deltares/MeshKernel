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

std::vector<double> meshkernel::algo::ComputeMeshEdgeLength(const Mesh& mesh)
{
    std::vector<double> length(mesh.GetNumEdges(), constants::missing::doubleValue);
    ComputeMeshEdgeLength(mesh, length);

    return length;
}

double meshkernel::algo::ComputeEdgeLength(const Mesh& mesh, const UInt edgeId)
{
    const auto [firstNode, secondNode] = mesh.GetEdge(edgeId);

    if (firstNode == constants::missing::uintValue ||
        secondNode == constants::missing::uintValue)
    {
        return constants::missing::doubleValue;
    }

    double val = ComputeDistance(mesh.Node(firstNode), mesh.Node(secondNode), mesh.m_projection);

    return val;
}

void meshkernel::algo::ComputeMeshEdgeLength(const Mesh& mesh, std::span<double> length)
{
    if (length.size() != mesh.GetNumEdges())
    {
        throw ConstraintError("array for length values is not the correct size");
    }

    const auto numEdges = static_cast<int>(mesh.GetNumEdges());

#pragma omp parallel for
    for (int e = 0; e < numEdges; e++)
    {
        length[e] = ComputeEdgeLength(mesh, static_cast<UInt>(e));
    }
}

double meshkernel::algo::MinEdgeLength(const Mesh& mesh, const Polygons& polygon, const std::span<const double> edgeLengths)
{
    const int numEdges = static_cast<int>(mesh.GetNumEdges());

    const auto isNodeInPolygon = mesh.IsLocationInPolygon(polygon, Location::Nodes);
    auto result = std::numeric_limits<double>::max();

    for (int e = 0; e < numEdges; e++)
    {
        const auto& [firstNode, secondNode] = mesh.GetEdge(static_cast<UInt>(e));

        if (isNodeInPolygon[firstNode] || isNodeInPolygon[secondNode])
        {
            result = std::min(result, edgeLengths[e]);
        }
    }

    return result;
}

double meshkernel::algo::MaxLengthSurroundingEdges(const Mesh& mesh,
                                                   const UInt nodeId,
                                                   const std::span<const double> edgeLengths)
{
    double maxEdgeLength = std::numeric_limits<double>::lowest();

    for (UInt e = 0; e < mesh.m_nodesNumEdges[nodeId]; ++e)
    {
        const auto edge = mesh.m_nodesEdges[nodeId][e];
        maxEdgeLength = std::max(maxEdgeLength, edgeLengths[edge]);
    }

    return maxEdgeLength;
}
