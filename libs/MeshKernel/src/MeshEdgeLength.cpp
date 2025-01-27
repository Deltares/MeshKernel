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

#include <algorithm>

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/MeshEdgeLength.hpp"
#include "MeshKernel/Operations.hpp"

std::vector<double> meshkernel::MeshEdgeLength::Compute(const Mesh& mesh)
{
    std::vector<double> edgeLengths(mesh.GetNumEdges(), constants::missing::doubleValue);
    Compute(mesh, edgeLengths);

    return edgeLengths;
}

double meshkernel::MeshEdgeLength::ComputeValue(const Mesh& mesh, const UInt edgeId)
{
    const auto [firstNodeIndex, secondNodeIndex] = mesh.GetEdge(edgeId);

    if (firstNodeIndex == constants::missing::uintValue ||
        secondNodeIndex == constants::missing::uintValue)
    {
        return constants::missing::doubleValue;
    }

    const Point& firstNode = mesh.Node(firstNodeIndex);
    const Point& secondNode = mesh.Node(secondNodeIndex);

    if (!firstNode.IsValid() || !secondNode.IsValid())
    {
        return constants::missing::doubleValue;
    }

    double val = ComputeDistance(firstNode, secondNode, mesh.m_projection);

    return val;
}

void meshkernel::MeshEdgeLength::Compute(const Mesh& mesh, std::span<double> edgeLength)
{
    if (edgeLength.size() != mesh.GetNumEdges())
    {
        throw ConstraintError("array for edgeLength values is not the correct size");
    }

    const auto numEdges = mesh.GetNumEdges();

    for (UInt e = 0; e < numEdges; e++)
    {
        edgeLength[e] = ComputeValue(mesh, e);
    }
}

double meshkernel::MeshEdgeLength::MaxLengthSurroundEdges(const Mesh& mesh, const UInt nodeId)
{
    if (nodeId == constants::missing::uintValue)
    {
        throw ConstraintError("Invalid node id");
    }

    const std::vector<UInt>& edgeIds = mesh.m_nodesEdges[nodeId];
    double maxEdgeLength = constants::missing::doubleValue;

    for (UInt i = 0; i < edgeIds.size(); ++i)
    {
        double length = ComputeValue(mesh, edgeIds[i]);
        maxEdgeLength = maxEdgeLength == constants::missing::doubleValue ? length : std::max(maxEdgeLength, length);
    }

    return maxEdgeLength;
}

double meshkernel::MeshEdgeLength::MinEdgeLength(const Mesh& mesh, const std::span<const double> edgeLengths, const Polygons& polygon)
{

    auto const numEdges = mesh.GetNumEdges();
    auto result = std::numeric_limits<double>::max();

    const auto isNodeInPolygon = mesh.IsLocationInPolygon(polygon, Location::Nodes);

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
