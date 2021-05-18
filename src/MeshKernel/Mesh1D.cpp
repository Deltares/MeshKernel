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

#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh1D.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Polygons.hpp>
#include <vector>

meshkernel::Mesh1D::Mesh1D(const std::vector<Edge>& edges,
                           const std::vector<Point>& nodes,
                           Projection projection) : Mesh(edges, nodes, projection){};

meshkernel::Mesh1D::Mesh1D(std::vector<std::vector<Point>> const& polylines, double offset, double mergingLength, Projection projection)
{
    if (offset <= mergingLength)
    {
        throw std::invalid_argument("Mesh1D::Mesh1D: Offset cannot be smaller than merging length");
    }

    std::vector<Edge> edges;
    std::vector<Point> nodes;
    size_t numNodes = 0;
    for (auto const& branch : polylines)
    {
        auto const branchNodes = RefinePolyLine(branch, offset, projection);
        std::copy(branchNodes.begin(), branchNodes.end(), back_inserter(nodes));
        for (auto i = numNodes; i < nodes.size() - 1; ++i)
        {
            edges.emplace_back(i, i + 1);
        }
        // branches are separated. If the end of a polyline coincides with the start of another, the two nodes will be merged.
        numNodes = numNodes + nodes.size();
    }

    // Sets the edges, nodes and projections
    m_edges = edges;
    m_nodes = nodes;
    m_projection = projection;

    // Perform node administration to fill the internal arrays
    AdministrateNodesEdges();

    // If there are computational nodes at a distance smaller than  the threshold, these are eliminated
    const Polygons polygon{};
    MergeNodesInPolygon(polygon, mergingLength);
}
