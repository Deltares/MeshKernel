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

#if USE_LIBFMT
#include <fmt/format.h>
namespace fmt_ns = fmt;
#else
#include <format>
namespace fmt_ns = std;
#endif

#include "MeshKernel/NetlinkContourPolygons.hpp"

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/MeshFaceCenters.hpp"
#include "MeshKernel/Operations.hpp"

void meshkernel::algo::NetlinkContourPolygons::CheckEdges(const Mesh2D& mesh)
{

    // Should call mesh.GetHangingEdges, however this function does not catch the case where
    // an edge connects two nodes and these two nodes are each part of different, unconnected
    // elements and the edge itself has no connected elements.
    //
    //  +------+     +-----+
    //  |      |     |     |
    //  |      |     |     |
    //  +------+=====+-----+
    //  |      |     |     |
    //  |      |     |     |
    //  +------+     +-----+
    //
    //  The edge marked with ='s
    //
    std::vector<UInt> hangingEdges;

    const auto edges = mesh.Edges();

    for (UInt e = 0; e < edges.size(); e++)
    {
        const auto firstNode = edges[e].first;
        const auto secondNode = edges[e].second;

        if (firstNode != constants::missing::uintValue && secondNode != constants::missing::uintValue)
        {
            // if one of the nodes has no other attached edges or attached faces, the current edge is an hanging edge
            if (mesh.m_nodesNumEdges[firstNode] <= 1 || mesh.m_nodesNumEdges[secondNode] <= 1 || mesh.m_edgesNumFaces[e] == 0)
            {
                hangingEdges.emplace_back(e);
            }
        }
    }

    if (hangingEdges.size() > 0)
    {
        std::string buffer;

        bool firstEdge = true;

        for (UInt edgeId : hangingEdges)
        {

            if (!firstEdge)
            {
                buffer += ", ";
            }

            const Edge& edge = mesh.GetEdge(edgeId);
            // It is not necessary to check if the node-id is valid.
            // Edges with invalid node-ids are not added to the hanging-edges list.
            const Point node = mesh.Node(edge.first);

            fmt_ns::format_to(std::back_inserter(buffer), "edge={}, node={{x={:5f}, y={:5f}}}", edgeId, node.x, node.y);
            firstEdge = false;
        }

        throw MeshKernelError("Deriving netlink contours failed due to {} hanging edge(s): {}", hangingEdges.size(), buffer);
    }
}

std::vector<meshkernel::Point> meshkernel::algo::NetlinkContourPolygons::Compute(const Mesh2D& mesh)
{

    CheckEdges(mesh);

    std::vector<Point> netlinkPolygons(4 * mesh.GetNumEdges(), Point{constants::missing::doubleValue, constants::missing::doubleValue});
    std::vector<Point> circumcentres(algo::ComputeFaceCircumcenters(mesh));

#pragma omp parallel for
    for (int edge = 0; edge < static_cast<int>(mesh.GetNumEdges()); ++edge)
    {
        if (mesh.IsValidEdge(edge))
        {
            const UInt pointCount = 4 * static_cast<UInt>(edge);

            const Edge currentEdge = mesh.GetEdge(edge);
            const auto edgeFaces = mesh.m_edgesFaces[edge];

            std::span<Point> edgePolygon(netlinkPolygons.data() + pointCount, netlinkPolygons.data() + pointCount + 4);
            Point startNode = mesh.Node(currentEdge.first);
            Point endNode = mesh.Node(currentEdge.second);
            Point circumcentreLeft = circumcentres[edgeFaces[0]];
            Point circumcentreRight;

            if (mesh.m_edgesNumFaces[edge] == 2)
            {
                circumcentreRight = circumcentres[edgeFaces[1]];
            }
            else
            {
                circumcentreRight.SetInvalid();
            }

            ComputePolygonForEdge(startNode, endNode, circumcentreLeft, circumcentreRight, mesh.m_projection, edgePolygon);
        }
    }

    return netlinkPolygons;
}

void meshkernel::algo::NetlinkContourPolygons::ComputePolygonForEdge(Point edgeStart, Point edgeEnd,
                                                                     const Point& circumcentre1, const Point& circumcentre2,
                                                                     const Projection projection,
                                                                     std::span<Point> polygon)
{

    // Project circumcenters onto the line passing through the edge
    // to find their lateral offsets relative to the edge line.

    // The polygon vertices are formed by shifting edgeStart and edgeEnd
    // along the normals pointing towards circumcentre1 and circumcentre2

    Point edgeNormal;
    bool normalReflected = false;
    // The normal computed here will be pointing to the same side of the edge as the circimcentre1 point
    NormalVectorInside(edgeStart, edgeEnd, circumcentre1, edgeNormal, normalReflected, projection);

    // Get distance of circumcentre from the edge
    auto [distanceToC1, np, rat] = DistanceFromLine(circumcentre1, edgeStart, edgeEnd, projection);

    if (normalReflected)
    {
        std::swap(edgeStart, edgeEnd);
    }

    polygon[0] = edgeEnd - ComputePointIncrement(edgeNormal, distanceToC1, edgeEnd, projection);
    polygon[1] = edgeStart - ComputePointIncrement(edgeNormal, distanceToC1, edgeStart, projection);

    // the second circumcentre is valid indicates the edge has two connecting elements
    if (circumcentre2.IsValid())
    {
        // Get distance of circumcentre from the edge
        auto [distanceToC2, np2, rat2] = DistanceFromLine(circumcentre2, edgeStart, edgeEnd, projection);

        polygon[2] = edgeStart + ComputePointIncrement(edgeNormal, distanceToC2, edgeStart, projection);
        polygon[3] = edgeEnd + ComputePointIncrement(edgeNormal, distanceToC2, edgeEnd, projection);
    }
    else
    {
        // the second circumcentre is invalid indicates the edge has only one connecting element
        polygon[2] = edgeStart;
        polygon[3] = edgeEnd;
    }
}
