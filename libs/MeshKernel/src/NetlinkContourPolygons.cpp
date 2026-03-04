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

#include "MeshKernel/NetlinkContourPolygons.hpp"

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/MeshFaceCenters.hpp"
#include "MeshKernel/Operations.hpp"

std::vector<meshkernel::Point> meshkernel::algo::NetlinkContourPolygons::Compute(const Mesh& mesh)
{
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
    // along the normals pointing towards circumcentre1 and circumcentre12

    Point edgeNormal;
    bool normalReflected = false;
    // The normal computed here will be pointing to the same side of the edge as the circimcentreLeft point
    NormalVectorInside(edgeStart, edgeEnd, circumcentre1, edgeNormal, normalReflected, projection);

    // Get distance of circumcentre from the edge
    auto [distanceToC1, np, rat] = DistanceFromLine(circumcentre1, edgeStart, edgeEnd, projection);

    if (normalReflected)
    {
        std::swap(edgeStart, edgeEnd);
    }

    polygon[0] = edgeEnd - ComputePointIncrement (edgeNormal, distanceToC1, edgeEnd, projection);
    polygon[1] = edgeStart - ComputePointIncrement (edgeNormal, distanceToC1, edgeStart, projection);

    // the second circumcentre is valid indicates the edge has two connecting elements
    if (circumcentre2.IsValid())
    {
        // Get distance of circumcentre from the edge
        auto [distanceToC2, np2, rat2] = DistanceFromLine(circumcentre2, edgeStart, edgeEnd, projection);

        polygon[2] = edgeStart + ComputePointIncrement (edgeNormal, distanceToC2, edgeStart, projection);
        polygon[3] = edgeEnd + ComputePointIncrement (edgeNormal, distanceToC2, edgeEnd, projection);
    }
    else
    {
        // the second circumcentre is invalid indicates the edge has only one connecting element
        polygon[2] = edgeStart;
        polygon[3] = edgeEnd;
    }
}
