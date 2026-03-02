#include "MeshKernel/NetlinkContourPolygons.hpp"

#include "MeshKernel/MeshFaceCenters.hpp"
#include "MeshKernel/Operations.hpp"

std::vector<meshkernel::Point> meshkernel::algo::NetlinkContourPolygons::Compute(const Mesh& mesh)
{
    std::vector<Point> netlinkPolygons(mesh.GetNumEdges() * 4);
    std::vector<Point> circumcentres(algo::ComputeFaceCircumcenters(mesh));

#pragma omp parallel
    for (int edge = 0; edge < static_cast<int>(mesh.GetNumEdges()); ++edge)
    {
        if (mesh.IsValidEdge(edge))
        {
            UInt pointCount = 4 * edge;

            const Edge currentEdge = mesh.GetEdge(edge);
            const auto edgeFaces = mesh.m_edgesFaces[edge];

            // Start and end are swapped in order to get the correct orientation of the polygon
            Point startNode = mesh.Node(currentEdge.second);
            Point endNode = mesh.Node(currentEdge.first);
            std::span<Point> edgePolygon(netlinkPolygons.data() + pointCount, netlinkPolygons.data() + pointCount + 4);
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

void meshkernel::algo::NetlinkContourPolygons::ComputePolygonForEdge(const Point& edgeStart, const Point& edgeEnd,
                                                                     const Point& circumcentreLeft, const Point& circumcentreRight,
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
    NormalVectorInside(edgeStart, edgeEnd, circumcentreLeft, edgeNormal, normalReflected, projection);

    // Distances from edge line to circumcenters
    auto getDist = [&](const Point& pnt)
    {
        return (pnt.x - edgeStart.x) * edgeNormal.x + (pnt.y - edgeStart.y) * edgeNormal.y;
    };

    double distanceToC1 = getDist(circumcentreLeft);

    if (projection == Projection::spherical)
    {
        distanceToC1 *= constants::conversion::mtsToDeg;
    }

    polygon[0] = {edgeStart.x + distanceToC1 * edgeNormal.x, edgeStart.y + distanceToC1 * edgeNormal.y};
    polygon[1] = {edgeEnd.x + distanceToC1 * edgeNormal.x, edgeEnd.y + distanceToC1 * edgeNormal.y};

    // the second circumcentre is valid indicates the edge has two connecting elements
    if (circumcentreRight.IsValid())
    {
        double distanceToC2 = getDist(circumcentreRight);

        if (projection == Projection::spherical)
        {
            distanceToC2 *= constants::conversion::mtsToDeg;
        }

        polygon[2] = {edgeEnd.x + distanceToC2 * edgeNormal.x, edgeEnd.y + distanceToC2 * edgeNormal.y};
        polygon[3] = {edgeStart.x + distanceToC2 * edgeNormal.x, edgeStart.y + distanceToC2 * edgeNormal.y};
    }
    else
    {
        // the second circumcentre is invalid indicates the edge has only one connecting element
        polygon[2] = edgeEnd;
        polygon[3] = edgeStart;
    }
}
