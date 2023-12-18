#include "MeshKernel/Mesh2DGenerateGlobalGrid.hpp"
#include <cmath>

using namespace meshkernel::constants;

double meshkernel::Mesh2DGenerateGlobalGrid::getDeltaLatitude(const double currentLatitude, const double longitudeDiscretization)
{

    double deltaLatitude = longitudeDiscretization * std::cos(conversion::degToRad * currentLatitude);
    constexpr UInt numIterations = 5;
    constexpr double tolerance = 1.0e-14;

    for (int i = 0; i < numIterations; ++i)
    {
        double phi = conversion::degToRad * (currentLatitude + 0.5 * deltaLatitude);
        double c = std::cos(phi);
        double s = std::sqrt(1.0 - c * c);
        double f = deltaLatitude - longitudeDiscretization * c;
        double df = 1.0 + 0.5 * conversion::degToRad * longitudeDiscretization * s;
        double yd = f / df;
        deltaLatitude = deltaLatitude - yd;

        if (yd < tolerance)
        {
            break;
        }
    }

    return deltaLatitude;
}

meshkernel::UInt meshkernel::Mesh2DGenerateGlobalGrid::getNodeIndexFromPosition(const Mesh& mesh, const Point& x)
{
    constexpr double tolerance = 1.0e-6;
    const auto nodeIndex = missing::uintValue;

    for (int i = static_cast<int>(mesh.m_nodes.size() - 1); i >= 0; --i)
    {
        if (IsEqual(x, mesh.m_nodes[i], tolerance))
        {
            return static_cast<UInt>(i);
        }
    }
    return nodeIndex;
}

void meshkernel::Mesh2DGenerateGlobalGrid::addFace(Mesh& mesh,
                                                   const std::array<Point, 8>& points,
                                                   const double latitudeDirection,
                                                   const UInt numNodes)
{
    std::array<UInt, 8> nodeIndices{};

    for (UInt i = 0; i < numNodes; ++i)
    {
        Point p = {points[i].x, latitudeDirection * points[i].y};
        nodeIndices[i] = getNodeIndexFromPosition(mesh, p);

        if (nodeIndices[i] == missing::uintValue)
        {
            nodeIndices[i] = mesh.InsertNode(p);
        }
    }

    for (UInt currentNode = 0; currentNode < numNodes; ++currentNode)
    {
        UInt nextNode = currentNode + 1;

        if (currentNode == numNodes - 1)
        {
            nextNode = 0;
        }
        mesh.ConnectNodes(nodeIndices[currentNode], nodeIndices[nextNode]);
    }
}

std::unique_ptr<meshkernel::Mesh2D> meshkernel::Mesh2DGenerateGlobalGrid::Compute(const UInt numX,
                                                                                  const UInt numY,
                                                                                  const Polygons& polygon)
{
    std::array<Point, 8> points;

    double deltaLongitude = 360.0 / static_cast<double>(numX);
    double currentLatitude = 0.0;

    bool pentagonFace = false;
    bool generationCompleted = false;
    UInt numberOfPoints;

    auto mesh2d = std::make_unique<Mesh2D>(Projection::sphericalAccurate);
    constexpr double minDiscretizationLength = 30000.0;

    for (UInt i = 0; i < numY; ++i)
    {
        double deltaLatitude = getDeltaLatitude(currentLatitude, deltaLongitude);

        if (currentLatitude + 1.5 * deltaLatitude > 90.0)
        {
            deltaLatitude = 90.0 - currentLatitude;
            generationCompleted = true;
            pentagonFace = false;
        }
        else
        {
            if (deltaLatitude * conversion::degToRad * geometric::earth_radius < minDiscretizationLength &&
                !pentagonFace)
            {
                deltaLongitude = 2.0 * deltaLongitude;
                pentagonFace = true;
                deltaLatitude = getDeltaLatitude(currentLatitude, deltaLongitude);
            }
            else
            {
                pentagonFace = false;
            }

            if (currentLatitude + 1.5 * deltaLatitude > 90.0)
            {
                deltaLatitude = 0.51 * (90.0 - currentLatitude);
            }
        }

        for (UInt j = 1; j <= numX; ++j)
        {
            double currentLongitude = static_cast<double>(j - 1) * deltaLongitude + -180.0;

            points[0] = {currentLongitude, currentLatitude};

            if (!pentagonFace)
            {
                points[1] = {currentLongitude + deltaLongitude, currentLatitude};
                points[2] = {currentLongitude + deltaLongitude, currentLatitude + deltaLatitude};
                points[3] = {currentLongitude, currentLatitude + deltaLatitude};
                numberOfPoints = 4;
            }
            else
            {
                points[1] = {currentLongitude + 0.5 * deltaLongitude, currentLatitude};
                points[2] = {currentLongitude + deltaLongitude, currentLatitude};
                points[3] = {currentLongitude + deltaLongitude, currentLatitude + deltaLatitude};
                points[4] = {currentLongitude, currentLatitude + deltaLatitude};
                numberOfPoints = 5;
            }

            if (points[2].x <= 180.0)
            {
                addFace(*mesh2d, points, 1.0, numberOfPoints);
                addFace(*mesh2d, points, -1.0, numberOfPoints);
            }
        }

        if (generationCompleted)
        {
            break;
        }

        currentLatitude += deltaLatitude;
    }

    mesh2d->MergeNodesInPolygon(polygon, 1e-3);

    constexpr double tolerance = 1.0e-6;
    for (UInt e = 0; e < mesh2d->GetNumEdges(); e++)
    {
        const auto& [firstNode, secondNode] = mesh2d->m_edges[e];

        if ((mesh2d->m_nodesNumEdges[firstNode] == geometric::numNodesInPentagon || mesh2d->m_nodesNumEdges[firstNode] == geometric::numNodesInhaxagon) &&
            (mesh2d->m_nodesNumEdges[secondNode] == geometric::numNodesInPentagon || mesh2d->m_nodesNumEdges[secondNode] == geometric::numNodesInhaxagon))
        {
            if (IsEqual(mesh2d->m_nodes[firstNode].y, mesh2d->m_nodes[secondNode].y, tolerance))
            {
                mesh2d->DeleteEdge(e);
            }
        }
    }

    mesh2d->Administrate();

    return mesh2d;
}
