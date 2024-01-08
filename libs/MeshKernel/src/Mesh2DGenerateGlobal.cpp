#include "MeshKernel/Mesh2DGenerateGlobal.hpp"
#include "MeshKernel/Operations.hpp"
#include "MeshKernel/Polygons.hpp"
#include <cmath>

using namespace meshkernel;

double Mesh2DGenerateGlobal::DeltaLatitude(const double currentLatitude, const double longitudeDiscretization)
{

    double deltaLatitude = longitudeDiscretization * std::cos(constants::conversion::degToRad * currentLatitude);

    for (UInt i = 0; i < numIterations; ++i)
    {
        const double phi = constants::conversion::degToRad * (currentLatitude + 0.5 * deltaLatitude);
        const double cosPhi = std::cos(phi);
        const double sinPhi = std::sqrt(1.0 - cosPhi * cosPhi);
        const double f = deltaLatitude - longitudeDiscretization * cosPhi;
        const double df = 1.0 + 0.5 * constants::conversion::degToRad * longitudeDiscretization * sinPhi;
        const double latitudeAdjustment = f / df;
        deltaLatitude = deltaLatitude - latitudeAdjustment;

        if (latitudeAdjustment < toleranceDeltaLatitude)
        {
            break;
        }
    }

    return deltaLatitude;
}

UInt Mesh2DGenerateGlobal::NodeIndexFromPosition(const Mesh& mesh, const Point& x)
{
    constexpr double tolerance = 1.0e-6;

    for (auto i = static_cast<int>(mesh.m_nodes.size() - 1); i >= 0; --i)
    {
        if (IsEqual(x, mesh.m_nodes[i], tolerance))
        {
            return static_cast<UInt>(i);
        }
    }
    return constants::missing::uintValue;
}

void Mesh2DGenerateGlobal::AddFace(Mesh& mesh,
                                   const std::array<Point, 5>& points,
                                   const GridExpansionDirection growingDirection,
                                   const UInt numNodes)
{
    std::array<UInt, 5> nodeIndices{};

    for (UInt n = 0; n < numNodes; ++n)
    {
        const auto expansionMultiplier = static_cast<double>(growingDirection);
        Point p = {points[n].x, expansionMultiplier * points[n].y};
        nodeIndices[n] = NodeIndexFromPosition(mesh, p);

        if (nodeIndices[n] == constants::missing::uintValue)
        {
            nodeIndices[n] = mesh.InsertNode(p);
        }
    }

    for (UInt n = 0; n < numNodes; ++n)
    {
        UInt nextNode = n + 1;

        if (n == numNodes - 1)
        {
            nextNode = 0;
        }
        const auto& firstNodeIndex = nodeIndices[n];
        const auto& secondNodeIndex = nodeIndices[nextNode];

        if (mesh.FindEdgeWithLinearSearch(firstNodeIndex, secondNodeIndex) == constants::missing::uintValue)
        {
            mesh.ConnectNodes(firstNodeIndex, secondNodeIndex);
        }
    }
}

std::unique_ptr<Mesh2D> Mesh2DGenerateGlobal::Compute(const UInt numLongitudeNodes,
                                                      const UInt numLatitudeNodes,
                                                      const Projection projection)
{
    if (numLongitudeNodes == 0)
    {
        throw MeshKernelError("The number of longitude nodes cannot be 0");
    }
    if (numLatitudeNodes == 0)
    {
        throw MeshKernelError("The number of latitude nodes cannot be 0");
    }
    if (projection != Projection::spherical && projection != Projection::sphericalAccurate)
    {
        throw MeshKernelError("Unsupported projection. The projection is not spherical nor sphericalAccurate");
    }

    std::array<Point, 5> points;
    double deltaLongitude = 360.0 / static_cast<double>(numLongitudeNodes);
    double currentLatitude = 0.0;
    bool pentagonFace = false;
    bool generationCompleted = false;

    auto mesh2d = std::make_unique<Mesh2D>(projection);
    constexpr double minDiscretizationLength = 30000.0;

    for (UInt i = 0; i < numLatitudeNodes; ++i)
    {
        double deltaLatitude = DeltaLatitude(currentLatitude, deltaLongitude);

        if (currentLatitude + 1.5 * deltaLatitude > 90.0)
        {
            deltaLatitude = 90.0 - currentLatitude;
            generationCompleted = true;
            pentagonFace = false;
        }
        else
        {
            if (deltaLatitude * constants::conversion::degToRad * constants::geometric::earth_radius < minDiscretizationLength &&
                !pentagonFace)
            {
                deltaLongitude = 2.0 * deltaLongitude;
                pentagonFace = true;
                deltaLatitude = DeltaLatitude(currentLatitude, deltaLongitude);
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

        for (UInt j = 0; j < numLongitudeNodes; ++j)
        {
            double currentLongitude = static_cast<double>(j) * deltaLongitude - 180.0;

            points[0] = {currentLongitude, currentLatitude};
            UInt numberOfPoints;
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
                pentagonFace = false;
                AddFace(*mesh2d, points, GridExpansionDirection::Northwards, numberOfPoints);
                AddFace(*mesh2d, points, GridExpansionDirection::Southwards, numberOfPoints);
            }
        }

        if (generationCompleted)
        {
            break;
        }

        currentLatitude += deltaLatitude;
    }

    constexpr double mergingDistance = 1e-3;
    const std::vector<Point> polygon;
    Polygons polygons(polygon, projection);
    mesh2d->MergeNodesInPolygon(polygons, mergingDistance);

    constexpr double tolerance = 1.0e-6;
    for (UInt e = 0; e < mesh2d->GetNumEdges(); e++)
    {
        const auto& [firstNode, secondNode] = mesh2d->m_edges[e];
        const auto numEdgesFirstNode = mesh2d->m_nodesNumEdges[firstNode];
        const auto numEdgesSecondNode = mesh2d->m_nodesNumEdges[secondNode];
        if ((numEdgesFirstNode == constants::geometric::numNodesInPentagon ||
             numEdgesFirstNode == constants::geometric::numNodesInhaxagon) &&
            (numEdgesSecondNode == constants::geometric::numNodesInPentagon ||
             numEdgesSecondNode == constants::geometric::numNodesInhaxagon) &&
            IsEqual(mesh2d->m_nodes[firstNode].y, mesh2d->m_nodes[secondNode].y, tolerance))
        {
            mesh2d->DeleteEdge(e);
        }
    }

    mesh2d->AdministrateNodesEdges();

    return mesh2d;
}
