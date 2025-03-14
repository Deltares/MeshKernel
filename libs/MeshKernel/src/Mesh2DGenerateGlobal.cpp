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

#include "MeshKernel/Mesh2DGenerateGlobal.hpp"
#include "MeshKernel/Operations.hpp"
#include "MeshKernel/Polygons.hpp"
#include "MeshKernel/UndoActions/CompoundUndoAction.hpp"
#include "MeshKernel/Utilities/RTreeFactory.hpp"
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

    for (auto i = static_cast<int>(mesh.GetNumNodes() - 1); i >= 0; --i)
    {
        if (IsEqual(x, mesh.Node(i), tolerance))
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
            auto [edgeId, nodeInsertionAction] = mesh.InsertNode(p);
            nodeIndices[n] = edgeId;
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
            auto [edgeId, connectionAction] = mesh.ConnectNodes(firstNodeIndex, secondNodeIndex);
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

            // TODO before merging with master is it possible to change which points get deleted in the Mesh::MergePointsInPolygon
            if (points[2].x < 180.0)
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

    // The merge action can be ignored in this case because we will not need to undo any merge operation
    [[maybe_unused]] auto mergeAction = mesh2d->MergeNodesInPolygon(polygons, mergingDistance);

    constexpr double tolerance = 1.0e-6;
    for (UInt e = 0; e < mesh2d->GetNumEdges(); e++)
    {
        const auto& [firstNode, secondNode] = mesh2d->GetEdge(e);

        if (firstNode == constants::missing::uintValue || secondNode == constants::missing::uintValue)
        {
            continue;
        }

        const auto numEdgesFirstNode = mesh2d->GetNumNodesEdges(firstNode);
        const auto numEdgesSecondNode = mesh2d->GetNumNodesEdges(secondNode);
        if ((numEdgesFirstNode == constants::geometric::numNodesInPentagon ||
             numEdgesFirstNode == constants::geometric::numNodesInHexagon) &&
            (numEdgesSecondNode == constants::geometric::numNodesInPentagon ||
             numEdgesSecondNode == constants::geometric::numNodesInHexagon) &&
            IsEqual(mesh2d->Node(firstNode).y, mesh2d->Node(secondNode).y, tolerance))
        {
            [[maybe_unused]] auto action = mesh2d->DeleteEdge(e);
        }
    }

    // A newly created grid should have no invalid nodes nor edges.
    // Delete any invalid node and edges that have been generated during calculation of the grid.
    mesh2d->DeleteInvalidNodesAndEdges();
    mesh2d->AdministrateNodesEdges();

    return mesh2d;
}
