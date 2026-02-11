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

#include "MeshKernel/MeshFaceCenters.hpp"
#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Operations.hpp"

meshkernel::Point meshkernel::algo::CircumcenterOfTriangle(const Point& firstNode, const Point& secondNode, const Point& thirdNode, const Projection projection)
{
    const double dx2 = GetDx(firstNode, secondNode, projection);
    const double dy2 = GetDy(firstNode, secondNode, projection);

    const double dx3 = GetDx(firstNode, thirdNode, projection);
    const double dy3 = GetDy(firstNode, thirdNode, projection);

    const double den = dy2 * dx3 - dy3 * dx2;
    double z = 0.0;

    if (std::abs(den) > 0.0)
    {
        z = (dx2 * (dx2 - dx3) + dy2 * (dy2 - dy3)) / den;
    }

    Point circumcenter;

    if (projection == Projection::cartesian)
    {
        circumcenter.x = firstNode.x + 0.5 * (dx3 - z * dy3);
        circumcenter.y = firstNode.y + 0.5 * (dy3 + z * dx3);
    }
    else if (projection == Projection::spherical)
    {
        const double phi = (firstNode.y + secondNode.y + thirdNode.y) * constants::numeric::oneThird;
        const double xf = 1.0 / cos(constants::conversion::degToRad * phi);
        circumcenter.x = firstNode.x + xf * 0.5 * (dx3 - z * dy3) * constants::conversion::radToDeg / constants::geometric::earth_radius;
        circumcenter.y = firstNode.y + 0.5 * (dy3 + z * dx3) * constants::conversion::radToDeg / constants::geometric::earth_radius;
    }
    else if (projection == Projection::sphericalAccurate)
    {
        // compute in case of spherical accurate (comp_circumcenter3D)
    }
    return circumcenter;
}

meshkernel::Point meshkernel::algo::ComputeCircumCenter(const Point& centerOfMass,
                                                        const UInt pointCount,
                                                        const std::array<Point, constants::geometric::maximumNumberOfNodesPerFace>& middlePoints,
                                                        const std::array<Point, constants::geometric::maximumNumberOfNodesPerFace>& normals,
                                                        const Projection projection)
{
    const double eps = constants::geometric::circumcentreTolerance * (projection == Projection::cartesian ? 1.0 : 1.0 / (constants::geometric::earth_radius * constants::conversion::degToRad));

    Point estimatedCircumCenter = centerOfMass;

    for (UInt iter = 0; iter < constants::numeric::MaximumNumberOfCircumcentreIterations; ++iter)
    {
        const Point previousCircumCenter = estimatedCircumCenter;
        for (UInt n = 0; n < pointCount; n++)
        {
            const Vector delta{GetDelta(middlePoints[n], previousCircumCenter, projection)};
            const auto increment = -0.1 * dot(delta, normals[n]);
            AddIncrementToPoint(normals[n], increment, centerOfMass, projection, estimatedCircumCenter);
        }
        if (iter > 0 &&
            abs(estimatedCircumCenter.x - previousCircumCenter.x) < eps &&
            abs(estimatedCircumCenter.y - previousCircumCenter.y) < eps)
        {
            break;
        }
    }

    return estimatedCircumCenter;
}

meshkernel::Point meshkernel::algo::ComputeFaceCircumenter(std::vector<Point>& polygon,
                                                           const std::vector<UInt>& edgesNumFaces,
                                                           const Projection projection)
{
    static constexpr double weightCircumCenter = 1.0; ///< Weight circum center

    std::array<Point, constants::geometric::maximumNumberOfNodesPerFace> middlePoints;
    std::array<Point, constants::geometric::maximumNumberOfNodesPerFace> normals;
    UInt pointCount = 0;

    const auto numNodes = static_cast<UInt>(polygon.size()) - 1;

    Point centerOfMass{0.0, 0.0};

    for (UInt n = 0; n < numNodes; ++n)
    {
        centerOfMass.x += polygon[n].x;
        centerOfMass.y += polygon[n].y;
    }

    centerOfMass /= static_cast<double>(numNodes);

    auto result = centerOfMass;
    if (numNodes == constants::geometric::numNodesInTriangle)
    {
        result = algo::CircumcenterOfTriangle(polygon[0], polygon[1], polygon[2], projection);
    }
    else if (!edgesNumFaces.empty())
    {
        UInt numValidEdges = CountNumberOfValidEdges(edgesNumFaces, numNodes);

        if (numValidEdges > 1)
        {
            ComputeMidPointsAndNormals(polygon, edgesNumFaces, numNodes, middlePoints, normals, pointCount, projection);
            result = algo::ComputeCircumCenter(centerOfMass, pointCount, middlePoints, normals, projection);
        }
    }

    if (weightCircumCenter != 1.0)
    {
        for (UInt n = 0; n < numNodes; ++n)
        {
            polygon[n] = weightCircumCenter * polygon[n] + (1.0 - weightCircumCenter) * centerOfMass;
        }
    }

    // The circumcenter is included in the face, then return the calculated circumcenter
    if (IsPointInPolygonNodes(result, polygon, projection))
    {
        return result;
    }

    // If the circumcenter is not included in the face,
    // the circumcenter will be placed at the intersection between an edge and the segment connecting the mass center with the circumcenter.
    for (UInt n = 0; n < numNodes; ++n)
    {
        const auto nextNode = NextCircularForwardIndex(n, numNodes);

        const auto [areLineCrossing,
                    intersection,
                    crossProduct,
                    intersectionAngle,
                    firstRatio,
                    secondRatio] = AreSegmentsCrossing(centerOfMass, result, polygon[n], polygon[nextNode], false, projection);

        if (areLineCrossing)
        {
            result = intersection;
            break;
        }
    }

    return result;
}

std::vector<meshkernel::Point> meshkernel::algo::ComputeFaceCircumcenters(const Mesh& mesh)
{
    std::vector<Point> faceCenters(mesh.GetNumFaces());
    ComputeFaceCircumcenters(mesh, faceCenters);

    return faceCenters;
}

void meshkernel::algo::ComputeFaceCircumcenters(const Mesh& mesh, std::span<Point> faceCenters)
{
    if (faceCenters.size() != mesh.GetNumFaces())
    {
        throw ConstraintError("array for faceCenters values is not the correct size");
    }

    const auto numFaces = static_cast<int>(mesh.GetNumFaces());

    std::vector<UInt> numEdgeFacesCache;
    numEdgeFacesCache.reserve(constants::geometric::maximumNumberOfEdgesPerFace);
    std::vector<Point> polygonNodesCache;
    polygonNodesCache.reserve(constants::geometric::maximumNumberOfEdgesPerFace);

#pragma omp parallel for private(numEdgeFacesCache, polygonNodesCache)
    for (int f = 0; f < numFaces; f++)
    {

        UInt numberOfInteriorEdges = 0;
        const auto numberOfFaceNodes = mesh.GetNumFaceEdges(f);

        for (UInt n = 0; n < numberOfFaceNodes; ++n)
        {
            if (!mesh.IsEdgeOnBoundary(mesh.m_facesEdges[f][n]))
            {
                numberOfInteriorEdges += 1;
            }
        }

        if (numberOfInteriorEdges == 0)
        {
            faceCenters[f] = mesh.m_facesMassCenters[f];
        }
        else
        {
            // need to account for spherical coordinates. Build a polygon around a face
            mesh.ComputeFaceClosedPolygon(f, polygonNodesCache);
            numEdgeFacesCache.clear();

            for (UInt n = 0; n < numberOfFaceNodes; ++n)
            {
                numEdgeFacesCache.emplace_back(mesh.m_edgesNumFaces[mesh.m_facesEdges[f][n]]);
            }

            faceCenters[f] = algo::ComputeFaceCircumenter(polygonNodesCache, numEdgeFacesCache, mesh.m_projection);
        }
    }
}
