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

// std::tuple<double, meshkernel::Point, meshkernel::TraversalDirection> meshkernel::Polygon::FaceAreaAndCenterOfMass(const std::vector<Point>& polygon, const Projection projection)
// {

//     if (polygon.size() < constants::geometric::numNodesInTriangle)
//     {
//         throw std::invalid_argument("FaceAreaAndCenterOfMass: The polygon has less than 3 unique nodes.");
//     }

//     double area = 0.0;

//     const double minArea = 1e-8;
//     const auto numberOfPointsOpenedPolygon = static_cast<UInt>(polygon.size()) - 1;

//     const double updateStepSize = 0.1;

//     Point centreOfMass(0.0, 0.0);

//     for (UInt n = 0; n < numberOfPointsOpenedPolygon; ++n)
//     {
//         centreOfMass += polygon[n];
//     }

//     centreOfMass *= 1.0 / static_cast<double>(numberOfPointsOpenedPolygon);
//     // Will be non-unity for spherical coordinates only
//     const double xTransformation = projection == Projection::cartesian ? 1.0 : 1.0 / std::cos(centreOfMass.y * constants::conversion::degToRad);
//     const double circumcentreTolerance = constants::geometric::circumcentreTolerance * (projection == Projection::cartesian ? 1.0 : 1.0 / (constants::geometric::earth_radius * constants::conversion::degToRad));

//     if (numberOfPointsOpenedPolygon == constants::geometric::numNodesInTriangle)
//     {

//         // Point midPoint1 = 0.5 * (polygon[0] + polygon[1]);
//         // Point midPoint2 = 0.5 * (polygon[1] + polygon[2]);
//         // Point midPoint3 = 0.5 * (polygon[2] + polygon[0]);

//         // Vector edgeVector1 = static_cast<Vector>(NormalVector(polygon[0], polygon[1], midPoint1, projection));
//         // Vector edgeVector2 = static_cast<Vector>(NormalVector(polygon[1], polygon[2], midPoint2, projection));
//         // Vector edgeVector3 = static_cast<Vector>(NormalVector(polygon[2], polygon[0], midPoint3, projection));

//         // edgeVector1.normalise();
//         // edgeVector2.normalise();
//         // edgeVector3.normalise();

//         // Vector edgeVectorSum = edgeVector1 + edgeVector2 + edgeVector3;

//         // double edgeVectorSumLength = edgeVectorSum.length();

//         // for (UInt i = 1; i <= MaximumNumberOfCircumcentreIterations; ++i)
//         // {
//         //     Vector delta1 = GetDelta(midPoint1, centreOfMass, projection);
//         //     Vector delta2 = GetDelta(midPoint2, centreOfMass, projection);
//         //     Vector delta3 = GetDelta(midPoint3, centreOfMass, projection);

//         //     double ds = dot(delta1, edgeVector1) + dot(delta2, edgeVector2) + dot(delta3, edgeVector3);

//         //     if (projection != Projection::cartesian)
//         //     {
//         //         ds *= constants::conversion::radToDeg * constants::geometric::inverse_earth_radius;
//         //     }

//         //     centreOfMass.x -= updateStepSize * ds * edgeVectorSum.x() * xTransformation;
//         //     centreOfMass.y -= updateStepSize * ds * edgeVectorSum.y();

//         //     if (ds * edgeVectorSumLength < circumcentreTolerance || i == MaximumNumberOfCircumcentreIterations)
//         //     {
//         //         break;
//         //     }
//         // }

//         Vector delta2 = GetDelta(polygon[0], polygon[1], projection);
//         Vector delta3 = GetDelta(polygon[0], polygon[2], projection);

//         double den = delta2.y() * delta3.x() - delta3.y() * delta2.x();
//         double correction = 0.0;

//         std::cout << "points: " << polygon[0].x << ", " << polygon[0].y << " -- "
//                   << polygon[1].x << ", " << polygon[1].y << " -- "
//                   << polygon[2].x << ", " << polygon[2].y << " -- "
//                   << std::endl;

//         if (den != 0.0)
//         {
//             correction = (delta2.x() * (delta2.x() - delta3.x()) + delta2.y() * (delta2.y() - delta3.y())) / den;
//         }

//         std::cout << "triangle average centre: " << centreOfMass.x << ", " << centreOfMass.y << "  " << den << "  " << correction << std::endl;

//         if (projection == Projection::cartesian)
//         {
//             centreOfMass.x = polygon[0].x + 0.5 * (delta3.x() - correction * delta3.y());
//             centreOfMass.y = polygon[0].y + 0.5 * (delta3.y() + correction * delta3.x());

//             // xz = x(1) + 0.5d0 * (dx3 - z * dy3)
//             // yz = y(1) + 0.5d0 * (dy3 + z * dx3)
//         }
//         else
//         {
//             double angle = (polygon[0].y + polygon[1].y + polygon[2].y) / 3.0;
//             double xf = 1.0 / std::cos(angle * constants::conversion::degToRad);

//             centreOfMass.x = polygon[0].x + xf * 0.5 * (delta3.x() - correction * delta3.y()) * constants::conversion::radToDeg * constants::geometric::inverse_earth_radius;
//             centreOfMass.y = polygon[0].y + 0.5 * (delta3.y() + correction * delta3.x()) * constants::conversion::radToDeg * constants::geometric::inverse_earth_radius;
//         }

//         std::cout << "triangle circum centre: " << centreOfMass.x << ", " << centreOfMass.y << std::endl;
//     }
//     else if (numberOfPointsOpenedPolygon == constants::geometric::numNodesInQuadrilateral)
//     {

//         Point midPoint1 = 0.5 * (polygon[0] + polygon[1]);
//         Point midPoint2 = 0.5 * (polygon[1] + polygon[2]);
//         Point midPoint3 = 0.5 * (polygon[2] + polygon[3]);
//         Point midPoint4 = 0.5 * (polygon[3] + polygon[0]);

//         Vector edgeVector1 = static_cast<Vector>(NormalVector(polygon[0], polygon[1], midPoint1, projection));
//         Vector edgeVector2 = static_cast<Vector>(NormalVector(polygon[1], polygon[2], midPoint2, projection));
//         Vector edgeVector3 = static_cast<Vector>(NormalVector(polygon[2], polygon[3], midPoint3, projection));
//         Vector edgeVector4 = static_cast<Vector>(NormalVector(polygon[3], polygon[0], midPoint4, projection));

//         edgeVector1.normalise();
//         edgeVector2.normalise();
//         edgeVector3.normalise();
//         edgeVector4.normalise();

//         Vector edgeVectorSum = edgeVector1 + edgeVector2 + edgeVector3 + edgeVector4;

//         if (projection != Projection::cartesian)
//         {
//             edgeVectorSum.x() *= constants::conversion::radToDeg * constants::geometric::inverse_earth_radius;
//         }

//         double edgeVectorSumLength = edgeVectorSum.length();

//         for (UInt i = 1; i <= constants::numeric::MaximumNumberOfCircumcentreIterations; ++i)
//         {
//             Vector delta1 = GetDelta(midPoint1, centreOfMass, projection);
//             Vector delta2 = GetDelta(midPoint2, centreOfMass, projection);
//             Vector delta3 = GetDelta(midPoint3, centreOfMass, projection);
//             Vector delta4 = GetDelta(midPoint4, centreOfMass, projection);

//             double ds = dot(delta1, edgeVector1) + dot(delta2, edgeVector2) + dot(delta3, edgeVector3) + dot(delta4, edgeVector4);

//             if (projection != Projection::cartesian)
//             {
//                 ds *= constants::conversion::radToDeg * constants::geometric::inverse_earth_radius;
//             }

//             centreOfMass.x -= updateStepSize * ds * edgeVectorSum.x() * xTransformation;
//             centreOfMass.y -= updateStepSize * ds * edgeVectorSum.y();

//             if (ds * edgeVectorSumLength < circumcentreTolerance)
//             {
//                 break;
//             }
//         }

//         std::cout << "quadrilateral circum centre: " << centreOfMass.x << ", " << centreOfMass.y << std::endl;
//     }
//     else
//     {
//         for (UInt j = 1; j <= constants::numeric::MaximumNumberOfCircumcentreIterations; ++j)
//         {
//             Vector edgeVectorSum(0.0, 0.0);
//             double ds = 0.0;

//             for (UInt i = 0; i < numberOfPointsOpenedPolygon; ++i)
//             {
//                 const auto nextNode = NextCircularForwardIndex(i, numberOfPointsOpenedPolygon);

//                 Point midPoint = 0.5 * (polygon[i] + polygon[nextNode]);
//                 Vector edgeVector = static_cast<Vector>(NormalVector(polygon[i], polygon[nextNode], midPoint, projection));
//                 Vector delta = GetDelta(midPoint, centreOfMass, projection);

//                 edgeVector.normalise();
//                 ds += dot(delta, edgeVector);
//                 edgeVectorSum += edgeVector;
//             }

//             if (projection != Projection::cartesian)
//             {
//                 ds *= constants::conversion::radToDeg * constants::geometric::inverse_earth_radius;
//             }

//             centreOfMass.x -= updateStepSize * ds * edgeVectorSum.x() * xTransformation;
//             centreOfMass.y -= updateStepSize * ds * edgeVectorSum.y();

//             if (j > 1 && ds * edgeVectorSum.length() < circumcentreTolerance)
//             {
//                 break;
//             }
//         }
//     }

//     area = ComputeArea(polygon, projection);
//     TraversalDirection direction = area > 0.0 ? TraversalDirection::AntiClockwise : TraversalDirection::Clockwise;

//     area = std::abs(area) < minArea ? minArea : area;

//     return {std::abs(area), centreOfMass, direction};
// }
