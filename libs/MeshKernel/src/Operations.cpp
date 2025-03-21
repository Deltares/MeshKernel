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

#include "MeshKernel/Operations.hpp"
#include "MeshKernel/Cartesian3DPoint.hpp"
#include "MeshKernel/Mesh.hpp"

namespace meshkernel::impl
{

    /// @brief Checks if a point is in polygonNodes for Cartesian and spherical coordinate systems, using the winding number method
    bool IsPointInPolygonNodesCartesian(const Point& point,
                                        const std::vector<Point>& polygonNodes,
                                        UInt startNode,
                                        UInt endNode)
    {
        int windingNumber = 0;
        for (auto n = startNode; n < endNode; n++)
        {

            const auto crossProductValue = crossProduct(polygonNodes[n], polygonNodes[n + 1], polygonNodes[n], point, Projection::cartesian);

            if (IsEqual(crossProductValue, 0.0))
            {
                // point on the line
                return true;
            }

            if (polygonNodes[n].y <= point.y) // an upward crossing
            {
                if (polygonNodes[n + 1].y > point.y && crossProductValue > 0.0)

                {
                    ++windingNumber; // have  a valid up intersect
                }
            }
            else
            {
                if (polygonNodes[n + 1].y <= point.y && crossProductValue < 0.0) // a downward crossing
                {

                    --windingNumber; // have  a valid down intersect
                }
            }
        }

        return windingNumber != 0;
    }

    /// @brief Checks if a point is in polygonNodes for accurate spherical coordinate system, using the winding number method
    bool IsPointInPolygonNodesSphericalAccurate(const Point& point,
                                                const std::vector<Point>& polygonNodes,
                                                Point polygonCenter,
                                                UInt startNode,
                                                UInt endNode)
    {
        const auto currentPolygonSize = endNode - startNode + 1;

        // get 3D polygon coordinates
        std::vector<Cartesian3DPoint> cartesian3DPoints;
        cartesian3DPoints.reserve(currentPolygonSize);
        for (UInt i = 0; i < currentPolygonSize; ++i)
        {
            cartesian3DPoints.emplace_back(SphericalToCartesian3D(polygonNodes[startNode + i]));
        }

        // enlarge around polygon
        const double enlargementFactor = 1.000001;
        const Cartesian3DPoint polygonCenterCartesian3D{SphericalToCartesian3D(polygonCenter)};
        for (UInt i = 0; i < currentPolygonSize; ++i)
        {
            cartesian3DPoints[i] = polygonCenterCartesian3D + enlargementFactor * (cartesian3DPoints[i] - polygonCenterCartesian3D);
        }

        // convert point
        const Cartesian3DPoint pointCartesian3D{SphericalToCartesian3D(point)};

        // get test direction: e_lambda
        const double lambda = point.x * constants::conversion::degToRad;
        const Cartesian3DPoint ee{-std::sin(lambda), std::cos(lambda), 0.0};
        int inside = 0;

        // loop over the polygon nodes
        for (UInt i = 0; i < currentPolygonSize - 1; ++i)
        {
            const auto nextNode = NextCircularForwardIndex(i, currentPolygonSize);
            const auto xiXxip1 = VectorProduct(cartesian3DPoints[i], cartesian3DPoints[nextNode]);
            const auto xpXe = VectorProduct(pointCartesian3D, ee);

            const auto D = InnerProduct(xiXxip1, ee);
            double zeta = 0.0;
            double xi = 0.0;
            double eta = 0.0;
            if (std::abs(D) > 0.0)
            {

                xi = -InnerProduct(xpXe, cartesian3DPoints[nextNode]) / D;
                eta = InnerProduct(xpXe, cartesian3DPoints[i]) / D;
                zeta = -InnerProduct(xiXxip1, pointCartesian3D) / D;
            }

            if (IsEqual(zeta, 0.0))
            {
                return true;
            }

            if (xi >= 0.0 && eta >= 0.0 && zeta >= 0.0)
            {
                inside = 1 - inside;
            }
        }

        return inside == 1;
    }

} // namespace meshkernel::impl

namespace meshkernel
{

    std::vector<std::pair<UInt, UInt>> FindIndices(const std::vector<Point>& vec,
                                                   size_t start,
                                                   size_t end,
                                                   double separator)
    {
        std::vector<std::pair<UInt, UInt>> result;

        if (vec.empty())
        {
            return result;
        }

        // set an invalid index
        if (start > vec.size() || end > vec.size())
        {
            return result;
        }

        bool inRange = false;
        UInt startRange = 0;
        for (auto n = start; n < end; n++)
        {
            if (!IsEqual(vec[n].x, separator) && !inRange)
            {
                startRange = static_cast<UInt>(n);
                inRange = true;
            }
            if (IsEqual(vec[n].x, separator) && inRange)
            {
                result.emplace_back(startRange, static_cast<UInt>(n - 1));
                inRange = false;
            }
        }

        // in case no separator was found
        if (inRange)
        {
            result.emplace_back(startRange, static_cast<UInt>(vec.size() - 1));
        }

        return result;
    }

    UInt InvalidPointCount(const std::vector<Point>& points)
    {
        if (points.size() > 0)
        {
            return InvalidPointCount(points, 0, points.size() - 1);
        }
        else
        {
            return 0;
        }
    }

    UInt InvalidPointCount(const std::vector<Point>& points, size_t start, size_t end)
    {
        UInt count = 0;

        for (size_t i = start; i <= end; ++i)
        {
            if (!points[i].IsValid())
            {
                ++count;
            }
        }

        return count;
    }

    UInt NextCircularForwardIndex(UInt currentIndex, UInt size)
    {
        if (size == 0)
        {
            throw ConstraintError("Invalid rotation range ");
        }

        if (currentIndex >= size)
        {
            throw ConstraintError("Index is out of range: {} not in [0 .. {}]", currentIndex, size - 1);
        }

        return currentIndex == size - 1 ? 0 : currentIndex + 1;
    }

    UInt NextCircularBackwardIndex(UInt currentIndex, UInt size)
    {
        if (size == 0)
        {
            throw ConstraintError("Invalid rotation range ");
        }

        if (currentIndex >= size)
        {
            throw ConstraintError("Index is out of range: {} not in [0 .. {}]", currentIndex, size - 1);
        }

        return currentIndex == 0 ? size - 1 : currentIndex - 1;
    }

    bool IsPointOnPole(const Point& point)
    {
        return std::abs(std::abs(point.y) - 90.0) < constants::geometric::absLatitudeAtPoles;
    }

    double crossProduct(const Point& firstSegmentFirstPoint, const Point& firstSegmentSecondPoint, const Point& secondSegmentFirstPoint, const Point& secondSegmentSecondPoint, const Projection& projection)
    {
        const auto dx1 = GetDx(firstSegmentFirstPoint, firstSegmentSecondPoint, projection);
        const auto dy1 = GetDy(firstSegmentFirstPoint, firstSegmentSecondPoint, projection);
        const auto dx2 = GetDx(secondSegmentFirstPoint, secondSegmentSecondPoint, projection);
        const auto dy2 = GetDy(secondSegmentFirstPoint, secondSegmentSecondPoint, projection);
        return dx1 * dy2 - dy1 * dx2;
    }

    bool IsPointInTriangle(const Point& point,
                           const std::span<const Point> triangleNodes,
                           const Projection& projection)
    {
        if (triangleNodes.empty())
        {
            return true;
        }

        bool isInTriangle = false;

        if (projection == Projection::cartesian || projection == Projection::spherical)
        {
            int windingNumber = 0;

            for (UInt n = 0; n < constants::geometric::numNodesInTriangle; n++)
            {
                UInt endIndex = n == 2 ? 0 : n + 1;

                const auto crossProductValue = crossProduct(triangleNodes[n], triangleNodes[endIndex], triangleNodes[n], point, Projection::cartesian);

                if (IsEqual(crossProductValue, 0.0))
                {
                    // point on the line
                    return true;
                }

                if (triangleNodes[n].y <= point.y) // an upward crossing
                {
                    if (triangleNodes[endIndex].y > point.y && crossProductValue > 0.0)
                    {
                        ++windingNumber; // have  a valid up intersect
                    }
                }
                else
                {
                    if (triangleNodes[endIndex].y <= point.y && crossProductValue < 0.0) // a downward crossing
                    {
                        --windingNumber; // have  a valid down intersect
                    }
                }
            }

            isInTriangle = windingNumber == 0 ? false : true;
        }

        if (projection == Projection::sphericalAccurate)
        {
            std::vector<Point> triangleCopy;
            triangleCopy.reserve(constants::geometric::numNodesInTriangle + 1);
            Point centre(0.0, 0.0);

            for (UInt i = 0; i < constants::geometric::numNodesInTriangle; ++i)
            {
                centre += triangleNodes[i];
                triangleCopy.emplace_back(triangleNodes[i]);
            }

            triangleCopy.emplace_back(triangleNodes[0]);

            centre *= 1.0 / 3.0;

            isInTriangle = IsPointInPolygonNodes(point, triangleCopy, projection, centre);
        }

        return isInTriangle;
    }

    bool IsPointInPolygonNodes(const Point& point,
                               const std::vector<Point>& polygonNodes,
                               const Projection& projection,
                               Point polygonCenter,
                               UInt startNode,
                               UInt endNode)
    {

        if (polygonNodes.empty())
        {
            return true;
        }

        if (startNode == constants::missing::uintValue && endNode == constants::missing::uintValue)
        {
            startNode = 0;
            endNode = static_cast<UInt>(polygonNodes.size()) - 1; // closed polygon
        }

        if (endNode <= startNode)
        {
            return true;
        }

        const auto currentPolygonSize = endNode - startNode + 1;
        if (currentPolygonSize < constants::geometric::numNodesInTriangle || polygonNodes.size() < currentPolygonSize)
        {
            return false;
        }
        if (polygonNodes[startNode] != polygonNodes[endNode])
        {
            return false;
        }

        if (const auto boundingBox = BoundingBox(polygonNodes); !boundingBox.Contains(point))
        {
            return false;
        }

        if (projection == Projection::cartesian || projection == Projection::spherical)
        {
            return impl::IsPointInPolygonNodesCartesian(point, polygonNodes, startNode, endNode);
        }
        else
        {
            return impl::IsPointInPolygonNodesSphericalAccurate(point, polygonNodes, polygonCenter, startNode, endNode);
        }
    }

    bool IsPointInPolygonNodes(const Point& point,
                               const std::vector<Point>& polygonNodes,
                               const Projection& projection,
                               const BoundingBox& boundingBox,
                               Point polygonCenter,
                               UInt startNode,
                               UInt endNode)
    {

        if (polygonNodes.empty())
        {
            return true;
        }

        if (startNode == constants::missing::uintValue && endNode == constants::missing::uintValue)
        {
            startNode = 0;
            endNode = static_cast<UInt>(polygonNodes.size()) - 1; // closed polygon
        }

        if (endNode <= startNode)
        {
            return true;
        }

        if (const auto currentPolygonSize = endNode - startNode + 1; currentPolygonSize < constants::geometric::numNodesInTriangle || polygonNodes.size() < currentPolygonSize)
        {
            return false;
        }

        if (polygonNodes[startNode] != polygonNodes[endNode])
        {
            return false;
        }

        if (!boundingBox.Contains(point))
        {
            return false;
        }

        if (projection == Projection::cartesian || projection == Projection::spherical)
        {
            return impl::IsPointInPolygonNodesCartesian(point, polygonNodes, startNode, endNode);
        }
        else
        {
            return impl::IsPointInPolygonNodesSphericalAccurate(point, polygonNodes, polygonCenter, startNode, endNode);
        }
    }

    void ComputeThreeBaseComponents(const Point& point, std::array<double, 3>& exxp, std::array<double, 3>& eyyp, std::array<double, 3>& ezzp)
    {
        const double phi0 = point.y * constants::conversion::degToRad;
        const double lambda0 = point.x * constants::conversion::degToRad;

        exxp[0] = cos(phi0) * cos(lambda0);
        exxp[1] = cos(phi0) * sin(lambda0);
        exxp[2] = sin(phi0);

        eyyp[0] = -sin(lambda0);
        eyyp[1] = cos(lambda0);
        eyyp[2] = 0.0;

        ezzp[0] = -sin(phi0) * cos(lambda0);
        ezzp[1] = -sin(phi0) * sin(lambda0);
        ezzp[2] = cos(phi0);
    }

    void ComputeTwoBaseComponents(const Point& point, std::array<double, 3>& elambda, std::array<double, 3>& ephi)
    {
        const double phi0 = point.y * constants::conversion::degToRad;
        const double lambda0 = point.x * constants::conversion::degToRad;

        elambda[0] = -sin(lambda0);
        elambda[1] = cos(lambda0);
        elambda[2] = 0.0;

        ephi[0] = -sin(phi0) * cos(lambda0);
        ephi[1] = -sin(phi0) * sin(lambda0);
        ephi[2] = cos(phi0);
    }

    Vector GetDelta(const Point& firstPoint, const Point& secondPoint, const Projection& projection)
    {
        if (projection == Projection::cartesian)
        {
            return GetDeltaCartesian(firstPoint, secondPoint);
        }

        // TODO some performance can be gained here, by combining the computing of dx and dy
        return Vector(GetDx(firstPoint, secondPoint, projection), GetDy(firstPoint, secondPoint, projection));
    }

    Vector ComputeNormalToline(const Point& start, const Point& end, const Projection projection)
    {
        Vector direction = GetDelta(start, end, projection);
        direction.normalise();
        Vector normal(-direction.y(), direction.x());
        return normal;
    }

    double GetDx(const Point& firstPoint, const Point& secondPoint, const Projection& projection)
    {

        if (projection == Projection::cartesian)
        {
            const double delta = secondPoint.x - firstPoint.x;
            return delta;
        }
        if (projection == Projection::spherical || projection == Projection::sphericalAccurate)
        {
            const bool isFirstPointOnPole = IsPointOnPole(firstPoint);
            const bool isSecondPointOnPole = IsPointOnPole(secondPoint);
            if ((isFirstPointOnPole && !isSecondPointOnPole) || (!isFirstPointOnPole && isSecondPointOnPole))
            {
                return 0.0;
            }
            double firstPointX = firstPoint.x;
            double secondPointX = secondPoint.x;
            if (firstPointX - secondPointX > 180.0)
            {
                firstPointX -= 360.0;
            }
            else if (firstPointX - secondPointX < -180.0)
            {
                firstPointX += 360.0;
            }

            firstPointX = firstPointX * constants::conversion::degToRad;
            secondPointX = secondPointX * constants::conversion::degToRad;
            const double firstPointY = firstPoint.y * constants::conversion::degToRad;
            const double secondPointY = secondPoint.y * constants::conversion::degToRad;
            const double cosPhi = cos(0.5 * (firstPointY + secondPointY));
            const double dx = constants::geometric::earth_radius * cosPhi * (secondPointX - firstPointX);
            return dx;
        }
        return constants::missing::doubleValue;
    }

    double GetDy(const Point& firstPoint, const Point& secondPoint, const Projection& projection)
    {

        if (projection == Projection::cartesian)
        {
            const double delta = secondPoint.y - firstPoint.y;
            return delta;
        }
        if (projection == Projection::spherical || projection == Projection::sphericalAccurate)
        {
            const double firstPointY = firstPoint.y * constants::conversion::degToRad;
            const double secondPointY = secondPoint.y * constants::conversion::degToRad;
            const double dy = constants::geometric::earth_radius * (secondPointY - firstPointY);
            return dy;
        }
        return constants::missing::doubleValue;
    }

    double OuterProductTwoSegments(const Point& firstPointFirstSegment, const Point& secondPointFirstSegment,
                                   const Point& firstPointSecondSegment, const Point& secondPointSecondSegment, const Projection& projection)
    {
        if (projection == Projection::sphericalAccurate)
        {
            const Cartesian3DPoint firstPointFirstSegmentCartesian{SphericalToCartesian3D(firstPointFirstSegment)};
            const auto xx1 = firstPointFirstSegmentCartesian.x;
            const auto yy1 = firstPointFirstSegmentCartesian.y;
            const auto zz1 = firstPointFirstSegmentCartesian.z;

            const Cartesian3DPoint secondPointFirstSegmentCartesian{SphericalToCartesian3D(secondPointFirstSegment)};
            const auto xx2 = secondPointFirstSegmentCartesian.x;
            const auto yy2 = secondPointFirstSegmentCartesian.y;
            const auto zz2 = secondPointFirstSegmentCartesian.z;

            const Cartesian3DPoint firstPointSecondSegmentCartesian{SphericalToCartesian3D(firstPointSecondSegment)};
            const auto xx3 = firstPointSecondSegmentCartesian.x;
            const auto yy3 = firstPointSecondSegmentCartesian.y;
            const auto zz3 = firstPointSecondSegmentCartesian.z;

            const Cartesian3DPoint secondPointSecondSegmentCartesian{SphericalToCartesian3D(secondPointSecondSegment)};
            const auto xx4 = secondPointSecondSegmentCartesian.x;
            const auto yy4 = secondPointSecondSegmentCartesian.y;
            const auto zz4 = secondPointSecondSegmentCartesian.z;

            const double vxx = (yy2 - yy1) * (zz4 - zz3) - (zz2 - zz1) * (yy4 - yy3);
            const double vyy = (zz2 - zz1) * (xx4 - xx3) - (xx2 - xx1) * (zz4 - zz3);
            const double vzz = (xx2 - xx1) * (yy4 - yy3) - (yy2 - yy1) * (xx4 - xx3);

            double result = std::sqrt(vxx * vxx + vyy * vyy + vzz * vzz);

            // check if vector is pointing outwards of earth
            if (vxx * xx1 + vyy * yy1 + vzz * zz1 < 0.0)
            {
                result = -result;
            }
            return result;
        }

        // cartesian and spherical
        if (projection == Projection::cartesian || projection == Projection::spherical)
        {
            const double dx1 = GetDx(firstPointFirstSegment, secondPointFirstSegment, projection);
            const double dx2 = GetDx(firstPointSecondSegment, secondPointSecondSegment, projection);

            const double dy1 = GetDy(firstPointFirstSegment, secondPointFirstSegment, projection);
            const double dy2 = GetDy(firstPointSecondSegment, secondPointSecondSegment, projection);

            return dx1 * dy2 - dy1 * dx2;
        }

        return constants::missing::doubleValue;
    }

    Point ComputeMiddlePoint(const Point& firstPoint, const Point& secondPoint, const Projection& projection)
    {
        if (!firstPoint.IsValid() || !secondPoint.IsValid())
        {
            return {constants::missing::doubleValue, constants::missing::doubleValue};
        }
        if (projection == Projection::sphericalAccurate)
        {
            const Cartesian3DPoint firstPointCartesianCoordinates{SphericalToCartesian3D(firstPoint)};
            const Cartesian3DPoint secondPointCartesianCoordinates{SphericalToCartesian3D(secondPoint)};

            Cartesian3DPoint middleCartesianPointCoordinate{constants::missing::doubleValue, constants::missing::doubleValue, constants::missing::doubleValue};
            middleCartesianPointCoordinate.x = 0.5 * (firstPointCartesianCoordinates.x + secondPointCartesianCoordinates.x);
            middleCartesianPointCoordinate.y = 0.5 * (firstPointCartesianCoordinates.y + secondPointCartesianCoordinates.y);
            const double referenceLongitude = std::max(firstPoint.x, secondPoint.x);
            const auto result = Cartesian3DToSpherical(middleCartesianPointCoordinate, referenceLongitude);
            return result;
        }

        // cartesian and spherical
        if (projection == Projection::cartesian || projection == Projection::spherical)
        {
            return {(firstPoint + secondPoint) * 0.5};
        }
        return {constants::missing::doubleValue, constants::missing::doubleValue};
    }

    Point ComputeMiddlePointAccountingForPoles(const Point& firstPoint, const Point& secondPoint, const Projection& projection)
    {
        if (!firstPoint.IsValid() || !secondPoint.IsValid())
        {
            return {constants::missing::doubleValue, constants::missing::doubleValue};
        }

        if (projection == Projection::spherical || projection == Projection::sphericalAccurate)
        {
            Point result((firstPoint + secondPoint) * 0.5);
            const auto isFirstNodeOnPole = IsPointOnPole(firstPoint);
            const auto isSecondNodeOnPole = IsPointOnPole(secondPoint);

            if (isFirstNodeOnPole && !isSecondNodeOnPole)
            {
                result.x = secondPoint.x;
            }
            else if (!isFirstNodeOnPole && isSecondNodeOnPole)
            {
                result.x = firstPoint.x;
            }
            else
            {
                const auto maxx = std::max(firstPoint.x, secondPoint.x);
                const auto minx = std::min(firstPoint.x, secondPoint.x);

                if (maxx - minx > 180.0)
                {
                    result.x = result.x + 180.0;
                }
            }
            return result;
        }

        // cartesian
        if (projection == Projection::cartesian)
        {
            return {(firstPoint + secondPoint) * 0.5};
        }
        return {constants::missing::doubleValue, constants::missing::doubleValue};
    }

    Point NormalVector(const Point& firstPoint, const Point& secondPoint, const Point& insidePoint, const Projection& projection)
    {
        if (!firstPoint.IsValid() || !secondPoint.IsValid())
        {
            return {constants::missing::doubleValue, constants::missing::doubleValue};
        }

        if (projection == Projection::sphericalAccurate)
        {
            const Cartesian3DPoint firstPointCartesianCoordinates{SphericalToCartesian3D(firstPoint)};
            const Cartesian3DPoint secondPointCartesianCoordinates{SphericalToCartesian3D(secondPoint)};

            std::array<double, 3> elambda{0.0, 0.0, 0.0};
            std::array<double, 3> ephi{0.0, 0.0, 0.0};
            ComputeTwoBaseComponents(insidePoint, elambda, ephi);

            const double dx = (secondPointCartesianCoordinates.x - firstPointCartesianCoordinates.x) * elambda[0] +
                              (secondPointCartesianCoordinates.y - firstPointCartesianCoordinates.y) * elambda[1] +
                              (secondPointCartesianCoordinates.z - firstPointCartesianCoordinates.z) * elambda[2];

            const double dy = (secondPointCartesianCoordinates.x - firstPointCartesianCoordinates.x) * ephi[0] +
                              (secondPointCartesianCoordinates.y - firstPointCartesianCoordinates.y) * ephi[1] +
                              (secondPointCartesianCoordinates.z - firstPointCartesianCoordinates.z) * ephi[2];

            const double squaredDistance = dx * dx + dy * dy;
            Point result{constants::missing::doubleValue, constants::missing::doubleValue};
            if (squaredDistance > 0.0)
            {
                const double distance = sqrt(squaredDistance);
                result.x = dx / distance;
                result.y = dy / distance;
            }
            return result;
        }

        // cartesian and spherical
        if (projection == Projection::cartesian || projection == Projection::spherical)
        {
            const double dx = GetDx(firstPoint, secondPoint, projection);
            const double dy = GetDy(firstPoint, secondPoint, projection);
            const double squaredDistance = dx * dx + dy * dy;
            Point result{constants::missing::doubleValue, constants::missing::doubleValue};
            if (squaredDistance > 0.0)
            {
                const double distance = sqrt(squaredDistance);
                result.x = dx / distance;
                result.y = dy / distance;
            }
            return result;
        }
        return {constants::missing::doubleValue, constants::missing::doubleValue};
    }

    void TransformGlobalVectorToLocal(const Point& reference, const Point& globalCoordinates, const Point& globalComponents, const Projection& projection, Point& localComponents)
    {
        if (projection == Projection::sphericalAccurate)
        {
            std::array<double, 3> exxp{0.0, 0.0, 0.0};
            std::array<double, 3> eyyp{0.0, 0.0, 0.0};
            std::array<double, 3> ezzp{0.0, 0.0, 0.0};
            ComputeThreeBaseComponents(reference, exxp, eyyp, ezzp);

            // get the 3D coordinate
            const Cartesian3DPoint globalCoordinatesCartesian{SphericalToCartesian3D(globalCoordinates)};

            // project to rotated frame
            Cartesian3DPoint globalCoordinatesCartesianRotated;
            globalCoordinatesCartesianRotated.x = exxp[0] * globalCoordinatesCartesian.x + exxp[1] * globalCoordinatesCartesian.y + exxp[2] * globalCoordinatesCartesian.z;
            globalCoordinatesCartesianRotated.y = eyyp[0] * globalCoordinatesCartesian.x + eyyp[1] * globalCoordinatesCartesian.y + eyyp[2] * globalCoordinatesCartesian.z;
            globalCoordinatesCartesianRotated.z = ezzp[0] * globalCoordinatesCartesian.x + ezzp[1] * globalCoordinatesCartesian.y + ezzp[2] * globalCoordinatesCartesian.z;

            // Compute global base vectors at other point in 3D(xx, yy, zz) frame
            std::array<double, 3> elambda{0.0, 0.0, 0.0};
            std::array<double, 3> ephi{0.0, 0.0, 0.0};
            ComputeTwoBaseComponents(globalCoordinates, elambda, ephi);

            const double vxx = globalComponents.x * elambda[0] + globalComponents.y * ephi[0];
            const double vyy = globalComponents.x * elambda[1] + globalComponents.y * ephi[1];
            const double vzz = globalComponents.x * elambda[2] + globalComponents.y * ephi[2];

            // transform to local spherical coordinates
            const auto globalCoordinatesToLocal = Cartesian3DToSpherical(globalCoordinatesCartesianRotated, reference.x);

            // compute base vectors at other point in rotated 3D(xxp, yyp, zzp) frame
            std::array<double, 3> elambdap{0.0, 0.0, 0.0};
            std::array<double, 3> ephip{0.0, 0.0, 0.0};
            ComputeTwoBaseComponents(globalCoordinatesToLocal, elambdap, ephip);

            // compute local base vectors in(xx, yy, zz) frame
            std::array<double, 3> elambdaloc;
            elambdaloc[0] = exxp[0] * elambdap[0] + eyyp[0] * elambdap[1] + ezzp[0] * elambda[2];
            elambdaloc[1] = exxp[1] * elambdap[0] + eyyp[1] * elambdap[1] + ezzp[1] * elambda[2];
            elambdaloc[2] = exxp[2] * elambdap[0] + eyyp[2] * elambdap[1] + ezzp[2] * elambda[2];

            std::array<double, 3> ephiloc;
            ephiloc[0] = exxp[0] * ephip[0] + eyyp[0] * ephip[1] + ezzp[0] * ephip[2];
            ephiloc[1] = exxp[1] * ephip[0] + eyyp[1] * ephip[1] + ezzp[1] * ephip[2];
            ephiloc[2] = exxp[2] * ephip[0] + eyyp[2] * ephip[1] + ezzp[2] * ephip[2];

            // compute vectors in other point in local base(elambdaloc, ephiloc)
            localComponents.x = elambdaloc[0] * vxx + elambdaloc[1] * vyy + elambdaloc[2] * vzz;
            localComponents.y = ephiloc[0] * vxx + ephiloc[1] * vyy + ephiloc[2] * vzz;
        }
        else
        {
            // cartesian and spherical
            if (projection == Projection::cartesian || projection == Projection::spherical)
            {
                localComponents = globalComponents;
            }
        }
    }

    Point NormalVectorOutside(const Point& firstPoint, const Point& secondPoint, const Projection& projection)
    {
        if (!firstPoint.IsValid() || !secondPoint.IsValid())
        {
            return {constants::missing::doubleValue, constants::missing::doubleValue};
        }

        if (projection == Projection::sphericalAccurate)
        {
            const auto middlePoint = ComputeMiddlePoint(firstPoint, secondPoint, projection);

            const Cartesian3DPoint firstPointCartesianCoordinates{SphericalToCartesian3D(firstPoint)};
            const Cartesian3DPoint secondPointCartesianCoordinates{SphericalToCartesian3D(secondPoint)};

            // compute the base vectors at middle point
            std::array<double, 3> elambda{0.0, 0.0, 0.0};
            std::array<double, 3> ephi{0.0, 0.0, 0.0};
            ComputeTwoBaseComponents(middlePoint, elambda, ephi);

            // project vector in local base
            const double dx = (secondPointCartesianCoordinates.x - firstPointCartesianCoordinates.x) * elambda[0] +
                              (secondPointCartesianCoordinates.y - firstPointCartesianCoordinates.y) * elambda[1] +
                              (secondPointCartesianCoordinates.z - firstPointCartesianCoordinates.z) * elambda[2];

            const double dy = (secondPointCartesianCoordinates.x - firstPointCartesianCoordinates.x) * ephi[0] +
                              (secondPointCartesianCoordinates.y - firstPointCartesianCoordinates.y) * ephi[1] +
                              (secondPointCartesianCoordinates.z - firstPointCartesianCoordinates.z) * ephi[2];

            const double squaredDistance = dx * dx + dy * dy;
            Point result{constants::missing::doubleValue, constants::missing::doubleValue};
            if (squaredDistance > 0.0)
            {
                const double distance = sqrt(squaredDistance);
                result.x = dy / distance;
                result.y = -dx / distance;
            }
            return result;
        }

        // cartesian and spherical
        if (projection == Projection::cartesian || projection == Projection::spherical)
        {
            const double dx = GetDx(firstPoint, secondPoint, projection);
            const double dy = GetDy(firstPoint, secondPoint, projection);

            const double squaredDistance = dx * dx + dy * dy;
            Point result{constants::missing::doubleValue, constants::missing::doubleValue};
            if (squaredDistance > 0.0)
            {
                const double distance = sqrt(squaredDistance);
                result.x = dy / distance;
                result.y = -dx / distance;
            }

            if (projection == Projection::spherical)
            {
                result.x = result.x / cos(constants::conversion::degToRad * 0.5 * (firstPoint.y + secondPoint.y));
                result.y = result.y;
            }
            return result;
        }
        return {constants::missing::doubleValue, constants::missing::doubleValue};
    }

    void NormalVectorInside(const Point& firstPoint, const Point& secondPoint, const Point& insidePoint, Point& normal, bool& flippedNormal, const Projection& projection)
    {
        normal = NormalVectorOutside(firstPoint, secondPoint, projection);
        Point thirdPoint;
        if (projection == Projection::cartesian || projection == Projection::spherical)
        {
            flippedNormal = false;
            thirdPoint.x = firstPoint.x + normal.x;
            thirdPoint.y = firstPoint.y + normal.y;
        }
        if (projection == Projection::sphericalAccurate)
        {
            const auto middle = ComputeMiddlePoint(firstPoint, secondPoint, projection);
            Point localComponents;
            TransformGlobalVectorToLocal(firstPoint, middle, normal, projection, localComponents);

            std::array<double, 3> elambda{0.0, 0.0, 0.0};
            std::array<double, 3> ephi{0.0, 0.0, 0.0};
            ComputeTwoBaseComponents(firstPoint, elambda, ephi);

            const double vxx = localComponents.x * elambda[0] + localComponents.y * ephi[0];
            const double vyy = localComponents.x * elambda[1] + localComponents.y * ephi[1];
            const double vzz = localComponents.x * elambda[2] + localComponents.y * ephi[2];

            const Cartesian3DPoint firstPointCartesian{SphericalToCartesian3D(firstPoint)};

            Cartesian3DPoint rotatedPoint;
            const double alpha = 0.0;
            rotatedPoint.x = firstPointCartesian.x + alpha * vxx;
            rotatedPoint.y = firstPointCartesian.y + alpha * vyy;
            rotatedPoint.z = firstPointCartesian.z + alpha * vzz;

            thirdPoint = Cartesian3DToSpherical(rotatedPoint, firstPoint.x);
        }

        if (OuterProductTwoSegments(firstPoint, thirdPoint, firstPoint, secondPoint, projection) * OuterProductTwoSegments(firstPoint, insidePoint, firstPoint, secondPoint, projection) > 0.0)
        {
            normal.x = -normal.x;
            normal.y = -normal.y;
            flippedNormal = true;
        }
        else
        {
            flippedNormal = false;
        }
    }

    void AddIncrementToPoint(const Point& normal, double increment, const Point& referencePoint, const Projection& projection, Point& point)
    {
        if (projection == Projection::cartesian)
        {
            point.x = point.x + normal.x * increment;
            point.y = point.y + normal.y * increment;
        }
        if (projection == Projection::spherical || projection == Projection::sphericalAccurate)
        {
            const double xf = 1.0 / std::cos(constants::conversion::degToRad * referencePoint.y);
            const double convertedIncrement = constants::conversion::radToDeg * increment / constants::geometric::earth_radius;
            point.x = point.x + normal.x * convertedIncrement * xf;
            point.y = point.y + normal.y * convertedIncrement;
        }
    }

    void TranslateSphericalCoordinates(std::vector<Point>& polygon)
    {
        double minX = std::numeric_limits<double>::max();
        double maxX = std::numeric_limits<double>::lowest();

        for (UInt i = 0; i < polygon.size(); ++i)
        {
            minX = std::min(polygon[i].x, minX);
            maxX = std::max(polygon[i].x, maxX);
        }

        if (maxX - minX > 180.0)
        {
            const double deltaX = maxX - 180.0;

            for (UInt i = 0; i < polygon.size(); ++i)
            {
                if (polygon[i].x < deltaX)
                {
                    polygon[i].x = polygon[i].x + 360.0;
                }
            }
        }
    }

    Point ReferencePoint(const std::vector<Point>& polygon, const Projection& projection)
    {
        double minX = std::numeric_limits<double>::max();
        // Used only in spherical coordinate system, but quicker to compute at the same time as the minX
        double maxX = std::numeric_limits<double>::lowest();
        double minY = std::numeric_limits<double>::max();
        const auto numPoints = static_cast<UInt>(polygon.size());

        for (UInt i = 0; i < numPoints; ++i)
        {
            minX = std::min(polygon[i].x, minX);
            maxX = std::max(polygon[i].x, maxX);

            if (abs(polygon[i].y) < abs(minY))
            {
                minY = polygon[i].y;
            }
        }

        if (projection == Projection::spherical)
        {
            if (maxX - minX > 180.0)
            {
                minX += 360.0;
            }
        }

        return Point{minX, minY};
    }

    Point ReferencePoint(const std::vector<Point>& nodes,
                         const std::vector<UInt>& polygonIndices,
                         const Projection& projection)
    {
        double minX = std::numeric_limits<double>::max();
        // Used only in spherical coordinate system, but quicker to compute at the same time as the minX
        double maxX = std::numeric_limits<double>::lowest();
        double minY = std::numeric_limits<double>::max();
        const auto numPoints = static_cast<UInt>(polygonIndices.size());

        for (UInt i = 0; i < numPoints; ++i)
        {
            minX = std::min(nodes[polygonIndices[i]].x, minX);
            maxX = std::max(nodes[polygonIndices[i]].x, maxX);

            if (abs(nodes[polygonIndices[i]].y) < abs(minY))
            {
                minY = nodes[polygonIndices[i]].y;
            }
        }

        if (projection == Projection::spherical)
        {
            if (maxX - minX > 180.0)
            {
                minX += 360.0;
            }
        }

        return Point{minX, minY};
    }

    double ComputeSquaredDistance(const Point& firstPoint, const Point& secondPoint, const Projection& projection)
    {

        if (!firstPoint.IsValid() || !secondPoint.IsValid())
        {
            return 0.0;
        }

        if (projection == Projection::sphericalAccurate)
        {
            const Cartesian3DPoint firstPointCartesian{SphericalToCartesian3D(firstPoint)};
            const auto xx1 = firstPointCartesian.x;
            const auto yy1 = firstPointCartesian.y;
            const auto zz1 = firstPointCartesian.z;

            const Cartesian3DPoint secondPointCartesian{SphericalToCartesian3D(secondPoint)};
            const auto xx2 = secondPointCartesian.x;
            const auto yy2 = secondPointCartesian.y;
            const auto zz2 = secondPointCartesian.z;

            return (xx2 - xx1) * (xx2 - xx1) + (yy2 - yy1) * (yy2 - yy1) + (zz2 - zz1) * (zz2 - zz1);
        }

        // cartesian and spherical
        if (projection == Projection::cartesian || projection == Projection::spherical)
        {
            const double dx = GetDx(firstPoint, secondPoint, projection);
            const double dy = GetDy(firstPoint, secondPoint, projection);
            return dx * dx + dy * dy;
        }

        return constants::missing::doubleValue;
    }

    double ComputeDistance(const Point& firstPoint, const Point& secondPoint, const Projection& projection)
    {
        double distance = ComputeSquaredDistance(firstPoint, secondPoint, projection);
        if (distance >= 0.0)
        {
            distance = sqrt(distance);
        }
        return distance;
    }

    std::tuple<double, Point, double> DistanceFromLine(const Point& point, const Point& firstNode, const Point& secondNode, const Projection& projection)
    {
        double distance = constants::missing::doubleValue;
        Point normalPoint;
        double ratio = constants::missing::doubleValue;
        if (projection == Projection::cartesian || projection == Projection::spherical)
        {
            const auto squaredDistance = ComputeSquaredDistance(secondNode, firstNode, projection);
            if (squaredDistance != 0.0)
            {
                ratio = (GetDx(firstNode, point, projection) * GetDx(firstNode, secondNode, projection) +
                         GetDy(firstNode, point, projection) * GetDy(firstNode, secondNode, projection)) /
                        squaredDistance;
                const auto correctedRatio = std::max(std::min(1.0, ratio), 0.0);

                normalPoint = firstNode + correctedRatio * (secondNode - firstNode);
                distance = ComputeDistance(point, normalPoint, projection);
            }
        }

        if (projection == Projection::sphericalAccurate)
        {
            const Cartesian3DPoint firstNodeCartesian{SphericalToCartesian3D(firstNode)};
            const auto xx1 = firstNodeCartesian.x;
            const auto yy1 = firstNodeCartesian.y;
            const auto zz1 = firstNodeCartesian.z;

            const Cartesian3DPoint secondNodeCartesian{SphericalToCartesian3D(secondNode)};
            const auto xx2 = secondNodeCartesian.x;
            const auto yy2 = secondNodeCartesian.y;
            const auto zz2 = secondNodeCartesian.z;

            const Cartesian3DPoint pointCartesian{SphericalToCartesian3D(point)};
            const auto xx3 = pointCartesian.x;
            const auto yy3 = pointCartesian.y;
            const auto zz3 = pointCartesian.z;

            const double x21 = xx2 - xx1;
            const double y21 = yy2 - yy1;
            const double z21 = zz2 - zz1;
            const double x31 = xx3 - xx1;
            const double y31 = yy3 - yy1;
            const double z31 = zz3 - zz1;

            const double r2 = x21 * x21 + y21 * y21 + z21 * z21;

            ratio = 0.0;
            if (r2 >= 0.0)
            {

                ratio = (x31 * x21 + y31 * y21 + z31 * z21) / r2;
                const double correctedRatio = std::max(std::min(1.0, ratio), 0.0);

                Cartesian3DPoint cartesianNormal3DPoint;
                cartesianNormal3DPoint.x = firstNodeCartesian.x + correctedRatio * x21;
                cartesianNormal3DPoint.y = firstNodeCartesian.y + correctedRatio * y21;
                cartesianNormal3DPoint.z = firstNodeCartesian.z + correctedRatio * z21;

                cartesianNormal3DPoint.x = cartesianNormal3DPoint.x - xx3;
                cartesianNormal3DPoint.y = cartesianNormal3DPoint.y - yy3;
                cartesianNormal3DPoint.z = cartesianNormal3DPoint.z - zz3;

                distance = std::sqrt(cartesianNormal3DPoint.x * cartesianNormal3DPoint.x +
                                     cartesianNormal3DPoint.y * cartesianNormal3DPoint.y +
                                     cartesianNormal3DPoint.z * cartesianNormal3DPoint.z);

                const double referenceLongitude = std::max({firstNode.x, secondNode.x, point.x});
                normalPoint = Cartesian3DToSpherical(cartesianNormal3DPoint, referenceLongitude);
            }
        }

        return {distance, normalPoint, ratio};
    }

    double InnerProductTwoSegments(const Point& firstPointFirstSegment, const Point& secondPointFirstSegment, const Point& firstPointSecondSegment, const Point& secondPointSecondSegment, const Projection& projection)
    {
        if (projection == Projection::sphericalAccurate)
        {
            const Cartesian3DPoint firstPointFirstSegment3D{SphericalToCartesian3D(firstPointFirstSegment)};
            const Cartesian3DPoint secondPointFirstSegment3D{SphericalToCartesian3D(secondPointFirstSegment)};
            const Cartesian3DPoint firstPointSecondSegment3D{SphericalToCartesian3D(firstPointSecondSegment)};
            const Cartesian3DPoint secondPointSecondSegment3D{SphericalToCartesian3D(secondPointSecondSegment)};

            const Cartesian3DPoint delta1 = secondPointFirstSegment3D - firstPointFirstSegment3D;
            const Cartesian3DPoint delta2 = secondPointSecondSegment3D - firstPointSecondSegment3D;

            return InnerProduct(delta1, delta2);
        }

        // cartesian and spherical
        if (projection == Projection::cartesian || projection == Projection::spherical)
        {
            const Vector delta1{GetDx(firstPointFirstSegment, secondPointFirstSegment, projection), GetDy(firstPointFirstSegment, secondPointFirstSegment, projection)};
            const Vector delta2{GetDx(firstPointSecondSegment, secondPointSecondSegment, projection), GetDy(firstPointSecondSegment, secondPointSecondSegment, projection)};

            return dot(delta1, delta2);
        }

        return constants::missing::doubleValue;
    }

    double NormalizedInnerProductTwoSegments(const Point& firstPointFirstSegment, const Point& secondPointFirstSegment, const Point& firstPointSecondSegment, const Point& secondPointSecondSegment, const Projection& projection)
    {
        if (projection == Projection::sphericalAccurate)
        {
            const Cartesian3DPoint firstPointFirstSegmentCartesian{SphericalToCartesian3D(firstPointFirstSegment)};
            const auto xx1 = firstPointFirstSegmentCartesian.x;
            const auto yy1 = firstPointFirstSegmentCartesian.y;
            const auto zz1 = firstPointFirstSegmentCartesian.z;

            const Cartesian3DPoint secondPointFirstSegmentCartesian{SphericalToCartesian3D(secondPointFirstSegment)};
            const auto xx2 = secondPointFirstSegmentCartesian.x;
            const auto yy2 = secondPointFirstSegmentCartesian.y;
            const auto zz2 = secondPointFirstSegmentCartesian.z;

            const Cartesian3DPoint firstPointSecondSegmentCartesian{SphericalToCartesian3D(firstPointSecondSegment)};
            const auto xx3 = firstPointSecondSegmentCartesian.x;
            const auto yy3 = firstPointSecondSegmentCartesian.y;
            const auto zz3 = firstPointSecondSegmentCartesian.z;

            const Cartesian3DPoint secondPointSecondSegmentCartesian{SphericalToCartesian3D(secondPointSecondSegment)};
            const auto xx4 = secondPointSecondSegmentCartesian.x;
            const auto yy4 = secondPointSecondSegmentCartesian.y;
            const auto zz4 = secondPointSecondSegmentCartesian.z;

            const auto dx1 = xx2 - xx1;
            const auto dy1 = yy2 - yy1;
            const auto dz1 = zz2 - zz1;
            const auto firstSegmentDistance = dx1 * dx1 + dy1 * dy1 + dz1 * dz1;

            const auto dx2 = xx4 - xx3;
            const auto dy2 = yy4 - yy3;
            const auto dz2 = zz4 - zz3;
            const auto secondSegmentDistance = dx2 * dx2 + dy2 * dy2 + dz2 * dz2;

            double cosphi;
            if (firstSegmentDistance <= 0.0 || secondSegmentDistance <= 0.0)
            {
                cosphi = constants::missing::doubleValue;
            }
            else
            {
                cosphi = (dx1 * dx2 + dy1 * dy2 + dz1 * dz2) / sqrt(firstSegmentDistance * secondSegmentDistance);
            }
            return cosphi;
        }

        // cartesian and spherical
        if (projection == Projection::cartesian || projection == Projection::spherical)
        {
            const auto dx1 = GetDx(firstPointFirstSegment, secondPointFirstSegment, projection);
            const auto dx2 = GetDx(firstPointSecondSegment, secondPointSecondSegment, projection);

            const auto dy1 = GetDy(firstPointFirstSegment, secondPointFirstSegment, projection);
            const auto dy2 = GetDy(firstPointSecondSegment, secondPointSecondSegment, projection);

            const auto r1 = dx1 * dx1 + dy1 * dy1;
            const auto r2 = dx2 * dx2 + dy2 * dy2;

            double cosphi;
            if (r1 <= 0.0 || r2 <= 0.0)
            {
                return constants::missing::doubleValue;
            }

            cosphi = (dx1 * dx2 + dy1 * dy2) / std::sqrt(r1 * r2);
            cosphi = std::max(std::min(cosphi, 1.0), -1.0);
            return cosphi;
        }

        return constants::missing::doubleValue;
    }

    Point CircumcenterOfTriangle(const Point& firstNode, const Point& secondNode, const Point& thirdNode, const Projection& projection)
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
        if (projection == Projection::spherical)
        {
            const double phi = (firstNode.y + secondNode.y + thirdNode.y) * constants::numeric::oneThird;
            const double xf = 1.0 / cos(constants::conversion::degToRad * phi);
            circumcenter.x = firstNode.x + xf * 0.5 * (dx3 - z * dy3) * constants::conversion::radToDeg / constants::geometric::earth_radius;
            circumcenter.y = firstNode.y + 0.5 * (dy3 + z * dx3) * constants::conversion::radToDeg / constants::geometric::earth_radius;
        }
        if (projection == Projection::sphericalAccurate)
        {
            // TODO: compute in case of spherical accurate (comp_circumcenter3D)
        }
        return circumcenter;
    }

    UInt CountNumberOfValidEdges(const std::vector<UInt>& edgesNumFaces, UInt numEdges)
    {
        if (numEdges > edgesNumFaces.size())
        {
            throw ConstraintError("Invalid range for array: {} > {}", numEdges, edgesNumFaces.size());
        }

        return static_cast<UInt>(std::count(edgesNumFaces.begin(), edgesNumFaces.begin() + numEdges, 2));
    }

    void ComputeMidPointsAndNormals(const std::vector<Point>& polygon,
                                    const std::vector<UInt>& edgesNumFaces,
                                    const UInt numNodes,
                                    std::array<Point, constants::geometric::maximumNumberOfNodesPerFace>& middlePoints,
                                    std::array<Point, constants::geometric::maximumNumberOfNodesPerFace>& normals,
                                    UInt& pointCount,
                                    const Projection& projection)
    {
        for (UInt n = 0; n < numNodes; n++)
        {
            if (edgesNumFaces[n] != 2)
            {
                continue;
            }

            const auto nextNode = NextCircularForwardIndex(n, numNodes);

            middlePoints[pointCount] = ((polygon[n] + polygon[nextNode]) * 0.5);
            normals[pointCount] = NormalVector(polygon[n], polygon[nextNode], middlePoints[pointCount], projection);
            ++pointCount;
        }
    }

    Point ComputeCircumCenter(const Point& centerOfMass,
                              const UInt pointCount,
                              const std::array<Point, constants::geometric::maximumNumberOfNodesPerFace>& middlePoints,
                              const std::array<Point, constants::geometric::maximumNumberOfNodesPerFace>& normals,
                              const Projection& projection)
    {
        const UInt maximumNumberCircumcenterIterations = 100;
        const double eps = projection == Projection::cartesian ? 1e-3 : 9e-10; // 111km = 0-e digit.

        Point estimatedCircumCenter = centerOfMass;

        for (UInt iter = 0; iter < maximumNumberCircumcenterIterations; ++iter)
        {
            const Point previousCircumCenter = estimatedCircumCenter;
            for (UInt n = 0; n < pointCount; n++)
            {
                const Point delta{GetDx(middlePoints[n], estimatedCircumCenter, projection), GetDy(middlePoints[n], estimatedCircumCenter, projection)};
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

    Point ComputeFaceCircumenter(std::vector<Point>& polygon,
                                 const std::vector<UInt>& edgesNumFaces,
                                 const Projection& projection)
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
            result = CircumcenterOfTriangle(polygon[0], polygon[1], polygon[2], projection);
        }
        else if (!edgesNumFaces.empty())
        {
            UInt numValidEdges = CountNumberOfValidEdges(edgesNumFaces, numNodes);

            if (numValidEdges > 1)
            {
                ComputeMidPointsAndNormals(polygon, edgesNumFaces, numNodes, middlePoints, normals, pointCount, projection);
                result = ComputeCircumCenter(centerOfMass, pointCount, middlePoints, normals, projection);
            }
        }

        for (UInt n = 0; n < numNodes; ++n)
        {
            polygon[n] = weightCircumCenter * polygon[n] + (1.0 - weightCircumCenter) * centerOfMass;
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

    std::tuple<bool, Point, double, double, double> AreSegmentsCrossing(const Point& firstSegmentFirstPoint,
                                                                        const Point& firstSegmentSecondPoint,
                                                                        const Point& secondSegmentFirstPoint,
                                                                        const Point& secondSegmentSecondPoint,
                                                                        bool adimensionalCrossProduct,
                                                                        const Projection& projection)
    {
        bool isCrossing = false;
        Point intersectionPoint;
        double ratioFirstSegment = constants::missing::doubleValue;
        double ratioSecondSegment = constants::missing::doubleValue;
        double crossProduct = constants::missing::doubleValue;

        if (projection == Projection::cartesian || projection == Projection::spherical)
        {

            auto const x21 = GetDx(firstSegmentFirstPoint, firstSegmentSecondPoint, projection);
            auto const y21 = GetDy(firstSegmentFirstPoint, firstSegmentSecondPoint, projection);

            auto const x43 = GetDx(secondSegmentFirstPoint, secondSegmentSecondPoint, projection);
            auto const y43 = GetDy(secondSegmentFirstPoint, secondSegmentSecondPoint, projection);

            auto const x31 = GetDx(firstSegmentFirstPoint, secondSegmentFirstPoint, projection);
            auto const y31 = GetDy(firstSegmentFirstPoint, secondSegmentFirstPoint, projection);

            auto const det = x43 * y21 - y43 * x21;

            double maxValue = std::max(std::max(std::abs(x21), std::abs(y21)),
                                       std::max(std::abs(x43), std::abs(y43)));
            const double eps = std::max(0.00001 * maxValue, std::numeric_limits<double>::denorm_min());

            if (std::abs(det) < eps)
            {
                return {isCrossing, intersectionPoint, crossProduct, ratioFirstSegment, ratioSecondSegment};
            }

            ratioSecondSegment = (y31 * x21 - x31 * y21) / det;
            ratioFirstSegment = (y31 * x43 - x31 * y43) / det;
            if (ratioFirstSegment >= 0.0 && ratioFirstSegment <= 1.0 && ratioSecondSegment >= 0.0 && ratioSecondSegment <= 1.0)
            {
                isCrossing = true;
            }
            intersectionPoint.x = firstSegmentFirstPoint.x + ratioFirstSegment * (firstSegmentSecondPoint.x - firstSegmentFirstPoint.x);
            intersectionPoint.y = firstSegmentFirstPoint.y + ratioFirstSegment * (firstSegmentSecondPoint.y - firstSegmentFirstPoint.y);
            crossProduct = -det;
            if (adimensionalCrossProduct)
            {
                crossProduct = -det / (std::sqrt(x21 * x21 + y21 * y21) * std::sqrt(x43 * x43 + y43 * y43) + 1e-8);
            }
        }

        if (projection == Projection::sphericalAccurate)
        {
            const Cartesian3DPoint firstSegmentFirstCartesian3DPoint{SphericalToCartesian3D(firstSegmentFirstPoint)};

            const Cartesian3DPoint firstSegmentSecondCartesian3DPoint{SphericalToCartesian3D(firstSegmentSecondPoint)};

            const Cartesian3DPoint secondSegmentFirstCartesian3DPoint{SphericalToCartesian3D(secondSegmentFirstPoint)};

            const Cartesian3DPoint secondSegmentSecondCartesian3DPoint{SphericalToCartesian3D(secondSegmentSecondPoint)};

            auto n12 = VectorProduct(firstSegmentFirstCartesian3DPoint, firstSegmentSecondCartesian3DPoint);
            const auto n12InnerProduct = std::sqrt(InnerProduct(n12, n12));
            n12.x = n12.x / n12InnerProduct;
            n12.y = n12.y / n12InnerProduct;
            n12.z = n12.z / n12InnerProduct;

            auto n34 = VectorProduct(secondSegmentFirstCartesian3DPoint, secondSegmentSecondCartesian3DPoint);
            const auto n34InnerProduct = std::sqrt(InnerProduct(n34, n34));
            n34.x = n34.x / n34InnerProduct;
            n34.y = n34.y / n34InnerProduct;
            n34.z = n34.z / n34InnerProduct;

            const auto n12n34InnerProduct = std::sqrt(std::abs(InnerProduct(n12, n34)));

            const double tolerance = 1e-12;
            if (n12n34InnerProduct > tolerance)
            {
                Cartesian3DPoint firstSegmentDifference;
                firstSegmentDifference.x = firstSegmentSecondCartesian3DPoint.x - firstSegmentFirstCartesian3DPoint.x;
                firstSegmentDifference.y = firstSegmentSecondCartesian3DPoint.y - firstSegmentFirstCartesian3DPoint.y;
                firstSegmentDifference.z = firstSegmentSecondCartesian3DPoint.z - firstSegmentFirstCartesian3DPoint.z;

                Cartesian3DPoint secondSegmentDifference;
                secondSegmentDifference.x = secondSegmentSecondCartesian3DPoint.x - secondSegmentFirstCartesian3DPoint.x;
                secondSegmentDifference.y = secondSegmentSecondCartesian3DPoint.y - secondSegmentFirstCartesian3DPoint.y;
                secondSegmentDifference.z = secondSegmentSecondCartesian3DPoint.z - secondSegmentFirstCartesian3DPoint.z;

                const auto Det12 = InnerProduct(firstSegmentDifference, n34);
                const auto Det34 = InnerProduct(secondSegmentDifference, n12);

                if (std::abs(Det12) > tolerance && std::abs(Det34) > tolerance)
                {
                    ratioFirstSegment = -InnerProduct(firstSegmentFirstCartesian3DPoint, n34) / Det12;
                    ratioSecondSegment = -InnerProduct(secondSegmentFirstCartesian3DPoint, n12) / Det34;
                }
            }

            if (ratioSecondSegment >= 0.0 && ratioSecondSegment <= 1.0 &&
                ratioFirstSegment >= 0.0 && ratioFirstSegment <= 1.0)
            {
                // check if segments are crossing
                isCrossing = true;

                // compute intersection
                Cartesian3DPoint intersectionCartesian3DPoint;
                intersectionCartesian3DPoint.x = firstSegmentFirstCartesian3DPoint.x + ratioFirstSegment * (firstSegmentSecondCartesian3DPoint.x - firstSegmentFirstCartesian3DPoint.x);
                intersectionCartesian3DPoint.y = firstSegmentFirstCartesian3DPoint.y + ratioFirstSegment * (firstSegmentSecondCartesian3DPoint.y - firstSegmentFirstCartesian3DPoint.y);
                intersectionCartesian3DPoint.z = firstSegmentFirstCartesian3DPoint.z + ratioFirstSegment * (firstSegmentSecondCartesian3DPoint.z - firstSegmentFirstCartesian3DPoint.z);
                intersectionPoint = Cartesian3DToSpherical(intersectionCartesian3DPoint, std::max(firstSegmentFirstPoint.x, firstSegmentSecondPoint.x));
            }
        }

        return {isCrossing, intersectionPoint, crossProduct, ratioFirstSegment, ratioSecondSegment};
    }

    std::tuple<std::vector<double>, double> ComputeAdimensionalDistancesFromPointSerie(const std::vector<Point>& v, const Projection& projection)
    {
        std::vector<double> result(v.size());

        result[0] = 0;
        for (UInt i = 1; i < v.size(); ++i)
        {
            result[i] = result[i - 1] + ComputeDistance(v[i - 1], v[i], projection);
        }
        auto totalDistance = result.back();
        if (IsEqual(totalDistance, 0.0))
        {
            return {result, totalDistance};
        }
        const double inverseTotalDistance = 1.0 / totalDistance;
        for (UInt i = 1; i < v.size(); ++i)
        {
            result[i] = result[i] * inverseTotalDistance;
        }

        return {result, totalDistance};
    }

    lin_alg::Matrix<Point> DiscretizeTransfinite(const std::vector<Point>& leftDiscretization,
                                                 const std::vector<Point>& rightDiscretization,
                                                 const std::vector<Point>& bottomDiscretization,
                                                 const std::vector<Point>& upperDiscretization,
                                                 const Projection& projection,
                                                 UInt numM,
                                                 UInt numN)
    {
        const auto [leftAdimensional, totalLengthLeft] = ComputeAdimensionalDistancesFromPointSerie(leftDiscretization, projection);
        const auto [rightAdimensional, totalLengthRight] = ComputeAdimensionalDistancesFromPointSerie(rightDiscretization, projection);
        const auto [bottomAdimensional, totalLengthBottom] = ComputeAdimensionalDistancesFromPointSerie(bottomDiscretization, projection);
        const auto [upperAdimensional, totalLengthUpper] = ComputeAdimensionalDistancesFromPointSerie(upperDiscretization, projection);

        // now compute the adimensional distance of each point to be filled
        const auto numNPoints = numN + 1;
        const auto numMPoints = numM + 1;

        lin_alg::Matrix<double> iWeightFactor(numNPoints, numMPoints);
        lin_alg::Matrix<double> jWeightFactor(numNPoints, numMPoints);

        for (UInt n = 0; n < numNPoints; n++)
        {
            for (UInt m = 0; m < numMPoints; m++)
            {
                const double nWeight = double(n) / double(numN);
                const double mWeight = double(m) / double(numM);

                const auto ifactor = (1.0 - mWeight) * bottomAdimensional[n] + mWeight * upperAdimensional[n];
                const auto jfactor = (1.0 - nWeight) * leftAdimensional[m] + nWeight * rightAdimensional[m];

                iWeightFactor(n, m) = ifactor;
                jWeightFactor(n, m) = jfactor;
            }
        }

        lin_alg::Matrix<double> ones = lin_alg::Matrix<double>::Ones(numNPoints, numMPoints);
        lin_alg::Matrix<double> weightOne = (ones - jWeightFactor) * totalLengthBottom + jWeightFactor * totalLengthUpper;
        lin_alg::Matrix<double> weightTwo = (ones - iWeightFactor) * totalLengthLeft + iWeightFactor * totalLengthRight;
        lin_alg::Matrix<double> weightThree = weightTwo.cwiseQuotient(weightOne);
        lin_alg::Matrix<double> weightFour = weightThree.cwiseInverse();
        lin_alg::Matrix<double> const wa = (weightThree + weightFour).cwiseInverse();
        weightOne = wa.cwiseProduct(weightThree);
        weightTwo = wa.cwiseProduct(weightFour);

        // border points
        lin_alg::Matrix<Point> result(numNPoints, numMPoints);
        for (UInt n = 0; n < numNPoints; n++)
        {
            result(n, 0) = bottomDiscretization[n];
            result(n, numM) = upperDiscretization[n];
        }
        for (UInt m = 0; m < numMPoints; m++)
        {
            result(0, m) = leftDiscretization[m];
            result(numN, m) = rightDiscretization[m];
        }

        // first interpolation
        for (UInt n = 1; n < numN; n++)
        {
            for (UInt m = 1; m < numM; m++)
            {
                const auto x_coord = (leftDiscretization[m].x * (1.0 - iWeightFactor(n, m)) +
                                      rightDiscretization[m].x * iWeightFactor(n, m)) *
                                         weightOne(n, m) +
                                     (bottomDiscretization[n].x * (1.0 - jWeightFactor(n, m)) +
                                      upperDiscretization[n].x * jWeightFactor(n, m)) *
                                         weightTwo(n, m);

                result(n, m).x = x_coord;

                const auto y_coord = (leftDiscretization[m].y * (1.0 - iWeightFactor(n, m)) +
                                      rightDiscretization[m].y * iWeightFactor(n, m)) *
                                         weightOne(n, m) +
                                     (bottomDiscretization[n].y * (1.0 - jWeightFactor(n, m)) +
                                      upperDiscretization[n].y * jWeightFactor(n, m)) *
                                         weightTwo(n, m);

                result(n, m).y = y_coord;
            }
        }

        // update weights
        for (UInt n = 0; n < numNPoints; n++)
        {
            for (UInt m = 0; m < numMPoints; m++)
            {
                const auto valOne = (1.0 - jWeightFactor(n, m)) * bottomAdimensional[n] * totalLengthBottom +
                                    jWeightFactor(n, m) * upperAdimensional[n] * totalLengthUpper;

                weightOne(n, m) = valOne;

                const auto valTwo = (1.0 - iWeightFactor(n, m)) * leftAdimensional[m] * totalLengthLeft +
                                    iWeightFactor(n, m) * rightAdimensional[m] * totalLengthRight;

                weightTwo(n, m) = valTwo;
            }
        }

        for (UInt n = 1; n < numNPoints; n++)
        {
            for (UInt m = 0; m < numMPoints; m++)
            {
                const auto val = weightOne(n, m) - weightOne(n - 1, m);
                weightThree(n, m) = val;
            }
        }

        for (UInt n = 0; n < numNPoints; n++)
        {
            for (UInt m = 1; m < numMPoints; m++)
            {

                const auto val = weightTwo(n, m) - weightTwo(n, m - 1);
                weightFour(n, m) = val;
            }
        }

        for (UInt n = 1; n < numNPoints; n++)
        {
            for (UInt m = 1; m < numMPoints - 1; m++)
            {
                const auto val = 0.25 *
                                 (weightFour(n, m) +
                                  weightFour(n, m + 1) +
                                  weightFour(n - 1, m) +
                                  weightFour(n - 1, m + 1)) /
                                 weightThree(n, m);

                weightOne(n, m) = val;
            }
        }

        for (UInt n = 1; n < numNPoints - 1; n++)
        {
            for (UInt m = 1; m < numMPoints; m++)
            {
                const auto val = 0.25 *
                                 (weightThree(n, m) +
                                  weightThree(n, m - 1) +
                                  weightThree(n + 1, m) +
                                  weightThree(n + 1, m - 1)) /
                                 weightFour(n, m);

                weightTwo(n, m) = val;
            }
        }

        // Iterate several times over
        // Move to constants
        static UInt constexpr numIterations = 25;

        for (UInt iter = 0; iter < numIterations; iter++)
        {
            // re-assign the weights
            for (UInt n = 0; n < numNPoints; n++)
            {
                for (UInt m = 0; m < numMPoints; m++)
                {
                    weightThree(n, m) = result(n, m).x;
                    weightFour(n, m) = result(n, m).y;
                }
            }

            for (UInt n = 1; n < numN; n++)
            {
                for (UInt m = 1; m < numM; m++)
                {

                    const double wa = 1.0 /
                                      (weightOne(n, m) +
                                       weightOne(n + 1, m) +
                                       weightTwo(n, m) +
                                       weightTwo(n, m + 1));

                    const auto x_coord = wa *
                                         (weightThree(n - 1, m) * weightOne(n, m) +
                                          weightThree(n + 1, m) * weightOne(n + 1, m) +
                                          weightThree(n, m - 1) * weightTwo(n, m) +
                                          weightThree(n, m + 1) * weightTwo(n, m + 1));

                    const auto y_coord = wa * (weightFour(n - 1, m) * weightOne(n, m) +
                                               weightFour(n + 1, m) * weightOne(n + 1, m) +
                                               weightFour(n, m - 1) * weightTwo(n, m) +
                                               weightFour(n, m + 1) * weightTwo(n, m + 1));

                    result(n, m).x = x_coord;

                    result(n, m).y = y_coord;
                }
            }
        }

        return result;
    }

    std::vector<Point> ComputeEdgeCenters(const std::vector<Point>& nodes, const std::vector<Edge>& edges)
    {
        std::vector<Point> edgesCenters(edges.size());

        for (UInt i = 0; i < edges.size(); ++i)
        {
            auto const first = edges[i].first;
            auto const second = edges[i].second;

            if (first == constants::missing::uintValue || second == constants::missing::uintValue)
            {
                // Initialise with invalid data.
                edgesCenters[i] = Point(constants::missing::doubleValue, constants::missing::doubleValue);
            }
            else
            {
                edgesCenters[i] = (nodes[first] + nodes[second]) * 0.5;
            }
        }
        return edgesCenters;
    }

    double LinearInterpolationInTriangle(const Point& interpolationPoint, const std::vector<Point>& polygon, const std::vector<double>& values, const Projection& projection)
    {
        double result = constants::missing::doubleValue;

        const auto a11 = GetDx(polygon[0], polygon[1], projection);
        const auto a21 = GetDy(polygon[0], polygon[1], projection);

        const auto a12 = GetDx(polygon[0], polygon[2], projection);
        const auto a22 = GetDy(polygon[0], polygon[2], projection);

        const auto b1 = GetDx(polygon[0], interpolationPoint, projection);
        const auto b2 = GetDy(polygon[0], interpolationPoint, projection);

        const double det = a11 * a22 - a12 * a21;
        if (std::abs(det) < 1e-12)
        {
            return result;
        }

        const double rlam = (a22 * b1 - a12 * b2) / det;
        const double rmhu = (-a21 * b1 + a11 * b2) / det;

        result = values[0] + rlam * (values[1] - values[0]) + rmhu * (values[2] - values[0]);

        return result;
    }

    meshkernel::Point ComputeAverageCoordinate(const std::span<const Point> points, const Projection& projection)
    {
        size_t validCount = std::ranges::count_if(points, [](const Point& p)
                                                  { return p.IsValid(); });

        if (projection == Projection::sphericalAccurate)
        {

            UInt firstValidPoint = 0;

            if (validCount != points.size())
            {
                auto iterator = std::ranges::find_if(points, [](const Point& p)
                                                     { return p.IsValid(); });
                firstValidPoint = static_cast<UInt>(iterator - points.begin());
            }

            Cartesian3DPoint averagePoint3D{0.0, 0.0, 0.0};
            for (const auto& point : points)
            {
                if (!point.IsValid())
                {
                    continue;
                }

                const Cartesian3DPoint point3D{SphericalToCartesian3D(point)};
                averagePoint3D.x += point3D.x;
                averagePoint3D.y += point3D.y;
                averagePoint3D.z += point3D.z;
            }
            averagePoint3D.x = averagePoint3D.x / static_cast<double>(validCount);
            averagePoint3D.y = averagePoint3D.y / static_cast<double>(validCount);
            averagePoint3D.z = averagePoint3D.z / static_cast<double>(validCount);

            return Cartesian3DToSpherical(averagePoint3D, points[firstValidPoint].x);
        }

        auto result = std::accumulate(points.begin(), points.end(), Point{0.0, 0.0}, [](const Point& sum, const Point& current)
                                      { return current.IsValid() ? sum + current : sum; });
        result.x = result.x / static_cast<double>(validCount);
        result.y = result.y / static_cast<double>(validCount);
        return result;
    }

    std::tuple<Point, double, bool> OrthogonalProjectionOnSegment(Point const& firstNode,
                                                                  Point const& secondNode,
                                                                  Point const& pointToProject)
    {
        const auto delta = secondNode - firstNode;
        const auto squaredDelta = delta.x * delta.x + delta.y * delta.y;
        auto segmentRatio = (pointToProject.x * delta.x + pointToProject.y * delta.y - firstNode.x * delta.x - firstNode.y * delta.y) / squaredDelta;

        Point projectedPoint = firstNode + delta * segmentRatio;
        auto projectionOnSegment = segmentRatio >= 0.0 && segmentRatio <= 1.0 ? true : false;
        segmentRatio = segmentRatio * std::sqrt(squaredDelta);
        return {projectedPoint, segmentRatio, projectionOnSegment};
    }

    std::vector<double> ComputePolyLineEdgesLengths(const std::vector<Point>& polyline, Projection projection)
    {
        std::vector<double> edgeLengths;
        if (polyline.empty())
        {
            return edgeLengths;
        }

        edgeLengths.reserve(polyline.size());

        for (UInt p = 0; p < polyline.size() - 1; ++p)
        {
            const auto firstNode = p;
            auto secondNode = p + 1;
            edgeLengths.emplace_back(ComputeDistance(polyline[firstNode], polyline[secondNode], projection));
        }
        return edgeLengths;
    }

    std::vector<double> ComputePolyLineNodalChainages(std::vector<Point> const& polyline, Projection projection)
    {
        // Compute the edge lengths and the edge coordinates
        auto const edgeLengths = ComputePolyLineEdgesLengths(polyline, projection);
        if (edgeLengths.empty())
        {
            return edgeLengths;
        }

        std::vector<double> chainages(polyline.size());
        chainages[0] = 0.0;
        for (UInt i = 0; i < edgeLengths.size(); ++i)
        {
            chainages[i + 1] = chainages[i] + edgeLengths[i];
        }
        return chainages;
    }

    std::vector<Point> ComputePolyLineDiscretization(std::vector<Point> const& polyline, std::vector<double>& chainages, Projection projection)
    {
        if (polyline.size() < 2)
        {
            throw std::invalid_argument("ComputePolyLineDiscretization polyline with less than 2 points");
        }

        // Compute the edge lengths and the edge coordinates
        auto const edgeLengths = ComputePolyLineEdgesLengths(polyline, projection);
        std::vector<double> const polylineNodalCoordinate = ComputePolyLineNodalChainages(polyline, projection);

        std::vector<Point> discretization;
        discretization.reserve(chainages.size());
        UInt curentNodalIndex = 0;
        std::sort(chainages.begin(), chainages.end());
        for (auto const& chainage : chainages)
        {
            if (chainage > polylineNodalCoordinate[curentNodalIndex + 1])
            {
                curentNodalIndex++;
            }

            double distanceFromLastNode = chainage - polylineNodalCoordinate[curentNodalIndex];
            Point p = polyline[curentNodalIndex] + (polyline[curentNodalIndex + 1] - polyline[curentNodalIndex]) * distanceFromLastNode / edgeLengths[curentNodalIndex];

            discretization.emplace_back(p);
        }

        return discretization;
    }

    double MatrixNorm(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& matCoefficients)
    {
        return (matCoefficients[0] * x[0] + matCoefficients[1] * x[1]) * y[0] + (matCoefficients[2] * x[0] + matCoefficients[3] * x[1]) * y[1];
    }

} // namespace meshkernel
