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
#include "MeshKernel/Mesh.hpp"

namespace meshkernel
{
    Cartesian3DPoint VectorProduct(const Cartesian3DPoint& a, const Cartesian3DPoint& b)
    {
        return Cartesian3DPoint{
            a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x};
    }

    double InnerProduct(const Cartesian3DPoint& a, const Cartesian3DPoint& b)
    {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

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
        UInt startRange;
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
        UInt index = currentIndex + 1;
        if (index >= size)
        {
            index = index - size;
        }
        return index;
    }

    UInt NextCircularBackwardIndex(UInt currentIndex, UInt size)
    {
        if (currentIndex == 0)
        {
            return currentIndex + size - 1;
        }
        return currentIndex - 1;
    }

    bool IsPointOnPole(const Point& point)
    {
        return std::abs(std::abs(point.y) - 90.0) < constants::geometric::absLatitudeAtPoles;
    }

    Cartesian3DPoint SphericalToCartesian3D(const Point& sphericalPoint)
    {
        Cartesian3DPoint result;
        result.z = constants::geometric::earth_radius * sin(sphericalPoint.y * constants::conversion::degToRad);
        const double rr = constants::geometric::earth_radius * cos(sphericalPoint.y * constants::conversion::degToRad);
        result.x = rr * cos(sphericalPoint.x * constants::conversion::degToRad);
        result.y = rr * sin(sphericalPoint.x * constants::conversion::degToRad);
        return result;
    }

    Point Cartesian3DToSpherical(const Cartesian3DPoint& cartesianPoint, double referenceLongitude)
    {
        Point sphericalPoint;
        const double angle = atan2(cartesianPoint.y, cartesianPoint.x) * constants::conversion::radToDeg;
        sphericalPoint.y = atan2(cartesianPoint.z, sqrt(cartesianPoint.x * cartesianPoint.x + cartesianPoint.y * cartesianPoint.y)) * constants::conversion::radToDeg;
        sphericalPoint.x = angle + std::lround((referenceLongitude - angle) / 360.0) * 360.0;
        return sphericalPoint;
    }

    double crossProduct(const Point& firstSegmentFirstPoint, const Point& firstSegmentSecondPoint, const Point& secondSegmentFirstPoint, const Point& secondSegmentSecondPoint, const Projection::Type& projection)
    {
        const auto dx1 = GetDx(firstSegmentFirstPoint, firstSegmentSecondPoint, projection);
        const auto dy1 = GetDy(firstSegmentFirstPoint, firstSegmentSecondPoint, projection);
        const auto dx2 = GetDx(secondSegmentFirstPoint, secondSegmentSecondPoint, projection);
        const auto dy2 = GetDy(secondSegmentFirstPoint, secondSegmentSecondPoint, projection);
        return dx1 * dy2 - dy1 * dx2;
    }

    bool IsPointInPolygonNodes(const Point& point,
                               const std::vector<Point>& polygonNodes,
                               const Projection::Type& projection,
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
        if (currentPolygonSize < Mesh::m_numNodesInTriangle || polygonNodes.size() < currentPolygonSize)
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

        bool isInPolygon = false;

        if (projection == Projection::Type::Cartesian || projection == Projection::Type::Spherical)
        {

            int windingNumber = 0;
            for (auto n = startNode; n < endNode; n++)
            {

                const auto crossProductValue = crossProduct(polygonNodes[n], polygonNodes[n + 1], polygonNodes[n], point, Projection::Type::Cartesian);

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

            isInPolygon = windingNumber == 0 ? false : true;
        }

        if (projection == Projection::Type::SphericalAccurate)
        {
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
                cartesian3DPoints[i].x = polygonCenterCartesian3D.x + enlargementFactor * (cartesian3DPoints[i].x - polygonCenterCartesian3D.x);
                cartesian3DPoints[i].y = polygonCenterCartesian3D.y + enlargementFactor * (cartesian3DPoints[i].y - polygonCenterCartesian3D.y);
                cartesian3DPoints[i].z = polygonCenterCartesian3D.z + enlargementFactor * (cartesian3DPoints[i].z - polygonCenterCartesian3D.z);
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

            if (inside == 1)
            {
                isInPolygon = true;
            }
        }
        return isInPolygon;
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

    Vector GetDelta(const Point& firstPoint, const Point& secondPoint, const Projection::Type& projection)
    {
        if (projection == Projection::Type::Cartesian)
        {
            return GetDeltaCartesian(firstPoint, secondPoint);
        }

        // TODO some performance can be gained here, by combining the computing of dx and dy
        return Vector(GetDx(firstPoint, secondPoint, projection), GetDy(firstPoint, secondPoint, projection));
    }

    Vector ComputeNormalToline(const Point& start, const Point& end, const Projection::Type projection)
    {
        Vector direction = GetDelta(start, end, projection);
        direction.normalise();
        Vector normal(-direction.y(), direction.x());
        return normal;
    }

    double GetDx(const Point& firstPoint, const Point& secondPoint, const Projection::Type& projection)
    {
        const double delta = secondPoint.x - firstPoint.x;

        if (projection == Projection::Type::Cartesian)
        {
            return delta;
        }
        if (projection == Projection::Type::Spherical || projection == Projection::Type::SphericalAccurate)
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

    double GetDy(const Point& firstPoint, const Point& secondPoint, const Projection::Type& projection)
    {
        const double delta = secondPoint.y - firstPoint.y;

        if (projection == Projection::Type::Cartesian)
        {
            return delta;
        }
        if (projection == Projection::Type::Spherical || projection == Projection::Type::SphericalAccurate)
        {
            const double firstPointY = firstPoint.y * constants::conversion::degToRad;
            const double secondPointY = secondPoint.y * constants::conversion::degToRad;
            const double dy = constants::geometric::earth_radius * (secondPointY - firstPointY);
            return dy;
        }
        return constants::missing::doubleValue;
    }

    double OuterProductTwoSegments(const Point& firstPointFirstSegment, const Point& secondPointFirstSegment,
                                   const Point& firstPointSecondSegment, const Point& secondPointSecondSegment, const Projection::Type& projection)
    {
        if (projection == Projection::Type::SphericalAccurate)
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
        if (projection == Projection::Type::Cartesian || projection == Projection::Type::Spherical)
        {
            const double dx1 = GetDx(firstPointFirstSegment, secondPointFirstSegment, projection);
            const double dx2 = GetDx(firstPointSecondSegment, secondPointSecondSegment, projection);

            const double dy1 = GetDy(firstPointFirstSegment, secondPointFirstSegment, projection);
            const double dy2 = GetDy(firstPointSecondSegment, secondPointSecondSegment, projection);

            return dx1 * dy2 - dy1 * dx2;
        }

        return constants::missing::doubleValue;
    }

    Point ComputeMiddlePoint(const Point& firstPoint, const Point& secondPoint, const Projection::Type& projection)
    {
        if (!firstPoint.IsValid() || !secondPoint.IsValid())
        {
            return {constants::missing::doubleValue, constants::missing::doubleValue};
        }
        if (projection == Projection::Type::SphericalAccurate)
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
        if (projection == Projection::Type::Cartesian || projection == Projection::Type::Spherical)
        {
            return {(firstPoint + secondPoint) * 0.5};
        }
        return {constants::missing::doubleValue, constants::missing::doubleValue};
    }

    Point ComputeMiddlePointAccountingForPoles(const Point& firstPoint, const Point& secondPoint, const Projection::Type& projection)
    {
        if (!firstPoint.IsValid() || !secondPoint.IsValid())
        {
            return {constants::missing::doubleValue, constants::missing::doubleValue};
        }

        if (projection == Projection::Type::Spherical || projection == Projection::Type::SphericalAccurate)
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
        if (projection == Projection::Type::Cartesian)
        {
            return {(firstPoint + secondPoint) * 0.5};
        }
        return {constants::missing::doubleValue, constants::missing::doubleValue};
    }

    Point NormalVector(const Point& firstPoint, const Point& secondPoint, const Point& insidePoint, const Projection::Type& projection)
    {
        if (!firstPoint.IsValid() || !secondPoint.IsValid())
        {
            return {constants::missing::doubleValue, constants::missing::doubleValue};
        }

        if (projection == Projection::Type::SphericalAccurate)
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
        if (projection == Projection::Type::Cartesian || projection == Projection::Type::Spherical)
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

    void TransformGlobalVectorToLocal(const Point& reference, const Point& globalCoordinates, const Point& globalComponents, const Projection::Type& projection, Point& localComponents)
    {
        if (projection == Projection::Type::SphericalAccurate)
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
            if (projection == Projection::Type::Cartesian || projection == Projection::Type::Spherical)
            {
                localComponents = globalComponents;
            }
        }
    }

    Point NormalVectorOutside(const Point& firstPoint, const Point& secondPoint, const Projection::Type& projection)
    {
        if (!firstPoint.IsValid() || !secondPoint.IsValid())
        {
            return {constants::missing::doubleValue, constants::missing::doubleValue};
        }

        if (projection == Projection::Type::SphericalAccurate)
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
        if (projection == Projection::Type::Cartesian || projection == Projection::Type::Spherical)
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

            if (projection == Projection::Type::Spherical)
            {
                result.x = result.x / cos(constants::conversion::degToRad * 0.5 * (firstPoint.y + secondPoint.y));
                result.y = result.y;
            }
            return result;
        }
        return {constants::missing::doubleValue, constants::missing::doubleValue};
    }

    void NormalVectorInside(const Point& firstPoint, const Point& secondPoint, const Point& insidePoint, Point& normal, bool& flippedNormal, const Projection::Type& projection)
    {
        normal = NormalVectorOutside(firstPoint, secondPoint, projection);
        Point thirdPoint;
        if (projection == Projection::Type::Cartesian || projection == Projection::Type::Spherical)
        {
            flippedNormal = false;
            thirdPoint.x = firstPoint.x + normal.x;
            thirdPoint.y = firstPoint.y + normal.y;
        }
        if (projection == Projection::Type::SphericalAccurate)
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

    void AddIncrementToPoint(const Point& normal, double increment, const Point& referencePoint, const Projection::Type& projection, Point& point)
    {
        if (projection == Projection::Type::Cartesian)
        {
            point.x = point.x + normal.x * increment;
            point.y = point.y + normal.y * increment;
        }
        if (projection == Projection::Type::Spherical || projection == Projection::Type::SphericalAccurate)
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

    Point ReferencePoint(const std::vector<Point>& polygon, const Projection::Type& projection)
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

        if (projection == Projection::Type::Spherical)
        {
            if (maxX - minX > 180.0)
            {
                minX += 360.0;
            }
        }

        return Point{minX, minY};
    }

    double ComputeSquaredDistance(const Point& firstPoint, const Point& secondPoint, const Projection::Type& projection)
    {

        if (!firstPoint.IsValid() || !secondPoint.IsValid())
        {
            return 0.0;
        }

        if (projection == Projection::Type::SphericalAccurate)
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
        if (projection == Projection::Type::Cartesian || projection == Projection::Type::Spherical)
        {
            const double dx = GetDx(firstPoint, secondPoint, projection);
            const double dy = GetDy(firstPoint, secondPoint, projection);
            return dx * dx + dy * dy;
        }

        return constants::missing::doubleValue;
    }

    double ComputeDistance(const Point& firstPoint, const Point& secondPoint, const Projection::Type& projection)
    {
        double distance = ComputeSquaredDistance(firstPoint, secondPoint, projection);
        if (distance >= 0.0)
        {
            distance = sqrt(distance);
        }
        return distance;
    }

    std::tuple<double, Point, double> DistanceFromLine(const Point& point, const Point& firstNode, const Point& secondNode, const Projection::Type& projection)
    {
        double distance = constants::missing::doubleValue;
        Point normalPoint;
        double ratio = constants::missing::doubleValue;
        if (projection == Projection::Type::Cartesian || projection == Projection::Type::Spherical)
        {
            const auto squaredDistance = ComputeSquaredDistance(secondNode, firstNode, projection);
            if (squaredDistance != 0.0)
            {
                ratio = (GetDx(firstNode, point, projection) * GetDx(firstNode, secondNode, projection) +
                         GetDy(firstNode, point, projection) * GetDy(firstNode, secondNode, projection)) /
                        squaredDistance;
                const auto correctedRatio = std::max(std::min(1.0, ratio), 0.0);
                normalPoint.x = firstNode.x + correctedRatio * (secondNode.x - firstNode.x);
                normalPoint.y = firstNode.y + correctedRatio * (secondNode.y - firstNode.y);
                distance = ComputeDistance(point, normalPoint, projection);
            }
        }

        if (projection == Projection::Type::SphericalAccurate)
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

    double InnerProductTwoSegments(const Point& firstPointFirstSegment, const Point& secondPointFirstSegment, const Point& firstPointSecondSegment, const Point& secondPointSecondSegment, const Projection::Type& projection)
    {
        if (projection == Projection::Type::SphericalAccurate)
        {
            const Cartesian3DPoint firstPointFirstSegment3D{SphericalToCartesian3D(firstPointFirstSegment)};
            const Cartesian3DPoint secondPointFirstSegment3D{SphericalToCartesian3D(secondPointFirstSegment)};
            const Cartesian3DPoint firstPointSecondSegment3D{SphericalToCartesian3D(firstPointSecondSegment)};
            const Cartesian3DPoint secondPointSecondSegment3D{SphericalToCartesian3D(secondPointSecondSegment)};

            const double dx1 = secondPointFirstSegment3D.x - firstPointFirstSegment3D.x;
            const double dy1 = secondPointFirstSegment3D.y - firstPointFirstSegment3D.y;
            const double dz1 = secondPointFirstSegment3D.z - firstPointFirstSegment3D.z;

            const double dx2 = secondPointSecondSegment3D.x - firstPointSecondSegment3D.x;
            const double dy2 = secondPointSecondSegment3D.y - firstPointSecondSegment3D.y;
            const double dz2 = secondPointSecondSegment3D.z - firstPointSecondSegment3D.z;

            return DotProduct(dx1, dx2, dy1, dy2, dz1, dz2);
        }

        // cartesian and spherical
        if (projection == Projection::Type::Cartesian || projection == Projection::Type::Spherical)
        {
            const double dx1 = GetDx(firstPointFirstSegment, secondPointFirstSegment, projection);
            const double dx2 = GetDx(firstPointSecondSegment, secondPointSecondSegment, projection);

            double dy1 = GetDy(firstPointFirstSegment, secondPointFirstSegment, projection);
            double dy2 = GetDy(firstPointSecondSegment, secondPointSecondSegment, projection);

            return DotProduct(dx1, dx2, dy1, dy2);
        }

        return constants::missing::doubleValue;
    }

    double NormalizedInnerProductTwoSegments(const Point& firstPointFirstSegment, const Point& secondPointFirstSegment, const Point& firstPointSecondSegment, const Point& secondPointSecondSegment, const Projection::Type& projection)
    {
        if (projection == Projection::Type::SphericalAccurate)
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
        if (projection == Projection::Type::Cartesian || projection == Projection::Type::Spherical)
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

    Point CircumcenterOfTriangle(const Point& firstNode, const Point& secondNode, const Point& thirdNode, const Projection::Type& projection)
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
        if (projection == Projection::Type::Cartesian)
        {
            circumcenter.x = firstNode.x + 0.5 * (dx3 - z * dy3);
            circumcenter.y = firstNode.y + 0.5 * (dy3 + z * dx3);
        }
        if (projection == Projection::Type::Spherical)
        {
            const double phi = (firstNode.y + secondNode.y + thirdNode.y) * constants::numeric::oneThird;
            const double xf = 1.0 / cos(constants::conversion::degToRad * phi);
            circumcenter.x = firstNode.x + xf * 0.5 * (dx3 - z * dy3) * constants::conversion::radToDeg / constants::geometric::earth_radius;
            circumcenter.y = firstNode.y + 0.5 * (dy3 + z * dx3) * constants::conversion::radToDeg / constants::geometric::earth_radius;
        }
        if (projection == Projection::Type::SphericalAccurate)
        {
            // TODO: compute in case of spherical accurate (comp_circumcenter3D)
        }
        return circumcenter;
    }

    bool AreSegmentsCrossing(const Point& firstSegmentFirstPoint,
                             const Point& firstSegmentSecondPoint,
                             const Point& secondSegmentFirstPoint,
                             const Point& secondSegmentSecondPoint,
                             bool adimensionalCrossProduct,
                             const Projection::Type& projection,
                             Point& intersectionPoint,
                             double& crossProduct,
                             double& ratioFirstSegment,
                             double& ratioSecondSegment)
    {
        bool isCrossing = false;
        ratioFirstSegment = constants::missing::doubleValue;
        ratioSecondSegment = constants::missing::doubleValue;
        crossProduct = constants::missing::doubleValue;

        if (projection == Projection::Type::Cartesian || projection == Projection::Type::Spherical)
        {

            auto const x21 = GetDx(firstSegmentFirstPoint, firstSegmentSecondPoint, projection);
            auto const y21 = GetDy(firstSegmentFirstPoint, firstSegmentSecondPoint, projection);

            auto const x43 = GetDx(secondSegmentFirstPoint, secondSegmentSecondPoint, projection);
            auto const y43 = GetDy(secondSegmentFirstPoint, secondSegmentSecondPoint, projection);

            auto const x31 = GetDx(firstSegmentFirstPoint, secondSegmentFirstPoint, projection);
            auto const y31 = GetDy(firstSegmentFirstPoint, secondSegmentFirstPoint, projection);

            auto const det = x43 * y21 - y43 * x21;

            std::vector<double> values{std::abs(x21), std::abs(y21), std::abs(x43), std::abs(y43)};
            const double eps = std::max(0.00001 * *std::max_element(values.begin(), values.end()), std::numeric_limits<double>::denorm_min());

            if (std::abs(det) < eps)
            {
                return isCrossing;
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

        if (projection == Projection::Type::SphericalAccurate)
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

        return isCrossing;
    }

    std::tuple<std::vector<double>, double> ComputeAdimensionalDistancesFromPointSerie(const std::vector<Point>& v, const Projection::Type& projection)
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
                                                 const Projection::Type& projection,
                                                 UInt numM,
                                                 UInt numN)
    {
        const auto [sideOneAdimensional, totalLengthOne] = ComputeAdimensionalDistancesFromPointSerie(leftDiscretization, projection);
        const auto [sideTwoAdimensional, totalLengthTwo] = ComputeAdimensionalDistancesFromPointSerie(rightDiscretization, projection);
        const auto [sideThreeAdimensional, totalLengthThree] = ComputeAdimensionalDistancesFromPointSerie(bottomDiscretization, projection);
        const auto [sideFourAdimensional, totalLengthFour] = ComputeAdimensionalDistancesFromPointSerie(upperDiscretization, projection);

        // now compute the adimensional distance of each point to be filled
        const auto numMPoints = numM + 1;
        const auto numNPoints = numN + 1;

        lin_alg::Matrix<double> iWeightFactor(numMPoints, numNPoints);
        lin_alg::Matrix<double> jWeightFactor(numMPoints, numNPoints);

        for (UInt m = 0; m < numMPoints; m++)
        {
            for (UInt n = 0; n < numNPoints; n++)
            {
                const double mWeight = double(m) / double(numM);
                const double nWeight = double(n) / double(numN);

                iWeightFactor(m, n) = (1.0 - nWeight) * sideThreeAdimensional[m] + nWeight * sideFourAdimensional[m];
                jWeightFactor(m, n) = (1.0 - mWeight) * sideOneAdimensional[n] + mWeight * sideTwoAdimensional[n];
            }
        }

        lin_alg::Matrix<double> ones = lin_alg::Matrix<double>::Ones(numMPoints, numNPoints);
        lin_alg::Matrix<double> weightOne = (ones - jWeightFactor) * totalLengthThree + jWeightFactor * totalLengthFour;
        lin_alg::Matrix<double> weightTwo = (ones - iWeightFactor) * totalLengthOne + iWeightFactor * totalLengthTwo;
        lin_alg::Matrix<double> weightThree = weightTwo.cwiseQuotient(weightOne);
        lin_alg::Matrix<double> weightFour = weightThree.cwiseInverse();
        lin_alg::Matrix<double> const wa = (weightThree + weightFour).cwiseInverse();
        weightOne = wa.cwiseProduct(weightThree);
        weightTwo = wa.cwiseProduct(weightFour);

        // border points
        lin_alg::Matrix<Point> result(numMPoints, numNPoints);
        for (UInt m = 0; m < numMPoints; m++)
        {
            result(m, 0) = bottomDiscretization[m];
            result(m, numN) = upperDiscretization[m];
        }
        for (UInt n = 0; n < numNPoints; n++)
        {
            result(0, n) = leftDiscretization[n];
            result(numM, n) = rightDiscretization[n];
        }

        // first interpolation
        for (UInt m = 1; m < numM; m++)
        {
            for (UInt n = 1; n < numN; n++)
            {

                result(m, n).x = (leftDiscretization[n].x * (1.0 - iWeightFactor(m, n)) +
                                  rightDiscretization[n].x * iWeightFactor(m, n)) *
                                     weightOne(m, n) +
                                 (bottomDiscretization[m].x * (1.0 - jWeightFactor(m, n)) +
                                  upperDiscretization[m].x * jWeightFactor(m, n)) *
                                     weightTwo(m, n);

                result(m, n).y = (leftDiscretization[n].y * (1.0 - iWeightFactor(m, n)) +
                                  rightDiscretization[n].y * iWeightFactor(m, n)) *
                                     weightOne(m, n) +
                                 (bottomDiscretization[m].y * (1.0 - jWeightFactor(m, n)) +
                                  upperDiscretization[m].y * jWeightFactor(m, n)) *
                                     weightTwo(m, n);
            }
        }

        // update weights
        for (UInt m = 0; m < numMPoints; m++)
        {
            for (UInt n = 0; n < numNPoints; n++)
            {
                weightOne(m, n) = (1.0 - jWeightFactor(m, n)) * sideThreeAdimensional[m] * totalLengthThree +
                                  jWeightFactor(m, n) * sideFourAdimensional[m] * totalLengthFour;
                weightTwo(m, n) = (1.0 - iWeightFactor(m, n)) * sideOneAdimensional[n] * totalLengthOne +
                                  iWeightFactor(m, n) * sideTwoAdimensional[n] * totalLengthTwo;
            }
        }

        for (UInt m = 1; m < numMPoints; m++)
        {
            for (UInt n = 0; n < numNPoints; n++)
            {
                weightThree(m, n) = weightOne(m, n) - weightOne(m - 1, n);
            }
        }

        for (UInt m = 0; m < numMPoints; m++)
        {
            for (UInt n = 1; n < numNPoints; n++)
            {
                weightFour(m, n) = weightTwo(m, n) - weightTwo(m, n - 1);
            }
        }

        for (UInt m = 1; m < numMPoints; m++)
        {
            for (UInt n = 1; n < numNPoints - 1; n++)
            {
                weightOne(m, n) = 0.25 *
                                  (weightFour(m, n) +
                                   weightFour(m, n + 1) +
                                   weightFour(m - 1, n) +
                                   weightFour(m - 1, n + 1)) /
                                  weightThree(m, n);
            }
        }

        for (UInt m = 1; m < numMPoints - 1; m++)
        {
            for (UInt n = 1; n < numNPoints; n++)
            {
                weightTwo(m, n) = 0.25 *
                                  (weightThree(m, n) +
                                   weightThree(m, n - 1) +
                                   weightThree(m + 1, n) +
                                   weightThree(m + 1, n - 1)) /
                                  weightFour(m, n);
            }
        }

        // Iterate several times over
        // Move to constants
        static UInt constexpr numIterations = 25;

        for (UInt iter = 0; iter < numIterations; iter++)
        {
            // re-assign the weights
            for (UInt m = 0; m < numMPoints; m++)
            {
                for (UInt n = 0; n < numNPoints; n++)
                {
                    weightThree(m, n) = result(m, n).x;
                    weightFour(m, n) = result(m, n).y;
                }
            }

            for (UInt m = 1; m < numM; m++)
            {
                for (UInt n = 1; n < numN; n++)
                {

                    const double wa = 1.0 /
                                      (weightOne(m, n) +
                                       weightOne(m + 1, n) +
                                       weightTwo(m, n) +
                                       weightTwo(m, n + 1));

                    result(m, n).x = wa *
                                     (weightThree(m - 1, n) * weightOne(m, n) +
                                      weightThree(m + 1, n) * weightOne(m + 1, n) +
                                      weightThree(m, n - 1) * weightTwo(m, n) +
                                      weightThree(m, n + 1) * weightTwo(m, n + 1));

                    result(m, n).y = wa * (weightFour(m - 1, n) * weightOne(m, n) +
                                           weightFour(m + 1, n) * weightOne(m + 1, n) +
                                           weightFour(m, n - 1) * weightTwo(m, n) +
                                           weightFour(m, n + 1) * weightTwo(m, n + 1));
                }
            }
        }

        return result;
    }

    std::vector<Point> ComputeEdgeCenters(const std::vector<Point>& nodes, const std::vector<Edge>& edges)
    {
        std::vector<Point> edgesCenters;
        edgesCenters.reserve(edges.size());

        for (const auto& edge : edges)
        {
            auto const first = edge.first;
            auto const second = edge.second;
            if (first == constants::missing::uintValue || second == constants::missing::uintValue)
            {
                continue;
            }
            edgesCenters.emplace_back((nodes[first] + nodes[second]) * 0.5);
        }
        return edgesCenters;
    }

    double LinearInterpolationInTriangle(const Point& interpolationPoint, const std::vector<Point>& polygon, const std::vector<double>& values, const Projection::Type& projection)
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

    Point ComputeAverageCoordinate(const std::vector<Point>& points, const Projection::Type& projection)
    {
        std::vector<Point> validPoints;
        validPoints.reserve(points.size());
        for (const auto& point : points)
        {
            if (!point.IsValid())
            {
                continue;
            }
            validPoints.emplace_back(point);
        }

        if (projection == Projection::Type::SphericalAccurate)
        {

            Cartesian3DPoint averagePoint3D{0.0, 0.0, 0.0};
            for (const auto& point : validPoints)
            {
                const Cartesian3DPoint point3D{SphericalToCartesian3D(point)};
                averagePoint3D.x += point3D.x;
                averagePoint3D.y += point3D.y;
                averagePoint3D.z += point3D.z;
            }
            averagePoint3D.x = averagePoint3D.x / static_cast<double>(validPoints.size());
            averagePoint3D.y = averagePoint3D.y / static_cast<double>(validPoints.size());
            averagePoint3D.z = averagePoint3D.z / static_cast<double>(validPoints.size());

            return Cartesian3DToSpherical(averagePoint3D, points[0].x);
        }

        auto result = std::accumulate(validPoints.begin(), validPoints.end(), Point{0.0, 0.0});
        result.x = result.x / static_cast<double>(validPoints.size());
        result.y = result.y / static_cast<double>(validPoints.size());
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

    std::vector<double> ComputePolyLineEdgesLengths(const std::vector<Point>& polyline, Projection::Type projection)
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

    std::vector<double> ComputePolyLineNodalChainages(std::vector<Point> const& polyline, Projection::Type projection)
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

    std::vector<Point> ComputePolyLineDiscretization(std::vector<Point> const& polyline, std::vector<double>& chainages, Projection::Type projection)
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
