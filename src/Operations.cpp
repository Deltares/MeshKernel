//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2020.
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

#include <algorithm>
#include <cmath>

#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/SpatialTrees.hpp>

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

    std::vector<std::vector<size_t>> FindIndices(const std::vector<Point>& vec,
                                                 size_t start,
                                                 size_t end,
                                                 double separator)
    {
        std::vector<std::vector<size_t>> result;

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
        size_t startRange;
        for (auto n = start; n < end; n++)
        {

            if (!IsEqual(vec[n].x, separator) && !inRange)
            {
                startRange = n;
                inRange = true;
            }
            if (IsEqual(vec[n].x, separator) && inRange)
            {
                result.emplace_back(std::initializer_list<size_t>{startRange, n - 1});
                inRange = false;
            }
        }

        //in case no separator have been found
        if (inRange)
        {
            result.emplace_back(std::initializer_list<size_t>{startRange, vec.size() - 1});
        }

        return result;
    }

    size_t NextCircularForwardIndex(size_t currentIndex, size_t size)
    {
        size_t index = currentIndex + 1;
        if (index >= size)
        {
            index = index - size;
        }
        return index;
    }

    size_t NextCircularBackwardIndex(size_t currentIndex, size_t size)
    {
        if (currentIndex == 0)
        {
            return currentIndex + size - 1;
        }
        return currentIndex - 1;
    }

    bool IsPointOnPole(const Point& point)
    {
        return std::abs(std::abs(point.y) - 90.0) < absLatitudeAtPoles;
    }

    Cartesian3DPoint SphericalToCartesian3D(const Point& sphericalPoint)
    {
        Cartesian3DPoint result;
        result.z = earth_radius * sin(sphericalPoint.y * degrad_hp);
        double rr = earth_radius * cos(sphericalPoint.y * degrad_hp);
        result.x = rr * cos(sphericalPoint.x * degrad_hp);
        result.y = rr * sin(sphericalPoint.x * degrad_hp);
        return result;
    }

    Point Cartesian3DToSpherical(const Cartesian3DPoint& cartesianPoint, double referenceLongitude)
    {
        Point sphericalPoint;
        double angle = atan2(cartesianPoint.y, cartesianPoint.x) * raddeg_hp;
        sphericalPoint.y = atan2(cartesianPoint.z, sqrt(cartesianPoint.x * cartesianPoint.x + cartesianPoint.y * cartesianPoint.y)) * raddeg_hp;
        sphericalPoint.x = angle + std::lround((referenceLongitude - angle) / 360.0) * 360.0;
        return sphericalPoint;
    }

    double IsLeft(const Point& leftPoint, const Point& rightPoint, const Point& point)
    {
        double left = (rightPoint.x - leftPoint.x) * (point.y - leftPoint.y) - (point.x - leftPoint.x) * (rightPoint.y - leftPoint.y);
        return left;
    }

    bool IsPointInPolygonNodes(const Point& point,
                               const std::vector<Point>& polygonNodes,
                               int startNode,
                               int endNode,
                               const Projection& projection,
                               Point polygonCenter)
    {

        if (endNode <= startNode)
        {
            return true;
        }
        const int currentPolygonSize = endNode - startNode + 1;
        if (currentPolygonSize < numNodesInTriangle || polygonNodes.size() < currentPolygonSize)
        {
            return false;
        }
        if (polygonNodes[startNode] != polygonNodes[endNode])
        {
            return false;
        }

        bool isInPolygon = false;

        if (projection == Projection::cartesian || projection == Projection::spherical)
        {

            int windingNumber = 0;
            for (auto n = startNode; n < endNode; n++)
            {
                const auto leftDifference = IsLeft(polygonNodes[n], polygonNodes[n + 1], point);
                if (IsEqual(leftDifference, 0.0))
                {
                    // point on the line
                    return true;
                }

                if (polygonNodes[n].y <= point.y) // an upward crossing
                {
                    if (polygonNodes[n + 1].y > point.y && leftDifference > 0.0)

                    {
                        ++windingNumber; // have  a valid up intersect
                    }
                }
                else
                {
                    if (polygonNodes[n + 1].y <= point.y && leftDifference < 0.0) // a downward crossing
                    {

                        --windingNumber; // have  a valid down intersect
                    }
                }
            }

            isInPolygon = windingNumber == 0 ? false : true;
        }

        if (projection == Projection::sphericalAccurate)
        {
            // get 3D polygon coordinates
            std::vector<Cartesian3DPoint> cartesian3DPoints;
            cartesian3DPoints.reserve(currentPolygonSize);
            for (auto i = 0; i < currentPolygonSize; i++)
            {
                cartesian3DPoints.emplace_back(SphericalToCartesian3D(polygonNodes[startNode + i]));
            }

            // enlarge around polygon
            const double enlargementFactor = 1.000001;
            const Cartesian3DPoint polygonCenterCartesian3D{SphericalToCartesian3D(polygonCenter)};
            for (auto i = 0; i < currentPolygonSize; i++)
            {
                cartesian3DPoints[i].x = polygonCenterCartesian3D.x + enlargementFactor * (cartesian3DPoints[i].x - polygonCenterCartesian3D.x);
                cartesian3DPoints[i].y = polygonCenterCartesian3D.y + enlargementFactor * (cartesian3DPoints[i].y - polygonCenterCartesian3D.y);
                cartesian3DPoints[i].z = polygonCenterCartesian3D.z + enlargementFactor * (cartesian3DPoints[i].z - polygonCenterCartesian3D.z);
            }

            // convert point
            const Cartesian3DPoint pointCartesian3D{SphericalToCartesian3D(point)};

            //get test direction: e_lambda
            const double lambda = point.x * degrad_hp;
            Cartesian3DPoint ee{-std::sin(lambda), std::cos(lambda), 0.0};
            int inside = 0;

            // loop over the polygon nodes
            for (auto i = 0; i < currentPolygonSize - 1; i++)
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
        double phi0 = point.y * degrad_hp;
        double lambda0 = point.x * degrad_hp;

        exxp[0] = cos(phi0) * cos(lambda0);
        exxp[1] = cos(phi0) * sin(lambda0);
        exxp[2] = sin(phi0);

        eyyp[0] = -sin(lambda0);
        eyyp[1] = cos(lambda0);
        eyyp[2] = 0.0;

        ezzp[0] = -sin(phi0) * cos(lambda0);
        ezzp[1] = -sin(phi0) * sin(lambda0);
        ezzp[2] = cos(phi0);
    };

    void ComputeTwoBaseComponents(const Point& point, std::array<double, 3>& elambda, std::array<double, 3>& ephi)
    {
        double phi0 = point.y * degrad_hp;
        double lambda0 = point.x * degrad_hp;

        elambda[0] = -sin(lambda0);
        elambda[1] = cos(lambda0);
        elambda[2] = 0.0;

        ephi[0] = -sin(phi0) * cos(lambda0);
        ephi[1] = -sin(phi0) * sin(lambda0);
        ephi[2] = cos(phi0);
    };

    double GetDx(const Point& firstPoint, const Point& secondPoint, const Projection& projection)
    {
        double delta = secondPoint.x - firstPoint.x;
        if (std::abs(delta) <= nearlyZero)
        {
            return 0.0;
        }

        if (projection == Projection::cartesian)
        {
            return delta;
        }
        if (projection == Projection::spherical || projection == Projection::sphericalAccurate)
        {
            bool isFirstPointOnPole = IsPointOnPole(firstPoint);
            bool isSecondPointOnPole = IsPointOnPole(secondPoint);
            if (isFirstPointOnPole && !isSecondPointOnPole || !isFirstPointOnPole && isSecondPointOnPole)
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

            firstPointX = firstPointX * degrad_hp;
            secondPointX = secondPointX * degrad_hp;
            const double firstPointY = firstPoint.y * degrad_hp;
            const double secondPointY = secondPoint.y * degrad_hp;
            const double cosPhi = cos(0.5 * (firstPointY + secondPointY));
            const double dx = earth_radius * cosPhi * (secondPointX - firstPointX);
            return dx;
        }
        return doubleMissingValue;
    }

    double GetDy(const Point& firstPoint, const Point& secondPoint, const Projection& projection)
    {
        double delta = secondPoint.y - firstPoint.y;
        if (std::abs(delta) <= nearlyZero)
        {
            return 0.0;
        }

        if (projection == Projection::cartesian)
        {
            return delta;
        }
        if (projection == Projection::spherical || projection == Projection::sphericalAccurate)
        {
            double firstPointY = firstPoint.y * degrad_hp;
            double secondPointY = secondPoint.y * degrad_hp;
            double dy = earth_radius * (secondPointY - firstPointY);
            return dy;
        }
        return doubleMissingValue;
    }

    double OuterProductTwoSegments(const Point& firstPointFirstSegment, const Point& secondPointFirstSegment,
                                   const Point& firstPointSecondSegment, const Point& secondPointSecondSegment, const Projection& projection)
    {
        if (projection == Projection::sphericalAccurate)
        {
            const Cartesian3DPoint firstPointFirstSegmentCartesian{SphericalToCartesian3D(firstPointFirstSegment)};
            auto xx1 = firstPointFirstSegmentCartesian.x;
            auto yy1 = firstPointFirstSegmentCartesian.y;
            auto zz1 = firstPointFirstSegmentCartesian.z;

            const Cartesian3DPoint secondPointFirstSegmentCartesian{SphericalToCartesian3D(secondPointFirstSegment)};
            auto xx2 = secondPointFirstSegmentCartesian.x;
            auto yy2 = secondPointFirstSegmentCartesian.y;
            auto zz2 = secondPointFirstSegmentCartesian.z;

            const Cartesian3DPoint firstPointSecondSegmentCartesian{SphericalToCartesian3D(firstPointSecondSegment)};
            auto xx3 = firstPointSecondSegmentCartesian.x;
            auto yy3 = firstPointSecondSegmentCartesian.y;
            auto zz3 = firstPointSecondSegmentCartesian.z;

            const Cartesian3DPoint secondPointSecondSegmentCartesian{SphericalToCartesian3D(secondPointSecondSegment)};
            auto xx4 = secondPointSecondSegmentCartesian.x;
            auto yy4 = secondPointSecondSegmentCartesian.y;
            auto zz4 = secondPointSecondSegmentCartesian.z;

            double vxx = (yy2 - yy1) * (zz4 - zz3) - (zz2 - zz1) * (yy4 - yy3);
            double vyy = (zz2 - zz1) * (xx4 - xx3) - (xx2 - xx1) * (zz4 - zz3);
            double vzz = (xx2 - xx1) * (yy4 - yy3) - (yy2 - yy1) * (xx4 - xx3);

            double result = std::sqrt(vxx * vxx + vyy * vyy + vzz * vzz);

            //check if vector is pointing outwards of earth
            if (vxx * xx1 + vyy * yy1 + vzz * zz1 < 0.0)
            {
                result = -result;
            }
            return result;
        }

        // cartesian and spherical
        if (projection == Projection::cartesian || projection == Projection::spherical)
        {
            double dx1 = GetDx(firstPointFirstSegment, secondPointFirstSegment, projection);
            double dx2 = GetDx(firstPointSecondSegment, secondPointSecondSegment, projection);

            double dy1 = GetDy(firstPointFirstSegment, secondPointFirstSegment, projection);
            double dy2 = GetDy(firstPointSecondSegment, secondPointSecondSegment, projection);

            return dx1 * dy2 - dy1 * dx2;
        }

        return doubleMissingValue;
    }

    Point ComputeMiddlePoint(const Point& firstPoint, const Point& secondPoint, const Projection& projection)
    {
        if (!firstPoint.IsValid() || !secondPoint.IsValid())
        {
            return {doubleMissingValue, doubleMissingValue};
        }
        if (projection == Projection::sphericalAccurate)
        {
            const Cartesian3DPoint firstPointCartesianCoordinates{SphericalToCartesian3D(firstPoint)};
            const Cartesian3DPoint secondPointCartesianCoordinates{SphericalToCartesian3D(secondPoint)};

            Cartesian3DPoint middleCartesianPointCoordinate{doubleMissingValue, doubleMissingValue};
            middleCartesianPointCoordinate.x = 0.5 * (firstPointCartesianCoordinates.x + secondPointCartesianCoordinates.x);
            middleCartesianPointCoordinate.y = 0.5 * (firstPointCartesianCoordinates.y + secondPointCartesianCoordinates.y);
            double referenceLongitude = std::max(firstPoint.x, secondPoint.x);
            const auto result = Cartesian3DToSpherical(middleCartesianPointCoordinate, referenceLongitude);
            return result;
        }

        // cartesian and spherical
        if (projection == Projection::cartesian || projection == Projection::spherical)
        {
            return {(firstPoint + secondPoint) * 0.5};
        }
        return {doubleMissingValue, doubleMissingValue};
    }

    Point ComputeMiddlePointAccountingForPoles(const Point& firstPoint, const Point& secondPoint, const Projection& projection)
    {
        if (!firstPoint.IsValid() || !secondPoint.IsValid())
        {
            return {doubleMissingValue, doubleMissingValue};
        }

        if (projection == Projection::spherical || projection == Projection::sphericalAccurate)
        {
            Point result(0.0, (firstPoint.y + secondPoint.y) * 0.5);
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

                if (std::abs(maxx - minx) > std::numeric_limits<double>::epsilon())
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
        return {doubleMissingValue, doubleMissingValue};
    }

    Point NormalVector(const Point& firstPoint, const Point& secondPoint, const Point& insidePoint, const Projection& projection)
    {
        if (!firstPoint.IsValid() || !secondPoint.IsValid())
        {
            return {doubleMissingValue, doubleMissingValue};
        }

        if (projection == Projection::sphericalAccurate)
        {
            const Cartesian3DPoint firstPointCartesianCoordinates{SphericalToCartesian3D(firstPoint)};
            const Cartesian3DPoint secondPointCartesianCoordinates{SphericalToCartesian3D(secondPoint)};

            std::array<double, 3> elambda{0.0, 0.0, 0.0};
            std::array<double, 3> ephi{0.0, 0.0, 0.0};
            ComputeTwoBaseComponents(insidePoint, elambda, ephi);

            double dx = (secondPointCartesianCoordinates.x - firstPointCartesianCoordinates.x) * elambda[0] +
                        (secondPointCartesianCoordinates.y - firstPointCartesianCoordinates.y) * elambda[1] +
                        (secondPointCartesianCoordinates.z - firstPointCartesianCoordinates.z) * elambda[2];

            double dy = (secondPointCartesianCoordinates.x - firstPointCartesianCoordinates.x) * ephi[0] +
                        (secondPointCartesianCoordinates.y - firstPointCartesianCoordinates.y) * ephi[1] +
                        (secondPointCartesianCoordinates.z - firstPointCartesianCoordinates.z) * ephi[2];

            double squaredDistance = dx * dx + dy * dy;
            Point result{doubleMissingValue, doubleMissingValue};
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
            double dx = GetDx(firstPoint, secondPoint, projection);
            double dy = GetDy(firstPoint, secondPoint, projection);
            const double squaredDistance = dx * dx + dy * dy;
            Point result{doubleMissingValue, doubleMissingValue};
            if (squaredDistance > 0.0)
            {
                const double distance = sqrt(squaredDistance);
                result.x = dx / distance;
                result.y = dy / distance;
            }
            return result;
        }
        return {doubleMissingValue, doubleMissingValue};
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

            //project to rotated frame
            Cartesian3DPoint globalCoordinatesCartesianRotated;
            globalCoordinatesCartesianRotated.x = exxp[0] * globalCoordinatesCartesian.x + exxp[1] * globalCoordinatesCartesian.y + exxp[2] * globalCoordinatesCartesian.z;
            globalCoordinatesCartesianRotated.y = eyyp[0] * globalCoordinatesCartesian.x + eyyp[1] * globalCoordinatesCartesian.y + eyyp[2] * globalCoordinatesCartesian.z;
            globalCoordinatesCartesianRotated.z = ezzp[0] * globalCoordinatesCartesian.x + ezzp[1] * globalCoordinatesCartesian.y + ezzp[2] * globalCoordinatesCartesian.z;

            // Compute global base vectors at other point in 3D(xx, yy, zz) frame
            std::array<double, 3> elambda{0.0, 0.0, 0.0};
            std::array<double, 3> ephi{0.0, 0.0, 0.0};
            ComputeTwoBaseComponents(globalCoordinates, elambda, ephi);

            double vxx = globalComponents.x * elambda[0] + globalComponents.y * ephi[0];
            double vyy = globalComponents.x * elambda[1] + globalComponents.y * ephi[1];
            double vzz = globalComponents.x * elambda[2] + globalComponents.y * ephi[2];

            //transform to local spherical coordinates
            const auto globalCoordinatesToLocal = Cartesian3DToSpherical(globalCoordinatesCartesianRotated, reference.x);

            //compute base vectors at other point in rotated 3D(xxp, yyp, zzp) frame
            std::array<double, 3> elambdap{0.0, 0.0, 0.0};
            std::array<double, 3> ephip{0.0, 0.0, 0.0};
            ComputeTwoBaseComponents(globalCoordinatesToLocal, elambdap, ephip);

            //compute local base vectors in(xx, yy, zz) frame
            std::array<double, 3> elambdaloc;
            elambdaloc[0] = exxp[0] * elambdap[0] + eyyp[0] * elambdap[1] + ezzp[0] * elambda[2];
            elambdaloc[1] = exxp[1] * elambdap[0] + eyyp[1] * elambdap[1] + ezzp[1] * elambda[2];
            elambdaloc[2] = exxp[2] * elambdap[0] + eyyp[2] * elambdap[1] + ezzp[2] * elambda[2];

            std::array<double, 3> ephiloc;
            ephiloc[0] = exxp[0] * ephip[0] + eyyp[0] * ephip[1] + ezzp[0] * ephip[2];
            ephiloc[1] = exxp[1] * ephip[0] + eyyp[1] * ephip[1] + ezzp[1] * ephip[2];
            ephiloc[2] = exxp[2] * ephip[0] + eyyp[2] * ephip[1] + ezzp[2] * ephip[2];

            //compute vectors in other point in local base(elambdaloc, ephiloc)
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
            return {doubleMissingValue, doubleMissingValue};
        }

        if (projection == Projection::sphericalAccurate)
        {
            const auto middlePoint = ComputeMiddlePoint(firstPoint, secondPoint, projection);

            const Cartesian3DPoint firstPointCartesianCoordinates{SphericalToCartesian3D(firstPoint)};
            const Cartesian3DPoint secondPointCartesianCoordinates{SphericalToCartesian3D(secondPoint)};

            //compute the base vectors at middle point
            std::array<double, 3> elambda{0.0, 0.0, 0.0};
            std::array<double, 3> ephi{0.0, 0.0, 0.0};
            ComputeTwoBaseComponents(middlePoint, elambda, ephi);

            // project vector in local base
            double dx = (secondPointCartesianCoordinates.x - firstPointCartesianCoordinates.x) * elambda[0] +
                        (secondPointCartesianCoordinates.y - firstPointCartesianCoordinates.y) * elambda[1] +
                        (secondPointCartesianCoordinates.z - firstPointCartesianCoordinates.z) * elambda[2];

            double dy = (secondPointCartesianCoordinates.x - firstPointCartesianCoordinates.x) * ephi[0] +
                        (secondPointCartesianCoordinates.y - firstPointCartesianCoordinates.y) * ephi[1] +
                        (secondPointCartesianCoordinates.z - firstPointCartesianCoordinates.z) * ephi[2];

            const double squaredDistance = dx * dx + dy * dy;
            Point result{doubleMissingValue, doubleMissingValue};
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
            double dx = GetDx(firstPoint, secondPoint, projection);
            double dy = GetDy(firstPoint, secondPoint, projection);

            const double squaredDistance = dx * dx + dy * dy;
            Point result{doubleMissingValue, doubleMissingValue};
            if (squaredDistance > 0.0)
            {
                const double distance = sqrt(squaredDistance);
                result.x = dy / distance;
                result.y = -dx / distance;
            }

            if (projection == Projection::spherical)
            {
                const auto middlePoint = ComputeMiddlePoint(firstPoint, secondPoint, projection);
                result.x = result.x / cos(degrad_hp * 0.5 * (firstPoint.y + secondPoint.y));
                result.y = result.y;
            }
            return result;
        }
        return {doubleMissingValue, doubleMissingValue};
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

            double vxx = localComponents.x * elambda[0] + localComponents.y * ephi[0];
            double vyy = localComponents.x * elambda[1] + localComponents.y * ephi[1];
            double vzz = localComponents.x * elambda[2] + localComponents.y * ephi[2];

            const Cartesian3DPoint firstPointCartesian{SphericalToCartesian3D(firstPoint)};

            Cartesian3DPoint rotatedPoint;
            double alpha = 0.0;
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
            const double xf = 1.0 / std::cos(degrad_hp * referencePoint.y);
            const double convertedIncrement = raddeg_hp * increment / earth_radius;
            point.x = point.x + normal.x * convertedIncrement * xf;
            point.y = point.y + normal.y * convertedIncrement;
        }
    }

    Point ReferencePoint(std::vector<Point>& polygon, const Projection& projection)
    {
        auto minX = std::numeric_limits<double>::max();
        auto minY = std::numeric_limits<double>::max();
        const auto numPoints = polygon.size();
        for (auto i = 0; i < numPoints; i++)
        {
            minX = std::min(polygon[i].x, minX);
            if (abs(polygon[i].y) < abs(minY))
            {
                minY = polygon[i].y;
            }
        }

        if (projection == Projection::spherical)
        {
            double maxX = std::numeric_limits<double>::lowest();
            for (auto i = 0; i < numPoints; i++)
            {
                maxX = std::max(polygon[i].x, maxX);
            }

            if (maxX - minX > 180.0)
            {
                double deltaX = maxX - 180.0;
                for (auto i = 0; i < numPoints; i++)
                {
                    if (polygon[i].x < deltaX)
                    {
                        polygon[i].x = polygon[i].x + 360.0;
                    }
                }
                minX = minX + 360.0;
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
            auto xx1 = firstPointCartesian.x;
            auto yy1 = firstPointCartesian.y;
            auto zz1 = firstPointCartesian.z;

            const Cartesian3DPoint secondPointCartesian{SphericalToCartesian3D(secondPoint)};
            auto xx2 = secondPointCartesian.x;
            auto yy2 = secondPointCartesian.y;
            auto zz2 = secondPointCartesian.z;

            return (xx2 - xx1) * (xx2 - xx1) + (yy2 - yy1) * (yy2 - yy1) + (zz2 - zz1) * (zz2 - zz1);
        }

        //cartesian and spherical
        if (projection == Projection::cartesian || projection == Projection::spherical)
        {
            double dx = GetDx(firstPoint, secondPoint, projection);
            double dy = GetDy(firstPoint, secondPoint, projection);
            return dx * dx + dy * dy;
        }

        return doubleMissingValue;
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

    double DistanceFromLine(const Point& point, const Point& firstNode, const Point& secondNode, const Projection& projection, Point& normalPoint, double& ratio)
    {
        if (projection == Projection::cartesian || projection == Projection::spherical)
        {
            double dis = 0.0;
            double squaredDistance = ComputeSquaredDistance(secondNode, firstNode, projection);
            if (squaredDistance != 0.0)
            {
                ratio = (GetDx(firstNode, point, projection) * GetDx(firstNode, secondNode, projection) +
                         GetDy(firstNode, point, projection) * GetDy(firstNode, secondNode, projection)) /
                        squaredDistance;
                const auto correctedRatio = std::max(std::min(1.0, ratio), 0.0);
                normalPoint.x = firstNode.x + correctedRatio * (secondNode.x - firstNode.x);
                normalPoint.y = firstNode.y + correctedRatio * (secondNode.y - firstNode.y);
                dis = ComputeDistance(point, normalPoint, projection);
            }
            return dis;
        }

        if (projection == Projection::sphericalAccurate)
        {
            const Cartesian3DPoint firstNodeCartesian{SphericalToCartesian3D(firstNode)};
            auto xx1 = firstNodeCartesian.x;
            auto yy1 = firstNodeCartesian.y;
            auto zz1 = firstNodeCartesian.z;

            const Cartesian3DPoint secondNodeCartesian{SphericalToCartesian3D(secondNode)};
            auto xx2 = secondNodeCartesian.x;
            auto yy2 = secondNodeCartesian.y;
            auto zz2 = secondNodeCartesian.z;

            const Cartesian3DPoint pointCartesian{SphericalToCartesian3D(point)};
            auto xx3 = pointCartesian.x;
            auto yy3 = pointCartesian.y;
            auto zz3 = pointCartesian.z;

            double x21 = xx2 - xx1;
            double y21 = yy2 - yy1;
            double z21 = zz2 - zz1;
            double x31 = xx3 - xx1;
            double y31 = yy3 - yy1;
            double z31 = zz3 - zz1;

            double r2 = x21 * x21 + y21 * y21 + z21 * z21;

            ratio = 0.0;
            if (r2 >= 0.0)
            {

                ratio = (x31 * x21 + y31 * y21 + z31 * z21) / r2;
                double correctedRatio = std::max(std::min(1.0, ratio), 0.0);

                Cartesian3DPoint cartesianNormal3DPoint;
                cartesianNormal3DPoint.x = firstNodeCartesian.x + correctedRatio * x21;
                cartesianNormal3DPoint.y = firstNodeCartesian.y + correctedRatio * y21;
                cartesianNormal3DPoint.z = firstNodeCartesian.z + correctedRatio * z21;

                cartesianNormal3DPoint.x = cartesianNormal3DPoint.x - xx3;
                cartesianNormal3DPoint.y = cartesianNormal3DPoint.y - yy3;
                cartesianNormal3DPoint.z = cartesianNormal3DPoint.z - zz3;

                double dis = std::sqrt(cartesianNormal3DPoint.x * cartesianNormal3DPoint.x +
                                       cartesianNormal3DPoint.y * cartesianNormal3DPoint.y +
                                       cartesianNormal3DPoint.z * cartesianNormal3DPoint.z);

                double referenceLongitude = std::max({firstNode.x, secondNode.x, point.x});
                normalPoint = Cartesian3DToSpherical(cartesianNormal3DPoint, referenceLongitude);

                return dis;
            }
        }

        return -1.0;
    }

    double InnerProductTwoSegments(const Point& firstPointFirstSegment, const Point& secondPointFirstSegment, const Point& firstPointSecondSegment, const Point& secondPointSecondSegment, const Projection& projection)
    {
        if (projection == Projection::sphericalAccurate)
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
        if (projection == Projection::cartesian || projection == Projection::spherical)
        {
            double dx1 = GetDx(firstPointFirstSegment, secondPointFirstSegment, projection);
            double dx2 = GetDx(firstPointSecondSegment, secondPointSecondSegment, projection);

            double dy1 = GetDy(firstPointFirstSegment, secondPointFirstSegment, projection);
            double dy2 = GetDy(firstPointSecondSegment, secondPointSecondSegment, projection);

            return DotProduct(dx1, dx2, dy1, dy2);
        }

        return doubleMissingValue;
    }

    double NormalizedInnerProductTwoSegments(const Point& firstPointFirstSegment, const Point& secondPointFirstSegment, const Point& firstPointSecondSegment, const Point& secondPointSecondSegment, const Projection& projection)
    {
        if (projection == Projection::sphericalAccurate)
        {
            const Cartesian3DPoint firstPointFirstSegmentCartesian{SphericalToCartesian3D(firstPointFirstSegment)};
            auto xx1 = firstPointFirstSegmentCartesian.x;
            auto yy1 = firstPointFirstSegmentCartesian.y;
            auto zz1 = firstPointFirstSegmentCartesian.z;

            const Cartesian3DPoint secondPointFirstSegmentCartesian{SphericalToCartesian3D(secondPointFirstSegment)};
            auto xx2 = secondPointFirstSegmentCartesian.x;
            auto yy2 = secondPointFirstSegmentCartesian.y;
            auto zz2 = secondPointFirstSegmentCartesian.z;

            const Cartesian3DPoint firstPointSecondSegmentCartesian{SphericalToCartesian3D(firstPointSecondSegment)};
            auto xx3 = firstPointSecondSegmentCartesian.x;
            auto yy3 = firstPointSecondSegmentCartesian.y;
            auto zz3 = firstPointSecondSegmentCartesian.z;

            const Cartesian3DPoint secondPointSecondSegmentCartesian{SphericalToCartesian3D(secondPointSecondSegment)};
            auto xx4 = secondPointSecondSegmentCartesian.x;
            auto yy4 = secondPointSecondSegmentCartesian.y;
            auto zz4 = secondPointSecondSegmentCartesian.z;

            auto dx1 = xx2 - xx1;
            auto dy1 = yy2 - yy1;
            auto dz1 = zz2 - zz1;
            auto firstSegmentDistance = dx1 * dx1 + dy1 * dy1 + dz1 * dz1;

            auto dx2 = xx4 - xx3;
            auto dy2 = yy4 - yy3;
            auto dz2 = zz4 - zz3;
            auto secondSegmentDistance = dx2 * dx2 + dy2 * dy2 + dz2 * dz2;

            double cosphi;
            if (firstSegmentDistance <= 0.0 || secondSegmentDistance <= 0.0)
            {
                cosphi = doubleMissingValue;
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
                return doubleMissingValue;
            }

            cosphi = (dx1 * dx2 + dy1 * dy2) / std::sqrt(r1 * r2);
            cosphi = std::max(std::min(cosphi, 1.0), -1.0);
            return cosphi;
        }

        return doubleMissingValue;
    }

    Point CircumcenterOfTriangle(const Point& firstNode, const Point& secondNode, const Point& thirdNode, const Projection& projection)
    {

        double dx2 = GetDx(firstNode, secondNode, projection);
        double dy2 = GetDy(firstNode, secondNode, projection);

        double dx3 = GetDx(firstNode, thirdNode, projection);
        double dy3 = GetDy(firstNode, thirdNode, projection);

        double den = dy2 * dx3 - dy3 * dx2;
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
            double phi = (firstNode.y + secondNode.y + thirdNode.y) * oneThird;
            double xf = 1.0 / cos(degrad_hp * phi);
            circumcenter.x = firstNode.x + xf * 0.5 * (dx3 - z * dy3) * raddeg_hp / earth_radius;
            circumcenter.y = firstNode.y + 0.5 * (dy3 + z * dx3) * raddeg_hp / earth_radius;
        }
        if (projection == Projection::sphericalAccurate)
        {
            //TODO: compute in case of spherical accurate
        }
        return circumcenter;
    }

    bool AreSegmentsCrossing(const Point& firstSegmentFistPoint,
                             const Point& firstSegmentSecondPoint,
                             const Point& secondSegmentFistPoint,
                             const Point& secondSegmentSecondPoint,
                             bool adimensionalCrossProduct,
                             const Projection& projection,
                             Point& intersectionPoint,
                             double& crossProduct,
                             double& ratioFirstSegment,
                             double& ratioSecondSegment)
    {
        bool isCrossing = false;
        ratioFirstSegment = doubleMissingValue;
        ratioSecondSegment = doubleMissingValue;
        crossProduct = doubleMissingValue;

        if (projection == Projection::cartesian || projection == Projection::spherical)
        {

            auto const x21 = GetDx(firstSegmentFistPoint, firstSegmentSecondPoint, projection);
            auto const y21 = GetDy(firstSegmentFistPoint, firstSegmentSecondPoint, projection);

            auto const x43 = GetDx(secondSegmentFistPoint, secondSegmentSecondPoint, projection);
            auto const y43 = GetDy(secondSegmentFistPoint, secondSegmentSecondPoint, projection);

            auto const x31 = GetDx(firstSegmentFistPoint, secondSegmentFistPoint, projection);
            auto const y31 = GetDy(firstSegmentFistPoint, secondSegmentFistPoint, projection);

            auto const det = x43 * y21 - y43 * x21;

            std::vector<double> values{x21, y21, x43, y43};
            double eps = std::max(0.00001 * (*std::max_element(values.begin(), values.end())), std::numeric_limits<double>::denorm_min());

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
            intersectionPoint.x = firstSegmentFistPoint.x + ratioFirstSegment * (firstSegmentSecondPoint.x - firstSegmentFistPoint.x);
            intersectionPoint.y = firstSegmentFistPoint.y + ratioFirstSegment * (firstSegmentSecondPoint.y - firstSegmentFistPoint.y);
            crossProduct = -det;
            if (adimensionalCrossProduct)
            {
                crossProduct = -det / (std::sqrt(x21 * x21 + y21 * y21) * std::sqrt(x43 * x43 + y43 * y43) + 1e-8);
            }
        }

        if (projection == Projection::sphericalAccurate)
        {
            const Cartesian3DPoint firstSegmentFistCartesian3DPoint{SphericalToCartesian3D(firstSegmentFistPoint)};

            const Cartesian3DPoint firstSegmentSecondCartesian3DPoint{SphericalToCartesian3D(firstSegmentSecondPoint)};

            const Cartesian3DPoint secondSegmentFistCartesian3DPoint{SphericalToCartesian3D(secondSegmentFistPoint)};

            const Cartesian3DPoint secondSegmentSecondCartesian3DPoint{SphericalToCartesian3D(secondSegmentSecondPoint)};

            auto n12 = VectorProduct(firstSegmentFistCartesian3DPoint, firstSegmentSecondCartesian3DPoint);
            const auto n12InnerProduct = std::sqrt(InnerProduct(n12, n12));
            n12.x = n12.x / n12InnerProduct;
            n12.y = n12.y / n12InnerProduct;
            n12.z = n12.z / n12InnerProduct;

            auto n34 = VectorProduct(secondSegmentFistCartesian3DPoint, secondSegmentSecondCartesian3DPoint);
            const auto n34InnerProduct = std::sqrt(InnerProduct(n34, n34));
            n34.x = n34.x / n34InnerProduct;
            n34.y = n34.y / n34InnerProduct;
            n34.z = n34.z / n34InnerProduct;

            const auto n12n34InnerProduct = std::sqrt(std::abs(InnerProduct(n12, n34)));

            const double tolerance = 1e-12;
            if (n12n34InnerProduct > tolerance)
            {
                Cartesian3DPoint firstSegmentDifference;
                firstSegmentDifference.x = firstSegmentSecondCartesian3DPoint.x - firstSegmentFistCartesian3DPoint.x;
                firstSegmentDifference.y = firstSegmentSecondCartesian3DPoint.y - firstSegmentFistCartesian3DPoint.y;
                firstSegmentDifference.z = firstSegmentSecondCartesian3DPoint.z - firstSegmentFistCartesian3DPoint.z;

                Cartesian3DPoint secondSegmentDifference;
                secondSegmentDifference.x = secondSegmentSecondCartesian3DPoint.x - secondSegmentFistCartesian3DPoint.x;
                secondSegmentDifference.y = secondSegmentSecondCartesian3DPoint.y - secondSegmentFistCartesian3DPoint.y;
                secondSegmentDifference.z = secondSegmentSecondCartesian3DPoint.z - secondSegmentFistCartesian3DPoint.z;

                const auto Det12 = InnerProduct(firstSegmentDifference, n34);
                const auto Det34 = InnerProduct(secondSegmentDifference, n12);

                if (std::abs(Det12) > tolerance && std::abs(Det34) > tolerance)
                {
                    ratioFirstSegment = -InnerProduct(firstSegmentFistCartesian3DPoint, n34) / Det12;
                    ratioSecondSegment = -InnerProduct(secondSegmentFistCartesian3DPoint, n12) / Det34;
                }
            }

            if (ratioSecondSegment >= 0.0 && ratioSecondSegment <= 1.0 &&
                ratioFirstSegment >= 0.0 && ratioFirstSegment <= 1.0)
            {
                // check if segments are crossing
                isCrossing = true;

                // compute intersection
                Cartesian3DPoint intersectionCartesian3DPoint;
                intersectionCartesian3DPoint.x = firstSegmentFistCartesian3DPoint.x + ratioFirstSegment * (firstSegmentSecondCartesian3DPoint.x - firstSegmentFistCartesian3DPoint.x);
                intersectionCartesian3DPoint.y = firstSegmentFistCartesian3DPoint.y + ratioFirstSegment * (firstSegmentSecondCartesian3DPoint.y - firstSegmentFistCartesian3DPoint.y);
                intersectionCartesian3DPoint.z = firstSegmentFistCartesian3DPoint.z + ratioFirstSegment * (firstSegmentSecondCartesian3DPoint.z - firstSegmentFistCartesian3DPoint.z);
                intersectionPoint = Cartesian3DToSpherical(intersectionCartesian3DPoint, std::max(firstSegmentFistPoint.x, firstSegmentSecondPoint.x));
            }
        }

        return isCrossing;
    }

    int CrossProductSign(const Point& firstSegmentFistPoint, const Point& firstSegmentSecondPoint, const Point& secondSegmentFistPoint, const Point& secondSegmentSecondPoint, const Projection& projection)
    {

        const auto dx1 = GetDx(firstSegmentFistPoint, firstSegmentSecondPoint, projection);
        const auto dy1 = GetDy(firstSegmentFistPoint, firstSegmentSecondPoint, projection);
        const auto dx2 = GetDx(secondSegmentFistPoint, secondSegmentSecondPoint, projection);
        const auto dy2 = GetDy(secondSegmentFistPoint, secondSegmentSecondPoint, projection);
        const auto val = dx1 * dy2 - dy1 * dx2;
        return sgn(val);
    }

    void FaceAreaAndCenterOfMass(std::vector<Point>& polygon, const Projection& projection, double& area, Point& centerOfMass, bool& isCounterClockWise)
    {
        if (polygon.empty())
        {
            throw std::invalid_argument("FaceAreaAndCenterOfMass: The polygon contains no nodes.");
        }

        if (polygon.size() - 1 < numNodesInTriangle)
        {
            throw std::invalid_argument("FaceAreaAndCenterOfMass: The polygon has less than 3 unique nodes.");
        }

        area = 0.0;
        double xCenterOfMass = 0.0;
        double yCenterOfMass = 0.0;
        const double minArea = 1e-8;
        const Point reference = ReferencePoint(polygon, projection);
        const auto numberOfPointsOpenedPolygon = polygon.size() - 1;
        for (auto n = 0; n < numberOfPointsOpenedPolygon; n++)
        {
            const auto nextNode = NextCircularForwardIndex(n, numberOfPointsOpenedPolygon);
            double dx0 = GetDx(reference, polygon[n], projection);
            double dy0 = GetDy(reference, polygon[n], projection);
            double dx1 = GetDx(reference, polygon[nextNode], projection);
            double dy1 = GetDy(reference, polygon[nextNode], projection);

            double xc = 0.5 * (dx0 + dx1);
            double yc = 0.5 * (dy0 + dy1);

            dx0 = GetDx(polygon[n], polygon[nextNode], projection);
            dy0 = GetDy(polygon[n], polygon[nextNode], projection);
            double dsx = dy0;
            double dsy = -dx0;
            double xds = xc * dsx + yc * dsy;
            area = area + 0.5 * xds;

            xCenterOfMass = xCenterOfMass + xds * xc;
            yCenterOfMass = yCenterOfMass + xds * yc;
        }

        isCounterClockWise = area > 0.0;

        area = std::abs(area) < minArea ? minArea : area;

        const double fac = 1.0 / (3.0 * area);
        xCenterOfMass = fac * xCenterOfMass;
        yCenterOfMass = fac * yCenterOfMass;

        if (projection == Projection::spherical)
        {
            yCenterOfMass = yCenterOfMass / (earth_radius * degrad_hp);
            xCenterOfMass = xCenterOfMass / (earth_radius * degrad_hp * std::cos((yCenterOfMass + reference.y) * degrad_hp));
        }

        centerOfMass.x = xCenterOfMass + reference.x;
        centerOfMass.y = yCenterOfMass + reference.y;

        area = std::abs(area);
    }

    void ComputeAdimensionalDistancesFromPointSerie(const std::vector<Point>& v, const Projection& projection, std::vector<double>& result, double& totalDistance)
    {
        result[0] = 0;
        for (auto i = 1; i < v.size(); i++)
        {
            result[i] = result[i - 1] + ComputeDistance(v[i - 1], v[i], projection);
        }
        totalDistance = result.back();
        if (IsEqual(totalDistance, 0.0))
        {
            return;
        }
        const double inverseTotalDistance = 1.0 / totalDistance;
        for (auto i = 1; i < v.size(); i++)
        {
            result[i] = result[i] * inverseTotalDistance;
        }
    }

    std::vector<std::vector<Point>> DiscretizeTransfinite(const std::vector<Point>& sideOne,
                                                          const std::vector<Point>& sideTwo,
                                                          const std::vector<Point>& sideThree,
                                                          const std::vector<Point>& sideFour,
                                                          const Projection& projection,
                                                          int numM,
                                                          int numN)
    {
        double totalLengthOne;
        std::vector<double> sideOneAdimensional(sideOne.size());
        ComputeAdimensionalDistancesFromPointSerie(sideOne, projection, sideOneAdimensional, totalLengthOne);

        double totalLengthTwo;
        std::vector<double> sideTwoAdimensional(sideTwo.size());
        ComputeAdimensionalDistancesFromPointSerie(sideTwo, projection, sideTwoAdimensional, totalLengthTwo);

        double totalLengthThree;
        std::vector<double> sideThreeAdimensional(sideThree.size());
        ComputeAdimensionalDistancesFromPointSerie(sideThree, projection, sideThreeAdimensional, totalLengthThree);

        double totalLengthFour;
        std::vector<double> sideFourAdimensional(sideFour.size());
        ComputeAdimensionalDistancesFromPointSerie(sideFour, projection, sideFourAdimensional, totalLengthFour);

        // now compute the adimensional distance of each point to be filled
        const int numMPoints = numM + 1;
        const int numNPoints = numN + 1;

        std::vector<std::vector<double>> iWeightFactor(numMPoints, std::vector<double>(numNPoints));
        std::vector<std::vector<double>> jWeightFactor(numMPoints, std::vector<double>(numNPoints));
        for (auto i = 0; i < numMPoints; i++)
        {
            for (auto j = 0; j < numNPoints; j++)
            {
                const double mWeight = double(i) / double(numM);
                const double nWeight = double(j) / double(numN);

                iWeightFactor[i][j] = (1.0 - nWeight) * sideThreeAdimensional[i] + nWeight * sideFourAdimensional[i];
                jWeightFactor[i][j] = (1.0 - mWeight) * sideOneAdimensional[j] + mWeight * sideTwoAdimensional[j];
            }
        }

        std::vector<std::vector<double>> weightOne(numMPoints, std::vector<double>(numNPoints));
        std::vector<std::vector<double>> weightTwo(numMPoints, std::vector<double>(numNPoints));
        std::vector<std::vector<double>> weightThree(numMPoints, std::vector<double>(numNPoints));
        std::vector<std::vector<double>> weightFour(numMPoints, std::vector<double>(numNPoints));
        for (auto i = 0; i < numMPoints; i++)
        {
            for (auto j = 0; j < numNPoints; j++)
            {

                weightOne[i][j] = (1.0 - jWeightFactor[i][j]) * totalLengthThree + jWeightFactor[i][j] * totalLengthFour;
                weightTwo[i][j] = (1.0 - iWeightFactor[i][j]) * totalLengthOne + iWeightFactor[i][j] * totalLengthTwo;
                weightThree[i][j] = weightTwo[i][j] / weightOne[i][j];
                weightFour[i][j] = weightOne[i][j] / weightTwo[i][j];
                const double wa = 1.0 / (weightThree[i][j] + weightFour[i][j]);
                weightOne[i][j] = wa * weightThree[i][j];
                weightTwo[i][j] = wa * weightFour[i][j];
            }
        }

        //border points
        std::vector<std::vector<Point>> result(numMPoints, std::vector<Point>(numNPoints));
        for (auto i = 0; i < numMPoints; i++)
        {
            result[i][0] = sideThree[i];
            result[i][numN] = sideFour[i];
        }
        for (auto i = 0; i < numNPoints; i++)
        {
            result[0][i] = sideOne[i];
            result[numM][i] = sideTwo[i];
        }

        // first interpolation
        for (auto i = 1; i < numM; i++)
        {
            for (auto j = 1; j < numN; j++)
            {

                result[i][j].x = (sideOne[j].x * (1.0 - iWeightFactor[i][j]) + sideTwo[j].x * iWeightFactor[i][j]) * weightOne[i][j] +
                                 (sideThree[i].x * (1.0 - jWeightFactor[i][j]) + sideFour[i].x * jWeightFactor[i][j]) * weightTwo[i][j];

                result[i][j].y = (sideOne[j].y * (1.0 - iWeightFactor[i][j]) + sideTwo[j].y * iWeightFactor[i][j]) * weightOne[i][j] +
                                 (sideThree[i].y * (1.0 - jWeightFactor[i][j]) + sideFour[i].y * jWeightFactor[i][j]) * weightTwo[i][j];
            }
        }

        // update weights
        for (auto i = 0; i < numMPoints; i++)
        {
            for (auto j = 0; j < numNPoints; j++)
            {
                weightOne[i][j] = (1.0 - jWeightFactor[i][j]) * sideThreeAdimensional[i] * totalLengthThree +
                                  jWeightFactor[i][j] * sideFourAdimensional[i] * totalLengthFour;
                weightTwo[i][j] = (1.0 - iWeightFactor[i][j]) * sideOneAdimensional[j] * totalLengthOne +
                                  iWeightFactor[i][j] * sideTwoAdimensional[j] * totalLengthTwo;
            }
        }

        for (auto i = 1; i < numMPoints; i++)
        {
            for (auto j = 0; j < numNPoints; j++)
            {
                weightThree[i][j] = weightOne[i][j] - weightOne[i - 1][j];
            }
        }

        for (auto i = 0; i < numMPoints; i++)
        {
            for (auto j = 1; j < numNPoints; j++)
            {
                weightFour[i][j] = weightTwo[i][j] - weightTwo[i][j - 1];
            }
        }

        for (auto i = 1; i < numMPoints; i++)
        {
            for (auto j = 1; j < numNPoints - 1; j++)
            {
                weightOne[i][j] = 0.25 * (weightFour[i][j] + weightFour[i][j + 1] + weightFour[i - 1][j] + weightFour[i - 1][j + 1]) / weightThree[i][j];
            }
        }

        for (auto i = 1; i < numMPoints - 1; i++)
        {
            for (auto j = 1; j < numNPoints; j++)
            {
                weightTwo[i][j] = 0.25 * (weightThree[i][j] + weightThree[i][j - 1] + weightThree[i + 1][j] + weightThree[i + 1][j - 1]) / weightFour[i][j];
            }
        }

        // Iterate several times over
        const int numIterations = 25;
        for (auto iter = 0; iter < numIterations; iter++)
        {
            // re-assign the weights
            for (auto i = 0; i < numMPoints; i++)
            {
                for (auto j = 0; j < numNPoints; j++)
                {
                    weightThree[i][j] = result[i][j].x;
                    weightFour[i][j] = result[i][j].y;
                }
            }

            for (auto i = 1; i < numM; i++)
            {
                for (auto j = 1; j < numN; j++)
                {

                    const double wa = 1.0 / (weightOne[i][j] + weightOne[i + 1][j] + weightTwo[i][j] + weightTwo[i][j + 1]);

                    result[i][j].x = wa * (weightThree[i - 1][j] * weightOne[i][j] + weightThree[i + 1][j] * weightOne[i + 1][j] +
                                           weightThree[i][j - 1] * weightTwo[i][j] + weightThree[i][j + 1] * weightTwo[i][j + 1]);

                    result[i][j].y = wa * (weightFour[i - 1][j] * weightOne[i][j] + weightFour[i + 1][j] * weightOne[i + 1][j] +
                                           weightFour[i][j - 1] * weightTwo[i][j] + weightFour[i][j + 1] * weightTwo[i][j + 1]);
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
            if (first == sizetMissingValue || second == sizetMissingValue)
            {
                continue;
            }

            edgesCenters.emplace_back((nodes[first] + nodes[second]) * 0.5);
        }
        return edgesCenters;
    }

    double LinearInterpolationInTriangle(const Point& interpolationPoint, const std::vector<Point>& polygon, const std::vector<double>& values, const Projection& projection)
    {
        double result = doubleMissingValue;

        const auto a11 = GetDx(polygon[0], polygon[1], projection);
        const auto a21 = GetDy(polygon[0], polygon[1], projection);

        const auto a12 = GetDx(polygon[0], polygon[2], projection);
        const auto a22 = GetDy(polygon[0], polygon[2], projection);

        const auto b1 = GetDx(polygon[0], interpolationPoint, projection);
        const auto b2 = GetDy(polygon[0], interpolationPoint, projection);

        double det = a11 * a22 - a12 * a21;
        if (std::abs(det) < 1e-12)
        {
            return result;
        }

        const double rlam = (a22 * b1 - a12 * b2) / det;
        const double rmhu = (-a21 * b1 + a11 * b2) / det;

        result = values[0] + rlam * (values[1] - values[0]) + rmhu * (values[2] - values[0]);

        return result;
    }

    Point ComputeAverageCoordinate(const std::vector<Point>& points, const Projection& projection)
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

        if (projection == Projection::sphericalAccurate)
        {

            Cartesian3DPoint averagePoint3D{0.0, 0.0, 0.0};
            for (const auto& point : validPoints)
            {
                const Cartesian3DPoint point3D{SphericalToCartesian3D(point)};
                averagePoint3D.x += point3D.x;
                averagePoint3D.y += point3D.y;
                averagePoint3D.z += point3D.z;
            }
            averagePoint3D.x = averagePoint3D.x / validPoints.size();
            averagePoint3D.y = averagePoint3D.y / validPoints.size();
            averagePoint3D.z = averagePoint3D.z / validPoints.size();

            return Cartesian3DToSpherical(averagePoint3D, points[0].x);
        }

        auto result = std::accumulate(validPoints.begin(), validPoints.end(), Point{0.0, 0.0});
        result.x = result.x / validPoints.size();
        result.y = result.y / validPoints.size();
        return result;
    }

} // namespace meshkernel
