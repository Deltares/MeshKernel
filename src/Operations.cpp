#pragma once

#include <cmath>
#include <algorithm>
#include <numeric>
#include "Entities.hpp"
#include "Constants.cpp"

namespace GridGeom
{
    // coordinate reference independent operations
    template<typename T>
    T DotProduct(T& dx1, T& dx2)
    {
        return dx1 * dx2;
    }

    template<typename T, typename... Args>
    T DotProduct(T& dx1, T& dx2, Args&... args)
    {
        return dx1 * dx2 + DotProduct(args...);
    }


    template<typename T>
    bool ResizeVector(int newSize, std::vector<T>& vectorToResize, T fillValue = T())
    {
        const int currentSize = vectorToResize.size();
        if (newSize > currentSize)
        {
            newSize = std::max(newSize, int(currentSize * 1.2));
            vectorToResize.resize(newSize, fillValue);
        }
        return true;
    }

    template<typename T>
    bool AllocateVector(int newSize, std::vector<T>& vectorToResize, int minSize, T fillValue = T())
    {
        const int currentSize = vectorToResize.size();
        if (currentSize < newSize)
        {
            newSize = std::max(minSize, int(5 * newSize));
            vectorToResize.resize(newSize, fillValue);
        }
        return true;
    }

    template <typename T>
    T FindIndex(const std::vector<T>& vec, const T& el)
    {
        T index = 0;
        for (int n = 0; n < vec.size(); n++)
        {
            if (vec[n] == el)
            {
                index = n;
                break;
            }
        }
        return index;
    }

    static int FindIndexes(const std::vector<Point>& vec, const double& separator, std::vector<std::vector<int>>& indexes)
    {
        int pos = 0;

        // set an invalid index
        for (int n = 0; n < indexes.size(); n++)
        {
            indexes[n][0] = -1;
            indexes[n][1] = -1;
        }

        for (int n = 0; n < vec.size(); n++)
        {
            if (vec[n].x != separator && indexes[pos][0] == -1)
            {
                indexes[pos][0] = n;
            }
            else if (vec[n].x == separator && indexes[pos][1] == -1)
            {
                indexes[pos][1] = n - 1;
                pos++;
            }
        }

        if(indexes[pos][1] == -1)
        {
            indexes[pos][1] = vec.size() - 1;
            pos++;
        }

        return pos;
    }

    template <typename T>
    std::vector<int> SortedIndexes(const std::vector<T> &v) 
    {
        std::vector<int> idx(v.size());
        iota(idx.begin(), idx.end(), 0);
        std::stable_sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });
        return std::move(idx);
    }

    //chmike's algorithm
    template< class T >
    void ReorderVector(std::vector<T> &v, std::vector<int> const &order) 
    {
        for (int s = 1, d; s < order.size(); ++s) 
        {
            for (d = order[s]; d < s; d = order[d]);
            if (d == s)
            {
                while (d = order[d], d != s)
                {
                    std::swap(v[s], v[d]);
                } 
            }
        }
    }

    template<typename T>
    bool MakeMonothonic(std::vector<T>& v)
    {

        bool isMonotonic = false;
        int maxIter = 10;
        int iter = 0;
        while (!isMonotonic && iter < maxIter)
        {
            isMonotonic = true;
            maxIter++;
                for (int n = 0; n < v.size(); ++n)
                {
                    if (v[n + 1] - v[n] < 0.0)
                    {
                        isMonotonic = false;
                        break;
                    }
                }
            if (!isMonotonic)
            {
                for (int n = 1; n < v.size() - 1; ++n)
                {
                    v[n] = 0.5 *(v[n - 1] + v[n + 1]);
                }
            }
        }
        return true;
    }


    // algorithm performing the zero's search using the golden section algorithm's
    template <typename F>
    double FindFunctionRootWithGoldenSearch(F func, double min, double max)
    {
        //golden distance factors
        const double c = 0.38196602;
        const double r = 0.61803399;
        const double tolerance = 0.00001;

        double left = min;
        double middle = (min + max)*0.5;
        double right = max;

        double x0 = left;
        double x1 = middle - c * (middle - left);
        double x2 = middle;
        double x3 = right;
        if (std::abs(right - middle)>std::abs(middle - left))
        {
            x1 = middle;
            x2 = middle + c * (middle - left);
        }

        double f1 = func(x1);
        double f2 = func(x2);

        while (std::abs(x3 - x0)>tolerance*std::max(std::abs(x1) + std::abs(x2), 1e-8))
        {
            if (f2<f1)
            {
                x0 = x1;
                x1 = x2;
                x2 = r*x1 + c*x3;

                f1 = f2;
                f2 = func(x2);
            }
            else
            {
                x3 = x2;
                x2 = x1;
                x1 = r*x2 + c*x0;

                f2 = f1;
                f1 = func(x1);
            }
        }
           
        return f1<f2? x1: x2;
    }

    ///sphertocart3D transform 2d spherical to 3d cartesian
    static void SphericalToCartesian(const Point& sphericalPoint, cartesian3DPoint& cartesianPoint)
    {
        cartesianPoint.z = earth_radius * sin(sphericalPoint.y * degrad_hp);
        double rr = earth_radius * cos(sphericalPoint.y * degrad_hp);
        cartesianPoint.x = rr * cos(sphericalPoint.x * degrad_hp);
        cartesianPoint.y = rr * sin(sphericalPoint.x * degrad_hp);
    }

    ///Cart3Dtospher Transform 3d cartesian coordinates to 2d spherical
    static void CartesianToSpherical(const cartesian3DPoint& cartesianPoint, const double referenceLongitude, Point& sphericalPoint)
    {
        double angle = atan2(cartesianPoint.y, cartesianPoint.x) * raddeg_hp;
        sphericalPoint.y = atan2(cartesianPoint.z, sqrt(cartesianPoint.x * cartesianPoint.x + cartesianPoint.y * cartesianPoint.y)) * raddeg_hp;
        sphericalPoint.x = angle + std::lround((referenceLongitude - angle) / 360.0) * 360.0;
    }

    /// IsLeft(): tests if a point is Left|On|Right of an infinite line.
    ///    Input:  three points leftPoint, rightPoint, and point
    ///    Return: >0 for point left of the line through leftPoint and rightPoint
    ///            =0 for point  on the line
    ///            <0 for point  right of the line
    static double IsLeft(const Point& leftPoint, const Point& rightPoint, const Point& point)
    {
        double left = (rightPoint.x - leftPoint.x) * (point.y - leftPoint.y) - (point.x - leftPoint.x) * (rightPoint.y - leftPoint.y);
        return left;
    }


    /// Check if a point is in polygon using the winding number method
    /// polygon: a closed polygon consisting f a vector of numberOfPolygonPoints + 1 in counter clockwise order
    static bool IsPointInPolygon(const Point& point, const std::vector<Point>& polygon, const int numberOfPolygonPoints)
    {
        if (numberOfPolygonPoints <= 0)
        {
            return true;
        }

        int windingNumber = 0;
        for (int n = 0; n < numberOfPolygonPoints; n++)
        {
            if (polygon[n].y <= point.y) // an upward crossing
            {
                if (polygon[n + 1].y > point.y)
                {
                    if (IsLeft(polygon[n], polygon[n + 1], point) > 0.0)
                    {
                        ++windingNumber; // have  a valid up intersect
                    }
                }
            }
            else
            {
                if (polygon[n + 1].y <= point.y) // a downward crossing
                {
                    if (IsLeft(polygon[n], polygon[n + 1], point) < 0.0)
                    {
                        --windingNumber; // have  a valid down intersect
                    }
                }
            }
        }
        return windingNumber == 0 ? false : true;
    }

    // dbpinpol_optinside_perpol
    static bool IsPointInPolygons(const Point& point, const std::vector<Point>& polygon, const int numberOfPolygonPoints)
    {
        if (numberOfPolygonPoints <= 0)
        {
            return true;
        }

        // Split input polygons in parts 
        std::vector<std::vector<int>> indexes(numberOfPolygonPoints, std::vector<int>(2));
        std::vector<Point> subpolygon(numberOfPolygonPoints);
        int numPolygons = FindIndexes(polygon, doubleMissingValue, indexes);

        bool inPolygon = false;
        auto compareX = [](const auto&p1, const auto&p2) { return p1->x > p2->x; };
        auto compareY = [](const auto&p1, const auto&p2) { return p1->y > p2->y; };
        for (int p = 0; p <= numPolygons; p++) 
        {
            // Calculate the bounding box
            auto startPolygon = polygon.begin() + indexes[p][0];
            auto endPolygon = polygon.begin() + indexes[p][1];
            const int numberOfSubPolygonPoints = indexes[p][1] - indexes[p][0] + 2;
            
            for (int pp = 0; pp < numberOfSubPolygonPoints; pp++) 
            {
                subpolygon[pp] = polygon[indexes[p][0] + pp];
            }
            subpolygon[numberOfSubPolygonPoints] = polygon[indexes[p][0]];
            
            // create bounding box before actual checking
            double XMin = std::min(startPolygon, endPolygon, compareX)->x;
            double XMax = std::max(startPolygon, endPolygon, compareX)->x;
            double YMin = std::min(startPolygon, endPolygon, compareY)->y;
            double YMax = std::max(startPolygon, endPolygon, compareY)->y;

            if ((point.x >= XMin && point.x <= XMax) && (point.y >= YMin && point.y <= YMax)) 
            {
                 inPolygon = IsPointInPolygon(point, subpolygon, numberOfSubPolygonPoints);
            }

            if (inPolygon) 
            {
                return true;
            }
        }

        return false;
    }



    static double getDx(const Point& firstPoint, const Point& secondPoint, const Projections& projection)
    {
        if (projection == Projections::cartesian)
        {
            return secondPoint.x - firstPoint.x;
        }
        if (projection == Projections::spherical)
        {
            double  firstPointYDiff = abs(abs(firstPoint.y) - 90.0);
            double  secondPointYDiff = abs(abs(secondPoint.y) - 90.0);
            if (firstPointYDiff <= dtol_pole && secondPointYDiff > dtol_pole || firstPointYDiff > dtol_pole && secondPointYDiff <= dtol_pole)
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
            double firstPointY = firstPoint.y * degrad_hp;
            double secondPointY = secondPoint.y * degrad_hp;
            double cosPhi = cos(0.5 * (firstPointY + secondPointY));
            double dx = earth_radius * cosPhi * (secondPointX - firstPointX);
            return dx;
        }
        return doubleMissingValue;
    }

    static double getDy(const Point& firstPoint, const Point& secondPoint, const Projections& projection)
    {
        if (projection == Projections::cartesian)
        {
            return secondPoint.y - firstPoint.y;
        }
        if (projection == Projections::spherical)
        {
            double firstPointY = firstPoint.y * degrad_hp;
            double secondPointY = secondPoint.y * degrad_hp;
            double dy = earth_radius * (secondPointY - firstPointY);
            return dy;
        }
        return doubleMissingValue;
    }

    //out product of two segments
    static double outerProductTwoSegments(const Point& firstPointFirstSegment, const Point& secondPointFirstSegment, const Point& firstPointSecondSegment, const Point& secondPointSecondSegment, const Projections& projection)
    {
        if (projection == Projections::cartesian)
        {
            double dx1 = getDx(firstPointFirstSegment, secondPointFirstSegment, projection);
            double dx2 = getDx(firstPointSecondSegment, secondPointSecondSegment, projection);

            double dy1 = getDy(firstPointFirstSegment, secondPointFirstSegment, projection);
            double dy2 = getDy(firstPointSecondSegment, secondPointSecondSegment, projection);

            return dx1 * dy2 - dy1 * dx2;
        }
        if (projection == Projections::spherical)
        {
            //TODO: IMPLEMENTATION IS MISSING
        }
        return doubleMissingValue;
    }

    static void MiddlePoint(const Point& firstPoint, const Point& secondPoint, Point& result, const Projections& projection)
    {
        if (projection == Projections::cartesian)
        {
            cartesian3DPoint firstPointCartesianCoordinates;
            SphericalToCartesian(firstPoint, firstPointCartesianCoordinates);
            cartesian3DPoint secondPointCartesianCoordinates;
            SphericalToCartesian(secondPoint, secondPointCartesianCoordinates);

            cartesian3DPoint middleCartesianPointCoordinate;
            middleCartesianPointCoordinate.x = 0.5*(firstPointCartesianCoordinates.x + secondPointCartesianCoordinates.x);
            middleCartesianPointCoordinate.y = 0.5*(firstPointCartesianCoordinates.y + secondPointCartesianCoordinates.y);
            double referenceLongitude = std::max(firstPoint.x, secondPoint.x);

            CartesianToSpherical(middleCartesianPointCoordinate, referenceLongitude, result);
        }
        if (projection == Projections::spherical)
        {
            result = { (firstPoint.x + secondPoint.x) / 2.0, (firstPoint.y + secondPoint.y) / 2.0 };
        }
    }

    ///normaloutchk, Computes the normal vector to a line 1-2, which is *outward* w.r.t. an 'inside' point 3. Similar to normalout, except that the normal
    ///vector may be flipped based on the 'inside' point.
    static void NormalVector(const Point& firstPoint, const Point& secondPoint, const Point& insidePoint, Point& result, const Projections& projection)
    {
        if (projection == Projections::cartesian)
        {
            double dx = getDx(firstPoint, secondPoint, projection);
            double dy = getDy(firstPoint, secondPoint, projection);
            const double squaredDistance = dx * dx + dy * dy;
            if (squaredDistance != 0.0)
            {
                const double distance = sqrt(squaredDistance);
                result.x = dy / distance;
                result.y = -dx / distance;
            }
        }
        if (projection == Projections::spherical)
        {
            cartesian3DPoint firstPointCartesianCoordinates;
            cartesian3DPoint secondPointCartesianCoordinates;
            SphericalToCartesian(firstPoint, firstPointCartesianCoordinates);
            SphericalToCartesian(secondPoint, secondPointCartesianCoordinates);

            double lambda = insidePoint.x * degrad_hp;
            double phi = insidePoint.y * degrad_hp;
            double elambda[3] = { -sin(lambda), cos(lambda), 0.0 };
            double ephi[3] = { -sin(lambda), cos(lambda), 0.0 };

            double dx = (secondPointCartesianCoordinates.x - firstPointCartesianCoordinates.x) * elambda[0] +
                (secondPointCartesianCoordinates.y - firstPointCartesianCoordinates.y) * elambda[1] +
                (secondPointCartesianCoordinates.z - firstPointCartesianCoordinates.z) * elambda[2];

            double dy = (secondPointCartesianCoordinates.x - firstPointCartesianCoordinates.x) * ephi[0] +
                (secondPointCartesianCoordinates.y - firstPointCartesianCoordinates.y) * ephi[1] +
                (secondPointCartesianCoordinates.z - firstPointCartesianCoordinates.z) * ephi[2];

            double squaredDistance = dx * dx + dy * dy;
            if (squaredDistance != 0.0)
            {
                const double distance = sqrt(squaredDistance);
                result.x = dx / distance;
                result.y = dy / distance;
            }
        }
    }

    ///normaloutchk
    static void NormalVectorInside(const Point& firstPoint, const Point& secondPoint, const Point& insidePoint, Point& result, bool& flippedNormal, const Projections& projection)
    {
        if (projection == Projections::cartesian)
        {
            NormalVector(firstPoint, secondPoint, insidePoint, result, projection);
            flippedNormal = false;
            Point thirdPoint{ firstPoint.x + result.x, firstPoint.y + result.y };

            if (outerProductTwoSegments(firstPoint, thirdPoint, firstPoint, secondPoint, projection) * outerProductTwoSegments(firstPoint, insidePoint, firstPoint, secondPoint, projection) > 0.0)
            {
                result.x = -result.x;
                result.y = -result.y;
                flippedNormal = true;
            }
            else
            {
                flippedNormal = false;
            }
        }
        if (projection == Projections::spherical)
        {
            //if (JSFERIC.eq.1 . and .jasfer3D.eq.0) xn = xn * cos(dg2rd * 0.5d0 * (y0 + y1)) !normal vector needs to be in Cartesian coordinates
        }
    }

    ///normalout
    static void NormalVectorOutside(const Point& firstPoint, const Point& secondPoint, Point& result, const Projections& projection)
    {
        if (projection == Projections::cartesian)
        {
            double dx = getDx(firstPoint, secondPoint, projection);
            double dy = getDy(firstPoint, secondPoint, projection);
            const double squaredDistance = dx * dx + dy * dy;
            if (squaredDistance != 0.0)
            {
                const double distance = sqrt(squaredDistance);
                result.x = dy / distance;
                result.y = -dx / distance;
            }
        }
        if (projection == Projections::spherical)
        {
            Point middlePoint;
            MiddlePoint(firstPoint, secondPoint, middlePoint, projection);

            cartesian3DPoint firstPointCartesianCoordinates;
            SphericalToCartesian(firstPoint, firstPointCartesianCoordinates);
            cartesian3DPoint secondPointCartesianCoordinates;
            SphericalToCartesian(secondPoint, secondPointCartesianCoordinates);

            //compute the base vectors at middle point
            double lambda = middlePoint.x * degrad_hp;
            double phi = middlePoint.y * degrad_hp;
            double elambda[3] = { -sin(lambda), cos(lambda), 0.0 };
            double ephi[3] = { -sin(phi)*cos(lambda), -sin(phi)*sin(lambda), cos(phi) };

            // project vector in local base
            double dx = (secondPointCartesianCoordinates.x - firstPointCartesianCoordinates.x) * elambda[0] +
                (secondPointCartesianCoordinates.y - firstPointCartesianCoordinates.y) * elambda[1] +
                (secondPointCartesianCoordinates.z - firstPointCartesianCoordinates.z) * elambda[2];

            double dy = (secondPointCartesianCoordinates.x - firstPointCartesianCoordinates.x) * ephi[0] +
                (secondPointCartesianCoordinates.y - firstPointCartesianCoordinates.y) * ephi[1] +
                (secondPointCartesianCoordinates.z - firstPointCartesianCoordinates.z) * ephi[2];

            result.x = result.x / cos(degrad_hp*0.5*(firstPoint.y + secondPoint.y));
            result.y = middlePoint.y;
        }

    }

    static void add(Point& point, const Point& normal, const double increment, const Projections& projection)
    {
        if (projection == Projections::cartesian)
        {
            point.x = point.x + normal.x * increment;
            point.y = point.y + normal.y * increment;
        }
        if (projection == Projections::spherical)
        {
            double convertedIncrement = raddeg_hp * increment / earth_radius;
            double xf = 1.0 / cos(degrad_hp * point.y);
            point.x = point.x + normal.x * convertedIncrement * xf;
            point.y = point.y + normal.y * convertedIncrement;
        }
    }

    static void referencePoint(std::vector<Point>& polygon, const int numPoints, double& minX, double& minY, const Projections& projection)
    {
        if (projection == Projections::cartesian)
        {
            minX = std::numeric_limits<double>::max();
            minY = std::numeric_limits<double>::max();
            for (int i = 0; i < numPoints; i++)
            {
                if (polygon[i].x < minX)
                {
                    minX = polygon[i].x;
                }
                if (abs(polygon[i].y) < abs(minY))
                {
                    minY = polygon[i].y;
                }
            }
        }
        if (projection == Projections::spherical)
        {
            minX = std::numeric_limits<double>::max();
            minY = std::numeric_limits<double>::max();
            double xmax = std::numeric_limits<double>::min();
            for (int i = 0; i < numPoints; i++)
            {
                if (polygon[i].x < minX)
                {
                    minX = polygon[i].x;
                }
                if (abs(polygon[i].y) < abs(minY))
                {
                    minY = polygon[i].y;
                }
                if (polygon[i].x > xmax)
                {
                    xmax = polygon[i].x;
                }
            }

            if (xmax - minX > 180.0)
            {
                double deltaX = xmax - 180.0;
                for (int i = 0; i < numPoints; i++)
                {
                    if (polygon[i].x < deltaX)
                    {
                        polygon[i].x = polygon[i].x + 360.0;
                    }
                }
                minX = minX + 360.0;
            }
            //TODO: check result
            minX = std::min_element(polygon.begin(), polygon.end(), [](const Point& p1, const Point& p2) { return p1.x < p2.x; })->x;
        }
    }

    //dbdistance
    static double Distance(const Point& firstPoint, const Point& secondPoint, const Projections& projection)
    {

        if (firstPoint.x == doubleMissingValue || firstPoint.y == doubleMissingValue ||
            secondPoint.x == doubleMissingValue || secondPoint.y == doubleMissingValue)
            return 0.0;

        if (projection == Projections::cartesian)
        {
            double dx = getDx(firstPoint, secondPoint, projection);
            double dy = getDy(firstPoint, secondPoint, projection);
            const double squaredDistance = dx * dx + dy * dy;
            double distance = 0.0;
            if (squaredDistance != 0.0)
            {
                distance = sqrt(squaredDistance);
            }
            return distance;
        }
        if (projection == Projections::spherical)
        {
        }
        return doubleMissingValue;
    }

    // dLINEDIS3
    // Computes the perpendicular distance from point to a line firstNode - secondNode.
    // normalPoint: coordinates of the projected point from point onto the line
    static double DistanceFromLine(const Point& point, const Point& firstNode, const Point& secondNode, Point& normalPoint, double& ratio, const Projections& projection)
    {
        if (projection == Projections::cartesian)
        {
            double dis = 0.0;
            double r2 = Distance(secondNode, firstNode, projection);
            if (r2 != 0.0)
            {
                ratio = (getDx(firstNode, point, projection) * getDx(firstNode, secondNode, projection) + getDy(firstNode, point, projection) * getDy(firstNode, secondNode, projection)) / (r2 * r2);
                double correctedRatio = std::max(std::min(1.0, ratio), 0.0);
                normalPoint.x = firstNode.x + correctedRatio * (secondNode.x - firstNode.x);
                normalPoint.y = firstNode.y + correctedRatio * (secondNode.y - firstNode.y);
                dis = Distance(point, normalPoint, projection);
            }
            return dis;
        }
        if (projection == Projections::spherical)
        {
            //TODO: implement me
        }
        return -1.0;
    }

    //inner product of two segments
    static double innerProductTwoSegments(const Point& firstPointFirstSegment, const Point& secondPointFirstSegment, const Point& firstPointSecondSegment, const Point& secondPointSecondSegment, const Projections& projection)
    {
        if (projection == Projections::cartesian)
        {
            double dx1 = getDx(firstPointFirstSegment, secondPointFirstSegment, projection);
            double dx2 = getDx(firstPointSecondSegment, secondPointSecondSegment, projection);

            double dy1 = getDy(firstPointFirstSegment, secondPointFirstSegment, projection);
            double dy2 = getDy(firstPointSecondSegment, secondPointSecondSegment, projection);

            return dx1 * dx2 + dy1 * dy2;
        }
        if (projection == Projections::spherical)
        {
            cartesian3DPoint firstPointFirstSegment3D;
            cartesian3DPoint secondPointFirstSegment3D;
            cartesian3DPoint firstPointSecondSegment3D;
            cartesian3DPoint secondPointSecondSegment3D;

            SphericalToCartesian(firstPointFirstSegment, firstPointFirstSegment3D);
            SphericalToCartesian(secondPointFirstSegment, secondPointFirstSegment3D);
            SphericalToCartesian(firstPointSecondSegment, firstPointSecondSegment3D);
            SphericalToCartesian(secondPointSecondSegment, secondPointSecondSegment3D);

            double dx1 = secondPointFirstSegment3D.x - firstPointFirstSegment3D.x;
            double dy1 = secondPointFirstSegment3D.y - firstPointFirstSegment3D.y;
            double dz1 = secondPointFirstSegment3D.z - firstPointFirstSegment3D.z;

            double dx2 = secondPointSecondSegment3D.x - firstPointSecondSegment3D.x;
            double dy2 = secondPointSecondSegment3D.y - firstPointSecondSegment3D.y;
            double dz2 = secondPointSecondSegment3D.z - firstPointSecondSegment3D.z;

            return DotProduct(dx1, dx2, dy1, dy2, dz1, dz2);
        }

        return doubleMissingValue;
    }

    // dcosphi
    static double normalizedInnerProductTwoSegments(const Point& firstPointFirstSegment, const Point& secondPointFirstSegment, const Point& firstPointSecondSegment, const Point& secondPointSecondSegment, const Projections& projection)
    {
        if (projection == Projections::cartesian)
        {
            double dx1 = getDx(firstPointFirstSegment, secondPointFirstSegment, projection);
            double dx2 = getDx(firstPointSecondSegment, secondPointSecondSegment, projection);

            double dy1 = getDy(firstPointFirstSegment, secondPointFirstSegment, projection);
            double dy2 = getDy(firstPointSecondSegment, secondPointSecondSegment, projection);

            double r1 = dx1 * dx1 + dy1 * dy1;
            double r2 = dx2 * dx2 + dy2 * dy2;

            double cosphi;
            if (r1 == 0.0 || r2 == 0.0)
            {
                cosphi = doubleMissingValue;
            }

            cosphi = (dx1 * dx2 + dy1 * dy2) / std::sqrt(r1 * r2);
            cosphi = std::max(std::min(cosphi, 1.0), -1.0);
            return cosphi;
        }
        if (projection == Projections::spherical)
        {
            //TODO: IMPLEMENTATION IS MISSING
        }
        return doubleMissingValue;
    }

    static bool orthogonalizationComputeLocalCoordinates(const std::vector<std::size_t>& m_nodesNumEdges, const std::vector<std::size_t>& numConnectedNodes, std::vector<int>& localCoordinates, const Projections& projection)
    {
        if (projection == Projections::cartesian)
        {
            //do nothing
        }
        if (projection == Projections::spherical)
        {
            localCoordinates.resize(m_nodesNumEdges.size(), 0);
            localCoordinates[0] = 1;
            for (int i = 0; i < m_nodesNumEdges.size(); i++)
            {
                localCoordinates[i + 1] = localCoordinates[i] + std::max(m_nodesNumEdges[i] + 1, numConnectedNodes[i]);
            }
        }
        return true;
    }

    static bool orthogonalizationComputeJacobian(const int currentNode, const std::vector<double>& Jxi, const std::vector<double>& Jeta, const std::vector<std::size_t>& connectedNodes, const int numNodes, const std::vector<Point>& nodes, std::vector<double>& J, const Projections& projection)
    {
        if (projection == Projections::cartesian)
        {
            J[0] = 0.0;
            J[1] = 0.0;
            J[2] = 0.0;
            J[3] = 0.0;
            for (int i = 0; i < numNodes; i++)
            {
                J[0] += Jxi[i] * nodes[connectedNodes[i]].x;
                J[1] += Jxi[i] * nodes[connectedNodes[i]].y;
                J[2] += Jeta[i] * nodes[connectedNodes[i]].x;
                J[3] += Jeta[i] * nodes[connectedNodes[i]].y;
            }
        }
        if (projection == Projections::spherical)
        {
            double factor = std::cos(nodes[currentNode].y) * degrad_hp;
            J[0] = 0.0;
            J[1] = 0.0;
            J[2] = 0.0;
            J[3] = 0.0;
            for (int i = 0; i < numNodes; i++)
            {
                J[0] += Jxi[i] * nodes[connectedNodes[i]].x;
                J[1] += Jxi[i] * nodes[connectedNodes[i]].y;
                J[2] += Jeta[i] * nodes[connectedNodes[i]].x;
                J[3] += Jeta[i] * nodes[connectedNodes[i]].y;
            }
        }
        return true;
    }

    static bool orthogonalizationComputeDeltasDxDy(int firstNode, int secondNode, double wwx, double wwy, const std::vector<Point>& nodes, double& dx0, double& dy0, const Projections& projection)
    {
        if (projection == Projections::cartesian)
        {
            dx0 += wwx * (nodes[firstNode].x - nodes[secondNode].x);
            dy0 += wwy * (nodes[firstNode].y - nodes[secondNode].y);
        }
        if (projection == Projections::spherical)
        {
            double wwxTransformed = wwx * earth_radius * degrad_hp;
            double wwyTransformed = wwy * earth_radius * degrad_hp;

            dx0 += dx0 + wwxTransformed * (nodes[firstNode].x - nodes[secondNode].x);
            dy0 += dy0 + wwyTransformed * (nodes[firstNode].y - nodes[secondNode].y);
        }
        return true;
    }

    static bool orthogonalizationComputeIncrements(double wwx, double wwy, double* increments, const Projections& projection)
    {
        if (projection == Projections::cartesian)
        {
            increments[0] += wwx;
            increments[1] += wwy;

        }
        if (projection == Projections::spherical)
        {
            double wwxTransformed = wwx * earth_radius * degrad_hp;
            double wwyTransformed = wwy * earth_radius * degrad_hp;

            increments[0] += wwxTransformed;
            increments[1] += wwyTransformed;
        }
        return true;
    }

    static bool orthogonalizationComputeDeltas(int firstNode, int secondNode, double wwx, double wwy, const std::vector<Point>& nodes, double& dx0, double& dy0, double* increments, const Projections& projection)
    {
        if (projection == Projections::cartesian)
        {
            increments[0] += wwx;
            increments[1] += wwy;

            dx0 = dx0 + wwx * (nodes[firstNode].x - nodes[secondNode].x);
            dy0 = dy0 + wwy * (nodes[firstNode].y - nodes[secondNode].y);
        }
        if (projection == Projections::spherical)
        {
            double wwxTransformed = wwx * earth_radius * degrad_hp;
            double wwyTransformed = wwy * earth_radius * degrad_hp;

            increments[0] += wwxTransformed;
            increments[1] += wwyTransformed;

            dx0 = dx0 + wwxTransformed * (nodes[firstNode].x - nodes[secondNode].x);
            dy0 = dy0 + wwxTransformed * (nodes[firstNode].y - nodes[secondNode].y);
        }
        return true;
    }

    static bool orthogonalizationComputeCoordinates(double dx0, double dy0, const Point& point, Point& updatedPoint, const Projections& projection)
    {
        if (projection == Projections::cartesian)
        {
            double x0 = point.x + dx0;
            double y0 = point.y + dy0;
            static constexpr double relaxationFactorCoordinates = 1.0 - relaxationFactorOrthogonalizationUpdate;

            updatedPoint.x = relaxationFactorOrthogonalizationUpdate * x0 + relaxationFactorCoordinates * point.x;
            updatedPoint.y = relaxationFactorOrthogonalizationUpdate * y0 + relaxationFactorCoordinates * point.y;
        }
        if (projection == Projections::spherical)
        {
            //TODO: implement
            //if (jsferic.eq.1 . and .jasfer3D.eq.1) then
            //    dumx(1) = relaxin * Dx0
            //    dumy(1) = relaxin * Dy0
            //    call loc2spher(xk(k), yk(k), 1, dumx, dumy, xk1(k), yk1(k))
            //else
        }
        return true;
    }

    static bool circumcenterOfTriangle(const Point& p1, const Point& p2, const Point& p3, Point& circumcenter, const Projections& projection)
    {
        if (projection == Projections::cartesian)
        {
            double dx2 = getDx(p1, p2, projection);
            double dy2 = getDy(p1, p2, projection);

            double dx3 = getDx(p1, p3, projection);
            double dy3 = getDy(p1, p3, projection);

            double den = dy2 * dx3 - dy3 * dx2;
            double z = 0.0;
            if (den != 0.0)
            {
                z = (dx2 * (dx2 - dx3) + dy2 * (dy2 - dy3)) / den;
            }

            circumcenter.x = p1.x + 0.5 * (dx3 - z * dy3);
            circumcenter.y = p1.y + 0.5 * (dy3 + z * dx3);
        }
        if (projection == Projections::spherical)
        {
            double dx2 = getDx(p1, p2, projection);
            double dy2 = getDy(p1, p2, projection);

            double dx3 = getDx(p1, p3, projection);
            double dy3 = getDy(p1, p3, projection);

            double den = dy2 * dx3 - dy3 * dx2;
            double z = 0.0;
            if (den > 1e-16)
            {
                z = (dx2 * (dx2 - dx3) + dy2 * (dy2 - dy3)) / den;
            }

            //TODO circumcenter3 FINISH
            //phi = (y(1) + y(2) + y(3)) / 3d0
            //    xf = 1d0 / dcos(degrad_hp * phi)
            //    xz = x(1) + xf * 0.5d0 * (dx3 - z * dy3) * raddeg_hp / earth_radius
            //    yz = y(1) + 0.5d0 * (dy3 + z * dx3) * raddeg_hp / earth_radius
        }
        return true;
    }

    //CROSS
    static bool AreLinesCrossing(const Point& firstSegmentFistPoint, 
        const Point& firstSegmentSecondPoint, 
        const Point& secondSegmentFistPoint, 
        const Point& secondSegmentSecondPoint, 
        bool adimensional, 
        Point& intersection, 
        double& crossProduct,
        double& firstRatio,
        double& secondRatio,
        const Projections& projection)
    {
        bool isCrossing = false;
        if (projection == Projections::cartesian)
        {
            firstRatio = doubleMissingValue;
            secondRatio = doubleMissingValue;
            double x21 = getDx(firstSegmentFistPoint, firstSegmentSecondPoint, projection);
            double y21 = getDy(firstSegmentFistPoint, firstSegmentSecondPoint, projection);

            double x43 = getDx(secondSegmentFistPoint, secondSegmentSecondPoint, projection);
            double y43 = getDy(secondSegmentFistPoint, secondSegmentSecondPoint, projection);

            double x31 = getDx(firstSegmentFistPoint, secondSegmentFistPoint, projection);
            double y31 = getDy(firstSegmentFistPoint, secondSegmentFistPoint, projection);

            double det = x43 * y21 - y43 * x21;

            std::vector<double> values{ x21, y21, x43, y43 };
            double eps = std::max(0.00001 * (*std::max_element(values.begin(), values.end())), std::numeric_limits<double>::denorm_min());

            if (std::abs(det) < eps)
            {
                return isCrossing;
            }

            secondRatio = (y31 * x21 - x31 * y21) / det;
            firstRatio = (y31 * x43 - x31 * y43) / det;
            if (firstRatio >= 0.0 && firstRatio <= 1.0 && secondRatio >= 0.0 && secondRatio <= 1.0)
            {
                isCrossing = true;
            }
            intersection.x = firstSegmentFistPoint.x + firstRatio * (firstSegmentSecondPoint.x - firstSegmentFistPoint.x);
            intersection.y = firstSegmentFistPoint.y + firstRatio * (firstSegmentSecondPoint.y - firstSegmentFistPoint.y);
            crossProduct = -det;
            if (adimensional)
            {
                crossProduct = -det / (std::sqrt(x21 * x21 + y21 * y21) * std::sqrt(x43 * x43 + y43 * y43) + 1e-8);
            }
        }

        if (projection == Projections::spherical)
        {
            double x21 = getDx(firstSegmentFistPoint, firstSegmentSecondPoint, projection);
            double y21 = getDy(firstSegmentFistPoint, firstSegmentSecondPoint, projection);

            double x43 = getDx(secondSegmentFistPoint, secondSegmentSecondPoint, projection);
            double y43 = getDy(secondSegmentFistPoint, secondSegmentSecondPoint, projection);

            double x31 = getDx(firstSegmentFistPoint, secondSegmentFistPoint, projection);
            double y31 = getDy(firstSegmentFistPoint, secondSegmentFistPoint, projection);

            double det = x43 * y21 - y43 * x21;

            std::vector<double> values{ x21, y21, x43, y43 };
            double eps = std::max(0.00001 * (*std::max_element(values.begin(), values.end())), std::numeric_limits<double>::denorm_min());

            if (det < eps)
            {
                return isCrossing;
            }

            double sm = (y31 * y21 - x31 * y21) / det;
            double sl = (y31 * x43 - x31 * y43) / det;
            if (sm >= 0.0 && sm <= 1.0 && sl >= 0.0 && sl <= 1.0)
            {
                isCrossing = true;
            }
            intersection.x = firstSegmentFistPoint.x + sl * (firstSegmentSecondPoint.x - firstSegmentFistPoint.x);
            intersection.y = firstSegmentFistPoint.y + sl * (firstSegmentSecondPoint.x - firstSegmentFistPoint.y);
            crossProduct = -det;
            if (adimensional)
            {
                crossProduct = -det / (std::sqrt(x21 * x21 + y21 * y21) * std::sqrt(x43 * x43 + y43 * y43) + 1e-8);
            }
        }
        return isCrossing;
    }

    //faceAreaAndCenterOfMass: for cartesian, spherical point and spherical3dPoint
    static bool faceAreaAndCenterOfMass(std::vector<Point>& polygon, const int numberOfPolygonPoints, double& area, Point& centerOfMass, const Projections& projection)
    {
        if (numberOfPolygonPoints < 1)
        {
            return false;
        }

        double minX = std::numeric_limits<double>::max();
        double minY = std::numeric_limits<double>::max();
        referencePoint(polygon, numberOfPolygonPoints, minX, minY, projection);

        Point reference{ minX, minY };
        area = 0.0;
        double xCenterOfMass = 0.0;
        double yCenterOfMass = 0.0;
        for (int p = 0; p < numberOfPolygonPoints; p++)
        {
            double dx0 = getDx(reference, polygon[p], projection);
            double dy0 = getDy(reference, polygon[p], projection);
            double dx1 = getDx(reference, polygon[p + 1], projection);
            double dy1 = getDy(reference, polygon[p + 1], projection);

            double xc = 0.5 * (dx0 + dx1);
            double yc = 0.5 * (dy0 + dy1);

            dx0 = getDx(polygon[p], polygon[p + 1], projection);
            dy0 = getDy(polygon[p], polygon[p + 1], projection);
            double dsx = dy0;
            double dsy = -dx0;
            double xds = xc * dsx + yc * dsy;
            area = area + 0.5 * xds;

            xCenterOfMass = xCenterOfMass + xds * xc;
            yCenterOfMass = yCenterOfMass + xds * yc;
        }

        double fac = 1.0 / (numberOfPolygonPoints * area);
        xCenterOfMass = fac * xCenterOfMass;
        yCenterOfMass = fac * yCenterOfMass;

        //if constexpr (operationType == sphericalOperations)
        //{
        //    yCenterOfMass = yCenterOfMass / (earth_radius * degrad_hp);
        //    xCenterOfMass = xCenterOfMass / (earth_radius * degrad_hp * cos((yCenterOfMass + minY) * degrad_hp));
        //}

        centerOfMass.x = xCenterOfMass + minX;
        centerOfMass.y = yCenterOfMass + minY;

        area = std::abs(area);

        return true;
    }




}