#ifndef OPERTIONS_CPP
#define OPERTIONS_CPP
#include <cmath>
#include "Entities.hpp"
#include "Constants.cpp"

#include <boost/geometry/geometries/segment.hpp> 
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>
#include <boost/geometry/geometries/register/ring.hpp>

#include <boost/geometry.hpp>
// register node so we can use boost geometry algorithms
BOOST_GEOMETRY_REGISTER_POINT_2D(GridGeom::cartesianPoint, double, boost::geometry::cs::cartesian, x, y);
BOOST_GEOMETRY_REGISTER_RING(std::vector<GridGeom::cartesianPoint>);

namespace GridGeom
{

    // functions that depends of the point type
    template<typename T>
    struct Operations;

    // coordinate reference indipendent operations
    template<typename T>
    static double dotProduct(const T& dx1, const T& dx2)
    {
        return dx1 * dx2;
    }

    template<typename T, typename... Args>
    static T dotProduct(T dx1, T dx2, Args... args) 
    {
        return dx1 * dx2 + dotProduct(args...);
    }

    // transform 2d spherical to 3d cartesian
    static void sphericalToCartesian(const sphericalPoint& sphericalPoint, cartesian3DPoint& cartesianPoint)
    {
        cartesianPoint.z = earth_radius * sin(sphericalPoint.y * degrad_hp);
        double rr = earth_radius * cos(sphericalPoint.y * degrad_hp);
        cartesianPoint.x = rr * cos(sphericalPoint.x * degrad_hp);
        cartesianPoint.y = rr * sin(sphericalPoint.x * degrad_hp);
    }

    //  transform 3d cartesian coordinates to 2d spherical
    static void cartesianToSpherical(const cartesian3DPoint& cartesianPoint, const double referenceLongitude, sphericalPoint& sphericalPoint)
    {
        double angle = atan2(cartesianPoint.y, cartesianPoint.x) * raddeg_hp;
        sphericalPoint.y = atan2(cartesianPoint.z, sqrt(cartesianPoint.x * cartesianPoint.x + cartesianPoint.y * cartesianPoint.y)) * raddeg_hp;
        sphericalPoint.x = angle + std::lround((referenceLongitude - angle) / 360.0) * 360.0;
    }

    // isLeft(): tests if a point is Left|On|Right of an infinite line.
    //    Input:  three points leftPoint, rightPoint, and point
    //    Return: >0 for point left of the line through leftPoint and rightPoint
    //            =0 for point  on the line
    //            <0 for point  right of the line
    static double isLeft(const cartesianPoint& leftPoint, const cartesianPoint& rightPoint, const cartesianPoint& point)
    {
        double left = (rightPoint.x - leftPoint.x) * (point.y - leftPoint.y) - (point.x - leftPoint.x) * (rightPoint.y - leftPoint.y);
        return left;
    }

    // check if a point is in polygon using the winding number method
    // polygon: vector of points in counter clockwise order
    static bool pointInPolygon(const cartesianPoint& point, const std::vector<cartesianPoint>& polygon, const int numberOfPolygonPoints)
    {
        int windingNumber = 0;
        for (int n = 0; n < numberOfPolygonPoints; n++)
        {
            if (polygon[n].y <= point.y) // an upward crossing
            {
                if (polygon[n + 1].y > point.y)
                {
                    if (isLeft(polygon[n], polygon[n + 1], point) > 0.0)
                    {
                        ++windingNumber; // have  a valid up intersect
                    }
                }
            }
            else
            {
                if (polygon[n + 1].y <= point.y) // a downward crossing
                {
                    if (isLeft(polygon[n], polygon[n + 1], point) < 0.0)
                    {
                        --windingNumber; // have  a valid down intersect
                    }
                }
            }
        }
        return windingNumber == 0 ? false : true;
    }


    //faceAreaAndCenterOfMass: for cartesian, spherical point and spherical3dPoint

    template<typename Point>
    static bool faceAreaAndCenterOfMass(const std::vector<Point>& polygon, const int numberOfPolygonPoints, double& area, Point& centerOfMass)
    {
        //double area = boost::geometry::area(pol);
        //return area;
        if (numberOfPolygonPoints < 1)
        {
            return false;
        }
        
        double minX = std::numeric_limits<double>::max();
        double minY = std::numeric_limits<double>::max();
        Operations<Point>::referencePoint(polygon, minX, minY);

        Point reference{ minX, minY };
        area = 0.0;
        double xCenterOfMass = 0.0;
        double yCenterOfMass = 0.0;
        for (int p = 0; p < numberOfPolygonPoints; p++)
        {
            double dx0 = Operations<Point>::getDx(reference, polygon[p]);
            double dy0 = Operations<Point>::getDy(reference, polygon[p]);
            double dx1 = Operations<Point>::getDx(reference, polygon[p + 1]);
            double dy1 = Operations<Point>::getDy(reference, polygon[p + 1]);

            double xc = 0.5 * (dx0 + dx1);
            double yc = 0.5 * (dy0 + dy1);

            dx0 = Operations<Point>::getDx(polygon[p], polygon[p + 1]);
            dy0 = Operations<Point>::getDy(polygon[p], polygon[p + 1]);
            double dsx = dy0;
            double dsy = -dx0;
            double xds = xc * dsx + yc * dsy;
            area = area + 0.5 * xds;

            xCenterOfMass = xCenterOfMass + xds * xc;
            yCenterOfMass = yCenterOfMass + xds * yc;
        }

        double fac = 1.0 / (3.0 * area);
        xCenterOfMass = fac * xCenterOfMass;
        yCenterOfMass = fac * yCenterOfMass;

        if constexpr (std::is_same<Point, cartesianPoint>::value)
        {
            yCenterOfMass = yCenterOfMass / (earth_radius * degrad_hp);
            xCenterOfMass = xCenterOfMass / (earth_radius * degrad_hp * cos((yCenterOfMass + minY) * degrad_hp));
        }

        centerOfMass.x = xCenterOfMass + minX;
        centerOfMass.y = yCenterOfMass + minY;

        return true;
    }

    static bool lineCrossing(const cartesianPoint& firstSegmentFistPoint, const cartesianPoint& firstSegmentSecondPoint, const cartesianPoint& secondSegmentFistPoint, const cartesianPoint& secondSegmentSecondPoint, cartesianPoint& intersection)
    {
        typedef boost::geometry::model::segment<cartesianPoint> Segment;
        Segment firstSegment(firstSegmentFistPoint, firstSegmentSecondPoint);
        Segment secondSegment(secondSegmentFistPoint, secondSegmentSecondPoint);

        std::vector<cartesianPoint> intersections;
        boost::geometry::intersection(firstSegment, secondSegment, intersections);

        if (!intersections.empty())
        {
            intersection = intersections[0];
            return true;
        }
        return false;
    }

    // cartesian points
    template <>
    struct Operations<cartesianPoint>
    {
        static void normalVector(const cartesianPoint& firstPoint, const cartesianPoint& secondPoint, const cartesianPoint& orientation, cartesianPoint& result)
        {
            double dx = getDx(firstPoint, secondPoint);
            double dy = getDy(firstPoint, secondPoint);
            const double squaredDistance = dx * dx + dy * dy;
            if (squaredDistance != 0.0)
            {
                const double distance = sqrt(squaredDistance);
                result.x = dx / distance;
                result.y = dy / distance;
            }
        }

        static double getDx(const cartesianPoint& firstPoint, const cartesianPoint& secondPoint)
        {
            return firstPoint.x - secondPoint.x;
        }

        static double getDy(const cartesianPoint& firstPoint, const cartesianPoint& secondPoint)
        {
            return firstPoint.y - secondPoint.y;
        }

        static void add(cartesianPoint& point, const cartesianPoint& normal, const double increment)
        {
            point.x = point.x + normal.x * increment;
            point.y = point.y + normal.y * increment;
        }

        static void referencePoint(const std::vector<cartesianPoint>& polygon, double& minX, double& minY)
        {
            minX = std::numeric_limits<double>::max();
            minY = std::numeric_limits<double>::max();
            for (const auto& point : polygon)
            {
                if (point.x < minX)
                {
                    minX = point.x;
                }
                if (abs(point.y) < abs(minY))
                {
                    minY = point.y;
                }
            }
        }

        static double distance(const cartesianPoint& firstPoint, const cartesianPoint& secondPoint)
        {
            double dx = getDx(firstPoint, secondPoint);
            double dy = getDy(firstPoint, secondPoint);
            const double squaredDistance = dx * dx + dy * dy;
            double distance = 0.0;
            if (squaredDistance != 0.0)
            {
                distance = sqrt(squaredDistance);
            }
            return distance;
        }

        static double innerProductTwoSegments(const cartesianPoint& firstPointFirstSegment, const cartesianPoint& secondPointFirstSegment, const cartesianPoint& firstPointSecondSegment, const cartesianPoint& secondPointSecondSegment)
        {
            double dx1 = getDx(firstPointFirstSegment, secondPointFirstSegment);
            double dx2 = getDx(firstPointSecondSegment, secondPointSecondSegment);

            double dy1 = getDy(firstPointFirstSegment, secondPointFirstSegment);
            double dy2 = getDy(firstPointSecondSegment, secondPointSecondSegment);
         
            return dotProduct(dx1, dx2, dy1, dy2);
        }
    };
    
    // spherical point
    template <>
    struct Operations<sphericalPoint>
    {
        static void normalVector(const sphericalPoint& firstPoint, const sphericalPoint& secondPoint, const sphericalPoint& orientation, sphericalPoint& result)
        {
            cartesian3DPoint firstPointCartesianCoordinates;
            cartesian3DPoint secondPointCartesianCoordinates;
            sphericalToCartesian(firstPoint, firstPointCartesianCoordinates);
            sphericalToCartesian(secondPoint, secondPointCartesianCoordinates);

            double lambda = orientation.x * degrad_hp;
            double phi = orientation.y * degrad_hp;
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

        static double getDx(const sphericalPoint& firstPoint, const sphericalPoint& secondPoint)
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

        static double getDy(const sphericalPoint& firstPoint, const sphericalPoint& secondPoint)
        {
            double firstPointY = firstPoint.y * degrad_hp;
            double secondPointY = secondPoint.y * degrad_hp;
            double dy = earth_radius * (secondPointY - firstPointY);
            return dy;
        }

        static void add(sphericalPoint& point, const sphericalPoint& normal, const double increment)
        {
            double convertedIncrement = raddeg_hp * increment / earth_radius;
            double xf = 1.0 / cos(degrad_hp * point.y);
            point.x = point.x + normal.x * convertedIncrement * xf;
            point.y = point.y + normal.y * convertedIncrement;
        }

        static void referencePoint(std::vector<sphericalPoint>& polygon, double& minX, double& minY)
        {
            minX = std::numeric_limits<double>::max();
            minY = std::numeric_limits<double>::max();
            double xmax = std::numeric_limits<double>::min();
            for (const auto& point : polygon)
            {
                if (point.x < minX)
                {
                    minX = point.x;
                }
                if (abs(point.y) < abs(minY))
                {
                    minY = point.y;
                }
                if (point.x> xmax)
                {
                    xmax = point.x;
                }
            }

            if (xmax - minX > 180.0)
            {
                double deltaX = xmax - 180.0;
                for (auto& point : polygon)
                {
                    if(point.x< deltaX)
                    {
                        point.x = point.x + 360.0;
                    } 
                }
                minX = minX + 360.0;
            }
            //TODO: check result
            minX = std::min_element(polygon.begin(), polygon.end(), [](const sphericalPoint& p1, const sphericalPoint& p2) { return p1.x < p2.x; })->x;
        }

        static double distance(const sphericalPoint& firstPoint, const sphericalPoint& secondPoint)
        {
            return -1.0;
        }

        static double innerProductTwoSegments(const sphericalPoint& firstPointFirstSegment, const sphericalPoint& secondPointFirstSegment, const sphericalPoint& firstPointSecondSegment, const sphericalPoint& secondPointSecondSegment)
        {
            cartesian3DPoint firstPointFirstSegment3D;
            cartesian3DPoint secondPointFirstSegment3D;
            cartesian3DPoint firstPointSecondSegment3D;
            cartesian3DPoint secondPointSecondSegment3D;

            sphericalToCartesian(firstPointFirstSegment, firstPointFirstSegment3D);
            sphericalToCartesian(secondPointFirstSegment, secondPointFirstSegment3D);
            sphericalToCartesian(firstPointSecondSegment, firstPointSecondSegment3D);
            sphericalToCartesian(secondPointSecondSegment, secondPointSecondSegment3D);

            double dx1 = secondPointFirstSegment3D.x - firstPointFirstSegment3D.x;
            double dy1 = secondPointFirstSegment3D.y - firstPointFirstSegment3D.y;
            double dz1 = secondPointFirstSegment3D.z - firstPointFirstSegment3D.z;

            double dx2 = secondPointSecondSegment3D.x - firstPointSecondSegment3D.x;
            double dy2 = secondPointSecondSegment3D.y - firstPointSecondSegment3D.y;
            double dz2 = secondPointSecondSegment3D.z - firstPointSecondSegment3D.z;

            return dotProduct(dx1, dx2, dy1, dy2, dz1, dz2);
        }
    };
    

}
#endif