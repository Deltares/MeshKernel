#pragma once

#include <cmath>
#include "Entities.hpp"
#include "Constants.cpp"
#include "IOperations.hpp"

#ifdef USE_BOOST 
#include <boost/geometry/geometries/segment.hpp> 
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>
#include <boost/geometry.hpp>
BOOST_GEOMETRY_REGISTER_POINT_2D(GridGeom::Point, double, boost::geometry::cs::cartesian, x, y);
#endif 

namespace GridGeom
{
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
    static void sphericalToCartesian(const Point& sphericalPoint, cartesian3DPoint& cartesianPoint)
    {
        cartesianPoint.z = earth_radius * sin(sphericalPoint.y * degrad_hp);
        double rr = earth_radius * cos(sphericalPoint.y * degrad_hp);
        cartesianPoint.x = rr * cos(sphericalPoint.x * degrad_hp);
        cartesianPoint.y = rr * sin(sphericalPoint.x * degrad_hp);
    }

    //  transform 3d cartesian coordinates to 2d spherical
    static void cartesianToSpherical(const cartesian3DPoint& cartesianPoint, const double referenceLongitude, Point& sphericalPoint)
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
    static double isLeft(const Point& leftPoint, const Point& rightPoint, const Point& point)
    {
        double left = (rightPoint.x - leftPoint.x) * (point.y - leftPoint.y) - (point.x - leftPoint.x) * (rightPoint.y - leftPoint.y);
        return left;
    }

    // check if a point is in polygon using the winding number method
    // polygon: vector of points in counter clockwise order
    static bool pointInPolygon(const Point& point, const std::vector<Point>& polygon, const int numberOfPolygonPoints)
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

    template <typename T>
    T findIndex(const std::vector<T>& vec, const T& el)
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

    //faceAreaAndCenterOfMass: for cartesian, spherical point and spherical3dPoint
    static bool faceAreaAndCenterOfMass(std::vector<Point>& polygon, const int numberOfPolygonPoints, double& area, Point& centerOfMass, IOperations* operations)
    {
        if (numberOfPolygonPoints < 1)
        {
            return false;
        }
        
        double minX = std::numeric_limits<double>::max();
        double minY = std::numeric_limits<double>::max();
        operations->referencePoint(polygon, numberOfPolygonPoints, minX, minY);

        Point reference{ minX, minY };
        area = 0.0;
        double xCenterOfMass = 0.0;
        double yCenterOfMass = 0.0;
        for (int p = 0; p < numberOfPolygonPoints; p++)
        {
            double dx0 = operations->getDx(reference, polygon[p]);
            double dy0 = operations->getDy(reference, polygon[p]);
            double dx1 = operations->getDx(reference, polygon[p + 1]);
            double dy1 = operations->getDy(reference, polygon[p + 1]);

            double xc = 0.5 * (dx0 + dx1);
            double yc = 0.5 * (dy0 + dy1);

            dx0 = operations->getDx(polygon[p], polygon[p + 1]);
            dy0 = operations->getDy(polygon[p], polygon[p + 1]);
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