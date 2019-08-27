#ifndef OPERTIONS_CPP
#define OPERTIONS_CPP
#include <cmath>
#include "Entities.cpp"
#include "Constants.cpp"


namespace GridGeom
{

    // coordinate reference indipendent operations
    static double dotProduct(const double& dx1, const double& dx2, const double& dy1, const double& dy2)
    {
        return dx1 * dx2 + dy1 * dy2;
    }

    // transform 2d spherical to 3d cartesian
    static void sphericalToCartesian(const Point& sphericalPoint, Point3D& cartesianPoint)
    {
        cartesianPoint.z = earth_radius * sin(sphericalPoint.y * degrad_hp);
        double rr = earth_radius * cos(sphericalPoint.y * degrad_hp);
        cartesianPoint.x = rr * cos(sphericalPoint.x * degrad_hp);
        cartesianPoint.y = rr * sin(sphericalPoint.x * degrad_hp);
    }

    //  transform 3d cartesian coordinates to 2d spherical
    static void cartesianToSpherical(const Point3D& cartesianPoint, const double referenceLongitude, Point& sphericalPoint)
    {
        double angle = atan2(cartesianPoint.y, cartesianPoint.x) * raddeg_hp;
        sphericalPoint.y = atan2(cartesianPoint.z, sqrt(cartesianPoint.x * cartesianPoint.x + cartesianPoint.y * cartesianPoint.y)) * raddeg_hp;
        sphericalPoint.x = angle + std::lround((referenceLongitude - angle) / 360.0) * 360.0;
    }


    enum class CoordinateSystems
    {
        cartesian, 
        spheric
    };

    // generic template
    template <CoordinateSystems CoordinateSystems>
    struct Operations
    {
        static void normalVector(const Point& firstPoint, const Point& secondPoint, const Point& orientation, Point& result);
        static double getDx(const Point& firstPoint, const Point& secondPoint);
        static double getDy(const Point& firstPoint, const Point& secondPoint); 
        static void add(Point& point, const Point& normal, const double increment);
    };

    // cartesian system
    template <>
    struct Operations<CoordinateSystems::cartesian>
    {
        static void normalVector(const Point& firstPoint, const Point& secondPoint, const Point& orientation, Point& result)
        {
            double dx = getDx(firstPoint, secondPoint);
            double dy = getDy(firstPoint, secondPoint);
            double squaredDistance = dx * dx + dy * dy;
            if (squaredDistance != 0.0)
            {
                const double distance = sqrt(squaredDistance);
                result.x = dx / distance;
                result.y = dy / distance;
            }
        }
        static double getDx(const Point& firstPoint, const Point& secondPoint)
        {
            return firstPoint.x - secondPoint.x;
        }


        static double getDy(const Point& firstPoint, const Point& secondPoint)
        {
            return firstPoint.y - secondPoint.y;
        }


        static void add(Point& point, const Point& normal, const double increment)
        {
            point.x = point.x + normal.x * increment;
            point.y = point.y + normal.y * increment;
        }
    };

    // spheric system
    template <>
    struct Operations<CoordinateSystems::spheric>
    {

        static void normalVector(const Point& firstPoint, const Point& secondPoint, const Point& orientation, Point& result)
        {
            Point3D firstPointCartesianCoordinates;
            Point3D secondPointCartesianCoordinates;
            sphericalToCartesian(firstPoint, firstPointCartesianCoordinates);
            sphericalToCartesian(secondPoint, secondPointCartesianCoordinates);

            double lambda = orientation.x * degrad_hp;
            double phi = orientation.y * degrad_hp;
            double elambda[3] = { -sin(lambda), cos(lambda), 0.0 };
            double ephi[3] = { -sin(lambda), cos(lambda), 0.0 };

            //VectorED elambda={ -sin(lambda), cos(lambda), 0.0 };
            //VectorED ephi = { -sin(lambda), cos(lambda), 0.0 };
            //double dx = (elambda * ephi).sum();

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

        static double getDx(const Point& firstPoint, const Point& secondPoint)
        {
            double  firstPointYDiff = abs(abs(firstPoint.y) - 90.0);
            double  secondPointYDiff = abs(abs(secondPoint.y) - 90.0);
            if (firstPointYDiff <= dtol_pole && secondPointYDiff > dtol_pole || firstPointYDiff > dtol_pole && secondPointYDiff <= dtol_pole)
            {
                return 0.0;
            }
            double firstPointX = firstPoint.x;
            double secondPointX = secondPoint.x;
            if(firstPointX- secondPointX > 180.0)
            {
                firstPointX -= 360.0;
            }
            else if(firstPointX - secondPointX < - 180.0)
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

        static double getDy(const Point& firstPoint, const Point& secondPoint)
        {
            double firstPointY = firstPoint.y * degrad_hp;
            double secondPointY = secondPoint.y * degrad_hp;
            double dy = earth_radius * (secondPointY - firstPointY);
            return dy;
        }

        static void add(Point& point, const Point& normal, const double increment)
        {

            double convertedIncrement = raddeg_hp * increment / earth_radius;
            double xf = 1.0 / cos(degrad_hp * point.y);
            point.x = point.x + normal.x * convertedIncrement* xf;
            point.y = point.y + normal.y * convertedIncrement;
        }
    };


}

#endif