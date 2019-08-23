#include <cmath> 
#include <utility>

struct SphericalPoint
{


    // constructor
    SphericalPoint(double x, double y)
        : x(x), y(y)
    {
    };

    // default point
    SphericalPoint()
        : x(999.0), y(999.0)
    {
    };

    //Normalized vector in direction orientation
    static SphericalPoint normalVector(const SphericalPoint& firstPoint, const SphericalPoint& secondPoint, const SphericalPoint& orientation)
    {
        double dx = getDx(firstPoint, secondPoint);
        double dy = getDy(firstPoint, secondPoint);
        double squaredDistance = dx * dx + dy * dy;
        if (squaredDistance != 0.0)
        {
            const double distance = sqrt(squaredDistance);
            SphericalPoint point{ dx / distance , dy / distance };
            return std::move(point);
        }
    }

    static double getDx(const SphericalPoint& firstPoint, const SphericalPoint& secondPoint)
    {
        return firstPoint.x - secondPoint.x;
    }

    static double getDy(const SphericalPoint& firstPoint, const SphericalPoint& secondPoint)
    {
        return firstPoint.y - secondPoint.y;
    }

    static void add(SphericalPoint& point, const SphericalPoint& normal, const double increment)
    {
        point.x = point.x + normal.x * increment;
        point.y = point.y + normal.y * increment;
    }

    // private members
    double x;
    double y;

};
