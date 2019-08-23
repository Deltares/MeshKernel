#include <cmath> 
#include <utility>

struct CartesianPoint
{
    // constructor
    CartesianPoint(double x, double y)
        : x(x), y(y)
    {
    };

    // default point
    CartesianPoint()
        : x(999.0), y(999.0)
    {
    };

    //Normalized vector in direction orientation
    static CartesianPoint normalVector(const CartesianPoint& firstPoint, const CartesianPoint& secondPoint, const CartesianPoint& orientation)
    {
        double dx = getDx(firstPoint, secondPoint);
        double dy = getDy(firstPoint, secondPoint);
        double squaredDistance = dx * dx + dy * dy;
        if (squaredDistance != 0.0)
        {
            const double distance = sqrt(squaredDistance);
            CartesianPoint point{ dx / distance , dy / distance };
            return std::move(point);
        }
    }

    static double getDx(const CartesianPoint& firstPoint, const CartesianPoint& secondPoint)
    {
        return firstPoint.x - secondPoint.x;
    }

    static double getDy(const CartesianPoint& firstPoint, const CartesianPoint& secondPoint)
    {
        return firstPoint.y - secondPoint.y;
    }

    static void add(CartesianPoint& point, const CartesianPoint& normal, const double increment)
    {
        point.x = point.x + normal.x * increment;
        point.y = point.y + normal.y * increment;
    }

    // private members
    double x;
    double y;

};
