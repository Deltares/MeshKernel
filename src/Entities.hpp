#pragma once

#include <utility>
#include <vector>
#include "Constants.cpp"
#include <cmath>

namespace GridGeom
{
    enum OperationTypes 
    {
        cartesianOperations,
        sphericalOperations
    };

    struct Point
    {
        double x;
        double y;

        Point operator+(Point const& rhs) const 
        {
            Point point
            {
                x + rhs.x,
                y + rhs.y
            };
            return std::move(point);
        }

        Point operator+(double const& rhs) const
        {
            Point point
            {
                x + rhs,
                y + rhs
            };
            return std::move(point);
        }

        Point operator-(Point const& rhs) const
        {
            Point point
            {
                x - rhs.x,
                y - rhs.y
            };
            return std::move(point);
        }

        Point operator-(double const& rhs) const
        {
            Point point
            {
                x - rhs,
                y - rhs
            };
            return std::move(point);
        }

        Point operator*(Point const& rhs) const
        {
            Point point
            {
                x * rhs.x,
                y * rhs.y
            };
            return std::move(point);
        }

        Point operator*(double const& rhs) const
        {
            Point point
            {
                x * rhs,
                y * rhs
            };
            return std::move(point);
        }

        Point operator/(Point const& rhs) const
        {
            Point point
            {
                x / rhs.x,
                y / rhs.y
            };
            return std::move(point);
        }

        Point operator/(double const& rhs) const
        {
            Point point
            {
                x / rhs,
                y / rhs
            };
            return std::move(point);
        }

        bool operator==(const Point& rhs) const
        {
            return x==rhs.x && y == rhs.y;
        }

        bool operator!=(const Point& rhs) const
        {
            return x != rhs.x || y != rhs.y;
        }

        void TransformToSpherical() 
        {
            x = x * degrad_hp *earth_radius * std::cos(degrad_hp*y);
            y = y * degrad_hp *earth_radius;
        }

        bool IsValid(const double missingValue = doubleMissingValue) const
        {
            return x != missingValue && y != missingValue ? true : false;
        }
    };

    struct Vector
    {
        double x;
        double y;
    };

    struct cartesian3DPoint
    {
        double x;
        double y;
        double z;
    };

    struct Sample 
    {
        double x;
        double y;
        double value;
    };

    struct Nodes
    {
        std::vector<double> x;
        std::vector<double> y;
    };

    enum class Projections
    {
        cartesian,
        spherical
    };

    typedef std::pair<int, int> Edge;

}