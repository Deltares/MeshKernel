#pragma once

#include <utility>
#include <vector>

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

    };

    struct cartesian3DPoint
    {
        double x;
        double y;
        double z;
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

    typedef std::pair<std::size_t, std::size_t> Edge;

}