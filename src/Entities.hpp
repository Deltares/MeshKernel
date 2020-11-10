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

#pragma once

#include <utility>
#include <vector>
#include "Constants.cpp"
#include <cmath>
#include <limits>

namespace meshkernel
{
    template <typename T>
    static bool IsDifferenceLessThanEpsilon(T firstValue, T secondValue)
    {
        return std::abs(firstValue - secondValue) < std::numeric_limits<T>::epsilon();
    }

    enum class Projections
    {
        cartesian,        // jsferic  = 0
        spherical,        // jsferic  = 1
        sphericalAccurate // jasfer3D = 1
    };

    struct Point
    {
        double x;
        double y;

        Point operator+(Point const& rhs) const
        {
            Point point{
                x + rhs.x,
                y + rhs.y};
            return point;
        }

        Point operator+(double const& rhs) const
        {
            Point point{
                x + rhs,
                y + rhs};
            return point;
        }

        Point operator-(Point const& rhs) const
        {
            Point point{
                x - rhs.x,
                y - rhs.y};
            return point;
        }

        Point operator-(double const& rhs) const
        {
            Point point{
                x - rhs,
                y - rhs};
            return point;
        }

        Point operator*(Point const& rhs) const
        {
            Point point{
                x * rhs.x,
                y * rhs.y};
            return point;
        }

        Point operator*(double const& rhs) const
        {
            Point point{
                x * rhs,
                y * rhs};
            return point;
        }

        Point operator/(Point const& rhs) const
        {
            Point point{
                x / rhs.x,
                y / rhs.y};
            return point;
        }

        Point operator/(double const& rhs) const
        {
            Point point{
                x / rhs,
                y / rhs};
            return point;
        }

        bool operator==(const Point& rhs) const
        {
            bool isEqual = IsDifferenceLessThanEpsilon(x, rhs.x) &&
                           IsDifferenceLessThanEpsilon(y, rhs.y);

            return isEqual;
        }

        bool operator!=(const Point& rhs) const
        {
            bool isEqual = IsDifferenceLessThanEpsilon(x, rhs.x) &&
                           IsDifferenceLessThanEpsilon(y, rhs.y);
            return !isEqual;
        }

        void TransformSphericalToCartesian(double referenceLatitude)
        {
            x = x * degrad_hp * earth_radius * std::cos(degrad_hp * referenceLatitude);
            y = y * degrad_hp * earth_radius;
        }

        [[nodiscard]] bool IsValid(const double missingValue = doubleMissingValue) const
        {
            bool isInvalid = IsDifferenceLessThanEpsilon(x, missingValue) ||
                             IsDifferenceLessThanEpsilon(y, missingValue);

            return !isInvalid;
        }
    };

    typedef std::pair<int, int> Edge;

    struct Cartesian3DPoint
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

    static std::vector<Edge> ConvertToEdgeNodesVector(int numEdges, const int* edge_nodes)
    {
        std::vector<Edge> edges(numEdges);

        int ei = 0;
        for (int e = 0; e < numEdges; e++)
        {
            edges[e].first = edge_nodes[ei];
            ei++;
            edges[e].second = edge_nodes[ei];
            ei++;
        }
        return edges;
    }

    static std::vector<Point> ConvertToNodesVector(int numNodes, const double* nodex, const double* nodey)
    {
        std::vector<Point> nodes(numNodes);
        for (int n = 0; n < numNodes; n++)
        {
            nodes[n].x = nodex[n];
            nodes[n].y = nodey[n];
        }
        return nodes;
    }

    static std::vector<Point> ConvertToFaceCentersVector(int numFaces, const double* facex, const double* facey)
    {
        std::vector<Point> faceCenters(numFaces);
        for (int n = 0; n < numFaces; n++)
        {
            faceCenters[n].x = facex[n];
            faceCenters[n].y = facey[n];
        }
        return faceCenters;
    };

} // namespace meshkernel
