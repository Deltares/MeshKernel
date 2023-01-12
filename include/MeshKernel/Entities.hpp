//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2021.
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

#include <type_traits>

#include <cmath>

#include <MeshKernel/Constants.hpp>

namespace meshkernel
{
    // to re-enable when compiling with c++20 support
    // template <typename T>
    // concept IsCoordinate = requires(T t)
    //{
    //    t.x;
    //    t.y;
    //    t.IsValid
    //};

    /// @brief Generic function determining if two values are equal
    ///
    /// This is especially useful for floating point values.
    template <typename T>
    static bool IsEqual(T value, T referenceValue)
    {
        if constexpr (std::is_integral<T>::value)
        {
            return value == referenceValue;
        }
        else
        {
            if (value == referenceValue)
            {
                return true;
            }

            const auto absDiff = std::abs(value - referenceValue);
            const auto absValue = std::abs(value);
            const auto absReferenceValue = std::abs(referenceValue);
            const auto tol = 10 * std::numeric_limits<T>::epsilon();

            return absDiff < tol * std::min(absValue, absReferenceValue);
            
        }
    }

    /// @brief Enumerator describing the supported projections
    enum class Projection
    {
        cartesian = 0,        // jsferic  = 0
        spherical = 1,        // jsferic  = 1
        sphericalAccurate = 2 // jasfer3D = 1
    };

    /// @brief A struct describing a point in a two-dimensional space
    struct Point
    {
        double x; ///< X-coordinate
        double y; ///< Y-coordinate

        /// @brief Constructor initializing with missing values
        Point() : x(doubleMissingValue), y(doubleMissingValue) {}

        /// @brief Constructor initializing with given arguments
        /// @param[in] x
        /// @param[in] y
        Point(double x, double y) : x(x), y(y) {}

        /// @brief Overloads addition with another Point
        Point operator+(Point const& rhs) const
        {
            Point point{
                x + rhs.x,
                y + rhs.y};
            return point;
        }

        /// @brief Overloads addition with a double
        Point operator+(double const& rhs) const
        {
            Point point{
                x + rhs,
                y + rhs};
            return point;
        }

        /// @brief Overloads subtraction with another Point
        Point operator-(Point const& rhs) const
        {
            Point point{
                x - rhs.x,
                y - rhs.y};
            return point;
        }

        /// @brief Overloads subtraction with a double
        Point operator-(double const& rhs) const
        {
            Point point{
                x - rhs,
                y - rhs};
            return point;
        }

        /// @brief Overloads multiplication with another Point
        Point operator*(Point const& rhs) const
        {
            Point point{
                x * rhs.x,
                y * rhs.y};
            return point;
        }

        /// @brief Overloads multiplication with a double
        Point operator*(double const& rhs) const
        {
            Point point{
                x * rhs,
                y * rhs};
            return point;
        }

        /// @brief Overloads multiplication with a double
        Point operator*(int const& rhs) const
        {
            Point point{
                x * rhs,
                y * rhs};
            return point;
        }

        /// @brief Overloads division with another Point
        Point operator/(Point const& rhs) const
        {
            Point point{
                x / rhs.x,
                y / rhs.y};
            return point;
        }

        /// @brief Overloads division with a double
        Point operator/(double const& rhs) const
        {
            Point point{
                x / rhs,
                y / rhs};
            return point;
        }

        /// @brief Overloads equality with another Point
        bool operator==(const Point& rhs) const
        {
            const bool isEqual = IsEqual(x, rhs.x) &&
                                 IsEqual(y, rhs.y);

            return isEqual;
        }

        /// @brief Overloads inequality with another Point
        bool operator!=(const Point& rhs) const
        {
            const bool isEqual = IsEqual(x, rhs.x) &&
                                 IsEqual(y, rhs.y);
            return !isEqual;
        }

        /// @brief Transforms spherical coordinates to cartesian
        void TransformSphericalToCartesian(double referenceLatitude)
        {
            x = x * degrad_hp * earth_radius * std::cos(degrad_hp * referenceLatitude);
            y = y * degrad_hp * earth_radius;
        }

        /// @brief Determines if one of the point coordinates equals to \p missingValue
        [[nodiscard]] bool IsValid(const double missingValue = doubleMissingValue) const
        {
            const bool isInvalid = IsEqual(x, missingValue) ||
                                   IsEqual(y, missingValue) ||
                                   IsEqual(x, innerOuterSeparator) ||
                                   IsEqual(y, innerOuterSeparator);

            return !isInvalid;
        }
    };

    /// @brief Describes an edge with two indices
    typedef std::pair<size_t, size_t> Edge;

    /// @brief Get the index of the node on the other side of the edge
    /// @param[in] edge The given edge
    /// @param[in] node The node where we want the other one
    /// @returns Node index of other node of the edge
    size_t static OtherNodeOfEdge(const Edge& edge, size_t node)
    {
        return node == edge.first ? edge.second : edge.first;
    }

    /// @brief A struct describing the three coordinates in a cartesian projection.
    struct Cartesian3DPoint
    {
        double x; ///< X-coordinate
        double y; ///< Y-coordinate
        double z; ///< Z-coordinate
    };

    /// @brief A struct describing a sample with two coordinates and a value
    struct Sample
    {
        double x;     ///< X-coordinate
        double y;     ///< Y-coordinate
        double value; ///< Value

        /// @brief Convert double arrays to std::vector<Sample>
        /// @param[in] numSamples Number of samples
        /// @param[in] samplesXCoordinate X-coordinates of the samples
        /// @param[in] samplesYCoordinate Y-coordinates of the samples
        /// @param[in] samplesValue Values of the samples
        static auto ConvertToSamples(int numSamples, const double** samplesXCoordinate, const double** samplesYCoordinate, const double** samplesValue)
        {
            // Build the samples
            std::vector<Sample> samples(numSamples);
            for (size_t i = 0; i < samples.size(); ++i)
            {
                samples[i].x = (*samplesXCoordinate)[i];
                samples[i].y = (*samplesYCoordinate)[i];
                samples[i].value = (*samplesValue)[i];
            }
            return samples;
        }

        /// @brief Determines if the sample instance has valid coordinates
        [[nodiscard]] bool IsValid(const double missingValue = doubleMissingValue) const
        {
            bool isInvalid = IsEqual(x, missingValue) ||
                             IsEqual(y, missingValue);

            return !isInvalid;
        }
    };

    /// @brief Converts array of edge nodes to corresponding vector
    static std::vector<Edge> ConvertToEdgeNodesVector(int numEdges, const int* edge_nodes)
    {
        std::vector<Edge> edges(numEdges);

        int ei = 0;
        for (auto e = 0; e < numEdges; e++)
        {
            edges[e].first = edge_nodes[ei];
            ei++;
            edges[e].second = edge_nodes[ei];
            ei++;
        }
        return edges;
    }

    /// @brief Converts array of nodes to corresponding vector
    static std::vector<Point> ConvertToNodesVector(int numNodes, const double* node_x, const double* node_y)
    {
        std::vector<Point> nodes(numNodes);
        for (auto n = 0; n < numNodes; n++)
        {
            nodes[n].x = node_x[n];
            nodes[n].y = node_y[n];
        }
        return nodes;
    }

    /// @brief Converts array of face centers to corresponding vector
    static std::vector<Point> ConvertToFaceCentersVector(int numFaces, const double* facex, const double* facey)
    {
        std::vector<Point> faceCenters(numFaces);
        for (auto n = 0; n < numFaces; n++)
        {
            faceCenters[n].x = facex[n];
            faceCenters[n].y = facey[n];
        }
        return faceCenters;
    }

    /// @brief Converts array of face centers to corresponding vector
    static std::vector<std::vector<size_t>> ConvertToFaceNodesVector(int num_faces, const int* face_nodes, const int* nodes_per_face)
    {
        std::vector<std::vector<size_t>> result;
        result.reserve(num_faces);

        std::vector<size_t> nodes;
        size_t index = 0;
        for (auto f = 0; f < num_faces; f++)
        {
            nodes.clear();
            for (auto n = 0; n < nodes_per_face[f]; n++)
            {
                nodes.emplace_back(face_nodes[index]);
                index++;
            }
            result.emplace_back(nodes);
        }
        return result;
    }

} // namespace meshkernel
