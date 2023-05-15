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

#include <cmath>
#include <type_traits>
#include <vector>

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

    /// @brief Generic function for determining if two floating point values are equal
    /// @param[value] The value to compare
    /// @param[ref_value] The reference value to compare to
    /// @param[eps_mutilpier] Multiplier of machine precision
    /// @return Boolean indicating whether the value and reference value are equal within machine precision multiplied by the multiplier
    template <std::floating_point T>
    static bool IsEqual(const T value, T ref_value, T eps_mutilpier = 10.0)
    {
        if (value == ref_value)
        {
            return true;
        }
        const T abs_diff = std::abs(value - ref_value);
        const T abs_value = std::abs(value);
        const T abs_ref_value = std::abs(ref_value);
        static const T tol = eps_mutilpier * std::numeric_limits<T>::epsilon();
        return abs_diff < tol * std::min(abs_value, abs_ref_value);
    }

    /// @brief Enumerator describing the supported projections
    enum class Projection
    {
        cartesian = 0,        // jsferic  = 0
        spherical = 1,        // jsferic  = 1
        sphericalAccurate = 2 // jasfer3D = 1
    };

    /// @brief A struct describing a point in a two-dimensional space
    class Point
    {
    public:
        double x; ///< X-coordinate
        double y; ///< Y-coordinate

        /// @brief Constructor initializing with missing values
        Point()
            : x(constants::missing::doubleValue),
              y(constants::missing::doubleValue)
        {
        }

        /// @brief Constructor initializing with given coordinates
        /// @param[in] x coordinate
        /// @param[in] y coordinate
        Point(double x, double y)
            : x(x),
              y(y)
        {
        }

        /// @brief Overloads addition with another Point
        [[nodiscard]] Point operator+(Point const& rhs) const { return Point(x + rhs.x, y + rhs.y); }

        /// @brief Overloads addition with a double
        [[nodiscard]] Point operator+(double const& rhs) const { return Point(x + rhs, y + rhs); }

        /// @brief Overloads subtraction with another Point
        [[nodiscard]] Point operator-(Point const& rhs) const { return Point(x - rhs.x, y - rhs.y); }

        /// @brief Overloads subtraction with a double
        [[nodiscard]] Point operator-(double const& rhs) const { return Point(x - rhs, y - rhs); }

        /// @brief Overloads multiplication with another Point
        [[nodiscard]] Point operator*(Point const& rhs) const { return Point(x * rhs.x, y * rhs.y); }

        /// @brief Overloads multiplication with a double
        [[nodiscard]] Point operator*(double const& rhs) const { return Point(x * rhs, y * rhs); }

        /// @brief Overloads multiplication with a double
        [[nodiscard]] Point operator*(int const& rhs) const { return Point(x * rhs, y * rhs); }

        /// @brief Overloads division with another Point
        [[nodiscard]] Point operator/(Point const& rhs) const { return Point(x / rhs.x, y / rhs.y); }

        /// @brief Overloads division with a double
        [[nodiscard]] Point operator/(double const& rhs) const { return Point(x / rhs, y / rhs); }

        /// @brief Overloads equality with another Point
        bool operator==(const Point& rhs) const { return IsEqual(x, rhs.x) && IsEqual(y, rhs.y); }

        /// @brief Transforms spherical coordinates to cartesian
        void TransformSphericalToCartesian(double referenceLatitude)
        {
            x = x * m_trans_factor * std::cos(constants::conversion::degToRad * referenceLatitude);
            y = y * m_trans_factor;
        }

        /// @brief Determines if one of the point coordinates equals to \p missingValue
        [[nodiscard]] bool IsValid(const double missingValue = constants::missing::doubleValue) const
        {
            const bool isInvalid = IsEqual(x, missingValue) ||
                                   IsEqual(y, missingValue) ||
                                   IsEqual(x, constants::missing::innerOuterSeparator) ||
                                   IsEqual(y, constants::missing::innerOuterSeparator);

            return !isInvalid;
        }

    private:
        static double constexpr m_trans_factor =
            constants::conversion::degToRad *
            constants::geometric::earth_radius; ///< Factor used in the transformation from spherical to Cartesian coordinates
    };

    /// @brief Describes an edge with two indices
    typedef std::pair<size_t, size_t> Edge;

    /// @brief Get the index of the node on the other node of the edge
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
    class Sample : public Point
    {
    public:
        double value = constants::missing::doubleValue; ///< Value

        /// @brief Default constructor
        Sample() = default;

        /// @brief Constructor taking coordinates and values
        Sample(double x, double y, double value)
            : Point(x, y),
              value(value)
        {
        }

        /// @brief Convert double arrays to std::vector<Sample>
        /// @param[in] numSamples Number of samples
        /// @param[in] samplesXCoordinate X-coordinates of the samples
        /// @param[in] samplesYCoordinate Y-coordinates of the samples
        /// @param[in] samplesValue Values of the samples
        static auto ConvertToSamples(int numSamples, const double** samplesXCoordinate,
                                     const double** samplesYCoordinate,
                                     const double** samplesValue)
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
        [[nodiscard]] bool IsValid(const double missingValue = constants::missing::doubleValue) const
        {
            bool isInvalid = IsEqual(x, missingValue) ||
                             IsEqual(y, missingValue);

            return !isInvalid;
        }
    };

    /// @brief Converts array of edge nodes to corresponding vector
    static std::vector<Edge> ConvertToEdgeNodesVector(int numEdges, const int* const edge_nodes)
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
    static std::vector<Point> ConvertToNodesVector(int numNodes,
                                                   const double* const node_x,
                                                   const double* const node_y)
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
    static std::vector<Point> ConvertToFaceCentersVector(int numFaces,
                                                         const double* const facex,
                                                         const double* const facey)
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
    static std::vector<std::vector<size_t>> ConvertToFaceNodesVector(int num_faces,
                                                                     const int* const face_nodes,
                                                                     const int* const nodes_per_face)
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
