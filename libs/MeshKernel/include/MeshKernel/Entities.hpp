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

#include <vector>

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Point.hpp"

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

    /// @brief Describes an edge with two indices
    using Edge = std::pair<UInt, UInt>;

    /// @brief Get the index of the node on the other node of the edge
    /// @param[in] edge The given edge
    /// @param[in] node The node where we want the other one
    /// @returns Node index of other node of the edge
    UInt static OtherNodeOfEdge(const Edge& edge, UInt node)
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
            for (UInt i = 0; i < samples.size(); ++i)
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

    /// @brief Separate array of nodes to raw pointers to arrays of x-, y-coordinates and the size
    static std::tuple<int, double*, double*> ConvertFromNodesVector(const std::vector<Point>& nodes);

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
    static std::vector<std::vector<UInt>> ConvertToFaceNodesVector(int num_faces,
                                                                   const int* const face_nodes,
                                                                   const int* const nodes_per_face)
    {
        std::vector<std::vector<UInt>> result;
        result.reserve(num_faces);

        std::vector<UInt> nodes;
        UInt index = 0;
        for (auto f = 0; f < num_faces; f++)
        {
            nodes.clear();
            for (auto n = 0; n < nodes_per_face[f]; n++)
            {
                nodes.emplace_back(static_cast<UInt>(face_nodes[index]));
                index++;
            }
            result.emplace_back(nodes);
        }
        return result;
    }

} // namespace meshkernel
