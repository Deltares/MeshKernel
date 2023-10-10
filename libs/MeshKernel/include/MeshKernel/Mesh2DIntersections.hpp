//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2023.
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
#include <algorithm>
#include <vector>

#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/Polygons.hpp"

namespace meshkernel
{
    /// An intersection with a mesh edge
    struct EdgeMeshPolylineIntersection
    {
        int polylineSegmentIndex{constants::missing::intValue};                      ///< The intersected segment index (a polyline can formed by several segments)
        double polylineDistance{constants::missing::doubleValue};                    ///< The location of the intersection expressed as distance from the polyline start
        double adimensionalPolylineSegmentDistance{constants::missing::doubleValue}; ///< The location of the intersection expressed as an adimensional distance from the segment start
        UInt edgeIndex{constants::missing::uintValue};                               ///< The edge index
        UInt edgeFirstNode{constants::missing::uintValue};                           ///< The first node of the edge is on the left (the virtual node)
        UInt edgeSecondNode{constants::missing::uintValue};                          ///< The second node of the edge is on the right (the inner node)
        double edgeDistance{constants::missing::doubleValue};                        ///< The location of the intersection expressed as an adimensional distance from the edge start
    };

    /// An intersection with a mesh face
    struct FaceMeshPolylineIntersection
    {
        double polylineDistance{constants::missing::doubleValue}; ///< The location of the intersection expressed as an adimensional distance from the polyline start
        UInt faceIndex{constants::missing::uintValue};            ///< The face index
        std::vector<UInt> edgeIndexses;                           ///< The indexes of crossed edges
        std::vector<UInt> edgeNodes;                              ///< The indexes of the nodes defining the crossed edges
    };

    /// @brief Compute the intersections of polygon inner and outer perimeters
    ///
    /// @note Uses a breadth first algorithm to reduce runtime complexity
    class Mesh2DIntersections final
    {
    public:
        /// @brief Constructor
        Mesh2DIntersections(Mesh2D& mesh);

        /// @brief Compute intersection with a polygon (multiple polylines)
        /// @param[in] polyLine An input polygon
        void Compute(const Polygons& polygon);

        /// @brief Compute intersection with a single polyline
        /// @param[in] polyLine An input polyline
        void Compute(const std::vector<Point>& polyLine);

        /// @brief Gets the edge intersections
        /// @returns The edges intersections
        [[nodiscard]] const auto& EdgeIntersections() const { return m_edgesIntersections; }

        /// @brief  Gets the face intersections
        /// @returns The faces intersections
        [[nodiscard]] const auto& FaceIntersections() const { return m_faceIntersections; }

        static void updateEdgeIntersections(const UInt segmentIndex,
                                            const UInt edgeIndex,
                                            const UInt edgeFirstNode,
                                            const UInt edgeSecondNode,
                                            const std::vector<double>& cumulativeLength,
                                            const double crossProductValue,
                                            const double adimensionalEdgeDistance,
                                            const double adimensionalPolylineSegmentDistance,
                                            std::vector<EdgeMeshPolylineIntersection>& intersections)
        {
            intersections[edgeIndex].polylineSegmentIndex = static_cast<int>(segmentIndex);
            intersections[edgeIndex].polylineDistance = cumulativeLength[segmentIndex] +
                                                        adimensionalPolylineSegmentDistance * (cumulativeLength[segmentIndex + 1] - cumulativeLength[segmentIndex]);
            intersections[edgeIndex].adimensionalPolylineSegmentDistance = adimensionalPolylineSegmentDistance;
            intersections[edgeIndex].edgeFirstNode = crossProductValue < 0 ? edgeSecondNode : edgeFirstNode;
            intersections[edgeIndex].edgeSecondNode = crossProductValue < 0 ? edgeFirstNode : edgeSecondNode;
            intersections[edgeIndex].edgeDistance = adimensionalEdgeDistance;
            intersections[edgeIndex].edgeIndex = edgeIndex;
        }

        static void updateFaceIntersections(const UInt faceIndex,
                                            const UInt edgeIndex,
                                            std::vector<FaceMeshPolylineIntersection>& intersections)
        {
            intersections[faceIndex].faceIndex = faceIndex;
            intersections[faceIndex].edgeIndexses.emplace_back(edgeIndex);
        }

        template <typename T>
        static void sortAndEraseIntersections(std::vector<T>& intersections)
        {
            std::ranges::sort(intersections,
                              [](const T& first, const T& second)
                              { return first.polylineDistance < second.polylineDistance; });

            std::erase_if(intersections, [](const T& v)
                          { return v.polylineDistance < 0; });
        }

    private:
        std::tuple<bool, UInt, UInt> GetIntersectionSeed(const Mesh2D& mesh, const std::vector<Point>& polyLine, const std::vector<bool>& vistedEdges) const;

        Mesh2D& m_mesh;
        std::vector<EdgeMeshPolylineIntersection> m_edgesIntersectionsCache; ///< A cache for saving the edge intersections of one inner or outer
        std::vector<FaceMeshPolylineIntersection> m_facesIntersectionsCache; ///< A cache for saving the local face intersections of one inner or outer
        std::vector<EdgeMeshPolylineIntersection> m_edgesIntersections;      ///< A vector collecting all edge intersection results
        std::vector<FaceMeshPolylineIntersection> m_faceIntersections;       ///< A vector collecting all face intersection results
        static constexpr UInt maxSearchSegments = 10;                        ///< mex number of steps in polyline intersection algorithm
    };

} // namespace meshkernel
