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

#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Polygons.hpp"

namespace meshkernel
{

    /// @brief A class containing the land boundary polylines
    class LandBoundary
    {
    public:
        /// @brief Construct with vector of points defining the land boundary
        explicit LandBoundary(const std::vector<Point>& landBoundary);

        /// @brief Find the nearest point on the land boundary (toland)
        void FindNearestPoint(const Point& samplePoint,
                              const Projection& projection,
                              Point& nearestPoint,
                              double& minimumDistance,
                              UInt& segmentStartIndex,
                              double& scaledDistanceToStart) const;

        /// @brief Find the nearest point on the land boundary (toland)
        Point FindNearestPoint(const Point& samplePoint,
                               const Projection& projection) const;

        /// @brief Gets the number of land boundary nodes.
        size_t GetNumNodes() const;

        ///@ Determine if the land boundary object is empty
        bool IsEmpty() const;

        /// @brief Get the node at position i.
        const Point& Node(const size_t i) const;

        /// @brief Get vector containing all land boundary nodes.
        const std::vector<Point>& GetNodes() const;

        /// @brief Add a new land boundary polyline segment
        void AddSegment(const Point& leftNode, const Point& rightNode);

        /// @brief Find the closest of two points to a given point.
        Point ClosestPoint(const Point& point, const size_t point1Index, const size_t point2Index, const Projection projection) const;

        /// @brief Find all start-end positions of the individual poly-lines that make up the land boundary
        std::vector<std::pair<UInt, UInt>> FindPolylineIndices() const;

        /// @brief Get vector of Boolean values indicating a valid node
        std::vector<bool> GetNodeMask(const Polygons& polygons) const;

    private:
        /// @brief The nodes making up the land boundary (XLAN, YLAN)
        std::vector<Point> m_nodes;
    };

} // namespace meshkernel

inline size_t meshkernel::LandBoundary::GetNumNodes() const
{
    return m_nodes.size();
}

inline bool meshkernel::LandBoundary::IsEmpty() const
{
    return m_nodes.empty();
}

inline const meshkernel::Point& meshkernel::LandBoundary::Node(const size_t i) const
{
    return m_nodes[i];
}

inline const std::vector<meshkernel::Point>& meshkernel::LandBoundary::GetNodes() const
{
    return m_nodes;
}
