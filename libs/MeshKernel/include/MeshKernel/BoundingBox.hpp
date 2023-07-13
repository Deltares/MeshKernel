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

#include "Entities.hpp"

#include <algorithm>

namespace meshkernel
{
    class BoundingBox
    {
    public:
        BoundingBox() = default;

        BoundingBox(const Point& m_lower_left, const Point& m_upper_right)
            : m_lowerLeft(m_lower_left),
              m_upperRight(m_upper_right)
        {
        }

        /// @brief Given a vector of coordinates, get the lowest upper and right points
        /// @tparam T Requires IsCoordinate<T>
        /// @param[in] points The point values
        /// @returns A tuple with bottom left and upper right corners of the bounding box
        template <typename T>
        BoundingBox(const std::vector<T>& points)
        {
            double minx = std::numeric_limits<double>::max();
            double maxx = std::numeric_limits<double>::lowest();
            double miny = std::numeric_limits<double>::max();
            double maxy = std::numeric_limits<double>::lowest();

            for (const auto& point : points)
            {
                if (point.IsValid())
                {
                    minx = std::min(minx, point.x);
                    maxx = std::max(maxx, point.x);
                    miny = std::min(miny, point.y);
                    maxy = std::max(maxy, point.y);
                }
            }
            m_lowerLeft = Point(minx, miny);
            m_upperRight = Point(maxx, maxy);
        };

        /// @brief Checks if one bounding box is contained in onother
        /// @tparam    T          Requires IsCoordinate<T>
        /// @param[in] other      The other bounding box
        /// @param[in] lowerLeft  The lower left corner of the bounding box
        /// @param[in] upperRight The upper right corner of the bounding box
        /// @returns If the point is in the bounding box
        bool
        IsContained(const BoundingBox& other) const
        {
            return IsContained(other.m_lowerLeft) && IsContained(other.m_upperRight);
        }

        bool IsEqual(const BoundingBox& other) const
        {
            return other.m_lowerLeft == m_lowerLeft && other.m_upperRight == m_upperRight;
        }

        /// @brief Checks if value is inside a bounding box
        /// @tparam    T          Requires IsCoordinate<T>
        /// @param[in] point      The point to inquire
        template <typename T>
        bool IsContained(T point) const
        {

            return point.x >= m_lowerLeft.x && point.x <= m_upperRight.x &&
                   point.y >= m_lowerLeft.y && point.y <= m_upperRight.y;
        }

        const auto lowerLeft() const { return m_lowerLeft; }
        const auto upperRight() const { return m_upperRight; }

    private:
        Point m_lowerLeft;
        Point m_upperRight;
    };
} // namespace meshkernel
