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
#include <limits>

namespace meshkernel
{

    /// @brief A class defining a bounding box
    class BoundingBox
    {
    public:
        /// @brief Default constructor
        BoundingBox() : m_lowerLeft(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()),
                        m_upperRight(std::numeric_limits<double>::max(), std::numeric_limits<double>::max()) {}

        /// @brief Constructor taking the corner points of the bounding box
        /// @param[in] lowerLeft The lower left corner of the bounding box
        /// @param[in] upperRight The upper right corner of the bounding box
        BoundingBox(const Point& lowerLeft, const Point& upperRight)
            : m_lowerLeft(lowerLeft),
              m_upperRight(upperRight)
        {
        }

        /// @brief Constructor taking a vector of coordinates types
        /// @tparam T Requires IsCoordinate<T>
        /// @param[in] points The point values
        template <typename T>
        BoundingBox(const std::vector<T>& points) : BoundingBox(points, 0, points.size() - 1) {}

        // @brief Constructor taking a vector of coordinates types
        // @tparam T Requires IsCoordinate<T>
        // @param[in] points The point values
        // @param[in] start The start index for the array slice
        // @param[in] end The end index for the array slice
        template <typename T>
        BoundingBox(const std::vector<T>& points, size_t start, size_t end)
        {
            double minx = std::numeric_limits<double>::max();
            double maxx = std::numeric_limits<double>::lowest();
            double miny = std::numeric_limits<double>::max();
            double maxy = std::numeric_limits<double>::lowest();

            for (size_t i = start; i <= end; ++i)
            {
                const auto& point = points[i];

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
        }

        /// @brief Not equal operator
        /// @param[in] other The other bounding box to compare
        /// @return True if the other bounding box is not equal
        bool operator!=(const BoundingBox& other) const
        {
            return other.m_lowerLeft != m_lowerLeft || other.m_upperRight != m_upperRight;
        }

        /// @brief Checks if a point is inside a bounding box
        /// @tparam    T          Requires IsCoordinate<T>
        /// @param[in] point      The point to inquire
        /// @return True if the point is contained, false otherwise
        template <typename T>
        bool Contains(const T& point) const
        {

            return point.x >= m_lowerLeft.x && point.x <= m_upperRight.x &&
                   point.y >= m_lowerLeft.y && point.y <= m_upperRight.y;
        }

        /// @brief Returns the lower left corner of the bounding box
        /// @return The lower left corner of the bounding box
        [[nodiscard]] auto& lowerLeft() const { return m_lowerLeft; }

        /// @brief Returns the upper right corner
        /// @return The upper right corner of the bounding box
        [[nodiscard]] auto& upperRight() const { return m_upperRight; }

        /// @brief Returns the mass centre
        /// @return The upper right corner of the bounding box
        Point MassCentre() const { return (m_lowerLeft + m_upperRight) * 0.5; }

        /// @brief Returns the bounding box width
        /// @return The bounding box width
        double Width() const { return m_upperRight.x - m_lowerLeft.x; }

        /// @brief Returns the bounding box height
        /// @return The bounding box height
        double Height() const { return m_upperRight.y - m_lowerLeft.y; }

        /// @brief Extends the bounding box by a factor
        void ExtendBoundingBox(double factor)
        {
            const double width = Width();
            const double height = Height();
            m_lowerLeft.x = m_lowerLeft.x - width * factor;
            m_lowerLeft.y = m_lowerLeft.y - height * factor;
            m_upperRight.x = m_upperRight.x + width * factor;
            m_upperRight.y = m_upperRight.y + height * factor;
        }

        /// @brief Return the centre of the bounding box.
        Point centre() const;

        /// @brief Return the delta of the bounding box.
        Point delta() const;

    private:
        Point m_lowerLeft;  ///< The lower left corner of the bounding box
        Point m_upperRight; ///< The upper right corner of the bounding box
    };

    /// @brief Merge two bounding boxes into a single bounding box that will contain both of the original.
    BoundingBox merge(const BoundingBox& b1, const BoundingBox& b2);

} // namespace meshkernel

inline meshkernel::BoundingBox meshkernel::merge(const BoundingBox& b1, const BoundingBox& b2)
{
    Point lowerLeft{std::min(b1.lowerLeft().x, b2.lowerLeft().x), std::min(b1.lowerLeft().y, b2.lowerLeft().y)};
    Point upperRight{std::max(b1.upperRight().x, b2.upperRight().x), std::max(b1.upperRight().y, b2.upperRight().y)};

    return BoundingBox(lowerLeft, upperRight);
}

inline meshkernel::Point meshkernel::BoundingBox::centre() const
{
    return 0.5 * (m_lowerLeft + m_upperRight);
}

inline meshkernel::Point meshkernel::BoundingBox::delta() const
{
    return m_upperRight - m_lowerLeft;
}
