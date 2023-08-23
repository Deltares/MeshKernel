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

#include <vector>

#include "MeshKernel/BoundingBox.hpp"
#include "MeshKernel/Entities.hpp"

namespace meshkernel
{

    /// @brief Forward declaration of the LandBoundary
    class LandBoundary;

    /// @brief A Polygon
    class Polygon
    {
    public:
        /// @brief Default constructor.
        Polygon() = default;

        /// @brief Default copy constructor.
        Polygon(const Polygon& copy) = default;

        /// @brief Default move constructor.
        Polygon(Polygon&& copy) = default;

        /// @brief Constructor
        Polygon(const std::vector<Point>& points,
                Projection projection);

        /// @brief Constructor
        /// @note Points are moved to the polygon
        Polygon(std::vector<Point>&& points,
                Projection projection);

        ///  @brief Reset the polygon
        void Reset(const std::vector<Point>& points,
                   Projection projection);

        Polygon& operator=(const Polygon& copy);

        Polygon& operator=(Polygon&& copy);

        // From here
        // Could be PointSequence?
        // Then polygon, spline and land boundary could use this.

        /// @brief Return the number of points in the polygon
        size_t Size() const;

        // TODO be consistent in naming
        /// @brief Return vector of points of the polygon
        const std::vector<Point>& Points() const;

        /// @brief Return single point at position
        Point& GetPoint(const size_t i);

        /// @brief Return single point at position
        const Point& GetPoint(const size_t i) const;

        /// @brief Determine if the polygon is closed or not.
        bool IsClosed() const;

        /// @brief Get the bounding box of the polygon
        const BoundingBox& GetBoundingBox() const;

        // To here

        /// @brief Determine if the polygon contains the point
        bool Contains(const Point& point) const;

        /// @brief Snap the section of the polygon defined by start- and end-index to the land boundary
        ///
        /// @note The bounding box may be changed
        void SnapToLandBoundary(const size_t startIndex, const size_t endIndex, const LandBoundary& boundary);

        /// @brief Refine the polygon
        /// @return The points for the refined polygon
        std::vector<Point> Refine(size_t startIndex, size_t endIndex, double refinementDistance) const;

        /// @brief Compute the area of the polygon, its centre of mass and the direction (true is anti-clockwise)
        std::tuple<double, Point, bool> FaceAreaAndCenterOfMass() const;

        // PerimeterClosed
        /// @brief Compute the perimiter length of the closed polygon
        double ClosedPerimeterLength() const;

        /// @brief Computes the edge lengths of the polygon
        /// @return edgeLengths The length of each polygon edge
        std::vector<double> EdgeLengths() const;

        /// @brief Compute the displaced poygon
        Polygon Displace(double displacement) const;

    private:
        /// @brief Determine if the polygon contains the point for Cartesian coordinate system
        ///
        /// Also for spherical coordinates
        bool ContainsCartesian(const Point& point) const;

        /// @brief Determine if the polygon contains the point for accurate spherical coordinate system
        bool ContainsSphericalAccurate(const Point& point) const;

        /// @brief The point sequence making up the corners of the polygon
        std::vector<Point> m_points;

        /// @brief The current projection
        Projection m_projection = Projection::cartesian;

        /// @brief The bounding box containing the polygon
        BoundingBox m_boundingBox;
    };

} // namespace meshkernel

inline size_t meshkernel::Polygon::Size() const
{
    return m_points.size();
}

inline const std::vector<meshkernel::Point>& meshkernel::Polygon::Points() const
{
    return m_points;
}

inline meshkernel::Point& meshkernel::Polygon::GetPoint(const size_t i)
{
    return m_points[i];
}

inline const meshkernel::Point& meshkernel::Polygon::GetPoint(const size_t i) const
{
    return m_points[i];
}

inline bool meshkernel::Polygon::IsClosed() const
{
    // m_points.back?
    return m_points[0] == m_points[m_points.size() - 1];
}

inline const meshkernel::BoundingBox& meshkernel::Polygon::GetBoundingBox() const
{
    return m_boundingBox;
}
