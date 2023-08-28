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
#include "MeshKernel/Polygon.hpp"

namespace meshkernel
{

    /// @brief A region enclosed by a polygonal permieter.
    ///
    /// Each region is described by an outer perimeter polygon.
    /// It may contain 0 or more holes, each described by an inner polygon.
    class PolygonalEnclosure
    {
    public:
        /// @brief The part of the enclosure a point is found
        enum class Region
        {
            None,     ///< The point is not contained with the enclosure
            Exterior, ///< The point is contained within the outer perimeter of the enclosure
            Interior  ///< The point is contained within one of the island, interior perimeters of the enclosure
        };

        /// @brief Constructor
        PolygonalEnclosure(const std::vector<Point>& points,
                           Projection projection);

        /// @brief The outer perimeter polygon
        const Polygon& Outer() const;

        /// @brief The number of inner hole regions
        size_t NumberOfInner() const;

        /// @brief Get an inner polygon
        const Polygon& Inner(size_t i) const;

        /// @brief Determine if the point lies in the polygon
        ///
        /// If the point lies within the outer polygon but outside any inner polygons
        bool Contains(const Point& pnt) const;

        /// @brief Determine in which part of the enclosure the point is.
        Region ContainsRegion(const Point& pnt) const;

        /// @brief Get the number of points making up the polygon, including interior if requested
        UInt NumberOfPoints(const bool includeInterior) const;

        /// @brief Snap all or part of the outer perimeter polygon to the land boundary
        void SnapToLandBoundary(size_t startIndex, size_t endIndex, const LandBoundary& landBoundary);

        /// @brief Refine the polygon.
        std::vector<Point> Refine(size_t startIndex, size_t endIndex, double refinementDistance);

    private:
        /// @typedef IndexRange
        /// @brief Contains the start and end of a section from the point array
        using IndexRange = std::pair<UInt, UInt>;

        /// @typedef IndexRangeArray
        /// @brief An array of IndexRange
        using IndexRangeArray = std::vector<IndexRange>;

        /// @brief Construct a polygon from a (sub) range of points.
        static Polygon ConstructPolygon(const std::vector<Point>& points,
                                        size_t start, size_t end,
                                        Projection projection);

        /// @brief Construct the outer perimiter polygon from the points.
        void ConstructOuterPolygon(const std::vector<Point>& points,
                                   size_t start, size_t end,
                                   const IndexRangeArray& innerIndices,
                                   Projection projection);

        /// @brief Construct all, if any, inner polygons from the points.
        void ConstructInnerPolygons(const std::vector<Point>& points,
                                    const IndexRangeArray& innerIndices,
                                    Projection projection);

        /// @brief The outer perimeter polygon of the enclosure
        Polygon m_outer;

        /// @brief The set of inner, island, polygons for the enclosure.
        std::vector<Polygon> m_inner;
    };

} // namespace meshkernel

inline const meshkernel::Polygon& meshkernel::PolygonalEnclosure::Outer() const
{
    return m_outer;
}

inline size_t meshkernel::PolygonalEnclosure::NumberOfInner() const
{
    return m_inner.size();
}

inline const meshkernel::Polygon& meshkernel::PolygonalEnclosure::Inner(const size_t i) const
{
    return m_inner[i];
}
