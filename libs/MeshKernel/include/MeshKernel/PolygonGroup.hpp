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

    class PolygonGroup {
    public :


        // TODO which constructor is better?
        PolygonGroup (const std::vector<Point>& points,
                      Projection projection);

        PolygonGroup (const std::vector<Point>& points,
                      size_t start, size_t end,
                      Projection projection);

        const Polygon& Outer () const;

        size_t NumberOfInner () const;

        const Polygon& Inner (size_t i) const;

        /// @brief Determine if the point lies in the polygon
        ///
        /// If the point lies within the outer polygon but outside any inner polygons
        bool Contains (const Point& pnt) const;

        // TODO
        // 0 none, 1 in exterior, 2 in interior
        int ContainsRegion (const Point& pnt) const;

        void SnapToLandBoundary(size_t startIndex, size_t endIndex, const LandBoundary& landBoundary);

    private :

        using IndexRange = std::pair<UInt, UInt>;

        using IndexRangeArray = std::vector<IndexRange>;

        static Polygon ConstructPolygon (const std::vector<Point>& points,
                                         size_t start, size_t end,
                                         Projection projection);

        void ConstructOuterPolygon ( const std::vector<Point>& points,
                                     size_t start, size_t end,
                                     const IndexRangeArray& innerIndices,
                                     Projection projection);

        void ConstructInnerPolygons ( const std::vector<Point>& points,
                                      const IndexRangeArray& innerIndices,
                                      Projection projection);


        Polygon m_outer;

        std::vector<Polygon> m_inner;

    };

} // namespace meshkernel

inline const meshkernel::Polygon& meshkernel::PolygonGroup::Outer () const
{
    return m_outer;
}

inline size_t meshkernel::PolygonGroup::NumberOfInner () const
{
    return m_inner.size ();
}

inline const meshkernel::Polygon& meshkernel::PolygonGroup::Inner (const size_t i) const
{
    return m_inner [i];
}
