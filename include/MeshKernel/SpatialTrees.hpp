//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2020.
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
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Constants.hpp>

// include boost
#define BOOST_ALLOW_DEPRECATED_HEADERS
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>
#undef BOOST_ALLOW_DEPRECATED_HEADERS

#include <vector>
#include <utility>
#include <stdexcept>

// r-tree
// https://gist.github.com/logc/10272165
// boost queries
// https://www.boost.org/doc/libs/1_66_0/libs/geometry/doc/html/geometry/spatial_indexes/queries.html

namespace meshkernel
{
    namespace SpatialTrees
    {
        namespace bg = boost::geometry;
        namespace bgi = boost::geometry::index;
        constexpr int QueryVectorCapacity = 100;

        class RTree
        {

            typedef bg::model::point<double, 2, bg::cs::cartesian> Point2D;
            typedef bg::model::box<Point2D> Box2D;
            typedef std::pair<Point2D, int> value2D;
            typedef bgi::rtree<value2D, bgi::linear<16>> RTree2D;

            typedef bg::model::point<double, 3, bg::cs::cartesian> Point3D;
            typedef bg::model::box<Point2D> Box3D;
            typedef std::pair<Point3D, int> value3D;
            typedef bgi::rtree<value3D, bgi::linear<16>> RTree3D;

        public:
            template <typename T>
            void BuildTree(std::vector<T>& nodes) //requires IsCoordinate<T>
            {
                m_points.reserve(m_points.size());
                m_points.clear();
                m_rtree2D.clear();

                for (int n = 0; n < nodes.size(); ++n)
                {
                    if (nodes[n].x != doubleMissingValue && nodes[n].y != doubleMissingValue)
                    {
                        m_points.emplace_back(Point2D{nodes[n].x, nodes[n].y}, n);
                    }
                }
                m_rtree2D = RTree2D(m_points.begin(), m_points.end());
            }

            void NearestNeighboursOnSquaredDistance(Point node, double searchRadiusSquared)
            {
                double searchRadius = std::sqrt(searchRadiusSquared);

                Box2D box(Point2D(node.x - searchRadius, node.y - searchRadius), Point2D(node.x + searchRadius, node.y + searchRadius));
                Point2D nodeSought = Point2D(node.x, node.y);

                m_queryCache.reserve(QueryVectorCapacity);
                m_queryCache.clear();
                m_rtree2D.query(
                    bgi::within(box) &&
                        bgi::satisfies([&nodeSought, &searchRadiusSquared](value2D const& v) { return bg::comparable_distance(v.first, nodeSought) <= searchRadiusSquared; }),
                    std::back_inserter(m_queryCache));

                m_queryIndices.reserve(m_queryCache.size());
                m_queryIndices.clear();
                for (const auto& v : m_queryCache)
                {
                    m_queryIndices.emplace_back(v.second);
                }
            }

            void NearestNeighbour(Point node)
            {

                m_queryCache.reserve(QueryVectorCapacity);
                m_queryCache.clear();
                Point2D nodeSought = Point2D(node.x, node.y);
                m_rtree2D.query(bgi::nearest(nodeSought, 1), std::back_inserter(m_queryCache));

                if (!m_queryCache.empty())
                {
                    m_queryIndices.clear();
                    m_queryIndices.emplace_back(m_queryCache[0].second);
                }
            }

            void RemoveNode(int position)
            {
                const auto numberRemoved = m_rtree2D.remove(m_points[position]);
                if (numberRemoved != 1)
                {
                    throw std::invalid_argument("SpatialTrees::RemoveNode: Could not remove node at given position.");
                }
                m_points[position] = {Point2D{doubleMissingValue, doubleMissingValue}, std::numeric_limits<size_t>::max()};
            }

            void InsertNode(const Point& node)
            {
                m_points.emplace_back(Point2D{node.x, node.y}, m_points.size());
                m_rtree2D.insert(m_points.end() - 1, m_points.end());
            }

            [[nodiscard]] auto Size() const
            {
                return m_rtree2D.size();
            }

            [[nodiscard]] auto Empty() const
            {
                return m_rtree2D.empty();
            }

            [[nodiscard]] auto GetQueryResultSize() const
            {
                return m_queryCache.size();
            }

            [[nodiscard]] auto GetQuerySampleIndex(int index) const
            {
                return m_queryIndices[index];
            }

        private:
            RTree2D m_rtree2D;
            std::vector<std::pair<Point2D, size_t>> m_points;
            std::vector<value2D> m_queryCache;
            std::vector<int> m_queryIndices;
        };

    } // namespace SpatialTrees
} // namespace meshkernel
