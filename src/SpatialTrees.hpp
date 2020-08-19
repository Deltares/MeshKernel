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
#include "Entities.hpp"
#include "Constants.cpp"

// include boost
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <vector>
#include <utility>

// r-tree
// https://gist.github.com/logc/10272165
// boost queries
// https://www.boost.org/doc/libs/1_66_0/libs/geometry/doc/html/geometry/spatial_indexes/queries.html

namespace GridGeom
{
    namespace SpatialTrees 
    {
        namespace bg = boost::geometry;
        namespace bgi = boost::geometry::index;

        class RTree
        {

            typedef bg::model::point<double, 2, bg::cs::cartesian> Point2D;
            typedef bg::model::box<Point2D> Box2D;
            typedef std::pair<Point2D, int> value2D;
            typedef bgi::rtree<value2D, bgi::quadratic<16>> RTree2D;

            typedef bg::model::point<double, 3, bg::cs::cartesian> Point3D;
            typedef bg::model::box<Point2D> Box3D;
            typedef std::pair<Point3D, int> value3D;
            typedef bgi::rtree<value3D, bgi::quadratic<16>> RTree3D;

        public:

            bool BuildTree(std::vector<Point>& nodes, const Projections projection)
            {
                m_points.resize(nodes.size());
                for (int n = 0; n < nodes.size(); ++n)
                {
                    m_points[n] = std::make_pair(Point2D{ nodes[n].x, nodes[n].y }, n);
                }
                m_rtree2D = RTree2D(m_points.begin(), m_points.end());
                
                m_queryCache.reserve(queryCapacity);
                m_queryIndexses.reserve(queryCapacity);

                return true;
            }

            bool NearestNeighboursOnSquaredDistance(Point node, double searchRadiusSquared)
            {
                m_queryCache.resize(0);
                double searchRadius = std::sqrt(searchRadiusSquared);
                
                Box2D box(Point2D(node.x - searchRadius, node.y - searchRadius), Point2D(node.x + searchRadius, node.y + searchRadius));
                Point2D nodeSought = Point2D(node.x, node.y);

                m_rtree2D.query(
                    bgi::within(box) &&
                    bgi::satisfies([&](value2D const& v) {return bg::comparable_distance(v.first, nodeSought) < searchRadiusSquared; }),
                    std::back_inserter(m_queryCache));

                m_queryIndexses.resize(m_queryCache.size());
                for (size_t i = 0; i < m_queryCache.size(); i++)
                {
                    m_queryIndexses[i] = m_queryCache[i].second;
                }

                return true;
            }

            bool NearestNeighbour(Point node)
            {
                m_queryCache.resize(0);

                Point2D nodeSought = Point2D(node.x, node.y);
                m_rtree2D.query(bgi::nearest(nodeSought,1),std::back_inserter(m_queryCache));

                if(!m_queryCache.empty())
                {
                    m_queryIndexses.resize(1);
                    m_queryIndexses[0] = m_queryCache[0].second;
                }

                return true;
            }

            bool RemoveNode(int position) 
            {
                int numberRemoved = m_rtree2D.remove(m_points[position]);
                if (numberRemoved != 1) 
                {
                    return false;
                }
                m_points[position] = std::make_pair(Point2D{ doubleMissingValue,doubleMissingValue }, -1);
                return true;
            }

            bool InsertNode(const Point& node)
            {
                const auto pointToInsert = std::make_pair(Point2D{ node.x, node.y }, m_points.size());
                m_points.push_back(pointToInsert);
                m_rtree2D.insert(m_points.end()-1, m_points.end());
                return true;
            }


            int Size() const
            {
                return m_rtree2D.size();
            }

            bool Empty() const
            {
                return m_rtree2D.empty();

            }

            bool Clear()
            {
                m_rtree2D.clear();
                return true;
            }

            int GetQueryResultSize() const
            {
                return m_queryCache.size();
            }

            int GetQuerySampleIndex(int index) const
            {
                return m_queryIndexses[index];
            }

        private:

            Projections m_projection;
            RTree2D m_rtree2D;
            std::vector<std::pair<Point2D, int>> m_points;
            std::vector<value2D> m_queryCache;
            std::vector<int> m_queryIndexses;
            int m_querySize = 0;

            // Rtree preallocated capacity of the query array 
            int queryCapacity = 100;
        };   

    }
}


