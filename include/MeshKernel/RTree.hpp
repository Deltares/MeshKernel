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

#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Entities.hpp>

// include boost
#define BOOST_ALLOW_DEPRECATED_HEADERS
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>
#undef BOOST_ALLOW_DEPRECATED_HEADERS

// r-tree
// https://gist.github.com/logc/10272165

namespace meshkernel
{
    // using a namespace alias
    namespace bg = boost::geometry;
    namespace bgi = boost::geometry::index;

    /// @brief Class used for inquiring adjacent nodes.
    ///
    /// The mesh class stores two RTree class instances, used for inquiring the closest mesh nodes and edge to a point.
    /// RTree is a class wrapping the boost::geometry::index::rtree code,
    /// adding an interface for performing common queries
    /// such as inquiring the nearest neighbors inside a specified distance(`meshkernel::RTree::SearchPoints`)
    /// or a vector of the nearest neighbors (`meshkernel::RTree::SearchNearestPoint`).
    /// RTee has a `m_queryCache`, a vector used for collecting all query results
    /// and avoid frequent re-allocations when the number of results changes.
    class RTree
    {

        typedef bg::model::point<double, 2, bg::cs::cartesian> Point2D; ///< Typedef for Point2D
        typedef bg::model::box<Point2D> Box2D;                          ///< Typedef for box of Point2D
        typedef std::pair<Point2D, size_t> value2D;                     ///< Typedef of pair of Point2D and size_t
        typedef bgi::rtree<value2D, bgi::linear<16>> RTree2D;           ///< Typedef for a 2D RTree

        typedef bg::model::point<double, 3, bg::cs::cartesian> Point3D; ///< Typedef for Point3D
        typedef std::pair<Point3D, size_t> value3D;                     ///< Typedef of pair of Point3D and size_t
        typedef bgi::rtree<value3D, bgi::linear<16>> RTree3D;           ///< Typedef for a 3D RTree

    public:
        /// @brief Builds the tree
        /// @tparam T Requires IsCoordinate<T>
        template <typename T>
        void BuildTree(std::vector<T>& nodes)
        {
            m_points.reserve(m_points.size());
            m_points.clear();
            m_rtree2D.clear();

            for (size_t n = 0; n < nodes.size(); ++n)
            {
                if (nodes[n].x != doubleMissingValue && nodes[n].y != doubleMissingValue)
                {
                    m_points.emplace_back(Point2D{nodes[n].x, nodes[n].y}, n);
                }
            }
            m_rtree2D = RTree2D(m_points.begin(), m_points.end());
        }

        /// @brief Finds all nodes in the search radius and stores the results in the query cache, to be inquired later
        /// @param[in] node The node
        /// @param[in] searchRadiusSquared The squared search radius around the node
        void SearchPoints(Point node, double searchRadiusSquared);

        /// @brief Finds the nearest node in the search radius and stores the results in the query cache, to be inquired later
        /// @param[in] node The node
        /// @param[in] searchRadiusSquared The squared search radius around the node
        void SearchNearestPoint(Point node, double searchRadiusSquared);

        /// @brief Gets the nearest of all nodes
        /// @param[in] node The node
        void SearchNearestPoint(Point node);

        /// @brief Deletes a node
        /// @param[in] position The index of the point to remove in m_points
        void DeleteNode(size_t position);

        /// @brief Determines size of the RTree
        [[nodiscard]] size_t Size() const { return m_rtree2D.size(); };

        /// @brief Determines if the RTree is empty
        [[nodiscard]] bool Empty() const { return m_rtree2D.empty(); }

        /// @brief Gets the size of the query
        [[nodiscard]] size_t GetQueryResultSize() const { return m_queryCache.size(); }

        /// @brief Gets the index of a sample in the query
        [[nodiscard]] size_t GetQueryResult(size_t index) const { return m_queryIndices[index]; }

        /// @brief True if a query has results, false otherwise
        [[nodiscard]] bool HasQueryResults() const { return GetQueryResultSize() > 0; }

    private:
        RTree2D m_rtree2D;                                ///< The 2D RTree
        std::vector<std::pair<Point2D, size_t>> m_points; ///< The points
        std::vector<value2D> m_queryCache;                ///< The query cache
        std::vector<size_t> m_queryIndices;               ///< The query indices
        int m_queryVectorCapacity = 100;                  ///< Capacity of the query vector
    };

} // namespace meshkernel
