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

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Entities.hpp"
// #include "MeshKernel/Exceptions.hpp"

// include boost
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <concepts>

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
    template <typename T = bgi::linear<16>>
    class RTree
    {
    public:
        /// @brief Default constructor
        RTree() = default;

        /// @brief Class constructor
        /// @param queryVectorCapacity Optional capacity of the query vector, the default is 100
        RTree(size_t queryVectorCapacity)
            : m_queryVectorCapacity(queryVectorCapacity)
        {
        }

        /// @brief Builds the tree
        /// @tparam T Requires IsCoordinate<T>
        template <std::derived_from<Point> T>
        void BuildTree(std::vector<T> const& nodes)
        {
            m_points.reserve(m_points.size());
            m_points.clear();
            m_rtree.clear();

            for (size_t n = 0; n < nodes.size(); ++n)
            {
                if (nodes[n].IsValid())
                {
                    m_points.emplace_back(Point2D(nodes[n].x, nodes[n].y), n);
                }
            }
            m_rtree = RTree2D(m_points);
        }

        /// @brief Finds all nodes in the search radius and stores the results in the query cache, to be inquired later
        /// @param[in] node The node
        /// @param[in] searchRadiusSquared The squared search radius around the node
        void SearchPoints(Point const& node, double searchRadiusSquared)
        {
            Point2D const nodeSought(node.x, node.y);
            double const searchRadius = std::sqrt(searchRadiusSquared);
            Box2D const box(Point2D(node.x - searchRadius, node.y - searchRadius),
                            Point2D(node.x + searchRadius, node.y + searchRadius));

            m_queryCache.reserve(m_queryVectorCapacity);
            m_queryCache.clear();
            m_rtree.query(bgi::within(box) &&
                              bgi::satisfies([&nodeSought, searchRadiusSquared](Value2D const& v)
                                             { return bg::comparable_distance(v.first, nodeSought) <= searchRadiusSquared; }),
                          std::back_inserter(m_queryCache));

            if (!m_queryCache.empty())
            {
                m_queryIndices.reserve(m_queryCache.size());
                m_queryIndices.clear();
                for (const auto& [first, second] : m_queryCache)
                {
                    m_queryIndices.emplace_back(second);
                }
            }
        }

        /// @brief Gets the nearest of all nodes
        /// @param[in] node The node
        void SearchNearestPoint(Point const& node)
        {
            const Point2D nodeSought(node.x, node.y);

            m_queryCache.reserve(m_queryVectorCapacity);
            m_queryCache.clear();
            m_rtree.query(bgi::nearest(nodeSought, 1),
                          std::back_inserter(m_queryCache));

            if (!m_queryCache.empty())
            {
                m_queryIndices.clear();
                m_queryIndices.emplace_back(m_queryCache.front().second);
            }
        }

        /// @brief Finds the nearest node in the search radius and stores the results in the query cache, to be inquired later
        /// @param[in] node The node
        /// @param[in] searchRadiusSquared The squared search radius around the node
        void SearchNearestPoint(Point const& node, double searchRadiusSquared)
        {
            Point2D const nodeSought(node.x, node.y);
            double const searchRadius = std::sqrt(searchRadiusSquared);
            Box2D const box(Point2D(node.x - searchRadius, node.y - searchRadius),
                            Point2D(node.x + searchRadius, node.y + searchRadius));

            m_queryCache.reserve(m_queryVectorCapacity);
            m_queryCache.clear();
            m_rtree.query(bgi::within(box) &&
                              bgi::satisfies([&nodeSought, searchRadiusSquared](Value2D const& v)
                                             { return bg::comparable_distance(v.first, nodeSought) <= searchRadiusSquared; }) &&
                              bgi::nearest(nodeSought, 1),
                          std::back_inserter(m_queryCache));

            if (!m_queryCache.empty())
            {
                m_queryIndices.clear();
                m_queryIndices.emplace_back(m_queryCache.front().second);
            }
        }

        /// @brief Deletes a node
        /// @param[in] position The index of the point to remove in m_points
        void DeleteNode(size_t position)
        {
            if (m_rtree.remove(m_points[position]) != 1)
            {
                throw std::invalid_argument("Could not remove node at given position.");
                // MeshKernelError("Could not remove node at given position.");
            }
            m_points[position] = {Point2D(constants::missing::doubleValue,
                                          constants::missing::doubleValue),
                                  constants::missing::sizetValue};
        }

        /// @brief Determines size of the RTree
        [[nodiscard]] size_t Size() const { return m_rtree.size(); };

        /// @brief Determines if the RTree is empty
        [[nodiscard]] bool Empty() const { return m_rtree.empty(); }

        /// @brief Gets the size of the query
        [[nodiscard]] size_t GetQueryResultSize() const { return m_queryCache.size(); }

        /// @brief Gets the index of a sample in the query
        [[nodiscard]] size_t GetQueryResult(size_t index) const { return m_queryIndices[index]; }

        /// @brief True if a query has results, false otherwise
        [[nodiscard]] bool HasQueryResults() const { return GetQueryResultSize() > 0; }

    private:
        using Point2D = bg::model::point<double, 2, bg::cs::cartesian>; ///< Typedef for Point2D
        using Box2D = bg::model::box<Point2D>;                          ///< Typedef for box of Point2D
        using Value2D = std::pair<Point2D, size_t>;                     ///< Typedef of pair of Point2D and size_t
        using RTree2D = bgi::rtree<Value2D, T>;                         ///< Typedef for a 2D RTree

        RTree2D m_rtree;                                  ///< The 2D RTree
        std::vector<std::pair<Point2D, size_t>> m_points; ///< The points
        std::vector<Value2D> m_queryCache;                ///< The query cache
        std::vector<size_t> m_queryIndices;               ///< The query indices
        size_t m_queryVectorCapacity = 100;               ///< Capacity of the query vector
    };

} // namespace meshkernel
