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
#include <MeshKernel/Exceptions.hpp>

// include boost
#define BOOST_ALLOW_DEPRECATED_HEADERS
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>
#undef BOOST_ALLOW_DEPRECATED_HEADERS

#include "MeshKernel/BoundingBox.hpp"
#include "MeshKernel/Utilities/RTreeBase.hpp"

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
    template <typename projection = bg::cs::cartesian>
    class RTree : public RTreeBase
    {
        using Point2D = bg::model::point<double, 2, projection>; ///< Typedef for Point2D
        using Box2D = bg::model::box<Point2D>;                   ///< Typedef for box of Point2D
        using Value2D = std::pair<Point2D, UInt>;                ///< Typedef of pair of Point2D and size_t
        using RTree2D = bgi::rtree<Value2D, bgi::linear<16>>;    ///< Typedef for a 2D RTree

    public:
        void BuildTree(const std::vector<Point>& nodes) override
        {
            BuildTreeFromVector(nodes);
        }

        void BuildTree(const std::vector<Sample>& samples) override
        {
            BuildTreeFromVector(samples);
        }

        void BuildTree(const std::vector<Point>& nodes, const BoundingBox& boundingBox) override
        {
            BuildTreeFromVectorWithBoundingBox(nodes, boundingBox);
        }

        void BuildTree(const std::vector<Sample>& samples, const BoundingBox& boundingBox) override
        {
            BuildTreeFromVectorWithBoundingBox(samples, boundingBox);
        }

        /// @brief Finds all nodes in the search radius and stores the results in the query cache, to be inquired later
        /// @param[in] node The node
        /// @param[in] searchRadiusSquared The squared search radius around the node
        void SearchPoints(Point const& node, double searchRadiusSquared) override
        {
            if (Empty())
            {
                throw AlgorithmError("RTree is empty, search cannot be performed");
            }

            const auto searchRadius = std::sqrt(searchRadiusSquared);

            Box2D const box(Point2D(node.x - searchRadius, node.y - searchRadius),
                            Point2D(node.x + searchRadius, node.y + searchRadius));
            Point2D nodeSought = Point2D(node.x, node.y);

            m_queryCache.reserve(m_queryVectorCapacity);
            m_queryCache.clear();
            m_rtree2D.query(bgi::within(box) &&
                                bgi::satisfies([&nodeSought, &searchRadiusSquared](Value2D const& v)
                                               { return bg::comparable_distance(v.first, nodeSought) <= searchRadiusSquared; }),
                            std::back_inserter(m_queryCache));

            m_queryIndices.reserve(m_queryCache.size());
            m_queryIndices.clear();

            for (size_t i = 0; i < m_queryCache.size(); ++i)
            {
                auto const index = std::get<1>(m_queryCache[i]);
                m_queryIndices.emplace_back(index);
            }
        }

        /// @brief Finds the nearest node in the search radius and stores the results in the query cache, to be inquired later
        /// @param[in] node The node
        /// @param[in] searchRadiusSquared The squared search radius around the node
        void SearchNearestPoint(Point const& node, double searchRadiusSquared) override
        {
            if (Empty())
            {
                throw AlgorithmError("RTree is empty, search cannot be performed");
            }

            m_queryCache.reserve(m_queryVectorCapacity);
            m_queryCache.clear();
            const Point2D nodeSought = Point2D(node.x, node.y);
            const auto searchRadius = std::sqrt(searchRadiusSquared);
            Box2D const box(Point2D(node.x - searchRadius, node.y - searchRadius),
                            Point2D(node.x + searchRadius, node.y + searchRadius));
            m_rtree2D.query(bgi::within(box) &&
                                bgi::satisfies([&nodeSought, &searchRadiusSquared](Value2D const& v)
                                               { return bg::comparable_distance(v.first, nodeSought) <= searchRadiusSquared; }) &&
                                bgi::nearest(nodeSought, 1),
                            std::back_inserter(m_queryCache));

            if (!m_queryCache.empty())
            {
                m_queryIndices.clear();
                m_queryIndices.emplace_back(m_queryCache[0].second);
            }
        }

        /// @brief Gets the nearest of all nodes
        /// @param[in] node The node
        void SearchNearestPoint(Point const& node) override
        {
            if (Empty())
            {
                throw AlgorithmError("RTree is empty, search cannot be performed");
            }

            m_queryCache.reserve(m_queryVectorCapacity);
            m_queryCache.clear();
            const Point2D nodeSought = Point2D(node.x, node.y);
            m_rtree2D.query(bgi::nearest(nodeSought, 1), std::back_inserter(m_queryCache));

            if (!m_queryCache.empty())
            {
                m_queryIndices.clear();
                m_queryIndices.emplace_back(m_queryCache[0].second);
            }
        }

        /// @brief Deletes a node
        /// @param[in] position The index of the point to remove in m_points
        void DeleteNode(UInt position) override
        {
            if (Empty())
            {
                throw AlgorithmError("RTree is empty, deletion cannot performed");
            }

            const auto numberRemoved = m_rtree2D.remove(m_points[position]);
            if (numberRemoved != 1)
            {
                return;
            }
            m_points[position] = {Point2D{constants::missing::doubleValue, constants::missing::doubleValue}, std::numeric_limits<UInt>::max()};
        }

        /// @brief Determines size of the RTree
        [[nodiscard]] UInt Size() const override { return static_cast<UInt>(m_rtree2D.size()); };

        /// @brief Determines if the RTree is empty
        [[nodiscard]] bool Empty() const override { return m_rtree2D.empty(); }

        /// @brief Gets the size of the query
        [[nodiscard]] UInt GetQueryResultSize() const override { return static_cast<UInt>(m_queryCache.size()); }

        /// @brief Gets the index of a sample in the query
        [[nodiscard]] UInt GetQueryResult(UInt index) const override { return m_queryIndices[index]; }

        /// @brief True if a query has results, false otherwise
        [[nodiscard]] bool HasQueryResults() const override { return !m_queryCache.empty(); }

    private:
        /// @brief Builds the tree
        /// @param[in] nodes The nodes
        template <std::derived_from<Point> T>
        void BuildTreeFromVector(const std::vector<T>& nodes)
        {
            m_points.reserve(m_points.size());
            m_points.clear();
            m_rtree2D.clear();

            for (UInt n = 0; n < nodes.size(); ++n)
            {
                if (nodes[n].x != constants::missing::doubleValue && nodes[n].y != constants::missing::doubleValue)
                {
                    m_points.emplace_back(Point2D{nodes[n].x, nodes[n].y}, n);
                }
            }
            m_rtree2D = RTree2D(m_points);
        }

        /// @brief Builds the tree with nodes
        /// @param[in] nodes The nodes
        template <std::derived_from<Point> T>
        void BuildTreeFromVectorWithBoundingBox(const std::vector<T>& nodes, const BoundingBox& boundingBox)
        {
            m_points.reserve(m_points.size());
            m_points.clear();
            m_rtree2D.clear();

            for (UInt n = 0; n < nodes.size(); ++n)
            {
                if (!boundingBox.Contains(nodes[n]))
                {
                    continue;
                }

                if (nodes[n].x != constants::missing::doubleValue && nodes[n].y != constants::missing::doubleValue)
                {
                    m_points.emplace_back(Point2D{nodes[n].x, nodes[n].y}, n);
                }
            }
            m_rtree2D = RTree2D(m_points);
        }

        RTree2D m_rtree2D;                              ///< The 2D RTree
        std::vector<std::pair<Point2D, UInt>> m_points; ///< The points
        std::vector<Value2D> m_queryCache;              ///< The query cache
        std::vector<UInt> m_queryIndices;               ///< The query indices
        UInt m_queryVectorCapacity = 100;               ///< Capacity of the query vector
    };
} // namespace meshkernel
