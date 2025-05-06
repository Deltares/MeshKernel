//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2025.
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

#include "MeshKernel/BoundingBox.hpp"
#include "MeshKernel/Utilities/RTreeBase.hpp"

#include <concepts>
#include <utility>

// r-tree
// https://gist.github.com/logc/10272165

namespace meshkernel
{
    // using a namespace alias
    namespace bg = boost::geometry;
    namespace bgi = boost::geometry::index;

    /// @brief Class used for inquiring adjacent nodes in a mesh.
    ///
    /// The RTree class is designed for querying adjacent nodes within a mesh.
    /// It encapsulates the boost::geometry::index::rtree functionality and extends it with a user-friendly
    /// interface to perform common spatial queries. This class is templated on the projection type,
    /// allowing flexibility in the geometric coordinate system used (default is bg::cs::cartesian).
    ///
    /// @tparam projection The geometric coordinate system projection (default is bg::cs::cartesian).
    ///
    /// The RTree class is primarily utilized within the mesh library for efficiently querying
    /// the closest mesh nodes and edges to a specified point. It employs the RTreeBase class as its base.
    ///
    /// Internally, the RTree class maintains a query cache (`m_queryCache`) which is a vector used
    /// to collect and store query results. This design helps optimize performance by avoiding frequent
    /// reallocations when the number of results changes between queries.
    ///
    /// Example usage:
    /// @code
    /// // Create an RTree instance with the default Cartesian projection.
    /// RTree<> cartesianRTree;
    ///
    /// // Perform a search for the nearest neighbors within a specified distance.
    /// auto resultPoints = cartesianRTree.SearchPoints(queryPoint, searchDistance);
    ///
    /// // Perform a search for the single nearest neighbor.
    /// auto nearestPoint = cartesianRTree.SearchNearestPoint(queryPoint);
    /// @endcode
    ///
    /// Note: For advanced use cases and different geometric coordinate systems, users can provide
    /// a custom projection template parameter when instantiating the RTree class.
    ///
    /// For more details on available query methods, refer to the base class documentation: meshkernel::RTreeBase.
    class RTreeSphericalToCartesian : public RTreeBase
    {
        using Point3D = bg::model::point<double, 3, bg::cs::cartesian>; ///< Typedef for Point3D
        using Box3D = bg::model::box<Point3D>;                          ///< Typedef for box of Point3D
        using Value3D = std::pair<Point3D, UInt>;                       ///< Typedef of pair of Point3D and UInt
        using RTree3D = bgi::rtree<Value3D, bgi::linear<16>>;           ///< Typedef for a 3D RTree

        /// @brief Ninety degrees
        static constexpr double NinetyDegrees = 90.0;

    public:
        /// @brief Builds the tree from a vector of Points
        /// @param[in] nodes The vector of nodes
        void BuildTree(const std::vector<Point>& nodes) override
        {
            auto conversion = [](const Point& p)
            { return convert(p); };
            BuildTreeFromVector(nodes, m_points3D, conversion);
            m_rtree3D = RTree3D(m_points3D);
        }

        /// @brief Builds the tree from a vector of samples
        /// @param[in] samples The vector of samples
        void BuildTree(const std::vector<Sample>& samples) override
        {
            auto conversion = [](const Point& p)
            { return convert(p); };
            BuildTreeFromVector(samples, m_points3D, conversion);
            m_rtree3D = RTree3D(m_points3D);
        }

        /// @brief Builds the tree from a vector of points within a bounding box
        /// @param[in] nodes The vector of nodes
        /// @param[in] boundingBox The vector bounding box
        void BuildTree(const std::vector<Point>& nodes, const BoundingBox& boundingBox) override
        {
            auto conversion = [](const Point& p)
            { return convert(p); };
            BuildTreeFromVectorWithinBoundingBox(nodes, m_points3D, conversion, boundingBox);
            m_rtree3D = RTree3D(m_points3D);
        }

        /// @brief Builds the tree from a vector of samples within a bounding box
        /// @param[in] samples The vector of samples
        /// @param[in] boundingBox The vector bounding box
        void BuildTree(const std::vector<Sample>& samples, const BoundingBox& boundingBox) override
        {
            auto conversion = [](const Point& p)
            { return convert(p); };
            BuildTreeFromVectorWithinBoundingBox(samples, m_points3D, conversion, boundingBox);
            m_rtree3D = RTree3D(m_points3D);
        }

        /// @brief Finds all nodes in the search radius and stores the results in the query cache, to be inquired later
        /// @param[in] node The node
        /// @param[in] searchRadiusSquared The squared search radius around the node
        void SearchPoints(Point const& node, double searchRadiusSquared) override;

        /// @brief Finds the nearest node in the search radius and stores the results in the query cache, to be inquired later
        /// @param[in] node The node
        /// @param[in] searchRadiusSquared The squared search radius around the node
        void SearchNearestPoint(Point const& node, double searchRadiusSquared) override;

        /// @brief Gets the nearest of all nodes
        /// @param[in] node The node
        void SearchNearestPoint(Point const& node) override;

        /// @brief Deletes a node
        /// @param[in] position The index of the point to remove in m_points3D
        void DeleteNode(UInt position) override;

        /// @brief Determines size of the RTree
        [[nodiscard]] UInt Size() const override { return static_cast<UInt>(m_rtree3D.size()); };

        /// @brief Determines if the RTree is empty
        [[nodiscard]] bool Empty() const override { return m_rtree3D.empty(); }

        /// @brief Gets the size of the query
        [[nodiscard]] UInt GetQueryResultSize() const override { return static_cast<UInt>(m_queryCache.size()); }

        /// @brief Gets the index of a sample in the query
        [[nodiscard]] UInt GetQueryResult(UInt index) const override { return m_queryIndices[index]; }

        /// @brief True if a query has results, false otherwise
        [[nodiscard]] bool HasQueryResults() const override { return !m_queryCache.empty(); }

    private:
        /// @brief Convert 2d point in spherical coordinates to 3d point in Cartesian coordinates.
        static Point3D convert(const Point& p);

        /// @brief Performs a spatial search within a search radius
        /// @param[in] node The reference point for the search.
        /// @param[in] searchRadiusSquared The squared search radius.
        /// @param[in] findNearest If true, finds the nearest point; otherwise, finds all points within the radius.
        void Search(Point const& node, double searchRadiusSquared, bool findNearest);

        RTree3D m_rtree3D;                              ///< The 3D RTree
        std::vector<std::pair<Point3D, UInt>> m_points3D; ///< The points
        std::vector<Value3D> m_queryCache;              ///< The query cache
        std::vector<UInt> m_queryIndices;               ///< The query indices
        UInt m_queryVectorCapacity = 100;               ///< Capacity of the query vector
    };

} // namespace meshkernel
