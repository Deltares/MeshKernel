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

#include "MeshKernel/BoundingBox.hpp"
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Exceptions.hpp>

// r-tree
// https://gist.github.com/logc/10272165

namespace meshkernel
{

    class RTreeBase
    {
    public:
        virtual ~RTreeBase() = default;

        virtual void BuildTree(const std::vector<Point>& nodes) = 0;

        virtual void BuildTree(const std::vector<Sample>& samples) = 0;

        virtual void BuildTree(const std::vector<Point>& nodes, const BoundingBox& boundingBox) = 0;

        virtual void BuildTree(const std::vector<Sample>& samples, const BoundingBox& boundingBox) = 0;

        /// @brief Finds all nodes in the search radius and stores the results in the query cache, to be inquired later
        /// @param[in] node The node
        /// @param[in] searchRadiusSquared The squared search radius around the node
        virtual void SearchPoints(Point const& node, double searchRadiusSquared) = 0;

        /// @brief Finds the nearest node in the search radius and stores the results in the query cache, to be inquired later
        /// @param[in] node The node
        /// @param[in] searchRadiusSquared The squared search radius around the node
        virtual void SearchNearestPoint(Point const& node, double searchRadiusSquared) = 0;

        /// @brief Gets the nearest of all nodes
        /// @param[in] node The node
        virtual void SearchNearestPoint(Point const& node) = 0;

        /// @brief Deletes a node
        /// @param[in] position The index of the point to remove in m_points
        virtual void DeleteNode(UInt position) = 0;

        /// @brief Determines size of the RTree
        virtual UInt Size() const = 0;

        /// @brief Determines if the RTree is empty
        virtual bool Empty() const = 0;

        /// @brief Gets the size of the query
        virtual UInt GetQueryResultSize() const = 0;

        /// @brief Gets the index of a sample in the query
        virtual UInt GetQueryResult(UInt index) const = 0;

        /// @brief True if a query has results, false otherwise
        virtual bool HasQueryResults() const = 0;
    };

} // namespace meshkernel
