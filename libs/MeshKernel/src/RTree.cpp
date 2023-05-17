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

#include "MeshKernel/RTree.hpp"
#include "MeshKernel/Exceptions.hpp"

#include <cmath>

namespace meshkernel
{

    void RTree::SearchPoints(Point const& node, double searchRadiusSquared)
    {
        Point2D const nodeSought(node.x, node.y);
        Box2D const box = MakeBox(nodeSought, std::sqrt(searchRadiusSquared));

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

    void RTree::SearchNearestPoint(Point const& node)
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

    void RTree::SearchNearestPoint(Point const& node, double searchRadiusSquared)
    {
        Point2D const nodeSought(node.x, node.y);
        Box2D const box = MakeBox(nodeSought, std::sqrt(searchRadiusSquared));

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

    void RTree::DeleteNode(size_t position)
    {
        if (m_rtree.remove(m_points[position]) != 1)
        {
            throw MeshKernelError(VariadicErrorMessage("Could not remove node at position {}.", position));
        }
        m_points[position] = {Point2D(constants::missing::doubleValue,
                                      constants::missing::doubleValue),
                              constants::missing::sizetValue};
    }

} // namespace meshkernel