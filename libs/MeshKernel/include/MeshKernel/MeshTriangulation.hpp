//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2024.
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

#include <array>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Operations.hpp"
#include "MeshKernel/Point.hpp"
#include "MeshKernel/Utilities/RTreeFactory.hpp"

namespace meshkernel
{

    /// @brief A simple bounded array
    template <const UInt Dimension>
    class BoundedArray
    {
    public:
        /// @brief Constructor
        BoundedArray() : m_size(0) {}

        /// @brief Number of elements in the array
        UInt size() const
        {
            return m_size;
        }

        /// @brief Add an element to the end of the array
        void push_back(const UInt index)
        {
            if (m_size == Dimension - 1)
            {
                throw ConstraintError("array already at size limit: {}", Dimension);
            }

            // m_indices.push_back(index);
            m_indices[m_size] = index;
            ++m_size;
        }

        /// @brief Get the element at the position
        UInt operator[](const UInt index) const
        {
            return m_indices[index];
        }

        /// @brief Get the element at the position
        UInt& operator[](const UInt index)
        {
            return m_indices[index];
        }

        /// @brief The iterator at the start of the array
        std::array<UInt, Dimension>::const_iterator begin() const
        {
            return m_indices.begin();
        }

        /// @brief The iterator at the end of the array
        std::array<UInt, Dimension>::const_iterator end() const
        {
            return m_indices.begin() + m_size;
        }

        /// @brief Does the array contain the element value or not.
        bool contains(const UInt index) const
        {
            return std::find(begin(), end(), index) != end();
        }

    private:
        /// @brief stack based array containing the values.
        std::array<UInt, Dimension> m_indices;

        /// @brief The current number of elements in the array
        UInt m_size = 0;
    };

    /// @brief Contains a mesh triangulated from a set of points.
    ///
    /// Contains the original set of nodes, the edges connecting nodes
    /// and the set of elements/faces making up the triangulation.
    class MeshTriangulation
    {
    public:
        /// @brief Constructor with array of points
        template <meshkernel::ValidConstPointArray PointVector>
        MeshTriangulation(const PointVector& nodes,
                          const Projection projection);

        /// @brief Constructor with separate arrays of x- and y-coordinates
        template <meshkernel::ValidConstDoubleArray VectorType>
        MeshTriangulation(const VectorType& xNodes,
                          const VectorType& yNodes,
                          const Projection projection);

        /// @brief Get the projection type used in the triangulation
        Projection GetProjection() const;

        /// @brief Get the number of nodes in the triangulation
        UInt NumberOfNodes() const;

        /// @brief Get the number of edges in the triangulation
        UInt NumberOfEdges() const;

        /// @brief Get the number of faces/elements in the triangulation
        UInt NumberOfFaces() const;

        /// @brief Get node as the position
        Point GetNode(const UInt nodeId) const;

        /// @brief Get edge as the position
        Edge GetEdge(const UInt nodeId) const;

        /// @brief Get the node values of the element
        std::array<Point, 3> GetNodes(const UInt faceId) const;

        /// @brief Get the node id's of the element
        std::array<UInt, 3> GetNodeIds(const UInt faceId) const;

        /// @brief Get the edge id's of the element
        std::array<UInt, 3> GetEdgeIds(const UInt faceId) const;

        /// @brief Get the id's of faces either side of the edge.
        ///
        /// May return invalid identifier in one or both values
        const std::array<UInt, 2>& GetFaceIds(const UInt edgeId) const;

        /// @brief Find the nearest face to the point
        UInt FindNearestFace(const Point& pnt) const;

        /// @brief Determine if the point lies within the element
        bool PointIsInElement(const Point& pnt, const UInt faceId) const;

        /// @brief Print the mesh graph
        void Print(std::ostream& out = std::cout) const;

    private:
        static constexpr UInt MaximumNumberOfEdgesPerNode = 16; ///< Maximum number of edges per node

        /// @brief Compute the triangulation.
        void Compute(const std::span<const double>& xNodes,
                     const std::span<const double>& yNodes);

        std::vector<Point> m_nodes;                    ///< x-node values
        std::vector<UInt> m_faceNodes;                 ///< Face nodes flat array passed to the triangulation library
        std::vector<UInt> m_edgeNodes;                 ///< Edge nodes flat array passed to the triangulation library
        std::vector<UInt> m_faceEdges;                 ///< Face edges flat array passed to the triangulation library
        std::vector<std::array<UInt, 2>> m_edgesFaces; ///< edge-face connectivity, generated from triangulation data
        std::vector<std::vector<UInt>> m_nodesEdges;   ///< node-edge connectivity, generated from triangulation data

        std::vector<Point> m_elementCentres; ///< Array of the centres of the elements

        UInt m_numEdges{0}; ///< number of triangulated edges
        UInt m_numFaces{0}; ///< number of triangulated faces

        Projection m_projection = Projection::cartesian; ///< The projection used
        std::unique_ptr<RTreeBase> m_elementCentreRTree; ///< RTree of element centres
        std::unique_ptr<RTreeBase> m_nodeRTree;          ///< RTree of mesh nods
    };

} // namespace meshkernel

inline meshkernel::Projection meshkernel::MeshTriangulation::GetProjection() const
{
    return m_projection;
}

template <meshkernel::ValidConstPointArray PointVector>
meshkernel::MeshTriangulation::MeshTriangulation(const PointVector& nodes,
                                                 const Projection projection)
    : m_nodes(nodes.begin(), nodes.end()),
      m_projection(projection)
{
    if (nodes.size() < constants::geometric::numNodesInTriangle)
    {
        throw ConstraintError("Not enough nodes for a single triangle: {}", nodes.size());
    }

    std::vector<double> xNodes(nodes.size());
    std::vector<double> yNodes(nodes.size());

    std::transform(nodes.begin(), nodes.end(), xNodes.begin(),
                   [](const Point& p)
                   { return p.x; });
    std::transform(nodes.begin(), nodes.end(), yNodes.begin(),
                   [](const Point& p)
                   { return p.y; });

    std::span<double> xNodesSpan(xNodes.data(), xNodes.size());
    std::span<double> yNodesSpan(yNodes.data(), yNodes.size());

    Compute(xNodesSpan, yNodesSpan);
}

template <meshkernel::ValidConstDoubleArray VectorType>
meshkernel::MeshTriangulation::MeshTriangulation(const VectorType& xNodes,
                                                 const VectorType& yNodes,
                                                 const Projection projection)
    : m_nodes(xNodes.size()),
      m_projection(projection)
{
    if (xNodes.size() < constants::geometric::numNodesInTriangle)
    {
        throw ConstraintError("Not enough nodes for a single triangle: {}", xNodes.size());
    }

    if (xNodes.size() != yNodes.size())
    {
        throw ConstraintError("Size mismatch between x- and y-node values: {} /= {}",
                              xNodes.size(), yNodes.size());
    }

    for (size_t i = 0; i < xNodes.size(); ++i)
    {
        m_nodes[i] = Point{xNodes[i], yNodes[i]};
    }

    std::span<const double> xNodesSpan(xNodes.data(), xNodes.size());
    std::span<const double> yNodesSpan(yNodes.data(), yNodes.size());

    Compute(xNodesSpan, yNodesSpan);
}

inline meshkernel::UInt meshkernel::MeshTriangulation::NumberOfNodes() const
{
    return static_cast<UInt>(m_nodes.size());
}

inline meshkernel::UInt meshkernel::MeshTriangulation::NumberOfEdges() const
{
    return m_numEdges;
}

inline meshkernel::UInt meshkernel::MeshTriangulation::NumberOfFaces() const
{
    return m_numFaces;
}

inline meshkernel::Point meshkernel::MeshTriangulation::GetNode(const UInt nodeId) const
{
    if (nodeId == constants::missing::uintValue)
    {
        throw ConstraintError("Invalid node id");
    }

    if (nodeId >= NumberOfNodes())
    {
        throw ConstraintError("node id out of range: {} >= {}", nodeId, NumberOfNodes());
    }

    return m_nodes[nodeId];
}

inline meshkernel::Edge meshkernel::MeshTriangulation::GetEdge(const UInt edgeId) const
{
    if (edgeId == constants::missing::uintValue)
    {
        throw ConstraintError("Invalid edge id");
    }

    if (edgeId >= m_numEdges)
    {
        throw ConstraintError("edge id out of range: {} >= {}", edgeId, m_numEdges);
    }

    return {m_edgeNodes[2 * edgeId], m_edgeNodes[2 * edgeId + 1]};
}

inline const std::array<meshkernel::UInt, 2>& meshkernel::MeshTriangulation::GetFaceIds(const UInt edgeId) const
{
    if (edgeId == constants::missing::uintValue)
    {
        throw ConstraintError("Invalid edge id");
    }

    if (edgeId >= m_numEdges)
    {
        throw ConstraintError("edge id out of range: {} >= {}", edgeId, m_numEdges);
    }

    return m_edgesFaces[edgeId];
}
