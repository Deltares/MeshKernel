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
#include <MeshKernel/Mesh.hpp>
#include <MeshKernel/Network1D.hpp>

/// \namespace meshkernel
/// @brief Contains the logic of the C++ static library
namespace meshkernel
{
    /// @brief A class derived from Mesh, which describes 1d meshes.
    ///
    /// A 1d mesh is composed of a series of connected edges
    /// representing 1d real word features, such as pipes or a sewage network.
    class Mesh1D final : public Mesh
    {

    public:
        /// @brief Default destructor
        ~Mesh1D() override = default;

        /// @brief Default constructor
        Mesh1D() = default;

        /// @brief Construct a mesh1d using only the projection
        /// @param[in] projection The projection to use
        explicit Mesh1D(Projection projection);

        /// @brief Construct a mesh1d starting from the edges and nodes
        /// @param[in] edges The input edges
        /// @param[in] nodes The input nodes
        /// @param[in] projection  The projection to use
        Mesh1D(std::vector<Edge> const& edges,
               std::vector<Point> const& nodes,
               Projection projection);

        /// @brief Constructs a mesh 1d from a network 1d. The network contains the chainages where the discratization points will be computed.
        /// @param[in] network1d The input network
        /// @param[in] minFaceSize The minimum face size below which two nodes will be merged
        Mesh1D(Network1D& network1d, double minFaceSize);

        /// @brief Inquire if a mesh 1d-node is on boundary
        /// @param[in] node The node index
        /// @return If the node is on boundary
        [[nodiscard]] bool IsNodeOnBoundary(UInt node) const { return m_nodesNumEdges[node] == 1; }

        /// @brief Compute a projected node along a line normal to the edges connected to the node.
        /// @param node [in] The node
        /// @param distanceFactor [in] The distance factor
        /// @return The projected node
        [[nodiscard]] Point ComputeProjectedNode(UInt node, double distanceFactor) const;
    };

} // namespace meshkernel
