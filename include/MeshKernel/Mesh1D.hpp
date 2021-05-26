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
#include <vector>

#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh.hpp>

/// \namespace meshkernel
/// @brief Contains the logic of the C++ static library
namespace meshkernel
{
    /// @brief A class derived from Mesh, which describes 1d meshes.
    ///
    /// A 1d mesh is composed of a series of connected edges
    /// representing 1d real word features, such as pipes or a sewage network.
    class Mesh1D : public Mesh
    {
    public:
        /// @brief Default constructor
        Mesh1D() = default;

        /// @brief Construct a mesh1d starting from the edges and nodes
        /// @param[in] edges The input edges
        /// @param[in] nodes The input nodes
        /// @param[in] projection  The projection to use
        explicit Mesh1D(std::vector<Edge> const& edges,
                        std::vector<Point> const& nodes,
                        Projection projection);

        /// @brief Construct a mesh1d by discretizing polyLines
        /// @param[in] polyLines The polylines to be discretize
        /// @param[in] fixedChainages The fixed chainages. These are locations where the first discretization points before and after should be at a distance equal to \ref offsetFixedChainages.
        /// @param[in] offset The regular offset between points
        /// @param[in] minFaceSize The minimum face size. The distance between two discratization point must be no less than this length.
        /// @param[in] offsetFixedChainages The offset to use for fixed chainages
        /// @param[in] projection The projection to use
        explicit Mesh1D(std::vector<std::vector<Point>> const& polyLines,
                        std::vector<std::vector<double>> const& fixedChainages,
                        double offset,
                        double minFaceSize,
                        double offsetFixedChainages,
                        Projection projection);

        /// @brief Inquire if a mesh 1d-node is on boundary
        /// @param[in] node The node index
        /// @return If the node is on boundary
        [[nodiscard]] bool IsNodeOnBoundary(size_t node) const { return m_nodesNumEdges[node] == 1; }

        /// @brief Compute the chainages from fixed point locations
        /// @param[in] fixedChainages The fixed chainages. These are locations where the first discretization points before and after should be at a distance equal to \ref offsetFixedChainages.
        /// @param[in] minFaceSize  The minimum face size. The distance between two discratization point must be no less than this length.
        /// @param[in] offsetFromFixedChainages  The offset to use for fixed chainages
        /// @param[in, out] chainages The vector containing the computed chainages
        void ComputeFixedChainages(std::vector<double> const& fixedChainages,
                                   double minFaceSize,
                                   double offsetFromFixedChainages,
                                   std::vector<double>& chainages);

        /// @brief Compute the chainages at a regular offset
        /// @param[in] offset The regular offset between points
        /// @param[in, out] chainages The vector containing the computed chainages
        void ComputeChainagesAtOffset(double offset, std::vector<double>& chainages);
    };

} // namespace meshkernel
