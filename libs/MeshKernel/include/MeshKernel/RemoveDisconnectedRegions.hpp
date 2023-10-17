//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2023.
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
#include <utility>
#include <vector>

#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Mesh2D.hpp"

namespace meshkernel
{
    /// @brief Removes any small regions that are disconnected from the main part of the domain.
    ///
    /// @note The assumption is that the main region of interest will contain the largest number of elements.
    /// The element of any other disconnected region (having fewer elements) will be removed.
    /// Element connectivity is across shared edges, so elements that share only a single corner node are not
    /// considered to be connected.
    class RemoveDisconnectedRegions final
    {
    public:
        /// @brief Remove the disconnected regions.
        void Compute(Mesh2D& mesh) const;

    private:
        /// @brief Get the neighbour element along an edge.
        ///
        /// Returns the null value if elementId is the null value or does not match any element connected to the edge.
        UInt GetNeighbour(const std::array<UInt, 2>& edge, const UInt elementId) const;

        /// @brief Label the elements in the connected (across faces only) region
        ///
        /// The faces of the region will be labeled with the regionId.
        /// The connected region is labeled using a flood fill algorithm
        void LabelConnectedRegion(const Mesh2D& mesh, const UInt regionId, std::vector<UInt>& elementRegionId, const UInt unlabledElementId, UInt& elementCount) const;

        /// @brief Label the elements in a single connected region.
        ///
        /// If no elements are left unlabeled then the elementCount is zero
        void LabelSingleDomainRegion(const Mesh2D& mesh, const UInt regionId, std::vector<UInt>& elementRegionId, UInt& elementCount) const;

        /// @brief Label the elements of all regions in the mesh.
        ///
        /// Each region will be assigned a unique identifier
        /// All elements in a single region will be assigned the same unique identifier (each region will have a different identifier)
        void LabelAllDomainRegions(const Mesh2D& mesh, std::vector<UInt>& elementRegionId, std::vector<std::pair<UInt, UInt>>& regionCount) const;

        /// @brief Remove elements from regions that do not have the main region identifier.
        void RemoveDetachedRegions(Mesh2D& mesh, const UInt regionId, std::vector<UInt>& elementRegionId, UInt& numberOfElementsRemoved) const;
    };

} // namespace meshkernel
