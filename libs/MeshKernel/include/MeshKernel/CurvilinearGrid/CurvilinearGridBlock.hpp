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

#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

#include "MeshKernel/CurvilinearGrid/CurvilinearGridNodeIndices.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridUtilities.hpp"
#include "MeshKernel/Point.hpp"
#include "MeshKernel/Utilities/LinearAlgebra.hpp"

namespace meshkernel
{

    /// @brief Forward declaration of the curvilinear mesh
    class CurvilinearGrid;

    /// @brief Block of grid node with start and end offsets.
    class CurvilinearGridBlock
    {
    public:
        /// @brief Constructor
        CurvilinearGridBlock(const CurvilinearGridNodeIndices& bottomLeft, const CurvilinearGridNodeIndices& topRight);

        /// @brief Copy block of nodes from the curvilinear grid.
        void CopyFrom(const CurvilinearGrid& grid);

        /// @brief Swap the saved grid nodes with those from the mesh.
        void Swap(CurvilinearGrid& grid);

        /// @brief Get the start offset indices
        const CurvilinearGridNodeIndices& StartOffset() const;

        /// @brief Get the end offset indices
        const CurvilinearGridNodeIndices& EndOffset() const;

        /// @brief Get the approximate number of bytes used.
        std::uint64_t MemorySize() const;

    private:
        /// @brief Member variable storing the grid nodes
        lin_alg::Matrix<Point> m_gridNodes;

        /// @brief Bottom left indices of the block
        CurvilinearGridNodeIndices m_bottomLeft;

        /// @brief Top right indices of the block (1 past end)
        CurvilinearGridNodeIndices m_topRight;
    };

} // namespace meshkernel

inline const meshkernel::CurvilinearGridNodeIndices& meshkernel::CurvilinearGridBlock::StartOffset() const
{
    return m_bottomLeft;
}

inline const meshkernel::CurvilinearGridNodeIndices& meshkernel::CurvilinearGridBlock::EndOffset() const
{
    return m_topRight;
}
