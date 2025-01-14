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

#include "MeshKernel/CurvilinearGrid/CurvilinearGridBlock.hpp"
#include "MeshKernel/Constants.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"

meshkernel::CurvilinearGridBlock::CurvilinearGridBlock(const CurvilinearGridNodeIndices& bottomLeft, const CurvilinearGridNodeIndices& topRight)
    : m_bottomLeft(bottomLeft), m_topRight(topRight)
{
    const UInt rows = m_topRight.m_n - m_bottomLeft.m_n;
    const UInt cols = m_topRight.m_m - m_bottomLeft.m_m;

    lin_alg::ResizeAndFillMatrix(m_gridNodes, rows, cols, false, {constants::missing::doubleValue, constants::missing::doubleValue});
}

void meshkernel::CurvilinearGridBlock::CopyFrom(const CurvilinearGrid& grid)
{
    const UInt rows = m_topRight.m_n - m_bottomLeft.m_n;
    const UInt cols = m_topRight.m_m - m_bottomLeft.m_m;

    for (UInt r = 0; r < rows; ++r)
    {
        for (UInt c = 0; c < cols; ++c)
        {
            m_gridNodes(r, c) = grid.GetNode(r + m_bottomLeft.m_n, c + m_bottomLeft.m_m);
        }
    }
}

void meshkernel::CurvilinearGridBlock::Swap(CurvilinearGrid& grid)
{
    const UInt rows = m_topRight.m_n - m_bottomLeft.m_n;
    const UInt cols = m_topRight.m_m - m_bottomLeft.m_m;

    for (UInt r = 0; r < rows; ++r)
    {
        for (UInt c = 0; c < cols; ++c)
        {
            std::swap(m_gridNodes(r, c), grid.GetNode(r + m_bottomLeft.m_n, c + m_bottomLeft.m_m));
        }
    }
}

std::uint64_t meshkernel::CurvilinearGridBlock::MemorySize() const
{
    std::uint64_t result = 0;

    result += sizeof(*this);
    result += static_cast<std::uint64_t>(m_gridNodes.rows() * m_gridNodes.cols() * sizeof(Point));

    return result;
}
