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

#include <memory>
#include <utility>
#include <vector>

#include "MeshKernel/CurvilinearGrid/CurvilinearGridNodeIndices.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridUtilities.hpp"
#include "MeshKernel/Point.hpp"
#include "MeshKernel/Utilities/LinearAlgebra.hpp"

namespace meshkernel
{

    class CurvilinearGrid;

    class CurvilinearGridBlock
    {
    public:
        CurvilinearGridBlock(const CurvilinearGridNodeIndices& bottomLeft, const CurvilinearGridNodeIndices& topRight);

        // void Extract(const lin_alg::Matrix<Point>& nodes,
        //              const lin_alg::Matrix<NodeType>& nodeTypes,
        //              const CurvilinearGridNodeIndices& start);

        // void Swap(lin_alg::Matrix<Point>& nodes,
        //           lin_alg::Matrix<NodeType>& nodeTypes,
        //           const CurvilinearGridNodeIndices& start);

        void Extract(const CurvilinearGrid& grid);

        void Swap(CurvilinearGrid& grid);

        const CurvilinearGridNodeIndices& StartOffset() const;

        const CurvilinearGridNodeIndices& EndOffset() const;

    private:
        lin_alg::Matrix<Point> m_gridNodes;         ///< Member variable storing the grid
        lin_alg::Matrix<NodeType> m_gridNodesTypes; ///< The grid node types

        CurvilinearGridNodeIndices m_bottomLeft;
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
