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
#include <vector>

#include "MeshKernel/CurvilinearGrid/CurvilinearGridNodeIndices.hpp"
#include "MeshKernel/Point.hpp"
#include "MeshKernel/UndoActions/BaseMeshUndoAction.hpp"
#include <MeshKernel/Utilities/LinearAlgebra.hpp>

namespace meshkernel
{
    class CurvilinearGrid;

    class CurvilinearGridRefinementUndoAction : public BaseMeshUndoAction<CurvilinearGridRefinementUndoAction, CurvilinearGrid>
    {
    public:
        static std::unique_ptr<CurvilinearGridRefinementUndoAction> Create(CurvilinearGrid& grid);

        CurvilinearGridRefinementUndoAction(CurvilinearGrid& grid);

        void Swap(lin_alg::Matrix<Point>& nodes, CurvilinearGridNodeIndices& lower, CurvilinearGridNodeIndices& upper);

    private:
        lin_alg::Matrix<Point> m_nodes; ///< Grid nodes
        CurvilinearGridNodeIndices m_lower;
        CurvilinearGridNodeIndices m_upper;
    };

} // namespace meshkernel
