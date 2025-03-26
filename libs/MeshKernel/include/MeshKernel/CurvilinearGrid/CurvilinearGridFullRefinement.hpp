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

#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"
#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Splines.hpp"
#include "MeshKernel/UndoActions/UndoAction.hpp"

namespace meshkernel
{
    /// @brief Refines the entire curvilinear grid.
    ///
    /// This refinement can be different in each m- and n-direction.
    /// Refinement factor of 2 indicates split the element into two halves
    /// Refinement factor of 4 indicates split the element into four quarters
    /// Refinement factor of 1 indicates do not split the element.
    class CurvilinearGridFullRefinement final
    {
    public:
        /// @brief Refine the entire grid by positive or negative refinement factors.
        ///
        /// @param grid [in out] The grid to be refined
        /// @param mRefinement [in] How much refinement (positive) or de-refinement (negative) is required in the m-direction
        /// @param nRefinement [in] How much refinement (positive) or de-refinement (negative) is required in the n-direction
        /// @returns the unto action, should the mesh want to be reverted to the original state.
        [[nodiscard]] UndoActionPtr Compute(CurvilinearGrid& grid,
                                            const int mRefinement,
                                            const int nRefinement) const;

    private:
        /// @brief Refine the entire grid by positive refinement factors.
        ///
        /// @param grid [in out] The grid to be refined
        /// @param mRefinement [in] How much refinement required in the m-direction
        /// @param nRefinement [in] How much refinement required in the n-direction
        /// @returns the unto action, should the mesh want to be reverted to the original state.
        [[nodiscard]] UndoActionPtr ComputeRefinement(CurvilinearGrid& grid,
                                                      const UInt mRefinement,
                                                      const UInt nRefinement) const;
        /// @brief Compute the discretisation of two opposite edges of a face
        ///
        /// This can be for either top and bottom or left and right edges.
        void ComputeRefinedElementEdges(const Splines& splines,
                                        const UInt splineIndex,
                                        const UInt currentIndex,
                                        const UInt refinement,
                                        std::vector<Point>& elementSide1,
                                        std::vector<Point>& elementSide2) const;

        /// @brief Determine if the all nodes of a face are valid
        bool ValidFace(const CurvilinearGrid& grid,
                       const UInt m, const UInt n) const;
    };

} // namespace meshkernel
