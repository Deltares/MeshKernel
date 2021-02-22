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

#pragma once

#include <memory>

#include <MeshKernel/CurvilinearGrid.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Splines.hpp>

namespace meshkernel
{
    class CurvilinearGrid;
    /// @brief A class implementing the curvilinear grid refinement algorithm
    class CurvilinearGridRefinement
    {
    public:
        CurvilinearGridRefinement(const std::shared_ptr<CurvilinearGrid>& m_grid, const Point& m_first_point, const Point& m_second_point, size_t m_refinement, size_t n_refinement)
            : m_grid(m_grid),
              m_firstPoint(m_first_point),
              m_secondPoint(m_second_point),
              m_refinement(m_refinement),
              n_refinement(n_refinement)
        {
        }

        void Compute()
        {
            auto [mFirstNode, nFirstNode] = m_grid->GetNodeIndices(m_firstPoint);
            auto [mSecondNode, nSecondNode] = m_grid->GetNodeIndices(m_secondPoint);

            if (mSecondNode - mFirstNode != 0 || nSecondNode - nFirstNode != 0)
            {
                throw std::invalid_argument("CurvilinearGridRefinement::Compute: The selected curvilinear grid nodes are not on the same grid-line");
            }

            const auto numMToRefine = mSecondNode - mFirstNode;
            const auto numNToRefine = nSecondNode - nFirstNode;

            size_t m = m_grid->m_numM - numMToRefine + numMToRefine * m_refinement;
            size_t n = m_grid->m_numN - numNToRefine + numNToRefine * n_refinement;
        }

    private:
        void ComputeSplineDerivates()
        {

            // second order derivates horizontal grid lines
            std::vector<Point> horizontalPoints(m_grid->m_numM);
            std::vector<std::vector<Point>> horizontalPointsDerivates(m_grid->m_numN);
            horizontalPoints.clear();
            for (auto n = 0; n < m_grid->m_numN; ++n)
            {
                for (auto m = 0; m < m_grid->m_numM; ++m)
                {
                    horizontalPoints.emplace_back(m_grid->m_nodes[m][n]);
                }
                horizontalPoints.clear();
                horizontalPointsDerivates[n] = Splines::SecondOrderDerivative(horizontalPoints);
            }

            // second order derivates vertical grid lines
            for (auto m = 0; m < m_grid->m_numM; ++m)
            {
                for (auto n = 0; n < m_grid->m_numN; ++n)
                {
                }
            }
        }

        std::shared_ptr<CurvilinearGrid> m_grid; ///< A pointer to mesh
        Point m_firstPoint;
        Point m_secondPoint;
        size_t m_refinement;
        size_t n_refinement;
    }
} // namespace meshkernel