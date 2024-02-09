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

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridRefinement.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Operations.hpp>

using meshkernel::CurvilinearGrid;
using meshkernel::CurvilinearGridRefinement;

CurvilinearGridRefinement::CurvilinearGridRefinement(CurvilinearGrid& grid, UInt refinement)
    : CurvilinearGridAlgorithm(grid),
      m_refinement(refinement)
{
    // Store the m and n grid lines in separate vectors, and compute the spline derivatives for each gridline
    m_splines = Splines(m_grid);
}

void CurvilinearGridRefinement::Compute()
{

    if (!m_lowerLeft.IsValid() || !m_upperRight.IsValid())
    {
        throw std::invalid_argument("CurvilinearGridRefinement::Compute: lower left and upper right corners defining the curvilinear grid block are not set");
    }

    // The points must lie on the same grid line
    if (m_lowerLeft.m_m != m_upperRight.m_m && m_lowerLeft.m_n != m_upperRight.m_n)
    {
        throw std::invalid_argument("CurvilinearGridRefinement::Compute: The selected curvilinear grid nodes are not on the same grid-line");
    }

    // Estimate the dimension of the refined grid
    const auto numMToRefine = m_upperRight.m_m - m_lowerLeft.m_m;
    const auto numNToRefine = m_upperRight.m_n - m_lowerLeft.m_n;
    const UInt maxM = m_grid.NumM() + numMToRefine * (m_refinement - 1);
    const UInt maxN = m_grid.NumN() + numNToRefine * (m_refinement - 1);

    // Local vector for each curvilinear grid face
    std::vector<Point> bottomRefinement(m_refinement);
    std::vector<Point> topRefinement(m_refinement);
    std::vector<Point> leftRefinement(m_refinement);
    std::vector<Point> rightRefinement(m_refinement);

    // The refined grid
    lin_alg::Matrix<Point> refinedGrid(maxM, maxN);

    UInt refinedM = 0;
    for (UInt currentM = 0; currentM < m_grid.NumM() - 1; ++currentM)
    {
        UInt localMRefinement = 1;
        if (currentM >= m_lowerLeft.m_m && currentM < m_upperRight.m_m)
        {
            localMRefinement = m_refinement;
        }

        UInt refinedN = 0;
        for (UInt currentN = 0; currentN < m_grid.NumN() - 1; ++currentN)
        {

            UInt localNRefinement = 1;
            if (currentN >= m_lowerLeft.m_n && currentN < m_upperRight.m_n)
            {
                localNRefinement = m_refinement;
            }

            // Only if all grid nodes of the face are valid, perform transfinite interpolation
            if (m_grid.GetNode(currentM, currentN).IsValid() &&
                m_grid.GetNode(currentM + 1, currentN).IsValid() &&
                m_grid.GetNode(currentM, currentN + 1).IsValid() &&
                m_grid.GetNode(currentM + 1, currentN + 1).IsValid())
            {
                // Calculate m-direction spline points
                bottomRefinement.clear();
                topRefinement.clear();
                for (UInt m = 0; m < localMRefinement + 1; ++m)
                {
                    const auto splineIndex = currentN;
                    const auto interpolationPoint = static_cast<double>(currentM) + static_cast<double>(m) / static_cast<double>(localMRefinement);
                    bottomRefinement.emplace_back(ComputePointOnSplineAtAdimensionalDistance(m_splines.m_splineNodes[splineIndex], m_splines.m_splineDerivatives[splineIndex], interpolationPoint));
                    topRefinement.emplace_back(ComputePointOnSplineAtAdimensionalDistance(m_splines.m_splineNodes[splineIndex + 1], m_splines.m_splineDerivatives[splineIndex + 1], interpolationPoint));
                }

                // Calculate m-direction spline points
                leftRefinement.clear();
                rightRefinement.clear();
                for (UInt n = 0; n < localNRefinement + 1; ++n)
                {
                    const auto splineIndex = m_grid.NumN() + currentM;
                    const auto interpolationPoint = static_cast<double>(currentN) + static_cast<double>(n) / static_cast<double>(localNRefinement);
                    leftRefinement.emplace_back(ComputePointOnSplineAtAdimensionalDistance(m_splines.m_splineNodes[splineIndex], m_splines.m_splineDerivatives[splineIndex], interpolationPoint));
                    rightRefinement.emplace_back(ComputePointOnSplineAtAdimensionalDistance(m_splines.m_splineNodes[splineIndex + 1], m_splines.m_splineDerivatives[splineIndex + 1], interpolationPoint));
                }

                // Perform transfinite interpolation on the current curvilinear face
                const auto localGrid = DiscretizeTransfinite(leftRefinement,
                                                             rightRefinement,
                                                             bottomRefinement,
                                                             topRefinement,
                                                             m_grid.m_projection,
                                                             localMRefinement,
                                                             localNRefinement);
                // Copy the local grid into the refined grid
                for (Eigen::Index m = 0; m < localMRefinement + 1; ++m)
                {
                    for (Eigen::Index n = 0; n < localNRefinement + 1; ++n)
                    {
                        refinedGrid(refinedM + m, refinedN + n) = localGrid(m, n);
                    }
                }
            }
            refinedN += localNRefinement;
        }
        refinedM += localMRefinement;
    }

    // Substitute original grid with the refined one
    m_grid.SetGridNodes(refinedGrid);
}
