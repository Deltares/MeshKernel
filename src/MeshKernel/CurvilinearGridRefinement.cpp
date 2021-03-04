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

#include <MeshKernel/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGridRefinement.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Splines.hpp>

meshkernel::CurvilinearGridRefinement::CurvilinearGridRefinement(const std::shared_ptr<CurvilinearGrid>& grid, const Point& firstPoint, const Point& secondPoint, size_t refinement)
    : m_grid(grid),
      m_firstPoint(firstPoint),
      m_secondPoint(secondPoint),
      m_refinement(refinement)
{
}

void meshkernel::CurvilinearGridRefinement::Compute() const
{
    // Get the m and n indices from the point coordinates
    auto [mFirstNode, nFirstNode] = m_grid->GetNodeIndices(m_firstPoint);
    auto [mSecondNode, nSecondNode] = m_grid->GetNodeIndices(m_secondPoint);

    // The points must lie on the same gridline
    if (mSecondNode - mFirstNode != 0 && nSecondNode - nFirstNode != 0)
    {
        throw std::invalid_argument("CurvilinearGridRefinement::Compute: The selected curvilinear grid nodes are not on the same grid-line");
    }

    // Estimate the dimension of the refined grid
    const auto numMToRefine = mSecondNode - mFirstNode;
    const auto numNToRefine = nSecondNode - nFirstNode;
    size_t maxM = m_grid->m_numM + numMToRefine * (m_refinement - 1);
    size_t maxN = m_grid->m_numN + numNToRefine * (m_refinement - 1);

    // Store the m and n grid lines in separate vectors, and compute the spline derivatives for each gridline
    const auto [mGridLines, mGridLineDerivates, nGridLines, nGridLinesDerivatives] = m_grid->ComputeGridLinesAndSplinesDerivatives();

    // Local vector for each curvilinear grid face
    std::vector<Point> bottomRefinement(m_refinement);
    std::vector<Point> topRefinement(m_refinement);
    std::vector<Point> leftRefinement(m_refinement);
    std::vector<Point> rightRefinement(m_refinement);

    // The refined grid
    std::vector<std::vector<Point>> refinedGrid(maxM, std::vector<Point>(maxN));

    size_t refinedM = 0;
    for (auto currentM = 0; currentM < m_grid->m_numM - 1; ++currentM)
    {
        size_t localMRefinement = 1;
        if (currentM >= mFirstNode && currentM < mSecondNode)
        {
            localMRefinement = m_refinement;
        }

        size_t refinedN = 0;
        for (auto currentN = 0; currentN < m_grid->m_numN - 1; ++currentN)
        {

            size_t localNRefinement = 1;
            if (currentN >= nFirstNode && currentN < nSecondNode)
            {
                localNRefinement = m_refinement;
            }

            // Only if all grid nodes of the face are valid, perform transfinite interpolation
            if (m_grid->m_gridNodes[currentM][currentN].IsValid() &&
                m_grid->m_gridNodes[currentM + 1][currentN].IsValid() &&
                m_grid->m_gridNodes[currentM][currentN + 1].IsValid() &&
                m_grid->m_gridNodes[currentM + 1][currentN + 1].IsValid())
            {
                // Calculate m-direction spline points
                bottomRefinement.clear();
                topRefinement.clear();
                for (auto m = 0; m < localMRefinement + 1; ++m)
                {
                    const auto interpolationPoint = static_cast<double>(currentM) + static_cast<double>(m) / static_cast<double>(localMRefinement);
                    bottomRefinement.emplace_back(InterpolateSplinePoint(mGridLines[currentN], mGridLineDerivates[currentN], interpolationPoint));
                    topRefinement.emplace_back(InterpolateSplinePoint(mGridLines[currentN + 1], mGridLineDerivates[currentN + 1], interpolationPoint));
                }

                // Calculate n-direction spline points
                leftRefinement.clear();
                rightRefinement.clear();
                for (auto n = 0; n < localNRefinement + 1; ++n)
                {
                    const auto interpolationPoint = static_cast<double>(currentN) + static_cast<double>(n) / static_cast<double>(localNRefinement);
                    leftRefinement.emplace_back(InterpolateSplinePoint(nGridLines[currentM], nGridLinesDerivatives[currentM], interpolationPoint));
                    rightRefinement.emplace_back(InterpolateSplinePoint(nGridLines[currentM + 1], nGridLinesDerivatives[currentM + 1], interpolationPoint));
                }

                // Perform transfinite interpolation on the current curvilinear face
                const auto localGrid = DiscretizeTransfinite(leftRefinement,
                                                             rightRefinement,
                                                             bottomRefinement,
                                                             topRefinement,
                                                             m_grid->m_projection,
                                                             localMRefinement,
                                                             localNRefinement);
                // Copy the local grid into the refined grid
                for (auto m = 0; m < localMRefinement + 1; ++m)
                {
                    for (auto n = 0; n < localNRefinement + 1; ++n)
                    {
                        refinedGrid[refinedM + m][refinedN + n] = localGrid[m][n];
                    }
                }
            }
            refinedN += localNRefinement;
        }
        refinedM += localMRefinement;
    }

    // Substitute original grid with the refined one
    *m_grid = CurvilinearGrid(std::move(refinedGrid), m_grid->m_projection);
}
