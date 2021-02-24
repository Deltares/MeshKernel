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
#include <MeshKernel/Operations.hpp>
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

            size_t maxM = m_grid->m_numM - numMToRefine + numMToRefine * m_refinement;
            size_t maxN = m_grid->m_numN - numNToRefine + numNToRefine * n_refinement;

            const auto [mGridLines, mGridLineDerivates, nGridLines, nGridLinesDerivatives] = ComputeGridLinesAndSplinesDerivatives();

            std::vector<Point> bottomRefinement(m_refinement);
            std::vector<Point> upRefinementUp(m_refinement);
            std::vector<Point> leftRefinement(n_refinement);
            std::vector<Point> rightRefinement(n_refinement);
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
                        localNRefinement = n_refinement;
                    }

                    // calculate m-direction spline points
                    bottomRefinement.clear();
                    upRefinementUp.clear();
                    for (auto m = 0; m < localMRefinement; ++m)
                    {
                        const double interpolationPoint = currentM + m / localMRefinement;
                        bottomRefinement.emplace_back(InterpolateSplinePoint(mGridLines[currentM], mGridLineDerivates[currentM], interpolationPoint));
                        upRefinementUp.emplace_back(InterpolateSplinePoint(mGridLines[currentM + 1], mGridLineDerivates[currentM + 1], interpolationPoint));
                    }
                    // calculate n-direction spline points
                    leftRefinement.clear();
                    rightRefinement.clear();
                    for (auto n = 0; n < localNRefinement; ++n)
                    {
                        const double interpolationPoint = currentN + n / localNRefinement;
                        leftRefinement.emplace_back(InterpolateSplinePoint(nGridLines[currentM], nGridLinesDerivatives[currentM], interpolationPoint));
                        rightRefinement.emplace_back(InterpolateSplinePoint(nGridLines[currentM + 1], nGridLinesDerivatives[currentM + 1], interpolationPoint));
                    }

                    // perform transfinite interpolation
                    const auto localGrid = DiscretizeTransfinite(bottomRefinement,
                                                                 upRefinementUp,
                                                                 leftRefinement,
                                                                 rightRefinement,
                                                                 m_grid->m_projection,
                                                                 localMRefinement,
                                                                 localNRefinement);
                    // copy the local grid into the refined grid
                    for (auto m = 0; m < localMRefinement; ++m)
                    {
                        for (auto n = 0; n < localNRefinement; ++n)
                        {
                            refinedGrid[refinedM + m][refinedN + n] = localGrid[m][n];
                        }
                    }
                    refinedN += localNRefinement;
                }
                refinedM += localMRefinement;
            }
        }

    private:
        /// @brief
        /// @return
        [[nodiscard]] std::tuple<std::vector<std::vector<Point>>,
                                 std::vector<std::vector<Point>>,
                                 std::vector<std::vector<Point>>,
                                 std::vector<std::vector<Point>>>
        ComputeGridLinesAndSplinesDerivatives() const
        {

            // m-gridlines
            std::vector<std::vector<Point>> mGridLines(m_grid->m_numN, std::vector<Point>(m_grid->m_numM));
            std::vector<std::vector<Point>> mGridLineDerivates(m_grid->m_numN, std::vector<Point>(m_grid->m_numM));
            for (auto n = 0; n < m_grid->m_numN; ++n)
            {
                for (auto m = 0; m < m_grid->m_numM; ++m)
                {
                    mGridLines[n][m] = m_grid->m_nodes[m][n];
                }
                mGridLineDerivates[n] = ComputeSplineDerivatesAlongGridLine(mGridLines[n]);
            }

            // n-gridlines
            std::vector<std::vector<Point>> nGridLines(m_grid->m_numM, std::vector<Point>(m_grid->m_numN));
            std::vector<std::vector<Point>> nGridLinesDerivatives(m_grid->m_numM, std::vector<Point>(m_grid->m_numN));
            for (auto m = 0; m < m_grid->m_numM; ++m)
            {
                nGridLines[m] = m_grid->m_nodes[m];
                nGridLinesDerivatives[m] = ComputeSplineDerivatesAlongGridLine(nGridLines[m]);
            }

            return {mGridLines, mGridLineDerivates, nGridLines, nGridLinesDerivatives};
        }

        /// @brief Compute spline derivatives along a gridline, accounting for missing values
        /// @param gridLine The input gridline
        /// @return the spline derivatives
        std::vector<Point> ComputeSplineDerivatesAlongGridLine(const std::vector<Point>& gridLine) const
        {
            std::vector<Point> gridlineDerivatives(gridLine.size());
            const auto indices = FindIndices(gridLine, 0, gridLine.size(), doubleMissingValue);
            for (auto i = 0; i < indices.size(); ++i)
            {
                const auto derivatives = Splines::SecondOrderDerivative(gridLine, indices[i][0], indices[i][1]);
                for (size_t j = indices[i][0], ind = 0; j < indices[i][1]; ++j, ind++)
                {
                    gridlineDerivatives[j] = derivatives[ind];
                }
            }
            return gridlineDerivatives;
        }

        std::shared_ptr<CurvilinearGrid> m_grid; ///< A pointer to mesh
        Point m_firstPoint;
        Point m_secondPoint;
        size_t m_refinement;
        size_t n_refinement;
    }
} // namespace meshkernel