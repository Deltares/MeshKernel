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
#include <MeshKernel/CurvilinearGridDeRefinement.hpp>
#include <MeshKernel/CurvilinearGridOrthogonalization.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Splines.hpp>

meshkernel::CurvilinearGridOrthogonalization::CurvilinearGridOrthogonalization(std::shared_ptr<CurvilinearGrid> grid,
                                                                               const meshkernelapi::OrthogonalizationParameters& orthogonalizationParameters,
                                                                               const Point& firstPoint,
                                                                               const Point& secondPoint)
    : m_grid(grid),
      m_orthogonalizationParameters(orthogonalizationParameters),
      m_firstPoint(firstPoint),
      m_secondPoint(secondPoint)

{
}

void meshkernel::CurvilinearGridOrthogonalization::Compute()
{

    // Get the m and n indices from the point coordinates
    auto [mFirstNode, nFirstNode] = m_grid->GetNodeIndices(m_firstPoint);
    auto [mSecondNode, nSecondNode] = m_grid->GetNodeIndices(m_secondPoint);

    /// The points must lie on the same GridLine
    if (mFirstNode == mSecondNode && nFirstNode == nSecondNode)
    {
        throw std::invalid_argument("CurvilinearGridDeRefinement::Compute: bottom left and upper right selection corners coincides");
    }

    /// define the orthogonalization bounding box
    m_minM = std::max(0, mFirstNode - 1);
    m_minN = std::max(0, nFirstNode - 1);
    m_maxM = std::min(static_cast<int>(m_grid->m_numM - 1), mSecondNode);
    m_maxN = std::min(static_cast<int>(m_grid->m_numN - 1), nSecondNode);

    /// Compute the masks
    m_grid->ComputeGridMasks();
    /// store the m and n grid lines in separate vectors, and compute the spline derivatives for each gridline
    const auto [mGridLines, mGridLineDerivates, nGridLines, nGridLinesDerivatives] = m_grid->ComputeGridLinesAndSplinesDerivatives();

    /// allocate matrix coefficients
    m_a.resize(m_grid->m_numM, std::vector<double>(m_grid->m_numN, doubleMissingValue));
    m_b.resize(m_grid->m_numM, std::vector<double>(m_grid->m_numN, doubleMissingValue));
    m_c.resize(m_grid->m_numM, std::vector<double>(m_grid->m_numN, doubleMissingValue));
    m_d.resize(m_grid->m_numM, std::vector<double>(m_grid->m_numN, doubleMissingValue));
    m_e.resize(m_grid->m_numM, std::vector<double>(m_grid->m_numN, doubleMissingValue));
    m_atp.resize(m_grid->m_numM, std::vector<double>(m_grid->m_numN, doubleMissingValue));

    for (auto outerIterations = 0; outerIterations < m_orthogonalizationParameters.OuterIterations; ++outerIterations)
    {
        for (auto m = m_minM; m < m_maxM; ++m)
        {
            for (auto n = m_minN; n < m_maxN; ++n)
            {
                if (!m_grid->m_gridFacesMask[m][n])
                {
                    continue;
                }

                const auto bottom = ComputeDistance(m_grid->m_gridNodes[m][n], m_grid->m_gridNodes[m + 1][n], Projection::cartesian);
                const auto upper = ComputeDistance(m_grid->m_gridNodes[m][n + 1], m_grid->m_gridNodes[m + 1][n + 1], Projection::cartesian);
                const auto left = ComputeDistance(m_grid->m_gridNodes[m][n], m_grid->m_gridNodes[m][n + 1], Projection::cartesian);
                const auto right = ComputeDistance(m_grid->m_gridNodes[m + 1][n], m_grid->m_gridNodes[m + 1][n + 1], Projection::cartesian);

                m_a[m][n] = (bottom + upper) * 0.5;
                m_b[m][n] = (left + right) * 0.5;

                m_atp[m][n] = m_a[m][n];
                m_e[m][n] = m_b[m][n];

                m_c[m][n] = m_atp[m][n];
                m_d[m][n] = m_e[m][n];
            }
        }
    }

    ComputeHorizontalMatrixCoefficients();
    ComputeHorizontalVerticalCoefficients();
}

void meshkernel::CurvilinearGridOrthogonalization::ComputeHorizontalVerticalCoefficients()
{
    for (auto m = m_minM + 1; m < m_maxM; ++m)
    {
        for (auto n = m_minN + 1; n < m_maxN; ++n)
        {
            int step = 0;
            if (m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeTypes::BottomLeft)
            {
                step = -1;
            }
            if (m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeTypes::BottomRight)
            {
                step = -1;
            }
            if (m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeTypes::UpperRight)
            {
                step = 1;
            }
            if (m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeTypes::UpperLeft)
            {
                step = 1;
            }
            if (step == 0)
            {
                continue;
            }
            auto lastValidN = n;
            while (m_grid->m_gridNodesMask[m][lastValidN] == CurvilinearGrid::NodeTypes::Valid)
            {
                lastValidN += step;
            }
            const auto start = std::min(n, lastValidN);
            const auto end = std::max(n, lastValidN);
            for (auto k = start; k < end; k += step)
            {
                m_grid->m_gridNodesMask[m][k] = CurvilinearGrid::NodeTypes::OrthogonalizationNodes;
            }
        }
    }

    // store the counter
    std::vector<std::vector<size_t>> counter(m_grid->m_numM, std::vector<size_t>(m_grid->m_numN, 0));

    // Perform left sum
    for (auto n = m_minN; n < m_maxN; ++n)
    {
        for (auto m = m_minM + 1; m < m_maxM; ++m)
        {

            if (m_grid->IsValidFace(m, n) &&
                !IsEqual(m_a[m][n], doubleMissingValue) &&
                !IsEqual(m_a[m - 1][n], doubleMissingValue) &&
                m_grid->m_gridNodesMask[m][n] != CurvilinearGrid::NodeTypes::OrthogonalizationNodes)
            {
                m_a[m][n] += m_a[m - 1][n];
                m_c[m][n] += m_c[m - 1][n];
                counter[m][n] = counter[m - 1][n] + 1;
            }
        }
    }

    // Perform right sum
    for (auto n = m_minN; n < m_maxN; ++n)
    {
        for (auto m = m_minM - 1; m >= m_maxM; --m)
        {
            if (m_grid->IsValidFace(m, n) &&
                !IsEqual(m_a[m][n], doubleMissingValue) &&
                !IsEqual(m_a[m + 1][n], doubleMissingValue) &&
                m_grid->m_gridNodesMask[m + 1][n] != CurvilinearGrid::NodeTypes::OrthogonalizationNodes)
            {
                m_a[m][n] = m_a[m + 1][n];
                m_c[m][n] = m_c[m + 1][n];
                counter[m][n] = counter[m + 1][n];
            }
        }
    }
    // make the average
    for (auto m = m_minM; m < m_maxM; ++m)
    {
        for (auto n = m_minN; n < m_maxN; ++n)
        {
            if (m_grid->IsValidFace(m, n))
            {
                m_a[m][n] /= static_cast<double>(counter[m][n] + 1);
                m_c[m][n] /= static_cast<double>(counter[m][n] + 1);
            }
        }
    }
}

void meshkernel::CurvilinearGridOrthogonalization::ComputeHorizontalMatrixCoefficients()
{

    for (auto m = m_minM + 1; m < m_maxM; ++m)
    {
        for (auto n = m_minN + 1; n < m_maxN; ++n)
        {
            int step = 0;
            if (m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeTypes::BottomLeft)
            {
                step = -1;
            }
            if (m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeTypes::BottomRight)
            {
                step = 1;
            }
            if (m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeTypes::UpperRight)
            {
                step = 1;
            }
            if (m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeTypes::UpperLeft)
            {
                step = -1;
            }
            if (step == 0)
            {
                continue;
            }
            auto lastValidM = m;
            while (m_grid->m_gridNodesMask[lastValidM][n] == CurvilinearGrid::NodeTypes::Valid)
            {
                lastValidM += step;
            }
            const auto start = std::min(m, lastValidM);
            const auto end = std::max(m, lastValidM);
            for (auto k = start; k < end; k += step)
            {
                m_grid->m_gridNodesMask[k][n] = CurvilinearGrid::NodeTypes::OrthogonalizationNodes;
            }
        }
    }

    std::vector<std::vector<size_t>> counter(m_grid->m_numM, std::vector<size_t>(m_grid->m_numN, 0));

    // Perform bottom sum
    for (auto m = m_minM; m < m_maxM; ++m)
    {
        for (auto n = m_minN + 1; n < m_maxN; ++n)
        {
            if (m_grid->IsValidFace(m, n) &&
                !IsEqual(m_b[m][n], doubleMissingValue) &&
                !IsEqual(m_b[m][n - 1], doubleMissingValue) &&
                m_grid->m_gridNodesMask[m][n] != CurvilinearGrid::NodeTypes::OrthogonalizationNodes)
            {
                m_b[m][n] += m_b[m][n - 1];
                m_d[m][n] += m_d[m][n - 1];
                counter[m][n] = counter[m][n - 1] + 1;
            }
        }
    }

    // Perform upper sum
    for (auto m = m_minM; m < m_maxM; ++m)
    {
        for (auto n = m_maxN - 1; n >= m_minN; --n)
        {
            if (m_grid->IsValidFace(m, n) &&
                !IsEqual(m_b[m][n], doubleMissingValue) &&
                !IsEqual(m_b[m][n + 1], doubleMissingValue) &&
                m_grid->m_gridNodesMask[m][n + 1] != CurvilinearGrid::NodeTypes::OrthogonalizationNodes)
            {
                m_b[m][n] = m_b[m][n + 1];
                m_d[m][n] = m_d[m][n + 1];
                counter[m][n] = counter[m][n + 1];
            }
        }
    }
    // make the average
    for (auto m = m_minM; m < m_maxM; ++m)
    {
        for (auto n = m_minN; n < m_maxN; ++n)
        {
            if (m_grid->IsValidFace(m, n))
            {
                m_b[m][n] /= static_cast<double>(counter[m][n] + 1);
                m_d[m][n] /= static_cast<double>(counter[m][n] + 1);
            }
        }
    }
}