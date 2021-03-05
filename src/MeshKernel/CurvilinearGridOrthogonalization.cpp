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
    auto [mGridLines, mGridLineDerivates, nGridLines, nGridLinesDerivatives] = m_grid->ComputeGridLinesAndSplinesDerivatives();

    m_mGridLines = std::move(mGridLines);
    m_mGridLineDerivates = std::move(mGridLineDerivates);
    m_nGridLines = std::move(nGridLines);
    m_nGridLinesDerivatives = std::move(nGridLinesDerivatives);

    // Compute the matrix coefficients
    for (auto outerIterations = 0; outerIterations < m_orthogonalizationParameters.OuterIterations; ++outerIterations)
    {
        ComputeHorizontalMatrixCoefficients();
        ApplyBoundaryConditions();
        for (auto boundaryIterations = 0; boundaryIterations < m_orthogonalizationParameters.BoundaryIterations; ++boundaryIterations)
        {
            Solve();
            TreatBoundaryConditions();
        }
    }
}

void meshkernel::CurvilinearGridOrthogonalization::TreatBoundaryConditions()
{
    // horizontal lines
    for (auto n = 0; n < m_grid->m_numN; ++n)
    {
        int in = 0;
        int init = 0;
        int j = n;
        size_t leftIndex;
        size_t rightIndex;
        size_t nextVertical = 0;
        for (auto m = 0; m < m_grid->m_numM; ++m)
        {

            const auto nodeType = m_grid->m_gridNodesMask[m][n];
            if (nodeType == CurvilinearGrid::NodeTypes::BottomLeft || nodeType == CurvilinearGrid::NodeTypes::UpperLeft)
            {
                leftIndex = m;
                continue;
            }
            if (nodeType == CurvilinearGrid::NodeTypes::Bottom)
            {
                nextVertical = 1;
                continue;
            }
            if (nodeType == CurvilinearGrid::NodeTypes::Up)
            {
                nextVertical = -1;
                continue;
            }
            if ((nodeType == CurvilinearGrid::NodeTypes::BottomRight || nodeType == CurvilinearGrid::NodeTypes::UpperRight) && nextVertical != 0)
            {
                rightIndex = m;

                for (int mm = leftIndex + 1; mm < rightIndex - 1; ++mm)
                {

                    if (mm < m_minM || mm > m_maxM || n < m_minN || n > m_maxN)
                    {
                        continue;
                    }
                    if (m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeTypes::Invalid)
                    {
                        continue;
                    }
                    const auto leftNode = m_grid->m_gridNodes[m - 1][n];
                    const auto upperNode = m_grid->m_gridNodes[m][n + nextVertical];
                    const auto rightNode = m_grid->m_gridNodes[m + 1][n];
                    double qb;
                    double qc;

                    if (nextVertical == 1)
                    {
                        qb = m_atp[m - 1][n];
                        qc = m_atp[m][n];
                    }

                    if (nextVertical == -1)
                    {
                        qb = m_atp[m - 1][n - 1];
                        qc = m_atp[m][n - 1];
                    }
                    const auto qbc = 1.0 / qb + 1.0 / qc;
                    const auto rn = qb + qc + qbc;
                    Point searhNode = (leftNode * qb +
                                       upperNode * qbc +
                                       rightNode * qc + rightNode - leftNode) /
                                      rn;
                }

                //    m_mGridLines
                //m_mGridLineDerivates
                //m_nGridLines
                //m_nGridLinesDerivatives

                continue;
            }
        }
    }
}

void meshkernel::CurvilinearGridOrthogonalization::Solve()
{

    double relaxationFactor = 1.0;
    const double relaxationCoefficientOne = 1.0;
    const double relaxationCoefficientTwo = 0.5;
    const double relaxationCoefficientThree = 0.9 * 0.9;
    const double relaxationCoefficientFour = 0.25;
    for (auto innerIterations = 0; innerIterations < m_orthogonalizationParameters.InnerIterations; ++innerIterations)
    {
        for (auto m = m_minM; m < m_maxM; ++m)
        {
            for (auto n = m_minN; n < m_maxN; ++n)
            {
                if (m_grid->m_gridNodesMask[m][n] != CurvilinearGrid::NodeTypes::Valid)
                {
                    continue;
                }

                const auto residual =
                    m_grid->m_gridNodes[m][n] * m_a[m][n] +
                    m_grid->m_gridNodes[m][n] * m_b[m][n] +
                    m_grid->m_gridNodes[m][n] * m_c[m][n] +
                    m_grid->m_gridNodes[m][n] * m_d[m][n] +
                    m_grid->m_gridNodes[m][n] * m_e[m][n];
                m_grid->m_gridNodes[m][n] = m_grid->m_gridNodes[m][n] - residual / m_e[m][n] * relaxationFactor;
            }
        }

        if (innerIterations == 1)
        {
            relaxationFactor = relaxationCoefficientOne / (relaxationCoefficientOne - relaxationCoefficientTwo * relaxationCoefficientThree);
        }
        else
        {
            relaxationFactor = relaxationCoefficientOne / (relaxationCoefficientOne - relaxationFactor * relaxationCoefficientFour * relaxationCoefficientThree);
        }
    }
}

void meshkernel::CurvilinearGridOrthogonalization::ApplyBoundaryConditions()
{
    // To verify
}

void meshkernel::CurvilinearGridOrthogonalization::ComputeHorizontalMatrixCoefficients()
{
    /// allocate matrix coefficients
    m_a.resize(m_grid->m_numM, std::vector<double>(m_grid->m_numN, doubleMissingValue));
    m_b.resize(m_grid->m_numM, std::vector<double>(m_grid->m_numN, doubleMissingValue));
    m_c.resize(m_grid->m_numM, std::vector<double>(m_grid->m_numN, doubleMissingValue));
    m_d.resize(m_grid->m_numM, std::vector<double>(m_grid->m_numN, doubleMissingValue));
    m_e.resize(m_grid->m_numM, std::vector<double>(m_grid->m_numN, doubleMissingValue));
    m_atp.resize(m_grid->m_numM, std::vector<double>(m_grid->m_numN, doubleMissingValue));

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

    ComputeHorizontalCoefficients();
    ComputeVerticalCoefficients();

    // Normalize
    const auto smoothingFactor = 1.0 - m_orthogonalizationParameters.OrthogonalizationToSmoothingFactor;
    for (auto m = m_minM; m < m_maxM; ++m)
    {
        for (auto n = m_minN; n < m_maxN; ++n)
        {
            if (!m_grid->m_gridFacesMask[m][n])
            {
                continue;
            }
            m_atp[m][n] = m_atp[m][n] * m_a[m][n] / m_c[m][n];
            m_atp[m][n] = m_orthogonalizationParameters.OrthogonalizationToSmoothingFactor * m_atp[m][n] + smoothingFactor * m_a[m][n];
            m_e[m][n] = m_e[m][n] * m_b[m][n] / m_d[m][n];
            m_e[m][n] = m_orthogonalizationParameters.OrthogonalizationToSmoothingFactor * m_e[m][n] + smoothingFactor * m_b[m][n];

            m_a[m][n] = m_atp[m][n];
            m_b[m][n] = m_e[m][n];
        }
    }

    // calculate m_atp
    for (auto m = m_minM; m < m_maxM; ++m)
    {
        for (auto n = m_minN; n < m_maxN; ++n)
        {
            if (!m_grid->m_gridFacesMask[m][n])
            {
                m_atp[m][n] = doubleMissingValue;
                continue;
            }
            m_atp[m][n] = m_b[m][n] / m_a[m][n];
        }
    }

    // re-set coefficents

    for (auto m = m_minM + 1; m < m_maxM; ++m)
    {
        for (auto n = m_minN + 1; n < m_maxN; ++n)
        {
            if (!m_grid->m_gridFacesMask[m][n])
            {
                continue;
            }
            m_a[m][n] = m_atp[m][n - 1] + m_atp[m][n];
            m_b[m][n] = m_atp[m - 1][n - 1] + m_atp[m - 1][n];
            m_c[m][n] = 1.0 / m_atp[m - 1][n] + 1.0 / m_atp[m][n];
            m_d[m][n] = 1.0 / m_atp[m - 1][n - 1] + 1.0 / m_atp[m][n - 1];

            m_e[m][n] = -m_a[m][n] - m_b[m][n] - m_c[m][n] - m_d[m][n];
        }
    }
}

void meshkernel::CurvilinearGridOrthogonalization::ComputeVerticalCoefficients()
{
    std::vector<std::vector<bool>> invalidBoundaryNodes(m_grid->m_numM, std::vector<bool>(m_grid->m_numN, false));
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
            auto lastValidN = n + step;
            while (m_grid->m_gridNodesMask[m][lastValidN] == CurvilinearGrid::NodeTypes::Valid)
            {
                lastValidN += step;
            }
            const auto start = std::min(n, lastValidN);
            const auto end = std::max(n, lastValidN);
            for (auto k = start; k < end; k += step)
            {
                invalidBoundaryNodes[m][k] = true;
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
                !invalidBoundaryNodes[m][n])
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
                !invalidBoundaryNodes[m + 1][n])
            {
                m_a[m][n] = m_a[m + 1][n];
                m_c[m][n] = m_c[m + 1][n];
                counter[m][n] = counter[m + 1][n];
            }
        }
    }
    // average contributions
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

void meshkernel::CurvilinearGridOrthogonalization::ComputeHorizontalCoefficients()
{
    std::vector<std::vector<bool>> invalidBoundaryNodes(m_grid->m_numM, std::vector<bool>(m_grid->m_numN, false));
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
            auto lastValidM = m + step;
            while (m_grid->m_gridNodesMask[lastValidM][n] == CurvilinearGrid::NodeTypes::Valid)
            {
                lastValidM += step;
            }
            const auto start = std::min(m, lastValidM);
            const auto end = std::max(m, lastValidM);
            for (auto k = start; k < end; k += step)
            {
                invalidBoundaryNodes[k][n] = true;
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
                !invalidBoundaryNodes[m][n])
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
                !invalidBoundaryNodes[m][n + 1])
            {
                m_b[m][n] = m_b[m][n + 1];
                m_d[m][n] = m_d[m][n + 1];
                counter[m][n] = counter[m][n + 1];
            }
        }
    }

    // average contributions
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