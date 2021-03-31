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

#include <MeshKernel/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGridDeRefinement.hpp>
#include <MeshKernel/CurvilinearGridOrthogonalization.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Splines.hpp>

using meshkernel::CurvilinearGrid;
using meshkernel::CurvilinearGridOrthogonalization;

CurvilinearGridOrthogonalization::CurvilinearGridOrthogonalization(std::shared_ptr<CurvilinearGrid> grid,
                                                                   const meshkernelapi::OrthogonalizationParameters& orthogonalizationParameters)
    : m_grid(grid),
      m_orthogonalizationParameters(orthogonalizationParameters)

{
    /// Store the grid lines of the curvilinear grid as splines
    m_splines = Splines(m_grid);

    /// allocate matrix coefficients
    m_a.resize(m_grid->m_numM, std::vector<double>(m_grid->m_numN, doubleMissingValue));
    m_b.resize(m_grid->m_numM, std::vector<double>(m_grid->m_numN, doubleMissingValue));
    m_c.resize(m_grid->m_numM, std::vector<double>(m_grid->m_numN, doubleMissingValue));
    m_d.resize(m_grid->m_numM, std::vector<double>(m_grid->m_numN, doubleMissingValue));
    m_e.resize(m_grid->m_numM, std::vector<double>(m_grid->m_numN, doubleMissingValue));
    m_atp.resize(m_grid->m_numM, std::vector<double>(m_grid->m_numN, doubleMissingValue));
    m_isGridNodeFrozen.resize(m_grid->m_numM, std::vector<bool>(m_grid->m_numN, false));
}

void meshkernel::CurvilinearGridOrthogonalization::SetBlock(Point const& firstCornerPoint, Point const& secondCornerPoint)
{
    // Get the m and n indices from the point coordinates
    auto [mFirstNode, nFirstNode] = m_grid->GetNodeIndices(firstCornerPoint);
    auto [mSecondNode, nSecondNode] = m_grid->GetNodeIndices(secondCornerPoint);

    // Coinciding corner nodes, no valid area, nothing to do
    if (mFirstNode == mSecondNode && nFirstNode == nSecondNode)
    {
        return;
    }

    // Compute orthogonalization bounding box
    const auto [minM, minN, maxM, maxN] = OrderCoordinates(mFirstNode, nFirstNode, mSecondNode, nSecondNode);
    m_minM = minM;
    m_minN = minN;
    m_maxM = maxM;
    m_maxN = maxN;
}

std::tuple<size_t, size_t, size_t, size_t> meshkernel::CurvilinearGridOrthogonalization::OrderCoordinates(size_t firstPointM, size_t firstPointN, size_t secondPointM, size_t secondPointN) const
{
    return {std::min(firstPointM, secondPointM), std::min(firstPointN, secondPointN), std::max(firstPointM, secondPointM), std::max(firstPointN, secondPointN)};
}

void meshkernel::CurvilinearGridOrthogonalization::SetFrozenLine(Point const& firstPoint, Point const& secondPoint)
{
    // Get the m and n indices from the point coordinates
    auto [mFirstNode, nFirstNode] = m_grid->GetNodeIndices(firstPoint);
    auto [mSecondNode, nSecondNode] = m_grid->GetNodeIndices(secondPoint);

    // Coinciding nodes, no valid line, nothing to do
    if (mFirstNode == mSecondNode && nFirstNode == nSecondNode)
    {
        return;
    }

    // The selected nodes must be on the vertical or horizontal line
    const auto [minM, minN, maxM, maxN] = OrderCoordinates(mFirstNode, nFirstNode, mSecondNode, nSecondNode);
    auto const deltaMNewLine = maxM - minM;
    if (auto const deltaNNewLine = maxN - minN; deltaMNewLine != 0 && deltaNNewLine != 0)
    {
        throw std::invalid_argument("CurvilinearGridOrthogonalization::SetFrozenLine the points of the line to freeze must lie on the same grid line");
    }

    // The start-end coordinates of the new line, along m or n direction
    auto const startNewLine = deltaMNewLine == 0 ? minN : minM;
    auto const endNewLine = deltaMNewLine == 0 ? maxN : maxM;
    auto const constantCoordinateLine = deltaMNewLine == 0 ? minM : minN;

    // Frozen lines cannot cross existing frozen lines
    for (auto const& frozenLine : m_frozenLines)
    {
        auto const& [minMCurrentLine, minNCurrentLine, maxMCurrentLine, maxNCurrentLine] = frozenLine;

        auto const deltaMCurrentLine = maxMCurrentLine - minMCurrentLine;
        auto const startCurrentLine = deltaMCurrentLine == 0 ? minNCurrentLine : minMCurrentLine;
        auto const endCurrentLine = deltaMCurrentLine == 0 ? maxNCurrentLine : maxMCurrentLine;
        auto const constantCoordinateCurrentLine = deltaMCurrentLine == 0 ? minMCurrentLine : minNCurrentLine;
        for (auto i = startCurrentLine; i <= endCurrentLine; ++i)
        {

            for (auto j = startNewLine; j <= endNewLine; ++j)
            {
                if (j == constantCoordinateCurrentLine && i == constantCoordinateLine)
                {
                    throw AlgorithmError("CurvilinearGridOrthogonalization::SetFrozenLine the new line to freeze is crossing an existing line");
                }
            }
        }
    }

    // Now a frozen line can be stored
    m_frozenLines.emplace_back(minM, minN, maxM, maxN);
}

void meshkernel::CurvilinearGridOrthogonalization::ComputeFrozenGridPoints()
{
    for (auto const& frozenLine : m_frozenLines)
    {
        auto const& [minMCurrentLine, minNCurrentLine, maxMCurrentLine, maxNCurrentLine] = frozenLine;

        for (auto m = minMCurrentLine; m <= maxMCurrentLine; ++m)
        {
            for (auto n = minNCurrentLine; n <= maxNCurrentLine; ++n)
            {
                m_isGridNodeFrozen[m][n] = true;
            }
        }
    }
}

void meshkernel::CurvilinearGridOrthogonalization::Compute()
{
    // Compute the grid node types
    m_grid->ComputeGridNodeTypes();

    // Set the frozen node mask
    ComputeFrozenGridPoints();

    // Compute the matrix coefficients
    for (auto outerIterations = 0; outerIterations < m_orthogonalizationParameters.OuterIterations; ++outerIterations)
    {
        ComputeCoefficients();
        for (auto boundaryIterations = 0; boundaryIterations < m_orthogonalizationParameters.BoundaryIterations; ++boundaryIterations)
        {
            Solve();
            ProjectHorizontalBoundaryGridNodes();
            ProjectVerticalBoundariesGridNodes();
        }
    }
}

void CurvilinearGridOrthogonalization::ProjectHorizontalBoundaryGridNodes()
{
    // m grid lines (horizontal)
    for (auto n = 0; n < m_grid->m_numN; ++n)
    {
        size_t startM = sizetMissingValue;
        int nextVertical = 0;
        for (auto m = 0; m < m_grid->m_numM; ++m)
        {
            const auto nodeType = m_grid->m_gridNodesMask[m][n];
            if (nodeType == CurvilinearGrid::NodeType::BottomLeft || nodeType == CurvilinearGrid::NodeType::UpperLeft)
            {
                startM = m;
                continue;
            }
            if (nodeType == CurvilinearGrid::NodeType::Bottom)
            {
                nextVertical = 1;
                continue;
            }
            if (nodeType == CurvilinearGrid::NodeType::Up)
            {
                nextVertical = -1;
                continue;
            }

            // Project the nodes at the boundary (Bottom and Up node types) if a valid interval has been found.
            // The interval ranges from startM to the next BottomRight or UpperRight node.
            if (startM != sizetMissingValue &&
                (nodeType == CurvilinearGrid::NodeType::BottomRight || nodeType == CurvilinearGrid::NodeType::UpperRight) &&
                nextVertical != 0)
            {
                for (auto mm = startM + 1; mm < m; ++mm)
                {

                    if (mm < m_minM || mm > m_maxM || n < m_minN || n > m_maxN)
                    {
                        continue;
                    }
                    if (m_grid->m_gridNodesMask[mm][n] == CurvilinearGrid::NodeType::Invalid)
                    {
                        continue;
                    }

                    const auto leftNode = m_grid->m_gridNodes[mm - 1][n];
                    const auto verticalNode = m_grid->m_gridNodes[mm][n + nextVertical];
                    const auto rightNode = m_grid->m_gridNodes[mm + 1][n];

                    Point boundaryNode;
                    if (nextVertical == 1)
                    {
                        const double qb = m_atp[mm - 1][n];
                        const double qc = m_atp[mm][n];
                        const auto qbc = 1.0 / qb + 1.0 / qc;
                        const auto rn = qb + qc + qbc;
                        boundaryNode.x = (leftNode.x * qb + verticalNode.x * qbc + rightNode.x * qc + rightNode.y - leftNode.y) / rn;
                        boundaryNode.y = (leftNode.y * qb + verticalNode.y * qbc + rightNode.y * qc + leftNode.x - rightNode.x) / rn;
                    }

                    if (nextVertical == -1)
                    {
                        const double qb = m_atp[mm - 1][n - 1];
                        const double qc = m_atp[mm][n - 1];
                        const auto qbc = 1.0 / qb + 1.0 / qc;
                        const auto rn = qb + qc + qbc;
                        boundaryNode.x = (leftNode.x * qb + verticalNode.x * qbc + rightNode.x * qc + leftNode.y - rightNode.y) / rn;
                        boundaryNode.y = (leftNode.y * qb + verticalNode.y * qbc + rightNode.y * qc + rightNode.x - leftNode.x) / rn;
                    }

                    m_grid->m_gridNodes[mm][n] = m_splines.ComputeClosestPointOnSplineSegment(n,
                                                                                              static_cast<double>(startM),
                                                                                              static_cast<double>(m),
                                                                                              boundaryNode);
                }
            }
        }
    }
}

void CurvilinearGridOrthogonalization::ProjectVerticalBoundariesGridNodes()
{
    // n gridlines (vertical)
    for (auto m = 0; m < m_grid->m_numM; ++m)
    {
        size_t startN = sizetMissingValue;
        int nextHorizontal = 0;
        for (auto n = 0; n < m_grid->m_numN; ++n)
        {
            const auto nodeType = m_grid->m_gridNodesMask[m][n];
            if (nodeType == CurvilinearGrid::NodeType::BottomLeft || nodeType == CurvilinearGrid::NodeType::BottomRight)
            {
                startN = n;
                continue;
            }
            if (nodeType == CurvilinearGrid::NodeType::Left)
            {
                nextHorizontal = 1;
                continue;
            }
            if (nodeType == CurvilinearGrid::NodeType::Right)
            {
                nextHorizontal = -1;
                continue;
            }

            // Project the nodes at the boundary (Left and Left node types) if a valid interval has been found.
            // The interval ranges from startN to the next UpperLeft or UpperRight node.
            if ((nodeType == CurvilinearGrid::NodeType::UpperLeft || nodeType == CurvilinearGrid::NodeType::UpperRight) &&
                nextHorizontal != 0 &&
                startN != sizetMissingValue)
            {
                for (auto nn = startN + 1; nn < n; ++nn)
                {

                    if (m < m_minM || m > m_maxM || nn < m_minN || nn > m_maxN)
                    {
                        continue;
                    }
                    if (m_grid->m_gridNodesMask[m][nn] == CurvilinearGrid::NodeType::Invalid)
                    {
                        continue;
                    }
                    const auto bottomNode = m_grid->m_gridNodes[m][nn - 1];
                    const auto horizontalNode = m_grid->m_gridNodes[m + nextHorizontal][nn];
                    const auto upperNode = m_grid->m_gridNodes[m][nn + 1];

                    Point boundaryNode;
                    if (nextHorizontal == 1)
                    {
                        const auto qb = 1.0 / m_atp[m][nn - 1];
                        const auto qc = 1.0 / m_atp[m][nn];
                        const auto qbc = m_atp[m][nn - 1] + m_atp[m][nn];
                        const auto rn = qb + qc + qbc;
                        boundaryNode.x = (bottomNode.x * qb + horizontalNode.x * qbc + upperNode.x * qc + bottomNode.y - upperNode.y) / rn;
                        boundaryNode.y = (bottomNode.y * qb + horizontalNode.y * qbc + upperNode.y * qc + upperNode.x - bottomNode.x) / rn;
                    }

                    if (nextHorizontal == -1)
                    {
                        const auto qb = 1.0 / m_atp[m - 1][nn - 1];
                        const auto qc = 1.0 / m_atp[m - 1][nn];
                        const auto qbc = m_atp[m - 1][nn - 1] + m_atp[m - 1][nn];
                        const auto rn = qb + qc + qbc;
                        boundaryNode.x = (bottomNode.x * qb + horizontalNode.x * qbc + upperNode.x * qc + upperNode.y - bottomNode.y) / rn;
                        boundaryNode.y = (bottomNode.y * qb + horizontalNode.y * qbc + upperNode.y * qc + bottomNode.x - upperNode.x) / rn;
                    }

                    // Vertical spline index
                    const auto splineIndex = m_grid->m_numN + m;
                    m_grid->m_gridNodes[m][nn] = m_splines.ComputeClosestPointOnSplineSegment(splineIndex,
                                                                                              static_cast<double>(startN),
                                                                                              static_cast<double>(n),
                                                                                              boundaryNode);
                }
            }
        }
    }
}

void CurvilinearGridOrthogonalization::Solve()
{

    double omega = 1.0;
    const double factor = 0.9 * 0.9;

    // Only the internal nodes of the orthogonalization box
    const auto minMInternal = std::max(static_cast<size_t>(1), m_minM);
    const auto minNInternal = std::max(static_cast<size_t>(1), m_minN);

    const auto maxMInternal = std::min(m_maxM, m_grid->m_numM - 1);
    const auto maxNInternal = std::min(m_maxN, m_grid->m_numN - 1);

    for (auto innerIterations = 0; innerIterations < m_orthogonalizationParameters.InnerIterations; ++innerIterations)
    {
        for (auto m = minMInternal; m < maxMInternal; ++m)
        {
            for (auto n = minNInternal; n < maxNInternal; ++n)
            {
                if (m_grid->m_gridNodesMask[m][n] != CurvilinearGrid::NodeType::InternalValid)
                {
                    continue;
                }

                if (m_isGridNodeFrozen[m][n])
                {
                    continue;
                }

                const auto residual =
                    m_grid->m_gridNodes[m + 1][n] * m_a[m][n] +
                    m_grid->m_gridNodes[m - 1][n] * m_b[m][n] +
                    m_grid->m_gridNodes[m][n + 1] * m_c[m][n] +
                    m_grid->m_gridNodes[m][n - 1] * m_d[m][n] +
                    m_grid->m_gridNodes[m][n] * m_e[m][n];

                m_grid->m_gridNodes[m][n] = m_grid->m_gridNodes[m][n] - residual / m_e[m][n] * omega;
            }
        }

        if (innerIterations == 0)
        {
            omega = 1.0 / (1.0 - 0.5 * factor);
        }
        else
        {
            omega = 1.0 / (1.0 - omega * 0.25 * factor);
        }
    }
}

void meshkernel::CurvilinearGridOrthogonalization::ComputeCoefficients()
{
    /// allocate matrix coefficients
    std::fill(m_a.begin(), m_a.end(), std::vector<double>(m_grid->m_numN, doubleMissingValue));
    std::fill(m_b.begin(), m_b.end(), std::vector<double>(m_grid->m_numN, doubleMissingValue));
    std::fill(m_c.begin(), m_c.end(), std::vector<double>(m_grid->m_numN, doubleMissingValue));
    std::fill(m_d.begin(), m_d.end(), std::vector<double>(m_grid->m_numN, doubleMissingValue));
    std::fill(m_e.begin(), m_e.end(), std::vector<double>(m_grid->m_numN, doubleMissingValue));
    std::fill(m_atp.begin(), m_atp.end(), std::vector<double>(m_grid->m_numN, doubleMissingValue));

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

    // Calculate m_atp
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

    // Re-set coefficients
    std::fill(m_a.begin(), m_a.end(), std::vector<double>(m_grid->m_numN, 0.0));
    std::fill(m_b.begin(), m_b.end(), std::vector<double>(m_grid->m_numN, 0.0));
    std::fill(m_c.begin(), m_c.end(), std::vector<double>(m_grid->m_numN, 0.0));
    std::fill(m_d.begin(), m_d.end(), std::vector<double>(m_grid->m_numN, 0.0));
    std::fill(m_e.begin(), m_e.end(), std::vector<double>(m_grid->m_numN, 0.0));

    for (auto m = m_minM + 1; m < m_maxM; ++m)
    {
        for (auto n = m_minN + 1; n < m_maxN; ++n)
        {
            if (m_grid->m_gridNodesMask[m][n] != CurvilinearGrid::NodeType::InternalValid)
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

void CurvilinearGridOrthogonalization::ComputeVerticalCoefficients()
{
    const auto invalidBoundaryNodes = ComputeInvalidVerticalBoundaryNodes();
    // Store the counter
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
        for (auto m = int(m_maxM) - 1; m >= int(m_minM); --m)
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
    // Average contributions
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

void CurvilinearGridOrthogonalization::ComputeHorizontalCoefficients()
{
    const auto invalidBoundaryNodes = ComputeInvalidHorizontalBoundaryNodes();
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
        for (auto n = int(m_maxN) - 1; n >= int(m_minN); --n)
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

    // Average contributions
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

std::vector<std::vector<bool>>
CurvilinearGridOrthogonalization::ComputeInvalidHorizontalBoundaryNodes() const
{
    std::vector<std::vector<bool>> invalidBoundaryNodes(m_grid->m_numM, std::vector<bool>(m_grid->m_numN, false));
    for (auto m = m_minM + 1; m < m_maxM; ++m)
    {
        for (auto n = m_minN + 1; n < m_maxN; ++n)
        {
            int step = 0;
            if (m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeType::BottomLeft)
            {
                step = -1;
            }
            if (m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeType::BottomRight)
            {
                step = 1;
            }
            if (m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeType::UpperRight)
            {
                step = 1;
            }
            if (m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeType::UpperLeft)
            {
                step = -1;
            }
            if (step == 0)
            {
                continue;
            }

            auto lastValidM = m + step;
            while (lastValidM > 0 &&
                   lastValidM < m_grid->m_numM &&
                   m_grid->m_gridNodesMask[lastValidM][n] == CurvilinearGrid::NodeType::InternalValid)
            {
                lastValidM += step;
            }
            const auto start = std::min(m, lastValidM);
            const auto end = std::max(m, lastValidM);
            for (auto k = start; k <= end; ++k)
            {
                invalidBoundaryNodes[k][n] = true;
            }
        }
    }

    return invalidBoundaryNodes;
}

std::vector<std::vector<bool>>
CurvilinearGridOrthogonalization::ComputeInvalidVerticalBoundaryNodes() const
{
    std::vector<std::vector<bool>> invalidBoundaryNodes(m_grid->m_numM, std::vector<bool>(m_grid->m_numN, false));
    for (auto m = m_minM + 1; m < m_maxM; ++m)
    {
        for (auto n = m_minN + 1; n < m_maxN; ++n)
        {
            int step = 0;
            if (m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeType::BottomLeft)
            {
                step = -1;
            }
            if (m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeType::BottomRight)
            {
                step = -1;
            }
            if (m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeType::UpperRight)
            {
                step = 1;
            }
            if (m_grid->m_gridNodesMask[m][n] == CurvilinearGrid::NodeType::UpperLeft)
            {
                step = 1;
            }
            if (step == 0)
            {
                continue;
            }
            auto lastValidN = n + step;
            while (lastValidN > 0 &&
                   lastValidN < m_grid->m_numN &&
                   m_grid->m_gridNodesMask[m][lastValidN] == CurvilinearGrid::NodeType::InternalValid)
            {
                lastValidN += step;
            }
            const auto start = std::min(n, lastValidN);
            const auto end = std::max(n, lastValidN);
            for (auto k = start; k <= end; ++k)
            {
                invalidBoundaryNodes[m][k] = true;
            }
        }
    }

    return invalidBoundaryNodes;
}
