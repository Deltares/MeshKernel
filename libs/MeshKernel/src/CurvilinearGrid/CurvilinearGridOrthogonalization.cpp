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
#include <MeshKernel/CurvilinearGrid/CurvilinearGridDeRefinement.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridOrthogonalization.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Splines.hpp>

using meshkernel::CurvilinearGrid;
using meshkernel::CurvilinearGridOrthogonalization;

CurvilinearGridOrthogonalization::CurvilinearGridOrthogonalization(std::shared_ptr<CurvilinearGrid> grid,
                                                                   const OrthogonalizationParameters& orthogonalizationParameters)
    : CurvilinearGridAlgorithm(grid),
      m_orthogonalizationParameters(orthogonalizationParameters),
      m_terms(m_grid.m_numM, m_grid.m_numN)

{
    lin_alg::ResizeAndFillMatrix(m_isGridNodeFrozen, m_grid.m_numM, m_grid.m_numN, true, false);

    /// Store the grid lines of the grid as splines
    m_splines = Splines(m_grid);
}

void CurvilinearGridOrthogonalization::ComputeFrozenGridPoints()
{
    for (auto const& frozenLine : m_lines)
    {
        for (auto m = frozenLine.m_startNode.m_m; m <= frozenLine.m_endNode.m_m; ++m)
        {
            for (auto n = frozenLine.m_startNode.m_n; n <= frozenLine.m_endNode.m_n; ++n)
            {
                m_isGridNodeFrozen(m, n) = true;
            }
        }
    }
}

CurvilinearGrid CurvilinearGridOrthogonalization::Compute()
{
    if (!m_lowerLeft.IsValid() || !m_upperRight.IsValid())
    {
        throw std::invalid_argument("CurvilinearGridOrthogonalization::Compute: lower left and upper right corners defining the curvilinear grid block are not set");
    }

    // Compute the grid node types
    m_grid.ComputeGridNodeTypes();

    // Set the frozen node mask
    ComputeFrozenGridPoints();

    // Compute the matrix coefficients
    for (auto outerIterations = 0; outerIterations < m_orthogonalizationParameters.outer_iterations; ++outerIterations)
    {
        ComputeCoefficients();
        for (auto boundaryIterations = 0; boundaryIterations < m_orthogonalizationParameters.boundary_iterations; ++boundaryIterations)
        {
            Solve();
            ProjectHorizontalBoundaryGridNodes();
            ProjectVerticalBoundariesGridNodes();
        }
    }

    return m_grid;
}

void CurvilinearGridOrthogonalization::ProjectHorizontalBoundaryGridNodes()
{
    // m grid lines (horizontal)
    for (UInt n = 0; n < m_grid.m_numN; ++n)
    {
        UInt startM = constants::missing::uintValue;
        int nextVertical = 0;
        for (UInt m = 0; m < m_grid.m_numM; ++m)
        {
            const auto nodeType = m_grid.m_gridNodesTypes(m, n);
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
            if (startM != constants::missing::uintValue &&
                (nodeType == CurvilinearGrid::NodeType::BottomRight || nodeType == CurvilinearGrid::NodeType::UpperRight) &&
                nextVertical != 0)
            {
                for (auto mm = startM + 1; mm < m; ++mm)
                {

                    if (mm < m_lowerLeft.m_m || mm > m_upperRight.m_m || n < m_lowerLeft.m_n || n > m_upperRight.m_n)
                    {
                        continue;
                    }
                    if (m_grid.m_gridNodesTypes(mm, n) == CurvilinearGrid::NodeType::Invalid)
                    {
                        continue;
                    }

                    const auto leftNode = m_grid.m_gridNodes(mm - 1, n);
                    const auto verticalNode = m_grid.m_gridNodes(mm, n + nextVertical);
                    const auto rightNode = m_grid.m_gridNodes(mm + 1, n);

                    Point boundaryNode;
                    if (nextVertical == 1)
                    {
                        const double qb = m_terms.atp(mm - 1, n);
                        const double qc = m_terms.atp(mm, n);
                        const auto qbc = 1.0 / qb + 1.0 / qc;
                        const auto rn = qb + qc + qbc;
                        boundaryNode.x = (leftNode.x * qb + verticalNode.x * qbc + rightNode.x * qc + rightNode.y - leftNode.y) / rn;
                        boundaryNode.y = (leftNode.y * qb + verticalNode.y * qbc + rightNode.y * qc + leftNode.x - rightNode.x) / rn;
                    }

                    if (nextVertical == -1)
                    {
                        const double qb = m_terms.atp(mm - 1, n - 1);
                        const double qc = m_terms.atp(mm, n - 1);
                        const auto qbc = 1.0 / qb + 1.0 / qc;
                        const auto rn = qb + qc + qbc;
                        boundaryNode.x = (leftNode.x * qb + verticalNode.x * qbc + rightNode.x * qc + leftNode.y - rightNode.y) / rn;
                        boundaryNode.y = (leftNode.y * qb + verticalNode.y * qbc + rightNode.y * qc + rightNode.x - leftNode.x) / rn;
                    }

                    m_grid.m_gridNodes(mm, n) = m_splines.ComputeClosestPointOnSplineSegment(n,
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
    // m gridlines (vertical)
    for (UInt m = 0; m < m_grid.m_numM; ++m)
    {
        UInt startN = constants::missing::uintValue;
        int nextHorizontal = 0;
        for (UInt n = 0; n < m_grid.m_numN; ++n)
        {
            const auto nodeType = m_grid.m_gridNodesTypes(m, n);
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
                startN != constants::missing::uintValue)
            {
                for (auto nn = startN + 1; nn < n; ++nn)
                {

                    if (m < m_lowerLeft.m_m || m > m_upperRight.m_m || nn < m_lowerLeft.m_n || nn > m_upperRight.m_n)
                    {
                        continue;
                    }
                    if (m_grid.m_gridNodesTypes(m, nn) == CurvilinearGrid::NodeType::Invalid)
                    {
                        continue;
                    }
                    const auto bottomNode = m_grid.m_gridNodes(m, nn - 1);
                    const auto horizontalNode = m_grid.m_gridNodes(m + nextHorizontal, nn);
                    const auto upperNode = m_grid.m_gridNodes(m, nn + 1);

                    Point boundaryNode;
                    if (nextHorizontal == 1)
                    {
                        const auto qb = 1.0 / m_terms.atp(m, nn - 1);
                        const auto qc = 1.0 / m_terms.atp(m, nn);
                        const auto qbc = m_terms.atp(m, nn - 1) + m_terms.atp(m, nn);
                        const auto rn = qb + qc + qbc;
                        boundaryNode.x = (bottomNode.x * qb + horizontalNode.x * qbc + upperNode.x * qc + bottomNode.y - upperNode.y) / rn;
                        boundaryNode.y = (bottomNode.y * qb + horizontalNode.y * qbc + upperNode.y * qc + upperNode.x - bottomNode.x) / rn;
                    }

                    if (nextHorizontal == -1)
                    {
                        const auto qb = 1.0 / m_terms.atp(m - 1, nn - 1);
                        const auto qc = 1.0 / m_terms.atp(m - 1, nn);
                        const auto qbc = m_terms.atp(m - 1, nn - 1) + m_terms.atp(m - 1, nn);
                        const auto rn = qb + qc + qbc;
                        boundaryNode.x = (bottomNode.x * qb + horizontalNode.x * qbc + upperNode.x * qc + upperNode.y - bottomNode.y) / rn;
                        boundaryNode.y = (bottomNode.y * qb + horizontalNode.y * qbc + upperNode.y * qc + bottomNode.x - upperNode.x) / rn;
                    }

                    // Vertical spline index
                    const auto splineIndex = m_grid.m_numN + m;
                    m_grid.m_gridNodes(m, nn) = m_splines.ComputeClosestPointOnSplineSegment(splineIndex,
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
    const auto minMInternal = std::max(static_cast<UInt>(1), m_lowerLeft.m_m);
    const auto minNInternal = std::max(static_cast<UInt>(1), m_lowerLeft.m_n);

    const auto maxMInternal = std::min(m_upperRight.m_m, m_grid.m_numM - 1);
    const auto maxNInternal = std::min(m_upperRight.m_n, m_grid.m_numN - 1);

    for (auto innerIterations = 0; innerIterations < m_orthogonalizationParameters.inner_iterations; ++innerIterations)
    {
        for (auto m = minMInternal; m < maxMInternal; ++m)
        {
            for (auto n = minNInternal; n < maxNInternal; ++n)
            {
                if (m_grid.m_gridNodesTypes(m, n) != CurvilinearGrid::NodeType::InternalValid)
                {
                    continue;
                }

                if (m_isGridNodeFrozen(m, n))
                {
                    continue;
                }

                const auto residual =
                    m_grid.m_gridNodes(m + 1, n) * m_terms.a(m, n) +
                    m_grid.m_gridNodes(m - 1, n) * m_terms.b(m, n) +
                    m_grid.m_gridNodes(m, n + 1) * m_terms.c(m, n) +
                    m_grid.m_gridNodes(m, n - 1) * m_terms.d(m, n) +
                    m_grid.m_gridNodes(m, n) * m_terms.e(m, n);

                m_grid.m_gridNodes(m, n) = m_grid.m_gridNodes(m, n) - residual / m_terms.e(m, n) * omega;
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

void CurvilinearGridOrthogonalization::ComputeCoefficients()
{
    m_terms.Invalidate();

    for (auto m = m_lowerLeft.m_m; m < m_upperRight.m_m; ++m)
    {
        for (auto n = m_lowerLeft.m_n; n < m_upperRight.m_n; ++n)
        {
            if (!m_grid.m_gridFacesMask(m, n))
            {
                continue;
            }

            const auto bottom = ComputeDistance(m_grid.m_gridNodes(m, n), m_grid.m_gridNodes(m + 1, n), Projection::cartesian);
            const auto upper = ComputeDistance(m_grid.m_gridNodes(m, n + 1), m_grid.m_gridNodes(m + 1, n + 1), Projection::cartesian);
            const auto left = ComputeDistance(m_grid.m_gridNodes(m, n), m_grid.m_gridNodes(m, n + 1), Projection::cartesian);
            const auto right = ComputeDistance(m_grid.m_gridNodes(m + 1, n), m_grid.m_gridNodes(m + 1, n + 1), Projection::cartesian);

            m_terms.a(m, n) = (bottom + upper) * 0.5;
            m_terms.b(m, n) = (left + right) * 0.5;

            m_terms.atp(m, n) = m_terms.a(m, n);
            m_terms.e(m, n) = m_terms.b(m, n);

            m_terms.c(m, n) = m_terms.atp(m, n);
            m_terms.d(m, n) = m_terms.e(m, n);
        }
    }

    ComputeHorizontalCoefficients();
    ComputeVerticalCoefficients();

    // Normalize
    const auto smoothingFactor = 1.0 - m_orthogonalizationParameters.orthogonalization_to_smoothing_factor;
    for (auto m = m_lowerLeft.m_m; m < m_upperRight.m_m; ++m)
    {
        for (auto n = m_lowerLeft.m_n; n < m_upperRight.m_n; ++n)
        {
            if (!m_grid.m_gridFacesMask(m, n))
            {
                continue;
            }
            m_terms.atp(m, n) = m_terms.atp(m, n) * m_terms.a(m, n) / m_terms.c(m, n);
            m_terms.atp(m, n) = m_orthogonalizationParameters.orthogonalization_to_smoothing_factor * m_terms.atp(m, n) +
                                smoothingFactor * m_terms.a(m, n);
            m_terms.e(m, n) = m_terms.e(m, n) * m_terms.b(m, n) / m_terms.d(m, n);
            m_terms.e(m, n) = m_orthogonalizationParameters.orthogonalization_to_smoothing_factor * m_terms.e(m, n) +
                              smoothingFactor * m_terms.b(m, n);

            m_terms.a(m, n) = m_terms.atp(m, n);
            m_terms.b(m, n) = m_terms.e(m, n);
        }
    }

    // Calculate m_terms.atp
    for (auto m = m_lowerLeft.m_m; m < m_upperRight.m_m; ++m)
    {
        for (auto n = m_lowerLeft.m_n; n < m_upperRight.m_n; ++n)
        {
            if (!m_grid.m_gridFacesMask(m, n))
            {
                m_terms.atp(m, n) = constants::missing::doubleValue;
                continue;
            }
            m_terms.atp(m, n) = m_terms.b(m, n) / m_terms.a(m, n);
        }
    }

    // Re-set coefficients, preserve atp
    m_terms.Reset();

    for (auto m = m_lowerLeft.m_m + 1; m < m_upperRight.m_m; ++m)
    {
        for (auto n = m_lowerLeft.m_n + 1; n < m_upperRight.m_n; ++n)
        {
            if (m_grid.m_gridNodesTypes(m, n) != CurvilinearGrid::NodeType::InternalValid)
            {
                continue;
            }
            m_terms.a(m, n) = m_terms.atp(m, n - 1) + m_terms.atp(m, n);
            m_terms.b(m, n) = m_terms.atp(m - 1, n - 1) + m_terms.atp(m - 1, n);
            m_terms.c(m, n) = 1.0 / m_terms.atp(m - 1, n) + 1.0 / m_terms.atp(m, n);
            m_terms.d(m, n) = 1.0 / m_terms.atp(m - 1, n - 1) + 1.0 / m_terms.atp(m, n - 1);
            m_terms.e(m, n) = -m_terms.a(m, n) - m_terms.b(m, n) - m_terms.c(m, n) - m_terms.d(m, n);
        }
    }
}

void CurvilinearGridOrthogonalization::ComputeVerticalCoefficients()
{
    const auto invalidBoundaryNodes = ComputeInvalidVerticalBoundaryNodes();
    // Store the counter
    lin_alg::MatrixRowMajor<UInt> counter(m_grid.m_numM, m_grid.m_numN);
    counter.fill(0);

    // Perform left sum
    for (auto n = m_lowerLeft.m_n; n < m_upperRight.m_n; ++n)
    {
        for (auto m = m_lowerLeft.m_m + 1; m < m_upperRight.m_m; ++m)
        {

            if (m_grid.IsValidFace(m, n) &&
                !IsEqual(m_terms.a(m, n), constants::missing::doubleValue) &&
                !IsEqual(m_terms.a(m - 1, n), constants::missing::doubleValue) &&
                !invalidBoundaryNodes[m][n])
            {
                m_terms.a(m, n) += m_terms.a(m - 1, n);
                m_terms.c(m, n) += m_terms.c(m - 1, n);
                counter(m, n) = counter(m - 1, n) + 1;
            }
        }
    }

    // Perform right sum
    for (auto n = m_lowerLeft.m_n; n < m_upperRight.m_n; ++n)
    {
        for (auto m = int(m_upperRight.m_m) - 1; m >= int(m_lowerLeft.m_m); --m)
        {
            if (m_grid.IsValidFace(m, n) &&
                !IsEqual(m_terms.a(m, n), constants::missing::doubleValue) &&
                !IsEqual(m_terms.a(m + 1, n), constants::missing::doubleValue) &&
                !invalidBoundaryNodes[m + 1][n])
            {
                m_terms.a(m, n) = m_terms.a(m + 1, n);
                m_terms.c(m, n) = m_terms.c(m + 1, n);
                counter(m, n) = counter(m + 1, n);
            }
        }
    }
    // Average contributions
    for (auto m = m_lowerLeft.m_m; m < m_upperRight.m_m; ++m)
    {
        for (auto n = m_lowerLeft.m_n; n < m_upperRight.m_n; ++n)
        {
            if (m_grid.IsValidFace(m, n))
            {
                m_terms.a(m, n) /= static_cast<double>(counter(m, n) + 1);
                m_terms.c(m, n) /= static_cast<double>(counter(m, n) + 1);
            }
        }
    }
}

void CurvilinearGridOrthogonalization::ComputeHorizontalCoefficients()
{
    const auto invalidBoundaryNodes = ComputeInvalidHorizontalBoundaryNodes();
    lin_alg::MatrixRowMajor<UInt> counter(m_grid.m_numM, m_grid.m_numN);
    counter.fill(0);

    // Perform bottom sum
    for (auto m = m_lowerLeft.m_m; m < m_upperRight.m_m; ++m)
    {
        for (auto n = m_lowerLeft.m_n + 1; n < m_upperRight.m_n; ++n)
        {
            if (m_grid.IsValidFace(m, n) &&
                !IsEqual(m_terms.b(m, n), constants::missing::doubleValue) &&
                !IsEqual(m_terms.b(m, n - 1), constants::missing::doubleValue) &&
                !invalidBoundaryNodes[m][n])
            {
                m_terms.b(m, n) += m_terms.b(m, n - 1);
                m_terms.d(m, n) += m_terms.d(m, n - 1);
                counter(m, n) = counter(m, n - 1) + 1;
            }
        }
    }

    // Perform upper sum
    for (auto m = m_lowerLeft.m_m; m < m_upperRight.m_m; ++m)
    {
        for (auto n = static_cast<int>(m_upperRight.m_n) - 1; n >= static_cast<int>(m_lowerLeft.m_n); --n)
        {
            if (m_grid.IsValidFace(m, n) &&
                !IsEqual(m_terms.b(m, n), constants::missing::doubleValue) &&
                !IsEqual(m_terms.b(m, n + 1), constants::missing::doubleValue) &&
                !invalidBoundaryNodes[m][n + 1])
            {
                m_terms.b(m, n) = m_terms.b(m, n + 1);
                m_terms.d(m, n) = m_terms.d(m, n + 1);
                counter(m, n) = counter(m, n + 1);
            }
        }
    }

    // Average contributions
    for (auto m = m_lowerLeft.m_m; m < m_upperRight.m_m; ++m)
    {
        for (auto n = m_lowerLeft.m_n; n < m_upperRight.m_n; ++n)
        {
            if (m_grid.IsValidFace(m, n))
            {
                m_terms.b(m, n) /= static_cast<double>(counter(m, n) + 1);
                m_terms.d(m, n) /= static_cast<double>(counter(m, n) + 1);
            }
        }
    }
}

std::vector<std::vector<bool>>
CurvilinearGridOrthogonalization::ComputeInvalidHorizontalBoundaryNodes() const
{
    std::vector<std::vector<bool>> invalidBoundaryNodes(m_grid.m_numM, std::vector<bool>(m_grid.m_numN, false));
    for (auto m = m_lowerLeft.m_m + 1; m < m_upperRight.m_m; ++m)
    {
        for (auto n = m_lowerLeft.m_n + 1; n < m_upperRight.m_n; ++n)
        {
            int step = 0;
            if (m_grid.m_gridNodesTypes(m, n) == CurvilinearGrid::NodeType::BottomLeft)
            {
                step = -1;
            }
            if (m_grid.m_gridNodesTypes(m, n) == CurvilinearGrid::NodeType::BottomRight)
            {
                step = 1;
            }
            if (m_grid.m_gridNodesTypes(m, n) == CurvilinearGrid::NodeType::UpperRight)
            {
                step = 1;
            }
            if (m_grid.m_gridNodesTypes(m, n) == CurvilinearGrid::NodeType::UpperLeft)
            {
                step = -1;
            }
            if (step == 0)
            {
                continue;
            }

            auto lastValidM = m + step;
            while (lastValidM > 0 &&
                   lastValidM < m_grid.m_numM &&
                   m_grid.m_gridNodesTypes(lastValidM, n) == CurvilinearGrid::NodeType::InternalValid)
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
    std::vector<std::vector<bool>> invalidBoundaryNodes(m_grid.m_numM, std::vector<bool>(m_grid.m_numN, false));
    for (auto m = m_lowerLeft.m_m + 1; m < m_upperRight.m_m; ++m)
    {
        for (auto n = m_lowerLeft.m_n + 1; n < m_upperRight.m_n; ++n)
        {
            int step = 0;
            if (m_grid.m_gridNodesTypes(m, n) == CurvilinearGrid::NodeType::BottomLeft)
            {
                step = -1;
            }
            if (m_grid.m_gridNodesTypes(m, n) == CurvilinearGrid::NodeType::BottomRight)
            {
                step = -1;
            }
            if (m_grid.m_gridNodesTypes(m, n) == CurvilinearGrid::NodeType::UpperRight)
            {
                step = 1;
            }
            if (m_grid.m_gridNodesTypes(m, n) == CurvilinearGrid::NodeType::UpperLeft)
            {
                step = 1;
            }
            if (step == 0)
            {
                continue;
            }
            auto lastValidN = n + step;
            while (lastValidN > 0 &&
                   lastValidN < m_grid.m_numN &&
                   m_grid.m_gridNodesTypes(m, lastValidN) == CurvilinearGrid::NodeType::InternalValid)
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
