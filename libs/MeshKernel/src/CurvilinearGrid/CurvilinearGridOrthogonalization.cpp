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

CurvilinearGridOrthogonalization::CurvilinearGridOrthogonalization(CurvilinearGrid& grid,
                                                                   const OrthogonalizationParameters& orthogonalizationParameters)
    : CurvilinearGridAlgorithm(grid),
      m_orthoEqTerms(m_grid.NumN(), m_grid.NumM()),
      m_isGridNodeFrozen(m_grid.NumN(), m_grid.NumM()),
      m_splines(Splines(m_grid))
{
    CheckOrthogonalizationParameters(orthogonalizationParameters);
    m_orthogonalizationParameters = orthogonalizationParameters;
    m_isGridNodeFrozen.fill(false);
}

void CurvilinearGridOrthogonalization::ComputeFrozenGridPoints()
{
    for (auto const& frozenLine : m_lines)
    {
        for (auto n = frozenLine.m_startNode.m_n; n <= frozenLine.m_endNode.m_n; ++n)
        {
            for (auto m = frozenLine.m_startNode.m_m; m <= frozenLine.m_endNode.m_m; ++m)
            {
                m_isGridNodeFrozen(n, m) = true;
            }
        }
    }
}

void CurvilinearGridOrthogonalization::Compute()
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
}

void CurvilinearGridOrthogonalization::ProjectHorizontalBoundaryGridNodes()
{
    // m grid lines (horizontal)
    for (UInt m = 0; m < m_grid.NumM(); ++m)
    {
        UInt startN = constants::missing::uintValue;
        int nextVertical = 0;
        for (UInt n = 0; n < m_grid.NumN(); ++n)
        {
            const auto nodeType = m_grid.GetNodeType(n, m);
            if (nodeType == CurvilinearGrid::NodeType::BottomLeft || nodeType == CurvilinearGrid::NodeType::UpperLeft)
            {
                startN = n;
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
            if (startN != constants::missing::uintValue &&
                (nodeType == CurvilinearGrid::NodeType::BottomRight || nodeType == CurvilinearGrid::NodeType::UpperRight) &&
                nextVertical != 0)
            {
                for (auto nn = startN + 1; nn < n; ++nn)
                {

                    if (m < m_lowerLeft.m_m || m > m_upperRight.m_m || nn < m_lowerLeft.m_n || nn > m_upperRight.m_n)
                    {
                        continue;
                    }
                    if (m_grid.GetNodeType(n, m) == CurvilinearGrid::NodeType::Invalid)
                    {
                        continue;
                    }

                    const auto leftNode = m_grid.GetNode(nn - 1, m);
                    const auto verticalNode = m_grid.GetNode(nn, m + nextVertical);
                    const auto rightNode = m_grid.GetNode(nn + 1, m);

                    Point boundaryNode;
                    if (nextVertical == 1)
                    {
                        const double qb = m_orthoEqTerms.atp(nn - 1, m);
                        const double qc = m_orthoEqTerms.atp(nn, m);
                        const auto qbc = 1.0 / qb + 1.0 / qc;
                        const auto rn = qb + qc + qbc;
                        boundaryNode.x = (leftNode.x * qb + verticalNode.x * qbc + rightNode.x * qc + rightNode.y - leftNode.y) / rn;
                        boundaryNode.y = (leftNode.y * qb + verticalNode.y * qbc + rightNode.y * qc + leftNode.x - rightNode.x) / rn;
                    }

                    if (nextVertical == -1)
                    {
                        const double qb = m_orthoEqTerms.atp(nn - 1, m - 1);
                        const double qc = m_orthoEqTerms.atp(nn, m - 1);
                        const auto qbc = 1.0 / qb + 1.0 / qc;
                        const auto rn = qb + qc + qbc;
                        boundaryNode.x = (leftNode.x * qb + verticalNode.x * qbc + rightNode.x * qc + leftNode.y - rightNode.y) / rn;
                        boundaryNode.y = (leftNode.y * qb + verticalNode.y * qbc + rightNode.y * qc + rightNode.x - leftNode.x) / rn;
                    }

                    const auto grid_node_value = m_splines.ComputeClosestPointOnSplineSegment(m,
                                                                                              startN,
                                                                                              n,
                                                                                              boundaryNode);

                    m_grid.GetNode(nn, m) = grid_node_value;
                }
            }
        }
    }
}

void CurvilinearGridOrthogonalization::ProjectVerticalBoundariesGridNodes()
{
    // m gridlines (vertical)
    for (UInt n = 0; n < m_grid.NumN(); ++n)
    {
        UInt startM = constants::missing::uintValue;
        int nextHorizontal = 0;
        for (UInt m = 0; m < m_grid.NumM(); ++m)
        {
            const auto nodeType = m_grid.GetNodeType(n, m);
            if (nodeType == CurvilinearGrid::NodeType::BottomLeft || nodeType == CurvilinearGrid::NodeType::BottomRight)
            {
                startM = m;
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
                startM != constants::missing::uintValue)
            {
                for (auto mm = startM + 1; mm < m; ++mm)
                {

                    if (mm < m_lowerLeft.m_m || mm > m_upperRight.m_m || n < m_lowerLeft.m_n || n > m_upperRight.m_n)
                    {
                        continue;
                    }
                    if (m_grid.GetNodeType(n, m) == CurvilinearGrid::NodeType::Invalid)
                    {
                        continue;
                    }
                    const auto bottomNode = m_grid.GetNode(n, mm - 1);
                    const auto horizontalNode = m_grid.GetNode(n + nextHorizontal, mm);
                    const auto upperNode = m_grid.GetNode(n, mm + 1);

                    Point boundaryNode;
                    if (nextHorizontal == 1)
                    {
                        const auto qb = 1.0 / m_orthoEqTerms.atp(n, mm - 1);
                        const auto qc = 1.0 / m_orthoEqTerms.atp(n, mm);
                        const auto qbc = m_orthoEqTerms.atp(n, mm - 1) + m_orthoEqTerms.atp(n, mm);
                        const auto rn = qb + qc + qbc;
                        boundaryNode.x = (bottomNode.x * qb + horizontalNode.x * qbc + upperNode.x * qc + bottomNode.y - upperNode.y) / rn;
                        boundaryNode.y = (bottomNode.y * qb + horizontalNode.y * qbc + upperNode.y * qc + upperNode.x - bottomNode.x) / rn;
                    }

                    if (nextHorizontal == -1)
                    {
                        const auto qb = 1.0 / m_orthoEqTerms.atp(n - 1, mm - 1);
                        const auto qc = 1.0 / m_orthoEqTerms.atp(n - 1, mm);
                        const auto qbc = m_orthoEqTerms.atp(n - 1, mm - 1) + m_orthoEqTerms.atp(n - 1, mm);
                        const auto rn = qb + qc + qbc;
                        boundaryNode.x = (bottomNode.x * qb + horizontalNode.x * qbc + upperNode.x * qc + upperNode.y - bottomNode.y) / rn;
                        boundaryNode.y = (bottomNode.y * qb + horizontalNode.y * qbc + upperNode.y * qc + bottomNode.x - upperNode.x) / rn;
                    }

                    // Vertical spline index
                    const auto splineIndex = m_grid.NumM() + n;

                    const auto node_value = m_splines.ComputeClosestPointOnSplineSegment(splineIndex,
                                                                                         startM,
                                                                                         m,
                                                                                         boundaryNode);

                    m_grid.GetNode(n, mm) = node_value;
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

    const auto maxMInternal = std::min(m_upperRight.m_m, m_grid.NumM() - 1);
    const auto maxNInternal = std::min(m_upperRight.m_n, m_grid.NumN() - 1);

    for (auto innerIterations = 0; innerIterations < m_orthogonalizationParameters.inner_iterations; ++innerIterations)
    {
        for (auto n = minNInternal; n < maxNInternal; ++n)
        {
            for (auto m = minMInternal; m < maxMInternal; ++m)
            {
                if (m_grid.GetNodeType(n, m) != CurvilinearGrid::NodeType::InternalValid)
                {
                    continue;
                }

                if (m_isGridNodeFrozen(n, m))
                {
                    continue;
                }

                const auto residual =
                    m_grid.GetNode(n + 1, m) * m_orthoEqTerms.a(n, m) +
                    m_grid.GetNode(n - 1, m) * m_orthoEqTerms.b(n, m) +
                    m_grid.GetNode(n, m + 1) * m_orthoEqTerms.c(n, m) +
                    m_grid.GetNode(n, m - 1) * m_orthoEqTerms.d(n, m) +
                    m_grid.GetNode(n, m) * m_orthoEqTerms.e(n, m);

                m_grid.GetNode(n, m) = m_grid.GetNode(n, m) - residual / m_orthoEqTerms.e(n, m) * omega;
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
    // initialise matrix coefficients
    m_orthoEqTerms.a.fill(constants::missing::doubleValue);
    m_orthoEqTerms.b.fill(constants::missing::doubleValue);
    m_orthoEqTerms.c.fill(constants::missing::doubleValue);
    m_orthoEqTerms.d.fill(constants::missing::doubleValue);
    m_orthoEqTerms.e.fill(constants::missing::doubleValue);
    m_orthoEqTerms.atp.fill(constants::missing::doubleValue);

    for (auto n = m_lowerLeft.m_n; n < m_upperRight.m_n; ++n)
    {
        for (auto m = m_lowerLeft.m_m; m < m_upperRight.m_m; ++m)
        {

            if (!m_grid.IsFaceMaskValid(n, m))
            {
                continue;
            }

            const auto bottom = ComputeDistance(m_grid.GetNode(n, m), m_grid.GetNode(n + 1, m), Projection::cartesian);
            const auto upper = ComputeDistance(m_grid.GetNode(n, m + 1), m_grid.GetNode(n + 1, m + 1), Projection::cartesian);
            const auto left = ComputeDistance(m_grid.GetNode(n, m), m_grid.GetNode(n, m + 1), Projection::cartesian);
            const auto right = ComputeDistance(m_grid.GetNode(n + 1, m), m_grid.GetNode(n + 1, m + 1), Projection::cartesian);

            m_orthoEqTerms.a(n, m) = (bottom + upper) * 0.5;
            m_orthoEqTerms.b(n, m) = (left + right) * 0.5;

            m_orthoEqTerms.atp(n, m) = m_orthoEqTerms.a(n, m);
            m_orthoEqTerms.e(n, m) = m_orthoEqTerms.b(n, m);

            m_orthoEqTerms.c(n, m) = m_orthoEqTerms.atp(n, m);
            m_orthoEqTerms.d(n, m) = m_orthoEqTerms.e(n, m);
        }
    }

    ComputeHorizontalCoefficients();
    ComputeVerticalCoefficients();

    // Normalize
    const auto smoothingFactor = 1.0 - m_orthogonalizationParameters.orthogonalization_to_smoothing_factor;
    for (auto n = m_lowerLeft.m_n; n < m_upperRight.m_n; ++n)
    {
        for (auto m = m_lowerLeft.m_m; m < m_upperRight.m_m; ++m)
        {
            if (!m_grid.IsFaceMaskValid(n, m))
            {
                continue;
            }
            m_orthoEqTerms.atp(n, m) = m_orthoEqTerms.atp(n, m) * m_orthoEqTerms.a(n, m) / m_orthoEqTerms.c(n, m);
            m_orthoEqTerms.atp(n, m) = m_orthogonalizationParameters.orthogonalization_to_smoothing_factor * m_orthoEqTerms.atp(n, m) + smoothingFactor * m_orthoEqTerms.a(n, m);
            m_orthoEqTerms.e(n, m) = m_orthoEqTerms.e(n, m) * m_orthoEqTerms.b(n, m) / m_orthoEqTerms.d(n, m);
            m_orthoEqTerms.e(n, m) = m_orthogonalizationParameters.orthogonalization_to_smoothing_factor * m_orthoEqTerms.e(n, m) + smoothingFactor * m_orthoEqTerms.b(n, m);

            m_orthoEqTerms.a(n, m) = m_orthoEqTerms.atp(n, m);
            m_orthoEqTerms.b(n, m) = m_orthoEqTerms.e(n, m);
        }
    }

    // Calculate m_atp
    for (auto n = m_lowerLeft.m_n; n < m_upperRight.m_n; ++n)
    {
        for (auto m = m_lowerLeft.m_m; m < m_upperRight.m_m; ++m)
        {
            if (!m_grid.IsFaceMaskValid(n, m))
            {
                m_orthoEqTerms.atp(n, m) = constants::missing::doubleValue;
                continue;
            }
            m_orthoEqTerms.atp(n, m) = m_orthoEqTerms.b(n, m) / m_orthoEqTerms.a(n, m);
        }
    }

    // Re-set coefficients
    m_orthoEqTerms.a.fill(0.0);
    m_orthoEqTerms.b.fill(0.0);
    m_orthoEqTerms.c.fill(0.0);
    m_orthoEqTerms.d.fill(0.0);
    m_orthoEqTerms.e.fill(0.0);

    for (auto n = m_lowerLeft.m_n + 1; n < m_upperRight.m_n; ++n)
    {
        for (auto m = m_lowerLeft.m_m + 1; m < m_upperRight.m_m; ++m)
        {
            if (m_grid.GetNodeType(n, m) != CurvilinearGrid::NodeType::InternalValid)
            {
                continue;
            }
            m_orthoEqTerms.a(n, m) = m_orthoEqTerms.atp(n, m - 1) + m_orthoEqTerms.atp(n, m);
            m_orthoEqTerms.b(n, m) = m_orthoEqTerms.atp(n - 1, m - 1) + m_orthoEqTerms.atp(n - 1, m);
            m_orthoEqTerms.c(n, m) = 1.0 / m_orthoEqTerms.atp(n - 1, m) + 1.0 / m_orthoEqTerms.atp(n, m);
            m_orthoEqTerms.d(n, m) = 1.0 / m_orthoEqTerms.atp(n - 1, m - 1) + 1.0 / m_orthoEqTerms.atp(n, m - 1);

            m_orthoEqTerms.e(n, m) = -m_orthoEqTerms.a(n, m) - m_orthoEqTerms.b(n, m) - m_orthoEqTerms.c(n, m) - m_orthoEqTerms.d(n, m);
        }
    }
}

void CurvilinearGridOrthogonalization::ComputeVerticalCoefficients()
{
    const auto invalidBoundaryNodes = ComputeInvalidVerticalBoundaryNodes();
    // Store the counter
    lin_alg::Matrix<UInt> counter(m_grid.NumN(), m_grid.NumM());
    counter.fill(0);

    // Perform left sum
    for (auto n = m_lowerLeft.m_n + 1; n < m_upperRight.m_n; ++n)
    {
        for (auto m = m_lowerLeft.m_m; m < m_upperRight.m_m; ++m)
        {

            if (m_grid.AreFaceNodesValid(n, m) &&
                !IsEqual(m_orthoEqTerms.a(n, m), constants::missing::doubleValue) &&
                !IsEqual(m_orthoEqTerms.a(n - 1, m), constants::missing::doubleValue) &&
                !invalidBoundaryNodes(n, m))
            {
                m_orthoEqTerms.a(n, m) += m_orthoEqTerms.a(n - 1, m);
                m_orthoEqTerms.c(n, m) += m_orthoEqTerms.c(n - 1, m);
                counter(n, m) = counter(n - 1, m) + 1;
            }
        }
    }

    // Perform right sum
    for (auto n = static_cast<int>(m_upperRight.m_n) - 1; n >= static_cast<int>(m_lowerLeft.m_n); --n)
    {
        for (auto m = m_lowerLeft.m_m; m < m_upperRight.m_m; ++m)
        {
            if (m_grid.AreFaceNodesValid(n, m) &&
                !IsEqual(m_orthoEqTerms.a(n, m), constants::missing::doubleValue) &&
                !IsEqual(m_orthoEqTerms.a(n + 1, m), constants::missing::doubleValue) &&
                !invalidBoundaryNodes(n + 1, m))
            {
                m_orthoEqTerms.a(n, m) = m_orthoEqTerms.a(n + 1, m);
                m_orthoEqTerms.c(n, m) = m_orthoEqTerms.c(n + 1, m);
                counter(n, m) = counter(n + 1, m);
            }
        }
    }
    // Average contributions
    for (auto n = m_lowerLeft.m_n; n < m_upperRight.m_n; ++n)
    {
        for (auto m = m_lowerLeft.m_m; m < m_upperRight.m_m; ++m)
        {
            if (m_grid.AreFaceNodesValid(n, m))
            {
                double const inv = 1.0 / static_cast<double>(counter(n, m) + 1);
                m_orthoEqTerms.a(n, m) *= inv;
                m_orthoEqTerms.c(n, m) *= inv;
            }
        }
    }
}

void CurvilinearGridOrthogonalization::ComputeHorizontalCoefficients()
{
    const auto invalidBoundaryNodes = ComputeInvalidHorizontalBoundaryNodes();
    lin_alg::Matrix<UInt> counter(m_grid.NumN(), m_grid.NumM());
    counter.fill(0);

    // Perform bottom sum
    for (auto n = m_lowerLeft.m_n; n < m_upperRight.m_n; ++n)
    {
        for (auto m = m_lowerLeft.m_m + 1; m < m_upperRight.m_m; ++m)
        {
            if (m_grid.AreFaceNodesValid(n, m) &&
                !IsEqual(m_orthoEqTerms.b(n, m), constants::missing::doubleValue) &&
                !IsEqual(m_orthoEqTerms.b(n, m - 1), constants::missing::doubleValue) &&
                !invalidBoundaryNodes(n, m))
            {
                m_orthoEqTerms.b(n, m) += m_orthoEqTerms.b(n, m - 1);
                m_orthoEqTerms.d(n, m) += m_orthoEqTerms.d(n, m - 1);
                counter(n, m) = counter(n, m - 1) + 1;
            }
        }
    }

    // Perform upper sum
    for (auto n = m_lowerLeft.m_n; n < m_upperRight.m_n; ++n)
    {
        for (auto m = static_cast<int>(m_upperRight.m_m) - 1; m >= static_cast<int>(m_lowerLeft.m_m); --m)
        {
            if (m_grid.AreFaceNodesValid(n, m) &&
                !IsEqual(m_orthoEqTerms.b(n, m), constants::missing::doubleValue) &&
                !IsEqual(m_orthoEqTerms.b(n, m + 1), constants::missing::doubleValue) &&
                !invalidBoundaryNodes(n, m + 1))
            {
                m_orthoEqTerms.b(n, m) = m_orthoEqTerms.b(n, m + 1);
                m_orthoEqTerms.d(n, m) = m_orthoEqTerms.d(n, m + 1);
                counter(n, m) = counter(n, m + 1);
            }
        }
    }

    // Average contributions
    for (auto n = m_lowerLeft.m_n; n < m_upperRight.m_n; ++n)
    {
        for (auto m = m_lowerLeft.m_m; m < m_upperRight.m_m; ++m)
        {
            if (m_grid.AreFaceNodesValid(n, m))
            {
                double const inv = 1.0 / static_cast<double>(counter(n, m) + 1);
                m_orthoEqTerms.b(n, m) *= inv;
                m_orthoEqTerms.d(n, m) *= inv;
            }
        }
    }
}

lin_alg::Matrix<bool> CurvilinearGridOrthogonalization::ComputeInvalidHorizontalBoundaryNodes() const
{
    lin_alg::Matrix<bool> invalidBoundaryNodes(m_grid.NumN(), m_grid.NumM());
    invalidBoundaryNodes.fill(false);

    for (auto n = m_lowerLeft.m_n + 1; n < m_upperRight.m_n; ++n)
    {
        for (auto m = m_lowerLeft.m_m + 1; m < m_upperRight.m_m; ++m)
        {
            int step = 0;
            if (m_grid.GetNodeType(n, m) == CurvilinearGrid::NodeType::BottomLeft)
            {
                step = -1;
            }
            if (m_grid.GetNodeType(n, m) == CurvilinearGrid::NodeType::BottomRight)
            {
                step = 1;
            }
            if (m_grid.GetNodeType(n, m) == CurvilinearGrid::NodeType::UpperRight)
            {
                step = 1;
            }
            if (m_grid.GetNodeType(n, m) == CurvilinearGrid::NodeType::UpperLeft)
            {
                step = -1;
            }
            if (step == 0)
            {
                continue;
            }

            auto lastValidN = n + step;
            while (lastValidN > 0 &&
                   lastValidN < m_grid.NumN() &&
                   m_grid.GetNodeType(lastValidN, n) == CurvilinearGrid::NodeType::InternalValid)
            {
                lastValidN += step;
            }
            const auto start = std::min(lastValidN, n);
            const auto end = std::max(lastValidN, n);
            for (auto k = start; k <= end; ++k)
            {
                invalidBoundaryNodes(k, m) = true;
            }
        }
    }

    return invalidBoundaryNodes;
}

lin_alg::Matrix<bool> CurvilinearGridOrthogonalization::ComputeInvalidVerticalBoundaryNodes() const
{
    lin_alg::Matrix<bool> invalidBoundaryNodes(m_grid.NumN(), m_grid.NumM());
    invalidBoundaryNodes.fill(false);

    for (auto n = m_lowerLeft.m_n + 1; n < m_upperRight.m_n; ++n)
    {
        for (auto m = m_lowerLeft.m_m + 1; m < m_upperRight.m_m; ++m)
        {
            int step = 0;
            if (m_grid.GetNodeType(n, m) == CurvilinearGrid::NodeType::BottomLeft)
            {
                step = -1;
            }
            if (m_grid.GetNodeType(n, m) == CurvilinearGrid::NodeType::BottomRight)
            {
                step = -1;
            }
            if (m_grid.GetNodeType(n, m) == CurvilinearGrid::NodeType::UpperRight)
            {
                step = 1;
            }
            if (m_grid.GetNodeType(n, m) == CurvilinearGrid::NodeType::UpperLeft)
            {
                step = 1;
            }
            if (step == 0)
            {
                continue;
            }
            auto lastValidM = m + step;
            while (lastValidM > 0 &&
                   lastValidM < m_grid.NumM() &&
                   m_grid.GetNodeType(n, lastValidM) == CurvilinearGrid::NodeType::InternalValid)
            {
                lastValidM += step;
            }
            const auto start = std::min(lastValidM, m);
            const auto end = std::max(m, lastValidM);
            for (auto k = start; k <= end; ++k)
            {
                invalidBoundaryNodes(n, k) = true;
            }
        }
    }

    return invalidBoundaryNodes;
}
