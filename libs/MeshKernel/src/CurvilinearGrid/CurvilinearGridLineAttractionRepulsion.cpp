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
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLine.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLineAttractionRepulsion.hpp>
#include <MeshKernel/CurvilinearGrid/UndoActions/CurvilinearGridBlockUndoAction.hpp>
#include <MeshKernel/Entities.hpp>

using meshkernel::CurvilinearGrid;
using meshkernel::CurvilinearGridLineAttractionRepulsion;
using meshkernel::Point;

CurvilinearGridLineAttractionRepulsion::CurvilinearGridLineAttractionRepulsion(CurvilinearGrid& grid,
                                                                               double attractionFactor) : CurvilinearGridAlgorithm(grid),
                                                                                                          m_originalGrid(grid),
                                                                                                          m_attractionFactor(attractionFactor)

{
}

meshkernel::UndoActionPtr CurvilinearGridLineAttractionRepulsion::Compute()
{
    if (m_grid.IsEmpty())
    {
        throw std::invalid_argument("CurvilinearGridLineAttractionRepulsion::Compute Empty grid");
    }

    if (m_lines.empty())
    {
        throw std::invalid_argument("CurvilinearGridLineAttractionRepulsion::Compute No candidate line to shift has been selected");
    }

    // Points are coinciding, no attraction/repulsion zone defined
    if ((m_lines[0].IsMGridLine() && m_lowerLeft.m_n == m_upperRight.m_n) ||
        (m_lines[0].IsNGridLine() && m_lowerLeft.m_m == m_upperRight.m_m))
    {
        throw std::invalid_argument("CurvilinearGridLineAttractionRepulsion::Compute The points defining the attraction area have the same direction of the attraction line.");
    }

    auto const startN = m_lines[0].IsNGridLine() ? m_lines[0].m_startCoordinate : m_lowerLeft.m_n;
    auto const endN = m_lines[0].IsNGridLine() ? m_lines[0].m_endCoordinate : m_upperRight.m_n;
    auto const startM = m_lines[0].IsMGridLine() ? m_lines[0].m_startCoordinate : m_lowerLeft.m_m;
    auto const endM = m_lines[0].IsMGridLine() ? m_lines[0].m_endCoordinate : m_upperRight.m_m;

    std::unique_ptr<CurvilinearGridBlockUndoAction> undoAction = CurvilinearGridBlockUndoAction::Create(m_grid, {startN, startM}, {endN + 1, endM + 1});

    for (auto n = startN; n <= endN; ++n)
    {
        for (auto m = startM; m <= endM; ++m)
        {
            // Not a valid grid node
            if (!m_originalGrid.GetNode(n, m).IsValid())
            {
                continue;
            }

            CurvilinearGridNodeIndices const nodeIndex{n, m};

            // Nodes on the line should not be moved
            if (m_lines[0].IsNodeOnLine(nodeIndex))
            {
                continue;
            }

            const auto [mSmoothing, nSmoothing, mixedSmoothing] = CurvilinearGrid::ComputeDirectionalSmoothingFactors(nodeIndex, m_lines[0].m_startNode, m_lowerLeft, m_upperRight);

            auto const distance = m_originalGrid.ComputeAverageNodalDistance(nodeIndex, m_lines[0].m_gridLineType);
            auto displacement = Point{0.0, 0.0};

            if (m_lines[0].IsMGridLine())
            {
                double const direction = n < m_lines[0].m_constantCoordinate ? 1.0 : -1.0;
                displacement.x = distance * m_attractionFactor * nSmoothing * direction;
            }
            if (m_lines[0].IsNGridLine())
            {
                double const direction = m < m_lines[0].m_constantCoordinate ? -1.0 : 1.0;
                displacement.y = distance * m_attractionFactor * mSmoothing * direction;
            }

            // project transformation
            displacement = m_originalGrid.TransformDisplacement(displacement, nodeIndex, false);

            // adjust nodes
            m_grid.GetNode(n, m) += displacement;
        }
    }

    return undoAction;
}
