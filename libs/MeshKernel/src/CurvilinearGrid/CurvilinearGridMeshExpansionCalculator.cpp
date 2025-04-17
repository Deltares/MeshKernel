//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2025.
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

#include "MeshKernel/CurvilinearGrid/CurvilinearGridMeshExpansionCalculator.hpp"

#include <cmath>
#include <numbers>

meshkernel::UserDefinedRegionExpasionCalculator::UserDefinedRegionExpasionCalculator(const CurvilinearGridNodeIndices& lowerLeft,
                                                                                     const CurvilinearGridNodeIndices& upperRight,
                                                                                     const CurvilinearGridNodeIndices& regionIndicator) : m_indexBoxLowerLeft(lowerLeft),
                                                                                                                                          m_indexBoxUpperRight(upperRight),
                                                                                                                                          m_expansionRegionIndicator(regionIndicator) {}

double meshkernel::UserDefinedRegionExpasionCalculator::compute(const CurvilinearGridNodeIndices& snappedNodeIndex,
                                                                const CurvilinearGridNodeIndices& gridLinePointIndex) const
{
    auto [scalingH, scalingV, scalingMixed] = CurvilinearGrid::ComputeDirectionalSmoothingFactors(gridLinePointIndex, snappedNodeIndex, m_indexBoxLowerLeft, m_indexBoxUpperRight);

    double factor = 0.0;

    if (m_expansionRegionIndicator.m_m == 1 && m_expansionRegionIndicator.m_n == 1)
    {
        factor = scalingMixed;
    }
    else if (m_expansionRegionIndicator.m_n == 1)
    {
        factor = scalingV;
    }
    else if (m_expansionRegionIndicator.m_m == 1)
    {
        factor = scalingH;
    }

    return factor;
}

meshkernel::DefaultRegionExpasionCalculator::DefaultRegionExpasionCalculator(const CurvilinearGrid& originalGrid,
                                                                             const CurvilinearGrid& snappedGrid,
                                                                             const LandBoundary& landBoundary) : m_originalGrid(originalGrid),
                                                                                                                 m_snappedGrid(snappedGrid),
                                                                                                                 m_expansionRegionMinimum(CalculateExpansionRegion(m_originalGrid.GetBoundingBox(),
                                                                                                                                                                   landBoundary.GetBoundingBox())) {}

meshkernel::DefaultRegionExpasionCalculator::DefaultRegionExpasionCalculator(const CurvilinearGrid& originalGrid,
                                                                             const CurvilinearGrid& snappedGrid,
                                                                             const Splines& spline) : m_originalGrid(originalGrid),
                                                                                                      m_snappedGrid(snappedGrid),
                                                                                                      m_expansionRegionMinimum(CalculateExpansionRegion(m_originalGrid.GetBoundingBox(),
                                                                                                                                                        spline.GetBoundingBox(0))) {}

double meshkernel::DefaultRegionExpasionCalculator::CalculateExpansionRegion(const BoundingBox& gridBoundingBox,
                                                                             const BoundingBox& entityBoundingBox)
{
    Vector delta = Merge(gridBoundingBox, entityBoundingBox).Delta();

    delta *= expansionRegionEnlargementFactor;

    delta.y() = std::max(delta.y(), aspectRatio * delta.x());
    return delta.y() / (6.0 * aspectRatio);
}

double meshkernel::DefaultRegionExpasionCalculator::compute(const CurvilinearGridNodeIndices& snappedNodeIndex,
                                                            const CurvilinearGridNodeIndices& gridLinePointIndex) const
{
    const Point meshDelta = m_snappedGrid.GetNode(snappedNodeIndex) - m_originalGrid.GetNode(snappedNodeIndex);
    const Point pointDelta = m_originalGrid.GetNode(gridLinePointIndex) - m_originalGrid.GetNode(snappedNodeIndex);
    double factor = 0.0;

    double rsx = std::max(m_expansionRegionMinimum, std::hypot(meshDelta.x, meshDelta.y));

    if (double rn = std::hypot(pointDelta.x, pointDelta.y); rn < rsx)
    {
        factor = 0.5 * (1.0 + std::cos(std::numbers::pi * rn / rsx));
    }

    return factor;
}
