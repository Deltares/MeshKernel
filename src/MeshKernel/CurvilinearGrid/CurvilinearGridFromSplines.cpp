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

#include <algorithm>
#include <cassert>
#include <vector>

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridFromSplines.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Splines.hpp>
#include <MeshKernelApi/CurvilinearParameters.hpp>
#include <MeshKernelApi/SplinesToCurvilinearParameters.hpp>

using meshkernel::CurvilinearGrid;
using meshkernel::CurvilinearGridFromSplines;

CurvilinearGridFromSplines::CurvilinearGridFromSplines(std::shared_ptr<Splines> splines,
                                                       const meshkernelapi::CurvilinearParameters& curvilinearParameters,
                                                       const meshkernelapi::SplinesToCurvilinearParameters& splinesToCurvilinearParameters) : m_splines(splines),
                                                                                                                                              m_curvilinearParameters(curvilinearParameters),
                                                                                                                                              m_splinesToCurvilinearParameters(splinesToCurvilinearParameters)
{
    m_onTopOfEachOtherSquaredTolerance = m_splinesToCurvilinearParameters.nodes_on_top_of_each_other_tolerance *
                                         m_splinesToCurvilinearParameters.nodes_on_top_of_each_other_tolerance;
};

/// to be called after all splines have been stored
void CurvilinearGridFromSplines::AllocateSplinesProperties()
{
    const auto numSplines = m_splines->GetNumSplines();
    m_type.resize(numSplines);

    m_centralSplineIndex.resize(numSplines);
    std::fill(m_centralSplineIndex.begin(), m_centralSplineIndex.end(), intMissingValue);

    m_numCrossingSplines.resize(numSplines, 0);
    std::fill(m_numCrossingSplines.begin(), m_numCrossingSplines.end(), 0);

    m_maximumGridHeights.resize(numSplines);
    std::fill(m_maximumGridHeights.begin(), m_maximumGridHeights.end(), doubleMissingValue);

    // multi-dimensional arrays
    m_crossingSplinesIndices.resize(numSplines);
    m_isLeftOriented.resize(numSplines);
    m_crossSplineCoordinates.resize(numSplines);
    m_cosCrossingAngle.resize(numSplines);
    m_crossSplineLeftHeights.resize(numSplines);
    m_crossSplineRightHeights.resize(numSplines);
    m_numCrossSplineLeftHeights.resize(numSplines);
    m_numCrossSplineRightHeights.resize(numSplines);
    for (auto s = 0; s < numSplines; ++s)
    {
        m_crossingSplinesIndices[s].resize(numSplines);
        std::fill(m_crossingSplinesIndices[s].begin(), m_crossingSplinesIndices[s].end(), 0);

        m_isLeftOriented[s].resize(numSplines, true);
        std::fill(m_isLeftOriented[s].begin(), m_isLeftOriented[s].end(), true);

        m_crossSplineCoordinates[s].resize(numSplines);
        std::fill(m_crossSplineCoordinates[s].begin(), m_crossSplineCoordinates[s].end(), doubleMissingValue);

        m_cosCrossingAngle[s].resize(numSplines, doubleMissingValue);
        std::fill(m_cosCrossingAngle[s].begin(), m_cosCrossingAngle[s].end(), doubleMissingValue);

        m_crossSplineLeftHeights[s].resize(numSplines);
        std::fill(m_crossSplineLeftHeights[s].begin(), m_crossSplineLeftHeights[s].end(), std::vector<double>(m_maxNumCenterSplineHeights, doubleMissingValue));

        m_crossSplineRightHeights[s].resize(numSplines);
        std::fill(m_crossSplineRightHeights[s].begin(), m_crossSplineRightHeights[s].end(), std::vector<double>(m_maxNumCenterSplineHeights, doubleMissingValue));

        m_numCrossSplineLeftHeights[s].resize(numSplines);
        std::fill(m_numCrossSplineLeftHeights[s].begin(), m_numCrossSplineLeftHeights[s].end(), 0);

        m_numCrossSplineRightHeights[s].resize(numSplines, intMissingValue);
        std::fill(m_numCrossSplineRightHeights[s].begin(), m_numCrossSplineRightHeights[s].end(), 0);
    }

    m_numMSplines.resize(numSplines);
    std::fill(m_numMSplines.begin(), m_numMSplines.end(), 0);
    m_leftGridLineIndex.resize(numSplines);
    std::fill(m_leftGridLineIndex.begin(), m_leftGridLineIndex.end(), sizetMissingValue);
    m_rightGridLineIndex.resize(numSplines);
    std::fill(m_rightGridLineIndex.begin(), m_rightGridLineIndex.end(), sizetMissingValue);
}

CurvilinearGrid CurvilinearGridFromSplines::Compute()
{
    Initialize();

    // Grow grid, from the second layer
    for (auto layer = 1; layer <= m_curvilinearParameters.n_refinement; ++layer)
    {
        Iterate(layer);
    }

    const auto deleteSkinnyTriangles = m_splinesToCurvilinearParameters.remove_skinny_triangles == 1 ? true : false;
    if (deleteSkinnyTriangles)
    {
        DeleteSkinnyTriangles();
    }

    return ComputeCurvilinearGridFromGridPoints();
}

void CurvilinearGridFromSplines::DeleteSkinnyTriangles()
{
    const size_t numMaxIterations = 10;
    const auto numN = m_gridPoints.size() - 2;
    const double squaredDistanceTolerance = 1e-4;
    const double cosineTolerance = 1e-2;
    const double maxCosine = 0.93969;
    for (auto j = numN - 1; j >= 1; --j)
    {
        for (auto iter = 0; iter < numMaxIterations; ++iter)
        {
            size_t numChanged = 0;
            size_t firstRightIndex = 0;
            size_t i = 0;

            while (firstRightIndex != m_numM - 1 || i != m_numM - 1)
            {
                if (firstRightIndex > i)
                {
                    i = firstRightIndex;
                }
                else
                {
                    i++;
                    if (i >= m_numM - 1)
                    {
                        break;
                    }
                }

                if (!m_gridPoints[j][i].IsValid())
                {
                    continue;
                }

                auto [firstLeftIndex, firstRightIndex] = GetNeighbours(m_gridPoints[j], i);

                const auto squaredRightDistance = ComputeSquaredDistance(m_gridPoints[j][i], m_gridPoints[j][firstRightIndex], m_splines->m_projection);

                if (squaredRightDistance < squaredDistanceTolerance)
                {
                    continue;
                }

                // Detect triangular cell
                if (!m_gridPoints[j + 1][i].IsValid())
                {
                    continue;
                }

                const auto squaredLeftDistance = ComputeSquaredDistance(m_gridPoints[j][firstLeftIndex], m_gridPoints[j][i], m_splines->m_projection);
                if (squaredLeftDistance < squaredDistanceTolerance)
                {
                    firstLeftIndex = i;
                }

                if (m_gridPoints[j + 1][firstRightIndex].IsValid())
                {
                    const auto squaredCurrentDistance = ComputeSquaredDistance(m_gridPoints[j + 1][i], m_gridPoints[j + 1][firstRightIndex], m_splines->m_projection);
                    const auto currentCosPhi = NormalizedInnerProductTwoSegments(
                        m_gridPoints[j + 1][i],
                        m_gridPoints[j][i],
                        m_gridPoints[j + 1][i],
                        m_gridPoints[j][firstRightIndex],
                        m_splines->m_projection);
                    if (squaredCurrentDistance < squaredDistanceTolerance && currentCosPhi > maxCosine)
                    {

                        //determine persistent node
                        const auto leftCosPhi = NormalizedInnerProductTwoSegments(
                            m_gridPoints[j - 1][i],
                            m_gridPoints[j][i],
                            m_gridPoints[j][i],
                            m_gridPoints[j + 1][i],
                            m_splines->m_projection);

                        const auto rightCosPhi = NormalizedInnerProductTwoSegments(
                            m_gridPoints[j - 1][firstRightIndex],
                            m_gridPoints[j][firstRightIndex],
                            m_gridPoints[j][firstRightIndex],
                            m_gridPoints[j + 1][firstRightIndex],
                            m_splines->m_projection);

                        const auto [secondLeftIndex, secondRightIndex] = GetNeighbours(m_gridPoints[j], firstRightIndex);

                        if ((secondRightIndex == firstRightIndex || leftCosPhi - rightCosPhi < -cosineTolerance) && firstLeftIndex != i)
                        {
                            //move left node
                            for (auto k = i; k <= firstRightIndex - 1; ++k)
                            {
                                m_gridPoints[j][k] = m_gridPoints[j][firstRightIndex];
                            }
                            numChanged++;
                        }
                        else if ((firstLeftIndex == i || rightCosPhi - leftCosPhi < -cosineTolerance) && secondRightIndex != firstRightIndex)
                        {
                            //move right node
                            for (auto k = firstRightIndex; k <= secondRightIndex - 1; ++k)
                            {
                                m_gridPoints[j][k] = m_gridPoints[j][i];
                            }
                            numChanged++;
                        }
                        else
                        {
                            //move both nodes
                            const Point middle = (m_gridPoints[j][i] + m_gridPoints[j][firstRightIndex]) * 0.5;
                            for (auto k = i; k <= firstRightIndex - 1; ++k)
                            {
                                m_gridPoints[j][k] = middle;
                            }
                            for (auto k = firstRightIndex; k <= secondRightIndex - 1; ++k)
                            {
                                m_gridPoints[j][k] = middle;
                            }
                            numChanged++;
                        }
                    }
                }
            }

            if (numChanged == 0)
            {
                break;
            }
        }
    }
}

void CurvilinearGridFromSplines::Initialize()
{
    // with only 1 spline, can't continue
    if (m_splines->GetNumSplines() <= 1)
    {
        throw std::invalid_argument("CurvilinearGridFromSplines::Initialize: Not enough splines to create a curvilinear grid.");
    }

    // compute properties
    ComputeSplineProperties(false);

    // get the properties of the center splines
    MakeAllGridLines();

    // Store original number of splines
    std::vector<Point> newCrossSpline(2);
    m_numOriginalSplines = m_splines->GetNumSplines();
    for (auto s = 0; s < m_numOriginalSplines; ++s)
    {
        // mirror only center splines
        if (m_type[s] != SplineTypes::central)
        {
            continue;
        }

        // construct the cross splines through the edges, along m discretization
        for (auto i = m_leftGridLineIndex[s]; i < m_leftGridLineIndex[s] + m_numMSplines[s]; ++i)
        {
            const auto normal = NormalVectorOutside(m_gridLine[i], m_gridLine[i + 1], m_splines->m_projection);

            const double xMiddle = (m_gridLine[i].x + m_gridLine[i + 1].x) * 0.5;
            const double yMiddle = (m_gridLine[i].y + m_gridLine[i + 1].y) * 0.5;

            double xs1 = doubleMissingValue;
            double xs2 = doubleMissingValue;
            double ys1 = doubleMissingValue;
            double ys2 = doubleMissingValue;

            if (m_splines->m_projection == Projection::cartesian)
            {
                xs1 = xMiddle + 2.0 * m_maximumGridHeights[s] * -normal.x;
                xs2 = xMiddle + 2.0 * m_maximumGridHeights[s] * normal.x;
                ys1 = yMiddle + 2.0 * m_maximumGridHeights[s] * -normal.y;
                ys2 = yMiddle + 2.0 * m_maximumGridHeights[s] * normal.y;
            }
            if (m_splines->m_projection == Projection::spherical)
            {
                const double factor = 1.0 / (earth_radius * degrad_hp);
                xs1 = xMiddle + 2.0 * m_maximumGridHeights[s] * -normal.x * factor;
                xs2 = xMiddle + 2.0 * m_maximumGridHeights[s] * normal.x * factor;
                ys1 = yMiddle + 2.0 * m_maximumGridHeights[s] * -normal.y * factor;
                ys2 = yMiddle + 2.0 * m_maximumGridHeights[s] * normal.y * factor;
            }

            newCrossSpline[0] = {xs1, ys1};
            newCrossSpline[1] = {xs2, ys2};
            m_splines->AddSpline(newCrossSpline, 0, newCrossSpline.size());
            // flag the cross spline as artificially added
            m_type.emplace_back(SplineTypes::artificial);
        }
    }

    // Backup original spline properties
    m_leftGridLineIndexOriginal.resize(m_numOriginalSplines);
    m_rightGridLineIndexOriginal.resize(m_numOriginalSplines);
    m_mfacOriginal.resize(m_numOriginalSplines);
    m_maximumGridHeightsOriginal.resize(m_numOriginalSplines);
    m_originalTypes.resize(m_numOriginalSplines);
    for (auto s = 0; s < m_numOriginalSplines; ++s)
    {
        m_leftGridLineIndexOriginal[s] = m_leftGridLineIndex[s];
        m_rightGridLineIndexOriginal[s] = m_rightGridLineIndex[s];
        m_mfacOriginal[s] = m_numMSplines[s];
        m_maximumGridHeightsOriginal[s] = m_maximumGridHeights[s];
        m_originalTypes[s] = m_type[s];
    }

    // compute spline properties with artificial splines added
    ComputeSplineProperties(true);

    // artificial cross spline: remove the last part of the sub-intervals (since it makes no sense,
    // as the artificial cross spline has an arbitrary, but sufficiently large, length)
    for (auto s = 0; s < m_numOriginalSplines; ++s)
    {
        // Remove the last part of the sub-intervals
        if (m_type[s] != SplineTypes::central)
        {
            continue;
        }

        // For number of intersecting splines
        for (auto i = 0; i < m_numCrossingSplines[s]; ++i)
        {
            const auto crossingSplineIndex = m_crossingSplinesIndices[s][i];
            if (m_type[crossingSplineIndex] == SplineTypes::artificial)
            {
                m_numCrossSplineLeftHeights[s][i] = m_numCrossSplineLeftHeights[s][i] - 1;
                m_numCrossSplineRightHeights[s][i] = m_numCrossSplineRightHeights[s][i] - 1;
            }
        }
    }

    // Compute edge velocities
    ComputeEdgeVelocities();

    // Increase curvilinear grid
    const auto numGridLayers = m_curvilinearParameters.n_refinement + 1;
    // The layer by coordinate to grow
    ResizeAndFill2DVector(m_gridPoints, numGridLayers + 1, m_numM + 1);
    m_validFrontNodes.resize(m_numM, 1);

    // Copy the first m point in m_gridPoints
    for (auto n = 0; n < m_numM; ++n)
    {
        m_gridPoints[0][n] = m_gridLine[n];
        if (!m_gridLine[n].IsValid())
        {
            m_validFrontNodes[n] = 0;
        }
        size_t sumLeft = 0;
        size_t sumRight = 0;
        const auto leftColumn = n == 0 ? 0 : n - 1;
        const auto rightColumn = n <= m_numM - 2 ? n : m_numM - 2;

        for (const auto& numFaces : m_numPerpendicularFacesOnSubintervalAndEdge)
        {
            sumLeft += numFaces[leftColumn];
            sumRight += numFaces[rightColumn];
        }

        if (sumLeft == 0 && sumRight == 0)
        {
            m_validFrontNodes[n] = 0;
        }
    }

    //compute maximum mesh width and get dtolLR in the proper dimension
    double squaredMaximumGridWidth = 0.0;
    for (auto i = 0; i < m_gridPoints[0].size() - 1; i++)
    {
        if (!m_gridPoints[0][i].IsValid() || !m_gridPoints[0][i + 1].IsValid())
        {
            continue;
        }
        squaredMaximumGridWidth = std::max(squaredMaximumGridWidth, ComputeSquaredDistance(m_gridPoints[0][i], m_gridPoints[0][i + 1], m_splines->m_projection));
    }
    m_onTopOfEachOtherSquaredTolerance = m_onTopOfEachOtherSquaredTolerance * squaredMaximumGridWidth;

    m_subLayerGridPoints.resize(m_numPerpendicularFacesOnSubintervalAndEdge.size());
}

void CurvilinearGridFromSplines::Iterate(size_t layer)
{
    GrowLayer(layer);

    for (auto j = 0; j < m_subLayerGridPoints.size(); ++j)
    {
        m_subLayerGridPoints[j] = m_numPerpendicularFacesOnSubintervalAndEdge[j][0];
    }

    auto results = ComputeGridLayerAndSubLayer(layer);
    auto subLayerRightIndex = std::get<1>(results);

    for (auto i = 0; i < m_numM; i++)
    {
        const auto subLayerLeftIndex = subLayerRightIndex;
        const auto minRight = std::min(static_cast<size_t>(i), m_numPerpendicularFacesOnSubintervalAndEdge[0].size() - static_cast<size_t>(1));
        for (auto j = 0; j < m_subLayerGridPoints.size(); ++j)
        {
            m_subLayerGridPoints[j] = m_numPerpendicularFacesOnSubintervalAndEdge[j][minRight];
        }

        results = ComputeGridLayerAndSubLayer(layer);
        const auto gridLayer = std::get<0>(results);
        subLayerRightIndex = std::get<1>(results);

        if (subLayerRightIndex != sizetMissingValue && i < m_numM - 1 && gridLayer != sizetMissingValue)
        {
            m_edgeVelocities[i] = m_growFactorOnSubintervalAndEdge[subLayerRightIndex][i] * m_edgeVelocities[i];
        }

        if (subLayerLeftIndex == sizetMissingValue && subLayerRightIndex == sizetMissingValue)
        {
            m_validFrontNodes[i] = sizetMissingValue;
        }
    }

    assert(m_timeStep > 1e-8 && "time step is smaller than 1e-8!");
}

CurvilinearGrid CurvilinearGridFromSplines::ComputeCurvilinearGridFromGridPoints()
{
    std::vector<std::vector<size_t>> mIndicesOtherSide(1, std::vector<size_t>(2));
    std::vector<std::vector<size_t>> nIndicesThisSide(1, std::vector<size_t>(2));
    std::vector<std::vector<Point>> gridPointsNDirection(m_gridPoints[0].size(),
                                                         std::vector<Point>(m_gridPoints.size(), {doubleMissingValue, doubleMissingValue}));
    std::vector<std::vector<Point>> curvilinearMeshPoints;
    const double squaredDistanceTolerance = 1e-12;

    // get the grid sizes in j-direction
    for (auto i = 0; i < m_gridPoints[0].size(); i++)
    {
        for (auto j = 0; j < m_gridPoints.size(); j++)
        {
            gridPointsNDirection[i][j] = m_gridPoints[j][i];
        }
    }

    size_t startIndex = 0;
    size_t startGridLine = 0;
    while (startIndex < m_gridPoints[0].size())
    {
        auto mIndicesThisSide = FindIndices(m_gridPoints[0], startIndex, m_numM, doubleMissingValue);

        mIndicesOtherSide[0][0] = mIndicesThisSide[0][1] + 2;
        mIndicesOtherSide[0][1] = mIndicesOtherSide[0][0] + (mIndicesThisSide[0][1] - mIndicesThisSide[0][0]);
        bool isConnected = true;

        size_t minN = m_curvilinearParameters.n_refinement;
        size_t maxN = 0;
        size_t minNOther = m_curvilinearParameters.n_refinement;
        size_t maxNOther = 0;
        //check if this part is connected to another part
        for (auto i = mIndicesThisSide[0][0]; i < mIndicesThisSide[0][1] + 1; ++i)
        {
            nIndicesThisSide = FindIndices(gridPointsNDirection[i], 0, gridPointsNDirection[i].size(), doubleMissingValue);
            minN = std::min(minN, nIndicesThisSide[0][0]);
            maxN = std::max(maxN, nIndicesThisSide[0][1]);

            const size_t mOther = mIndicesThisSide[0][1] + 2 + (mIndicesThisSide[0][1] - i);

            if (mOther > m_numM - 1)
            {
                // no more grid available
                isConnected = false;
            }
            else
            {
                const double squaredDistance = ComputeSquaredDistance(m_gridPoints[0][i], m_gridPoints[0][mOther], m_splines->m_projection);
                if (squaredDistance > squaredDistanceTolerance)
                {
                    isConnected = false;
                }
                else
                {
                    const auto nIndicesOtherSide = FindIndices(gridPointsNDirection[mOther], 0, gridPointsNDirection[mOther].size(), doubleMissingValue);
                    minNOther = std::min(minNOther, nIndicesOtherSide[0][0]);
                    maxNOther = std::max(maxNOther, nIndicesOtherSide[0][1]);
                }
            }
        }

        const auto endGridlineIndex = startGridLine + mIndicesThisSide[0][1] - mIndicesThisSide[0][0];
        if (isConnected)
        {
            startIndex = mIndicesOtherSide[0][1] + 2;
        }
        else
        {
            maxNOther = 1;
            startIndex = mIndicesThisSide[0][1] + 2;
        }

        // increment points
        curvilinearMeshPoints.resize(endGridlineIndex + 1);
        const auto NSize = std::max(curvilinearMeshPoints[0].size(), maxN + maxNOther + 1);
        for (auto& element : curvilinearMeshPoints)
        {
            element.resize(NSize);
        }

        // fill first part
        size_t columnIncrement = 0;
        for (auto i = startGridLine; i < endGridlineIndex + 1; ++i)
        {
            for (auto j = 0; j < maxN + 1; ++j)
            {
                curvilinearMeshPoints[i][j + maxNOther] = m_gridPoints[j][mIndicesThisSide[0][0] + columnIncrement];
            }
            columnIncrement++;
        }

        columnIncrement = 0;
        for (auto i = startGridLine; i < endGridlineIndex + 1; ++i)
        {
            for (auto j = 0; j < maxNOther + 1; ++j)
            {
                curvilinearMeshPoints[i][maxNOther - j] = m_gridPoints[j][mIndicesOtherSide[0][1] - columnIncrement];
            }
            columnIncrement++;
        }

        startGridLine = endGridlineIndex + 2;
    }

    return CurvilinearGrid(std::move(curvilinearMeshPoints), m_splines->m_projection);
}

std::tuple<size_t, size_t>
CurvilinearGridFromSplines::ComputeGridLayerAndSubLayer(size_t layerIndex)
{

    if (layerIndex == 0)
    {
        return {sizetMissingValue, sizetMissingValue};
    }

    size_t gridLayer = layerIndex - 1;
    auto sum = std::accumulate(m_subLayerGridPoints.begin(), m_subLayerGridPoints.end(), size_t(0));

    size_t subLayerIndex;
    if (layerIndex >= sum)
    {
        subLayerIndex = sizetMissingValue;
    }
    else
    {
        subLayerIndex = 0;
        sum = m_subLayerGridPoints[0] + 1;
        while (sum <= layerIndex && subLayerIndex < m_maxNumCenterSplineHeights)
        {
            subLayerIndex = subLayerIndex + 1;
            sum += m_subLayerGridPoints[subLayerIndex];
        }
        gridLayer = layerIndex - sum + m_subLayerGridPoints[subLayerIndex];
    }

    return {gridLayer, subLayerIndex};
}

void CurvilinearGridFromSplines::GrowLayer(size_t layerIndex)
{
    auto velocityVectorAtGridPoints = ComputeVelocitiesAtGridPoints(layerIndex - 1);

    std::vector<Point> activeLayerPoints(m_gridPoints[layerIndex - 1]);
    for (auto m = 0; m < velocityVectorAtGridPoints.size(); ++m)
    {
        if (!velocityVectorAtGridPoints[m].IsValid())
        {
            m_gridPoints[layerIndex - 1][m] = {doubleMissingValue, doubleMissingValue};
            activeLayerPoints[m] = {doubleMissingValue, doubleMissingValue};
        }
    }

    auto frontVelocities = CopyVelocitiesToFront(layerIndex - 1, velocityVectorAtGridPoints);

    double totalTimeStep = 0.0;
    std::vector<Point> gridLine(m_gridPoints[layerIndex - 1]);
    double localTimeStep = 0.0;
    double otherTimeStep = std::numeric_limits<double>::max();
    const auto numGridPoints = m_gridPoints.size() * m_gridPoints[0].size();
    std::vector<size_t> newValidFrontNodes(numGridPoints);

    while (totalTimeStep < m_timeStep)
    {
        // Copy old front velocities
        newValidFrontNodes = m_validFrontNodes;

        for (auto i = 0; i < m_validFrontNodes.size(); ++i)
        {
            if (m_validFrontNodes[i] == sizetMissingValue)
            {
                activeLayerPoints[i] = {doubleMissingValue, doubleMissingValue};
            }
        }

        const auto maximumGridLayerGrowTime = ComputeMaximumEdgeGrowTime(activeLayerPoints, velocityVectorAtGridPoints);
        localTimeStep = std::min(m_timeStep - totalTimeStep, *std::min_element(maximumGridLayerGrowTime.begin(), maximumGridLayerGrowTime.end()));

        if (m_splinesToCurvilinearParameters.check_front_collisions)
        {
            // TODO: implement front collisions
            otherTimeStep = 0;
        }
        localTimeStep = std::min(localTimeStep, otherTimeStep);

        // remove isolated points at the start end end of the masl
        if (newValidFrontNodes[0] == 1 && newValidFrontNodes[1] == 0)
        {
            newValidFrontNodes[0] = 0;
        }

        if (newValidFrontNodes[m_numM - 1] == 1 && newValidFrontNodes[m_numM - 2] == 0)
        {
            newValidFrontNodes[m_numM - 1] = 0;
        }

        for (auto i = 0; i < newValidFrontNodes.size() - 2; ++i)
        {
            if (newValidFrontNodes[i + 1] == 1 && newValidFrontNodes[i] == 0 && newValidFrontNodes[i + 2] == 0)
            {
                newValidFrontNodes[i + 1] = 0;
            }
        }

        m_validFrontNodes = newValidFrontNodes;

        for (auto i = 0; i < velocityVectorAtGridPoints.size(); ++i)
        {
            if (m_validFrontNodes[i] == 1 && velocityVectorAtGridPoints[i].IsValid())
            {
                activeLayerPoints[i].x = activeLayerPoints[i].x + localTimeStep * velocityVectorAtGridPoints[i].x;
                activeLayerPoints[i].y = activeLayerPoints[i].y + localTimeStep * velocityVectorAtGridPoints[i].y;
            }
            else
            {
                activeLayerPoints[i].x = doubleMissingValue;
                activeLayerPoints[i].y = doubleMissingValue;
            }
        }

        // update the grid points
        m_gridPoints[layerIndex] = activeLayerPoints;

        // update the time step
        totalTimeStep += localTimeStep;

        if (totalTimeStep < m_timeStep)
        {
            velocityVectorAtGridPoints = ComputeVelocitiesAtGridPoints(layerIndex);

            for (auto i = 0; i < m_numM; ++i)
            {
                // Disable points that have no valid normal vector
                // Remove stationary points
                if (!frontVelocities[i].IsValid() || m_validFrontNodes[i] == 0)
                {
                    activeLayerPoints[i] = {doubleMissingValue, doubleMissingValue};
                }
            }
            frontVelocities = CopyVelocitiesToFront(layerIndex - 1, velocityVectorAtGridPoints);
        }
    }

    auto [gridPointsIndices, frontGridPoints, numFrontPoints] = FindFront();
    if (layerIndex >= 2)
    {
        for (auto i = 1; i < m_numM - 1; ++i)
        {

            if (!activeLayerPoints[i].IsValid())
            {
                continue;
            }
            const auto cosphi = NormalizedInnerProductTwoSegments(m_gridPoints[layerIndex - 2][i],
                                                                  m_gridPoints[layerIndex - 1][i],
                                                                  m_gridPoints[layerIndex - 1][i],
                                                                  activeLayerPoints[i],
                                                                  m_splines->m_projection);
            if (cosphi < -0.5)
            {
                const auto [currentLeftIndex, currentRightIndex] = GetNeighbours(frontGridPoints, i);
                for (auto j = currentLeftIndex + 1; j < currentRightIndex; ++j)
                {
                    newValidFrontNodes[j] = 0;
                    m_gridPoints[layerIndex - 1][j] = {doubleMissingValue, doubleMissingValue};
                }
            }
        }
    }

    m_validFrontNodes = newValidFrontNodes;
}

std::vector<double> CurvilinearGridFromSplines::ComputeMaximumEdgeGrowTime(const std::vector<Point>& coordinates,
                                                                           const std::vector<Point>& velocities) const
{
    std::vector<double> maximumGridLayerGrowTime(m_validFrontNodes.size(), std::numeric_limits<double>::max());
    std::vector<double> edgeWidth(coordinates.size() - 1);
    std::vector<double> edgeIncrement(coordinates.size() - 1);
    const double minEdgeWidth = 1e-8;
    const double dt = 1.0;
    for (auto i = 0; i < coordinates.size() - 1; ++i)
    {
        if (!coordinates[i].IsValid() || !coordinates[i + 1].IsValid())
        {
            continue;
        }

        edgeWidth[i] = ComputeDistance(coordinates[i], coordinates[i + 1], m_splines->m_projection);

        if (edgeWidth[i] < minEdgeWidth)
        {
            continue;
        }

        Point firstPointIncremented(coordinates[i] + velocities[i] * dt);
        Point secondPointIncremented(coordinates[i + 1] + velocities[i + 1] * dt);

        edgeIncrement[i] = InnerProductTwoSegments(coordinates[i], coordinates[i + 1], firstPointIncremented, secondPointIncremented, m_splines->m_projection) / edgeWidth[i] - edgeWidth[i];
        edgeIncrement[i] = edgeIncrement[i] / dt;
    }

    for (auto i = 0; i < coordinates.size() - 1; ++i)
    {
        if (edgeIncrement[i] < 0.0)
        {
            maximumGridLayerGrowTime[i] = -edgeWidth[i] / edgeIncrement[i];
        }
    }

    return maximumGridLayerGrowTime;
}

std::vector<meshkernel::Point> CurvilinearGridFromSplines::CopyVelocitiesToFront(size_t layerIndex,
                                                                                 const std::vector<Point>& previousFrontVelocities)
{
    const auto numGridPoints = m_gridPoints.size() * m_gridPoints[0].size();
    std::vector<Point> velocities(numGridPoints, {0.0, 0.0});

    size_t numCornerNodes = 0;
    size_t p = 0;
    auto [gridPointsIndices, frontGridPoints, numFrontPoints] = FindFront();
    while (p < numFrontPoints)
    {
        if (gridPointsIndices[p][1] == layerIndex && m_validFrontNodes[gridPointsIndices[p][0]] == 1)
        {
            velocities[p] = previousFrontVelocities[gridPointsIndices[p][0]];
            if (!velocities[p].IsValid())
            {
                velocities[p] = {0.0, 0.0};
            }

            // Check for corner nodes
            const auto previous = p == 0 ? 0 : p - 1;
            const auto previousIndices = gridPointsIndices[previous];
            const auto next = std::min(p + 1, numFrontPoints);
            const auto nextIndices = gridPointsIndices[next];

            // Check corner nodes
            bool ll = previousIndices[0] == gridPointsIndices[p][0] - 1 &&
                      previousIndices[1] == gridPointsIndices[p][1] &&
                      m_validFrontNodes[previousIndices[0]] == sizetMissingValue;

            bool lr = nextIndices[0] == gridPointsIndices[p][0] + 1 &&
                      nextIndices[1] == gridPointsIndices[p][1] &&
                      m_validFrontNodes[nextIndices[0]] == sizetMissingValue;

            ll = ll || previousIndices[0] == gridPointsIndices[p][0] && previousIndices[1] < gridPointsIndices[p][1];
            lr = lr || nextIndices[0] == gridPointsIndices[p][0] && nextIndices[1] < gridPointsIndices[p][1];
            if (ll || lr)
            {
                numCornerNodes++;
                if (numFrontPoints + 1 > frontGridPoints.size())
                {
                    continue;
                }
                for (auto i = numFrontPoints; i >= p; --i)
                {
                    frontGridPoints[i + 1] = frontGridPoints[i];
                    velocities[i + 1] = velocities[i];
                    gridPointsIndices[i + 1] = gridPointsIndices[i];
                }
                numFrontPoints++;

                if (ll)
                {
                    velocities[p] = {0.0, 0.0};
                }
                else
                {
                    velocities[p + 1] = {0.0, 0.0};
                }
                p = p + 1;
            }
        }
        p = p + 1;
    }

    return velocities;
}

std::tuple<std::vector<std::vector<size_t>>, std::vector<meshkernel::Point>, size_t>
CurvilinearGridFromSplines::FindFront()
{
    const auto numGridPoints = m_gridPoints.size() * m_gridPoints[0].size();
    std::vector<std::vector<size_t>> gridPointsIndices(numGridPoints, std::vector<size_t>(2, sizetMissingValue));
    std::vector<Point> frontGridPoints(numGridPoints, {0.0, 0.0});
    size_t numFrontPoints;

    std::vector<int> frontPosition(m_gridPoints[0].size() - 2, static_cast<int>(m_gridPoints.size()));
    for (auto m = 0; m < frontPosition.size(); ++m)
    {
        for (auto n = 0; n < m_gridPoints.size(); ++n)
        {
            if (!m_gridPoints[n][m].IsValid() || !m_gridPoints[n][m + 1].IsValid())
            {
                frontPosition[m] = n - 1;
                break;
            }
        }
    }

    numFrontPoints = 0;
    // check for circular connectivity
    int previousFrontPosition = 0;
    const auto [leftNode, rightNode] = GetNeighbours(m_gridPoints[0], 0);
    if (leftNode == 0)
    {
        frontGridPoints[0] = m_gridPoints[0][0];
        // store front index
        gridPointsIndices[numFrontPoints][0] = 0;
        gridPointsIndices[numFrontPoints][1] = 0;
        numFrontPoints++;
    }
    else
    {
        previousFrontPosition = frontPosition[leftNode];
        frontGridPoints[numFrontPoints] = m_gridPoints[0][frontPosition[0]];
        gridPointsIndices[numFrontPoints][0] = frontPosition[0];
        gridPointsIndices[numFrontPoints][1] = 0;
        numFrontPoints++;
    }

    for (auto m = 0; m < m_gridPoints[0].size() - 2; ++m)
    {
        const auto currentFrontPosition = frontPosition[m];
        if (currentFrontPosition >= 0)
        {
            if (previousFrontPosition == -1)
            {
                frontGridPoints[numFrontPoints] = m_gridPoints[0][m];
                gridPointsIndices[numFrontPoints][0] = m;
                gridPointsIndices[numFrontPoints][1] = 0;
                numFrontPoints++;
            }
            for (auto i = previousFrontPosition + 1; i <= currentFrontPosition; ++i)
            {
                frontGridPoints[numFrontPoints] = m_gridPoints[i][m];
                gridPointsIndices[numFrontPoints][0] = m;
                gridPointsIndices[numFrontPoints][1] = i;
                numFrontPoints++;
            }
            for (auto i = previousFrontPosition; i > currentFrontPosition; --i)
            {
                frontGridPoints[numFrontPoints] = m_gridPoints[i][m];
                gridPointsIndices[numFrontPoints][0] = m;
                gridPointsIndices[numFrontPoints][1] = i;
                numFrontPoints++;
            }

            frontGridPoints[numFrontPoints] = m_gridPoints[currentFrontPosition][m + 1];
            gridPointsIndices[numFrontPoints][0] = m + 1;
            gridPointsIndices[numFrontPoints][1] = currentFrontPosition;
            numFrontPoints++;
        }
        else if (previousFrontPosition >= 0)
        {
            for (auto i = previousFrontPosition - 1; i >= 0; --i)
            {
                frontGridPoints[numFrontPoints] = m_gridPoints[i][m];
                gridPointsIndices[numFrontPoints][0] = m;
                gridPointsIndices[numFrontPoints][1] = i;
                numFrontPoints++;
            }

            frontGridPoints[numFrontPoints] = {doubleMissingValue, doubleMissingValue};
            gridPointsIndices[numFrontPoints][0] = m;
            gridPointsIndices[numFrontPoints][1] = sizetMissingValue;
            numFrontPoints++;
        }

        previousFrontPosition = currentFrontPosition;
    }

    // add last j-edge, check for circular connectivity
    const auto lastPoint = m_gridPoints[0].size() - 2;
    const auto [currentLeftIndex, currentRightIndex] = GetNeighbours(m_gridPoints[0], lastPoint);
    if (currentRightIndex == m_gridPoints[0].size() - 2)
    {
        for (auto i = previousFrontPosition; i >= 0; --i)
        {
            frontGridPoints[numFrontPoints] = m_gridPoints[i][lastPoint];
            gridPointsIndices[numFrontPoints][0] = lastPoint;
            gridPointsIndices[numFrontPoints][1] = i;
            numFrontPoints++;
        }
    }

    return {gridPointsIndices, frontGridPoints, numFrontPoints};
}

std::vector<meshkernel::Point>
CurvilinearGridFromSplines::ComputeVelocitiesAtGridPoints(size_t layerIndex)
{
    std::vector<Point> velocityVector(m_numM);
    std::fill(velocityVector.begin(), velocityVector.end(), Point{doubleMissingValue, doubleMissingValue});
    Point normalVectorLeft;
    Point normalVectorRight;
    const double cosTolerance = 1e-8;
    const double eps = 1e-10;
    for (auto m = 0; m < velocityVector.size(); ++m)
    {
        if (!m_gridPoints[layerIndex][m].IsValid())
        {
            continue;
        }

        const auto [currentLeftIndex, currentRightIndex] = GetNeighbours(m_gridPoints[layerIndex], m);
        const auto squaredLeftRightDistance = ComputeSquaredDistance(m_gridPoints[layerIndex][currentLeftIndex], m_gridPoints[layerIndex][currentRightIndex], m_splines->m_projection);

        if (squaredLeftRightDistance <= m_onTopOfEachOtherSquaredTolerance)
        {
            continue;
        }
        const auto squaredLeftDistance = ComputeSquaredDistance(m_gridPoints[layerIndex][currentLeftIndex], m_gridPoints[layerIndex][m], m_splines->m_projection);
        const auto squaredRightDistance = ComputeSquaredDistance(m_gridPoints[layerIndex][currentRightIndex], m_gridPoints[layerIndex][m], m_splines->m_projection);

        if (squaredLeftDistance <= m_onTopOfEachOtherSquaredTolerance || squaredRightDistance <= m_onTopOfEachOtherSquaredTolerance)
        {
            normalVectorLeft = NormalVectorOutside(m_gridPoints[layerIndex][currentRightIndex], m_gridPoints[layerIndex][currentLeftIndex], m_splines->m_projection);
            if (m_splines->m_projection == Projection::spherical)
            {
                normalVectorLeft.x = normalVectorLeft.x * std::cos(degrad_hp * 0.5 * (m_gridPoints[layerIndex][currentLeftIndex].y + m_gridPoints[layerIndex][currentRightIndex].y));
            }
            normalVectorRight = normalVectorLeft;
        }
        else
        {
            normalVectorLeft = NormalVectorOutside(m_gridPoints[layerIndex][m], m_gridPoints[layerIndex][currentLeftIndex], m_splines->m_projection);
            normalVectorRight = NormalVectorOutside(m_gridPoints[layerIndex][currentRightIndex], m_gridPoints[layerIndex][m], m_splines->m_projection);

            if (m_splines->m_projection == Projection::spherical)
            {
                normalVectorLeft.x = normalVectorLeft.x * std::cos(degrad_hp * 0.5 * (m_gridPoints[layerIndex][currentLeftIndex].y + m_gridPoints[layerIndex][m].y));
                normalVectorRight.x = normalVectorRight.x * std::cos(degrad_hp * 0.5 * (m_gridPoints[layerIndex][currentRightIndex].y + m_gridPoints[layerIndex][m].y));
            }
        }

        if (currentLeftIndex == velocityVector.size() - 1)
        {
            continue;
        }

        const auto cosphi = DotProduct(normalVectorLeft.x, normalVectorRight.x, normalVectorLeft.y, normalVectorRight.y);
        if (cosphi < -1.0 + cosTolerance)
        {
            continue;
        }
        const auto leftVelocity = normalVectorLeft * m_edgeVelocities[currentLeftIndex];
        const auto rightVelocity = normalVectorRight * m_edgeVelocities[currentRightIndex - 1];
        double rightLeftVelocityRatio = m_edgeVelocities[currentRightIndex - 1] / m_edgeVelocities[currentLeftIndex];

        if (rightLeftVelocityRatio - cosphi > eps && 1.0 / rightLeftVelocityRatio - cosphi > eps || cosphi <= cosTolerance)
        {
            velocityVector[m] = (leftVelocity * (1.0 - rightLeftVelocityRatio * cosphi) +
                                 rightVelocity * (1.0 - 1.0 / rightLeftVelocityRatio * cosphi)) /
                                (1.0 - cosphi * cosphi);
        }
        else if (cosphi - rightLeftVelocityRatio > eps)
        {
            velocityVector[m] = leftVelocity * rightLeftVelocityRatio / cosphi;
        }
        else
        {
            velocityVector[m] = rightVelocity * 1.0 / (rightLeftVelocityRatio * cosphi);
        }

        if (m_splines->m_projection == Projection::spherical)
        {
            velocityVector[m].x = velocityVector[m].x * one_over_earth_radius * raddeg_hp / std::cos(degrad_hp * m_gridPoints[layerIndex][m].y);
            velocityVector[m].y = velocityVector[m].y * one_over_earth_radius * raddeg_hp;
        }
    }

    return velocityVector;
}

std::tuple<size_t, size_t> CurvilinearGridFromSplines::GetNeighbours(const std::vector<Point>& gridPoints,
                                                                     size_t index) const
{

    if (gridPoints.empty())
    {
        return {sizetMissingValue, sizetMissingValue};
    }

    bool circularConnection = false;
    auto localLeftIndex = static_cast<int>(index);
    auto localRightIndex = static_cast<int>(index);

    // left
    while (ComputeSquaredDistance(gridPoints[localLeftIndex], gridPoints[index], m_splines->m_projection) < m_onTopOfEachOtherSquaredTolerance)
    {
        if (!circularConnection)
        {
            if (localLeftIndex == 0)
            {
                break;
            }
        }
        else if (localLeftIndex == 0)
        {
            localLeftIndex = static_cast<int>(gridPoints.size());
            circularConnection = false;
        }

        if (localLeftIndex == 0 || !gridPoints[localLeftIndex - 1].IsValid())
        {
            break;
        }
        localLeftIndex--;
    }

    // right
    while (ComputeSquaredDistance(gridPoints[localRightIndex], gridPoints[index], m_splines->m_projection) < m_onTopOfEachOtherSquaredTolerance)
    {
        if (!circularConnection)
        {
            if (localRightIndex == static_cast<int>(gridPoints.size()))
            {
                break;
            }
        }
        else if (localRightIndex == static_cast<int>(gridPoints.size()))
        {
            localRightIndex = -1;
            circularConnection = false;
        }

        if (localRightIndex == static_cast<int>(gridPoints.size()) || !gridPoints[localRightIndex + 1].IsValid())
        {
            break;
        }
        localRightIndex++;
    }

    return {static_cast<size_t>(localLeftIndex), static_cast<size_t>(localRightIndex)};
}

void CurvilinearGridFromSplines::ComputeEdgeVelocities()
{
    m_edgeVelocities.resize(m_numM - 1, doubleMissingValue);
    ResizeAndFill2DVector(m_growFactorOnSubintervalAndEdge, m_maxNumCenterSplineHeights, m_numM - 1, true, 1.0);
    ResizeAndFill2DVector(m_numPerpendicularFacesOnSubintervalAndEdge, m_maxNumCenterSplineHeights, m_numM - 1, true, static_cast<size_t>(0));

    ComputeGridHeights();

    std::fill(m_numPerpendicularFacesOnSubintervalAndEdge[0].begin(), m_numPerpendicularFacesOnSubintervalAndEdge[0].end(), 1);

    for (auto s = 0; s < m_splines->GetNumSplines(); s++)
    {

        if (m_type[s] != SplineTypes::central)
        {
            continue;
        }

        // Get true crossing splines heights
        auto numLeftHeights = m_maxNumCenterSplineHeights;
        auto numRightHeights = m_maxNumCenterSplineHeights;
        size_t numTrueCrossings = 0;
        for (auto i = 0; i < m_numCrossingSplines[s]; ++i)
        {
            if (m_type[m_crossingSplinesIndices[s][i]] != SplineTypes::crossing)
            {
                // true crossing splines only
                continue;
            }
            numTrueCrossings++;
            numLeftHeights = std::min(numLeftHeights, m_numCrossSplineLeftHeights[s][i]);
            numRightHeights = std::min(numRightHeights, m_numCrossSplineRightHeights[s][i]);
        }

        // no true cross splines: exponentially growing grid only
        if (numTrueCrossings == 0)
        {
            numLeftHeights = 0;
            numRightHeights = 0;
        }

        const auto startGridLineLeft = m_leftGridLineIndex[s];
        const auto endGridLineLeft = startGridLineLeft + m_numMSplines[s];
        const auto startGridLineRight = m_rightGridLineIndex[s];
        const auto endGridLineRight = startGridLineRight + m_numMSplines[s];

        double hh0LeftMaxRatio;
        double hh0RightMaxRatio;
        const size_t numIterations = 2;

        double maxHeight = std::numeric_limits<double>::lowest();
        for (const auto& e : m_gridHeights[0])
        {
            if (!IsEqual(e, doubleMissingValue) && e > maxHeight)
            {
                maxHeight = e;
            }
        }
        const double firstHeight = std::min(maxHeight, m_splinesToCurvilinearParameters.aspect_ratio * m_splinesToCurvilinearParameters.average_width);

        for (auto iter = 0; iter < numIterations; ++iter)
        {
            ComputeVelocitiesSubIntervals(s, startGridLineLeft, endGridLineLeft, numLeftHeights, numRightHeights, firstHeight,
                                          m_leftGridLineIndex, m_rightGridLineIndex, m_numPerpendicularFacesOnSubintervalAndEdge, m_edgeVelocities, hh0LeftMaxRatio);

            ComputeVelocitiesSubIntervals(s, startGridLineRight, endGridLineRight, numLeftHeights, numRightHeights, firstHeight,
                                          m_rightGridLineIndex, m_leftGridLineIndex, m_numPerpendicularFacesOnSubintervalAndEdge, m_edgeVelocities, hh0RightMaxRatio);
        }

        // re-evaluate if growing grid outside is needed

        if (numLeftHeights == 0 && numRightHeights <= 1 ||
            numRightHeights == 0 && numLeftHeights <= 1 ||
            numLeftHeights == 1 && numRightHeights == 1)
        {
            m_splinesToCurvilinearParameters.grow_grid_outside = 1;
        }

        // left part
        size_t numNLeftExponential = 0;
        if (m_splinesToCurvilinearParameters.grow_grid_outside == 1)
        {
            numNLeftExponential = std::min(ComputeNumberExponentialLayers(hh0LeftMaxRatio), static_cast<size_t>(m_curvilinearParameters.n_refinement));
        }
        for (auto i = startGridLineLeft; i < endGridLineLeft; ++i)
        {
            m_numPerpendicularFacesOnSubintervalAndEdge[1][i] = numNLeftExponential;
        }

        // right part
        size_t numNRightExponential = 0;
        if (m_splinesToCurvilinearParameters.grow_grid_outside == 1)
        {
            numNRightExponential = std::min(ComputeNumberExponentialLayers(hh0RightMaxRatio), static_cast<size_t>(m_curvilinearParameters.n_refinement));
        }
        for (auto i = startGridLineRight; i < endGridLineRight; ++i)
        {
            m_numPerpendicularFacesOnSubintervalAndEdge[1][i] = numNRightExponential;
        }
    }

    // compute local grow factors
    for (size_t s = 0; s < m_splines->GetNumSplines(); s++)
    {
        if (m_numMSplines[s] < 1)
        {
            continue;
        }

        for (auto i = m_leftGridLineIndex[s]; i < m_rightGridLineIndex[s] + m_numMSplines[s]; ++i)
        {
            if (!m_gridLine[i].IsValid() || !m_gridLine[i + 1].IsValid() || m_numPerpendicularFacesOnSubintervalAndEdge[1][i] < 1)
            {
                continue;
            }
            m_growFactorOnSubintervalAndEdge[1][i] = ComputeGrowFactor(i);
        }
    }
}

double CurvilinearGridFromSplines::ComputeGrowFactor(size_t splineIndex) const
{

    // eheight m_gridHeights
    double aspectRatioGrowFactor = 1.0;
    auto heightDifference = ComputeTotalExponentialHeight(aspectRatioGrowFactor, m_edgeVelocities[splineIndex], m_numPerpendicularFacesOnSubintervalAndEdge[1][splineIndex]) - m_gridHeights[1][splineIndex];

    const double deps = 0.01;
    double aspectRatioGrowFactorIncremented = 1.0 + deps;
    auto heightDifferenceIncremented = ComputeTotalExponentialHeight(aspectRatioGrowFactorIncremented, m_edgeVelocities[splineIndex], m_numPerpendicularFacesOnSubintervalAndEdge[1][splineIndex]) - m_gridHeights[1][splineIndex];

    const double tolerance = 1e-8;
    const size_t numIterations = 1000;
    const double relaxationFactor = 0.5;
    double oldAspectRatio;
    double oldHeightDifference = heightDifference;

    if (std::abs(heightDifferenceIncremented) > tolerance && std::abs(heightDifferenceIncremented - heightDifference) > tolerance)
    {
        for (auto i = 0; i < numIterations; ++i)
        {
            oldAspectRatio = aspectRatioGrowFactor;
            oldHeightDifference = heightDifference;

            aspectRatioGrowFactor = aspectRatioGrowFactorIncremented;
            heightDifference = heightDifferenceIncremented;

            aspectRatioGrowFactorIncremented = aspectRatioGrowFactor - relaxationFactor * heightDifference /
                                                                           (heightDifference - oldHeightDifference + 1e-16) *
                                                                           (aspectRatioGrowFactor - oldAspectRatio);

            heightDifferenceIncremented = ComputeTotalExponentialHeight(aspectRatioGrowFactorIncremented,
                                                                        m_edgeVelocities[splineIndex],
                                                                        m_numPerpendicularFacesOnSubintervalAndEdge[1][splineIndex]) -
                                          m_gridHeights[1][splineIndex];

            if (std::abs(oldHeightDifference) < tolerance)
            {
                break;
            }
        }
    }

    if (oldHeightDifference < tolerance)
    {
        return aspectRatioGrowFactorIncremented;
    }
    return 1.0;
}

double CurvilinearGridFromSplines::ComputeTotalExponentialHeight(double aspectRatio, double height, size_t numLayers) const
{

    if (std::abs(aspectRatio - 1.0) > 1e-8)
    {
        return (std::pow(aspectRatio, numLayers) - 1.0) / (aspectRatio - 1.0) * height;
    }
    return height * static_cast<double>(numLayers);
}

size_t CurvilinearGridFromSplines::ComputeNumberExponentialLayers(const double heightRatio) const
{
    if (m_splinesToCurvilinearParameters.aspect_ratio_grow_factor - 1.0 > 1e-8)
    {
        return size_t(std::floor(std::log((m_splinesToCurvilinearParameters.aspect_ratio_grow_factor - 1.0) * heightRatio + 1.0) /
                                 log(m_splinesToCurvilinearParameters.aspect_ratio_grow_factor)));
    }
    return size_t(std::floor(0.999 + heightRatio));
}

void CurvilinearGridFromSplines::ComputeVelocitiesSubIntervals(size_t s,
                                                               size_t startGridLineIndex,
                                                               size_t endGridLineIndex,
                                                               size_t numHeights,
                                                               size_t numOtherSideHeights,
                                                               const double firstHeight,
                                                               const std::vector<size_t>& gridLineIndex,
                                                               const std::vector<size_t>& otherGridLineIndex,
                                                               std::vector<std::vector<size_t>>& numPerpendicularFacesOnSubintervalAndEdge,
                                                               std::vector<double>& edgeVelocities,
                                                               double& hh0MaxRatio)
{
    hh0MaxRatio = 0.0;
    if (numHeights > 1 && numHeights == numOtherSideHeights || numHeights > numOtherSideHeights)
    {
        const auto maxHeight = *std::max_element(m_gridHeights[0].begin() + startGridLineIndex, m_gridHeights[0].begin() + endGridLineIndex);

        auto numNUniformPart = static_cast<size_t>(std::floor(maxHeight / firstHeight + 0.99999));
        numNUniformPart = std::min(numNUniformPart, m_maxNUniformPart);

        for (auto i = startGridLineIndex; i < endGridLineIndex; ++i)
        {
            numPerpendicularFacesOnSubintervalAndEdge[0][i] = numNUniformPart;
            edgeVelocities[i] = m_gridHeights[0][i] / static_cast<double>(numNUniformPart);
            hh0MaxRatio = std::max(hh0MaxRatio, m_gridHeights[1][i] / edgeVelocities[i]);
        }
    }
    else
    {
        // only one subinterval: no uniform part
        const size_t numNUniformPart = 0;
        for (auto i = startGridLineIndex; i < endGridLineIndex; ++i)
        {
            numPerpendicularFacesOnSubintervalAndEdge[0][i] = numNUniformPart;
            edgeVelocities[i] = firstHeight;

            //compare with other side of spline
            const auto otherSideIndex = otherGridLineIndex[s] + m_numMSplines[s] - (i - gridLineIndex[s] + 1);

            if (edgeVelocities[otherSideIndex] != doubleMissingValue)
            {
                if (numPerpendicularFacesOnSubintervalAndEdge[0][otherSideIndex] == 0)
                {
                    edgeVelocities[i] = std::max(edgeVelocities[i], edgeVelocities[otherSideIndex]);
                }
                else
                {
                    edgeVelocities[i] = edgeVelocities[otherSideIndex];
                }
            }

            for (auto j = 1; j < m_maxNumCenterSplineHeights; ++j)
            {
                m_gridHeights[j][i] = m_gridHeights[j - 1][i];
            }

            for (auto j = startGridLineIndex; j < endGridLineIndex; ++j)
            {

                hh0MaxRatio = std::max(hh0MaxRatio, m_gridHeights[1][j] / edgeVelocities[j]);
            }
        }
    }
}

void CurvilinearGridFromSplines::ComputeGridHeights()
{
    const auto numSplines = m_splines->GetNumSplines();
    ResizeAndFill2DVector(m_gridHeights, m_maxNumCenterSplineHeights, m_numM - 1, true, doubleMissingValue);

    std::vector<std::vector<double>> heightsLeft(m_maxNumCenterSplineHeights, std::vector<double>(m_curvilinearParameters.m_refinement, 0.0));
    std::vector<std::vector<double>> heightsRight(m_maxNumCenterSplineHeights, std::vector<double>(m_curvilinearParameters.m_refinement, 0.0));
    std::vector<double> edgesCenterPoints(m_numM, 0.0);
    std::vector<double> crossingSplinesDimensionalCoordinates(numSplines, 0.0);
    std::vector<double> localSplineDerivatives(numSplines, 0.0);
    std::vector<size_t> localValidSplineIndices(numSplines, 0);

    for (auto s = 0; s < numSplines; s++)
    {
        if (m_type[s] != SplineTypes::central)
        {
            continue;
        }

        const auto numM = m_numMSplines[s];

        // Get the minimum number of sub-intervals in the cross splines for this center spline
        const auto minNumLeftIntervals = *std::min_element(m_numCrossSplineLeftHeights[s].begin(), m_numCrossSplineLeftHeights[s].end());
        const auto minNumRightIntervals = *std::min_element(m_numCrossSplineRightHeights[s].begin(), m_numCrossSplineRightHeights[s].end());

        std::fill(heightsLeft[0].begin(), heightsLeft[0].begin() + numM, m_maximumGridHeights[s]);
        std::fill(heightsRight[0].begin(), heightsRight[0].begin() + numM, m_maximumGridHeights[s]);

        if (m_numCrossingSplines[s] == 1)
        {
            // only one crossing spline present:
            for (auto i = 0; i < minNumLeftIntervals; ++i)
            {
                std::fill(heightsLeft[i].begin(), heightsLeft[i].begin() + numM, m_crossSplineRightHeights[s][i][0]);
            }
            for (auto i = 0; i < minNumRightIntervals; ++i)
            {
                std::fill(heightsRight[i].begin(), heightsRight[i].begin() + numM, m_crossSplineLeftHeights[s][i][0]);
            }
        }
        else
        {
            const auto leftGridLineIndex = m_leftGridLineIndex[s];
            edgesCenterPoints[0] = m_splines->ComputeSplineLength(s, 0, m_gridLineDimensionalCoordinates[leftGridLineIndex]);
            for (auto i = 0; i < numM; ++i)
            {
                edgesCenterPoints[i + 1] = edgesCenterPoints[i] +
                                           m_splines->ComputeSplineLength(s,
                                                                          m_gridLineDimensionalCoordinates[leftGridLineIndex + i],
                                                                          m_gridLineDimensionalCoordinates[leftGridLineIndex + i + 1]);
            }

            // compute at edge center points
            for (auto i = 0; i < numM; ++i)
            {
                edgesCenterPoints[i] = 0.5 * (edgesCenterPoints[i] + edgesCenterPoints[i + 1]);
            }
            edgesCenterPoints[numM] = doubleMissingValue;

            //compute center spline path length of cross splines
            crossingSplinesDimensionalCoordinates[0] = m_splines->ComputeSplineLength(s, 0.0, m_crossSplineCoordinates[s][0]);
            for (auto i = 0; i < m_numCrossingSplines[s] - 1; ++i)
            {
                crossingSplinesDimensionalCoordinates[i + 1] = crossingSplinesDimensionalCoordinates[i] +
                                                               m_splines->ComputeSplineLength(s, m_crossSplineCoordinates[s][i], m_crossSplineCoordinates[s][i + 1]);
            }

            for (auto j = 0; j < m_maxNumCenterSplineHeights; ++j)
            {

                FindNearestCrossSplines(s,
                                        j,
                                        m_numCrossSplineLeftHeights[s],
                                        m_crossSplineLeftHeights[s],
                                        edgesCenterPoints,
                                        localValidSplineIndices,
                                        localSplineDerivatives,
                                        crossingSplinesDimensionalCoordinates,
                                        heightsLeft);

                FindNearestCrossSplines(s,
                                        j,
                                        m_numCrossSplineRightHeights[s],
                                        m_crossSplineRightHeights[s],
                                        edgesCenterPoints,
                                        localValidSplineIndices,
                                        localSplineDerivatives,
                                        crossingSplinesDimensionalCoordinates,
                                        heightsRight);
            }
        }

        // store grid height
        for (auto j = 0; j < m_maxNumCenterSplineHeights; ++j)
        {
            for (auto i = 0; i < m_numMSplines[s]; ++i)
            {
                m_gridHeights[j][m_leftGridLineIndex[s] + i] = heightsLeft[j][i];
                m_gridHeights[j][m_rightGridLineIndex[s] + m_numMSplines[s] - i - 1] = heightsRight[j][i];
            }
        }
    }
}

void CurvilinearGridFromSplines::FindNearestCrossSplines(size_t s,
                                                         size_t j,
                                                         const std::vector<size_t>& numHeightsLeft,
                                                         const std::vector<std::vector<double>>& crossSplineLeftHeights,
                                                         const std::vector<double>& edgesCenterPoints,
                                                         std::vector<size_t>& localValidSplineIndices,
                                                         std::vector<double>& localSplineDerivatives,
                                                         std::vector<double>& crossingSplinesDimensionalCoordinates,
                                                         std::vector<std::vector<double>>& heights)
{
    size_t numValid = 0;
    for (auto i = 0; i < m_numCrossingSplines[s]; ++i)
    {
        if (numHeightsLeft[i] != 0)
        {
            localValidSplineIndices[numValid] = i;
            numValid++;
        }
    }

    // no sub-heights to compute
    if (numValid == 0)
    {
        return;
    }

    const auto numM = m_numMSplines[s];
    std::vector<double> localCornerPoints(numValid);

    // TODO: strided memory access
    for (auto i = 0; i < numValid; ++i)
    {
        const auto index = localValidSplineIndices[i];
        localCornerPoints[i] = crossSplineLeftHeights[index][j];
    }

    localSplineDerivatives = Splines::SecondOrderDerivative(localCornerPoints, 0, localCornerPoints.size() - 1);

    crossingSplinesDimensionalCoordinates[0] = m_splines->ComputeSplineLength(s, 0.0, m_crossSplineCoordinates[s][0]);
    for (auto i = 0; i < numM; ++i)
    {
        size_t leftIndex = 0;
        double leftCoordinate = crossingSplinesDimensionalCoordinates[localValidSplineIndices[leftIndex]];
        auto rightIndex = std::min(size_t(1), numValid - 1);
        double rightCoordinate = crossingSplinesDimensionalCoordinates[localValidSplineIndices[rightIndex]];
        // Find two nearest cross splines
        while (rightCoordinate < edgesCenterPoints[i] && rightIndex < numValid)
        {
            leftIndex = rightIndex;
            leftCoordinate = rightCoordinate;
            rightIndex++;
            rightCoordinate = crossingSplinesDimensionalCoordinates[localValidSplineIndices[rightIndex]];
            if (rightIndex == numValid - 1)
            {
                break;
            }
        }

        double factor = 0.0;
        if (std::abs(rightCoordinate - leftCoordinate) > 1e-8)
        {
            factor = (edgesCenterPoints[i] - leftCoordinate) / (rightCoordinate - leftCoordinate);
        }

        factor = std::max(std::min(double(leftIndex + 1) + factor - 1.0, double(numValid - 1)), 0.0);

        heights[j][i] = ComputePointOnSplineAtAdimensionalDistance(localCornerPoints, localSplineDerivatives, factor);
    }
}

void CurvilinearGridFromSplines::GetSplineIntersections(size_t splineIndex)
{
    m_numCrossingSplines[splineIndex] = 0;
    const auto numSplines = m_splines->GetNumSplines();
    std::fill(m_crossingSplinesIndices[splineIndex].begin(), m_crossingSplinesIndices[splineIndex].end(), 0);
    std::fill(m_isLeftOriented[splineIndex].begin(), m_isLeftOriented[splineIndex].end(), true);
    std::fill(m_crossSplineCoordinates[splineIndex].begin(), m_crossSplineCoordinates[splineIndex].end(), std::numeric_limits<double>::max());
    std::fill(m_cosCrossingAngle[splineIndex].begin(), m_cosCrossingAngle[splineIndex].end(), doubleMissingValue);

    for (auto s = 0; s < numSplines; ++s)
    {
        // a crossing is a spline with 2 nodes and another with more than 2 nodes
        const auto numSplineNodesS = m_splines->m_splineNodes[s].size();
        const auto numSplineNodesI = m_splines->m_splineNodes[splineIndex].size();
        if (numSplineNodesS == 2 && numSplineNodesI == 2 ||
            numSplineNodesS > 2 && numSplineNodesI > 2)
        {
            continue;
        }

        double crossProductIntersection;
        Point intersectionPoint;
        double firstSplineRatio;
        double secondSplineRatio;
        bool crossing = m_splines->GetSplinesIntersection(splineIndex, s, crossProductIntersection, intersectionPoint, firstSplineRatio, secondSplineRatio);

        if (std::abs(crossProductIntersection) < m_splinesToCurvilinearParameters.min_cosine_crossing_angles)
        {
            crossing = false;
        }

        if (crossing)
        {
            m_numCrossingSplines[splineIndex]++;
            m_crossingSplinesIndices[splineIndex][s] = s;
            if (crossProductIntersection > 0.0)
            {
                m_isLeftOriented[splineIndex][s] = false;
            }
            m_crossSplineCoordinates[splineIndex][s] = firstSplineRatio;
            m_cosCrossingAngle[splineIndex][s] = crossProductIntersection;
        }
    }

    const auto sortedIndices = SortedIndices(m_crossSplineCoordinates[splineIndex]);
    m_crossSplineCoordinates[splineIndex] = ReorderVector(m_crossSplineCoordinates[splineIndex], sortedIndices);
    m_crossingSplinesIndices[splineIndex] = ReorderVector(m_crossingSplinesIndices[splineIndex], sortedIndices);
    m_isLeftOriented[splineIndex] = ReorderVector(m_isLeftOriented[splineIndex], sortedIndices);
}

void CurvilinearGridFromSplines::MakeAllGridLines()
{
    m_numM = 0;
    size_t numCenterSplines = 0;
    for (auto s = 0; s < m_splines->GetNumSplines(); ++s)
    {
        //center splines only
        if (m_type[s] != SplineTypes::central)
        {
            continue;
        }
        numCenterSplines += 1;
    }

    if (numCenterSplines == 0)
    {
        throw std::invalid_argument("CurvilinearGridFromSplines::MakeAllGridLines: There are no center splines.");
    }

    size_t gridLineIndex = 0;
    for (size_t s = 0; s < m_splines->GetNumSplines(); ++s)
    {
        //center splines only
        if (m_type[s] != SplineTypes::central)
        {
            continue;
        }

        // upper bound of m_gridLine size, with two sides of spline and two missing values added
        const auto sizeGridLine = gridLineIndex + 1 + 2 * (m_curvilinearParameters.m_refinement + 1) + 2;
        m_gridLine.resize(sizeGridLine);
        m_gridLineDimensionalCoordinates.resize(sizeGridLine);

        if (gridLineIndex > 0)
        {
            m_gridLine[gridLineIndex] = {doubleMissingValue, doubleMissingValue};
            m_gridLineDimensionalCoordinates[gridLineIndex] = doubleMissingValue;
            gridLineIndex++;
        }

        m_leftGridLineIndex[s] = gridLineIndex;

        const auto numM = MakeGridLine(s, gridLineIndex);

        gridLineIndex = gridLineIndex + numM + 1;
        m_gridLine[gridLineIndex] = Point{doubleMissingValue, doubleMissingValue};
        m_gridLineDimensionalCoordinates[gridLineIndex] = doubleMissingValue;
        gridLineIndex++;

        //add other side of gridline
        m_rightGridLineIndex[s] = gridLineIndex;
        auto rightIndex = m_rightGridLineIndex[s] - 1;
        for (auto j = m_rightGridLineIndex[s] - 1; j >= m_leftGridLineIndex[s] && j != static_cast<size_t>(0) - 1; --j)
        {
            m_gridLine[rightIndex] = m_gridLine[j];
            m_gridLineDimensionalCoordinates[rightIndex] = m_gridLineDimensionalCoordinates[j];
            ++rightIndex;
        }

        gridLineIndex = rightIndex;

        m_numMSplines[s] = numM;
        m_numM = gridLineIndex;
    }
}

size_t CurvilinearGridFromSplines::MakeGridLine(size_t splineIndex,
                                                size_t startingIndex)
{
    // first estimation of nodes along m
    auto numM = 1 + static_cast<size_t>(std::floor(m_splines->m_splinesLength[splineIndex] / m_splinesToCurvilinearParameters.average_width));
    numM = std::min(numM, static_cast<size_t>(m_curvilinearParameters.m_refinement));

    const auto endSplineAdimensionalCoordinate = static_cast<double>(m_splines->m_splineNodes[splineIndex].size()) - 1;
    const auto splineLength = m_splines->ComputeSplineLength(splineIndex,
                                                             0.0,
                                                             endSplineAdimensionalCoordinate,
                                                             10,
                                                             m_splinesToCurvilinearParameters.curvature_adapted_grid_spacing,
                                                             m_maximumGridHeights[splineIndex]);

    m_gridLine[startingIndex] = m_splines->m_splineNodes[splineIndex][0];

    auto currentMaxWidth = std::numeric_limits<double>::max();
    std::vector<double> distances(numM);
    while (currentMaxWidth > m_splinesToCurvilinearParameters.average_width)
    {
        currentMaxWidth = 0.0;
        for (auto n = 0; n < numM; ++n)
        {
            distances[n] = splineLength * (n + 1.0) / static_cast<double>(numM);
        }

        auto [points, adimensionalDistances] = m_splines->ComputePointOnSplineFromAdimensionalDistance(splineIndex,
                                                                                                       m_maximumGridHeights[splineIndex],
                                                                                                       m_splinesToCurvilinearParameters.curvature_adapted_grid_spacing,
                                                                                                       distances);

        for (auto n = 0; n < numM; ++n)
        {
            const auto index = startingIndex + n + 1;
            m_gridLineDimensionalCoordinates[index] = adimensionalDistances[n];
            m_gridLine[index] = points[n];
            currentMaxWidth = std::max(currentMaxWidth, ComputeDistance(m_gridLine[index - 1], m_gridLine[index], m_splines->m_projection));
        }

        // a gridline is computed
        if (currentMaxWidth < m_splinesToCurvilinearParameters.average_width ||
            numM == static_cast<size_t>(m_curvilinearParameters.m_refinement))
        {
            break;
        }

        // room for sub-division
        if (currentMaxWidth > m_splinesToCurvilinearParameters.average_width)
        {
            numM = std::min(std::max(static_cast<size_t>(m_curvilinearParameters.m_refinement / m_maximumGridHeights[splineIndex] * static_cast<double>(numM)),
                                     numM + static_cast<size_t>(1)),
                            static_cast<size_t>(m_curvilinearParameters.m_refinement));

            distances.resize(numM);
            adimensionalDistances.resize(numM);
            points.resize(numM);
        }
    }

    return numM;
}

void CurvilinearGridFromSplines::ComputeSplineProperties(const bool restoreOriginalProperties)
{
    AllocateSplinesProperties();

    for (size_t s = 0; s < m_splines->GetNumSplines(); ++s)
    {
        GetSplineIntersections(s);
    }

    // select all non-cross splines only
    for (size_t s = 0; s < m_splines->GetNumSplines(); ++s)
    {
        m_type[s] = SplineTypes::crossing;
        // if more than 2 nodes, the spline must be a central spline
        if (m_splines->m_splineNodes[s].size() > 2)
        {
            m_type[s] = SplineTypes::central;
        }
    }
    // check the cross splines. The center spline is the middle spline that crosses the cross spline
    for (size_t s = 0; s < m_splines->GetNumSplines(); ++s)
    {
        // only crossing splines with one or more central spline
        if (m_splines->m_splineNodes[s].size() != 2 || m_numCrossingSplines[s] < 1)
        {
            continue;
        }

        auto middleCrossingSpline = std::min(m_numCrossingSplines[s] / 2, m_numCrossingSplines[s]);
        auto crossingSplineIndex = m_crossingSplinesIndices[s][middleCrossingSpline];

        // if m_numIntersectingSplines[s] is even, check if the middle spline has already been assigned as a bounding spline
        if (m_type[crossingSplineIndex] != SplineTypes::central && 2 * crossingSplineIndex == m_numCrossingSplines[s])
        {
            middleCrossingSpline = std::min(middleCrossingSpline + 1, m_numCrossingSplines[s] - 1);
            crossingSplineIndex = m_crossingSplinesIndices[s][middleCrossingSpline];
        }

        if (m_type[crossingSplineIndex] == SplineTypes::central)
        {
            // associate bounding splines with the middle spline
            for (auto i = 0; i < middleCrossingSpline; ++i)
            {
                const auto index = m_crossingSplinesIndices[s][i];
                m_type[index] = SplineTypes::lateral; // lateral spline
                m_centralSplineIndex[index] = -static_cast<int>(crossingSplineIndex);
            }
            for (auto i = middleCrossingSpline + 1; i < m_numCrossingSplines[s]; ++i)
            {
                const auto index = m_crossingSplinesIndices[s][i];
                m_type[index] = SplineTypes::lateral; // lateral spline
                m_centralSplineIndex[index] = -static_cast<int>(crossingSplineIndex);
            }
        }
    }

    if (restoreOriginalProperties)
    {
        // restore original spline properties
        for (size_t s = 0; s < m_numOriginalSplines; ++s)
        {
            m_leftGridLineIndex[s] = m_leftGridLineIndexOriginal[s];
            m_rightGridLineIndex[s] = m_rightGridLineIndexOriginal[s];
            m_numMSplines[s] = m_mfacOriginal[s];
            m_maximumGridHeights[s] = m_maximumGridHeightsOriginal[s];
            m_type[s] = m_originalTypes[s];
        }

        //mark new splines as artificial cross splines
        for (size_t s = m_numOriginalSplines; s < m_splines->GetNumSplines(); ++s)
        {
            m_type[s] = SplineTypes::artificial;
        }
    }

    ComputeHeights();
}

void CurvilinearGridFromSplines::ComputeHeights()
{
    for (size_t i = 0; i < m_splines->GetNumSplines(); ++i)
    {
        // Heights should be computed only for center splines
        if (m_splines->m_splineNodes[i].size() <= 2)
        {
            continue;
        }
        for (auto j = 0; j < m_numCrossingSplines[i]; ++j)
        {
            ComputeSubHeights(i, j);
        }
    }

    // compute m_maximumGridHeight
    for (size_t s = 0; s < m_splines->GetNumSplines(); ++s)
    {
        if (m_numCrossingSplines[s] == 0)
        {
            m_maximumGridHeights[s] = m_splinesToCurvilinearParameters.aspect_ratio * m_splines->m_splinesLength[s];
            continue;
        }
        double maximumHeight = 0.0;
        for (auto c = 0; c < m_numCrossingSplines[s]; ++c)
        {
            double sumLeftHeights = 0.0;
            for (auto ss = 0; ss < m_numCrossSplineLeftHeights[s][c]; ++ss)
            {
                sumLeftHeights += m_crossSplineLeftHeights[s][c][ss];
            }
            double sumRightHeights = 0.0;
            for (auto ss = 0; ss < m_numCrossSplineRightHeights[s][c]; ++ss)
            {
                sumRightHeights += m_crossSplineRightHeights[s][c][ss];
            }
            maximumHeight = std::max(maximumHeight, std::max(sumLeftHeights, sumRightHeights));
        }

        m_maximumGridHeights[s] = maximumHeight;
    }
}

void CurvilinearGridFromSplines::ComputeSubHeights(size_t centerSplineIndex, size_t crossingSplineLocalIndex)
{
    // find center spline index
    size_t centerSplineLocalIndex = 0;
    const auto crossingSplineIndex = m_crossingSplinesIndices[centerSplineIndex][crossingSplineLocalIndex]; //js
    for (auto s = 0; s < m_numCrossingSplines[crossingSplineIndex]; ++s)
    {
        if (m_crossingSplinesIndices[crossingSplineIndex][s] == centerSplineIndex)
        {
            centerSplineLocalIndex = s;
            break;
        }
    }

    // right part
    size_t numSubIntervalsRight = 0;
    size_t rightCenterSplineIndex = centerSplineLocalIndex;
    size_t leftCenterSplineIndex;
    m_crossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex].resize(m_maxNumCenterSplineHeights, 0);
    for (auto s = centerSplineLocalIndex; s < m_numCrossingSplines[crossingSplineIndex] - 1; ++s)
    {
        if (numSubIntervalsRight >= m_maxNumCenterSplineHeights)
        {
            break;
        }
        if (m_centralSplineIndex[m_crossingSplinesIndices[crossingSplineIndex][s + 1]] != -static_cast<int>(centerSplineIndex))
        {
            continue;
        }
        leftCenterSplineIndex = rightCenterSplineIndex;
        rightCenterSplineIndex = s + 1;
        m_crossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex][numSubIntervalsRight] =
            m_splines->ComputeSplineLength(crossingSplineIndex,
                                           m_crossSplineCoordinates[crossingSplineIndex][leftCenterSplineIndex],
                                           m_crossSplineCoordinates[crossingSplineIndex][rightCenterSplineIndex]);
        numSubIntervalsRight++;
    }

    const auto numSplineNodes = m_splines->m_splineNodes[crossingSplineIndex].size();
    m_crossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex][numSubIntervalsRight] =
        m_splines->ComputeSplineLength(crossingSplineIndex,
                                       m_crossSplineCoordinates[crossingSplineIndex][rightCenterSplineIndex],
                                       static_cast<double>(numSplineNodes) - 1.0);
    numSubIntervalsRight++;
    std::fill(m_crossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex].begin() + numSubIntervalsRight,
              m_crossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex].end(), 0.0);

    m_numCrossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex] = numSubIntervalsRight;

    // left part
    size_t numSubIntervalsLeft = 0;
    leftCenterSplineIndex = centerSplineLocalIndex;
    m_crossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex].resize(m_maxNumCenterSplineHeights, 0);
    for (auto s = centerSplineLocalIndex; s >= 1; --s)
    {
        if (numSubIntervalsLeft >= m_maxNumCenterSplineHeights)
        {
            break;
        }
        if (m_centralSplineIndex[m_crossingSplinesIndices[crossingSplineIndex][s - 1]] != -static_cast<int>(centerSplineIndex))
        {
            continue;
        }
        rightCenterSplineIndex = leftCenterSplineIndex;
        leftCenterSplineIndex = s - 1;
        m_crossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex][numSubIntervalsLeft] =
            m_splines->ComputeSplineLength(crossingSplineIndex,
                                           m_crossSplineCoordinates[crossingSplineIndex][leftCenterSplineIndex],
                                           m_crossSplineCoordinates[crossingSplineIndex][rightCenterSplineIndex]);
        numSubIntervalsLeft++;
    }

    m_crossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex][numSubIntervalsLeft] =
        m_splines->ComputeSplineLength(crossingSplineIndex,
                                       0.0,
                                       m_crossSplineCoordinates[crossingSplineIndex][leftCenterSplineIndex]);

    numSubIntervalsLeft++;
    std::fill(m_crossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex].begin() + numSubIntervalsLeft,
              m_crossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex].end(), 0.0);

    m_numCrossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex] = numSubIntervalsLeft;

    // if not left oriented, swap
    if (!m_isLeftOriented[centerSplineIndex][crossingSplineLocalIndex])
    {
        m_numCrossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex] = numSubIntervalsRight;
        m_numCrossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex] = numSubIntervalsLeft;

        const std::vector<double> leftSubIntervalsTemp(m_crossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex]);
        m_crossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex] = m_crossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex];
        m_crossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex] = leftSubIntervalsTemp;
    }
}
