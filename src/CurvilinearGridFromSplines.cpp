//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2020.
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

#include <vector>
#include <algorithm>
#include <cassert>
#include "Operations.cpp"
#include "Entities.hpp"
#include "CurvilinearParametersNative.hpp"
#include "SplinesToCurvilinearParametersNative.hpp"
#include "Splines.hpp"
#include "CurvilinearGridFromSplines.hpp"
#include "CurvilinearGrid.hpp"

meshkernel::CurvilinearGridFromSplines::CurvilinearGridFromSplines(std::shared_ptr<Splines> splines,
                                                                   const meshkernelapi::CurvilinearParametersNative& curvilinearParametersNative,
                                                                   const meshkernelapi::SplinesToCurvilinearParametersNative& splinesToCurvilinearParametersNative) :
    m_splines(splines),
    m_curvilinearParametersNative(curvilinearParametersNative),
    m_splinesToCurvilinearParametersNative(splinesToCurvilinearParametersNative)
{
    m_onTopOfEachOtherSquaredTolerance = m_splinesToCurvilinearParametersNative.GridsOnTopOfEachOtherTolerance *
                                         m_splinesToCurvilinearParametersNative.GridsOnTopOfEachOtherTolerance;
};

/// to be called after all splines have been stored
bool meshkernel::CurvilinearGridFromSplines::AllocateSplinesProperties()
{
    const auto numSplines = m_splines->m_numSplines;
    m_type.resize(numSplines);

    m_centralSplineIndex.resize(numSplines);
    std::fill(m_centralSplineIndex.begin(), m_centralSplineIndex.end(), intMissingValue);

    m_numCrossingSplines.resize(numSplines, 0);
    std::fill(m_numCrossingSplines.begin(), m_numCrossingSplines.end(), 0);

    m_maximumGridHeights.resize(numSplines);
    std::fill(m_maximumGridHeights.begin(), m_maximumGridHeights.end(), doubleMissingValue);

    // multi-dimensional arrays
    m_crossingSplinesIndexses.resize(numSplines);
    m_isLeftOriented.resize(numSplines);
    m_crossSplineCoordinates.resize(numSplines);
    m_cosCrossingAngle.resize(numSplines);
    m_crossSplineLeftHeights.resize(numSplines);
    m_crossSplineRightHeights.resize(numSplines);
    m_numCrossSplineLeftHeights.resize(numSplines);
    m_numCrossSplineRightHeights.resize(numSplines);
    m_nfacL.resize(numSplines);
    m_nfacR.resize(numSplines);
    for (int s = 0; s < numSplines; ++s)
    {
        m_crossingSplinesIndexses[s].resize(numSplines);
        std::fill(m_crossingSplinesIndexses[s].begin(), m_crossingSplinesIndexses[s].end(), -1);

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

        m_nfacL[s].resize(numSplines);
        std::fill(m_nfacL[s].begin(), m_nfacL[s].end(), 0);

        m_nfacR[s].resize(numSplines);
        std::fill(m_nfacR[s].begin(), m_nfacR[s].end(), 0);
    }

    m_numMSplines.resize(numSplines);
    std::fill(m_numMSplines.begin(), m_numMSplines.end(), 0);
    m_leftGridLineIndex.resize(numSplines);
    std::fill(m_leftGridLineIndex.begin(), m_leftGridLineIndex.end(), intMissingValue);
    m_rightGridLineIndex.resize(numSplines);
    std::fill(m_rightGridLineIndex.begin(), m_rightGridLineIndex.end(), intMissingValue);

    return true;
}



bool meshkernel::CurvilinearGridFromSplines::Compute(CurvilinearGrid& curvilinearGrid)
{

    bool successful = Initialize();
    if (!successful)
    {
        return false;
    }

    // Grow grid, from the second layer
    for (int layer = 1; layer <= m_curvilinearParametersNative.NRefinement; ++layer)
    {
        successful = Iterate(layer);
        if (!successful)
        {
            break;
        }
    }

    bool removeSkinnyTriangles = m_splinesToCurvilinearParametersNative.RemoveSkinnyTriangles == 1 ? true : false;
    if (removeSkinnyTriangles)
    {
        RemoveSkinnyTriangles();
    }

    successful = ComputeCurvilinearGrid(curvilinearGrid);
    return successful;
}


bool meshkernel::CurvilinearGridFromSplines::RemoveSkinnyTriangles()
{
    int numMaxIterations = 10;
    auto numN = m_gridPoints.size() - 2;
    const double squaredDistanceTolerance = 1e-4;
    const double cosineTolerance = 1e-2;
    const double maxCosine = 0.93969;
    for (auto j = numN - 1; j >= 1; --j)
    {
        for (int iter = 0; iter < numMaxIterations; ++iter)
        {
            int numChanged = 0;
            
            int firstLeftIndex;
            int firstRightIndex = 0;
            int i = 0;

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

                GetNeighbours(m_gridPoints[j], i, firstLeftIndex, firstRightIndex);

                double squaredRightDistance = ComputeSquaredDistance(m_gridPoints[j][i], m_gridPoints[j][firstRightIndex], m_splines->m_projection);

                if (squaredRightDistance < squaredDistanceTolerance)
                {
                    continue;
                }

                // Detect triangular cell
                if (!m_gridPoints[j + 1][i].IsValid())
                {
                    continue;
                }

                double squaredLeftDistance = ComputeSquaredDistance(m_gridPoints[j][firstLeftIndex], m_gridPoints[j][i], m_splines->m_projection);
                if (squaredLeftDistance < squaredDistanceTolerance)
                {
                    firstLeftIndex = i;
                }

                if (m_gridPoints[j + 1][firstRightIndex].IsValid())
                {
                    double squaredCurrentDistance = ComputeSquaredDistance(m_gridPoints[j + 1][i], m_gridPoints[j + 1][firstRightIndex], m_splines->m_projection);
                    double currentCosPhi = NormalizedInnerProductTwoSegments(
                        m_gridPoints[j + 1][i],
                        m_gridPoints[j][i],
                        m_gridPoints[j + 1][i],
                        m_gridPoints[j][firstRightIndex],
                        m_splines->m_projection);
                    if (squaredCurrentDistance < squaredDistanceTolerance && currentCosPhi > maxCosine)
                    {

                        //determine persistent node
                        double leftCosPhi = NormalizedInnerProductTwoSegments(
                            m_gridPoints[j - 1][i],
                            m_gridPoints[j][i],
                            m_gridPoints[j][i],
                            m_gridPoints[j + 1][i],
                            m_splines->m_projection);

                        double rightCosPhi = NormalizedInnerProductTwoSegments(
                            m_gridPoints[j - 1][firstRightIndex],
                            m_gridPoints[j][firstRightIndex],
                            m_gridPoints[j][firstRightIndex],
                            m_gridPoints[j + 1][firstRightIndex],
                            m_splines->m_projection);


                        int secondLeftIndex;
                        int secondRightIndex;
                        GetNeighbours(m_gridPoints[j], firstRightIndex, secondLeftIndex, secondRightIndex);


                        if ((secondRightIndex == firstRightIndex || leftCosPhi - rightCosPhi < -cosineTolerance) && firstLeftIndex != i)
                        {
                            //move left node
                            for (int k = i; k <= firstRightIndex - 1; ++k)
                            {
                                m_gridPoints[j][k] = m_gridPoints[j][firstRightIndex];
                            }
                            numChanged++;
                        }
                        else if ((firstLeftIndex == i || rightCosPhi - leftCosPhi < -cosineTolerance) && secondRightIndex != firstRightIndex)
                        {
                            //move right node
                            for (int k = firstRightIndex; k <= secondRightIndex - 1; ++k)
                            {
                                m_gridPoints[j][k] = m_gridPoints[j][i];
                            }
                            numChanged++;
                        }
                        else
                        {
                            //move both nodes
                            Point middle = (m_gridPoints[j][i] + m_gridPoints[j][firstRightIndex]) * 0.5;
                            for (int k = i; k <= firstRightIndex - 1; ++k)
                            {
                                m_gridPoints[j][k] = middle;
                            }
                            for (int k = firstRightIndex; k <= secondRightIndex - 1; ++k)
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

    return true;
}

bool meshkernel::CurvilinearGridFromSplines::Initialize()
{
    // no splines
    if (m_splines->m_numSplines < 2)
    {
        return false;
    }

    // compute properties
    bool successful = ComputeSplineProperties(false);
    if (!successful)
    {
        return false;
    }

    // get the properties of the center splines
    m_numM = 0;
    successful = MakeAllGridLines(true);
    if (!successful)
    {
        return false;
    }

    // Store original number of splines
    std::vector<Point> newCrossSpline(2);
    m_numOriginalSplines = m_splines->m_numSplines;
    for (int s = 0; s < m_numOriginalSplines; ++s)
    {
        // mirrow only center splines
        if (m_type[s] != SplineTypes::central)
        {
            continue;
        }

        // construct the cross splines through the edges, along m discretization
        for (int i = m_leftGridLineIndex[s]; i < m_leftGridLineIndex[s] + m_numMSplines[s]; ++i)
        {
            Point normal{ doubleMissingValue, doubleMissingValue };
            NormalVectorOutside(m_gridLine[i], m_gridLine[i + 1], normal, m_splines->m_projection);

            double xMiddle = (m_gridLine[i].x + m_gridLine[i + 1].x) * 0.5;
            double yMiddle = (m_gridLine[i].y + m_gridLine[i + 1].y) * 0.5;
            
            double xs1;
            double xs2;
            double ys1;
            double ys2;

            if (m_splines->m_projection == Projections::cartesian)
            {
                xs1 = xMiddle + 2.0 * m_maximumGridHeights[s] * -normal.x;
                xs2 = xMiddle + 2.0 * m_maximumGridHeights[s] * normal.x;
                ys1 = yMiddle + 2.0 * m_maximumGridHeights[s] * -normal.y;
                ys2 = yMiddle + 2.0 * m_maximumGridHeights[s] * normal.y;
            }
            if (m_splines->m_projection == Projections::spherical)
            {
                const double factor = 1.0 / (earth_radius * degrad_hp);
                xs1 = xMiddle + 2.0 * m_maximumGridHeights[s] * -normal.x * factor;
                xs2 = xMiddle + 2.0 * m_maximumGridHeights[s] * normal.x * factor;
                ys1 = yMiddle + 2.0 * m_maximumGridHeights[s] * -normal.y * factor;
                ys2 = yMiddle + 2.0 * m_maximumGridHeights[s] * normal.y * factor;
            }

            newCrossSpline[0] = { xs1, ys1 };
            newCrossSpline[1] = { xs2, ys2 };
            m_splines->AddSpline(newCrossSpline, 0, newCrossSpline.size());
            // flag the cross spline as artificially added
            ResizeVectorIfNeeded(m_splines->m_numAllocatedSplines, m_type);
            m_type[m_splines->m_numSplines - 1] = SplineTypes::arficial;
        }
    }

    // Backup original spline properties
    m_leftGridLineIndexOriginal.resize(m_numOriginalSplines);
    m_rightGridLineIndexOriginal.resize(m_numOriginalSplines);
    m_mfacOriginal.resize(m_numOriginalSplines);
    m_maximumGridHeightsOriginal.resize(m_numOriginalSplines);
    m_originalTypes.resize(m_numOriginalSplines);
    for (int s = 0; s < m_numOriginalSplines; ++s)
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
    for (int s = 0; s < m_numOriginalSplines; ++s)
    {
        // Remove the last part of the sub-intervals
        if (m_type[s] != SplineTypes::central)
        {
            continue;
        }

        // For number of intersecting splines
        for (int i = 0; i < m_numCrossingSplines[s]; ++i)
        {
            int crossingSplineIndex = m_crossingSplinesIndexses[s][i];
            if (m_type[crossingSplineIndex] == SplineTypes::arficial)
            {
                m_numCrossSplineLeftHeights[s][i] = m_numCrossSplineLeftHeights[s][i] - 1;
                m_numCrossSplineRightHeights[s][i] = m_numCrossSplineRightHeights[s][i] - 1;
            }
        }
    }

    // Compute edge velocities
    m_edgeVelocities.resize(m_numM - 1, doubleMissingValue);
    m_growFactorOnSubintervalAndEdge.resize(m_maxNumCenterSplineHeights, std::vector<double>(m_numM - 1, 1.0));
    m_numPerpendicularFacesOnSubintervalAndEdge.resize(m_maxNumCenterSplineHeights, std::vector<int>(m_numM - 1, 0));
    successful = ComputeEdgeVelocities(m_edgeVelocities, m_growFactorOnSubintervalAndEdge, m_numPerpendicularFacesOnSubintervalAndEdge);
    if (!successful)
    {
        return false;
    }

    // Increase curvilinear grid
    const int numGridLayers = m_curvilinearParametersNative.NRefinement + 1;
    // The layer by coordinate to grow
    m_gridPoints.resize(numGridLayers + 1, std::vector<Point>(m_numM + 1, { doubleMissingValue, doubleMissingValue }));
    m_validFrontNodes.resize(m_numM, 1);

    // Copy the first n in m_gridPoints
    for (int n = 0; n < m_numM; ++n)
    {
        m_gridPoints[0][n] = m_gridLine[n];
        if (m_gridLine[n].x == doubleMissingValue)
        {
            m_validFrontNodes[n] = 0;
        }
        int sumLeft = 0;
        int sumRight = 0;
        auto leftColumn = std::max(n - 1, 0);
        auto rightColumn = std::min(n, int(m_numM) - 2);
        for (int j = 0; j < m_numPerpendicularFacesOnSubintervalAndEdge.size(); ++j)
        {
            sumLeft += m_numPerpendicularFacesOnSubintervalAndEdge[j][leftColumn];
            sumRight += m_numPerpendicularFacesOnSubintervalAndEdge[j][rightColumn];
        }
        if (sumLeft == 0 && sumRight == 0)
        {
            m_validFrontNodes[n] = 0;
        }
    }

    //compute maximum mesh width and get dtolLR in the proper dimension
    double squaredMaximumGridWidth = 0.0;
    for (int i = 0; i < m_gridPoints[0].size() - 1; i++)
    {
        if (!m_gridPoints[0][i].IsValid() || !m_gridPoints[0][i + 1].IsValid())
        {
            continue;
        }
        squaredMaximumGridWidth = std::max(squaredMaximumGridWidth, ComputeSquaredDistance(m_gridPoints[0][i], m_gridPoints[0][i + 1], m_splines->m_projection));
    }
    m_onTopOfEachOtherSquaredTolerance = m_onTopOfEachOtherSquaredTolerance * squaredMaximumGridWidth;


    m_subLayerGridPoints.resize(m_numPerpendicularFacesOnSubintervalAndEdge.size());

    return true;
}

bool meshkernel::CurvilinearGridFromSplines::Iterate(int layer)
{
    bool successful = GrowLayer(layer);
    if (!successful)
    {
        return false;
    }

    for (int j = 0; j < m_subLayerGridPoints.size(); ++j)
    {
        m_subLayerGridPoints[j] = m_numPerpendicularFacesOnSubintervalAndEdge[j][0];
    }

    int gridLayer;
    int subLayerRightIndex;

    successful = GetSubIntervalAndGridLayer(layer, gridLayer, subLayerRightIndex);
    if (!successful)
    {
        return false;
    }

    for (int i = 0; i < m_numM; i++)
    {
        int subLayerLeftIndex = subLayerRightIndex;
        int minRight = std::min(i, int(m_numPerpendicularFacesOnSubintervalAndEdge[0].size() - 1));
        for (int j = 0; j < m_subLayerGridPoints.size(); ++j)
        {
            m_subLayerGridPoints[j] = m_numPerpendicularFacesOnSubintervalAndEdge[j][minRight];
        }

        successful = GetSubIntervalAndGridLayer(layer, gridLayer, subLayerRightIndex);
        if (!successful)
        {
            return false;
        }

        if (subLayerRightIndex >= 0 && i < m_numM - 1 && gridLayer >= 0)
        {
            m_edgeVelocities[i] = m_growFactorOnSubintervalAndEdge[subLayerRightIndex][i] * m_edgeVelocities[i];
        }

        if (subLayerLeftIndex < 0 && subLayerRightIndex < 0)
        {
            m_validFrontNodes[i] = -1;
        }
    }

    if (m_timeStep < 1e-8)
    {
        return false;
    }

    return true;
}

bool meshkernel::CurvilinearGridFromSplines::ComputeCurvilinearGrid(CurvilinearGrid& curvilinearGrid)
{
    std::vector<std::vector<size_t>> mIndexsesThisSide(1, std::vector<size_t>(2));
    std::vector<std::vector<size_t>> mIndexsesOtherSide(1, std::vector<size_t>(2));
    std::vector<std::vector<size_t>> nIndexsesThisSide(1, std::vector<size_t>(2));
    std::vector<std::vector<size_t>> nIndexsesOtherSide(1, std::vector<size_t>(2));
    std::vector<std::vector<Point>> gridPointsNDirection(m_gridPoints[0].size(), std::vector<Point>(m_gridPoints.size()));
    std::vector<std::vector<Point>> curvilinearMeshPoints;
    double squaredDistanceTolerance = 1e-12;

    // get the grid sizes in j-direction
    for (int i = 0; i < m_gridPoints[0].size(); i++)
    {
        for (int j = 0; j < m_gridPoints.size(); j++)
        {
            gridPointsNDirection[i][j] = m_gridPoints[j][i];

        }
    }

    size_t startIndex = 0;
    size_t startGridLine = 0;
    while (startIndex < m_gridPoints[0].size())
    {
        int pos = FindIndexes(m_gridPoints[0], startIndex, m_numM, doubleMissingValue, mIndexsesThisSide);

        mIndexsesOtherSide[0][0] = mIndexsesThisSide[0][1] + 2;
        mIndexsesOtherSide[0][1] = mIndexsesOtherSide[0][0] + (mIndexsesThisSide[0][1] - mIndexsesThisSide[0][0]);
        bool isConnected = true;

        size_t minN = m_curvilinearParametersNative.NRefinement;
        size_t maxN = 0;
        size_t minNOther = m_curvilinearParametersNative.NRefinement;
        size_t maxNOther = 0;
        //check if this part is connected to another part
        for (auto i = mIndexsesThisSide[0][0]; i < mIndexsesThisSide[0][1] + 1; ++i)
        {
            pos = FindIndexes(gridPointsNDirection[i], 0, gridPointsNDirection[i].size(), doubleMissingValue, nIndexsesThisSide);
            minN = std::min(minN, nIndexsesThisSide[0][0]);
            maxN = std::max(maxN, nIndexsesThisSide[0][1]);

            size_t mOther = mIndexsesThisSide[0][1] + 2 + (mIndexsesThisSide[0][1] - i);

            if (mOther > m_numM - 1)
            {
                // no more grid available
                isConnected = false;
            }
            else
            {
                double squaredDistance = ComputeSquaredDistance(m_gridPoints[0][i], m_gridPoints[0][mOther], m_splines->m_projection);
                if (squaredDistance > squaredDistanceTolerance)
                {
                    isConnected = false;
                }
                else
                {
                    pos = FindIndexes(gridPointsNDirection[mOther], 0, gridPointsNDirection[mOther].size(), doubleMissingValue, nIndexsesOtherSide);
                    minNOther = std::min(minNOther, nIndexsesOtherSide[0][0]);
                    maxNOther = std::max(maxNOther, nIndexsesOtherSide[0][1]);
                }
            }
        }

        const auto endGridlineIndex = startGridLine + mIndexsesThisSide[0][1] - mIndexsesThisSide[0][0];
        if (isConnected)
        {
            startIndex = mIndexsesOtherSide[0][1] + 2;
        }
        else
        {
            maxNOther = 1;
            startIndex = mIndexsesThisSide[0][1] + 2;
        }

        // increment points
        curvilinearMeshPoints.resize(endGridlineIndex + 1);
        const auto NSize = std::max(curvilinearMeshPoints[0].size(), maxN + maxNOther + 1);

        for (int i = 0; i < curvilinearMeshPoints.size(); i++)
        {
            curvilinearMeshPoints[i].resize(NSize);
        }

        // fill first part
        int columnIncrement = 0;
        for (auto i = startGridLine; i < endGridlineIndex + 1; ++i)
        {
            for (int j = 0; j < maxN + 1; ++j)
            {
                curvilinearMeshPoints[i][j + maxNOther] = m_gridPoints[j][mIndexsesThisSide[0][0] + columnIncrement];
            }
            columnIncrement++;
        }

        columnIncrement = 0;
        for (auto i = startGridLine; i < endGridlineIndex + 1; ++i)
        {
            for (int j = 0; j < maxNOther + 1; ++j)
            {
                curvilinearMeshPoints[i][maxNOther - j] = m_gridPoints[j][mIndexsesOtherSide[0][1] - columnIncrement];
            }
            columnIncrement++;
        }

        startGridLine = endGridlineIndex + 2;
    }

    curvilinearGrid.Set(int(curvilinearMeshPoints.size()), int(curvilinearMeshPoints[0].size()));
    curvilinearGrid.Set(curvilinearMeshPoints);

    return true;
}

bool meshkernel::CurvilinearGridFromSplines::GetSubIntervalAndGridLayer(int layer, int& gridLayer, int& subLayerIndex)
{
    gridLayer = layer - 1;
    int sum = std::accumulate(m_subLayerGridPoints.begin(), m_subLayerGridPoints.end(), 0);
    if (layer >= sum)
    {
        subLayerIndex = -1;
        return true;
    }

    subLayerIndex = 0;
    sum = m_subLayerGridPoints[0] + 1;
    while (sum <= layer && subLayerIndex < m_maxNumCenterSplineHeights)
    {
        subLayerIndex = subLayerIndex + 1;
        sum += m_subLayerGridPoints[subLayerIndex];
    }
    gridLayer = layer - sum + m_subLayerGridPoints[subLayerIndex];

    return true;
}

bool meshkernel::CurvilinearGridFromSplines::GrowLayer(int layerIndex)
{
    assert(layerIndex - 1 >= 0);
    std::vector<Point> velocityVectorAtGridPoints(m_numM);
    bool successful = ComputeVelocitiesAtGridPoints(layerIndex - 1, velocityVectorAtGridPoints);
    if (!successful)
    {
        return false;
    }

    std::vector<Point> activeLayerPoints(m_gridPoints[layerIndex - 1]);
    for (int m = 0; m < velocityVectorAtGridPoints.size(); ++m)
    {
        if (!velocityVectorAtGridPoints[m].IsValid())
        {
            m_gridPoints[layerIndex - 1][m] = { doubleMissingValue, doubleMissingValue };
            activeLayerPoints[m] = { doubleMissingValue, doubleMissingValue };
        }
    }

    const auto numGridPoints = m_gridPoints.size() * m_gridPoints[0].size();
    std::vector<std::vector<int>> gridPointsIndexses(numGridPoints, std::vector<int>(2, -1));
    std::vector<Point> frontGridPoints(numGridPoints);
    int numFrontPoints;
    successful = FindFront(gridPointsIndexses, frontGridPoints, numFrontPoints);
    if (!successful)
    {
        return false;
    }

    std::vector<Point> frontVelocities(numGridPoints);
    successful = CopyVelocitiesToFront(layerIndex - 1, 
                                    velocityVectorAtGridPoints, 
                                    numFrontPoints,
                                    gridPointsIndexses, 
                                    frontGridPoints, 
                                    frontVelocities);
    if (!successful)
    {
        return false;
    }

    double totalTimeStep = 0.0;
    std::vector<Point> gridLine(m_gridPoints[layerIndex - 1]);
    double localTimeStep = 0.0;
    double otherTimeStep = std::numeric_limits<double>::max();
    std::vector<int> newValidFrontNodes(numGridPoints);

    while (totalTimeStep < m_timeStep)
    {
        // Copy old front velocities
        newValidFrontNodes = m_validFrontNodes;

        for (int i = 0; i < m_validFrontNodes.size(); ++i)
        {
            if(m_validFrontNodes[i]<=0)
            {
                activeLayerPoints[i] = { doubleMissingValue, doubleMissingValue };
            }
        }

        std::vector<double> maximumGridLayerGrowTime(newValidFrontNodes.size(), std::numeric_limits<double>::max());
        successful = ComputeMaximumGridLayerGrowTime(activeLayerPoints, velocityVectorAtGridPoints, maximumGridLayerGrowTime);
        if (!successful)
        {
            return false;
        }
        localTimeStep = std::min(m_timeStep - totalTimeStep, *std::min_element(maximumGridLayerGrowTime.begin(), maximumGridLayerGrowTime.end()));

        if (m_splinesToCurvilinearParametersNative.CheckFrontCollisions)
        {
            //TODO: implement front collisions
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

        for (int i = 0; i < newValidFrontNodes.size() - 2; ++i)
        {
            if (newValidFrontNodes[i + 1] == 1 && newValidFrontNodes[i] == 0 && newValidFrontNodes[i + 2] == 0)
            {
                newValidFrontNodes[i + 1] = 0;
            }
        }

        m_validFrontNodes = newValidFrontNodes;

        for (int i = 0; i < velocityVectorAtGridPoints.size(); ++i)
        {
            if (m_validFrontNodes[i] == 1 && velocityVectorAtGridPoints[i].IsValid())
            {
                if (velocityVectorAtGridPoints[i].x == 0.0 && velocityVectorAtGridPoints[i].y == 0.0)
                {
                    continue;
                }
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
            successful = ComputeVelocitiesAtGridPoints(layerIndex, velocityVectorAtGridPoints);
            if (!successful)
            {
                return false;
            }

            for (int i = 0; i < m_numM; ++i)
            {
                // Disable points that have no valid normal vector
                // Remove stationary points
                if (!frontVelocities[i].IsValid() || m_validFrontNodes[i] == 0)
                {
                    activeLayerPoints[i] = { doubleMissingValue, doubleMissingValue };
                }
            }

            successful = FindFront(gridPointsIndexses, frontGridPoints, numFrontPoints);
            if (!successful)
            {
                return false;
            }

            successful = CopyVelocitiesToFront(layerIndex - 1, velocityVectorAtGridPoints, numFrontPoints,
                gridPointsIndexses, frontGridPoints, frontVelocities);

            if (!successful)
            {
                return false;
            }
        }
    }

    if (layerIndex >= 2)
    {
        for (int i = 1; i < m_numM - 1; ++i)
        {

            if (!activeLayerPoints[i].IsValid())
            {
                continue;
            }
            double cosphi = NormalizedInnerProductTwoSegments(m_gridPoints[layerIndex - 2][i],
                m_gridPoints[layerIndex - 1][i],
                m_gridPoints[layerIndex - 1][i],
                activeLayerPoints[i],
                m_splines->m_projection);
            if (cosphi < -0.5)
            {
                int currentLeftIndex;
                int currentRightIndex;
                GetNeighbours(frontGridPoints, i, currentLeftIndex, currentRightIndex);
                for (int j = currentLeftIndex + 1; j < currentRightIndex; ++j)
                {
                    newValidFrontNodes[j] = 0;
                    m_gridPoints[layerIndex - 1][j] = { doubleMissingValue, doubleMissingValue };
                }
            }
        }
    }

    m_validFrontNodes = newValidFrontNodes;


    return true;
}

bool meshkernel::CurvilinearGridFromSplines::ComputeMaximumGridLayerGrowTime( const std::vector<Point>& coordinates,
                                                                              const std::vector<Point>& velocities,
                                                                              std::vector<double>& maximumGridLayerGrowTime) const 
{
    std::vector<double> edgeWidth(coordinates.size() - 1);
    std::vector<double> edgeIncrement(coordinates.size() - 1);
    double minEdgeWidth = 1e-8;
    double dt = 1.0;
    for (int i = 0; i < coordinates.size() - 1; ++i)
    {
        if (!coordinates[i].IsValid() || !coordinates[i + 1].IsValid())
        {
            continue;
        }

        edgeWidth[i] = Distance(coordinates[i], coordinates[i + 1], m_splines->m_projection);

        if (edgeWidth[i] < minEdgeWidth)
        {
            continue;
        }

        Point firstPointIncremented(coordinates[i] + velocities[i] * dt);
        Point secondPointIncremented(coordinates[i + 1] + velocities[i + 1] * dt);;
        edgeIncrement[i] = InnerProductTwoSegments ( coordinates[i], coordinates[i + 1], firstPointIncremented, secondPointIncremented, m_splines->m_projection) / edgeWidth[i] - edgeWidth[i];
        edgeIncrement[i] = edgeIncrement[i] / dt;
    }

    for (int i = 0; i < coordinates.size() - 1; ++i)
    {
        if (edgeIncrement[i] < 0.0)
        {
            maximumGridLayerGrowTime[i] = -edgeWidth[i] / edgeIncrement[i];
        }
    }

    return true;
}

bool meshkernel::CurvilinearGridFromSplines::CopyVelocitiesToFront(const int layerIndex,
                                              const std::vector<Point>& previousVelocities,
                                              int& numFrontPoints,
                                              std::vector<std::vector<int>>& gridPointsIndexses,
                                              std::vector<Point>& frontGridPoints,
                                              std::vector<Point>& velocities)
{
    int numCornerNodes = 0;
    int p = -1;
    while (p < numFrontPoints)
    {
        p = p + 1;
        if (gridPointsIndexses[p][1] == layerIndex && m_validFrontNodes[gridPointsIndexses[p][0]] == 1)
        {
            velocities[p] = previousVelocities[gridPointsIndexses[p][0]];
            if (!velocities[p].IsValid())
            {
                velocities[p] = { 0.0,0.0 };
            }

            // Check for cornernodes
            int previous = std::max(p - 1, 0);
            std::vector<int> previousIndexses = gridPointsIndexses[previous];
            int next = std::min(p + 1, numFrontPoints);
            std::vector<int> nextIndexses = gridPointsIndexses[next];

            // Check corner nodes
            bool ll = previousIndexses[0] == gridPointsIndexses[p][0] - 1 &&
                previousIndexses[1] == gridPointsIndexses[p][1] &&
                m_validFrontNodes[previousIndexses[0]] == -1;

            bool lr = nextIndexses[0] == gridPointsIndexses[p][0] + 1 &&
                nextIndexses[1] == gridPointsIndexses[p][1] &&
                m_validFrontNodes[nextIndexses[0]] == -1;

            ll = ll || previousIndexses[0] == gridPointsIndexses[p][0] && previousIndexses[1] < gridPointsIndexses[p][1];
            lr = lr || nextIndexses[0] == gridPointsIndexses[p][0] && nextIndexses[1] < gridPointsIndexses[p][1];
            if (ll || lr)
            {
                numCornerNodes++;
                if (numFrontPoints + 1 > frontGridPoints.size())
                {
                    continue;
                }
                for (int i = numFrontPoints; i >= p; --i)
                {
                    frontGridPoints[i + 1] = frontGridPoints[i];
                    velocities[i + 1] = velocities[i];
                    gridPointsIndexses[i + 1] = gridPointsIndexses[i];
                }
                numFrontPoints++;

                if (ll)
                {
                    velocities[p] = { 0.0,0.0 };
                }
                else
                {
                    velocities[p + 1] = { 0.0,0.0 };
                }
                p = p + 1;
            }
        }
    }
    return true;
}


bool meshkernel::CurvilinearGridFromSplines::FindFront( std::vector<std::vector<int>>& gridPointsIndexses,
                                                        std::vector<Point>& frontGridPoints,
                                                        int& numFrontPoints )
{

    std::vector<int> frontPosition(m_gridPoints[0].size() - 2, int(m_gridPoints.size()));
    for (int m = 0; m < frontPosition.size(); ++m)
    {
        for (int n = 0; n < m_gridPoints.size(); ++n)
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
    int currentLeftIndex;
    int currentRightIndex;
    int previousFrontPosition = 0;
    GetNeighbours(m_gridPoints[0], 0, currentLeftIndex, currentRightIndex);
    if (currentLeftIndex == 0)
    {
        frontGridPoints[0] = m_gridPoints[0][0];
        // store front index
        gridPointsIndexses[numFrontPoints][0] = 0;
        gridPointsIndexses[numFrontPoints][1] = 0;
        numFrontPoints++;
    }
    else
    {
        previousFrontPosition = frontPosition[currentLeftIndex];
        frontGridPoints[numFrontPoints] = m_gridPoints[0][frontPosition[0]];
        gridPointsIndexses[numFrontPoints][0] = frontPosition[0];
        gridPointsIndexses[numFrontPoints][1] = 0;
        numFrontPoints++;
    }

    for (int m = 0; m < m_gridPoints[0].size() - 2; ++m)
    {
        GetNeighbours(m_gridPoints[0], m, currentLeftIndex, currentRightIndex);
        int currentFrontPosition = frontPosition[m];
        if (currentFrontPosition >= 0)
        {
            if (previousFrontPosition == -1)
            {
                frontGridPoints[numFrontPoints] = m_gridPoints[0][m];
                gridPointsIndexses[numFrontPoints][0] = m;
                gridPointsIndexses[numFrontPoints][1] = 0;
                numFrontPoints++;
            }
            for (int i = previousFrontPosition + 1; i <= currentFrontPosition; ++i)
            {
                frontGridPoints[numFrontPoints] = m_gridPoints[i][m];
                gridPointsIndexses[numFrontPoints][0] = m;
                gridPointsIndexses[numFrontPoints][1] = i;
                numFrontPoints++;
            }
            for (int i = previousFrontPosition; i > currentFrontPosition; --i)
            {
                frontGridPoints[numFrontPoints] = m_gridPoints[i][m];
                gridPointsIndexses[numFrontPoints][0] = m;
                gridPointsIndexses[numFrontPoints][1] = i;
                numFrontPoints++;
            }

            frontGridPoints[numFrontPoints] = m_gridPoints[currentFrontPosition][m + 1];
            gridPointsIndexses[numFrontPoints][0] = m + 1;
            gridPointsIndexses[numFrontPoints][1] = currentFrontPosition;
            numFrontPoints++;
        }
        else if (previousFrontPosition >= 0)
        {
            for (int i = previousFrontPosition - 1; i >= 0; --i)
            {
                frontGridPoints[numFrontPoints] = m_gridPoints[i][m];
                gridPointsIndexses[numFrontPoints][0] = m;
                gridPointsIndexses[numFrontPoints][1] = i;
                numFrontPoints++;
            }

            frontGridPoints[numFrontPoints] = { doubleMissingValue,doubleMissingValue };
            gridPointsIndexses[numFrontPoints][0] = m;
            gridPointsIndexses[numFrontPoints][1] = -1;
            numFrontPoints++;
        }

        previousFrontPosition = currentFrontPosition;
    }

    // add last j-edge, check for circular connectivity
    int lastPoint = int(m_gridPoints[0].size()) - 2;
    GetNeighbours(m_gridPoints[0], lastPoint, currentLeftIndex, currentRightIndex);
    if (currentRightIndex == m_gridPoints[0].size() - 2)
    {
        for (int i = previousFrontPosition; i >= 0; --i)
        {
            frontGridPoints[numFrontPoints] = m_gridPoints[i][lastPoint];
            gridPointsIndexses[numFrontPoints][0] = lastPoint;
            gridPointsIndexses[numFrontPoints][1] = i;
            numFrontPoints++;
        }

    }

    return true;
}


//comp_vel
bool meshkernel::CurvilinearGridFromSplines::ComputeVelocitiesAtGridPoints
(
    int layerIndex,
    std::vector<Point>& velocityVector
)
{
    std::fill(velocityVector.begin(), velocityVector.end(), Point{ doubleMissingValue,doubleMissingValue });
    Point normalVectorLeft;
    Point normalVectorRight;
    const double cosTolerance = 1e-8;
    double eps = 1e-10;
    for (int m = 0; m < velocityVector.size(); ++m)
    {
        if (!m_gridPoints[layerIndex][m].IsValid())
        {
            continue;
        }

        int currentLeftIndex;
        int currentRightIndex;
        GetNeighbours(m_gridPoints[layerIndex], m, currentLeftIndex, currentRightIndex);

        double squaredLeftRightDistance = ComputeSquaredDistance(m_gridPoints[layerIndex][currentLeftIndex], m_gridPoints[layerIndex][currentRightIndex], m_splines->m_projection);
        double squaredLeftDistance = ComputeSquaredDistance(m_gridPoints[layerIndex][currentLeftIndex], m_gridPoints[layerIndex][m], m_splines->m_projection);
        double squaredRightDistance = ComputeSquaredDistance(m_gridPoints[layerIndex][currentRightIndex], m_gridPoints[layerIndex][m], m_splines->m_projection);

        if (squaredLeftRightDistance <= m_onTopOfEachOtherSquaredTolerance)
        {
            continue;
        }

        if (squaredLeftDistance <= m_onTopOfEachOtherSquaredTolerance || squaredRightDistance <= m_onTopOfEachOtherSquaredTolerance)
        {
            NormalVectorOutside(m_gridPoints[layerIndex][currentRightIndex], m_gridPoints[layerIndex][currentLeftIndex], normalVectorLeft, m_splines->m_projection);
            if (m_splines->m_projection == Projections::spherical)
            {
                normalVectorLeft.x = normalVectorLeft.x * std::cos(degrad_hp * 0.5 * (m_gridPoints[layerIndex][currentLeftIndex].y + m_gridPoints[layerIndex][currentRightIndex].y));
            }
            normalVectorRight = normalVectorLeft;
        }
        else
        {
            NormalVectorOutside(m_gridPoints[layerIndex][m], m_gridPoints[layerIndex][currentLeftIndex], normalVectorLeft, m_splines->m_projection);
            NormalVectorOutside(m_gridPoints[layerIndex][currentRightIndex], m_gridPoints[layerIndex][m], normalVectorRight, m_splines->m_projection);

            if (m_splines->m_projection == Projections::spherical)
            {
                normalVectorLeft.x = normalVectorLeft.x * std::cos(degrad_hp * 0.5 * (m_gridPoints[layerIndex][currentLeftIndex].y + m_gridPoints[layerIndex][m].y));
                normalVectorRight.x = normalVectorRight.x * std::cos(degrad_hp * 0.5 * (m_gridPoints[layerIndex][currentRightIndex].y + m_gridPoints[layerIndex][m].y));
            }
        }

        if (currentLeftIndex == velocityVector.size() - 1)
        {
            continue;
        }

        double cosphi = DotProduct(normalVectorLeft.x, normalVectorRight.x, normalVectorLeft.y, normalVectorRight.y);
        Point leftVelocity = normalVectorLeft * m_edgeVelocities[currentLeftIndex];
        Point rightVelocity = normalVectorRight * m_edgeVelocities[currentRightIndex - 1];
        double rightLeftVelocityRatio = m_edgeVelocities[currentRightIndex - 1] / m_edgeVelocities[currentLeftIndex];

        if (cosphi < -1.0 + cosTolerance)
        {
            continue;
        }

        if (rightLeftVelocityRatio - cosphi > eps && 1.0 / rightLeftVelocityRatio - cosphi > eps || cosphi <= cosTolerance)
        {
            velocityVector[m] = (leftVelocity * (1.0 - rightLeftVelocityRatio * cosphi) +
                rightVelocity * (1.0 - 1.0 / rightLeftVelocityRatio * cosphi)) / (1.0 - cosphi * cosphi);
        }
        else if (cosphi - rightLeftVelocityRatio > eps)
        {
            velocityVector[m] = leftVelocity * rightLeftVelocityRatio / cosphi;
        }
        else
        {
            velocityVector[m] = rightVelocity * 1.0 / (rightLeftVelocityRatio * cosphi);
        }

        if (m_splines->m_projection == Projections::spherical)
        {
            velocityVector[m].x = velocityVector[m].x * one_over_earth_radius * raddeg_hp / std::cos(degrad_hp * m_gridPoints[layerIndex][m].y);
            velocityVector[m].y = velocityVector[m].y * one_over_earth_radius * raddeg_hp;
        }
    }
    return true;
}

bool meshkernel::CurvilinearGridFromSplines::GetNeighbours( const std::vector<Point>& gridPoints,
                                                            int index,
                                                            int& currentLeftIndex,
                                                            int& currentRightIndex) const
{
    bool circularConnection = false;
    currentLeftIndex = index;
    currentRightIndex = index;
    int start = 0;
    int end = int(gridPoints.size()) - 1;

    // left
    while (ComputeSquaredDistance(gridPoints[currentLeftIndex], gridPoints[index], m_splines->m_projection) < m_onTopOfEachOtherSquaredTolerance)
    {
        if (!circularConnection)
        {
            if (currentLeftIndex - 1 < 0)
            {
                break;
            }
        }
        else if (currentLeftIndex - 1 < 0)
        {
            currentLeftIndex = end + 1;
            circularConnection = false;
        }
        if (currentLeftIndex - 1 < 0 || !gridPoints[currentLeftIndex - 1].IsValid())
        {
            break;
        }
        currentLeftIndex--;
    }

    // right
    while (ComputeSquaredDistance(gridPoints[currentRightIndex], gridPoints[index], m_splines->m_projection) < m_onTopOfEachOtherSquaredTolerance)
    {
        if (!circularConnection)
        {
            if (currentRightIndex + 1 > gridPoints.size())
            {
                break;
            }
        }
        else if (currentRightIndex + 1 > gridPoints.size())
        {
            currentRightIndex = start - 1;
            circularConnection = false;
        }
        if (currentRightIndex + 1 >= gridPoints.size() || !gridPoints[currentRightIndex + 1].IsValid())
        {
            break;
        }
        currentRightIndex++;
    }
    return true;
}

bool meshkernel::CurvilinearGridFromSplines::ComputeEdgeVelocities( std::vector<double>& edgeVelocities, // edgevel
                                               std::vector<std::vector<double>>& growFactorOnSubintervalAndEdge, //dgrow1
                                               std::vector<std::vector<int>>& numPerpendicularFacesOnSubintervalAndEdge //nfac1
)
{

    bool successful = ComputeGridHeights();
    if (!successful)
    {
        return false;
    }

    for (int i = 0; i < numPerpendicularFacesOnSubintervalAndEdge[0].size(); ++i)
    {
        numPerpendicularFacesOnSubintervalAndEdge[0][i] = 1;
    }

    for (size_t s = 0; s < m_splines->m_numSplines; s++)
    {
        double maxHeight = std::numeric_limits<double>::min();
        for (int i = 0; i < m_gridHeights[0].size(); ++i)
        {
            if (m_gridHeights[0][i] != doubleMissingValue && m_gridHeights[0][i] > maxHeight)
            {
                maxHeight = m_gridHeights[0][i];
            }
        }

        double firstHeight = std::min(maxHeight, m_splinesToCurvilinearParametersNative.AspectRatio * m_splinesToCurvilinearParametersNative.AverageWidth);

        // Get true crossing splines heights
        int numLeftHeights = m_maxNumCenterSplineHeights;
        int numRightHeights = m_maxNumCenterSplineHeights;
        int numTrueCrossings = 0;
        for (int i = 0; i < m_numCrossingSplines[s]; ++i)
        {
            if (m_type[m_crossingSplinesIndexses[s][i]] != SplineTypes::crossing)
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


        int startGridLineLeft = m_leftGridLineIndex[s];
        int endGridLineLeft = startGridLineLeft + m_numMSplines[s];
        int startGridLineRight = m_rightGridLineIndex[s];
        int endGridLineRight = startGridLineRight + m_numMSplines[s];

        double hh0LeftMaxRatio;
        double hh0RightMaxRatio;
        const int numIterations = 2;
        for (int iter = 0; iter < numIterations; ++iter)
        {
            ComputeVelocitiesSubIntervals(s, startGridLineLeft, endGridLineLeft, numLeftHeights, numRightHeights, firstHeight,
                m_leftGridLineIndex, m_rightGridLineIndex, numPerpendicularFacesOnSubintervalAndEdge, edgeVelocities, hh0LeftMaxRatio);

            ComputeVelocitiesSubIntervals(s, startGridLineRight, endGridLineRight, numLeftHeights, numRightHeights, firstHeight,
                m_rightGridLineIndex, m_leftGridLineIndex, numPerpendicularFacesOnSubintervalAndEdge, edgeVelocities, hh0RightMaxRatio);
        }

        // re-evaluate if growing grid outside is needed

        if (numLeftHeights == 0 && numRightHeights <= 1 ||
            numRightHeights == 0 && numLeftHeights <= 1 ||
            numLeftHeights == numRightHeights == 1)
        {
            m_growGridOutside = true;
        }

        // left part
        int numNLeftExponential = 0;
        if (m_growGridOutside)
        {
            numNLeftExponential = std::min(ComputeNumberExponentialIntervals(hh0LeftMaxRatio), m_curvilinearParametersNative.NRefinement);
        }
        for (int i = startGridLineLeft; i < endGridLineLeft; ++i)
        {
            numPerpendicularFacesOnSubintervalAndEdge[1][i] = numNLeftExponential;
        }

        // right part
        int numNRightExponential = 0;
        if (m_growGridOutside)
        {
            numNRightExponential = std::min(ComputeNumberExponentialIntervals(hh0RightMaxRatio), m_curvilinearParametersNative.NRefinement);
        }
        for (int i = startGridLineRight; i < endGridLineRight; ++i)
        {
            numPerpendicularFacesOnSubintervalAndEdge[1][i] = numNRightExponential;
        }
    }

    // compute local grow factors
    for (size_t s = 0; s < m_splines->m_numSplines; s++)
    {
        if (m_numMSplines[s] < 1)
        {
            continue;
        }

        for (int i = m_leftGridLineIndex[s]; i < m_rightGridLineIndex[s] + m_numMSplines[s]; ++i)
        {
            if (m_gridLine[i].x == doubleMissingValue || m_gridLine[i + 1].x == doubleMissingValue || numPerpendicularFacesOnSubintervalAndEdge[1][i] < 1)
            {
                continue;
            }
            bool successfull = ComputeGrowFactor( m_gridHeights[1][i],
                                                  edgeVelocities[i],
                                                  numPerpendicularFacesOnSubintervalAndEdge[1][i],
                                                  growFactorOnSubintervalAndEdge[1][i] );
            if (!successfull)
            {
                growFactorOnSubintervalAndEdge[1][i] = 1.0;
            }
        }
    }

    return true;
}

///comp_dgrow: this is another root finding algorithm, could go in the general part
bool meshkernel::CurvilinearGridFromSplines::ComputeGrowFactor(
    double totalGridHeight,
    double firstGridLayerHeight,
    int numberOfGridLayers,
    double& result)
{
    // eheight m_gridHeights
    double aspectRatioGrowFactor = 1.0;
    double heightDifference = ComputeTotalExponentialHeight(aspectRatioGrowFactor, firstGridLayerHeight, numberOfGridLayers) - totalGridHeight;

    double deps = 0.01;
    double aspectRatioGrowFactorIncremented = 1.0 + deps;
    double heightDifferenceIncremented = ComputeTotalExponentialHeight(aspectRatioGrowFactorIncremented, firstGridLayerHeight, numberOfGridLayers) - totalGridHeight;

    const double tolerance = 1e-8;
    const int numIterations = 1000;
    const double relaxationFactor = 0.5;
    double oldAspectRatio;
    double oldHeightDifference = heightDifference;

    if (std::abs(heightDifferenceIncremented) > tolerance&& std::abs(heightDifferenceIncremented - heightDifference) > tolerance)
    {
        for (int i = 0; i < numIterations; ++i)
        {
            oldAspectRatio = aspectRatioGrowFactor;
            oldHeightDifference = heightDifference;

            aspectRatioGrowFactor = aspectRatioGrowFactorIncremented;
            heightDifference = heightDifferenceIncremented;

            aspectRatioGrowFactorIncremented = aspectRatioGrowFactor - relaxationFactor * heightDifference / (heightDifference - oldHeightDifference + 1e-16) * (aspectRatioGrowFactor - oldAspectRatio);
            heightDifferenceIncremented = ComputeTotalExponentialHeight(aspectRatioGrowFactorIncremented, firstGridLayerHeight, numberOfGridLayers) - totalGridHeight;

            if (std::abs(oldHeightDifference) < tolerance)
            {
                break;
            }
        }
    }

    if (oldHeightDifference > tolerance)
    {
        result = doubleMissingValue;
        return false;
    }

    result = aspectRatioGrowFactorIncremented;
    return true;
}

double meshkernel::CurvilinearGridFromSplines::ComputeTotalExponentialHeight(double aspectRatioGrowFactor, double firstGridLayerHeights, int numberOfGridLayers) const
{
    double height;
    if (std::abs(aspectRatioGrowFactor - 1.0) > 1e-8)
    {
        height = (std::pow(aspectRatioGrowFactor, numberOfGridLayers) - 1.0) / (aspectRatioGrowFactor - 1.0) * firstGridLayerHeights;
    }
    else
    {
        height = firstGridLayerHeights * numberOfGridLayers;
    }
    return height;
}


int meshkernel::CurvilinearGridFromSplines::ComputeNumberExponentialIntervals(const double hhMaxRatio) const
{
    int numIntervals = 0;
    if (m_splinesToCurvilinearParametersNative.AspectRatioGrowFactor - 1.0 > 1e-8)
    {
        numIntervals = int(std::floor(std::log((m_splinesToCurvilinearParametersNative.AspectRatioGrowFactor - 1.0) * hhMaxRatio + 1.0) / log(m_splinesToCurvilinearParametersNative.AspectRatioGrowFactor)));
    }
    else
    {
        numIntervals = int(std::floor(0.999 + hhMaxRatio));
    }
    return numIntervals;
}

bool meshkernel::CurvilinearGridFromSplines::ComputeVelocitiesSubIntervals( const int s,
                                                       const int startGridLineIndex,
                                                       const int endGridLineIndex,
                                                       const int numHeights,
                                                       const int numOtherSideHeights,
                                                       const double firstHeight,
                                                       const std::vector<int>& gridLineIndex,
                                                       const std::vector<int>& otherGridLineIndex,
                                                       std::vector<std::vector<int>>& numPerpendicularFacesOnSubintervalAndEdge,
                                                       std::vector<double>& edgeVelocities,
                                                       double& hh0MaxRatio)
{
    hh0MaxRatio = 0.0;
    if (numHeights > 1 && numHeights == numOtherSideHeights || numHeights > numOtherSideHeights)
    {
        double maxHeight = *std::max_element(m_gridHeights[0].begin() + startGridLineIndex, m_gridHeights[0].begin() + endGridLineIndex);

        int numNUniformPart = int(std::floor(maxHeight / firstHeight + 0.99999));
        numNUniformPart = std::min(numNUniformPart, m_maxNUniformPart);

        for (int i = startGridLineIndex; i < endGridLineIndex; ++i)
        {
            numPerpendicularFacesOnSubintervalAndEdge[0][i] = numNUniformPart;
            edgeVelocities[i] = m_gridHeights[0][i] / numNUniformPart;
            hh0MaxRatio = std::max(hh0MaxRatio, m_gridHeights[1][i] / edgeVelocities[i]);
        }
    }
    else
    {
        // only one subinterval: no uniform part
        int numNUniformPart = 0;
        for (int i = startGridLineIndex; i < endGridLineIndex; ++i)
        {
            numPerpendicularFacesOnSubintervalAndEdge[0][i] = numNUniformPart;
            edgeVelocities[i] = firstHeight;

            //compare with other side of spline
            int otherSideIndex = otherGridLineIndex[s] + m_numMSplines[s] - (i - gridLineIndex[s] + 1);

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

            for (int j = 1; j < m_maxNumCenterSplineHeights; ++j)
            {
                m_gridHeights[j][i] = m_gridHeights[j - 1][i];
            }

            for (int j = startGridLineIndex; j < endGridLineIndex; ++j)
            {

                hh0MaxRatio = std::max(hh0MaxRatio, m_gridHeights[1][j] / edgeVelocities[j]);
            }
        }
    }

    return true;
}

bool meshkernel::CurvilinearGridFromSplines::ComputeGridHeights()
{
    const auto numSplines = m_splines->m_numSplines;

    m_gridHeights.resize(m_maxNumCenterSplineHeights, std::vector<double>(m_numM - 1));
    std::fill(m_gridHeights.begin(), m_gridHeights.end(), std::vector<double>(m_numM - 1, doubleMissingValue));

    std::vector<std::vector<double>> heightsLeft(m_maxNumCenterSplineHeights, std::vector<double>(m_curvilinearParametersNative.MRefinement, 0.0));
    std::vector<std::vector<double>> heightsRight(m_maxNumCenterSplineHeights, std::vector<double>(m_curvilinearParametersNative.MRefinement, 0.0));
    std::vector<double> edgesCenterPoints(m_numM, 0.0);
    std::vector<double> crossingSplinesDimensionalCoordinates(numSplines, 0.0);
    std::vector<int> numHeightsLeft(numSplines, 0);
    std::vector<int> numHeightsRight(numSplines, 0);
    std::vector<double> localSplineDerivatives(numSplines, 0.0);
    std::vector<int> localValidSplineIndexes(numSplines, 0);

    for (auto s = 0; s < numSplines; s++)
    {
        if (m_type[s] != SplineTypes::central)
        {
            continue;
        }

        int numM = m_numMSplines[s];

        // Get the minimum number of sub-intervals in the cross splines for this center spline
        int minNumLeftIntervals = *std::min_element(m_numCrossSplineLeftHeights[s].begin(), m_numCrossSplineLeftHeights[s].end());
        int minNumRightIntervals = *std::min_element(m_numCrossSplineRightHeights[s].begin(), m_numCrossSplineRightHeights[s].end());

        std::fill(heightsLeft[0].begin(), heightsLeft[0].begin() + numM, m_maximumGridHeights[s]);
        std::fill(heightsRight[0].begin(), heightsRight[0].begin() + numM, m_maximumGridHeights[s]);

        if (m_numCrossingSplines[s] == 1)
        {
            // only one crossing spline present: 
            for (int i = 0; i < minNumLeftIntervals; ++i)
            {
                std::fill(heightsLeft[i].begin(), heightsLeft[i].begin() + numM, m_crossSplineRightHeights[s][i][0]);

            }
            for (int i = 0; i < minNumRightIntervals; ++i)
            {
                std::fill(heightsRight[i].begin(), heightsRight[i].begin() + numM, m_crossSplineLeftHeights[s][i][0]);
            }
        }
        else
        {
            int leftGridLineIndex = m_leftGridLineIndex[s];
            edgesCenterPoints[0] = m_splines->GetSplineLength(s, 0, m_gridLineDimensionalCoordinates[leftGridLineIndex]);
            for (int i = 0; i < numM; ++i)
            {
                edgesCenterPoints[i + 1] = edgesCenterPoints[i] + m_splines->GetSplineLength(s, m_gridLineDimensionalCoordinates[leftGridLineIndex + i], m_gridLineDimensionalCoordinates[leftGridLineIndex + i + 1]);
            }

            // compute at edge center points
            for (int i = 0; i < numM; ++i)
            {
                edgesCenterPoints[i] = 0.5 * (edgesCenterPoints[i] + edgesCenterPoints[i + 1]);
            }
            edgesCenterPoints[numM] = doubleMissingValue;

            //compute center spline path length of cross splines
            crossingSplinesDimensionalCoordinates[0] = m_splines->GetSplineLength(s, 0.0, m_crossSplineCoordinates[s][0]);
            for (int i = 0; i < m_numCrossingSplines[s] - 1; ++i)
            {
                crossingSplinesDimensionalCoordinates[i + 1] = crossingSplinesDimensionalCoordinates[i] + m_splines->GetSplineLength(s, m_crossSplineCoordinates[s][i], m_crossSplineCoordinates[s][i + 1]);
            }

            for (int i = 0; i < m_numCrossingSplines[s]; ++i)
            {
                numHeightsLeft[i] = m_numCrossSplineLeftHeights[s][i];
                numHeightsRight[i] = m_numCrossSplineRightHeights[s][i];
            }


            for (int j = 0; j < m_maxNumCenterSplineHeights; ++j)
            {
                AddValueToVector(numHeightsLeft, -1);
                AddValueToVector(numHeightsRight, -1);

                bool successful = FindNearestCrossSplines(s, j,
                    numHeightsLeft,
                    edgesCenterPoints,
                    m_crossSplineLeftHeights[s],
                    localValidSplineIndexes,
                    localSplineDerivatives,
                    crossingSplinesDimensionalCoordinates,
                    heightsLeft);

                if (!successful)
                {
                    return false;
                }

                successful = FindNearestCrossSplines(s, j,
                    numHeightsRight,
                    edgesCenterPoints,
                    m_crossSplineRightHeights[s],
                    localValidSplineIndexes,
                    localSplineDerivatives,
                    crossingSplinesDimensionalCoordinates,
                    heightsRight);

                if (!successful)
                {
                    return false;
                }
            }
        }


        // store grid height
        for (int j = 0; j < m_maxNumCenterSplineHeights; ++j)
        {
            for (int i = 0; i < m_numMSplines[s]; ++i)
            {
                m_gridHeights[j][m_leftGridLineIndex[s] + i] = heightsLeft[j][i];
                m_gridHeights[j][m_rightGridLineIndex[s] + m_numMSplines[s] - i - 1] = heightsRight[j][i];
            }
        }
    }

    return true;
}

bool meshkernel::CurvilinearGridFromSplines::FindNearestCrossSplines(const int s,
                                                const int j,
                                                const std::vector<int>& numHeightsLeft,
                                                const std::vector<double>& edgesCenterPoints,
                                                const std::vector<std::vector<double>>& crossSplineLeftHeights,
                                                std::vector<int>& localValidSplineIndexes,
                                                std::vector<double>& localSplineDerivatives,
                                                std::vector<double>& crossingSplinesDimensionalCoordinates,
                                                std::vector<std::vector<double>>& heights)
{
    int numValid;
    GetValidSplineIndexses(s, m_numCrossingSplines[s], numHeightsLeft, localValidSplineIndexes, numValid);

    // no sub-heights to compute 
    if (numValid == 0)
    {
        return true;
    }

    int numM = m_numMSplines[s];
    std::vector<double> localCornerPoints(numValid);

    // TODO: strided memory access
    for (int i = 0; i < numValid; ++i)
    {
        int index = localValidSplineIndexes[i];
        localCornerPoints[i] = crossSplineLeftHeights[index][j];
    }

    m_splines->SecondOrderDerivative(localCornerPoints, numValid, localSplineDerivatives);

    crossingSplinesDimensionalCoordinates[0] = m_splines->GetSplineLength(s, 0.0, m_crossSplineCoordinates[s][0]);
    for (int i = 0; i < numM; ++i)
    {
        int leftIndex = 0;
        double leftCoordinate = crossingSplinesDimensionalCoordinates[localValidSplineIndexes[leftIndex]];
        int rightIndex = std::min(1, numValid - 1);
        double rightCoordinate = crossingSplinesDimensionalCoordinates[localValidSplineIndexes[rightIndex]];
        // Find two nearest cross splines
        while (rightCoordinate < edgesCenterPoints[i] && rightIndex < numValid)
        {
            leftIndex = rightIndex;
            leftCoordinate = rightCoordinate;
            rightIndex++;
            rightCoordinate = crossingSplinesDimensionalCoordinates[localValidSplineIndexes[rightIndex]];
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
        else
        {
            rightIndex = leftIndex;
        }

        factor = std::max(std::min(double(leftIndex + 1) + factor - 1.0, double(numValid - 1)), 0.0);

        bool successful = InterpolateSplinePoint(localCornerPoints, localSplineDerivatives, factor, heights[j][i]);
        if (!successful)
        {
            return false;
        }
    }

    return true;
}


bool meshkernel::CurvilinearGridFromSplines::GetValidSplineIndexses(const int s, const int numValues, const std::vector<int>& v, std::vector<int>& validIndexses, int& numValid) const
{
    numValid = 0;
    for (int i = 0; i < numValues; ++i)
    {
        if (v[i] >= 0)
        {
            validIndexses[numValid] = i;
            numValid++;
        }
    }
    return true;
};

/// get_crosssplines
/// compute the intersection of two splines, one must have only two nodes
bool meshkernel::CurvilinearGridFromSplines::GetSplineIntersections(const int index)
{
    m_numCrossingSplines[index] = 0;
    const auto numSplines = m_splines->m_numSplines;
    std::fill(m_crossingSplinesIndexses[index].begin(), m_crossingSplinesIndexses[index].end(), -1);
    std::fill(m_isLeftOriented[index].begin(), m_isLeftOriented[index].end(), true);
    std::fill(m_crossSplineCoordinates[index].begin(), m_crossSplineCoordinates[index].end(), std::numeric_limits<double>::max());
    std::fill(m_cosCrossingAngle[index].begin(), m_cosCrossingAngle[index].end(), doubleMissingValue);

    for (auto s = 0; s < numSplines; ++s)
    {
        // a crossing is a spline with 2 nodes and another with more than 2 nodes
        if (m_splines->m_numSplineNodes[s] == 2 && m_splines->m_numSplineNodes[index] == 2 ||
            m_splines->m_numSplineNodes[s] > 2 && m_splines->m_numSplineNodes[index] > 2)
        {
            continue;
        }

        double crossProductIntersection;
        Point intersectionPoint;
        double firstSplineRatio;
        double secondSplineRatio;
        bool crossing = m_splines->GetSplinesIntersection(index, s, crossProductIntersection, intersectionPoint, firstSplineRatio, secondSplineRatio);

        if (std::abs(crossProductIntersection) < m_splinesToCurvilinearParametersNative.MinimumCosineOfCrossingAngles)
        {
            crossing = false;
        }


        if (crossing)
        {
            m_numCrossingSplines[index]++;
            m_crossingSplinesIndexses[index][s] = s;
            if (crossProductIntersection > 0.0)
            {
                m_isLeftOriented[index][s] = false;
            }
            m_crossSplineCoordinates[index][s] = firstSplineRatio;
            m_cosCrossingAngle[index][s] = crossProductIntersection;
        }
    }

    const auto sortedIndexses = SortedIndexes(m_crossSplineCoordinates[index]);
    ReorderVector(m_crossSplineCoordinates[index], sortedIndexses);
    ReorderVector(m_crossingSplinesIndexses[index], sortedIndexses);
    ReorderVector(m_isLeftOriented[index], sortedIndexses);

    return true;
}

bool meshkernel::CurvilinearGridFromSplines::MakeAllGridLines(bool isSpacingCurvatureAdapeted)
{

    int numCenterSplines = 0;
    for (auto s = 0; s < m_splines->m_numSplines; ++s)
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
        return false;
    }


    int gridLineIndex = 0;
    for (size_t s = 0; s < m_splines->m_numSplines; ++s)
    {
        //center splines only
        if (m_type[s] != SplineTypes::central)
        {
            continue;
        }

        // upper bound of m_gridLine, with two sides of spline and two missing values added
        int sizeGridLine = gridLineIndex + 1 + 2 * (m_curvilinearParametersNative.MRefinement + 1) + 2;
        // increase size
        ResizeVectorIfNeeded(sizeGridLine, m_gridLine, { doubleMissingValue,doubleMissingValue });
        m_gridLineDimensionalCoordinates.resize(sizeGridLine);

        if (gridLineIndex > 0)
        {
            gridLineIndex++;
            m_gridLine[gridLineIndex] = { doubleMissingValue,doubleMissingValue };
            m_gridLineDimensionalCoordinates[gridLineIndex] = doubleMissingValue;
        }

        m_leftGridLineIndex[s] = gridLineIndex;

        int numM = 0;
        bool successful = MakeGridLine(s, gridLineIndex, m_gridLine, m_gridLineDimensionalCoordinates, numM);
        if (!successful)
        {
            return false;
        }

        gridLineIndex = gridLineIndex + numM + 1;
        m_gridLine[gridLineIndex] = Point{ doubleMissingValue,  doubleMissingValue };
        m_gridLineDimensionalCoordinates[gridLineIndex] = doubleMissingValue;
        gridLineIndex++;

        //add other side of gridline
        m_rightGridLineIndex[s] = gridLineIndex;
        int i = m_rightGridLineIndex[s] - 1;
        for (int j = m_rightGridLineIndex[s] - 1; j >= m_leftGridLineIndex[s]; --j)
        {
            m_gridLine[i] = m_gridLine[j];
            m_gridLineDimensionalCoordinates[i] = m_gridLineDimensionalCoordinates[j];
            ++i;
        }

        //compute new (actual) grid size
        //new size   old size   both sides of spline   DMISS between both sides
        gridLineIndex = gridLineIndex + numM + 1;

        m_gridLine[gridLineIndex] = Point{ doubleMissingValue,  doubleMissingValue };
        m_gridLineDimensionalCoordinates[gridLineIndex] = doubleMissingValue;
        m_numMSplines[s] = numM;
        m_numM = gridLineIndex;
    }

    return true;
}

bool meshkernel::CurvilinearGridFromSplines::MakeGridLine( int splineIndex,
                                                         int startingIndex,
                                                         std::vector<Point>& gridLine,
                                                         std::vector<double>& adimensionalCoordinates,
                                                         int& numM )
{
    // first estimation of nodes along m
    numM = 1 + int(std::floor(m_splines->m_splinesLength[splineIndex] / m_splinesToCurvilinearParametersNative.AverageWidth));
    numM = std::min(numM, m_curvilinearParametersNative.MRefinement);

    double endSplineAdimensionalCoordinate = m_splines->m_numSplineNodes[splineIndex] - 1;
    double splineLength = m_splines->GetSplineLength(splineIndex, 0.0, endSplineAdimensionalCoordinate, 10, m_splinesToCurvilinearParametersNative.CurvatureAdapetedGridSpacing, m_maximumGridHeights[splineIndex]);

    gridLine[startingIndex] = m_splines->m_splineNodes[splineIndex][0];

    double currentMaxWidth = std::numeric_limits<double>::max();
    std::vector<double> distances(numM);
    std::vector<double> adimensionalDistances(numM);
    std::vector<Point> points(numM);
    while (currentMaxWidth > m_splinesToCurvilinearParametersNative.AverageWidth)
    {
        currentMaxWidth = 0.0;
        for (int n = 0; n < numM; ++n)
        {
            distances[n] = splineLength * double(n + 1.0) / double(numM);
        }
        
        m_splines->InterpolatePointsOnSpline(splineIndex,
                                             m_maximumGridHeights[splineIndex],
                                             m_splinesToCurvilinearParametersNative.CurvatureAdapetedGridSpacing,
                                             distances,
                                             points,
                                             adimensionalDistances);

        for (int n = 0; n < numM; ++n)
        {
            int index = startingIndex + n + 1;
            adimensionalCoordinates[index] = adimensionalDistances[n];
            gridLine[index] = points[n];
            currentMaxWidth = std::max(currentMaxWidth, Distance(gridLine[index-1], gridLine[index], m_splines->m_projection));
        }

        // a gridline is computed
        if (currentMaxWidth < m_splinesToCurvilinearParametersNative.AverageWidth || numM == m_curvilinearParametersNative.MRefinement)
        {
            break;
        }

        // room for sub-division
        if (currentMaxWidth > m_splinesToCurvilinearParametersNative.AverageWidth)
        {
            numM = std::min(std::max(int(m_curvilinearParametersNative.MRefinement / m_maximumGridHeights[splineIndex] * numM), numM + 1), m_curvilinearParametersNative.MRefinement);
            distances.resize(numM);
            adimensionalDistances.resize(numM);
            points.resize(numM);
        }
    }
    return true;
}

bool meshkernel::CurvilinearGridFromSplines::ComputeSplineProperties(const bool restoreOriginalProperties)
{
    bool successful = AllocateSplinesProperties();
    if (!successful)
    {
        return false;
    }

    for (size_t s = 0; s < m_splines->m_numSplines; ++s)
    {
        successful = GetSplineIntersections(s);

        if (!successful)
        {
            return false;
        }
    }
    // select all non-cross splines only
    for (size_t s = 0; s < m_splines->m_numSplines; ++s)
    {
        m_type[s] = SplineTypes::crossing;
        // select all non-cross splines for growing the grid
        if (m_splines->m_numSplineNodes[s] > 2)
        {
            m_type[s] = SplineTypes::central;
        }
    }
    // check the cross splines. The center spline is the middle spline that crosses the cross spline
    for (size_t s = 0; s < m_splines->m_numSplines; ++s)
    {
        // only crossing splines with one or more center spline
        if (m_splines->m_numSplineNodes[s] != 2 || m_numCrossingSplines[s] < 1)
        {
            continue;
        }

        int middleCrossingSpline = std::min(m_numCrossingSplines[s] / 2, m_numCrossingSplines[s]);
        int crossingSplineIndex = m_crossingSplinesIndexses[s][middleCrossingSpline];

        // if m_numIntersectingSplines[s] is even, check if the middle spline has already been assigned as a bounding spline
        if (m_type[crossingSplineIndex] != SplineTypes::central && 2 * crossingSplineIndex == m_numCrossingSplines[s])
        {
            middleCrossingSpline = std::min(middleCrossingSpline + 1, m_numCrossingSplines[s] - 1);
            crossingSplineIndex = m_crossingSplinesIndexses[s][middleCrossingSpline];
        }

        if (m_type[crossingSplineIndex] == SplineTypes::central)
        {
            // associate bounding splines with the middle spline
            for (int i = 0; i < middleCrossingSpline; ++i)
            {
                int index = m_crossingSplinesIndexses[s][i];
                m_type[index] = SplineTypes::lateral; // lateral spline
                m_centralSplineIndex[index] = -crossingSplineIndex;

            }
            for (int i = middleCrossingSpline + 1; i < m_numCrossingSplines[s]; ++i)
            {
                int index = m_crossingSplinesIndexses[s][i];
                m_type[index] = SplineTypes::lateral; // lateral spline
                m_centralSplineIndex[index] = -crossingSplineIndex;
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
        for (size_t s = m_numOriginalSplines; s < m_splines->m_numSplines; ++s)
        {
            m_type[s] = SplineTypes::arficial;
        }
    }

    bool successfull = ComputeHeights();
    return successfull;
}

bool meshkernel::CurvilinearGridFromSplines::ComputeHeights()
{
    for (size_t i = 0; i < m_splines->m_numSplines; ++i)
    {
        // Heights should be computed only for center splines
        if (m_splines->m_numSplineNodes[i] <= 2)
        {
            continue;
        }
        for (int j = 0; j < m_numCrossingSplines[i]; ++j)
        {
            int intersectingSplineIndex = m_crossingSplinesIndexses[i][j];
            bool successful = ComputeSubHeights(i, j);
            if (!successful)
            {
                return false;
            }
        }
    }

    // compute m_maximumGridHeight
    for (size_t s = 0; s < m_splines->m_numSplines; ++s)
    {
        if (m_numCrossingSplines[s] == 0)
        {
            m_maximumGridHeights[s] = m_splinesToCurvilinearParametersNative.AspectRatio * m_splines->m_splinesLength[s];
            continue;
        }
        double maximumHeight = 0.0;
        for (int c = 0; c < m_numCrossingSplines[s]; ++c)
        {
            double sumLeftHeights = 0.0;
            for (int ss = 0; ss < m_numCrossSplineLeftHeights[s][c]; ++ss)
            {
                sumLeftHeights += m_crossSplineLeftHeights[s][c][ss];
            }
            double sumRightHeights = 0.0;
            for (int ss = 0; ss < m_numCrossSplineRightHeights[s][c]; ++ss)
            {
                sumRightHeights += m_crossSplineRightHeights[s][c][ss];
            }
            maximumHeight = std::max(maximumHeight, std::max(sumLeftHeights, sumRightHeights));
        }

        m_maximumGridHeights[s] = maximumHeight;
    }
    return true;
}

bool meshkernel::CurvilinearGridFromSplines::ComputeSubHeights( int centerSplineIndex,  int crossingSplineLocalIndex)
{
    // find center spline index
    int centerSplineLocalIndex = 0;
    int crossingSplineIndex = m_crossingSplinesIndexses[centerSplineIndex][crossingSplineLocalIndex]; //js
    for (int s = 0; s < m_numCrossingSplines[crossingSplineIndex]; ++s)
    {
        if (m_crossingSplinesIndexses[crossingSplineIndex][s] == centerSplineIndex)
        {
            centerSplineLocalIndex = s;
            break;
        }
    }

    // right part
    int numSubIntervalsRight = 0;
    int rightCenterSplineIndex = centerSplineLocalIndex;
    int leftCenterSplineIndex;
    m_crossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex].resize(m_maxNumCenterSplineHeights, 0);
    for (int s = centerSplineLocalIndex; s < m_numCrossingSplines[crossingSplineIndex] - 1; ++s)
    {
        if (numSubIntervalsRight >= m_maxNumCenterSplineHeights)
        {
            break;
        }
        if (m_centralSplineIndex[m_crossingSplinesIndexses[crossingSplineIndex][s + 1]] != -centerSplineIndex)
        {
            continue;
        }
        leftCenterSplineIndex = rightCenterSplineIndex;
        rightCenterSplineIndex = s + 1;
        m_crossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex][numSubIntervalsRight] = m_splines->GetSplineLength(crossingSplineIndex,
            m_crossSplineCoordinates[crossingSplineIndex][leftCenterSplineIndex], m_crossSplineCoordinates[crossingSplineIndex][rightCenterSplineIndex]);
        numSubIntervalsRight++;
    }

    m_crossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex][numSubIntervalsRight] = m_splines->GetSplineLength(crossingSplineIndex,
        m_crossSplineCoordinates[crossingSplineIndex][rightCenterSplineIndex], m_splines->m_numSplineNodes[crossingSplineIndex] - 1);

    numSubIntervalsRight++;
    std::fill(m_crossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex].begin() + numSubIntervalsRight, m_crossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex].end(), 0.0);
    m_numCrossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex] = numSubIntervalsRight;

    // left part
    int numSubIntervalsLeft = 0;
    leftCenterSplineIndex = centerSplineLocalIndex;
    m_crossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex].resize(m_maxNumCenterSplineHeights, 0);
    for (int s = centerSplineLocalIndex; s >= 1; --s)
    {
        if (numSubIntervalsLeft >= m_maxNumCenterSplineHeights)
        {
            break;
        }
        if (m_centralSplineIndex[m_crossingSplinesIndexses[crossingSplineIndex][s - 1]] != -centerSplineIndex)
        {
            continue;
        }
        rightCenterSplineIndex = leftCenterSplineIndex;
        leftCenterSplineIndex = s - 1;
        m_crossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex][numSubIntervalsLeft] = m_splines->GetSplineLength(crossingSplineIndex,
            m_crossSplineCoordinates[crossingSplineIndex][leftCenterSplineIndex], m_crossSplineCoordinates[crossingSplineIndex][rightCenterSplineIndex]);
        numSubIntervalsLeft++;
    }

    m_crossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex][numSubIntervalsLeft] = m_splines->GetSplineLength(crossingSplineIndex,
        0.0, m_crossSplineCoordinates[crossingSplineIndex][leftCenterSplineIndex]);

    numSubIntervalsLeft++;
    std::fill(m_crossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex].begin() + numSubIntervalsLeft, m_crossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex].end(), 0.0);
    m_numCrossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex] = numSubIntervalsLeft;

    // if not left oriented, swap
    if (!m_isLeftOriented[centerSplineIndex][centerSplineLocalIndex])
    {
        m_numCrossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex] = numSubIntervalsRight;
        m_numCrossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex] = numSubIntervalsLeft;

        std::vector<double> leftSubIntervalsTemp(m_crossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex]);
        m_crossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex] = m_crossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex];
        m_crossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex] = leftSubIntervalsTemp;
    }

    return true;
}