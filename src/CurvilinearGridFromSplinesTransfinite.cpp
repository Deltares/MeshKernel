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
#include "CurvilinearGridFromSplinesTransfinite.hpp"
#include "Splines.hpp"
#include "CurvilinearGrid.hpp"

GridGeom::CurvilinearGridFromSplinesTransfinite::CurvilinearGridFromSplinesTransfinite(): m_splines(nullptr)
{
}

GridGeom::CurvilinearGridFromSplinesTransfinite::CurvilinearGridFromSplinesTransfinite(Splines* splines) : m_splines(splines)
{
    m_numN = 40;
    m_numM = 20;

};


bool GridGeom::CurvilinearGridFromSplinesTransfinite::Compute(CurvilinearGrid& curvilinearGrid)
{
    // compute the intersections
    if (m_numN == 0 || m_numM == 0)
    {
        return false;
    }

    // if the number of splines is less than 4 return false
    const auto numSplines = m_splines->m_numSplines;
    if (numSplines < 4)
    {
        return false;
    }

    // compute the intersections 
    bool successful = ComputeIntersections();
    if (!successful)
    {
        return false;
    }

    const int numMPoints = m_numM + 1;
    const int numNPoints = m_numN + 1;
    const int maxNumPoints = std::max(numMPoints, numNPoints);
    std::vector<double> distances;
    std::vector<double> adimensionalDistances;
    std::vector<double> intersectionDistances;
    std::vector<Point> points;

    // For each side of the plane reserve a vector
    std::vector<Point> sideOne;
    std::vector<Point> sideTwo;
    std::vector<Point> sideThree;
    std::vector<Point> sideFour;

    // The vector storing the transfinite interpolation results
    std::vector<std::vector<Point>> interpolationResult(numMPoints, std::vector<Point>(numNPoints));

    // allocate local vectors
    distances.reserve(maxNumPoints);
    adimensionalDistances.reserve(maxNumPoints);
    points.reserve(maxNumPoints);
    intersectionDistances.resize(numSplines);

    sideOne.reserve(maxNumPoints);
    sideTwo.reserve(maxNumPoints);
    sideThree.reserve(maxNumPoints);
    sideFour.reserve(maxNumPoints);

    // allocate the curvilinear grid
    curvilinearGrid.Set(m_numM, m_numN);

    int numMSplines = 0;
    int numNSplines = 0;
    for (int splineIndex = 0; splineIndex < numSplines; splineIndex++)
    {
        int numIntersections = 0;
        for (int i = 0; i < m_splineIntersectionRatios[splineIndex].size(); i++)
        {
            if (std::abs(m_splineIntersectionRatios[splineIndex][i]) > 0.0)
            {
                intersectionDistances[numIntersections] = m_splines->GetSplineLength(splineIndex, 0.0, m_splineIntersectionRatios[splineIndex][i], 10, false);
                numIntersections++;
            }
        }

        if (numIntersections < 0)
        {
            return false;
        }

        int numPoints;
        int position;
        int from;
        int to;
        if (splineIndex < m_firstMSplines)
        {
            numPoints = m_numM;
            position = (m_splineGroupIndexAndFromToIntersections[splineIndex][0] - 1) * m_numN;
            from = (m_splineGroupIndexAndFromToIntersections[splineIndex][1] - 1) * m_numM;
            to = (m_splineGroupIndexAndFromToIntersections[splineIndex][2] - 1) * m_numM + 1;
            numNSplines = std::max(numNSplines, m_splineGroupIndexAndFromToIntersections[splineIndex][0]);
        }
        else
        {
            numPoints = m_numN;
            position = (m_splineGroupIndexAndFromToIntersections[splineIndex][0] - 1) * m_numM;
            from = (m_splineGroupIndexAndFromToIntersections[splineIndex][1] - 1) * m_numN;
            to = (m_splineGroupIndexAndFromToIntersections[splineIndex][2] - 1) * m_numN + 1;
            numMSplines = std::max(numMSplines, m_splineGroupIndexAndFromToIntersections[splineIndex][0]);
        }

        distances.resize(numPoints + 1);
        adimensionalDistances.resize(numPoints + 1);
        points.resize(numPoints + 1);

        for (int i = 0; i < distances.size(); i++)
        {
            distances[i] = intersectionDistances[0] + (intersectionDistances[1] - intersectionDistances[0]) * i / numPoints;
        }

        m_splines->InterpolatePointsOnSpline(splineIndex,
            doubleMissingValue,
            false,
            distances,
            points,
            adimensionalDistances);

        // Start filling curvilinear grid
        int index = 0;
        for (int i = from; i < to; i++)
        {

            if (splineIndex < m_firstMSplines)
            {
                curvilinearGrid.m_grid[i][position] = points[index];
            }
            else
            {
                curvilinearGrid.m_grid[position][i] = points[index];
            }
            index++;
        }
    }

    sideOne.resize(numNPoints);
    sideTwo.resize(numNPoints);
    sideThree.resize(numMPoints);
    sideFour.resize(numMPoints);
    for (int i = 0; i < numMSplines - 1; i++)
    {
        for (int j = 0; j < numNSplines - 1; j++)
        {
            //Fill each block of the interpolation plane
            for (int k = 0; k < numMPoints; k++)
            {
                for (int l = 0; l < numNPoints; l++)
                {
                    const int m = i * m_numM + k;
                    const int n = j * m_numN + l;

                    // We are at the boundary 
                    if (!curvilinearGrid.m_grid[m][n].IsValid())
                    {
                        continue;
                    }

                    if (k == 0)
                    {
                        sideOne[l] = curvilinearGrid.m_grid[m][n];
                    }
                    if (k == m_numM)
                    {
                        sideTwo[l] = curvilinearGrid.m_grid[m][n];
                    }
                    if (l == 0)
                    {
                        sideThree[k] = curvilinearGrid.m_grid[m][n];
                    }
                    if (l == m_numN)
                    {
                        sideFour[k] = curvilinearGrid.m_grid[m][n];
                    }
                }
            }

            // call transfinite interpolation
            Interpolate(sideOne, sideTwo, sideThree, sideFour, interpolationResult);

            // assign the points
            for (int k = 0; k < numMPoints; k++)
            {
                for (int l = 0; l < numNPoints; l++)
                {
                    const int m = i * m_numM + k;
                    const int n = j * m_numN + l;

                    if (curvilinearGrid.m_grid[m][n].IsValid())
                    {
                        continue;
                    }

                    curvilinearGrid.m_grid[m][n] = interpolationResult[k][l];
                }
            }
        }
    }

    return true;
}


bool GridGeom::CurvilinearGridFromSplinesTransfinite::Interpolate(const std::vector<Point>& sideOne,
                                                                  const std::vector<Point>& sideTwo,
                                                                  const std::vector<Point>& sideThree,
                                                                  const std::vector<Point>& sideFour,
                                                                  std::vector<std::vector<Point>>& result)
{

    double totalLengthOne;
    std::vector<double> sideOneAdimensional(sideOne.size());
    ComputeAdimensionalDistancesFromPointSerie(sideOne, m_splines->m_projection, sideOneAdimensional, totalLengthOne);

    double totalLengthTwo;
    std::vector<double> sideTwoAdimensional(sideTwo.size());
    ComputeAdimensionalDistancesFromPointSerie(sideTwo, m_splines->m_projection, sideTwoAdimensional, totalLengthTwo);

    double totalLengthThree;
    std::vector<double> sideThreeAdimensional(sideThree.size());
    ComputeAdimensionalDistancesFromPointSerie(sideThree, m_splines->m_projection, sideThreeAdimensional, totalLengthThree);

    double totalLengthFour;
    std::vector<double> sideFourAdimensional(sideFour.size());
    ComputeAdimensionalDistancesFromPointSerie(sideFour, m_splines->m_projection, sideFourAdimensional, totalLengthFour);

    // now compute the adimensional distance of each point to be filled
    const int numMPoints = m_numM + 1;
    const int numNPoints = m_numN + 1;

    std::vector<std::vector<double>> iWeightFactor(numMPoints, std::vector<double>(numNPoints));
    std::vector<std::vector<double>> jWeightFactor(numMPoints, std::vector<double>(numNPoints));
    for (int i = 0; i < numMPoints; i++)
    {
        for (int j = 0; j < numNPoints; j++)
        {
            const double mWeight = i / m_numM;
            const double nWeight = j / m_numN;

            iWeightFactor[i][j] = (1.0 - nWeight) * sideThreeAdimensional[i] + nWeight * sideFourAdimensional[i];
            jWeightFactor[i][j] = (1.0 - mWeight) * sideOneAdimensional[j] + mWeight * sideTwoAdimensional[j];
        }
    }

    std::vector<std::vector<double>> weightOne(numMPoints, std::vector<double>(numNPoints));
    std::vector<std::vector<double>> weightTwo(numMPoints, std::vector<double>(numNPoints));
    std::vector<std::vector<double>> weightThree(numMPoints, std::vector<double>(numNPoints));
    std::vector<std::vector<double>> weightFour(numMPoints, std::vector<double>(numNPoints));
    for (int i = 0; i < numMPoints; i++)
    {
        for (int j = 0; j < numNPoints; j++)
        {

            weightOne[i][j] = (1.0 - jWeightFactor[i][j]) * totalLengthThree + jWeightFactor[i][j] * totalLengthFour;
            weightTwo[i][j] = (1.0 - iWeightFactor[i][j]) * totalLengthOne + iWeightFactor[i][j] * totalLengthTwo;
            weightThree[i][j] = weightTwo[i][j] / weightOne[i][j];
            weightFour[i][j] = weightOne[i][j] / weightTwo[i][j];
            const double wa = 1.0 / (weightThree[i][j] + weightFour[i][j]);
            weightOne[i][j] = wa * weightThree[i][j];
            weightTwo[i][j] = wa * weightFour[i][j];
        }
    }

    //border points
    result.resize(numMPoints, std::vector<Point>(numNPoints));
    for (int i = 0; i < numMPoints; i++)
    {
        result[i][0] = sideThree[i];
        result[i][m_numN] = sideFour[i];
    }
    for (int i = 0; i < numNPoints; i++)
    {
        result[0][i] = sideOne[i];
        result[m_numM][i] = sideTwo[i];
    }

    // first interpolation
    for (int i = 1; i < m_numM; i++)
    {
        for (int j = 1; j < m_numN; j++)
        {

            result[i][j].x = (sideOne[j].x * (1.0 - iWeightFactor[i][j]) + sideTwo[j].x * iWeightFactor[i][j]) * weightOne[i][j] +
                             (sideThree[i].x * (1.0 - jWeightFactor[i][j]) + sideFour[i].x * jWeightFactor[i][j]) * weightTwo[i][j];

            result[i][j].y = (sideOne[j].y * (1.0 - iWeightFactor[i][j]) + sideTwo[j].y * iWeightFactor[i][j]) * weightOne[i][j] +
                             (sideThree[i].y * (1.0 - jWeightFactor[i][j]) + sideFour[i].y * jWeightFactor[i][j]) * weightTwo[i][j];
        }
    }

    // update weights
    for (int i = 0; i < numMPoints; i++)
    {
        for (int j = 0; j < numNPoints; j++)
        {
            weightOne[i][j] = (1.0 - jWeightFactor[i][j]) * sideThreeAdimensional[i] * totalLengthThree +
                jWeightFactor[i][j] * sideFourAdimensional[i] * totalLengthFour;
            weightTwo[i][j] = (1.0 - iWeightFactor[i][j]) * sideOneAdimensional[j] * totalLengthOne +
                iWeightFactor[i][j] * sideTwoAdimensional[j] * totalLengthTwo;
        }
    }

    for (int i = 1; i < numMPoints; i++)
    {
        for (int j = 0; j < numNPoints; j++)
        {
            weightThree[i][j] = weightOne[i][j] - weightOne[i - 1][j];
        }
    }

    for (int i = 0; i < numMPoints; i++)
    {
        for (int j = 1; j < numNPoints; j++)
        {
            weightFour[i][j] = weightTwo[i][j] - weightTwo[i][j - 1];
        }
    }

    for (int i = 1; i < numMPoints; i++)
    {
        for (int j = 1; j < numNPoints - 1; j++)
        {
            weightOne[i][j] = 0.25 * (weightFour[i][j] + weightFour[i][j + 1] + weightFour[i - 1][j] + weightFour[i - 1][j + 1]) / weightThree[i][j];
        }
    }


    for (int i = 1; i < numMPoints - 1; i++)
    {
        for (int j = 1; j < numNPoints; j++)
        {
            weightTwo[i][j] = 0.25 * (weightThree[i][j] + weightThree[i][j - 1] + weightThree[i + 1][j] + weightThree[i + 1][j - 1]) / weightFour[i][j];
        }
    }

    // Iterate several times over
    const int numIterations = 25;
    for (int iter = 0; iter < numIterations; iter++)
    {
        // re-assign the weights
        for (int i = 0; i < numMPoints; i++)
        {
            for (int j = 0; j < numNPoints; j++)
            {
                weightThree[i][j] = result[i][j].x;
                weightFour[i][j] = result[i][j].y;
            }
        }


        for (int i = 1; i < m_numM; i++)
        {
            for (int j = 1; j < m_numN; j++)
            {

                const double wa = 1.0 / (weightOne[i][j] + weightOne[i + 1][j] + weightTwo[i][j] + weightTwo[i][j + 1]);

                result[i][j].x = wa * (weightThree[i][j] * weightOne[i][j] + weightThree[i][j] * weightOne[i + 1][j] +
                    weightThree[i][j] * weightTwo[i][j] + weightThree[i][j] * weightTwo[i][j + 1]);

                result[i][j].y = wa * (weightFour[i][j] * weightOne[i][j] + weightFour[i][j] * weightOne[i + 1][j] +
                    weightFour[i][j] * weightTwo[i][j] + weightFour[i][j] * weightTwo[i][j + 1]);
            }
        }
    }

    return true;
}


bool GridGeom::CurvilinearGridFromSplinesTransfinite::ComputeIntersections()
{
    const auto numSplines = m_splines->m_numSplines;

    // fill the splines with zeros
    m_splineType.resize(numSplines);
    std::fill(m_splineType.begin(), m_splineType.end(), 0);
    m_splineType[0] = 1;

    m_splineIntersectionRatios.resize(numSplines);
    std::fill(m_splineIntersectionRatios.begin(), m_splineIntersectionRatios.end(), std::vector<double>(numSplines, 0.0));

    for (int i = 0; i < numSplines; i++)
    {
        for (int j = i + 1; j < numSplines; j++)
        {
            double crossProductIntersection;
            Point intersectionPoint;
            double firstSplineRatio;
            double secondSplineRatio;

            // find intersections
            const auto areCrossing = m_splines->GetSplinesIntersection(i, j, crossProductIntersection, intersectionPoint, firstSplineRatio, secondSplineRatio);

            if (areCrossing) 
            {
                if (m_splineType[i] * m_splineType[j] == 1) 
                {
                    // intersecting twice
                    return false;
                }
                else if (m_splineType[i] == 0 && m_splineType[j] == 0)
                {
                    // both undefined
                }
                else if (m_splineType[j] == 0)
                {
                    m_splineType[j] = -m_splineType[i];
                    if (crossProductIntersection * m_splineType[i] < 0.0) 
                    {
                        // switch j
                        SwapVectorElements(m_splines->m_splineNodes[j], m_splines->m_numSplineNodes[j]);
                        secondSplineRatio = double(m_splines->m_numSplineNodes[j]) - 1.0 - secondSplineRatio;
                    }
                }
                else if (m_splineType[i] == 0)
                {
                    m_splineType[i] = -m_splineType[j];
                    if (crossProductIntersection * m_splineType[j] > 0.0)
                    {
                        // switch i
                        SwapVectorElements(m_splines->m_splineNodes[i], m_splines->m_numSplineNodes[i]);
                        firstSplineRatio = double(m_splines->m_numSplineNodes[i]) - 1.0 - firstSplineRatio;
                    }
                }
                m_splineIntersectionRatios[i][j] = firstSplineRatio;
                m_splineIntersectionRatios[j][i] = secondSplineRatio;
            }
        }
    }

    // find if one of the spline could not be classified
    // put the horizontal types on front first (done from swapping before)
    for (int i = 0; i < numSplines; i++)
    {
        if (m_splineType[i] == 0) 
        {
            return false;
        }
    }

    // find the first non m spline
    m_firstMSplines = FindIndex(m_splineType, -1);

    int maxExternalIterations = 10;
    for (int i = 0; i < maxExternalIterations; i++)
    {
        // sort along m 
        int maxInternalIterations = 100;
        for (int j = 0; j < maxInternalIterations; j++)
        {
            auto succeded = OrderSplines(0, m_firstMSplines, m_firstMSplines, numSplines);
            if (succeded)
            {
                break;
            }
        }

        // sort along n
        bool nSplineSortingHasNotChanged = true;
        for (int j = 0; j < maxInternalIterations; j++)
        {
            auto succeded = OrderSplines(m_firstMSplines, numSplines, 0, m_firstMSplines);
            if (succeded)
            {
                break;
            }
            else 
            {
                nSplineSortingHasNotChanged = false;
            }
        }

        if (nSplineSortingHasNotChanged) 
        {
            break;
        }
    }

    // Now determine the start and end spline corner points for each spline
    m_splineGroupIndexAndFromToIntersections.resize(numSplines, std::vector<int>(3, 0));

    // n direction
    for (int i = 0; i < m_firstMSplines; i++)
    {
        for (int j = m_firstMSplines; j < numSplines; j++)
        {
            int maxIndex = 0;
            int lastIndex = 0;
            for (int k = 0; k <= i; k++)
            {
                if (std::abs(m_splineIntersectionRatios[j][k]) > 0.0)
                {
                    maxIndex = m_splineGroupIndexAndFromToIntersections[lastIndex][0] + 1;
                    lastIndex = k;
                }
            }
            m_splineGroupIndexAndFromToIntersections[j][1] = maxIndex;
        }
        int maxIndex = 0;
        for (int j = m_firstMSplines; j < numSplines; j++)
        {
            if (std::abs(m_splineIntersectionRatios[j][i]) > 0.0)
            {
                maxIndex = std::max(maxIndex, m_splineGroupIndexAndFromToIntersections[j][1]);
            }
        }
        m_splineGroupIndexAndFromToIntersections[i][0] = maxIndex;
    }


    // m direction
    for (int i = m_firstMSplines; i < numSplines; i++)
    {
        for (int j = 0; j < m_firstMSplines; j++)
        {
            int maxIndex = 0;
            int lastIndex = m_firstMSplines;
            for (int k = m_firstMSplines; k <= i; k++)
            {
                if (std::abs(m_splineIntersectionRatios[j][k]) > 0.0)
                {
                    maxIndex = m_splineGroupIndexAndFromToIntersections[lastIndex][0] + 1;
                    lastIndex = k;
                }
            }
            m_splineGroupIndexAndFromToIntersections[j][2] = maxIndex;
        }
        int maxIndex = 0;
        for (int j = 0; j < m_firstMSplines; j++)
        {
            if (std::abs(m_splineIntersectionRatios[j][i]) > 0.0)
            {
                maxIndex = std::max(maxIndex, m_splineGroupIndexAndFromToIntersections[j][2]);
            }
        }
        m_splineGroupIndexAndFromToIntersections[i][0] = maxIndex;
    }

    for (int i = 0; i < numSplines; i++)
    {
        m_splineGroupIndexAndFromToIntersections[i][1] = 0;
        m_splineGroupIndexAndFromToIntersections[i][2] = 0;
    }

    // n constant, spline start end end
    for (int i = 0; i < m_firstMSplines; i++)
    {
        for (int j = m_firstMSplines; j < numSplines; j++)
        {
            if (std::abs(m_splineIntersectionRatios[i][j]) > 0.0)
            {
                if (m_splineGroupIndexAndFromToIntersections[i][1] == 0)
                {
                    m_splineGroupIndexAndFromToIntersections[i][1] = m_splineGroupIndexAndFromToIntersections[j][0];
                }
                m_splineGroupIndexAndFromToIntersections[i][2] = m_splineGroupIndexAndFromToIntersections[j][0];
            }
        }
    }

    // m constant, spline start end end
    for (int i = m_firstMSplines; i < numSplines; i++)
    {
        for (int j = 0; j < m_firstMSplines; j++)
        {
            if (std::abs(m_splineIntersectionRatios[i][j]) > 0.0)
            {
                if (m_splineGroupIndexAndFromToIntersections[i][1] == 0)
                {
                    m_splineGroupIndexAndFromToIntersections[i][1] = m_splineGroupIndexAndFromToIntersections[j][0];
                }
                m_splineGroupIndexAndFromToIntersections[i][2] = m_splineGroupIndexAndFromToIntersections[j][0];
            }
        }
    }

    return true;
}

bool GridGeom::CurvilinearGridFromSplinesTransfinite::OrderSplines(int startFirst,
    int endFirst,
    int startSecond,
    int endSecond)
{
    int numSwap = 0;
    const auto numSplines = m_splines->m_numSplines;

    for (int i = startFirst; i < endFirst; i++)
    {
        for (int j = startSecond; j < endSecond; j++)
        {
            const auto firstIntersectionRatio = m_splineIntersectionRatios[i][j];
            if (firstIntersectionRatio == 0)
            {
                continue;
            }

            for (int k = j + 1; k < endSecond; k++)
            {
                const auto secondIntersectionRatio = m_splineIntersectionRatios[i][k];

                // all fine nothing to do, they are already sorted
                if (secondIntersectionRatio == 0 || firstIntersectionRatio <= secondIntersectionRatio)
                {
                    continue;
                }
                //they must be swapped
                SwapRows(m_splines->m_splineNodes, j, k);
                SwapRows(m_splineIntersectionRatios, j, k);
                SwapColumns(m_splineIntersectionRatios, j, k);


                //repeat the entire procedure once more
                return false;
            }
        }
    }

    return true;
}

template<typename T>
bool  GridGeom::CurvilinearGridFromSplinesTransfinite::SwapRows(std::vector<std::vector<T>>& v, int firstRow, int secondRow)
{
    int minSize = std::min(v[firstRow].size(), v[secondRow].size());
    minSize = std::min(minSize, m_splines->m_numSplines);

    for (int i = 0; i < minSize; i++)
    {
        std::swap(v[firstRow][i], v[secondRow][i]);
    }
    return true;
}

template<typename T>
bool  GridGeom::CurvilinearGridFromSplinesTransfinite::SwapColumns(std::vector<std::vector<T>>& v, int firstColumn, int secondColumn)
{
    for (int i = 0; i < m_splines->m_numSplines; i++)
    {
        if (firstColumn >= v[i].size() || secondColumn >= v[i].size())
        {
            continue;
        }

        std::swap(v[i][firstColumn], v[i][secondColumn]);
    }
    return true;
}