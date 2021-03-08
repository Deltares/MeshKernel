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
#include <memory>
#include <vector>

#include <MeshKernel/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGridFromSplinesTransfinite.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Splines.hpp>

meshkernel::CurvilinearGridFromSplinesTransfinite::CurvilinearGridFromSplinesTransfinite(std::shared_ptr<Splines> splines,
                                                                                         const meshkernelapi::CurvilinearParameters& curvilinearParameters) : m_splines(splines),
                                                                                                                                                              m_numM(curvilinearParameters.MRefinement),
                                                                                                                                                              m_numN(curvilinearParameters.NRefinement)

{
}

meshkernel::CurvilinearGrid meshkernel::CurvilinearGridFromSplinesTransfinite::Compute()
{
    if (m_numN == 0 || m_numM == 0)
    {
        throw std::invalid_argument("CurvilinearGridFromSplinesTransfinite::Compute: There are no rows or no columns.");
    }

    const auto numSplines = m_splines->GetNumSplines();
    if (numSplines < 4)
    {
        throw std::invalid_argument("CurvilinearGridFromSplinesTransfinite::Compute: The number of splines is less than four.");
    }

    ComputeIntersections();

    const auto numMPoints = m_numM + 1;
    const auto numNPoints = m_numN + 1;
    const auto maxNumPoints = std::max(numMPoints, numNPoints);
    std::vector<double> distances;
    std::vector<double> adimensionalDistances;
    std::vector<double> intersectionDistances;
    std::vector<Point> points;

    // For each side of the plane reserve a vector
    std::vector<Point> sideOne;
    std::vector<Point> sideTwo;
    std::vector<Point> sideThree;
    std::vector<Point> sideFour;

    // allocate local vectors
    distances.reserve(maxNumPoints);
    adimensionalDistances.reserve(maxNumPoints);
    points.reserve(maxNumPoints);
    intersectionDistances.resize(numSplines);

    sideOne.reserve(maxNumPoints);
    sideTwo.reserve(maxNumPoints);
    sideThree.reserve(maxNumPoints);
    sideFour.reserve(maxNumPoints);

    // Allocate the curvilinear grid. We can have multiple divisions along N and M.
    const auto TotalMColumns = (m_numNSplines - 1) * m_numM;
    const auto TotalNRows = (m_numMSplines - 1) * m_numN;
    std::vector<std::vector<Point>> gridNodes(TotalMColumns + 1, std::vector<Point>(TotalNRows + 1));

    size_t numMSplines = 0;
    size_t numNSplines = 0;
    for (auto splineIndex = 0; splineIndex < numSplines; splineIndex++)
    {
        size_t numIntersections = 0;

        for (const auto& value : m_splineIntersectionRatios[splineIndex])
        {
            if (std::abs(value) > 0.0)
            {
                intersectionDistances[numIntersections] = m_splines->ComputeSplineLength(splineIndex, 0.0, value, 10, false);
                numIntersections++;
            }
        }

        if (numIntersections < 2)
        {
            throw std::invalid_argument("CurvilinearGridFromSplinesTransfinite::Compute: The number of intersections are less than two.");
        }

        size_t numPoints;
        size_t numDiscretizations;
        size_t position;
        size_t from;
        size_t to;
        if (splineIndex < m_numMSplines)
        {
            numPoints = TotalMColumns + 1;
            numDiscretizations = m_numM;
            position = (m_splineGroupIndexAndFromToIntersections[splineIndex][0] - 1) * m_numN;
            from = (m_splineGroupIndexAndFromToIntersections[splineIndex][1] - 1) * m_numM;
            to = (m_splineGroupIndexAndFromToIntersections[splineIndex][2] - 1) * m_numM + 1;
            numNSplines = std::max(numNSplines, m_splineGroupIndexAndFromToIntersections[splineIndex][0]);
        }
        else
        {
            numPoints = TotalNRows + 1;
            numDiscretizations = m_numN;
            position = (m_splineGroupIndexAndFromToIntersections[splineIndex][0] - 1) * m_numM;
            from = (m_splineGroupIndexAndFromToIntersections[splineIndex][1] - 1) * m_numN;
            to = (m_splineGroupIndexAndFromToIntersections[splineIndex][2] - 1) * m_numN + 1;
            numMSplines = std::max(numMSplines, m_splineGroupIndexAndFromToIntersections[splineIndex][0]);
        }

        distances.resize(numPoints);
        adimensionalDistances.resize(numPoints);
        points.resize(numPoints);

        ComputeDiscretizations(numIntersections,
                               numPoints,
                               numDiscretizations,
                               intersectionDistances,
                               distances);

        const auto [points, adimensionalDistances] = m_splines->InterpolatePointsOnSpline(splineIndex,
                                                                                          doubleMissingValue,
                                                                                          false,
                                                                                          distances);

        // Start filling curvilinear grid
        size_t index = 0;
        for (auto i = from; i < to; i++)
        {
            if (splineIndex < m_numMSplines)
            {
                gridNodes[i][position] = points[index];
            }
            else
            {
                gridNodes[position][i] = points[index];
            }
            index++;
        }
    }

    sideOne.resize(numNPoints);
    sideTwo.resize(numNPoints);
    sideThree.resize(numMPoints);
    sideFour.resize(numMPoints);
    for (auto i = 0; i < numMSplines - 1; i++)
    {
        for (auto j = 0; j < numNSplines - 1; j++)
        {
            //Fill each block of the interpolation plane
            for (auto k = 0; k < numMPoints; k++)
            {
                for (auto l = 0; l < numNPoints; l++)
                {
                    const auto m = i * m_numM + k;
                    const auto n = j * m_numN + l;

                    // We are at the boundary
                    if (!gridNodes[m][n].IsValid())
                    {
                        continue;
                    }

                    if (k == 0)
                    {
                        sideOne[l] = gridNodes[m][n];
                    }
                    if (k == m_numM)
                    {
                        sideTwo[l] = gridNodes[m][n];
                    }
                    if (l == 0)
                    {
                        sideThree[k] = gridNodes[m][n];
                    }
                    if (l == m_numN)
                    {
                        sideFour[k] = gridNodes[m][n];
                    }
                }
            }

            // call transfinite interpolation
            const auto interpolationResult = DiscretizeTransfinite(sideOne,
                                                                   sideTwo,
                                                                   sideThree,
                                                                   sideFour,
                                                                   m_splines->m_projection,
                                                                   m_numM,
                                                                   m_numN);

            // assign the points
            for (auto k = 0; k < numMPoints; k++)
            {
                for (auto l = 0; l < numNPoints; l++)
                {
                    const auto m = i * m_numM + k;
                    const auto n = j * m_numN + l;

                    if (gridNodes[m][n].IsValid())
                    {
                        continue;
                    }

                    gridNodes[m][n] = interpolationResult[k][l];
                }
            }
        }
    }

    return CurvilinearGrid(std::move(gridNodes), m_splines->m_projection);
}

void meshkernel::CurvilinearGridFromSplinesTransfinite::ComputeDiscretizations(size_t numIntersections,
                                                                               size_t numPoints,
                                                                               size_t numDiscretizations,
                                                                               const std::vector<double>& intersectionDistances,
                                                                               std::vector<double>& distances) const
{

    if (numIntersections == 2)
    {
        for (auto i = 0; i < numPoints; i++)
        {
            distances[i] = intersectionDistances[0] + (intersectionDistances[1] - intersectionDistances[0]) * i / static_cast<double>(numDiscretizations);
        }
    }
    if (numIntersections > 2)
    {
        std::vector<double> ratioSegments(numIntersections, 0.0);
        for (auto i = 1; i < numIntersections - 1; i++)
        {
            ratioSegments[i] = (intersectionDistances[i + 1] - intersectionDistances[i]) / (intersectionDistances[i] - intersectionDistances[i - 1]);
        }
        ratioSegments[0] = ratioSegments[1];
        ratioSegments[numIntersections - 1] = ratioSegments[numIntersections - 2];

        // the ratios
        std::vector<double> leftDiscretization(numDiscretizations + 1, 0.0);
        std::vector<double> rightDiscretization(numDiscretizations + 1, 0.0);
        for (auto i = 0; i < numIntersections - 1; i++)
        {
            const double rightRatio = std::pow(ratioSegments[i + 1], 1.0 / static_cast<double>(numDiscretizations));
            ComputeExponentialDistances(rightRatio,
                                        intersectionDistances[i],
                                        intersectionDistances[i + 1],
                                        rightDiscretization);

            const double leftRatio = std::pow(ratioSegments[i], 1.0 / static_cast<double>(numDiscretizations));
            ComputeExponentialDistances(leftRatio,
                                        intersectionDistances[i],
                                        intersectionDistances[i + 1],
                                        leftDiscretization);

            for (auto j = 0; j < numDiscretizations + 1; j++)
            {

                double ar = double(j) / static_cast<double>(numDiscretizations);
                double al = 1.0 - ar;

                const auto index = i * numDiscretizations + j;
                distances[index] = ar * rightDiscretization[j] + al * leftDiscretization[j];

                //adjust a second time
                ar = (distances[index] - intersectionDistances[i]) / (intersectionDistances[i + 1] - intersectionDistances[i]);
                al = 1.0 - ar;
                distances[index] = ar * rightDiscretization[j] + al * leftDiscretization[j];
            }
        }
    }
}

void meshkernel::CurvilinearGridFromSplinesTransfinite::ComputeExponentialDistances(double factor,
                                                                                    double leftDistance,
                                                                                    double rightDistance,
                                                                                    std::vector<double>& distances) const
{
    distances[0] = 0.0;
    double incrementRatio = 1.0;
    for (auto i = 0; i < distances.size() - 1; i++)
    {
        distances[i + 1] = distances[i] + incrementRatio;
        incrementRatio = incrementRatio * factor;
    }

    incrementRatio = (rightDistance - leftDistance) / distances.back();

    for (auto& value : distances)
    {
        value = leftDistance + incrementRatio * value;
    }
}

void meshkernel::CurvilinearGridFromSplinesTransfinite::ComputeIntersections()
{
    const auto numSplines = m_splines->GetNumSplines();

    // fill the splines with zeros
    m_splineType.resize(numSplines);
    std::fill(m_splineType.begin(), m_splineType.end(), 0);
    m_splineType[0] = 1;

    m_splineIntersectionRatios.resize(numSplines);
    std::fill(m_splineIntersectionRatios.begin(), m_splineIntersectionRatios.end(), std::vector<double>(numSplines, 0.0));

    for (auto i = 0; i < numSplines; i++)
    {
        for (auto j = i + 1; j < numSplines; j++)
        {
            double crossProductIntersection;
            Point intersectionPoint;
            double firstSplineRatio;
            double secondSplineRatio;
            const auto numNodesISpline = static_cast<int>(m_splines->m_splineNodes[i].size());
            const auto numNodesJSpline = static_cast<int>(m_splines->m_splineNodes[j].size());

            // find intersections
            const auto areCrossing = m_splines->GetSplinesIntersection(i, j, crossProductIntersection, intersectionPoint, firstSplineRatio, secondSplineRatio);

            if (areCrossing)
            {
                if (m_splineType[i] * m_splineType[j] == 1)
                {
                    throw std::invalid_argument("CurvilinearGridFromSplinesTransfinite::Compute: At least two splines are intersecting twice.");
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
                        SwapVectorElements(m_splines->m_splineNodes[j]);
                        secondSplineRatio = static_cast<double>(numNodesJSpline) - 1.0 - secondSplineRatio;
                    }
                }
                else if (m_splineType[i] == 0)
                {
                    m_splineType[i] = -m_splineType[j];
                    if (crossProductIntersection * m_splineType[j] > 0.0)
                    {
                        // switch i
                        SwapVectorElements(m_splines->m_splineNodes[i]);
                        firstSplineRatio = static_cast<double>(numNodesISpline) - 1.0 - firstSplineRatio;
                    }
                }
                m_splineIntersectionRatios[i][j] = firstSplineRatio;
                m_splineIntersectionRatios[j][i] = secondSplineRatio;
            }
        }
    }

    for (auto i = 0; i < numSplines; i++)
    {
        if (m_splineType[i] == 0)
        {
            throw std::invalid_argument("CurvilinearGridFromSplinesTransfinite::Compute: At least one of the splines could not be classified.");
        }
    }

    // find the first non m spline
    m_numMSplines = FindIndex(m_splineType, -1);
    m_numNSplines = numSplines - m_numMSplines;

    const size_t maxExternalIterations = 10;
    for (auto i = 0; i < maxExternalIterations; i++)
    {
        // sort along m
        const size_t maxInternalIterations = 100;
        for (auto j = 0; j < maxInternalIterations; j++)
        {
            const auto successful = OrderSplines(0, m_numMSplines, m_numMSplines, numSplines);
            if (successful)
            {
                break;
            }
        }

        // sort along n
        bool nSplineSortingHasNotChanged = true;
        for (auto j = 0; j < maxInternalIterations; j++)
        {
            const auto successful = OrderSplines(m_numMSplines, numSplines, 0, m_numMSplines);
            if (successful)
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
    m_splineGroupIndexAndFromToIntersections.resize(numSplines, std::vector<size_t>(3, 0));

    // n direction
    for (auto i = 0; i < m_numMSplines; i++)
    {
        for (auto j = m_numMSplines; j < numSplines; j++)
        {
            size_t maxIndex = 0;
            size_t lastIndex = 0;
            for (auto k = 0; k <= i; k++)
            {
                if (std::abs(m_splineIntersectionRatios[j][k]) > 0.0)
                {
                    maxIndex = m_splineGroupIndexAndFromToIntersections[lastIndex][0] + 1;
                    lastIndex = k;
                }
            }
            m_splineGroupIndexAndFromToIntersections[j][1] = maxIndex;
        }
        size_t maxIndex = 0;
        for (auto j = m_numMSplines; j < numSplines; j++)
        {
            if (std::abs(m_splineIntersectionRatios[j][i]) > 0.0)
            {
                maxIndex = std::max(maxIndex, m_splineGroupIndexAndFromToIntersections[j][1]);
            }
        }
        m_splineGroupIndexAndFromToIntersections[i][0] = maxIndex;
    }

    // m direction
    for (auto i = m_numMSplines; i < numSplines; i++)
    {
        for (auto j = 0; j < m_numMSplines; j++)
        {
            size_t maxIndex = 0;
            size_t lastIndex = m_numMSplines;
            for (auto k = m_numMSplines; k <= i; k++)
            {
                if (std::abs(m_splineIntersectionRatios[j][k]) > 0.0)
                {
                    maxIndex = m_splineGroupIndexAndFromToIntersections[lastIndex][0] + 1;
                    lastIndex = k;
                }
            }
            m_splineGroupIndexAndFromToIntersections[j][2] = maxIndex;
        }
        size_t maxIndex = 0;
        for (auto j = 0; j < m_numMSplines; j++)
        {
            if (std::abs(m_splineIntersectionRatios[j][i]) > 0.0)
            {
                maxIndex = std::max(maxIndex, m_splineGroupIndexAndFromToIntersections[j][2]);
            }
        }
        m_splineGroupIndexAndFromToIntersections[i][0] = maxIndex;
    }

    for (auto i = 0; i < numSplines; i++)
    {
        m_splineGroupIndexAndFromToIntersections[i][1] = 0;
        m_splineGroupIndexAndFromToIntersections[i][2] = 0;
    }

    // n constant, spline start end end
    for (auto i = 0; i < m_numMSplines; i++)
    {
        for (auto j = m_numMSplines; j < numSplines; j++)
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
    for (auto i = m_numMSplines; i < numSplines; i++)
    {
        for (auto j = 0; j < m_numMSplines; j++)
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
}

bool meshkernel::CurvilinearGridFromSplinesTransfinite::OrderSplines(size_t startFirst,
                                                                     size_t endFirst,
                                                                     size_t startSecond,
                                                                     size_t endSecond)
{
    for (auto i = startFirst; i < endFirst; i++)
    {
        for (auto j = startSecond; j < endSecond; j++)
        {
            const auto firstIntersectionRatio = m_splineIntersectionRatios[i][j];
            if (IsEqual(firstIntersectionRatio, 0.0))
            {
                continue;
            }

            for (auto k = j + 1; k < endSecond; k++)
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

template <typename T>
void meshkernel::CurvilinearGridFromSplinesTransfinite::SwapRows(std::vector<std::vector<T>>& v, size_t firstRow, size_t secondRow) const
{
    auto minSize = std::min(v[firstRow].size(), v[secondRow].size());
    minSize = std::min(minSize, m_splines->GetNumSplines());

    for (auto i = 0; i < minSize; i++)
    {
        std::swap(v[firstRow][i], v[secondRow][i]);
    }
}

template <typename T>
void meshkernel::CurvilinearGridFromSplinesTransfinite::SwapColumns(std::vector<std::vector<T>>& v, size_t firstColumn, size_t secondColumn) const
{
    for (auto i = 0; i < m_splines->GetNumSplines(); i++)
    {
        if (firstColumn >= v[i].size() || secondColumn >= v[i].size())
        {
            continue;
        }

        std::swap(v[i][firstColumn], v[i][secondColumn]);
    }
}
