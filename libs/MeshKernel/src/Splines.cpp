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

#include <iostream>

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/LandBoundaries.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/SplineAlgorithms.hpp>
#include <MeshKernel/Splines.hpp>

using meshkernel::Splines;

Splines::Splines(Projection projection) : m_projection(projection) {}

Splines::Splines(CurvilinearGrid const& grid)
{
    // first the m_n m_m-gridlines
    lin_alg::Matrix<Point> mGridLines(grid.m_numN, grid.m_numM);
    for (UInt n = 0; n < grid.m_numN; ++n)
    {
        for (UInt m = 0; m < grid.m_numM; ++m)
        {
            mGridLines(n, m) = grid.m_gridNodes(m, n);
        }
        AddSpline(lin_alg::MatrixRowToSTLVector(mGridLines, n));
    }

    // then the m_m m_n-gridlines
    for (UInt m = 0; m < grid.m_numM; ++m)
    {
        AddSpline(lin_alg::MatrixRowToSTLVector(grid.m_gridNodes, m));
    }

    m_projection = grid.m_projection;
}

/// add a new spline, return the index
void Splines::AddSpline(const std::vector<Point>& splines, UInt start, UInt size)
{
    if (size == 0)
    {
        return;
    }

    // copy the spline nodes from start to start + size
    UInt count = 0;
    std::vector<Point> splinesNodes(size);
    for (auto i = start; i < start + size; ++i)
    {
        splinesNodes[count] = splines[i];
        count++;
    }
    m_splineNodes.emplace_back(splinesNodes);

    // compute second order derivatives
    std::vector<Point> splineDerivatives(splinesNodes.size());
    const auto indices = FindIndices(splinesNodes, 0, static_cast<UInt>(splinesNodes.size()), constants::missing::doubleValue);
    for (auto index : indices)
    {
        const auto& [startIndex, endIndex] = index;
        const auto derivatives = SplineAlgorithms::SecondOrderDerivative(splinesNodes, startIndex, endIndex);
        for (auto j = startIndex; j <= endIndex; ++j)
        {
            splineDerivatives[j] = derivatives[j - startIndex];
        }
    }
    m_splineDerivatives.emplace_back(splineDerivatives);

    m_splinesLength.emplace_back(ComputeSplineLength(GetNumSplines() - 1, 0.0, static_cast<double>(size - 1)));
}

void Splines::AddSpline(const std::vector<Point>& splines)
{
    AddSpline(splines, 0, static_cast<UInt>(splines.size()));
}

void Splines::DeleteSpline(UInt splineIndex)
{
    m_splineNodes.erase(m_splineNodes.begin() + splineIndex);
    m_splineDerivatives.erase(m_splineDerivatives.begin() + splineIndex);
    m_splinesLength.erase(m_splinesLength.begin() + splineIndex);
}

void Splines::AddPointInExistingSpline(UInt splineIndex, const Point& point)
{
    if (splineIndex > GetNumSplines())
    {
        throw std::invalid_argument("Splines::AddPointInExistingSpline: Invalid spline index.");
    }
    m_splineNodes[splineIndex].emplace_back(point);
}

bool Splines::GetSplinesIntersection(UInt first,
                                     UInt second,
                                     double& crossProductIntersection,
                                     Point& intersectionPoint,
                                     double& firstSplineRatio,
                                     double& secondSplineRatio)
{
    double minimumCrossingDistance = std::numeric_limits<double>::max();
    double crossingDistance;
    UInt numCrossing = 0;
    double firstCrossingRatio = -1.0;
    double secondCrossingRatio = -1.0;
    UInt firstCrossingIndex = 0;
    UInt secondCrossingIndex = 0;
    Point closestIntersection;
    const auto numNodesFirstSpline = static_cast<UInt>(m_splineNodes[first].size());
    const auto numNodesSecondSpline = static_cast<UInt>(m_splineNodes[second].size());

    // First find a valid crossing, the closest to spline central point
    for (UInt n = 0; n < numNodesFirstSpline - 1; n++)
    {
        for (UInt nn = 0; nn < numNodesSecondSpline - 1; nn++)
        {
            Point intersection;
            double crossProduct;
            double firstRatio;
            double secondRatio;
            const bool areCrossing = AreSegmentsCrossing(m_splineNodes[first][n],
                                                         m_splineNodes[first][n + 1],
                                                         m_splineNodes[second][nn],
                                                         m_splineNodes[second][nn + 1],
                                                         false,
                                                         m_projection,
                                                         intersection,
                                                         crossProduct,
                                                         firstRatio,
                                                         secondRatio);

            if (areCrossing)
            {
                if (numNodesFirstSpline == 2)
                {
                    crossingDistance = std::min(minimumCrossingDistance, std::abs(firstRatio - 0.5));
                }
                else if (numNodesSecondSpline == 2)
                {
                    crossingDistance = std::abs(secondRatio - 0.5);
                }
                else
                {
                    crossingDistance = minimumCrossingDistance;
                }

                if (crossingDistance < minimumCrossingDistance || numCrossing == 0)
                {
                    minimumCrossingDistance = crossingDistance;
                    numCrossing = 1;
                    firstCrossingIndex = n;            // TI0
                    secondCrossingIndex = nn;          // TJ0
                    firstCrossingRatio = firstRatio;   // SL
                    secondCrossingRatio = secondRatio; // SM
                }
            }
            closestIntersection = intersection;
        }
    }

    // if no crossing found, return
    if (numCrossing == 0)
    {
        return false;
    }

    double firstCrossing = IsEqual(firstCrossingRatio, -1.0) ? 0.0 : static_cast<double>(firstCrossingIndex) + firstCrossingRatio;
    double secondCrossing = IsEqual(secondCrossingRatio, -1.0) ? 0.0 : static_cast<double>(secondCrossingIndex) + secondCrossingRatio;

    // use bisection to find the intersection
    double squaredDistanceBetweenCrossings = std::numeric_limits<double>::max();
    const double maxSquaredDistanceBetweenCrossings = 1e-12;
    const double maxDistanceBetweenNodes = 0.0001;
    double firstRatioIterations = 1.0;
    double secondRatioIterations = 1.0;
    UInt numIterations = 0;
    while (squaredDistanceBetweenCrossings > maxSquaredDistanceBetweenCrossings && numIterations < 20)
    {
        // increment counter
        numIterations++;

        if (firstCrossingRatio > 0 && firstCrossingRatio < 1.0)
        {
            firstRatioIterations = 0.5 * firstRatioIterations;
        }
        if (secondCrossingRatio > 0 && secondCrossingRatio < 1.0)
        {
            secondRatioIterations = 0.5 * secondRatioIterations;
        }

        firstCrossing = std::max(0.0, std::min(firstCrossing, double(numNodesFirstSpline)));
        secondCrossing = std::max(0.0, std::min(secondCrossing, double(numNodesSecondSpline)));

        const double firstLeft = std::max(0.0, std::min(double(numNodesFirstSpline - 1), firstCrossing - firstRatioIterations / 2.0));
        const double firstRight = std::max(0.0, std::min(double(numNodesFirstSpline - 1), firstCrossing + firstRatioIterations / 2.0));

        const double secondLeft = std::max(0.0, std::min(double(numNodesSecondSpline - 1), secondCrossing - secondRatioIterations / 2.0));
        const double secondRight = std::max(0.0, std::min(double(numNodesSecondSpline - 1), secondCrossing + secondRatioIterations / 2.0));

        firstRatioIterations = firstRight - firstLeft;
        secondRatioIterations = secondRight - secondLeft;

        const auto firstLeftSplinePoint = ComputePointOnSplineAtAdimensionalDistance(m_splineNodes[first], m_splineDerivatives[first], firstLeft);
        if (!firstLeftSplinePoint.IsValid())
        {
            return false;
        }

        const auto firstRightSplinePoint = ComputePointOnSplineAtAdimensionalDistance(m_splineNodes[first], m_splineDerivatives[first], firstRight);
        if (!firstRightSplinePoint.IsValid())
        {
            return false;
        }

        const auto secondLeftSplinePoint = ComputePointOnSplineAtAdimensionalDistance(m_splineNodes[second], m_splineDerivatives[second], secondLeft);
        if (!secondLeftSplinePoint.IsValid())
        {
            return false;
        }

        const auto secondRightSplinePoint = ComputePointOnSplineAtAdimensionalDistance(m_splineNodes[second], m_splineDerivatives[second], secondRight);
        if (!secondRightSplinePoint.IsValid())
        {
            return false;
        }

        Point oldIntersection = closestIntersection;

        double crossProduct;
        double firstRatio = constants::missing::doubleValue;
        double secondRatio = constants::missing::doubleValue;
        const bool areCrossing = AreSegmentsCrossing(firstLeftSplinePoint,
                                                     firstRightSplinePoint,
                                                     secondLeftSplinePoint,
                                                     secondRightSplinePoint,
                                                     true,
                                                     m_projection,
                                                     closestIntersection,
                                                     crossProduct,
                                                     firstRatio,
                                                     secondRatio);

        // search close by
        if (firstRatio > -2.0 && firstRatio < 3.0 && secondRatio > -2.0 && secondRatio < 3.0)
        {
            const double previousFirstCrossing = firstCrossing;
            const double previousSecondCrossing = secondCrossing;

            firstCrossing = firstLeft + firstRatio * (firstRight - firstLeft);
            secondCrossing = secondLeft + secondRatio * (secondRight - secondLeft);

            firstCrossing = std::max(0.0, std::min(static_cast<double>(numNodesFirstSpline) - 1.0, firstCrossing));
            secondCrossing = std::max(0.0, std::min(static_cast<double>(numNodesSecondSpline) - 1.0, secondCrossing));

            if (areCrossing)
            {
                numCrossing = 1;
                crossProductIntersection = crossProduct;
            }

            if (std::abs(firstCrossing - previousFirstCrossing) > maxDistanceBetweenNodes ||
                std::abs(secondCrossing - previousSecondCrossing) > maxDistanceBetweenNodes)
            {
                squaredDistanceBetweenCrossings = ComputeSquaredDistance(oldIntersection, closestIntersection, m_projection);
            }
            else
            {
                break;
            }
        }
    }

    if (numCrossing == 1)
    {
        intersectionPoint = closestIntersection;
        firstSplineRatio = firstCrossing;
        secondSplineRatio = secondCrossing;
        return true;
    }

    // not crossing
    return false;
}

double Splines::ComputeSplineLength(UInt index,
                                    double startAdimensionalCoordinate,
                                    double endAdimensionalCoordinate,
                                    UInt numSamples,
                                    bool accountForCurvature,
                                    double height,
                                    double assignedDelta) const
{
    if (m_splineNodes[index].empty())
    {
        return 0.0;
    }

    double delta = assignedDelta;
    UInt numPoints = static_cast<UInt>(endAdimensionalCoordinate / delta) + 1;
    if (delta < 0.0)
    {
        delta = 1.0 / static_cast<double>(numSamples);
        // TODO: Refactor or at least document the calculation of "numPoints"
        numPoints = static_cast<UInt>(std::max(std::floor(0.9999 + (endAdimensionalCoordinate - startAdimensionalCoordinate) / delta), 10.0));
        delta = (endAdimensionalCoordinate - startAdimensionalCoordinate) / static_cast<double>(numPoints);
    }

    // first point
    auto leftPoint = ComputePointOnSplineAtAdimensionalDistance(m_splineNodes[index], m_splineDerivatives[index], startAdimensionalCoordinate);
    double splineLength = 0.0;

    auto rightPointCoordinateOnSpline = startAdimensionalCoordinate;
    for (UInt p = 0; p < numPoints; ++p)
    {
        const double leftPointCoordinateOnSpline = rightPointCoordinateOnSpline;
        rightPointCoordinateOnSpline += delta;
        if (rightPointCoordinateOnSpline > endAdimensionalCoordinate)
        {
            rightPointCoordinateOnSpline = endAdimensionalCoordinate;
        }

        const auto rightPoint = ComputePointOnSplineAtAdimensionalDistance(m_splineNodes[index], m_splineDerivatives[index], rightPointCoordinateOnSpline);
        if (!rightPoint.IsValid())
        {
            continue;
        }

        double curvatureFactor = 0.0;
        if (accountForCurvature)
        {
            const auto [normalVector, tangentialVector, computedCurvatureFactor] = ComputeCurvatureOnSplinePoint(index, 0.5 * (rightPointCoordinateOnSpline + leftPointCoordinateOnSpline));
            curvatureFactor = computedCurvatureFactor;
        }
        splineLength = splineLength + ComputeDistance(leftPoint, rightPoint, m_projection) * (1.0 + curvatureFactor * height);
        leftPoint = rightPoint;
    }

    return splineLength;
}

std::tuple<meshkernel::Point, meshkernel::Point, double>
Splines::ComputeCurvatureOnSplinePoint(UInt splineIndex, double adimensionalPointCoordinate) const
{
    return SplineAlgorithms::ComputeCurvatureOnSplinePoint(m_splineNodes[splineIndex], m_splineDerivatives[splineIndex], adimensionalPointCoordinate, m_projection);
}

std::tuple<std::vector<meshkernel::Point>, std::vector<double>>
Splines::ComputePointOnSplineFromAdimensionalDistance(UInt index,
                                                      double maximumGridHeight,
                                                      bool isSpacingCurvatureAdapted,
                                                      const std::vector<double>& distances)
{

    std::vector<Point> points(distances.size());
    std::vector<double> adimensionalDistances(distances.size());

    FuncAdimensionalToDimensionalDistanceOnSpline func(this, index, isSpacingCurvatureAdapted, maximumGridHeight);
    const auto numNodes = static_cast<UInt>(m_splineNodes[index].size());
    for (UInt i = 0, size = static_cast<UInt>(distances.size()); i < size; ++i)
    {
        func.SetDimensionalDistance(distances[i]);
        adimensionalDistances[i] = FindFunctionRootWithGoldenSectionSearch(func, 0, static_cast<double>(numNodes) - 1.0);
        points[i] = ComputePointOnSplineAtAdimensionalDistance(m_splineNodes[index], m_splineDerivatives[index], adimensionalDistances[i]);
        if (!points[i].IsValid())
        {
            throw AlgorithmError("Could not interpolate spline points.");
        }
    }
    return {points, adimensionalDistances};
}

meshkernel::Point Splines::ComputeClosestPointOnSplineSegment(UInt index, double startSplineSegment, double endSplineSegment, Point point)
{
    FuncDistanceFromAPoint func(this, index, point);
    const auto adimensionalDistance = FindFunctionRootWithGoldenSectionSearch(func, startSplineSegment, endSplineSegment);
    return ComputePointOnSplineAtAdimensionalDistance(m_splineNodes[index], m_splineDerivatives[index], adimensionalDistance);
}

void Splines::SnapSpline(const size_t splineIndex,
                         const LandBoundary& landBoundary,
                         const int numberOfIterations)
{
    if (splineIndex >= GetNumSplines())
    {
        throw meshkernel::ConstraintError("Invalid spline index: {}, not in range 0 .. {}",
                                          splineIndex,
                                          GetNumSplines() - 1);
    }

    if (m_splineNodes[splineIndex].empty())
    {
        throw meshkernel::ConstraintError("Empty spline at index: {}", splineIndex);
    }

    const auto indices = FindIndices(m_splineNodes[splineIndex], 0, static_cast<UInt>(m_splineNodes[splineIndex].size()), constants::missing::doubleValue);

    for (auto [firstIndex, secondIndex] : indices)
    {
        std::vector<Point> splinePoints(m_splineNodes[splineIndex].begin() + firstIndex, m_splineNodes[splineIndex].begin() + secondIndex + 1);
        std::vector<Point> splineDerivative(m_splineDerivatives[splineIndex].begin() + firstIndex, m_splineDerivatives[splineIndex].begin() + secondIndex + 1);

        SplineAlgorithms::SnapSplineToBoundary(splinePoints, splineDerivative, landBoundary, m_projection, numberOfIterations);

        // Now that the spline has changed, need to update the second derivative.
        splineDerivative = SplineAlgorithms::SecondOrderDerivative(splinePoints, 0, splinePoints.size() - 1);

        // Copy updated spline and derivative back.
        for (auto j = firstIndex; j <= secondIndex; ++j)
        {
            m_splineNodes[splineIndex][j] = splinePoints[j - firstIndex];
            m_splineDerivatives[splineIndex][j] = splineDerivative[j - firstIndex];
        }
    }
}
