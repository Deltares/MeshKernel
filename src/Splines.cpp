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
#include <stdexcept>

#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/Operations.cpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Splines.hpp>

meshkernel::Splines::Splines(){};

meshkernel::Splines::Splines(Projections projection) : m_projection(projection){};

/// add a new spline, return the index
void meshkernel::Splines::AddSpline(const std::vector<Point>& splines, size_t start, size_t size)
{
    m_splineNodes.emplace_back();
    for (auto n = start; n < start + size; ++n)
    {
        m_splineNodes.back().emplace_back(splines[n]);
    }

    // compute basic properties
    m_splineDerivatives.emplace_back();
    SecondOrderDerivative(m_splineNodes.back(), size, m_splineDerivatives.back());
    m_splinesLength.emplace_back(GetSplineLength(static_cast<int>(GetNumSplines() - 1), 0, size - 1));
}

void meshkernel::Splines::DeleteSpline(int splineIndex)
{
    m_splineNodes.erase(m_splineNodes.begin() + splineIndex);
    m_splineDerivatives.erase(m_splineDerivatives.begin() + splineIndex);
    m_splinesLength.erase(m_splinesLength.begin() + splineIndex);
}

void meshkernel::Splines::AddPointInExistingSpline(int splineIndex, const Point& point)
{
    if (splineIndex > GetNumSplines())
    {
        throw std::invalid_argument("Splines::AddPointInExistingSpline: Invalid spline index.");
    }
    m_splineNodes[splineIndex].emplace_back(point);
}

bool meshkernel::Splines::GetSplinesIntersection(int first,
                                                 int second,
                                                 double& crossProductIntersection,
                                                 Point& intersectionPoint,
                                                 double& firstSplineRatio,
                                                 double& secondSplineRatio)
{
    double minimumCrossingDistance = std::numeric_limits<double>::max();
    double crossingDistance;
    int numCrossing = 0;
    double firstCrossingRatio;
    double secondCrossingRatio;
    int firstCrossingIndex = 0;
    int secondCrossingIndex = 0;
    Point closestIntersection;
    const auto numNodesFirstSpline = static_cast<int>(m_splineNodes[first].size());
    const auto numNodesSecondSpline = static_cast<int>(m_splineNodes[second].size());

    // First find a valid crossing, the closest to spline central point
    for (int n = 0; n < numNodesFirstSpline - 1; n++)
    {
        for (int nn = 0; nn < numNodesSecondSpline - 1; nn++)
        {
            Point intersection;
            double crossProduct;
            double firstRatio;
            double secondRatio;
            bool areCrossing = AreLinesCrossing(m_splineNodes[first][n],
                                                m_splineNodes[first][n + 1],
                                                m_splineNodes[second][nn],
                                                m_splineNodes[second][nn + 1],
                                                false,
                                                intersection,
                                                crossProduct,
                                                firstRatio,
                                                secondRatio,
                                                m_projection);

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
                    firstCrossingIndex = n;            //TI0
                    secondCrossingIndex = nn;          //TJ0
                    firstCrossingRatio = firstRatio;   //SL
                    secondCrossingRatio = secondRatio; //SM
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

    double firstCrossing = firstCrossingRatio == -1 ? 0 : firstCrossingIndex + firstCrossingRatio;
    double secondCrossing = secondCrossingRatio == -1 ? 0 : secondCrossingIndex + secondCrossingRatio;

    // use bisection to find the intersection
    double squaredDistanceBetweenCrossings = std::numeric_limits<double>::max();
    double maxSquaredDistanceBetweenCrossings = 1e-12;
    double maxDistanceBetweenVertices = 0.0001;
    double firstRatioIterations = 1.0;
    double secondRatioIterations = 1.0;
    double previousFirstCrossing;
    double previousSecondCrossing;
    int numIterations = 0;
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

        double firstLeft = std::max(0.0, std::min(double(numNodesFirstSpline - 1), firstCrossing - firstRatioIterations / 2.0));
        double firstRight = std::max(0.0, std::min(double(numNodesFirstSpline - 1), firstCrossing + firstRatioIterations / 2.0));

        double secondLeft = std::max(0.0, std::min(double(numNodesSecondSpline - 1), secondCrossing - secondRatioIterations / 2.0));
        double secondRight = std::max(0.0, std::min(double(numNodesSecondSpline - 1), secondCrossing + secondRatioIterations / 2.0));

        firstRatioIterations = firstRight - firstLeft;
        secondRatioIterations = secondRight - secondLeft;

        Point firstLeftSplinePoint{doubleMissingValue, doubleMissingValue};
        bool successful = InterpolateSplinePoint(m_splineNodes[first], m_splineDerivatives[first], firstLeft, firstLeftSplinePoint);
        if (!successful)
        {
            return false;
        }

        Point firstRightSplinePoint{doubleMissingValue, doubleMissingValue};
        successful = InterpolateSplinePoint(m_splineNodes[first], m_splineDerivatives[first], firstRight, firstRightSplinePoint);
        if (!successful)
        {
            return false;
        }

        Point secondLeftSplinePoint{doubleMissingValue, doubleMissingValue};
        successful = InterpolateSplinePoint(m_splineNodes[second], m_splineDerivatives[second], secondLeft, secondLeftSplinePoint);
        if (!successful)
        {
            return false;
        }

        Point secondRightSplinePoint{doubleMissingValue, doubleMissingValue};
        successful = InterpolateSplinePoint(m_splineNodes[second], m_splineDerivatives[second], secondRight, secondRightSplinePoint);
        if (!successful)
        {
            return false;
        }

        Point oldIntersection = closestIntersection;

        double crossProduct;
        double firstRatio = doubleMissingValue;
        double secondRatio = doubleMissingValue;
        bool areCrossing = AreLinesCrossing(firstLeftSplinePoint,
                                            firstRightSplinePoint,
                                            secondLeftSplinePoint,
                                            secondRightSplinePoint,
                                            true,
                                            closestIntersection,
                                            crossProduct,
                                            firstRatio,
                                            secondRatio,
                                            m_projection);

        // search close by
        if (firstRatio > -2.0 && firstRatio < 3.0 && secondRatio > -2.0 && secondRatio < 3.0)
        {
            previousFirstCrossing = firstCrossing;
            previousSecondCrossing = secondCrossing;

            firstCrossing = firstLeft + firstRatio * (firstRight - firstLeft);
            secondCrossing = secondLeft + secondRatio * (secondRight - secondLeft);

            firstCrossing = std::max(0.0, std::min(numNodesFirstSpline - 1.0, firstCrossing));
            secondCrossing = std::max(0.0, std::min(numNodesSecondSpline - 1.0, secondCrossing));

            if (areCrossing)
            {
                numCrossing = 1;
                crossProductIntersection = crossProduct;
            }

            if (std::abs(firstCrossing - previousFirstCrossing) > maxDistanceBetweenVertices ||
                std::abs(secondCrossing - previousSecondCrossing) > maxDistanceBetweenVertices)
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

    //not crossing
    return false;
}

double meshkernel::Splines::GetSplineLength(int index,
                                            double startIndex,
                                            double endIndex,
                                            int numSamples,
                                            bool accountForCurvature,
                                            double height,
                                            double assignedDelta)
{
    double splineLength = 0.0;
    double delta = assignedDelta;
    int numPoints = int(endIndex / delta) + 1;
    if (delta < 0.0)
    {
        delta = 1.0 / numSamples;
        // TODO: Refactor or at least document the calculation of "numPoints"
        numPoints = int(std::max(std::floor(0.9999 + (endIndex - startIndex) / delta), 10.0));
        delta = (endIndex - startIndex) / numPoints;
    }

    // first point
    Point leftPoint{doubleMissingValue, doubleMissingValue};
    bool successful = InterpolateSplinePoint(m_splineNodes[index], m_splineDerivatives[index], startIndex, leftPoint);
    if (!successful)
    {
        throw AlgorithmError("Splines::GetSplineLength: Could not interpolate spline points.");
    }

    double rightPointCoordinateOnSpline = startIndex;
    double leftPointCoordinateOnSpline;
    for (int p = 0; p < numPoints; ++p)
    {
        leftPointCoordinateOnSpline = rightPointCoordinateOnSpline;
        rightPointCoordinateOnSpline += delta;
        if (rightPointCoordinateOnSpline > endIndex)
        {
            rightPointCoordinateOnSpline = endIndex;
        }

        Point rightPoint{doubleMissingValue, doubleMissingValue};
        successful = InterpolateSplinePoint(m_splineNodes[index], m_splineDerivatives[index], rightPointCoordinateOnSpline, rightPoint);
        if (!successful)
        {
            throw AlgorithmError("Splines::GetSplineLength: Could not interpolate spline points.");
        }

        double curvatureFactor = 0.0;
        if (accountForCurvature)
        {
            Point normalVector;
            Point tangentialVector;
            ComputeCurvatureOnSplinePoint(index, 0.5 * (rightPointCoordinateOnSpline + leftPointCoordinateOnSpline), curvatureFactor, normalVector, tangentialVector);
        }
        splineLength = splineLength + ComputeDistance(leftPoint, rightPoint, m_projection) * (1.0 + curvatureFactor * height);
        leftPoint = rightPoint;
    }

    return splineLength;
}

void meshkernel::Splines::ComputeCurvatureOnSplinePoint(int splineIndex,
                                                        double adimensionalPointCoordinate,
                                                        double& curvatureFactor,
                                                        Point& normalVector,
                                                        Point& tangentialVector)
{
    const auto numNodesFirstSpline = static_cast<int>(m_splineNodes[splineIndex].size());
    auto const leftCornerPoint = int(std::max(std::min(double(std::floor(adimensionalPointCoordinate)), static_cast<double>(numNodesFirstSpline - 1)), 0.0));
    auto const rightCornerPoint = int(std::max(double(leftCornerPoint + 1.0), 0.0));

    double leftSegment = rightCornerPoint - adimensionalPointCoordinate;
    double rightSegment = adimensionalPointCoordinate - leftCornerPoint;

    Point pointCoordinate;
    bool successful = InterpolateSplinePoint(m_splineNodes[splineIndex], m_splineDerivatives[splineIndex], adimensionalPointCoordinate, pointCoordinate);
    if (!successful)
    {
        curvatureFactor = 0.0;
        throw AlgorithmError("Splines::ComputeCurvatureOnSplinePoint: Could not interpolate spline points.");
    }

    Point p = m_splineNodes[splineIndex][rightCornerPoint] - m_splineNodes[splineIndex][leftCornerPoint] +
              (m_splineDerivatives[splineIndex][leftCornerPoint] * (-3.0 * leftSegment * leftSegment + 1.0) +
               m_splineDerivatives[splineIndex][rightCornerPoint] * (3.0 * rightSegment * rightSegment - 1.0)) /
                  6.0;

    Point pp = m_splineDerivatives[splineIndex][leftCornerPoint] * leftSegment +
               m_splineDerivatives[splineIndex][rightCornerPoint] * rightSegment;

    if (m_projection == Projections::spherical)
    {
        p.TransformSphericalToCartesian(pointCoordinate.y);
        pp.TransformSphericalToCartesian(pointCoordinate.y);
    }

    curvatureFactor = std::abs(pp.x * p.y - pp.y * p.x) / std::pow((p.x * p.x + p.y * p.y + 1e-8), 1.5);

    Point incremenetedPointCoordinate = pointCoordinate + p * 1e-4;
    NormalVectorOutside(pointCoordinate, incremenetedPointCoordinate, normalVector, m_projection);

    double distance = ComputeDistance(pointCoordinate, incremenetedPointCoordinate, m_projection);
    double dx = GetDx(pointCoordinate, incremenetedPointCoordinate, m_projection);
    double dy = GetDy(pointCoordinate, incremenetedPointCoordinate, m_projection);

    tangentialVector.x = dx / distance;
    tangentialVector.y = dy / distance;
}

void meshkernel::Splines::SecondOrderDerivative(const std::vector<Point>& spline, size_t numNodes, std::vector<Point>& coordinatesDerivatives)
{
    std::vector<Point> u(numNodes);
    u[0] = {0.0, 0.0};
    coordinatesDerivatives.resize(spline.size(), {0.0, 0.0});
    coordinatesDerivatives[0] = {0.0, 0.0};

    for (int i = 1; i < numNodes - 1; i++)
    {
        const Point p = coordinatesDerivatives[i - 1] * 0.5 + 2.0;
        coordinatesDerivatives[i].x = -0.5 / p.x;
        coordinatesDerivatives[i].y = -0.5 / p.y;

        const Point delta = spline[i + 1] - spline[i] - (spline[i] - spline[i - 1]);
        u[i] = (delta * 6.0 / 2.0 - u[i - 1] * 0.5) / p;
    }

    coordinatesDerivatives[numNodes - 1] = {0.0, 0.0};
    for (int i = numNodes - 2; i >= 0; i--)
    {
        coordinatesDerivatives[i] = coordinatesDerivatives[i] * coordinatesDerivatives[i + 1] + u[i];
    }
}

void meshkernel::Splines::SecondOrderDerivative(const std::vector<double>& coordinates, size_t numNodes, std::vector<double>& coordinatesDerivatives)
{
    std::vector<double> u(numNodes);
    u[0] = 0.0;
    coordinatesDerivatives[0] = 0.0;

    for (int i = 1; i < numNodes - 1; i++)
    {
        const double p = coordinatesDerivatives[i - 1] * 0.5 + 2.0;
        coordinatesDerivatives[i] = -0.5 / p;

        const double delta = coordinates[i + 1] - coordinates[i] - (coordinates[i] - coordinates[i - 1]);
        u[i] = (delta * 6.0 / 2.0 - u[i - 1] * 0.5) / p;
    }

    coordinatesDerivatives[numNodes - 1] = 0.0;
    for (int i = numNodes - 2; i >= 0; i--)
    {
        coordinatesDerivatives[i] = coordinatesDerivatives[i] * coordinatesDerivatives[i + 1] + u[i];
    }
}

void meshkernel::Splines::InterpolatePointsOnSpline(int index,
                                                    double maximumGridHeight,
                                                    bool isSpacingCurvatureAdapted,
                                                    const std::vector<double>& distances,
                                                    std::vector<Point>& points,
                                                    std::vector<double>& adimensionalDistances)
{
    FuncDimensionalToAdimensionalDistance func(this, index, isSpacingCurvatureAdapted, maximumGridHeight);
    const auto numNodes = static_cast<int>(m_splineNodes[index].size());
    for (size_t i = 0, size = distances.size(); i < size; ++i)
    {
        func.SetDimensionalDistance(distances[i]);
        adimensionalDistances[i] = FindFunctionRootWithGoldenSectionSearch(func, 0, numNodes - 1);
        auto successful = InterpolateSplinePoint(m_splineNodes[index], m_splineDerivatives[index], adimensionalDistances[i], points[i]);
        if (!successful)
        {
            throw AlgorithmError("Splines::InterpolatePointsOnSpline: Could not interpolate spline points.");
        }
    }
}
