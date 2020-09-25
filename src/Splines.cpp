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
#include "CurvilinearGrid.hpp"

MeshKernel::Splines::Splines() : m_projection(Projections::cartesian)
{
};

MeshKernel::Splines::Splines(Projections projection) : m_projection(projection)
{
};

/// add a new spline, return the index
bool MeshKernel::Splines::AddSpline(const std::vector<Point>& splines, size_t start, size_t size)
{
    ResizeVectorIfNeededWithMinimumSize(m_numSplines + 1, m_splineNodes, m_allocationSize, std::vector<Point>(10, { doubleMissingValue, doubleMissingValue }));

    m_numAllocatedSplines = int(m_splineNodes.size());
    m_numAllocatedSplineNodes.resize(m_numAllocatedSplines, 10);

    m_numSplineNodes.resize(m_numAllocatedSplines, 0);
    m_numSplineNodes[m_numSplines] = size;

    m_splineDerivatives.resize(m_numAllocatedSplines);
    int index = 0;
    for (int n = start; n < start + size; ++n)
    {
        m_splineNodes[m_numSplines][index] = splines[n];
        index++;
    }
    m_splinesLength.resize(m_numAllocatedSplines);

    // compute basic properties
    SecondOrderDerivative(m_splineNodes[m_numSplines], m_numSplineNodes[m_numSplines], m_splineDerivatives[m_numSplines]);
    m_splinesLength[m_numSplines] = GetSplineLength(m_numSplines, 0, m_numSplineNodes[m_numSplines] - 1);
    m_numSplines++;

    return true;
}

bool MeshKernel::Splines::DeleteSpline(int splineIndex)
{
    m_splineNodes.erase(m_splineNodes.begin() + splineIndex);
    m_numSplineNodes.erase(m_numSplineNodes.begin() + splineIndex);
    m_splineDerivatives.erase(m_splineDerivatives.begin() + splineIndex);
    m_splinesLength.erase(m_splinesLength.begin() + splineIndex);
    m_numSplines--;
    return true;
}

/// add a new spline point in an existing spline
bool MeshKernel::Splines::AddPointInExistingSpline(int splineIndex, const Point& point)
{
    if (splineIndex > m_numSplines)
    {
        return false;
    }
    ResizeVectorIfNeededWithMinimumSize(m_numSplineNodes[splineIndex] + 1, m_splineNodes[splineIndex], m_allocationSize, { doubleMissingValue, doubleMissingValue });
    m_numAllocatedSplineNodes[splineIndex] = int(m_splineNodes[splineIndex].size());

    m_splineNodes[splineIndex][m_numSplineNodes[splineIndex]] = point;
    m_numSplineNodes[splineIndex]++;
    return true;
}

bool MeshKernel::Splines::GetSplinesIntersection(int first, 
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

    // First find a valid crossing, the closest to spline central point
    for (int n = 0; n < m_numSplineNodes[first] - 1; n++)
    {
        for (int nn = 0; nn < m_numSplineNodes[second] - 1; nn++)
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
                if (m_numSplineNodes[first] == 2)
                {
                    crossingDistance = std::min(minimumCrossingDistance, std::abs(firstRatio - 0.5));
                }
                else if (m_numSplineNodes[second] == 2)
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
                    firstCrossingIndex = n;             //TI0
                    secondCrossingIndex = nn;           //TJ0
                    firstCrossingRatio = firstRatio;    //SL
                    secondCrossingRatio = secondRatio;  //SM
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
    while (squaredDistanceBetweenCrossings > maxSquaredDistanceBetweenCrossings&& numIterations < 20)
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

        firstCrossing = std::max(0.0, std::min(firstCrossing, double(m_numSplineNodes[first])));
        secondCrossing = std::max(0.0, std::min(secondCrossing, double(m_numSplineNodes[second])));

        double firstLeft = std::max(0.0, std::min(double(m_numSplineNodes[first] - 1), firstCrossing - firstRatioIterations / 2.0));
        double firstRight = std::max(0.0, std::min(double(m_numSplineNodes[first] - 1), firstCrossing + firstRatioIterations / 2.0));

        double secondLeft = std::max(0.0, std::min(double(m_numSplineNodes[second] - 1), secondCrossing - secondRatioIterations / 2.0));
        double secondRight = std::max(0.0, std::min(double(m_numSplineNodes[second] - 1), secondCrossing + secondRatioIterations / 2.0));

        firstRatioIterations = firstRight - firstLeft;
        secondRatioIterations = secondRight - secondLeft;

        Point firstLeftSplinePoint{ doubleMissingValue, doubleMissingValue };
        bool successful = InterpolateSplinePoint( m_splineNodes[first], m_splineDerivatives[first],  firstLeft, firstLeftSplinePoint);
        if (!successful) 
        {
            return false;
        }

        Point firstRightSplinePoint{ doubleMissingValue, doubleMissingValue };
        successful = InterpolateSplinePoint(m_splineNodes[first], m_splineDerivatives[first], firstRight, firstRightSplinePoint);
        if (!successful)
        {
            return false;
        }

        Point secondLeftSplinePoint{ doubleMissingValue, doubleMissingValue };
        successful = InterpolateSplinePoint(m_splineNodes[second], m_splineDerivatives[second], secondLeft, secondLeftSplinePoint);
        if (!successful)
        {
            return false;
        }

        Point secondRightSplinePoint{ doubleMissingValue, doubleMissingValue };
        successful = InterpolateSplinePoint(m_splineNodes[second], m_splineDerivatives[second], secondRight, secondRightSplinePoint);
        if (!successful)
        {
            return false;
        }

        Point oldIntersection = closestIntersection;

        double crossProduct;
        double firstRatio = doubleMissingValue;
        double secondRatio = doubleMissingValue;
        bool areCrossing = AreLinesCrossing( firstLeftSplinePoint, 
                                             firstRightSplinePoint,
                                             secondLeftSplinePoint, 
                                             secondRightSplinePoint,
                                             true,
                                             closestIntersection,
                                             crossProduct,
                                             firstRatio,
                                             secondRatio,
                                             m_projection );

        // search close by
        if (-2.0 < firstRatio < 3.0 && -2.0 < secondRatio < 3.0)
        {
            previousFirstCrossing = firstCrossing;
            previousSecondCrossing = secondCrossing;

            firstCrossing = firstLeft + firstRatio * (firstRight - firstLeft);
            secondCrossing = secondLeft + secondRatio * (secondRight - secondLeft);

            firstCrossing = std::max(0.0, std::min(m_numSplineNodes[first] - 1.0, firstCrossing));
            secondCrossing = std::max(0.0, std::min(m_numSplineNodes[second] - 1.0, secondCrossing));

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

double MeshKernel::Splines::GetSplineLength(int index,
                                          double startIndex,
                                          double endIndex,
                                          int numSamples,
                                          bool accountForCurvature,
                                          double height,
                                          double assignedDelta)
{
    double delta = assignedDelta;
    int numPoints = int(endIndex / delta) + 1;
    if (delta < 0.0)
    {
        delta = 1.0 / numSamples;
        numPoints = std::max(std::floor(0.9999 + (endIndex - startIndex) / delta), 10.0);
        delta = (endIndex - startIndex) / numPoints;
    }

    // first point
    Point leftPoint;
    InterpolateSplinePoint(m_splineNodes[index], m_splineDerivatives[index], startIndex, leftPoint);

    double splineLength = 0.0;

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

        Point rightPoint;
        InterpolateSplinePoint(m_splineNodes[index], m_splineDerivatives[index], rightPointCoordinateOnSpline, rightPoint);
        double curvatureFactor = 0.0;
        if (accountForCurvature)
        {
            Point normalVector;
            Point tangentialVector;
            ComputeCurvatureOnSplinePoint(index, 0.5 * (rightPointCoordinateOnSpline + leftPointCoordinateOnSpline), curvatureFactor, normalVector, tangentialVector);
        }
        splineLength = splineLength + Distance(leftPoint, rightPoint, m_projection) * (1.0 + curvatureFactor * height);
        leftPoint = rightPoint;
    }

    return splineLength;
}

bool MeshKernel::Splines::ComputeCurvatureOnSplinePoint( int splineIndex,
                                                       double adimensionalPointCoordinate,
                                                       double& curvatureFactor,
                                                       Point& normalVector,
                                                       Point& tangentialVector)
{
    auto const leftCornerPoint = int(std::max(std::min(double(std::floor(adimensionalPointCoordinate)), double(m_numSplineNodes[splineIndex] - 1)), 0.0));
    auto const rightCornerPoint = int(std::max(double(leftCornerPoint + 1.0), 0.0));

    double leftSegment = rightCornerPoint - adimensionalPointCoordinate;
    double rightSegment = adimensionalPointCoordinate - leftCornerPoint;

    Point pointCoordinate;
    bool successful = InterpolateSplinePoint(m_splineNodes[splineIndex], m_splineDerivatives[splineIndex], adimensionalPointCoordinate, pointCoordinate);

    if (!successful) 
    {
        curvatureFactor = 0.0;
        return false;
    }


    Point p = m_splineNodes[splineIndex][rightCornerPoint] - m_splineNodes[splineIndex][leftCornerPoint] + 
              (m_splineDerivatives[splineIndex][leftCornerPoint] * (-3.0 * leftSegment * leftSegment + 1.0) + 
               m_splineDerivatives[splineIndex][rightCornerPoint] * (3.0 * rightSegment * rightSegment - 1.0)) / 6.0;

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

    double distance = Distance(pointCoordinate, incremenetedPointCoordinate, m_projection);
    double dx = GetDx(pointCoordinate, incremenetedPointCoordinate, m_projection);
    double dy = GetDy(pointCoordinate, incremenetedPointCoordinate, m_projection);

    tangentialVector.x = dx / distance;
    tangentialVector.y = dy / distance;

    return true;
}

bool MeshKernel::Splines::SecondOrderDerivative(const std::vector<Point>& spline, int numNodes, std::vector<Point>& coordinatesDerivatives)
{
    std::vector<Point> u(numNodes);
    u[0] = { 0.0, 0.0 };
    coordinatesDerivatives.resize(spline.size(), { 0.0, 0.0 });
    coordinatesDerivatives[0] = { 0.0, 0.0 };

    for (int i = 1; i < numNodes - 1; i++)
    {
        const Point p = coordinatesDerivatives[i - 1] * 0.5 + 2.0;
        coordinatesDerivatives[i].x = -0.5 / p.x;
        coordinatesDerivatives[i].y = -0.5 / p.y;

        const Point delta = spline[i + 1] - spline[i] - (spline[i] - spline[i - 1]);
        u[i] = (delta * 6.0 / 2.0 - u[i - 1] * 0.5) / p;
    }

    coordinatesDerivatives[numNodes - 1] = { 0.0, 0.0 };
    for (int i = numNodes - 2; i >= 0; i--)
    {
        coordinatesDerivatives[i] = coordinatesDerivatives[i] * coordinatesDerivatives[i + 1] + u[i];
    }

    return true;
}

bool MeshKernel::Splines::SecondOrderDerivative(const std::vector<double>& coordinates, int numNodes, std::vector<double>& coordinatesDerivatives)
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

    return true;
}

bool MeshKernel::Splines::InterpolatePointsOnSpline( int index,
                                                   double maximumGridHeight,
                                                   bool isSpacingCurvatureAdapted,
                                                   const std::vector<double>& distances,
                                                   std::vector<Point>& points,
                                                   std::vector<double>& adimensionalDistances)
{
    MeshKernel::FuncDimensionalToAdimensionalDistance func(this, index, isSpacingCurvatureAdapted, maximumGridHeight);

    for (int i = 0; i < distances.size(); i++)
    {
        func.SetDimensionalDistance(distances[i]);
        adimensionalDistances[i] = FindFunctionRootWithGoldenSectionSearch(func, 0, m_numSplineNodes[index] - 1);
        InterpolateSplinePoint(m_splineNodes[index], m_splineDerivatives[index], adimensionalDistances[i], points[i]);
    }

    return true;
}