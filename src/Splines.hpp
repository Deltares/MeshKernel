#pragma once

#include <vector>
#include <algorithm>
#include "Operations.cpp"
#include "Entities.hpp"
#include "Polygons.hpp"
#include "CurvilinearParametersNative.hpp"
#include "SplinesToCurvilinearParametersNative.hpp"

namespace GridGeom
{
    class Splines
    {
    public:

        Splines() : m_numAllocatedSplines(0), m_numSplines(0)
        {
            AllocateVector(m_numAllocatedSplines,
                m_splines,
                std::vector<Point>(10, { doubleMissingValue, doubleMissingValue }),
                5);

            m_numAllocatedSplineNodes.resize(m_numAllocatedSplines, 10);
            m_numSplineNodes.resize(m_numAllocatedSplines, 0);
        }

        /// add a new spline
        bool Set(const std::vector<Point>& splines)
        {
            AllocateVector(m_numAllocatedSplines, m_splines, std::vector<Point>(10, { doubleMissingValue, doubleMissingValue }), 5);

            m_numAllocatedSplineNodes.resize(m_numAllocatedSplines, 10);
            m_numSplineNodes.resize(m_numAllocatedSplines, 0);

            m_splines[m_numSplines] = splines;
            m_numSplineNodes[m_numSplines] = splines.size();
            m_numSplines++;
            return true;
        }

        /// add a new spline point in an existing spline
        bool Set(int splineIndex, const Point& point)
        {
            if (splineIndex >= m_numSplines)
            {
                return false;
            }
            AllocateVector(m_numAllocatedSplineNodes[splineIndex], m_splines[splineIndex], { doubleMissingValue, doubleMissingValue }, 10);

            m_splines[splineIndex][m_numSplineNodes[splineIndex]] = point;
            m_numSplineNodes[splineIndex]++;
            return true;
        }

        /// splint
        static bool Interpolate(const std::vector<Point>& coordinates, const std::vector<Point>& coordinatesDerivatives, double pointAdimensionalCoordinate, Point& pointCoordinate)
        {

            const double eps = 1e-5;
            const double splFac = 1.0;
            int intCoordinate = std::floor(pointAdimensionalCoordinate);
            if (pointAdimensionalCoordinate - intCoordinate < eps)
            {
                pointCoordinate = coordinates[intCoordinate];
                return true;
            }

            int low = intCoordinate;
            int high = low + 1;
            double a = high - pointAdimensionalCoordinate;
            double b = pointAdimensionalCoordinate - low;

            pointCoordinate = coordinates[low] * a + coordinates[high] * b +
                (coordinatesDerivatives[low] * (pow(a, 3) - a) + coordinatesDerivatives[high] * (pow(b, 3) - b)) / 6.0 * splFac;

            return true;
        }


        bool OrthogonalCurvilinearMeshFromSplines(
            const GridGeomApi::CurvilinearParametersNative& curvilinearParametersNative,
            const GridGeomApi::SplinesToCurvilinearParametersNative& splinesToCurvilinearParametersNative,
            const GridGeom::Polygons& polygon)
        {
            // no splines
            if (m_numSplines < 1)
            {
                return false;
            }

            // Deep copy of all members
            Splines oldSplines(*this);

            // Delete the splines that are not fully inside the polygons
            for (int s = 0; s < m_numSplines; s++)
            {
                for (int n = 0; n < m_numSplineNodes[s]; n++)
                {
                    bool isInPolygons = IsPointInPolygons(m_splines[s][n], polygon.m_nodes, polygon.m_numNodes);
                    if (!isInPolygons)
                    {
                        DeleteSpline(s);
                        break;
                    }
                }
            }

            // get the properties of the center splines (no angle check)


            return true;
        }

        /// get_crosssplines
        bool GetCrossedSplines(int position)
        {

            return true;
        }

        /// SECT3R
        /// compute the intersection of two splines
        bool GetSplinesIntersection(int first, int second,
            const Projections& projection,
            double& crossProductIntersection,
            Point& dimensionalIntersection,
            Point& adimensionalIntersection)
        {
            double minimumCrossingDistance = std::numeric_limits<double>::max();
            double crossingDistance;
            int numCrossing = 0;
            double firstCrossingRatio;
            double secondCrossingRatio;
            int firstCrossingIndex = 0;
            int secondCrossingIndex = 0;
            Point closestIntersection;
            for (int n = 0; n < m_numSplineNodes[first] - 1; n++)
            {
                for (int nn = 0; nn < m_numSplineNodes[second] - 1; nn++)
                {
                    Point intersection;
                    double crossProduct;
                    double firstRatio;
                    double secondRatio;
                    bool areCrossing = GridGeom::AreLinesCrossing(m_splines[first][n],
                        m_splines[first][n + 1],
                        m_splines[second][nn],
                        m_splines[second][nn + 1],
                        false,
                        intersection,
                        crossProduct,
                        firstRatio,
                        secondRatio,
                        projection);


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
                            firstCrossingIndex = n;
                            secondCrossingIndex = nn;
                            firstCrossingRatio = firstRatio;
                            secondCrossingRatio = secondRatio;
                        }
                    }
                    closestIntersection = intersection;
                }
            }

            if (numCrossing == 0)
            {
                return false;
            }

            std::vector<Point> firstSpline(m_numSplineNodes[first]);
            std::vector<Point> secondSpline(m_numSplineNodes[second]);

            bool successful = Derivative(m_splines[first], m_numSplineNodes[first], firstSpline);
            successful = Derivative(m_splines[second], m_numSplineNodes[second], secondSpline);

            double firstCrossing = firstCrossingRatio == -1 ? 0 : firstCrossingIndex + firstCrossingRatio;
            double secondCrossing = secondCrossingRatio == -1 ? 0 : secondCrossingIndex + secondCrossingRatio;

            // use bisection to find the intersection 
            double distanceBetweenCrossings = std::numeric_limits<double>::max();
            double maxDistanceBetweenCrossings = 0.000001;
            double maxDistanceBetweenVertices = 0.0001;
            double firstRatioIterations = 1.0;
            double secondRatioIterations = 1.0;
            double previousFirstCrossing = firstCrossing;
            double previousSecondCrossing = secondCrossing;
            int numIterations = 0;
            while (distanceBetweenCrossings > maxDistanceBetweenCrossings && numIterations < 20)
            {
                numIterations++;
                if (firstCrossingRatio > 0 && firstCrossingRatio < 1.0)
                {
                    firstRatioIterations = 0.5 * firstRatioIterations;
                }
                if (secondCrossingRatio > 0 && secondCrossingRatio < 1.0)
                {
                    secondRatioIterations = 0.5 * secondRatioIterations;
                }

                double firstLeft = std::min(m_numSplineNodes[first] - 1.0, firstCrossing - firstRatioIterations / 2.0);
                double firstRight = std::min(m_numSplineNodes[first] - 1.0, firstCrossing + firstRatioIterations / 2.0);

                double secondLeft = std::min(m_numSplineNodes[second] - 1.0, secondCrossing - secondRatioIterations / 2.0);
                double secondRight = std::min(m_numSplineNodes[second] - 1.0, secondCrossing + secondRatioIterations / 2.0);

                Point firstLeftSplinePoint;
                Interpolate(m_splines[first], firstSpline, firstLeft, firstLeftSplinePoint);
                Point firstRightSplinePoint;
                Interpolate(m_splines[first], firstSpline, firstRight, firstRightSplinePoint);

                Point secondLeftSplinePoint;
                Interpolate(m_splines[second], secondSpline, secondLeft, secondLeftSplinePoint);
                Point secondRightSplinePoint;
                Interpolate(m_splines[second], secondSpline, secondRight, secondRightSplinePoint);

                Point oldIntersection = closestIntersection;

                double crossProduct;
                double firstRatio;
                double secondRatio;
                bool areCrossing = AreLinesCrossing(firstLeftSplinePoint, firstRightSplinePoint,
                    secondLeftSplinePoint, secondRightSplinePoint,
                    true,
                    closestIntersection,
                    crossProduct,
                    firstRatio,
                    secondRatio,
                    projection);

                // outside the segments
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
                        distanceBetweenCrossings = Distance(oldIntersection, closestIntersection, projection);
                    }
                    else
                    {
                        break;
                    }
                }
            }

            if (numCrossing == 1)
            {
                dimensionalIntersection = closestIntersection;
                adimensionalIntersection = { firstCrossing, secondCrossing };
            }
            return true;
        }


        bool DeleteSpline(int position)
        {
            m_splines.erase(m_splines.begin() + position);
            m_numSplineNodes.erase(m_numSplineNodes.begin() + position);
            m_numSplines--;
            return true;
        }

        /// SPLINE
        static bool Derivative(const std::vector<Point>& coordinates, int numNodes, std::vector<Point>& coordinatesDerivatives)
        {
            std::vector<Point> u(numNodes);
            u[0] = { 0.0, 0.0 };
            coordinatesDerivatives[0] = { 0.0, 0.0 };

            for (int i = 1; i < numNodes - 1; i++)
            {
                const Point p = coordinatesDerivatives[i - 1] * 0.5 + 2.0;
                coordinatesDerivatives[i].x = -0.5 / p.x;
                coordinatesDerivatives[i].y = -0.5 / p.y;

                const Point delta = coordinates[i + 1] - coordinates[i] - (coordinates[i] - coordinates[i - 1]);
                u[i] = (delta *6.0 / 2.0 - u[i - 1] * 0.5) / p;
            }

            coordinatesDerivatives[numNodes - 1] = { 0.0, 0.0 };
            for (int i = numNodes - 2; i >= 0; i--)
            {
                coordinatesDerivatives[i] = coordinatesDerivatives[i] * coordinatesDerivatives[i + 1] + u[i];
            }

            return true;
        }

        int m_numSplines;
        int m_numAllocatedSplines;
        std::vector<int> m_numSplineNodes;
        std::vector<int> m_numAllocatedSplineNodes;
        std::vector<std::vector<Point>> m_splines;
    };

}