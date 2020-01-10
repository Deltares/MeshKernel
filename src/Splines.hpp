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

        Splines() : m_numAllocatedSplines(5), m_numSplines(0)
        {
            AllocateSplines(m_numAllocatedSplines);
        }

        /// add a new spline
        bool Set(const std::vector<Point>& splines)
        {
            AllocateSplines(m_numSplines + 1);
            m_splines[m_numSplines] = splines;
            m_numSplineNodes[m_numSplines] = splines.size();
            m_numSplines++;
            return true;
        }

        bool AllocateSplines(int numSplines)
        {
            AllocateVector(m_numAllocatedSplines, m_splines, std::vector<Point>(10, { doubleMissingValue, doubleMissingValue }), numSplines);
            m_numAllocatedSplineNodes.resize(m_numAllocatedSplines, 10);
            m_numSplineNodes.resize(m_numAllocatedSplines, 0);
            return true;
        }

        bool AllocateSplineProperties(int numSplines)
        {
            m_id.resize(numSplines);
            m_numCrossedSplines.resize(numSplines);
            m_splinesLength.resize(numSplines);
            m_maximumGridHeight.resize(numSplines);
            m_intersectedSplinesIndexes.resize(numSplines);
            m_isLeftOriented.resize(numSplines);
            m_splineCrossings.resize(numSplines);
            m_cosCrossingAngle.resize(numSplines);
            m_hL.resize(numSplines);
            m_hR.resize(numSplines);
            m_NsubL.resize(numSplines);
            m_NsubR.resize(numSplines);
            m_mfac.resize(numSplines);
            m_nfacL.resize(numSplines);
            m_nfacR.resize(numSplines);
            m_iL.resize(numSplines);
            m_iR.resize(numSplines);
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

            // if no crossing found, return
            if (numCrossing == 0)
            {
                return false;
            }

            std::vector<Point> firstSpline(m_numSplineNodes[first]);
            std::vector<Point> secondSpline(m_numSplineNodes[second]);
            bool successful = SecondOrderDerivative(m_splines[first], m_numSplineNodes[first], firstSpline);
            if (!successful)
                return false;
            successful = SecondOrderDerivative(m_splines[second], m_numSplineNodes[second], secondSpline);
            if (!successful)
                return false;

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
                intersectionPoint = closestIntersection;
                firstSplineRatio = firstCrossing;
                secondSplineRatio = secondCrossing;
                return true;
            }

            //not crossing
            return false;
        }

        /// get_crosssplines
        /// compute the intersection of two splines, one must have only two nodes
        bool GetSplineIntersections(int index,
            const Projections& projection,
            int& numCrossedSplines,
            std::vector<int>& intersectedSplines,
            std::vector<bool>& orientation,
            std::vector<double>& adimensionalCrossPoints,
            std::vector<double>& crossingAngles)
        {
            numCrossedSplines = 0;
            intersectedSplines.resize(m_numSplines, -1);
            orientation.resize(m_numSplines,false);
            adimensionalCrossPoints.resize(m_numSplines,std::numeric_limits<double>::max());
            crossingAngles.resize(m_numSplines,-1.0);

            for (int s = 0; s < m_numSplines; ++s)
            {
                if (m_numSplineNodes[s] == 2 && m_numSplineNodes[index] == 2 ||
                    m_numSplineNodes[s] > 2 && m_numSplineNodes[index] > 2)
                {
                    continue;
                }

                double crossProductIntersection;
                Point intersectionPoint;
                double firstSplineRatio;
                double secondSplineRatio;                
                bool crossing = GetSplinesIntersection(index, s, projection, crossProductIntersection, intersectionPoint, firstSplineRatio, secondSplineRatio);

                if(crossing)
                {
                    numCrossedSplines++;
                    intersectedSplines[s] = s;
                    if(crossProductIntersection>0.0)
                    {
                        orientation[s] = true;
                    }
                    adimensionalCrossPoints[s]= firstSplineRatio;
                    crossingAngles[s]= crossProductIntersection;
                }
            }

            const auto sortedIndexses = SortedIndexes(adimensionalCrossPoints);
            ReorderVector(intersectedSplines, sortedIndexses);
            ReorderVector(orientation, sortedIndexses);
            ReorderVector(crossingAngles, sortedIndexses);

            return true;
        }

        double GetSplineLength(int index)
        {
            std::vector<Point> coordinatesDerivatives(m_numSplineNodes[index]);
            SecondOrderDerivative(m_splines[index], m_numSplineNodes[index], coordinatesDerivatives);

            // first point
            double pointAdimensionalCoordinate = 0.0;
            Point leftPoint;
            Interpolate(m_splines[index], coordinatesDerivatives, pointAdimensionalCoordinate, leftPoint);

            const int numPoints = 100;
            const double delta = 1.0 / numPoints;
            double splineLength = 0.0;

            for (int n = 0; n < m_numSplineNodes[index]; ++n)
            {
                pointAdimensionalCoordinate = double(n);
                for (int p = 0; p < numPoints; ++p)
                {
                    pointAdimensionalCoordinate = pointAdimensionalCoordinate + delta;
                    Point rightPoint;
                    Interpolate(m_splines[index], coordinatesDerivatives, pointAdimensionalCoordinate, rightPoint);
                    splineLength = splineLength + Distance(leftPoint, rightPoint, m_projection);
                    leftPoint = rightPoint;
                }
            }

            return splineLength;
        }

        ///get_splineprops
        bool GetSplineProperties()
        {

            AllocateSplineProperties(m_numSplines);

            for (int s = 0; s < m_numSplines; ++s)
            {
                m_splinesLength[s] = GetSplineLength(s);
                bool successful = GetSplineIntersections(s, m_projection, m_numCrossedSplines[s], m_intersectedSplinesIndexes[s], m_isLeftOriented[s], m_splineCrossings[s], m_cosCrossingAngle[s]);


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
        /// second order derivative of spline coordinates
        static bool SecondOrderDerivative(const std::vector<Point>& coordinates, int numNodes, std::vector<Point>& coordinatesDerivatives)
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


        int m_numSplines;
        int m_numAllocatedSplines;
        std::vector<int> m_numSplineNodes;
        std::vector<int> m_numAllocatedSplineNodes;
        std::vector<std::vector<Point>> m_splines;
        Projections m_projection;

        // Spline properties (first index is the spline number)
        std::vector<int> m_id;                                     // m_id maximum number of subintervals of grid layers, each having their own exponential grow factor
        std::vector<int> m_numCrossedSplines;                      // m_numCrossedSplines center spline type that contains information derived from cross splines
        std::vector<double> m_splinesLength;                       // m_splinesLength spline path length
        std::vector<double> m_maximumGridHeight;                   // m_hmax maximum grid height
        std::vector<std::vector<int>> m_intersectedSplinesIndexes; // m_intersectedSplinesIndexes cross spline numbers
        std::vector<std::vector<bool>> m_isLeftOriented;           // m_isLeftOriented cross spline is left to right(.true.) or not (.false.) w.r.t.center spline
        std::vector<std::vector<double>>  m_splineCrossings;       // m_t center spline coordinates of cross splines
        std::vector<std::vector<double>> m_cosCrossingAngle;       // m_cosCrossingAngle cosine of crossing angle
        std::vector<std::vector<std::vector<double>>> m_hL;        // left - hand side grid heights at cross spline locations for each grid layer subinterval, hL(1, :) being the height of the first subinterval, etc.
        std::vector<std::vector<std::vector<double>>> m_hR;        // right - hand side grid heights at cross spline locations for each grid layer subinterval, hR(1, :) being the height of the first subinterval, etc.
        std::vector<std::vector<double>>   m_NsubL;                // number of subintervals of grid layers at cross spline locations at the left - hand side of the spline, each having their own exponential grow factor
        std::vector<std::vector<double>>  m_NsubR;                 // number of subintervals of grid layers at cross spline locations at the right - hand side of the spline, each having their own exponential grow factor
        std::vector<int> m_mfac;                                   // number of grid intervals on the spline
        std::vector<std::vector<int>> m_nfacL;                     // number of grid layers in each subinterval at the left - hand side of the spline * not used yet*
        std::vector<std::vector<int>> m_nfacR;                     // number of grid layers in each subinterval at the right - hand side of the spline * not used yet*
        std::vector<int> m_iL;                                     // index in the whole gridline array of the first grid point on the left - hand side of the spline
        std::vector<int> m_iR;                                     // index in the whole gridline array of the first grid point on the right - hand side of the spline
    };

}