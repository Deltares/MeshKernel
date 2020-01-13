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

        Splines& operator=(const Splines& other)
        {
            if (this == &other)
                return *this;
            m_numSplines = other.m_numSplines;
            m_numAllocatedSplines = other.m_numAllocatedSplines;
            m_numSplineNodes = other.m_numSplineNodes;
            m_numAllocatedSplineNodes = other.m_numAllocatedSplineNodes;
            m_splines = other.m_splines;
            m_projection = other.m_projection;
            m_numLayers = other.m_numLayers;
            m_numCrossingSplines = other.m_numCrossingSplines;
            m_centerSplinesIndexes = other.m_centerSplinesIndexes;
            m_splinesLength = other.m_splinesLength;
            m_maximumGridHeight = other.m_maximumGridHeight;
            m_isLeftOriented = other.m_isLeftOriented;
            m_centerSplineCrossings = other.m_centerSplineCrossings;
            m_cosCrossingAngle = other.m_cosCrossingAngle;
            m_crossSplineLeftHeights = other.m_crossSplineLeftHeights;
            m_crossSplineRightHeights = other.m_crossSplineRightHeights;
            m_numHeightsCrossSplineLeft = other.m_numHeightsCrossSplineLeft;
            m_numHeightsCrossSplineRight = other.m_numHeightsCrossSplineRight;
            m_mfac = other.m_mfac;
            m_nfacL = other.m_nfacL;
            m_nfacR = other.m_nfacR;
            m_iL = other.m_iL;
            m_iR = other.m_iR;
            return *this;
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
            m_numLayers.resize(numSplines);
            m_numCrossingSplines.resize(numSplines);
            m_splinesLength.resize(numSplines);
            m_maximumGridHeight.resize(numSplines);
            m_centerSplinesIndexes.resize(numSplines);
            m_isLeftOriented.resize(numSplines);
            m_centerSplineCrossings.resize(numSplines);
            m_cosCrossingAngle.resize(numSplines);
            m_crossSplineLeftHeights.resize(numSplines);
            m_crossSplineRightHeights.resize(numSplines);
            m_numHeightsCrossSplineLeft.resize(numSplines);
            m_numHeightsCrossSplineRight.resize(numSplines);
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
        bool GetSplinesIntersection(const int first, const int second,
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
        bool GetSplineIntersections(const int index,
            const Projections& projection,
            int& numCrossedSplines,
            std::vector<int>& intersectedSplines,
            std::vector<bool>& orientation,
            std::vector<double>& adimensionalCrossPoints,
            std::vector<double>& crossingAngles)
        {
            numCrossedSplines = 0;
            intersectedSplines.resize(m_numSplines, -1);
            orientation.resize(m_numSplines, false);
            adimensionalCrossPoints.resize(m_numSplines, std::numeric_limits<double>::max());
            crossingAngles.resize(m_numSplines, -1.0);

            for (int s = 0; s < m_numSplines; ++s)
            {
                // a crossing is a spline with 2 nodes and another with more than 2 nodes
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

                if (crossing)
                {
                    numCrossedSplines++;
                    intersectedSplines[s] = s;
                    if (crossProductIntersection > 0.0)
                    {
                        orientation[s] = true;
                    }
                    adimensionalCrossPoints[s] = firstSplineRatio;
                    crossingAngles[s] = crossProductIntersection;
                }
            }

            const auto sortedIndexses = SortedIndexes(adimensionalCrossPoints);
            ReorderVector(intersectedSplines, sortedIndexses);
            ReorderVector(orientation, sortedIndexses);
            ReorderVector(crossingAngles, sortedIndexses);

            return true;
        }

        ///get_splineprops
        bool ComputeSplineProperties()
        {
            // cross spline is a spline with only two points, others are non-cross splines
            AllocateSplineProperties(m_numSplines);
            for (int s = 0; s < m_numSplines; ++s)
            {
                m_splinesLength[s] = GetSplineLength(s, 0, m_numSplineNodes[s]);
                bool successful = GetSplineIntersections(s, m_projection, m_numCrossingSplines[s], 
                    m_centerSplinesIndexes[s], m_isLeftOriented[s], m_centerSplineCrossings[s], m_cosCrossingAngle[s]);
            }

            for (int s = 0; s < m_numSplines; ++s)
            {
                m_numLayers[s] = 1;
                // select all non-cross splines for growing the grid
                if (m_numSplineNodes[s] > 2)
                {
                    m_numLayers[s] = 0;
                }
            }

            for (int s = 0; s < m_numSplines; ++s)
            {
                if (m_numCrossingSplines[s] < 1)
                {
                    continue;
                }

                int middleIndex = std::min(m_numCrossingSplines[s] / 2 + 1, m_numCrossingSplines[s]);
                int middleSplineIndex = m_centerSplinesIndexes[s][middleIndex];

                // check if the middle spline has already been assigned as a bounding spline
                if (m_numLayers[middleSplineIndex] != 0 && 2 * middleSplineIndex == m_numCrossingSplines[s])
                {
                    middleIndex = std::min(middleIndex + 1, m_numCrossingSplines[s]);
                    middleSplineIndex = m_centerSplinesIndexes[s][middleIndex];
                }
                else if (m_numLayers[middleSplineIndex] == 0)
                {
                    // associate bounding splines with the middle spline
                    for (int i = 0; i < middleIndex; ++i)
                    {
                        int index = m_centerSplinesIndexes[s][i];
                        m_numLayers[index] = -middleSplineIndex;

                    }
                    for (int i = middleIndex; i < m_numCrossingSplines[s]; ++i)
                    {
                        int index = m_centerSplinesIndexes[s][i];
                        m_numLayers[index] = -middleSplineIndex;
                    }
                }
            }

            bool successfull = ComputeHeights();


            return true;
        }

        // get_heights
        // get the grid heights from the cross spline information
        // get_heights
        bool ComputeHeights()
        {
            for (int s = 0; s < m_numSplines; ++s)
            {
                if (m_numSplineNodes[s] <= 2)
                {
                    // not a center spline
                    continue;
                }
                for (int c = 0; c < m_numCrossingSplines[s]; ++c)
                {
                    int crossingSplineIndex = m_centerSplinesIndexes[s][c];

                    //for this cross spline, find the left and right center splines
                    int leftCenterSpline = 0;
                    int rightCenterSpline = 0;
                    for (int i = 0; i < m_numCrossingSplines[crossingSplineIndex]; ++i)
                    {
                        int centerSplineOfCrossingSpline = m_centerSplinesIndexes[crossingSplineIndex][i];

                        // we have found a couple!
                        if (centerSplineOfCrossingSpline == s)
                        {
                            for (int kk = c - 1; kk >= 0; --kk)
                            {
                                int index = m_centerSplinesIndexes[crossingSplineIndex][kk];
                                if (m_numLayers[index] == -centerSplineOfCrossingSpline)
                                {
                                    leftCenterSpline = index;
                                    break;
                                }
                            }

                            for (int kk = c + 1; kk < m_numCrossingSplines[crossingSplineIndex]; ++kk)
                            {
                                int index = m_centerSplinesIndexes[crossingSplineIndex][kk];
                                if (m_numLayers[index] == -centerSplineOfCrossingSpline)
                                {
                                    rightCenterSpline = index;
                                    break;
                                }
                            }

                            break;
                        }
                    }

                    bool success = ComputeSubHeights(crossingSplineIndex, s);

                    if(!success)
                    {
                        return false;
                    }
                }
            }


            return true;
        }

        ///comp_subheights, compute the lengths of a cross spline (left and right of the spline)
        bool ComputeSubHeights(const int crossSplineIndex,const int centerSplineIndex)
        {
            // find center spline index
            int centerSplineIndexInVector = 0;
            for (int s = 0; s < m_numCrossingSplines[centerSplineIndex]; ++s)
            {
                if (m_centerSplinesIndexes[crossSplineIndex][s] == centerSplineIndex)
                {
                    centerSplineIndexInVector = s;
                    break;
                }
            }

            // right part
            int numSubIntervalsRight = 0;
            int rightCenterSplineIndex = centerSplineIndexInVector;
            int leftCenterSplineIndex;
            m_crossSplineRightHeights[crossSplineIndex][centerSplineIndex].resize(m_maxNumHeights, 0);
            for (int s = centerSplineIndexInVector; s < m_numCrossingSplines[centerSplineIndex] - 1; ++s)
            {
                if (numSubIntervalsRight >= m_maxNumHeights)
                {
                    break;
                }
                if (m_numLayers[m_centerSplinesIndexes[crossSplineIndex][s + 1]] != -centerSplineIndex)
                {
                    continue;
                }
                leftCenterSplineIndex = rightCenterSplineIndex;
                rightCenterSplineIndex = s + 1;
                m_crossSplineRightHeights[crossSplineIndex][centerSplineIndex][numSubIntervalsRight] = GetSplineLength(crossSplineIndex,
                    m_centerSplineCrossings[centerSplineIndex][leftCenterSplineIndex], m_centerSplineCrossings[centerSplineIndex][rightCenterSplineIndex]);
                numSubIntervalsRight++;
            }

            m_crossSplineRightHeights[crossSplineIndex][centerSplineIndex][numSubIntervalsRight] = GetSplineLength(crossSplineIndex,
                m_centerSplineCrossings[centerSplineIndex][rightCenterSplineIndex], m_centerSplineCrossings[centerSplineIndex][m_numSplineNodes[crossSplineIndex] - 1]);

            numSubIntervalsRight++;
            std::fill(m_crossSplineRightHeights[crossSplineIndex][centerSplineIndex].begin() + numSubIntervalsRight, m_crossSplineRightHeights[crossSplineIndex][centerSplineIndex].end(), 0.0);
            m_numHeightsCrossSplineRight[crossSplineIndex][centerSplineIndex] = numSubIntervalsRight;

            // left part
            int numSubIntervalsLeft = 0;
            leftCenterSplineIndex = centerSplineIndexInVector;
            m_crossSplineLeftHeights[crossSplineIndex][centerSplineIndex].resize(m_maxNumHeights, 0);
            for (int s = centerSplineIndexInVector; s >= 1; --s)
            {
                if (numSubIntervalsLeft >= m_maxNumHeights)
                {
                    break;
                }
                if (m_numLayers[m_centerSplinesIndexes[crossSplineIndex][s - 1]] != -centerSplineIndex)
                {
                    continue;
                }
                rightCenterSplineIndex = leftCenterSplineIndex;
                leftCenterSplineIndex = s - 1;
                m_crossSplineLeftHeights[crossSplineIndex][centerSplineIndex][numSubIntervalsLeft] = GetSplineLength(crossSplineIndex,
                    m_centerSplineCrossings[centerSplineIndex][leftCenterSplineIndex], m_centerSplineCrossings[centerSplineIndex][rightCenterSplineIndex]);
                numSubIntervalsLeft++;
            }

            m_crossSplineLeftHeights[crossSplineIndex][centerSplineIndex][numSubIntervalsLeft] = GetSplineLength(crossSplineIndex,
                0.0, m_centerSplineCrossings[centerSplineIndex][leftCenterSplineIndex]);

            numSubIntervalsLeft++;
            std::fill(m_crossSplineLeftHeights[crossSplineIndex][centerSplineIndex].begin() + numSubIntervalsLeft, m_crossSplineLeftHeights[crossSplineIndex][centerSplineIndex].end(), 0.0);
            m_numHeightsCrossSplineLeft[crossSplineIndex][centerSplineIndex] = numSubIntervalsLeft;


            // if not left oriented, swap
            if (!m_isLeftOriented[crossSplineIndex][centerSplineIndex])
            {
                m_numHeightsCrossSplineLeft[crossSplineIndex][centerSplineIndex] = numSubIntervalsRight;
                m_numHeightsCrossSplineRight[crossSplineIndex][centerSplineIndex] = numSubIntervalsLeft;

                std::vector<double> leftSubIntervalsTemp(m_crossSplineLeftHeights[crossSplineIndex][centerSplineIndex]);
                m_crossSplineLeftHeights[crossSplineIndex][centerSplineIndex] = m_crossSplineRightHeights[crossSplineIndex][centerSplineIndex];
                m_crossSplineRightHeights[crossSplineIndex][centerSplineIndex] = leftSubIntervalsTemp;
            }

            return true;
        }

        ///splinelength
        double GetSplineLength(int index, double beginFactor, double endFactor)
        {
            std::vector<Point> coordinatesDerivatives(m_numSplineNodes[index]);
            SecondOrderDerivative(m_splines[index], m_numSplineNodes[index], coordinatesDerivatives);

            double delta = 1.0 / 100.0;
            double numPoints = std::max(std::floor(0.9999 + (endFactor - beginFactor) / delta), 10.0);
            delta = (endFactor - beginFactor) / numPoints;

            // first point
            Point leftPoint;
            Interpolate(m_splines[index], coordinatesDerivatives, beginFactor, leftPoint);

            double splineLength = 0.0;
            for (int n = 0; n < m_numSplineNodes[index]; ++n)
            {
                double pointAdimensionalCoordinate = double(n);
                if (pointAdimensionalCoordinate<beginFactor || pointAdimensionalCoordinate>endFactor)
                {
                    continue;
                }
                for (int p = 0; p < numPoints; ++p)
                {
                    pointAdimensionalCoordinate = beginFactor + delta;
                    Point rightPoint;
                    Interpolate(m_splines[index], coordinatesDerivatives, pointAdimensionalCoordinate, rightPoint);
                    splineLength = splineLength + Distance(leftPoint, rightPoint, m_projection);
                    leftPoint = rightPoint;
                }
            }

            return splineLength;
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
        const int m_maxNumHeights = 10;                                                 // Nsubmax maximum number of subintervals of grid layers, each having their own exponential grow factor
        std::vector<int> m_numLayers;                                                   // id maximum number of subintervals of grid layers, each having their own exponential grow factor
        std::vector<int> m_numCrossingSplines;                                          // ncs num of cross splines
        std::vector<std::vector<int>> m_centerSplinesIndexes;                           // ics for each cross spline, the indexses of the center splines
        std::vector<double> m_splinesLength;                                            // splinesLength spline path length
        std::vector<double> m_maximumGridHeight;                                        // hmax maximum grid height
        std::vector<std::vector<bool>> m_isLeftOriented;                                // isLeftOriented cross spline is left to right(.true.) or not (.false.) w.r.t.center spline
        std::vector<std::vector<double>>  m_centerSplineCrossings;                      // t center spline coordinates of cross splines
        std::vector<std::vector<double>> m_cosCrossingAngle;                            // cosPhi cosine of crossing angle
        std::vector<std::vector<std::vector<double>>> m_crossSplineLeftHeights;         // left - hand side grid heights at cross spline locations for each grid layer subinterval, hL(1, :) being the height of the first subinterval, etc.
        std::vector<std::vector<std::vector<double>>> m_crossSplineRightHeights;        // right - hand side grid heights at cross spline locations for each grid layer subinterval, hR(1, :) being the height of the first subinterval, etc.
        std::vector<std::vector<double>>   m_numHeightsCrossSplineLeft;                 // number of subintervals of grid layers at cross spline locations at the left - hand side of the spline, each having their own exponential grow factor
        std::vector<std::vector<double>>  m_numHeightsCrossSplineRight;                 // number of subintervals of grid layers at cross spline locations at the right - hand side of the spline, each having their own exponential grow factor
        std::vector<int> m_mfac;                                                        // number of grid intervals on the spline
        std::vector<std::vector<int>> m_nfacL;                                          // number of grid layers in each subinterval at the left - hand side of the spline * not used yet*
        std::vector<std::vector<int>> m_nfacR;                                          // number of grid layers in each subinterval at the right - hand side of the spline * not used yet*
        std::vector<int> m_iL;                                                          // index in the whole gridline array of the first grid point on the left - hand side of the spline
        std::vector<int> m_iR;                                                          // index in the whole gridline array of the first grid point on the right - hand side of the spline
    };

}