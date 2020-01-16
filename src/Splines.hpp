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


        Splines() : m_projection(Projections::cartesian), m_numAllocatedSplines(5), m_numSplines(0)
        {
        }

        Splines(Projections projection) : m_projection(projection), m_numAllocatedSplines(0), m_numSplines(0)
        {
            AllocateSplines(5);
            m_aspectRatioFirstLayer = 0.10;
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
            m_numIntersectingSplines = other.m_numIntersectingSplines;
            m_intersectingSplinesIndexses = other.m_intersectingSplinesIndexses;
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
            m_leftGridLineIndex = other.m_leftGridLineIndex;
            m_rightGridLineIndex = other.m_rightGridLineIndex;
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
            m_numIntersectingSplines.resize(numSplines);
            m_splinesLength.resize(numSplines);
            m_maximumGridHeight.resize(numSplines);
            m_intersectingSplinesIndexses.resize(numSplines);
            m_isLeftOriented.resize(numSplines);
            m_centerSplineCrossings.resize(numSplines);
            m_cosCrossingAngle.resize(numSplines);
            m_crossSplineLeftHeights.resize(numSplines, std::vector<std::vector<double>>(numSplines));
            m_crossSplineRightHeights.resize(numSplines, std::vector<std::vector<double>>(numSplines));
            m_numHeightsCrossSplineLeft.resize(numSplines, std::vector<int>(numSplines));
            m_numHeightsCrossSplineRight.resize(numSplines, std::vector<int>(numSplines));
            m_mfac.resize(numSplines);
            m_nfacL.resize(numSplines);
            m_nfacR.resize(numSplines);
            m_leftGridLineIndex.resize(numSplines);
            m_rightGridLineIndex.resize(numSplines);

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

        bool OrthogonalCurvilinearMeshFromSplines(
            const GridGeomApi::CurvilinearParametersNative& curvilinearParametersNative,
            const GridGeomApi::SplinesToCurvilinearParametersNative& splinesToCurvilinearParametersNative,
            const Polygons& polygon)
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

                firstCrossing = std::max(0.0, std::min(firstCrossing, double(m_numSplineNodes[first])));
                secondCrossing = std::max(0.0, std::min(secondCrossing, double(m_numSplineNodes[second])));

                double firstLeft = std::max(0.0,std::min(double(m_numSplineNodes[first]-1), firstCrossing - firstRatioIterations / 2.0));
                double firstRight = std::max(0.0, std::min(double(m_numSplineNodes[first]-1), firstCrossing + firstRatioIterations / 2.0));

                double secondLeft = std::max(0.0, std::min(double(m_numSplineNodes[second]-1), secondCrossing - secondRatioIterations / 2.0));
                double secondRight = std::max(0.0, std::min(double(m_numSplineNodes[second]-1), secondCrossing + secondRatioIterations / 2.0));

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
            std::vector<bool>& isLeftOriented,
            std::vector<double>& adimensionalCrossPoints,
            std::vector<double>& crossingAngles)
        {
            numCrossedSplines = 0;
            intersectedSplines.resize(m_numSplines, -1);
            isLeftOriented.resize(m_numSplines, false);
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
                    isLeftOriented[s] = true;
                    if (crossProductIntersection > 0.0)
                    {
                        isLeftOriented[s] = false;
                    }
                    adimensionalCrossPoints[s] = firstSplineRatio;
                    crossingAngles[s] = crossProductIntersection;
                }
            }

            const auto sortedIndexses = SortedIndexes(adimensionalCrossPoints);
            ReorderVector(intersectedSplines, sortedIndexses);
            ReorderVector(isLeftOriented, sortedIndexses);
            ReorderVector(crossingAngles, sortedIndexses);

            return true;
        }


        //make_wholegridline
        bool MakeAllGridLines()
        {

            int numCenterSplines = 0;
            for (int s = 0; s < m_numSplines; ++s)
            {
                //center splines only
                if (m_numLayers[s] != 0)
                {
                    continue;
                }
                numCenterSplines += 1;
            }

            if (numCenterSplines == 0)
            {
                return false;
            }

            std::vector<Point> gridLine(numCenterSplines*(m_maximumNumMeshNodesAlongM + 2));

            // allocate
            int gridLineIndex = 0;
            for (int s = 0; s < m_numSplines; ++s)
            {
                //center splines only
                if (m_numLayers[s] != 0)
                {
                    continue;
                }

                if (gridLineIndex > 0)
                {
                    gridLine[gridLineIndex] = { doubleMissingValue,doubleMissingValue };
                }


                m_leftGridLineIndex[s] = gridLineIndex;
                int numMLayers = 0;

                // add...


                gridLineIndex += numMLayers;
                m_leftGridLineIndex[s] = gridLineIndex;
                gridLineIndex++;
            }

            return false;
        }

        /// make_gridline, generate a gridline on a spline with a prescribed maximum mesh width
        /// generate a gridline on a spline with a prescribed maximum mesh width
        bool MakeGridLine(const int splineIndex,
            const double maximumMeshWidth,
            const double maximumMeshHeight,
            const bool isSpacingCurvatureAdapeted,
            std::vector<Point>& gridLine,
            std::vector<Point>& gridLineOnSplineCoordinate)
        {

            // Compute spline derivatives
            std::vector<Point> coordinatesDerivatives;
            bool successful = SecondOrderDerivative(m_splines[splineIndex], m_numSplineNodes[splineIndex], coordinatesDerivatives);
            
            double curentMaxWidth = std::numeric_limits<double>::max();
            double startEndIndexses[]{ 0,m_numSplineNodes[splineIndex] - 1 };
            double splineLength = GetSplineLength(splineIndex, 0.0, m_numSplineNodes[splineIndex] - 1);
            int numMeshNodesAlongM = 1 + std::floor(splineLength / m_averageMeshWidth);
            
            std::vector<Point> distancesNewSpline;
            std::vector<Point> coordinatesDerivativeNewSpline;
            bool success = ComputeSplineAlongDistance(splineIndex, isSpacingCurvatureAdapeted, distancesNewSpline, coordinatesDerivativeNewSpline);

            while (curentMaxWidth < maximumMeshWidth)
            {
                int maximumNumAlongM = numMeshNodesAlongM + 1;
                std::vector<double> ssq(numMeshNodesAlongM + 1);
                for (int n = 0; n < numMeshNodesAlongM + 1; ++n)
                {
                    // could be the coordinate derivative instead!
                    ssq[n] = distancesNewSpline[0].x + (distancesNewSpline[1].x - distancesNewSpline[1].x) * double(n) / double(numMeshNodesAlongM);
                }

                bool isMonotonic = false;
                while(!isMonotonic)
                {
                    isMonotonic = true;
                    for (int n = 0; n < maximumNumAlongM; ++n)
                    {
                        if(ssq[n + 1] - ssq[n]< 0.0)
                        {
                            isMonotonic = false;
                            break;
                        }
                    }
                    if ( !isMonotonic )
                    {
                        for (int n = 1; n < maximumNumAlongM - 1; ++n)
                        {
                            ssq[n] = 0.5 *(ssq[n - 1] + ssq[n + 1]);
                        }
                    }
                }
               //mc = numGridNodes
               //num = numSpine
               //makespl(startstop, xsp, ysp, max(mc, num), num, 2, mc - 1, xc, yc, kmax, sc, h)
               //MAKESPL(T,           X,   Y, imax,           N, NT, MNFAC, XH, YH, KMAX, TT, H)
                  
            }

            return false;

        }

        ///MAKES
        bool ComputeSplineAlongDistance(
            const int splineIndex, 
            const bool isSpacingCurvatureAdapted,
            std::vector<Point>& distances,
            std::vector<Point>& coordinatesDerivative)
        {
            // Compute new spline, only in the x-direction
            distances.resize(m_numSplineNodes[splineIndex],{0,0});
            for (int i = 0; i < m_numSplineNodes[splineIndex]; ++i)
            {
                distances[i].x = GetSplineLength(splineIndex, 0.0, i, 10, isSpacingCurvatureAdapted);
            }
            bool success = SecondOrderDerivative(distances, m_numSplineNodes[splineIndex], coordinatesDerivative);
            return success;
        }

        ///get_splineprops
        bool ComputeSplineProperties()
        {
            // cross spline is a spline with only two points, others are non-cross splines
            AllocateSplineProperties(m_numSplines);
            for (int s = 0; s < m_numSplines; ++s)
            {
                m_splinesLength[s] = GetSplineLength(s, 0, m_numSplineNodes[s] - 1);
                bool successful = GetSplineIntersections(s, m_projection, m_numIntersectingSplines[s],
                    m_intersectingSplinesIndexses[s], m_isLeftOriented[s], m_centerSplineCrossings[s], m_cosCrossingAngle[s]);

                if(!successful)
                {
                    return false;
                }
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
                // only crossing splines with one or more center spline
                if (m_numSplineNodes[s] != 2 || m_numIntersectingSplines[s] < 1)
                {
                    continue;
                }

                int middleSplineIndexInCenterSplines = std::min(m_numIntersectingSplines[s] / 2, m_numIntersectingSplines[s]);
                int middleSplineIndex = m_intersectingSplinesIndexses[s][middleSplineIndexInCenterSplines];

                // if m_numIntersectingSplines[s] is even, check if the middle spline has already been assigned as a bounding spline
                if (m_numLayers[middleSplineIndex] != 0 && 2 * middleSplineIndex == m_numIntersectingSplines[s])
                {
                    middleSplineIndexInCenterSplines = std::min(middleSplineIndexInCenterSplines + 1, m_numIntersectingSplines[s]);
                    middleSplineIndex = m_intersectingSplinesIndexses[s][middleSplineIndexInCenterSplines];
                }
                
                if (m_numLayers[middleSplineIndex] == 0)
                {
                    // associate bounding splines with the middle spline
                    for (int i = 0; i < middleSplineIndexInCenterSplines; ++i)
                    {
                        int index = m_intersectingSplinesIndexses[s][i];
                        m_numLayers[index] = -middleSplineIndex;

                    }
                    for (int i = middleSplineIndexInCenterSplines + 1; i < m_numIntersectingSplines[s]; ++i)
                    {
                        int index = m_intersectingSplinesIndexses[s][i];
                        m_numLayers[index] = -middleSplineIndex;
                    }
                }
            }

            bool successfull = ComputeHeights();
            return successfull;
        }

        // get_heights
        // get the grid heights from the cross spline information
        // get_heights
        bool ComputeHeights()
        {
            for (int i = 0; i < m_numSplines; ++i)
            {
                // Heights should be computed only for center splines
                if (m_numSplineNodes[i] <= 2)
                {
                    continue;
                }
                for (int j = 0; j < m_numIntersectingSplines[i]; ++j)
                {
                    int intersectingSplineIndex = m_intersectingSplinesIndexses[i][j];
                    bool success = ComputeSubHeights(i, j);
                    if (!success)
                    {
                        return false;
                    }
                }
            }

            // compute m_maximumGridHeight
            for (int s = 0; s < m_numSplines; ++s)
            {
                int numIntersectingSplines = m_numIntersectingSplines[s];
                double splineLength = m_splinesLength[s];
                if (numIntersectingSplines == 0)
                {
                    m_maximumGridHeight[s] = m_aspectRatioFirstLayer * splineLength;
                    continue;
                }
                double maximumHeight = 0.0;
                for (int c = 0; c < numIntersectingSplines; ++c)
                {
                    double sumLeftHeights = 0.0;
                    for (int ss = 0; ss < m_numHeightsCrossSplineLeft[s][c]; ++ss)
                    {
                        sumLeftHeights += m_crossSplineLeftHeights[s][c][ss];
                    }
                    double sumRightHeights = 0.0;
                    for (int ss = 0; ss < m_numHeightsCrossSplineRight[s][c]; ++ss)
                    {
                        sumRightHeights += m_crossSplineRightHeights[s][c][ss];
                    }
                    maximumHeight = std::max(maximumHeight, std::max(sumLeftHeights, sumRightHeights));
                }

                m_maximumGridHeight[s] = maximumHeight;
            }
            return true;
        }

        ///comp_subheights, compute the lengths of a cross spline (left and right of the spline)
        bool ComputeSubHeights(const int centerSplineLocalIndex, const int intersectingSplineLocalIndex)
        {
            // find center spline index
            int centerSplineLocalIndexInIntersectingSpline = 0;
            int intersectingSplineGlobalIndex = m_intersectingSplinesIndexses[centerSplineLocalIndex][intersectingSplineLocalIndex]; //js
            //m_intersectingSplinesIndexses[intersectingSplinesIndex] // ics
            for (int s = 0; s < m_numIntersectingSplines[intersectingSplineGlobalIndex]; ++s)
            {
                if (m_intersectingSplinesIndexses[intersectingSplineGlobalIndex][s]== centerSplineLocalIndex)
                {
                    centerSplineLocalIndexInIntersectingSpline = s;
                    break;
                }
            }

            // right part
            int numSubIntervalsRight = 0;
            int rightCenterSplineIndex = centerSplineLocalIndexInIntersectingSpline;
            int leftCenterSplineIndex;
            m_crossSplineRightHeights[centerSplineLocalIndex][intersectingSplineLocalIndex].resize(m_maxNumHeights, 0);
            for (int s = centerSplineLocalIndexInIntersectingSpline; s <  m_numIntersectingSplines[intersectingSplineGlobalIndex] - 1; ++s)
            {
                if (numSubIntervalsRight >= m_maxNumHeights)
                {
                    break;
                }
                if (m_numLayers[m_intersectingSplinesIndexses[intersectingSplineGlobalIndex][s + 1]] != -centerSplineLocalIndex)
                {
                    continue;
                }
                leftCenterSplineIndex = rightCenterSplineIndex;
                rightCenterSplineIndex = s + 1;
                m_crossSplineRightHeights[centerSplineLocalIndex][intersectingSplineLocalIndex][numSubIntervalsRight] = GetSplineLength(intersectingSplineGlobalIndex,
                    m_centerSplineCrossings[intersectingSplineGlobalIndex][leftCenterSplineIndex], m_centerSplineCrossings[intersectingSplineGlobalIndex][rightCenterSplineIndex]);
                numSubIntervalsRight++;
            }

            m_crossSplineRightHeights[centerSplineLocalIndex][intersectingSplineLocalIndex][numSubIntervalsRight] = GetSplineLength(intersectingSplineGlobalIndex,
                m_centerSplineCrossings[intersectingSplineGlobalIndex][rightCenterSplineIndex], m_numSplineNodes[intersectingSplineGlobalIndex] - 1);

            numSubIntervalsRight++;
            std::fill(m_crossSplineRightHeights[centerSplineLocalIndex][intersectingSplineLocalIndex].begin() + numSubIntervalsRight, m_crossSplineRightHeights[centerSplineLocalIndex][intersectingSplineLocalIndex].end(), 0.0);
            m_numHeightsCrossSplineRight[centerSplineLocalIndex][intersectingSplineLocalIndex] = numSubIntervalsRight;

            // left part
            int numSubIntervalsLeft = 0;
            leftCenterSplineIndex = centerSplineLocalIndexInIntersectingSpline;
            m_crossSplineLeftHeights[centerSplineLocalIndex][intersectingSplineLocalIndex].resize(m_maxNumHeights, 0);
            for (int s = centerSplineLocalIndexInIntersectingSpline; s >= 1; --s)
            {
                if (numSubIntervalsLeft >= m_maxNumHeights)
                {
                    break;
                }
                if (m_numLayers[m_intersectingSplinesIndexses[intersectingSplineGlobalIndex][s - 1]] != -centerSplineLocalIndex)
                {
                    continue;
                }
                rightCenterSplineIndex = leftCenterSplineIndex;
                leftCenterSplineIndex = s - 1;
                m_crossSplineLeftHeights[centerSplineLocalIndex][intersectingSplineLocalIndex][numSubIntervalsLeft] = GetSplineLength(intersectingSplineGlobalIndex,
                    m_centerSplineCrossings[intersectingSplineGlobalIndex][leftCenterSplineIndex], m_centerSplineCrossings[intersectingSplineGlobalIndex][rightCenterSplineIndex]);
                numSubIntervalsLeft++;
            }

            m_crossSplineLeftHeights[centerSplineLocalIndex][intersectingSplineLocalIndex][numSubIntervalsLeft] = GetSplineLength(intersectingSplineGlobalIndex,
                0.0, m_centerSplineCrossings[intersectingSplineGlobalIndex][leftCenterSplineIndex]);

            numSubIntervalsLeft++;
            std::fill(m_crossSplineLeftHeights[centerSplineLocalIndex][intersectingSplineLocalIndex].begin() + numSubIntervalsLeft, m_crossSplineLeftHeights[centerSplineLocalIndex][intersectingSplineLocalIndex].end(), 0.0);
            m_numHeightsCrossSplineLeft[centerSplineLocalIndex][intersectingSplineLocalIndex] = numSubIntervalsLeft;

            // if not left oriented, swap
            if (!m_isLeftOriented[centerSplineLocalIndex][centerSplineLocalIndexInIntersectingSpline])
            {
                m_numHeightsCrossSplineLeft[centerSplineLocalIndex][intersectingSplineLocalIndex] = numSubIntervalsRight;
                m_numHeightsCrossSplineRight[centerSplineLocalIndex][intersectingSplineLocalIndex] = numSubIntervalsLeft;

                std::vector<double> leftSubIntervalsTemp(m_crossSplineLeftHeights[centerSplineLocalIndex][intersectingSplineLocalIndex]);
                m_crossSplineLeftHeights[centerSplineLocalIndex][intersectingSplineLocalIndex] = m_crossSplineRightHeights[centerSplineLocalIndex][intersectingSplineLocalIndex];
                m_crossSplineRightHeights[centerSplineLocalIndex][intersectingSplineLocalIndex] = leftSubIntervalsTemp;
            }

            return true;
        }

        ///splinelength
        double GetSplineLength(int index, double beginFactor, double endFactor, int numSamples = 100, bool accountForCurvature = false )
        {
            std::vector<Point> coordinatesDerivatives(m_numSplineNodes[index]);
            SecondOrderDerivative(m_splines[index], m_numSplineNodes[index], coordinatesDerivatives);

            double delta = 1.0 / numSamples;
            double numPoints = std::max(std::floor( 0.9999 + (endFactor - beginFactor) / delta), 10.0);
            delta = (endFactor - beginFactor) / numPoints;

            // first point
            Point leftPoint;
            Interpolate(m_splines[index], coordinatesDerivatives, beginFactor, leftPoint);

            double splineLength = 0.0;

            double pointAdimensionalCoordinate = beginFactor;
            for (int p = 0; p < numPoints; ++p)
            {
                pointAdimensionalCoordinate += delta;
                Point rightPoint;
                Interpolate(m_splines[index], coordinatesDerivatives, pointAdimensionalCoordinate, rightPoint);
                double curvatureFactor = 1.0;
                if (accountForCurvature)
                {
                    //TODO, compute curvature factor using comp_curv
                }
                splineLength = splineLength + Distance(leftPoint, rightPoint, m_projection) * curvatureFactor;
                leftPoint = rightPoint;
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
            coordinatesDerivatives.resize(coordinates.size());
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
        double m_aspectRatioFirstLayer;
        int m_maximumNumMeshNodesAlongM = 2000;
        double m_averageMeshWidth = 500.0;

        // Spline properties (first index is the spline number)
        const int m_maxNumHeights = 10;                                                 // Nsubmax maximum number of subintervals of grid layers, each having their own exponential grow factor
        std::vector<int> m_numLayers;                                                   // id number of layers ( >0 only for center spline)
        std::vector<int> m_numIntersectingSplines;                                      // ncs num of cross splines
        std::vector<std::vector<int>> m_intersectingSplinesIndexses;                           // ics for each cross spline, the indexses of the center splines
        std::vector<double> m_splinesLength;                                            // splinesLength spline path length
        std::vector<double> m_maximumGridHeight;                                        // hmax maximum grid height
        std::vector<std::vector<bool>> m_isLeftOriented;                                // isLeftOriented cross spline is left to right(.true.) or not (.false.) w.r.t.center spline
        std::vector<std::vector<double>>  m_centerSplineCrossings;                      // t center spline coordinates of cross splines
        std::vector<std::vector<double>> m_cosCrossingAngle;                            // cosPhi cosine of crossing angle
        std::vector<std::vector<std::vector<double>>> m_crossSplineLeftHeights;         // hL left - hand side grid heights at cross spline locations for each grid layer subinterval, hL(1, :) being the height of the first subinterval, etc.
        std::vector<std::vector<std::vector<double>>> m_crossSplineRightHeights;        // hR right - hand side grid heights at cross spline locations for each grid layer subinterval, hR(1, :) being the height of the first subinterval, etc.
        std::vector<std::vector<int>>   m_numHeightsCrossSplineLeft;                    // NsubL number of subintervals of grid layers at cross spline locations at the left - hand side of the spline, each having their own exponential grow factor
        std::vector<std::vector<int>>  m_numHeightsCrossSplineRight;                    // NsubR number of subintervals of grid layers at cross spline locations at the right - hand side of the spline, each having their own exponential grow factor
        std::vector<int> m_mfac;                                                        // mfac number of grid intervals on the spline
        std::vector<std::vector<int>> m_nfacL;                                          // nfacL number of grid layers in each subinterval at the left - hand side of the spline * not used yet*
        std::vector<std::vector<int>> m_nfacR;                                          // nfacR number of grid layers in each subinterval at the right - hand side of the spline * not used yet*
        std::vector<int> m_leftGridLineIndex;                                           // iL index in the whole gridline array of the first grid point on the left - hand side of the spline
        std::vector<int> m_rightGridLineIndex;                                          // iR index in the whole gridline array of the first grid point on the right - hand side of the spline
    };

}