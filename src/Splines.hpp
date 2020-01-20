#pragma once

#include <vector>
#include <algorithm>
#include <functional>
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

        Splines(const Splines& other)
            : m_numSplines(other.m_numSplines),
              m_numAllocatedSplines(other.m_numAllocatedSplines),
              m_numSplineNodes(other.m_numSplineNodes),
              m_numAllocatedSplineNodes(other.m_numAllocatedSplineNodes),
              m_splineCornerPoints(other.m_splineCornerPoints),
              m_projection(other.m_projection),
              m_aspectRatioFirstLayer(other.m_aspectRatioFirstLayer),
              m_maximumNumMeshNodesAlongM(other.m_maximumNumMeshNodesAlongM),
              m_averageMeshWidth(other.m_averageMeshWidth),
              m_maxNumHeights(other.m_maxNumHeights),
              m_numLayers(other.m_numLayers),
              m_numIntersectingSplines(other.m_numIntersectingSplines),
              m_intersectingSplinesIndexses(other.m_intersectingSplinesIndexses),
              m_splinesLength(other.m_splinesLength),
              m_maximumGridHeight(other.m_maximumGridHeight),
              m_isLeftOriented(other.m_isLeftOriented),
              m_centerSplineCrossings(other.m_centerSplineCrossings),
              m_cosCrossingAngle(other.m_cosCrossingAngle),
              m_crossSplineLeftHeights(other.m_crossSplineLeftHeights),
              m_crossSplineRightHeights(other.m_crossSplineRightHeights),
              m_numHeightsCrossSplineLeft(other.m_numHeightsCrossSplineLeft),
              m_numHeightsCrossSplineRight(other.m_numHeightsCrossSplineRight),
              m_mfac(other.m_mfac),
              m_nfacL(other.m_nfacL),
              m_nfacR(other.m_nfacR),
              m_leftGridLineIndex(other.m_leftGridLineIndex),
              m_rightGridLineIndex(other.m_rightGridLineIndex)
        {
        }

        Splines& operator=(const Splines& other)
        {
            if (this == &other)
                return *this;
            m_numSplines = other.m_numSplines;
            m_numAllocatedSplines = other.m_numAllocatedSplines;
            m_numSplineNodes = other.m_numSplineNodes;
            m_numAllocatedSplineNodes = other.m_numAllocatedSplineNodes;
            m_splineCornerPoints = other.m_splineCornerPoints;
            m_projection = other.m_projection;
            m_aspectRatioFirstLayer = other.m_aspectRatioFirstLayer;
            m_maximumNumMeshNodesAlongM = other.m_maximumNumMeshNodesAlongM;
            m_averageMeshWidth = other.m_averageMeshWidth;
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

        Splines() : m_projection(Projections::cartesian), m_numAllocatedSplines(5), m_numSplines(0)
        {
        }

        Splines(Projections projection) : m_projection(projection), m_numAllocatedSplines(0), m_numSplines(0)
        {
            AllocateSplines(5);
            m_aspectRatioFirstLayer = 0.10;
        }

        /// add a new spline
        bool Set(const std::vector<Point>& splines)
        {
            AllocateSplines(m_numSplines + 1);
            m_splineCornerPoints[m_numSplines] = splines;
            m_numSplineNodes[m_numSplines] = splines.size();

            // compute basic properties
            SecondOrderDerivative(m_splineCornerPoints[m_numSplines], m_numSplineNodes[m_numSplines], m_splineDerivatives[m_numSplines]);
            m_splinesLength[m_numSplines] = GetSplineLength(m_numSplines, 0, m_numSplineNodes[m_numSplines] - 1);
            m_numSplines++;

            return true;
        }

        bool AllocateSplines(int numSplines)
        {
            AllocateVector(m_numAllocatedSplines, m_splineCornerPoints, std::vector<Point>(10, { doubleMissingValue, doubleMissingValue }), numSplines);
            m_numAllocatedSplineNodes.resize(m_numAllocatedSplines, 10);
            m_numSplineNodes.resize(m_numAllocatedSplines, 0);

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
            m_splineDerivatives.resize(m_numAllocatedSplines);

            return true;
        }



        bool DeleteSpline(int position)
        {
            m_splineCornerPoints.erase(m_splineCornerPoints.begin() + position);
            m_numSplineNodes.erase(m_numSplineNodes.begin() + position);
            m_numSplines--;
            return true;
        }

        /// add a new spline point in an existing spline
        bool Set(int splineIndex, const Point& point)
        {
            if (splineIndex >= m_numSplines)
            {
                return false;
            }
            AllocateVector(m_numAllocatedSplineNodes[splineIndex], m_splineCornerPoints[splineIndex], { doubleMissingValue, doubleMissingValue }, 10);

            m_splineCornerPoints[splineIndex][m_numSplineNodes[splineIndex]] = point;
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
                    bool isInPolygons = IsPointInPolygons(m_splineCornerPoints[s][n], polygon.m_nodes, polygon.m_numNodes);
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
                    bool areCrossing = GridGeom::AreLinesCrossing(m_splineCornerPoints[first][n],
                        m_splineCornerPoints[first][n + 1],
                        m_splineCornerPoints[second][nn],
                        m_splineCornerPoints[second][nn + 1],
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
                Interpolate(m_splineCornerPoints[first], m_splineDerivatives[first], firstLeft, firstLeftSplinePoint);
                Point firstRightSplinePoint;
                Interpolate(m_splineCornerPoints[first], m_splineDerivatives[first], firstRight, firstRightSplinePoint);

                Point secondLeftSplinePoint;
                Interpolate(m_splineCornerPoints[second], m_splineDerivatives[second], secondLeft, secondLeftSplinePoint);
                Point secondRightSplinePoint;
                Interpolate(m_splineCornerPoints[second], m_splineDerivatives[second], secondRight, secondRightSplinePoint);

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
        bool MakeAllGridLines(const bool isSpacingCurvatureAdapeted)
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
            std::vector<double> gridLineOnSplineCoordinate(numCenterSplines*(m_maximumNumMeshNodesAlongM + 2));

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
                    gridLineIndex++;
                }


                m_leftGridLineIndex[s] = gridLineIndex;
                int numMLayers = 0;

                // add...
                bool success = MakeGridLine(s, isSpacingCurvatureAdapeted, gridLineIndex, gridLine, gridLineOnSplineCoordinate);


                gridLineIndex += numMLayers;
                m_leftGridLineIndex[s] = gridLineIndex;
                gridLineIndex++;
            }

            return false;
        }

        /// make_gridline, generate a gridline on a spline with a prescribed maximum mesh width
        /// generate a gridline on a spline with a prescribed maximum mesh width
        /// call makespl(startstop, xsp, ysp, max(mc,num), num, 2, mc-1, xc, yc, kmax, sc, h)
        ///     MAKESPL(        T,   X,   Y,        imax,   N,NT,MNFAC,  XH,YH, KMAX, TT, H)
        bool MakeGridLine(const int splineIndex,
            const bool isSpacingCurvatureAdapeted,
            const int startingIndex,
            std::vector<Point>& gridLine,
            std::vector<double>& centerLineCoordinates)
        {
            int numNodesAlongM = 1 + std::floor(m_splinesLength[splineIndex] / m_averageMeshWidth);
            numNodesAlongM = std::min(numNodesAlongM, m_maximumNumMeshNodesAlongM);

            double endSplineAdimensionalCoordinate = m_numSplineNodes[splineIndex] - 1;
            double splineLength = GetSplineLength(splineIndex, 0.0, endSplineAdimensionalCoordinate, 10, isSpacingCurvatureAdapeted, m_maximumGridHeight[splineIndex]);
            
            gridLine[startingIndex] = m_splineCornerPoints[splineIndex][0];
            FuncDimensionalToAdimensionalDistance func(*this, splineIndex, isSpacingCurvatureAdapeted, m_maximumGridHeight[splineIndex]);
            
            double currentMaxWidth = std::numeric_limits<double>::max();
            while (currentMaxWidth > m_averageMeshWidth && numNodesAlongM < m_maximumNumMeshNodesAlongM)
            {
                currentMaxWidth = 0.0;
                for (int n = 1; n < numNodesAlongM + 1; ++n)
                {
                    int index = startingIndex + n;
                    centerLineCoordinates[index] = splineLength * double(n) / double(numNodesAlongM);
                    func.SetDimensionalDistance(centerLineCoordinates[index]);
                    double pointAdimensionalCoordinate = FindFunctionRootWithGoldenSearch(func, 0, endSplineAdimensionalCoordinate); //1.09068327269294;
                    Interpolate(m_splineCornerPoints[splineIndex], m_splineDerivatives[splineIndex], pointAdimensionalCoordinate, gridLine[index]);
                    currentMaxWidth = std::max(currentMaxWidth, Distance(gridLine[index - 1], gridLine[index], m_projection));
                }

                // room for sub-division
                if (currentMaxWidth > m_averageMeshWidth)
                {
                    numNodesAlongM = std::min(std::max(int(m_maximumNumMeshNodesAlongM / m_maximumGridHeight[splineIndex] *numNodesAlongM), numNodesAlongM + 1), m_maximumNumMeshNodesAlongM);
                }
            }
            return true;
        }

        // functor wrapping GetSplineLength to use generic root finders
        struct FuncDimensionalToAdimensionalDistance
        {
            FuncDimensionalToAdimensionalDistance(Splines& spline, int splineIndex, bool isSpacingCurvatureAdapted, double h) :
                m_spline(spline),
                m_splineIndex(splineIndex),
                m_isSpacingCurvatureAdapted(isSpacingCurvatureAdapted),
                m_h(h)
            {
            };

            void SetDimensionalDistance(const double distance)
            {
                m_DimensionalDistance = distance;
            }

            // this is the fancion we want to find the root
            double operator()(double const& adimensionalDistancereferencePoint)
            {
                double distanceFromReferencePoint = m_spline.GetSplineLength(m_splineIndex, 0, adimensionalDistancereferencePoint, m_numSamples, m_isSpacingCurvatureAdapted, m_h);
                distanceFromReferencePoint = std::abs(distanceFromReferencePoint - m_DimensionalDistance);
                return distanceFromReferencePoint;
            }

            Splines& m_spline;
            int m_splineIndex;
            bool m_isSpacingCurvatureAdapted;
            double m_h;
            int m_numSamples = 10;
            double m_DimensionalDistance;
        };

        ///get_splineprops
        bool ComputeSplineProperties()
        {
            // cross spline is a spline with only two points, others are non-cross splines
            for (int s = 0; s < m_numSplines; ++s)
            {
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
        double GetSplineLength(int index, 
            double beginFactor, 
            double endFactor, 
            int numSamples = 100, 
            bool accountForCurvature = false, 
            double height = 1.0 )
        {

            double delta = 1.0 / numSamples;
            double numPoints = std::max(std::floor( 0.9999 + (endFactor - beginFactor) / delta), 10.0);
            delta = (endFactor - beginFactor) / numPoints;

            // first point
            Point leftPoint;
            Interpolate(m_splineCornerPoints[index], m_splineDerivatives[index], beginFactor, leftPoint);

            double splineLength = 0.0;

            double rightPointCoordinateOnSpline = beginFactor;
            double leftPointCoordinateOnSpline;
            for (int p = 0; p < numPoints; ++p)
            {
                leftPointCoordinateOnSpline = rightPointCoordinateOnSpline;
                rightPointCoordinateOnSpline += delta;
                Point rightPoint;
                Interpolate(m_splineCornerPoints[index], m_splineDerivatives[index], rightPointCoordinateOnSpline, rightPoint);
                double curvatureFactor = 0.0;
                if (accountForCurvature)
                {
                    Point normalVector;
                    Point tangentialVector;
                    ComputeCurvatureOnSplinePoint(index, 0.5*(rightPointCoordinateOnSpline + leftPointCoordinateOnSpline), curvatureFactor, normalVector, tangentialVector);
                }
                splineLength = splineLength + Distance(leftPoint, rightPoint, m_projection) * (1.0 + curvatureFactor * height);
                leftPoint = rightPoint;
            }

            return splineLength;
        }

        /// comp_curv
        /// compute curvature in a point on a spline
        bool ComputeCurvatureOnSplinePoint(
            const int splineIndex,
            const double adimensionalPointCoordinate,
            double& curvatureFactor,
            Point& normalVector,
            Point& tangentialVector)
        {
            double leftCornerPoint = std::max(std::min(double(std::floor(adimensionalPointCoordinate)), double(m_numSplineNodes[splineIndex] - 1)), 0.0);
            double rightCornerPoint = std::max(double(leftCornerPoint + 1.0), 0.0);

            double leftSegment = rightCornerPoint - adimensionalPointCoordinate;
            double rightSegment = adimensionalPointCoordinate - leftCornerPoint;
       
            Point pointCoordinate;
            Interpolate(m_splineCornerPoints[splineIndex], m_splineDerivatives[splineIndex], adimensionalPointCoordinate, pointCoordinate);

            Point p = m_splineCornerPoints[splineIndex][rightCornerPoint] - m_splineCornerPoints[splineIndex][leftCornerPoint] 
                + (m_splineDerivatives[splineIndex][leftCornerPoint] * (-3.0 *leftSegment*leftSegment + 1.0)
                    + m_splineDerivatives[splineIndex][rightCornerPoint] * (3.0 *rightSegment*rightSegment - 1.0)) / 6.0;

            Point pp = m_splineDerivatives[splineIndex][leftCornerPoint] * leftSegment + m_splineDerivatives[splineIndex][rightCornerPoint] * rightSegment;

            if (m_projection == Projections::spherical)
            {
                p.TransformToSpherical();
                pp.TransformToSpherical();
            }

            curvatureFactor = std::abs(pp.x*p.y - pp.y*p.x) / std::pow((p.x*p.x + p.y*p.y + 1e-8), 1.5);

            Point incremenetedPointCoordinate = pointCoordinate + p * 1e-4;
            NormalVectorOutside(pointCoordinate, incremenetedPointCoordinate, normalVector, m_projection);

            double distance = Distance(pointCoordinate, incremenetedPointCoordinate, m_projection);
            double dx = getDx(pointCoordinate, incremenetedPointCoordinate, m_projection);
            double dy = getDy(pointCoordinate, incremenetedPointCoordinate, m_projection);

            tangentialVector.x = dx / distance;
            tangentialVector.y = dy / distance;

            return true;
        }


        /// these functions are used in the gridgeom api, made them static
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
        std::vector<std::vector<Point>> m_splineCornerPoints;
        std::vector<std::vector<Point>> m_splineDerivatives;
        Projections m_projection;
        double m_aspectRatioFirstLayer;
        int m_maximumNumMeshNodesAlongM = 20;
        double m_averageMeshWidth = 500.0;

        // Spline properties (first index is the spline number)
        const int m_maxNumHeights = 10;                                                 // Nsubmax maximum number of subintervals of grid layers, each having their own exponential grow factor
        std::vector<int> m_numLayers;                                                   // id number of layers ( >0 only for center spline)
        std::vector<int> m_numIntersectingSplines;                                      // ncs num of cross splines
        std::vector<std::vector<int>> m_intersectingSplinesIndexses;                    // ics for each cross spline, the indexses of the center splines
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