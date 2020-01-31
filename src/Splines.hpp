#pragma once

#include <vector>
#include <algorithm>
#include <functional>
#include "Operations.cpp"
#include "Entities.hpp"
#include "Polygons.hpp"
#include "CurvilinearParametersNative.hpp"
#include "SplinesToCurvilinearParametersNative.hpp"
#include <cassert>

namespace GridGeom
{
    class Splines
    {
    public:
        Splines() : m_projection(Projections::cartesian)
        {
        };

        Splines(const Projections projection, Polygons& polygon) : m_projection(projection) ,m_polygon(polygon)
        {
        };

        /// add a new spline, return the index
        bool AddSpline(const std::vector<Point>& splines)
        {
            AllocateVector(m_numSplines + 1, m_splineCornerPoints, m_allocationSize, std::vector<Point>(10, { doubleMissingValue, doubleMissingValue }));
           
            m_numAllocatedSplines = m_splineCornerPoints.size();
            m_numAllocatedSplineNodes.resize(m_numAllocatedSplines, 10);

            m_numSplineNodes.resize(m_numAllocatedSplines, 0);
            m_numSplineNodes[m_numSplines] = splines.size();

            m_splineDerivatives.resize(m_numAllocatedSplines);
            for (int n = 0; n < m_numSplineNodes[m_numSplines]; ++n)
            {
                m_splineCornerPoints[m_numSplines][n] = splines[n];
            }
            

            m_splinesLength.resize(m_numAllocatedSplines);
            m_numLayers.resize(m_numAllocatedSplines);

            // compute basic properties
            SecondOrderDerivative(m_splineCornerPoints[m_numSplines], m_numSplineNodes[m_numSplines], m_splineDerivatives[m_numSplines]);
            m_splinesLength[m_numSplines] = GetSplineLength(m_numSplines, 0, m_numSplineNodes[m_numSplines] - 1);
            m_numSplines++;

            return true;
        }

        //Run-time settings
        bool SetParameters(const GridGeomApi::CurvilinearParametersNative& curvilinearParametersNative,
                           const GridGeomApi::SplinesToCurvilinearParametersNative& splinesToCurvilinearParametersNative)
        {
            return true;
        }

        // set the reference to polygon
        bool SetPolygon(const Polygons& polygon)
        {
            m_polygon = polygon;
            return true;
        }

        bool DeleteSpline(const int position)
        {
            m_splineCornerPoints.erase(m_splineCornerPoints.begin() + position);
            m_numSplineNodes.erase(m_numSplineNodes.begin() + position);
            m_splineDerivatives.erase(m_splineDerivatives.begin() + position);
            m_splinesLength.erase(m_splinesLength.begin() + position);
            m_numSplines--;
            return true;
        }

        /// to be called after all splines have been stored
        bool AllocateSplinesProperties()
        {
            m_numLayers.resize(m_numSplines);
            std::fill(m_numLayers.begin(), m_numLayers.end(), intMissingValue);

            m_numCrossingSplines.resize(m_numSplines, 0);
            std::fill(m_numCrossingSplines.begin(), m_numCrossingSplines.end(), 0);

            m_maximumGridHeights.resize(m_numSplines);
            std::fill(m_maximumGridHeights.begin(), m_maximumGridHeights.end(), doubleMissingValue);

            // multi-dimensional arrays
            m_crossingSplinesIndexses.resize(m_numSplines);
            m_isLeftOriented.resize(m_numSplines);
            m_crossSplineCoordinates.resize(m_numSplines);
            m_cosCrossingAngle.resize(m_numSplines);
            m_crossSplineLeftHeights.resize(m_numSplines);
            m_crossSplineRightHeights.resize(m_numSplines);
            m_numCrossSplineLeftHeights.resize(m_numSplines);
            m_numCrossSplineRightHeights.resize(m_numSplines);
            m_nfacL.resize(m_numSplines);
            m_nfacR.resize(m_numSplines);
            for (int s = 0; s < m_numSplines; ++s)
            {
                m_crossingSplinesIndexses[s].resize(m_numSplines);
                std::fill(m_crossingSplinesIndexses[s].begin(), m_crossingSplinesIndexses[s].end(), -1);

                m_isLeftOriented[s].resize(m_numSplines, true);
                std::fill(m_isLeftOriented[s].begin(), m_isLeftOriented[s].end(), true);

                m_crossSplineCoordinates[s].resize(m_numSplines);
                std::fill(m_crossSplineCoordinates[s].begin(), m_crossSplineCoordinates[s].end(), doubleMissingValue);

                m_cosCrossingAngle[s].resize(m_numSplines, doubleMissingValue);
                std::fill(m_cosCrossingAngle[s].begin(), m_cosCrossingAngle[s].end(), doubleMissingValue);

                m_crossSplineLeftHeights[s].resize(m_numSplines);
                std::fill(m_crossSplineLeftHeights[s].begin(), m_crossSplineLeftHeights[s].end(), std::vector<double>(m_maxNumCenterSplineHeights, doubleMissingValue));

                m_crossSplineRightHeights[s].resize(m_numSplines);
                std::fill(m_crossSplineRightHeights[s].begin(), m_crossSplineRightHeights[s].end(), std::vector<double>(m_maxNumCenterSplineHeights, doubleMissingValue));

                m_numCrossSplineLeftHeights[s].resize(m_numSplines);
                std::fill(m_numCrossSplineLeftHeights[s].begin(), m_numCrossSplineLeftHeights[s].end(), 0);

                m_numCrossSplineRightHeights[s].resize(m_numSplines,doubleMissingValue);
                std::fill(m_numCrossSplineRightHeights[s].begin(), m_numCrossSplineRightHeights[s].end(), 0);

                m_nfacL[s].resize(m_numSplines);
                std::fill(m_nfacL[s].begin(), m_nfacL[s].end(), 0);

                m_nfacR[s].resize(m_numSplines);
                std::fill(m_nfacR[s].begin(), m_nfacR[s].end(), 0);
            }

            m_numMSpline.resize(m_numSplines);
            std::fill(m_numMSpline.begin(), m_numMSpline.end(), 0);
            m_leftGridLineIndex.resize(m_numSplines);
            std::fill(m_leftGridLineIndex.begin(), m_leftGridLineIndex.end(), intMissingValue);
            m_rightGridLineIndex.resize(m_numSplines);
            std::fill(m_rightGridLineIndex.begin(), m_rightGridLineIndex.end(), intMissingValue);

            return true;
        }

        /// add a new spline point in an existing spline
        bool AddPointInExistingSpline(const int splineIndex, const Point& point)
        {
            if (splineIndex >= m_numSplines)
            {
                return false;
            }
            AllocateVector(m_numSplineNodes[splineIndex] + 1, m_splineCornerPoints[splineIndex], m_allocationSize, { doubleMissingValue, doubleMissingValue });
            m_numAllocatedSplineNodes[splineIndex] = m_splineCornerPoints[splineIndex].size();

            m_splineCornerPoints[splineIndex][m_numSplineNodes[splineIndex]] = point;
            m_numSplineNodes[splineIndex]++;
            return true;
        }

        /// spline2curvi
        /// 1. Eliminate spline that are not in polygon
        /// 2. Compute the properties
        /// 3. Make all grid lines of the central spline
        /// 4. Add artificial splines
        /// 5. Compute properties with artificial spline added
        /// 6. Compute edge velocities
        /// 7. Grow layers
        bool OrthogonalCurvilinearMeshFromSplines()
        {
            // no splines
            if (m_numSplines < 1)
            {
                return false;
            }

            // Delete the splines that are not fully inside the polygons
            for (int s = 0; s < m_numSplines; s++)
            {
                for (int n = 0; n < m_numSplineNodes[s]; n++)
                {
                    bool isInPolygons = IsPointInPolygons(m_splineCornerPoints[s][n], m_polygon.m_nodes, m_polygon.m_numNodes);
                    if (!isInPolygons)
                    {
                        DeleteSpline(s);
                        break;
                    }
                }
            }

            // compute properties
            bool success = ComputeSplineProperties(false);
            if (!success)
            {
                return false;
            }

            // get the properties of the center splines
            m_numM = 0;
            success = MakeAllGridLines(true);
            if (!success)
            {
                return false;
            }

            // Store original number of splines
            std::vector<Point> newCrossSpline(2);
            m_numOriginalSplines = m_numSplines;
            for (int s = 0; s < m_numOriginalSplines; ++s)
            {
                // mirrow only center splines
                if (m_numLayers[s] != 0.0)
                {
                    continue;
                }

                // construct the cross splines through the edges, along m discretization
                for (int i = m_leftGridLineIndex[s]; i < m_leftGridLineIndex[s] + m_numMSpline[s]; ++i)
                {
                    Point normal;
                    NormalVectorOutside(m_gridLine[i], m_gridLine[i + 1], normal, m_projection);

                    double xMiddle = (m_gridLine[i].x + m_gridLine[i + 1].x)*0.5;
                    double yMiddle = (m_gridLine[i].y + m_gridLine[i + 1].y)*0.5;
                    double xs1 = xMiddle + 2.0 * m_maximumGridHeights[s] * -normal.x;
                    double xs2 = xMiddle + 2.0 * m_maximumGridHeights[s] * normal.x;
                    double ys1 = yMiddle + 2.0 * m_maximumGridHeights[s] * -normal.y;
                    double ys2 = yMiddle + 2.0 * m_maximumGridHeights[s] * normal.y;

                    if (m_projection == Projections::spherical)
                    {
                        const double factor = 1.0 / (earth_radius *degrad_hp);
                        xs1 = xs1 * factor;
                        ys1 = ys1 * factor;
                        xs2 = xs2 * factor;
                        ys2 = ys2 * factor;
                    }

                    newCrossSpline[0] = { xs1, ys1 };
                    newCrossSpline[1] = { xs2, ys2 };
                    AddSpline(newCrossSpline);
                    // flag the cross spline as artificially added
                    m_numLayers[m_numSplines - 1] = 3;
                }
            }

            // Backup original spline properties
            m_leftGridLineIndexOriginal.resize(m_numOriginalSplines);
            m_rightGridLineIndexOriginal.resize(m_numOriginalSplines);
            m_mfacOriginal.resize(m_numOriginalSplines);
            m_maximumGridHeightsOriginal.resize(m_numOriginalSplines);
            m_numLayersOriginal.resize(m_numOriginalSplines);
            for (int s = 0; s < m_numOriginalSplines; ++s)
            {
                m_leftGridLineIndexOriginal[s] = m_leftGridLineIndex[s];
                m_rightGridLineIndexOriginal[s] = m_rightGridLineIndex[s];
                m_mfacOriginal[s] = m_numMSpline[s];
                m_maximumGridHeightsOriginal[s] = m_maximumGridHeights[s];
                m_numLayersOriginal[s] = m_numLayers[s];
            }

            // compute spline properties with artificial splines added
            ComputeSplineProperties(true);

            // artificial cross spline: remove the last part of the sub-intervals (since it makes no sense, 
            // as the artificial cross spline has an arbitrary, but sufficiently large, length)
            for (int s = 0; s < m_numOriginalSplines; ++s)
            {
                // Remove the last part of the sub-intervals
                if (m_numLayers[s] != 0)
                {
                    continue;
                }

                // For number of intersecting splines
                for (int i = 0; i < m_numCrossingSplines[s]; ++i)
                {
                    int crossingSplineIndex = m_crossingSplinesIndexses[s][i];
                    if (m_numLayers[crossingSplineIndex] == 3)
                    {
                        m_numCrossSplineLeftHeights[s][i] = m_numCrossSplineLeftHeights[s][i] - 1;
                        m_numCrossSplineRightHeights[s][i] = m_numCrossSplineRightHeights[s][i] - 1;
                    }
                }
            }

            // Compute edge velocity
            std::vector<double> edgeVelocities(m_numM - 1, doubleMissingValue);
            std::vector<std::vector<double>> growFactorOnSubintervalAndEdge(m_maxNumCenterSplineHeights, std::vector<double>(m_numM - 1, doubleMissingValue));
            std::vector<std::vector<int>> numPerpendicularFacesOnSubintervalAndEdge(m_maxNumCenterSplineHeights, std::vector<int>(m_numM, 0));

            success = ComputeEdgeVelocities(edgeVelocities, growFactorOnSubintervalAndEdge, numPerpendicularFacesOnSubintervalAndEdge);
            if (!success)
            {
                return false;
            }

            // Increase curvilinear grid
            int maxNumPoints = std::max(m_numM + 1, m_maxNumN + 1);
            // The layer by coordinate to grow
            std::vector<std::vector<Point>> gridPoints(m_maxNumN + 1, std::vector<Point>(m_numM + 1, { doubleMissingValue, doubleMissingValue }));
            std::vector<int> validFrontNodes(m_numM, 1);

            // Copy the first n
            for (int n = 0; n < m_numM; ++n)
            {
                gridPoints[0][n] = m_gridLine[n];
                if (m_gridLine[n].x == doubleMissingValue)
                {
                    validFrontNodes[n] = 0;
                }
                int sumLeft = 0;
                int sumRight = 0;
                int leftColumn = std::max(n - 1, 0);
                int rightColumn = std::min(n, m_numM - 2);
                for (int j = 0; j < m_maxNumCenterSplineHeights; ++j)
                {
                    sumLeft += numPerpendicularFacesOnSubintervalAndEdge[j][leftColumn];
                    sumRight += numPerpendicularFacesOnSubintervalAndEdge[j][rightColumn];
                }
                if (sumLeft == 0 && sumRight == 0)
                {
                    validFrontNodes[n] = 0;
                }
            }

            //compute maximum mesh width and get dtolLR in the proper dimension
            double maximumGridWidth = 0.0;
            for (int i = 0; i < gridPoints[0].size() - 1; i++)
            {
                if (!gridPoints[0][i].IsValid() || !gridPoints[0][i + 1].IsValid())
                {
                    continue;
                }
                maximumGridWidth = std::max(maximumGridWidth, Distance(gridPoints[0][i], gridPoints[0][i + 1], m_projection));
            }
            m_onTopOfEachOtherTolerance = m_onTopOfEachOtherTolerance * maximumGridWidth;

            // grow grid, from second n
            for (int n = 1; n < m_maxNumN + 1; ++n)
            {
                success = GrowLayer(n, validFrontNodes, edgeVelocities, gridPoints);
                if (!success)
                {
                    return false;
                }
            }

            return true;
        }

        /// growlayer
        bool GrowLayer(const int layerIndex,
            const std::vector<int>& validFrontNodes,
            const std::vector<double>& edgeVelocities,
            std::vector<std::vector<Point>>& gridPoints)
        {
            //std::vector<int> currentValidFrontNodes(validFrontNodes);

            assert(layerIndex - 1 >= 0);
            
            std::vector<Point> velocityVector(validFrontNodes.size());
            bool success = ComputeVelocities(gridPoints[layerIndex - 1], edgeVelocities, velocityVector);
            if (!success) 
            {
                return false;
            }

            std::vector<Point> activeLayerPoints(gridPoints[layerIndex - 1]);
            for (int m = 0; m < velocityVector.size(); ++m)
            {
                if (!velocityVector[m].IsValid())
                {
                    gridPoints[layerIndex - 1][m] = {doubleMissingValue, doubleMissingValue};
                    activeLayerPoints[m] = { doubleMissingValue, doubleMissingValue };
                }
            }

            // FindFront
            int numGridPoints = gridPoints.size() * gridPoints[0].size();
            std::vector<std::vector<int>> gridPointsIndexses(numGridPoints,std::vector<int>(2, - 1));
            std::vector<Point> frontGridPoints(numGridPoints);
            int numFrontPoints;
            success = FindFront(gridPoints, gridPointsIndexses, frontGridPoints, numFrontPoints);
            if (!success)
            {
                return false;
            }

            std::vector<Point> frontVelocities(numGridPoints);
            success = CopyVelocitiesToFront(layerIndex - 1, validFrontNodes, velocityVector, numFrontPoints, 
                gridPointsIndexses, frontGridPoints, frontVelocities);

            if (!success)
            {
                return false;
            }

            return true;
        }


        /// copy growth velocities to the front, and add points in the front at corners
        /// copy_vel_to_front
        bool CopyVelocitiesToFront(
            const int layerIndex,
            const std::vector<int>& validFrontNodes,
            const std::vector<Point>& previousVelocities,
            int& numFrontPoints,
            std::vector<std::vector<int>>& gridPointsIndexses,
            std::vector<Point>& frontGridPoints,
            std::vector<Point>& velocities)
        {
            int numCornerNodes = 0;
            int p = -1;
            while(p<numFrontPoints)
            {
                p = p + 1;
                if (gridPointsIndexses[p][1] == layerIndex && validFrontNodes[gridPointsIndexses[p][0]] == 1)
                {
                    velocities[p] = previousVelocities[gridPointsIndexses[p][0]];
                    if (!velocities[p].IsValid())
                    {
                        velocities[p] = { 0.0,0.0 };
                    }

                    // Check for cornernodes
                    int previous = std::max(p - 1, 0);
                    std::vector<int> previousIndexses = gridPointsIndexses[previous];
                    int next = std::min(p + 1, numFrontPoints);
                    std::vector<int> nextIndexses = gridPointsIndexses[next];

                    // Check corner nodes
                    bool ll = previousIndexses[0] == gridPointsIndexses[p][0] - 1 &&
                        previousIndexses[1] == gridPointsIndexses[p][1] &&
                        validFrontNodes[previousIndexses[0]] == 0;

                    bool lr = nextIndexses[0] == gridPointsIndexses[p][0] + 1 &&
                        nextIndexses[1] == gridPointsIndexses[p][1] &&
                        validFrontNodes[nextIndexses[0]] == 0;

                    ll = ll || previousIndexses[0] == gridPointsIndexses[p][0] && previousIndexses[1] < gridPointsIndexses[p][1];
                    lr = lr || nextIndexses[0] == gridPointsIndexses[p][0] && nextIndexses[1] < gridPointsIndexses[p][1];
                    if (ll || lr)
                    {
                        numCornerNodes++;
                        if (numFrontPoints + 1 > frontGridPoints.size())
                        {
                            continue;
                        }
                        for (int i = numFrontPoints; i >= p; --i)
                        {
                            frontGridPoints[i + 1] = frontGridPoints[i];
                            velocities[i + 1] = velocities[i];
                            gridPointsIndexses[i + 1] = gridPointsIndexses[i];
                        }
                        numFrontPoints++;

                        if (ll)
                        {
                            velocities[p] = { 0.0,0.0 };
                        }
                        else
                        {
                            velocities[p + 1] = { 0.0,0.0 };
                        }
                        p = p + 1;
                    }
                }
            }
            return true;
        }


        ///find the frontline of the old (static) grid
        ///findfront
        bool FindFront(
            const std::vector<std::vector<Point>>& allGridPoints,
            std::vector<std::vector<int>>& gridPointsIndexses,
            std::vector<Point>& frontGridPoints,
            int& numFrontPoints)
        {

            std::vector<int> frontPosition(allGridPoints[0].size(), allGridPoints.size());
            for (int m = 0; m < allGridPoints[0].size(); ++m)
            {
                for (int n = 0; n < allGridPoints.size() - 1; ++n)
                {
                    if (!allGridPoints[n][m].IsValid() || !allGridPoints[n][m + 1].IsValid())
                    {
                        frontPosition[m] = n - 1;
                        break;
                    }
                }
            }

            numFrontPoints = 0;
            // check for circular connectivity
            int currentLeftIndex;
            int currentRightIndex;
            int previousFrontPosition = 0;
            GetNeighbours(allGridPoints[0], 0, currentLeftIndex, currentRightIndex);
            if (currentLeftIndex == 0)
            {
                frontGridPoints[0] = allGridPoints[0][0];
                // store front index
                gridPointsIndexses[numFrontPoints][0] = 0;
                gridPointsIndexses[numFrontPoints][1] = 0;
                numFrontPoints++;
            }
            else
            {
                previousFrontPosition = frontPosition[currentLeftIndex];
                frontGridPoints[numFrontPoints] = allGridPoints[0][frontPosition[0]];
                gridPointsIndexses[numFrontPoints][0] = frontPosition[0];
                gridPointsIndexses[numFrontPoints][1] = 0;
                numFrontPoints++;
            }

            for (int m = 0; m < allGridPoints[0].size() - 2; ++m)
            {
                GetNeighbours(allGridPoints[0], m, currentLeftIndex, currentRightIndex);
                int currentFrontPosition = frontPosition[m];
                if (currentFrontPosition >= 0)
                {
                    if (previousFrontPosition == -1)
                    {
                        frontGridPoints[numFrontPoints] = allGridPoints[0][m];
                        gridPointsIndexses[numFrontPoints][0] = m;
                        gridPointsIndexses[numFrontPoints][1] = 0;
                        numFrontPoints++;
                    }
                    for (int i = previousFrontPosition + 1; i <= currentFrontPosition; ++i)
                    {
                        frontGridPoints[numFrontPoints] = allGridPoints[i][m];
                        gridPointsIndexses[numFrontPoints][0] = m;
                        gridPointsIndexses[numFrontPoints][1] = i;
                        numFrontPoints++;
                    }
                    for (int i = previousFrontPosition; i > currentFrontPosition; --i)
                    {
                        frontGridPoints[numFrontPoints] = allGridPoints[i][m];
                        gridPointsIndexses[numFrontPoints][0] = m;
                        gridPointsIndexses[numFrontPoints][1] = i;
                        numFrontPoints++;
                    }

                    frontGridPoints[numFrontPoints] = allGridPoints[currentFrontPosition][m + 1];
                    gridPointsIndexses[numFrontPoints][0] = m + 1;
                    gridPointsIndexses[numFrontPoints][1] = currentFrontPosition;
                    numFrontPoints++;
                }
                else if (previousFrontPosition >= 0)
                {
                    for (int i = previousFrontPosition - 1; i >= 0; --i)
                    {
                        frontGridPoints[numFrontPoints] = allGridPoints[i][m];
                        gridPointsIndexses[numFrontPoints][0] = m;
                        gridPointsIndexses[numFrontPoints][1] = i;
                        numFrontPoints++;
                    }

                    frontGridPoints[numFrontPoints] = { doubleMissingValue,doubleMissingValue };
                    gridPointsIndexses[numFrontPoints][0] = m;
                    gridPointsIndexses[numFrontPoints][1] = -1;
                    numFrontPoints++;
                }

                previousFrontPosition = currentFrontPosition;
            }

            // add last j-edge, check for circular connectivity
            int lastPoint = allGridPoints[0].size() - 2;
            GetNeighbours(allGridPoints[0], lastPoint, currentLeftIndex, currentRightIndex);
            if(currentRightIndex == allGridPoints[0].size() - 2)
            {
                for (int i = previousFrontPosition; i >= 0; --i)
                {
                    frontGridPoints[numFrontPoints] = allGridPoints[i][lastPoint];
                    gridPointsIndexses[numFrontPoints][0] = lastPoint;
                    gridPointsIndexses[numFrontPoints][1] = i;
                    numFrontPoints++;
                }

            }

            return true;
        }


        //comp_vel
        bool ComputeVelocities
        (
            const std::vector<Point>& gridPoints, 
            const std::vector<double>& edgeVelocities,  
            std::vector<Point>& velocityVector
        )
        {
            std::fill(velocityVector.begin(), velocityVector.end(), Point{ doubleMissingValue,doubleMissingValue });
            Point normalVectorLeft;
            Point normalVectorRight;
            const double cosTolerance = 1e-8;
            double eps = std::numeric_limits<double>::min();
            for (int m = 0; m < velocityVector.size(); ++m)
            {
                if (!gridPoints[m].IsValid())
                {
                    continue;
                }
                int currentLeftIndex;
                int currentRightIndex;
                GetNeighbours(gridPoints, m, currentLeftIndex, currentRightIndex);

                double leftRightDistance = Distance(gridPoints[currentLeftIndex], gridPoints[currentRightIndex], m_projection);
                double leftDistance = Distance(gridPoints[currentLeftIndex], gridPoints[m], m_projection);
                double rightDistance = Distance(gridPoints[currentRightIndex], gridPoints[m], m_projection);

                if (leftRightDistance <= m_onTopOfEachOtherTolerance)
                {
                    continue;
                }

                if (leftDistance <= m_onTopOfEachOtherTolerance || rightDistance <= m_onTopOfEachOtherTolerance)
                {
                    NormalVectorOutside(gridPoints[currentRightIndex], gridPoints[currentLeftIndex], normalVectorLeft, m_projection);
                    if (m_projection == Projections::spherical)
                    {
                        normalVectorLeft.x = normalVectorLeft.x * std::cos(degrad_hp * 0.5 *(gridPoints[currentLeftIndex].y + gridPoints[currentRightIndex].y));
                    }
                    normalVectorRight = normalVectorLeft;
                }
                else
                {
                    NormalVectorOutside(gridPoints[m], gridPoints[currentLeftIndex], normalVectorLeft, m_projection);
                    NormalVectorOutside(gridPoints[currentRightIndex], gridPoints[m], normalVectorRight, m_projection);

                    if (m_projection == Projections::spherical)
                    {
                        normalVectorLeft.x = normalVectorLeft.x * std::cos(degrad_hp * 0.5 *(gridPoints[currentLeftIndex].y + gridPoints[m].y));
                        normalVectorRight.x = normalVectorRight.x * std::cos(degrad_hp * 0.5 *(gridPoints[currentRightIndex].y + gridPoints[m].y));
                    }
                }

                if (currentLeftIndex == velocityVector.size() - 1)
                {
                    continue;
                }

                double cosphi = DotProduct(normalVectorLeft.x, normalVectorRight.x, normalVectorLeft.y, normalVectorRight.y);
                Point leftVelocity = normalVectorLeft * edgeVelocities[currentLeftIndex];
                Point rightVelocity = normalVectorRight * edgeVelocities[currentRightIndex - 1];
                double rightLeftVelocityRatio = edgeVelocities[currentRightIndex - 1] / edgeVelocities[currentLeftIndex];

                if (cosphi -(cosTolerance - 1.0) < eps)
                {
                    continue;
                }

                if (rightLeftVelocityRatio - cosphi > eps  && 1.0 / rightLeftVelocityRatio - cosphi > eps || cosphi<= cosTolerance)
                {
                    velocityVector[m] = (leftVelocity * (1.0 - rightLeftVelocityRatio * cosphi) +
                        rightVelocity * (1.0 - 1.0 / rightLeftVelocityRatio*cosphi)) / (1.0 - cosphi*cosphi);
                }
                else if (rightLeftVelocityRatio - cosphi < eps)
                {
                    velocityVector[m] = leftVelocity * rightLeftVelocityRatio / cosphi;
                }
                else 
                {
                    velocityVector[m] = rightVelocity * 1.0 / (rightLeftVelocityRatio * cosphi);
                }

                if (m_projection == Projections::spherical)
                {
                    // use to spherical?
                    velocityVector[m].TransformToSpherical();
                    // TODO: check the results
                    //vel(1, i) = vel(1, i) * Rai*rd2dg / cos(dg2rd*yc(i));
                    //vel(2, i) = vel(2, i) * Rai*rd2dg; 
                }
            }
            return true;
        }


        bool GetNeighbours(
            const std::vector<Point>& gridPoints, 
            const int index,
            int& currentLeftIndex,
            int& currentRightIndex)
        {
            bool circularConnection = false;
            currentLeftIndex = index; 
            currentRightIndex = index;
            int start = 0;
            int end = gridPoints.size() - 1;

            // left
            while (Distance(gridPoints[currentLeftIndex], gridPoints[index], m_projection) < m_onTopOfEachOtherTolerance) 
            {
                if (!circularConnection)
                {
                    if (currentLeftIndex - 1 < 0) 
                    {
                        break;
                    }
                }
                else if (currentLeftIndex - 1 < 0)
                {
                    currentLeftIndex = end + 1;
                    circularConnection = false;
                }
                if(currentLeftIndex - 1 < 0 || !gridPoints[currentLeftIndex -1 ].IsValid())
                {
                    break;
                }
                currentLeftIndex--;
            }

            // right
            while (Distance(gridPoints[currentRightIndex], gridPoints[index], m_projection) < m_onTopOfEachOtherTolerance)
            {
                if (!circularConnection)
                {
                    if (currentRightIndex + 1 > gridPoints.size())
                    {
                        break;
                    }
                }
                else if (currentRightIndex + 1 > gridPoints.size())
                {
                    currentRightIndex = start - 1;
                    circularConnection = false;
                }
                if (currentRightIndex + 1>= gridPoints.size() || !gridPoints[currentRightIndex + 1].IsValid())
                {
                    break;
                }
                currentRightIndex++;
            }
            return true;
        }
    
        ///comp_edgevel
        bool ComputeEdgeVelocities(
            std::vector<double>& edgeVelocities, // edgevel
            std::vector<std::vector<double>>& growFactorOnSubintervalAndEdge, //dgrow1
            std::vector<std::vector<int>>& numPerpendicularFacesOnSubintervalAndEdge //nfac1
        )
        {

            bool success = ComputeGridHeights();
            if (!success)
            {
                return false;
            }

            for (int s = 0; s < m_numSplines; s++)
            {
                double maxHeight = std::numeric_limits<double>::min();
                for (int i = 0; i < m_gridHeights[0].size(); ++i)
                {
                    if (m_gridHeights[0][i] != doubleMissingValue && m_gridHeights[0][i] > maxHeight)
                    {
                        maxHeight = m_gridHeights[0][i];
                    }
                }

                double firstHeight = std::min(maxHeight, m_aspectRatio * m_averageMeshWidth);

                // Get true crossing splines heights
                int numLeftHeights = m_maxNumCenterSplineHeights;
                int numRightHeights = m_maxNumCenterSplineHeights;
                int numTrueCrossings = 0;
                for (int i = 0; i < m_numCrossingSplines[s]; ++i)
                {
                    if (m_numLayers[m_crossingSplinesIndexses[s][i]] != 1)
                    {
                        // true crossing splines only
                        continue;
                    }
                    numTrueCrossings++;
                    numLeftHeights = std::min(numLeftHeights, m_numCrossSplineLeftHeights[s][i]);
                    numRightHeights = std::min(numRightHeights, m_numCrossSplineRightHeights[s][i]);
                }

                // no true cross splines: exponentially growing grid only
                if (numTrueCrossings == 0)
                {
                    numLeftHeights = 0;
                    numRightHeights = 0;
                }


                int startGridLineLeft = m_leftGridLineIndex[s];
                int endGridLineLeft = startGridLineLeft + m_numMSpline[s];
                int startGridLineRight = m_rightGridLineIndex[s];
                int endGridLineRight = startGridLineRight + m_numMSpline[s];

                double hh0LeftMaxRatio;
                double hh0RightMaxRatio;
                const int numIterations = 2;
                for (int iter = 0; iter < numIterations; ++iter)
                {
                    ComputeVelocitiesSubIntervals(s, startGridLineLeft, endGridLineLeft, numLeftHeights, numRightHeights, firstHeight,
                        m_leftGridLineIndex, m_rightGridLineIndex, numPerpendicularFacesOnSubintervalAndEdge, edgeVelocities, hh0LeftMaxRatio);

                    ComputeVelocitiesSubIntervals(s, startGridLineRight, endGridLineRight, numLeftHeights, numRightHeights, firstHeight,
                        m_rightGridLineIndex, m_leftGridLineIndex, numPerpendicularFacesOnSubintervalAndEdge, edgeVelocities, hh0RightMaxRatio);
                }

                // re-evaluate if growing grid outside is needed

                if (numLeftHeights == 0 && numRightHeights <= 1 ||
                    numRightHeights == 0 && numLeftHeights <= 1 ||
                    numLeftHeights == numRightHeights == 1)
                {
                    m_growGridOutside = true;
                }

                // left part
                int numNLeftExponential = 0;
                if (m_growGridOutside)
                {
                    numNLeftExponential = std::min(ComputeNumberExponentialIntervals(hh0LeftMaxRatio), m_maxNumN);
                }
                for (int i = startGridLineLeft; i < endGridLineLeft; ++i)
                {
                    numPerpendicularFacesOnSubintervalAndEdge[1][i] = numNLeftExponential;
                }

                // right part
                int numNRightExponential = 0;
                if (m_growGridOutside)
                {
                    numNRightExponential = std::min(ComputeNumberExponentialIntervals(hh0RightMaxRatio), m_maxNumN);
                }
                for (int i = startGridLineRight; i < endGridLineRight; ++i)
                {
                    numPerpendicularFacesOnSubintervalAndEdge[1][i] = numNRightExponential;
                }
            }

            // compute local grow factors
            for (int s = 0; s < m_numSplines; s++)
            {
                if (m_numMSpline[s] < 1)
                {
                    continue;
                }

                for (int i = m_leftGridLineIndex[s]; i < m_rightGridLineIndex[s] + m_numMSpline[s]; ++i)
                {
                    if (m_gridLine[i].x == doubleMissingValue || m_gridLine[i + 1].x == doubleMissingValue || numPerpendicularFacesOnSubintervalAndEdge[1][i] < 1)
                    {
                        continue;
                    }
                    bool successfull = ComputeGrowFactor(m_gridHeights[1][i],
                        edgeVelocities[i],
                        numPerpendicularFacesOnSubintervalAndEdge[1][i],
                        growFactorOnSubintervalAndEdge[1][i]);
                    if (!successfull)
                    {
                        growFactorOnSubintervalAndEdge[1][i] = 1.0;
                    }
                }
            }

            return true;
        }

        ///comp_dgrow: this is another root finding algorithm, could go in the general part
        bool ComputeGrowFactor(
            const double totalGridHeight,
            const int numGridLayers,
            const double firstGridLayerHeight,
            double& result)
        {
            // eheight m_gridHeights
            double aspectRatioGrowFactor = 1.0;
            double heightDifference = ComputeTotalExponentialHeight(aspectRatioGrowFactor, numGridLayers, firstGridLayerHeight) - totalGridHeight;

            double deps = 0.01;
            double aspectRatioGrowFactorIncremented = 1.0 + deps;
            double heightDifferenceIncremented = ComputeTotalExponentialHeight(aspectRatioGrowFactorIncremented, numGridLayers, firstGridLayerHeight) - totalGridHeight;

            const double tolerance = 1e-8;
            const int numIterations = 1000;
            const double relaxationFactor = 0.5;
            double oldAspectRatio;
            double oldHeightDifference;

            if (std::abs(heightDifferenceIncremented) > tolerance && std::abs(heightDifferenceIncremented - heightDifference) > tolerance)
            {
                for (int i = 0; i < numIterations; ++i)
                {
                    oldAspectRatio = aspectRatioGrowFactor;
                    oldHeightDifference = heightDifference;

                    aspectRatioGrowFactor = aspectRatioGrowFactorIncremented;
                    heightDifference = heightDifferenceIncremented;

                    aspectRatioGrowFactorIncremented = aspectRatioGrowFactor - relaxationFactor * heightDifference / (heightDifference - oldHeightDifference + 1e-16) * (aspectRatioGrowFactor - oldAspectRatio);
                    heightDifferenceIncremented = ComputeTotalExponentialHeight(aspectRatioGrowFactorIncremented, numGridLayers, firstGridLayerHeight) - totalGridHeight;

                    if (std::abs(oldHeightDifference) < tolerance)
                    {
                        break;
                    }
                }
            }

            if (oldHeightDifference > tolerance)
            {
                result = doubleMissingValue;
                return false;
            }

            result = aspectRatioGrowFactorIncremented;
            return true;
        }

        double ComputeTotalExponentialHeight(const double aspectRatioGrowFactor, const int numberOfGridLayers, const double firstGridLayerHeight)
        {
            double height;
            if (aspectRatioGrowFactor - 1.0 > 1e-8)
            {
                height = (std::pow(aspectRatioGrowFactor, firstGridLayerHeight) - 1.0) / (aspectRatioGrowFactor - 1.0) * numberOfGridLayers;
            }
            else
            {
                height = numberOfGridLayers * firstGridLayerHeight;
            }
            return height;
        }

        ///comp_nfac
        ///compute the number of grid layers for a given grow factor, first grid layer height and total grid height
        int ComputeNumberExponentialIntervals(const double hhMaxRatio)
        {
            int numIntervals = 0;
            if( m_aspectRatioGrowFactor - 1.0 > 1e-8)
            {
                numIntervals = std::floor(std::log((m_aspectRatioGrowFactor - 1.0) * hhMaxRatio + 1.0) / log(m_aspectRatioGrowFactor));
            }
            else
            {
                numIntervals = std::floor(0.999 + hhMaxRatio);
            }
            return numIntervals;
        }



        bool ComputeVelocitiesSubIntervals(
            const int s,
            const int startGridLineIndex,
            const int endGridLineIndex,
            const int numHeights,
            const int numOtherSideHeights,
            const double firstHeight,
            const std::vector<int>& gridLineIndex,
            const std::vector<int>& otherGridLineIndex,
            std::vector<std::vector<int>>& numPerpendicularFacesOnSubintervalAndEdge,
            std::vector<double>& edgeVelocities,
            double& hh0MaxRatio)
        {
            hh0MaxRatio = 0.0;
            if (numHeights > 1 && numHeights == numOtherSideHeights || numHeights > numOtherSideHeights)
            {
                double maxHeight = *std::max_element(m_gridHeights[0].begin() + startGridLineIndex, m_gridHeights[0].begin() + endGridLineIndex);

                int numNUniformPart = std::floor(maxHeight / firstHeight + 0.99999);
                numNUniformPart = std::min(numNUniformPart, m_maxNUniformPart);

                for (int i = startGridLineIndex; i < endGridLineIndex; ++i)
                {
                    numPerpendicularFacesOnSubintervalAndEdge[0][i] = numNUniformPart;
                    edgeVelocities[i] = m_gridHeights[0][i] / numNUniformPart;
                    hh0MaxRatio = std::max(hh0MaxRatio, m_gridHeights[1][i] / edgeVelocities[i]);
                }
            }
            else
            {
                // only one subinterval: no uniform part
                int numNUniformPart = 0;
                for (int i = startGridLineIndex; i < endGridLineIndex; ++i)
                {
                    numPerpendicularFacesOnSubintervalAndEdge[0][i] = numNUniformPart;
                    edgeVelocities[i] = firstHeight;

                    //compare with other side of spline
                    int otherSideIndex = otherGridLineIndex[s] + m_numMSpline[s] - (i - gridLineIndex[s] + 1);

                    if (edgeVelocities[otherSideIndex] != doubleMissingValue)
                    {
                        if (numPerpendicularFacesOnSubintervalAndEdge[0][otherSideIndex] == 0)
                        {
                            edgeVelocities[i] = std::max(edgeVelocities[i], edgeVelocities[otherSideIndex]);
                        }
                        else
                        {
                            edgeVelocities[i] = edgeVelocities[otherSideIndex];
                        }
                    }

                    for (int j = 1; j < m_maxNumCenterSplineHeights; ++j)
                    {
                        m_gridHeights[j][i] = m_gridHeights[j - 1][i];
                    }

                    for (int j = startGridLineIndex; j < endGridLineIndex; ++j)
                    {

                        hh0MaxRatio = std::max(hh0MaxRatio, m_gridHeights[1][j] / edgeVelocities[j]);
                    }
                }
            }

            return true;
        }

        /// compute the grid heights at grid edges on the center spline
        ///comp_gridheights
        bool ComputeGridHeights()
        {
            m_gridHeights.resize(m_maxNumCenterSplineHeights, std::vector<double>(m_numM - 1 ));
            std::fill(m_gridHeights.begin(), m_gridHeights.end(), std::vector<double>(m_numM - 1, doubleMissingValue));

            std::vector<std::vector<double>> heightsLeft(m_maxNumCenterSplineHeights, std::vector<double>(m_maxNumM, 0.0));
            std::vector<std::vector<double>> heightsRight(m_maxNumCenterSplineHeights, std::vector<double>(m_maxNumM, 0.0));
            int maxNumCornerPoints = *std::max_element(m_numSplineNodes.begin(), m_numSplineNodes.end());
            std::vector<double> edgesCenterPoints(maxNumCornerPoints, 0.0);
            std::vector<double> crossingSplinesDimensionalCoordinates(m_numSplines, 0.0);
            std::vector<int> numHeightsLeft(m_numSplines, 0.0);
            std::vector<int> numHeightsRight(m_numSplines, 0.0);
            std::vector<double> localSplineDerivatives(m_numSplines, 0.0);
            std::vector<int> localValidSplineIndexes(m_numSplines, 0.0);

            for (int s = 0; s < m_numSplines; s++)
            {
                if (m_numLayers[s] != 0)
                {
                    continue;
                }

                int numM = m_numMSpline[s];

                // Get the minimum number of sub-intervals in the cross splines for this center spline
                int minNumLeftIntervals = *std::min_element(m_numCrossSplineLeftHeights[s].begin(), m_numCrossSplineLeftHeights[s].end());
                int minNumRightIntervals = *std::min_element(m_numCrossSplineRightHeights[s].begin(), m_numCrossSplineRightHeights[s].end());

                std::fill(heightsLeft[0].begin(), heightsLeft[0].begin() + numM, m_maximumGridHeights[s]);
                std::fill(heightsRight[0].begin(), heightsRight[0].begin() + numM, m_maximumGridHeights[s]);

                if (m_numCrossingSplines[s] == 1)
                {
                    // only one crossing spline present: 
                    for (int i = 0; i < minNumLeftIntervals; ++i)
                    {
                        std::fill(heightsLeft[i].begin(), heightsLeft[i].begin() + numM, m_crossSplineRightHeights[s][i][0]);

                    }
                    for (int i = 0; i < minNumRightIntervals; ++i)
                    {
                        std::fill(heightsRight[i].begin(), heightsRight[i].begin() + numM, m_crossSplineLeftHeights[s][i][0]);
                    }
                }
                else
                {
                    int leftGridLineIndex = m_leftGridLineIndex[s];
                    edgesCenterPoints[0] = GetSplineLength(s, 0, m_gridLineDimensionalCoordinates[leftGridLineIndex]);
                    for (int i = 0; i < numM; ++i)
                    {
                        edgesCenterPoints[i + 1] = edgesCenterPoints[i] + GetSplineLength(s, m_gridLineDimensionalCoordinates[leftGridLineIndex + i], m_gridLineDimensionalCoordinates[leftGridLineIndex + i + 1]);
                    }

                    // compute at edge center points
                    for (int i = 0; i < numM; ++i)
                    {
                        edgesCenterPoints[i] = 0.5 *(edgesCenterPoints[i] + edgesCenterPoints[i + 1]);
                    }
                    edgesCenterPoints[numM] = doubleMissingValue;

                    //compute center spline path length of cross splines
                    crossingSplinesDimensionalCoordinates[0] = GetSplineLength(s, 0.0, m_crossSplineCoordinates[s][0]);
                    for (int i = 0; i < m_numCrossingSplines[s] - 1; ++i)
                    {
                        crossingSplinesDimensionalCoordinates[i + 1] = crossingSplinesDimensionalCoordinates[i] + GetSplineLength(s, m_crossSplineCoordinates[s][i], m_crossSplineCoordinates[s][i + 1]);
                        numHeightsLeft[i] = m_numCrossSplineLeftHeights[s][i];
                        numHeightsRight[i] = m_numCrossSplineRightHeights[s][i];
                    }


                    for (int j = 0; j < m_maxNumCenterSplineHeights; ++j)
                    {
                        AddValueToVector(numHeightsLeft, -1);
                        AddValueToVector(numHeightsRight, -1);

                        FindNearestCrossSplines(s, j,
                            numHeightsLeft,
                            edgesCenterPoints,
                            m_crossSplineLeftHeights[s],
                            localValidSplineIndexes,
                            localSplineDerivatives,
                            crossingSplinesDimensionalCoordinates,
                            heightsLeft);

                        FindNearestCrossSplines(s, j,
                            numHeightsRight,
                            edgesCenterPoints,
                            m_crossSplineRightHeights[s],
                            localValidSplineIndexes,
                            localSplineDerivatives,
                            crossingSplinesDimensionalCoordinates,
                            heightsRight);
                    }
                }

                
                // store grid height
                for (int j = 0; j < m_maxNumCenterSplineHeights; ++j)
                {
                    for (int i = 0; i < m_numMSpline[s]; ++i)
                    {
                        m_gridHeights[j][m_leftGridLineIndex[s] + i] = heightsLeft[j][i];
                        m_gridHeights[j][m_rightGridLineIndex[s] + m_numMSpline[s] - i - 1] = heightsRight[j][i];
                    }
                }
            }

            return true;
        }

        bool FindNearestCrossSplines(const int s, 
            const int j,
            const std::vector<int>& numHeightsLeft,
            const std::vector<double>& edgesCenterPoints,
            const std::vector<std::vector<double>>& crossSplineLeftHeights,
            std::vector<int>& localValidSplineIndexes,
            std::vector<double>& localSplineDerivatives,
            std::vector<double>& crossingSplinesDimensionalCoordinates,
            std::vector<std::vector<double>>& heights)
        {
            int numValid;
            GetValidSplineIndexses(s, m_numCrossingSplines[s], numHeightsLeft, localValidSplineIndexes, numValid);

            // no sub-heights to compute 
            if (numValid == 0)
            {
                return false;
            }

            int numM = m_numMSpline[s];
            std::vector<double> localCornerPoints(numValid);
            
            // TODO: strided memory access
            for (int i = 0; i < numValid; ++i)
            {
                int index = localValidSplineIndexes[i];
                localCornerPoints[i] = crossSplineLeftHeights[index][j];
            }
            
            SecondOrderDerivative(localCornerPoints, numValid, localSplineDerivatives);

            crossingSplinesDimensionalCoordinates[0] = GetSplineLength(s, 0.0, m_crossSplineCoordinates[s][0]);
            for (int i = 0; i < numM; ++i)
            {
                double leftCoordinate = crossingSplinesDimensionalCoordinates[localValidSplineIndexes[0]];
                int rightIndex = std::min(1, numValid - 1);
                double rightCoordinate = crossingSplinesDimensionalCoordinates[localValidSplineIndexes[rightIndex]];
                int leftIndex = 1;
                // Find two nearest cross splines
                for (int k = leftIndex; k < numValid; k++)
                {
                    leftIndex = k;
                    leftCoordinate = crossingSplinesDimensionalCoordinates[localValidSplineIndexes[k]];
                    rightIndex++;
                    rightCoordinate = crossingSplinesDimensionalCoordinates[localValidSplineIndexes[rightIndex]];
                    if (rightIndex >= numValid || rightCoordinate >= edgesCenterPoints[i])
                    {
                        break;
                    }
                }

                double factor = 0.0;
                if (rightCoordinate - leftCoordinate > 1e-8)
                {
                    factor = (crossingSplinesDimensionalCoordinates[i] - leftCoordinate) / (rightCoordinate - leftCoordinate);
                }
                else
                {
                    rightIndex = leftIndex;
                }

                factor = std::max(std::min(double(leftIndex) + factor - 1.0, double(numValid - 1)), 0.0);

                Interpolate(localCornerPoints, localSplineDerivatives, factor, heights[j][i]);
            }

            return true;
        }


        // GetValidSplineIndexses
        bool GetValidSplineIndexses(const int s, const int numValues, const std::vector<int>& v, std::vector<int>& validIndexses, int& numValid)
        {
            numValid = 0;
            for (int i = 0; i < numValues; ++i)
            {
                if (v[i] >= 0)
                {
                    validIndexses[numValid] = i;
                    numValid++;
                }
            }
            return true;
        };

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

        bool GetSplineIntersections(const int index)
        {
            m_numCrossingSplines[index] = 0;
            std::fill(m_crossingSplinesIndexses[index].begin(), m_crossingSplinesIndexses[index].end(), -1);
            std::fill(m_isLeftOriented[index].begin(), m_isLeftOriented[index].end(), false);
            std::fill(m_crossSplineCoordinates[index].begin(), m_crossSplineCoordinates[index].end(), std::numeric_limits<double>::max());
            std::fill(m_cosCrossingAngle[index].begin(), m_cosCrossingAngle[index].end(), -1);

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
                bool crossing = GetSplinesIntersection(index, s, m_projection, crossProductIntersection, intersectionPoint, firstSplineRatio, secondSplineRatio);

                if(std::abs(crossProductIntersection)<m_dtolcos)
                {
                    crossing = false;
                }


                if (crossing)
                {
                    m_numCrossingSplines[index]++;
                    m_crossingSplinesIndexses[index][s] = s;
                    m_isLeftOriented[index][s] = true;
                    if (crossProductIntersection > 0.0)
                    {
                        m_isLeftOriented[index][s] = false;
                    }
                    m_crossSplineCoordinates[index][s] = firstSplineRatio;
                    m_cosCrossingAngle[index][s] = crossProductIntersection;
                }
            }

            const auto sortedIndexses = SortedIndexes(m_crossSplineCoordinates[index]);
            ReorderVector(m_crossSplineCoordinates[index], sortedIndexses);
            ReorderVector(m_crossingSplinesIndexses[index], sortedIndexses);
            ReorderVector(m_isLeftOriented[index], sortedIndexses);
            ReorderVector(m_cosCrossingAngle[index], sortedIndexses);

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

            // allocate
            m_gridLine.resize(numCenterSplines*(m_maxNumM + 2));
            m_gridLineDimensionalCoordinates.resize(numCenterSplines*(m_maxNumM + 2));

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
                    gridLineIndex++;
                    m_gridLine[gridLineIndex] = {doubleMissingValue,doubleMissingValue};
                    m_gridLineDimensionalCoordinates[gridLineIndex] = doubleMissingValue;
                }

                m_leftGridLineIndex[s] = gridLineIndex;

                int numM = 0;
                bool success = MakeGridLine(s, isSpacingCurvatureAdapeted, gridLineIndex, m_gridLine, m_gridLineDimensionalCoordinates, numM);
                if (!success)
                {
                    return false;
                }

                gridLineIndex = gridLineIndex + numM + 1;
                m_gridLine[gridLineIndex] = Point{doubleMissingValue,  doubleMissingValue };
                m_gridLineDimensionalCoordinates[gridLineIndex] = doubleMissingValue;
                gridLineIndex++;

                //add other side of gridline
                m_rightGridLineIndex[s] = gridLineIndex;
                for (int i = m_rightGridLineIndex[s] - 1, j = m_rightGridLineIndex[s] - 1; j >= m_leftGridLineIndex[s]; ++i, --j)
                {
                    m_gridLine[i] = m_gridLine[j];
                    m_gridLineDimensionalCoordinates[i] = m_gridLineDimensionalCoordinates[j];
                }

                //compute new (actual) grid size
                //new size   old size   both sides of spline   DMISS between both sides
                gridLineIndex = gridLineIndex + numM + 1;

                m_gridLine[gridLineIndex] = Point{ doubleMissingValue,  doubleMissingValue };
                m_gridLineDimensionalCoordinates[gridLineIndex] = doubleMissingValue;
                m_numMSpline[s] = numM;
                m_numM = gridLineIndex;
            }

            return true;
        }

        /// make_gridline, generate a gridline on a spline with a prescribed maximum mesh width
        /// generate a gridline on a spline with a prescribed maximum mesh width
        /// call makespl(startstop, xsp, ysp, max(mc,num), num, 2, mc-1, xc, yc, kmax, sc, h)
        ///     MAKESPL(        T,   X,   Y,        imax,   N,NT,MNFAC,  XH,YH, KMAX, TT, H)
        bool MakeGridLine(const int splineIndex,
            const bool isSpacingCurvatureAdapeted,
            const int startingIndex,
            std::vector<Point>& gridLine,
            std::vector<double>& adimensionalCoordinates,
            int& numM)
        {
            // first estimation of nodes along m
            numM = 1 + std::floor(m_splinesLength[splineIndex] / m_averageMeshWidth);
            numM = std::min(numM, m_maxNumM);

            double endSplineAdimensionalCoordinate = m_numSplineNodes[splineIndex] - 1;
            double splineLength = GetSplineLength(splineIndex, 0.0, endSplineAdimensionalCoordinate, 10, isSpacingCurvatureAdapeted, m_maximumGridHeights[splineIndex]);
            
            gridLine[startingIndex] = m_splineCornerPoints[splineIndex][0];
            FuncDimensionalToAdimensionalDistance func(*this, splineIndex, isSpacingCurvatureAdapeted, m_maximumGridHeights[splineIndex]);
            
            double currentMaxWidth = std::numeric_limits<double>::max();
            while (currentMaxWidth > m_averageMeshWidth && numM < m_maxNumM)
            {
                currentMaxWidth = 0.0;
                for (int n = 1; n < numM + 1; ++n)
                {
                    int index = startingIndex + n;
                    func.SetDimensionalDistance(splineLength * double(n) / double(numM));
                    adimensionalCoordinates[index] = FindFunctionRootWithGoldenSectionSearch(func, 0, endSplineAdimensionalCoordinate);
                    Interpolate(m_splineCornerPoints[splineIndex], m_splineDerivatives[splineIndex], adimensionalCoordinates[index], gridLine[index]);
                    currentMaxWidth = std::max(currentMaxWidth, Distance(gridLine[index - 1], gridLine[index], m_projection));
                }

                // room for sub-division
                if (currentMaxWidth > m_averageMeshWidth)
                {
                    numM = std::min(std::max(int(m_maxNumM / m_maximumGridHeights[splineIndex] *numM), numM + 1), m_maxNumM);
                }
            }
            return true;
        }

        /// functor wrapping GetSplineLength to use generic root finders
        /// https://www.fluentcpp.com/2017/01/23/stl-function-objects-stateless-is-stressless/
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
        bool ComputeSplineProperties(const bool restoreOriginalProperties)
        {
            // cross spline is a spline with only two points, others are non-cross splines
            bool successful = AllocateSplinesProperties();
            if (!successful)
            {
                return false;
            }

            for (int s = 0; s < m_numSplines; ++s)
            {
                successful = GetSplineIntersections(s);

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
                if (m_numSplineNodes[s] != 2 || m_numCrossingSplines[s] < 1)
                {
                    continue;
                }

                int middleCrossingSpline = std::min(m_numCrossingSplines[s] / 2, m_numCrossingSplines[s]);
                int crossingSplineIndex = m_crossingSplinesIndexses[s][middleCrossingSpline];

                // if m_numIntersectingSplines[s] is even, check if the middle spline has already been assigned as a bounding spline
                if (m_numLayers[crossingSplineIndex] != 0 && 2 * crossingSplineIndex == m_numCrossingSplines[s])
                {
                    middleCrossingSpline = std::min(middleCrossingSpline + 1, m_numCrossingSplines[s] - 1);
                    crossingSplineIndex = m_crossingSplinesIndexses[s][middleCrossingSpline];
                }
                
                if (m_numLayers[crossingSplineIndex] == 0)
                {
                    // associate bounding splines with the middle spline
                    for (int i = 0; i < middleCrossingSpline; ++i)
                    {
                        int index = m_crossingSplinesIndexses[s][i];
                        m_numLayers[index] = -crossingSplineIndex;

                    }
                    for (int i = middleCrossingSpline + 1; i < m_numCrossingSplines[s]; ++i)
                    {
                        int index = m_crossingSplinesIndexses[s][i];
                        m_numLayers[index] = -crossingSplineIndex;
                    }
                }
            }

            if(restoreOriginalProperties)
            {
                // restore original spline properties
                for (int s = 0; s < m_numOriginalSplines; ++s)
                {
                    m_leftGridLineIndex[s] = m_leftGridLineIndexOriginal[s];
                    m_rightGridLineIndex[s] = m_rightGridLineIndexOriginal[s];
                    m_numMSpline[s] = m_mfacOriginal[s];
                    m_maximumGridHeights[s] = m_maximumGridHeightsOriginal[s];
                    m_numLayers[s] = m_numLayersOriginal[s];
                }

                //mark new splines as artificial cross splines
                for (int s = m_numOriginalSplines; s < m_numSplines; ++s)
                {
                    m_numLayers[s] = 3;
                }
            }

            bool successfull = ComputeHeights();
            return successfull;
        }

        /// get the grid heights from the cross spline information
        /// get_heights
        bool ComputeHeights()
        {
            for (int i = 0; i < m_numSplines; ++i)
            {
                // Heights should be computed only for center splines
                if (m_numSplineNodes[i] <= 2)
                {
                    continue;
                }
                for (int j = 0; j < m_numCrossingSplines[i]; ++j)
                {
                    int intersectingSplineIndex = m_crossingSplinesIndexses[i][j];
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
                if (m_numCrossingSplines[s] == 0)
                {
                    m_maximumGridHeights[s] = m_aspectRatioFirstLayer * m_splinesLength[s];
                    continue;
                }
                double maximumHeight = 0.0;
                for (int c = 0; c <  m_numCrossingSplines[s]; ++c)
                {
                    double sumLeftHeights = 0.0;
                    for (int ss = 0; ss < m_numCrossSplineLeftHeights[s][c]; ++ss)
                    {
                        sumLeftHeights += m_crossSplineLeftHeights[s][c][ss];
                    }
                    double sumRightHeights = 0.0;
                    for (int ss = 0; ss < m_numCrossSplineRightHeights[s][c]; ++ss)
                    {
                        sumRightHeights += m_crossSplineRightHeights[s][c][ss];
                    }
                    maximumHeight = std::max(maximumHeight, std::max(sumLeftHeights, sumRightHeights));
                }

                m_maximumGridHeights[s] = maximumHeight;
            }
            return true;
        }

        ///comp_subheights, compute the height of the subintervals of grid layers on a cross spline, w.r.t. a center spline
        bool ComputeSubHeights(const int centerSplineIndex, const int crossingSplineLocalIndex)
        {
            // find center spline index
            int centerSplineLocalIndex = 0;
            int crossingSplineIndex = m_crossingSplinesIndexses[centerSplineIndex][crossingSplineLocalIndex]; //js
            //m_intersectingSplinesIndexses[intersectingSplinesIndex] // ics
            for (int s = 0; s < m_numCrossingSplines[crossingSplineIndex]; ++s)
            {
                if (m_crossingSplinesIndexses[crossingSplineIndex][s]== centerSplineIndex)
                {
                    centerSplineLocalIndex = s;
                    break;
                }
            }

            // right part
            int numSubIntervalsRight = 0;
            int rightCenterSplineIndex = centerSplineLocalIndex;
            int leftCenterSplineIndex;
            m_crossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex].resize(m_maxNumCenterSplineHeights, 0);
            for (int s = centerSplineLocalIndex; s <  m_numCrossingSplines[crossingSplineIndex] - 1; ++s)
            {
                if (numSubIntervalsRight >= m_maxNumCenterSplineHeights)
                {
                    break;
                }
                if (m_numLayers[m_crossingSplinesIndexses[crossingSplineIndex][s + 1]] != -centerSplineIndex)
                {
                    continue;
                }
                leftCenterSplineIndex = rightCenterSplineIndex;
                rightCenterSplineIndex = s + 1;
                m_crossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex][numSubIntervalsRight] = GetSplineLength(crossingSplineIndex,
                    m_crossSplineCoordinates[crossingSplineIndex][leftCenterSplineIndex], m_crossSplineCoordinates[crossingSplineIndex][rightCenterSplineIndex]);
                numSubIntervalsRight++;
            }

            m_crossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex][numSubIntervalsRight] = GetSplineLength(crossingSplineIndex,
                m_crossSplineCoordinates[crossingSplineIndex][rightCenterSplineIndex], m_numSplineNodes[crossingSplineIndex] - 1);

            numSubIntervalsRight++;
            std::fill(m_crossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex].begin() + numSubIntervalsRight, m_crossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex].end(), 0.0);
            m_numCrossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex] = numSubIntervalsRight;

            // left part
            int numSubIntervalsLeft = 0;
            leftCenterSplineIndex = centerSplineLocalIndex;
            m_crossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex].resize(m_maxNumCenterSplineHeights, 0);
            for (int s = centerSplineLocalIndex; s >= 1; --s)
            {
                if (numSubIntervalsLeft >= m_maxNumCenterSplineHeights)
                {
                    break;
                }
                if (m_numLayers[m_crossingSplinesIndexses[crossingSplineIndex][s - 1]] != -centerSplineIndex)
                {
                    continue;
                }
                rightCenterSplineIndex = leftCenterSplineIndex;
                leftCenterSplineIndex = s - 1;
                m_crossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex][numSubIntervalsLeft] = GetSplineLength(crossingSplineIndex,
                    m_crossSplineCoordinates[crossingSplineIndex][leftCenterSplineIndex], m_crossSplineCoordinates[crossingSplineIndex][rightCenterSplineIndex]);
                numSubIntervalsLeft++;
            }

            m_crossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex][numSubIntervalsLeft] = GetSplineLength(crossingSplineIndex,
                0.0, m_crossSplineCoordinates[crossingSplineIndex][leftCenterSplineIndex]);

            numSubIntervalsLeft++;
            std::fill(m_crossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex].begin() + numSubIntervalsLeft, m_crossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex].end(), 0.0);
            m_numCrossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex] = numSubIntervalsLeft;

            // if not left oriented, swap
            if (!m_isLeftOriented[centerSplineIndex][centerSplineLocalIndex])
            {
                m_numCrossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex] = numSubIntervalsRight;
                m_numCrossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex] = numSubIntervalsLeft;

                std::vector<double> leftSubIntervalsTemp(m_crossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex]);
                m_crossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex] = m_crossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex];
                m_crossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex] = leftSubIntervalsTemp;
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
            double dx = GetDx(pointCoordinate, incremenetedPointCoordinate, m_projection);
            double dy = GetDy(pointCoordinate, incremenetedPointCoordinate, m_projection);

            tangentialVector.x = dx / distance;
            tangentialVector.y = dy / distance;

            return true;
        }


        /// these functions are used in the gridgeom api, made them static
        /// SPLINE
        /// second order derivative of spline coordinates
        static inline bool SecondOrderDerivative(const std::vector<Point>& coordinates, int numNodes, std::vector<Point>& coordinatesDerivatives)
        {
            std::vector<Point> u(numNodes);
            u[0] = { 0.0, 0.0 };
            coordinatesDerivatives.resize(coordinates.size(), { 0.0, 0.0 });
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

        static inline bool SecondOrderDerivative(const std::vector<double>& coordinates, int numNodes, std::vector<double>& coordinatesDerivatives)
        {
            std::vector<double> u(numNodes);
            u[0] = 0.0;
            //coordinatesDerivatives.resize(coordinates.size());
            coordinatesDerivatives[0] = 0.0 ;

            for (int i = 1; i < numNodes - 1; i++)
            {
                const double p = coordinatesDerivatives[i - 1] * 0.5 + 2.0;
                coordinatesDerivatives[i] = -0.5 / p;

                const double delta = coordinates[i + 1] - coordinates[i] - (coordinates[i] - coordinates[i - 1]);
                u[i] = (delta *6.0 / 2.0 - u[i - 1] * 0.5) / p;
            }

            coordinatesDerivatives[numNodes - 1] =  0.0;
            for (int i = numNodes - 2; i >= 0; i--)
            {
                coordinatesDerivatives[i] = coordinatesDerivatives[i] * coordinatesDerivatives[i + 1] + u[i];
            }

            return true;
        }

        /// splint
        template<typename T>
        static bool Interpolate(const std::vector<T>& coordinates, const std::vector<T>& coordinatesDerivatives, double pointAdimensionalCoordinate, T& pointCoordinate)
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

        int m_numSplines = 0;
        int m_numAllocatedSplines = 0;
        std::vector<int> m_numSplineNodes;
        std::vector<int> m_numAllocatedSplineNodes;
        std::vector<std::vector<Point>> m_splineCornerPoints;
        std::vector<std::vector<Point>> m_splineDerivatives;
        Projections m_projection;
        double m_aspectRatioFirstLayer = 0.10;
        int m_maxNumM = 20;                                                             // mfacmax
        int m_numM = 0;                                                                 // mc
        int m_maxNumN = 40;                                                             //  N - refinement factor for regular grid generation.
        int m_maxNumCenterSplineHeights = 10;                                           // Nsubmax, naz number of different heights a cross spline can have (is determined by how many crossing spline the user can input).
        double m_averageMeshWidth = 500.0;
        double m_dtolcos = 0.95;                                                        // minimum allowed absolute value of crossing - angle cosine
        double m_aspectRatio = 0.1;                                                     // daspect aspect ratio
        double m_aspectRatioGrowFactor = 1.1;                                           // dgrow grow factor of aspect ratio
        double m_gridLayerHeight0 = 10.0;                                               // dheight0 grid layer height
        double m_maxaspect = 1.0;                                                       // maxaspect maximum cell aspect ratio *inoperative*
        int m_maxNUniformPart = 5;                                                      // maximum number of layers in the uniform part
        bool m_growGridOutside = true;                                                  // grow the grid outside the prescribed grid height
        double m_onTopOfEachOtherTolerance = 1e-4;                                      // On - top - of - each - other tolerance *IMPORTANT*

        Polygons m_polygon;                                                             // selecting polygon

        // Spline properties (first index is the spline number)                 
        std::vector<int> m_numLayers;                                                   // id number of layers ( >0 only for center spline)
        std::vector<int> m_numCrossingSplines;                                          // ncs num of cross splines
        std::vector<std::vector<int>> m_crossingSplinesIndexses;                        // ics for each cross spline, the indexses of the center splines
        std::vector<double> m_splinesLength;                                            // splinesLength spline path length
        std::vector<double> m_maximumGridHeights;                                       // hmax maximum grid height
        std::vector<std::vector<bool>> m_isLeftOriented;                                // isLeftOriented cross spline is left to right(.true.) or not (.false.) w.r.t.center spline
        std::vector<std::vector<double>>  m_crossSplineCoordinates;                     // t center spline coordinates of cross splines
        std::vector<std::vector<double>> m_cosCrossingAngle;                            // cosPhi cosine of crossing angle
        std::vector<std::vector<std::vector<double>>> m_crossSplineLeftHeights;         // hL left - hand side grid heights at cross spline locations for each grid layer subinterval, hL(1, :) being the height of the first subinterval, etc.
        std::vector<std::vector<std::vector<double>>> m_crossSplineRightHeights;        // hR right - hand side grid heights at cross spline locations for each grid layer subinterval, hR(1, :) being the height of the first subinterval, etc.
        std::vector<std::vector<int>>   m_numCrossSplineLeftHeights;                    // NsubL number of subintervals of grid layers at cross spline locations at the left - hand side of the spline, each having their own exponential grow factor
        std::vector<std::vector<int>>  m_numCrossSplineRightHeights;                    // NsubR number of subintervals of grid layers at cross spline locations at the right - hand side of the spline, each having their own exponential grow factor
        std::vector<int> m_numMSpline;                                                  // mfac number of grid intervals on the spline
        std::vector<std::vector<int>> m_nfacL;                                          // nfacL number of grid layers in each subinterval at the left - hand side of the spline * not used yet*
        std::vector<std::vector<int>> m_nfacR;                                          // nfacR number of grid layers in each subinterval at the right - hand side of the spline * not used yet*
        std::vector<int> m_leftGridLineIndex;                                           // iL index in the whole gridline array of the first grid point on the left - hand side of the spline
        std::vector<int> m_rightGridLineIndex;                                          // iR index in the whole gridline array of the first grid point on the right - hand side of the spline
        std::vector<Point> m_gridLine;                                                  // xg1, yg1 coordinates of the first gridline
        std::vector<double> m_gridLineDimensionalCoordinates;                           // sg1 center spline coordinates of the first gridline
        

        int m_allocationSize = 5;                                                       // allocation cache size

        //original spline chaches
        int m_numOriginalSplines = 0;
        
        std::vector<int> m_leftGridLineIndexOriginal;
        std::vector<int> m_rightGridLineIndexOriginal;
        std::vector<int> m_mfacOriginal;
        std::vector<double> m_maximumGridHeightsOriginal;
        std::vector<int> m_numLayersOriginal;
        std::vector<std::vector<double>> m_gridHeights;                                 //eheight
    };

}
