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
                std::fill(m_crossSplineLeftHeights[s].begin(), m_crossSplineLeftHeights[s].end(), std::vector<double>(m_maxNumHeights, doubleMissingValue));

                m_crossSplineRightHeights[s].resize(m_numSplines);
                std::fill(m_crossSplineRightHeights[s].begin(), m_crossSplineRightHeights[s].end(), std::vector<double>(m_maxNumHeights, doubleMissingValue));

                m_numCrossSplineLeftHeights[s].resize(m_numSplines);
                std::fill(m_numCrossSplineLeftHeights[s].begin(), m_numCrossSplineLeftHeights[s].end(), doubleMissingValue);

                m_numCrossSplineRightHeights[s].resize(m_numSplines,doubleMissingValue);
                std::fill(m_numCrossSplineRightHeights[s].begin(), m_numCrossSplineRightHeights[s].end(), doubleMissingValue);

                m_nfacL[s].resize(m_numSplines);
                std::fill(m_nfacL[s].begin(), m_nfacL[s].end(), 0);

                m_nfacR[s].resize(m_numSplines);
                std::fill(m_nfacR[s].begin(), m_nfacR[s].end(), 0);
            }

            m_mfac.resize(m_numSplines);
            std::fill(m_mfac.begin(), m_mfac.end(), 0);
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

            // get the properties of the center splines
            bool success = MakeAllGridLines(true);
            
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
                for (int i = m_leftGridLineIndex[s]; i < m_leftGridLineIndex[s] + m_mfac[s]; ++i)
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
                m_mfacOriginal[s] = m_mfac[s];
                m_maximumGridHeightsOriginal[s] = m_maximumGridHeights[s];
                m_numLayersOriginal[s] = m_numLayers[s];
            }

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
            //ComputeEdgeVelocities();

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
            if(!success)
            {
                return false;
            }

            for (int s = 0; s < m_numSplines; s++)
            {
               
                double maxHeight = std::numeric_limits<double>::min();
                for (int i = 0; i < m_maxNumHeights; ++i)
                {
                    if(m_gridHeights[0][i]!=doubleMissingValue && m_gridHeights[0][i]> maxHeight)
                    {
                        maxHeight = m_gridHeights[0][i];
                    }
                }

                double firstHeight = std::min(maxHeight, m_aspectRatio * m_averageMeshWidth);

                // Get true crossing splines heights
                int numLeftHeights = m_maxNumHeights;
                int numRightHeights = m_maxNumHeights;
                int numTrueCrossings = 0;
                for (int i = 0; i < m_numCrossingSplines[s]; ++i)
                {
                    if(m_numLayers[m_crossingSplinesIndexses[s][i]]!=1)
                    {
                        // true crossing splines only
                        continue;
                    }
                    numTrueCrossings++;
                    numLeftHeights = std::min(numLeftHeights, m_numCrossSplineLeftHeights[s][i]);
                    numRightHeights = std::min(numRightHeights, m_numCrossSplineRightHeights[s][i]);
                }

                // no true cross splines: exponentially growing grid only
                if(numTrueCrossings==0)
                {
                    numLeftHeights = 0;
                    numRightHeights = 0;
                }


                int startHeightsLeft = m_leftGridLineIndex[s];
                int endHeightsLeft = startHeightsLeft + m_mfac[s] - 1;
                double hh0LeftMaxRatio;
                
                ComputeVelocitiesSubIntervals(s, startHeightsLeft, endHeightsLeft, numLeftHeights, numRightHeights, firstHeight,
                    m_leftGridLineIndex, m_rightGridLineIndex, numPerpendicularFacesOnSubintervalAndEdge, edgeVelocities, hh0LeftMaxRatio);

                int startHeightsRight = m_rightGridLineIndex[s];
                int endHeightsRight = startHeightsRight + m_mfac[s] - 1;
                double hh0RightMaxRatio;

                ComputeVelocitiesSubIntervals(s, startHeightsRight, endHeightsRight, numLeftHeights, numRightHeights, firstHeight,
                    m_rightGridLineIndex, m_leftGridLineIndex, numPerpendicularFacesOnSubintervalAndEdge, edgeVelocities, hh0RightMaxRatio);

                // re-evaluate if growing grid outside is needed

                if(numLeftHeights ==0 && numRightHeights<=1 || 
                    numRightHeights == 0 && numLeftHeights <= 1 ||
                    numLeftHeights== numRightHeights==1)
                {
                    m_growGridOutside = true;
                }

                // left part
                int numNLeftExponential = 0;
                if(m_growGridOutside)
                {
                    numNLeftExponential = std::min(ComputeNumberExponentialIntervals(hh0LeftMaxRatio), m_maximumNumGridNodesAlongN);
                }
                for (int i = startHeightsLeft; i < endHeightsLeft; ++i)
                {
                    numPerpendicularFacesOnSubintervalAndEdge[2][i] = numNLeftExponential;
                }

                // right part
                int numNRightExponential = 0;
                if (m_growGridOutside)
                {
                    numNRightExponential = std::min(ComputeNumberExponentialIntervals(hh0RightMaxRatio), m_maximumNumGridNodesAlongN);
                }
                for (int i = startHeightsRight; i < endHeightsRight; ++i)
                {
                    numPerpendicularFacesOnSubintervalAndEdge[2][i] = numNRightExponential;
                }
            }

            return true;
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
            const int startHeights,
            const int endHeights,
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
                double maxHeight = *std::max_element(m_gridHeights[0].begin() + startHeights, m_gridHeights[0].begin() + endHeights);

                int numNUniformPart = std::floor(maxHeight / firstHeight + 0.99999);
                numNUniformPart = std::min(numNUniformPart, m_maxNUniformPart);

                for (int i = startHeights; i < endHeights; ++i)
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
                for (int i = startHeights; i < endHeights; ++i)
                {
                    numPerpendicularFacesOnSubintervalAndEdge[0][i] = numNUniformPart;
                    edgeVelocities[i] = firstHeight;

                    //compare with other side of spline
                    int otherSideIndex = gridLineIndex[s] + m_mfac[s] - (i - otherGridLineIndex[s] + 1);

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

                    for (int j = 1; j < m_maxNumHeights; ++j)
                    {
                        m_gridHeights[j][i] = m_gridHeights[j - 1][i];
                    }

                    for (int j = startHeights; j < endHeights; ++j)
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
            m_gridHeights.resize(m_maxNumHeights, std::vector<double>(m_maximumNumMeshNodesAlongM));
            std::fill(m_gridHeights.begin(), m_gridHeights.end(), std::vector<double>(m_maximumNumMeshNodesAlongM, 0.0));

            std::vector<std::vector<double>> heightsLeft(m_maxNumHeights, std::vector<double>(m_maximumNumMeshNodesAlongM, 0.0));
            std::vector<std::vector<double>> heightsRight(m_maxNumHeights, std::vector<double>(m_maximumNumMeshNodesAlongM, 0.0));
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

                int mfac = m_mfac[s];

                // Get the minimum number of sub-intervals in the cross splines for this center spline
                int minNumLeftIntervals = *std::min_element(m_numCrossSplineLeftHeights[s].begin(), m_numCrossSplineLeftHeights[s].end());
                int minNumRightIntervals = *std::min_element(m_numCrossSplineRightHeights[s].begin(), m_numCrossSplineRightHeights[s].end());

                std::fill(heightsLeft[0].begin(), heightsLeft[0].begin() + mfac, m_maximumGridHeights[s]);
                std::fill(heightsRight[0].begin(), heightsRight[0].begin() + mfac, m_maximumGridHeights[s]);

                if (m_numCrossingSplines[s] == 1)
                {
                    // only one crossing spline present: 
                    for (int i = 0; i < minNumLeftIntervals; ++i)
                    {
                        std::fill(heightsLeft[i].begin(), heightsLeft[i].begin() + mfac, m_crossSplineRightHeights[s][i][0]);

                    }
                    for (int i = 0; i < minNumRightIntervals; ++i)
                    {
                        std::fill(heightsRight[i].begin(), heightsRight[i].begin() + mfac, m_crossSplineLeftHeights[s][i][0]);
                    }
                }
                else
                {
                    int leftGridLineIndex = m_leftGridLineIndex[s];
                    edgesCenterPoints[0] = GetSplineLength(s, 0, m_gridLineDimensionalCoordinates[leftGridLineIndex]);
                    for (int i = 0; i < mfac; ++i)
                    {
                        edgesCenterPoints[i + 1] = edgesCenterPoints[i] + GetSplineLength(s, m_gridLineDimensionalCoordinates[leftGridLineIndex], m_gridLineDimensionalCoordinates[leftGridLineIndex + i]);
                    }

                    // compute at edge center points
                    for (int i = 0; i < mfac; ++i)
                    {
                        edgesCenterPoints[i] = 0.5 *(edgesCenterPoints[i] + edgesCenterPoints[i + 1]);
                    }
                    edgesCenterPoints[mfac] = doubleMissingValue;

                    //compute center spline path length of cross splines
                    crossingSplinesDimensionalCoordinates[0] = GetSplineLength(s, 0.0, m_crossSplineCoordinates[s][0]);
                    for (int i = 0; i < m_numCrossingSplines[s] - 1; ++i)
                    {
                        crossingSplinesDimensionalCoordinates[i + 1] = crossingSplinesDimensionalCoordinates[i] + GetSplineLength(s, m_crossSplineCoordinates[s][i], m_crossSplineCoordinates[s][i + 1]);
                        numHeightsLeft[i] = m_numCrossSplineLeftHeights[s][i];
                        numHeightsRight[i] = m_numCrossSplineRightHeights[s][i];
                    }


                    for (int j = 0; j < m_maxNumHeights; ++j)
                    {
                        AddValueToVector(numHeightsLeft, -1);
                        AddValueToVector(numHeightsRight, -1);

                        FindNearestCrossSplines(s, j,
                            numHeightsLeft,
                            edgesCenterPoints,
                            localValidSplineIndexes,
                            localSplineDerivatives,
                            crossingSplinesDimensionalCoordinates,
                            heightsLeft);

                        FindNearestCrossSplines(s, j,
                            numHeightsRight,
                            edgesCenterPoints,
                            localValidSplineIndexes,
                            localSplineDerivatives,
                            crossingSplinesDimensionalCoordinates,
                            heightsRight);
                    }
                }

                
                // store grid height
                for (int j = 0; j < m_maxNumHeights; ++j)
                {
                    for (int i = 0; i < m_mfac[s]; ++i)
                    {
                        m_gridHeights[j][m_leftGridLineIndex[s] + i - 1] = heightsLeft[j][i];
                    }
                }
            }

            return true;
        }

        bool FindNearestCrossSplines(const int s, 
            const int j,
            const std::vector<int>& numHeightsLeft,
            const std::vector<double>& edgesCenterPoints,
            std::vector<int>& localValidSplineIndexes,
            std::vector<double>& localSplineDerivatives,
            std::vector<double>& crossingSplinesDimensionalCoordinates,
            std::vector<std::vector<double>>& heights)
        {
            int numValid;
            GetValidSplineIndexses(s, m_numCrossingSplines[s], numHeightsLeft, localValidSplineIndexes, numValid);

            // no subheights to compute 
            if (numValid == 0)
            {
                return false;
            }

            int mfac = m_mfac[s];
            SecondOrderDerivative(m_crossSplineLeftHeights[s][j], m_numCrossingSplines[s], localSplineDerivatives);

            crossingSplinesDimensionalCoordinates[0] = GetSplineLength(s, 0.0, m_crossSplineCoordinates[s][0]);
            for (int i = 0; i < mfac; ++i)
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

                Interpolate(m_crossSplineLeftHeights[s][j], localSplineDerivatives, factor, heights[j][i]);
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
            m_gridLine.resize(numCenterSplines*(m_maximumNumMeshNodesAlongM + 2));
            m_gridLineDimensionalCoordinates.resize(numCenterSplines*(m_maximumNumMeshNodesAlongM + 2));

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

                gridLineIndex = gridLineIndex + numM + 1;
                m_gridLine[gridLineIndex] = Point{ doubleMissingValue,  doubleMissingValue };
                m_gridLineDimensionalCoordinates[gridLineIndex] = doubleMissingValue;
                m_mfac[s] = numM;
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
            numM = std::min(numM, m_maximumNumMeshNodesAlongM);

            double endSplineAdimensionalCoordinate = m_numSplineNodes[splineIndex] - 1;
            double splineLength = GetSplineLength(splineIndex, 0.0, endSplineAdimensionalCoordinate, 10, isSpacingCurvatureAdapeted, m_maximumGridHeights[splineIndex]);
            
            gridLine[startingIndex] = m_splineCornerPoints[splineIndex][0];
            FuncDimensionalToAdimensionalDistance func(*this, splineIndex, isSpacingCurvatureAdapeted, m_maximumGridHeights[splineIndex]);
            
            double currentMaxWidth = std::numeric_limits<double>::max();
            while (currentMaxWidth > m_averageMeshWidth && numM < m_maximumNumMeshNodesAlongM)
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
                    numM = std::min(std::max(int(m_maximumNumMeshNodesAlongM / m_maximumGridHeights[splineIndex] *numM), numM + 1), m_maximumNumMeshNodesAlongM);
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
                    m_mfac[s] = m_mfacOriginal[s];
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
            m_crossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex].resize(m_maxNumHeights, 0);
            for (int s = centerSplineLocalIndex; s <  m_numCrossingSplines[crossingSplineIndex] - 1; ++s)
            {
                if (numSubIntervalsRight >= m_maxNumHeights)
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
            m_crossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex].resize(m_maxNumHeights, 0);
            for (int s = centerSplineLocalIndex; s >= 1; --s)
            {
                if (numSubIntervalsLeft >= m_maxNumHeights)
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
        int m_maximumNumMeshNodesAlongM = 20;                                           // mfacmax
        int m_maximumNumGridNodesAlongN = 40;                                           //  N - refinement factor for regular grid generation.
        int m_maxNumHeights = 10;                                                       // Nsubmax maximum number of subintervals of grid layers, each having their own exponential grow factor
        double m_averageMeshWidth = 500.0;
        int m_allocationSize = 5;                                                       // allocation cache size
        double m_dtolcos = 0.95;                                                        // minimum allowed absolute value of crossing - angle cosine
        
        double m_aspectRatio = 0.1;                                                     // daspect aspect ratio
        double m_aspectRatioGrowFactor = 1.1;                                           // dgrow grow factor of aspect ratio
        double m_gridLayerHeight0 = 10.0;                                               // dheight0 grid layer height
        double m_maxaspect = 1.0;                                                       // maxaspect maximum cell aspect ratio *inoperative*
        int m_maxNUniformPart = 5;                                                      // maximum number of layers in the uniform part
        bool m_growGridOutside = true;                                                  // grow the grid outside the prescribed grid height
        

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
        std::vector<int> m_mfac;                                                        // mfac number of grid intervals on the spline
        std::vector<std::vector<int>> m_nfacL;                                          // nfacL number of grid layers in each subinterval at the left - hand side of the spline * not used yet*
        std::vector<std::vector<int>> m_nfacR;                                          // nfacR number of grid layers in each subinterval at the right - hand side of the spline * not used yet*
        std::vector<int> m_leftGridLineIndex;                                           // iL index in the whole gridline array of the first grid point on the left - hand side of the spline
        std::vector<int> m_rightGridLineIndex;                                          // iR index in the whole gridline array of the first grid point on the right - hand side of the spline
        std::vector<Point> m_gridLine;                                                  // xg1, yg1 coordinates of the first gridline
        std::vector<double> m_gridLineDimensionalCoordinates;                           // sg1 center spline coordinates of the first gridline
        std::vector<std::vector<double>> m_gridHeights;


        //original spline chaches
        int m_numOriginalSplines = 0;
        std::vector<int> m_leftGridLineIndexOriginal;
        std::vector<int> m_rightGridLineIndexOriginal;
        std::vector<int> m_mfacOriginal;
        std::vector<double> m_maximumGridHeightsOriginal;
        std::vector<int> m_numLayersOriginal;
    
    };

}