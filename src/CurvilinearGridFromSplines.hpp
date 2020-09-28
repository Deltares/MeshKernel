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
#include "Entities.hpp"
#include "CurvilinearParametersNative.hpp"
#include "SplinesToCurvilinearParametersNative.hpp"
#include "Splines.hpp"

namespace meshkernel
{
    class CurvilinearGrid;
    struct SplinesToCurvilinearParametersNative;

    class CurvilinearGridFromSplines
    {
    public:

        /// <summary>
        /// Ctor
        /// </summary>
        /// <param name="curvilinearParametersNative">The parameters for OrthogonalCurvilinearGridFromSplines algoritm</param>
        /// <param name="splinesToCurvilinearParametersNative">The parameters for OrthogonalCurvilinearGridFromSplines algoritm</param>
        /// <returns></returns>
        CurvilinearGridFromSplines( std::shared_ptr<Splines> splines,
                                    const meshkernelapi::CurvilinearParametersNative& curvilinearParametersNative,
                                    const meshkernelapi::SplinesToCurvilinearParametersNative& splinesToCurvilinearParametersNative);

        /// <summary>
        /// Computes the spline properties, such as cross splines (get_splineprops)
        /// </summary>
        /// <param name="restoreOriginalProperties"></param>
        /// <returns></returns>
        bool ComputeSplineProperties(bool restoreOriginalProperties);

        /// <summary>
        /// Computes a curvilinear grid using the growing front method (spline2curvi). OrthogonalCurvilinearGrid algorithm.
        /// 1. Compute spline properties (the crossings)
        /// 2. Make all grid lines of the central spline
        /// 3. Adds artificial splines
        /// 4. Compute properties with artificial splines added
        /// 5. Compute the edge velocities
        /// 6. Grow layers
        /// </summary>
        /// <param name="curvilinearGrid">The computed curvilinear grid</param>
        /// <returns></returns>
        bool Compute(CurvilinearGrid& curvilinearGrid);

        /// <summary>
        /// Initialize the OrthogonalCurvilinearGrid algorithm.
        /// </summary>
        /// <returns></returns>
        bool Initialize();

        /// <summary>
        /// Performs one iteration for generating another layer on the advancing fronts
        /// </summary>
        /// <param name="layer">The index of the layer to be generated</param>
        /// <returns></returns>
        bool Iterate(int layer);

        /// <summary>
        /// Get the curvilinear grid
        /// </summary>
        /// <param name="curvilinearGrid"></param>
        /// <returns></returns>
        bool ComputeCurvilinearGrid(CurvilinearGrid& curvilinearGrid);

        /// <summary>
        /// For the central spline, computes the spline subdivisions along the spline (make_wholegridline)
        /// </summary>
        /// <returns></returns>
        bool MakeAllGridLines();

        std::vector<Point> m_gridLine;                           // coordinates of the first gridline (xg1, yg1)
        std::vector<double> m_gridLineDimensionalCoordinates;    // center spline coordinates of the first gridline (sg1)
        std::vector<double> m_maximumGridHeights;                // maximum transversal grid height ()
        size_t m_numM = 0;                                       // Number of columns in the curvilinear grid
        
        std::shared_ptr<Splines> m_splines;                      // A pointer to spline

    private:

        /// <summary>
        /// From the layer index gets the previous grid layer where start growing and the transversal sublayer index (get_isub)
        /// </summary>
        /// <param name="layer">The current layer</param>
        /// <param name="gridLayer">The next grid layer</param>
        /// <param name="subLayerIndex">The transversal sub-layer index</param>
        /// <returns></returns>
        bool GetSubIntervalAndGridLayer(int layer, 
                                        int& gridLayer, 
                                        int& subLayerIndex);

        /// <summary>
        /// Grow layer at layer index
        /// </summary>
        /// <param name="layerIndex">The layer index to grow</param>
        /// <returns></returns>
        bool GrowLayer(int layerIndex);

        /// <summary>
        /// Compute the maximum allowable grid layer growth time self crossings (comp_tmax_self)
        /// </summary>
        /// <param name="coordinates">The coordinates to grow</param>
        /// <param name="velocities">The velocities</param>
        /// <param name="maximumGridLayerGrowTime">The maximum grow layer time</param>
        /// <returns></returns>
        bool ComputeMaximumGridLayerGrowTime(const std::vector<Point>& coordinates,
                                             const std::vector<Point>& velocities,
                                             std::vector<double>& maximumGridLayerGrowTime) const;

        /// <summary>
        /// Copy growth velocities to the advancing front, add points at front corners corners (copy_vel_to_front)
        /// </summary>
        /// <param name="layerIndex"></param>
        /// <param name="previousVelocities"></param>
        /// <param name="numFrontPoints"></param>
        /// <param name="gridPointsIndexses"></param>
        /// <param name="frontGridPoints"></param>
        /// <param name="velocities"></param>
        /// <returns></returns>
        bool CopyVelocitiesToFront( const int layerIndex,
                                    const std::vector<Point>& previousVelocities,
                                    int& numFrontPoints,
                                    std::vector<std::vector<int>>& gridPointsIndexses,
                                    std::vector<Point>& frontGridPoints,
                                    std::vector<Point>& velocities );

        /// <summary>
        /// Computes the points at front, which have to be moved.
        /// </summary>
        /// <param name="gridPointsIndexses"></param>
        /// <param name="frontGridPoints"></param>
        /// <param name="numFrontPoints"></param>
        /// <returns></returns>
        bool FindFront(std::vector<std::vector<int>>& gridPointsIndexses,
                       std::vector<Point>& frontGridPoints,
                       int& numFrontPoints);

        /// <summary>
        /// Compute growth velocity vectors at grid points (comp_vel)
        /// </summary>
        /// <param name="layerIndex"></param>
        /// <param name="velocityVector"></param>
        /// <returns></returns>
        bool ComputeVelocitiesAtGridPoints(int layerIndex, std::vector<Point>& velocityVector);

        /// <summary>
        /// Get left and right points at given layer for a given index (get_LR)
        /// </summary>
        /// <param name="gridPoints">The layer</param>
        /// <param name="index"></param>
        /// <param name="currentLeftIndex"></param>
        /// <param name="currentRightIndex"></param>
        /// <returns></returns>
        bool GetNeighbours(const std::vector<Point>& gridPoints, 
                           int index, 
                           int& currentLeftIndex, 
                           int& currentRightIndex) const;

        /// <summary>
        /// Compute the edge grow velocities (comp_edgevel)
        /// TODO: can this be splitted in compute heights and computeGrowFactors
        /// </summary>
        /// <param name="edgeVelocities"></param>
        /// <param name="growFactorOnSubintervalAndEdge"></param>
        /// <param name="numPerpendicularFacesOnSubintervalAndEdge"></param>
        /// <returns></returns>
        bool ComputeEdgeVelocities(std::vector<double>& edgeVelocities, 
                                   std::vector<std::vector<double>>& growFactorOnSubintervalAndEdge, 
                                   std::vector<std::vector<int>>& numPerpendicularFacesOnSubintervalAndEdge);

        /// <summary>
        /// Compute the grid grow factor for a given total grid height, first grid layer height and number of grid layers (comp_dgrow)
        /// </summary>
        /// <param name="totalGridHeight"></param>
        /// <param name="firstGridLayerHeight"></param>
        /// <param name="numberOfGridLayers"></param>
        /// <param name="result"></param>
        /// <returns></returns>
        bool ComputeGrowFactor(double totalGridHeight,
                               double firstGridLayerHeight,
                               int numberOfGridLayers,
                               double& result);

        /// <summary>
        /// Computes the exponential grid height
        /// </summary>
        /// <param name="aspectRatioGrowFactor"></param>
        /// <param name="firstGridLayerHeights"></param>
        /// <param name="numberOfGridLayers"></param>
        /// <returns></returns>
        double ComputeTotalExponentialHeight(double aspectRatioGrowFactor, 
                                             double firstGridLayerHeights, 
                                             int numberOfGridLayers) const;

        /// <summary>
        /// Compute the number of grid layers for a given grow factor, first grid layer height and total grid height (comp_nfac)
        /// </summary>
        /// <param name="hhMaxRatio"></param>
        /// <returns></returns>
        int ComputeNumberExponentialIntervals(double hhMaxRatio) const;

        /// <summary>
        /// Computes the sub-interval velocities (left and right)
        /// </summary>
        /// <param name="s"></param>
        /// <param name="startGridLineIndex"></param>
        /// <param name="endGridLineIndex"></param>
        /// <param name="numHeights"></param>
        /// <param name="numOtherSideHeights"></param>
        /// <param name="firstHeight"></param>
        /// <param name="gridLineIndex"></param>
        /// <param name="otherGridLineIndex"></param>
        /// <param name="numPerpendicularFacesOnSubintervalAndEdge"></param>
        /// <param name="edgeVelocities"></param>
        /// <param name="hh0MaxRatio"></param>
        /// <returns></returns>
        bool ComputeVelocitiesSubIntervals( int s,
                                            int startGridLineIndex,
                                            int endGridLineIndex,
                                            int numHeights,
                                            int numOtherSideHeights,
                                            double firstHeight,
                                            const std::vector<int>& gridLineIndex,
                                            const std::vector<int>& otherGridLineIndex,
                                            std::vector<std::vector<int>>& numPerpendicularFacesOnSubintervalAndEdge,
                                            std::vector<double>& edgeVelocities,
                                            double& hh0MaxRatio);

        /// <summary>
        /// Compute the grid heights at grid edges on the center spline (comp_gridheights)
        /// </summary>
        /// <returns></returns>
        bool ComputeGridHeights();

        /// <summary>
        /// Find the nearest cross spline
        /// </summary>
        /// <param name="s"></param>
        /// <param name="j"></param>
        /// <param name="numHeightsLeft"></param>
        /// <param name="edgesCenterPoints"></param>
        /// <param name="crossSplineLeftHeights"></param>
        /// <param name="localValidSplineIndexes"></param>
        /// <param name="localSplineDerivatives"></param>
        /// <param name="crossingSplinesDimensionalCoordinates"></param>
        /// <param name="heights"></param>
        /// <returns></returns>
        bool FindNearestCrossSplines(int s,
                                     int j,
                                     const std::vector<int>& numHeightsLeft,
                                     const std::vector<double>& edgesCenterPoints,
                                     const std::vector<std::vector<double>>& crossSplineLeftHeights,
                                     std::vector<int>& localValidSplineIndexes,
                                     std::vector<double>& localSplineDerivatives,
                                     std::vector<double>& crossingSplinesDimensionalCoordinates,
                                     std::vector<std::vector<double>>& heights);

        /// <summary>
        /// Gets the valid spline indexses
        /// </summary>
        /// <param name="numValues"></param>
        /// <param name="v"></param>
        /// <param name="validIndexses"></param>
        /// <param name="numValid"></param>
        /// <returns></returns>
        bool GetValidSplineIndexses(int numValues, 
                                    const std::vector<int>& v, 
                                    std::vector<int>& validIndexses, 
                                    int& numValid) const;

        /// <summary>
        /// Computes the intersection of two splines, one must have only two nodes (get_crosssplines)
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        bool GetSplineIntersections(int index);

        /// <summary>
        /// generate a gridline on a spline with a prescribed maximum mesh width (make_gridline)
        /// </summary>
        /// <param name="splineIndex"></param>
        /// <param name="startingIndex"></param>
        /// <param name="gridLine"></param>
        /// <param name="adimensionalCoordinates"></param>
        /// <param name="numM"></param>
        /// <returns></returns>
        bool MakeGridLine(int splineIndex,
                          int startingIndex,
                          std::vector<Point>& gridLine,
                          std::vector<double>& adimensionalCoordinates,
                          int& numM);



        /// <summary>
        /// Compute the grid heights using ComputeSubHeights and calculates the maximum sub height (get_heights)
        /// </summary>
        /// <returns></returns>
        bool ComputeHeights();

        /// <summary>
        /// Compute the heights at left and right of the center spline index. 
        /// Heights are determined by the lenghts of the left and right parts of the crossing splines (comp_subheights)
        /// </summary>
        /// <param name="centerSplineIndex"></param>
        /// <param name="crossingSplineLocalIndex"></param>
        /// <returns></returns>
        bool ComputeSubHeights(int centerSplineIndex, int crossingSplineLocalIndex);
        
        /// <summary>
        /// Computes curvature in a point on a spline (comp_curv)
        /// </summary>
        /// <param name="splineIndex"></param>
        /// <param name="adimensionalPointCoordinate"></param>
        /// <param name="curvatureFactor"></param>
        /// <param name="normalVector"></param>
        /// <param name="tangentialVector"></param>
        /// <returns></returns>
        bool ComputeCurvatureOnSplinePoint(int splineIndex,
                                           double adimensionalPointCoordinate,
                                           double& curvatureFactor,
                                           Point& normalVector,
                                           Point& tangentialVector);

        /// <summary>
        /// Remove skewed cells and cells whose aspect ratio exceeds a prescibed value (postgrid)
        /// </summary>
        /// <returns></returns>
        bool RemoveSkinnyTriangles();

        /// <summary>
        /// Delete a spline
        /// </summary>
        /// <param name="splineIndex">The spline index to delete</param>
        /// <returns></returns>
        bool DeleteSpline(int splineIndex);

        /// <summary>
        /// Allocate spline properties arrays 
        /// </summary>
        /// <returns></returns>
        bool AllocateSplinesProperties();

        /// <summary>
        /// The spline type
        /// </summary>
        enum class SplineTypes
        {
            central,
            crossing,
            arficial,
            lateral
        };


        // algorithm parameters
        meshkernelapi::CurvilinearParametersNative          m_curvilinearParametersNative;
        meshkernelapi::SplinesToCurvilinearParametersNative m_splinesToCurvilinearParametersNative;

        const int m_maxNumCenterSplineHeights = 10;                                     // Nsubmax, naz number of different heights a cross spline can have (is determined by how many crossing spline the user can input)
        const int m_maxNUniformPart = 5;                                                // maximum number of layers in the uniform part
        bool m_growGridOutside = true;                                                  // grow the grid outside the prescribed grid height
        double m_onTopOfEachOtherSquaredTolerance;                                      // On top of each other squared tolerance
        size_t m_numOriginalSplines = 0;                                                // The original number of splines
        int m_allocationSize = 5;                                                       // allocation cache size

        // Spline properties (first index is the spline number) 
        std::vector<SplineTypes> m_type;
        std::vector<int> m_centralSplineIndex;                                          // for each spline the index to its central
        std::vector<int> m_numCrossingSplines;                                          // ncs num of cross splines
        std::vector<std::vector<int>> m_crossingSplinesIndexses;                        // ics for each cross spline, the indexses of the center splines
  
        std::vector<std::vector<bool>> m_isLeftOriented;                                // isLeftOriented cross spline is left to right(.true.) or not (.false.) w.r.t.center spline
        std::vector<std::vector<double>>  m_crossSplineCoordinates;                     // t center spline coordinates of cross splines
        std::vector<std::vector<double>> m_cosCrossingAngle;                            // cosPhi cosine of crossing angle
        std::vector<std::vector<std::vector<double>>> m_crossSplineLeftHeights;         // hL left - hand side grid heights at cross spline locations for each grid layer subinterval, hL(1, :) being the height of the first subinterval, etc.
        std::vector<std::vector<std::vector<double>>> m_crossSplineRightHeights;        // hR right - hand side grid heights at cross spline locations for each grid layer subinterval, hR(1, :) being the height of the first subinterval, etc.
        std::vector<std::vector<int>>   m_numCrossSplineLeftHeights;                    // NsubL number of subintervals of grid layers at cross spline locations at the left - hand side of the spline, each having their own exponential grow factor
        std::vector<std::vector<int>>  m_numCrossSplineRightHeights;                    // NsubR number of subintervals of grid layers at cross spline locations at the right - hand side of the spline, each having their own exponential grow factor
        std::vector<int> m_numMSplines;                                                  // mfac number of grid intervals on the spline
        std::vector<std::vector<int>> m_nfacL;                                          // nfacL number of grid layers in each subinterval at the left - hand side of the spline * not used yet*
        std::vector<std::vector<int>> m_nfacR;                                          // nfacR number of grid layers in each subinterval at the right - hand side of the spline * not used yet*
        std::vector<int> m_leftGridLineIndex;                                           // iL index in the whole gridline array of the first grid point on the left - hand side of the spline
        std::vector<int> m_rightGridLineIndex;                                          // iR index in the whole gridline array of the first grid point on the right - hand side of the spline
        std::vector<std::vector<double>> m_gridHeights;                                 // heights of all grid elements        

        std::vector<int> m_leftGridLineIndexOriginal;
        std::vector<int> m_rightGridLineIndexOriginal;
        std::vector<int> m_mfacOriginal;
        std::vector<double> m_maximumGridHeightsOriginal;
        std::vector<SplineTypes> m_originalTypes;

        //cache variables during iterations
        std::vector<double> m_edgeVelocities;
        std::vector<int> m_validFrontNodes;
        std::vector<std::vector<Point>> m_gridPoints;
        double m_timeStep = 1.0;
        std::vector<int> m_subLayerGridPoints;
        std::vector<std::vector<int>> m_numPerpendicularFacesOnSubintervalAndEdge;
        std::vector<std::vector<double>> m_growFactorOnSubintervalAndEdge;

    };
}


