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
#include <memory>
#include <MeshKernel/SplinesToCurvilinearParametersNative.hpp>
#include <MeshKernel/CurvilinearParametersNative.hpp>
#include <MeshKernel/Entities.hpp>

namespace meshkernel
{
    class CurvilinearGrid;
    class Splines;

    class CurvilinearGridFromSplines
    {
    public:
        /// @brief Ctor
        /// @param curvilinearParametersNative The parameters for OrthogonalCurvilinearGridFromSplines algorithm
        /// @param splinesToCurvilinearParametersNative The parameters for OrthogonalCurvilinearGridFromSplines algorithm
        CurvilinearGridFromSplines(std::shared_ptr<Splines> splines,
                                   const meshkernelapi::CurvilinearParametersNative& curvilinearParametersNative,
                                   const meshkernelapi::SplinesToCurvilinearParametersNative& splinesToCurvilinearParametersNative);

        /// Computes the spline properties, such as cross splines (get_splineprops)
        /// @param restoreOriginalProperties
        void ComputeSplineProperties(bool restoreOriginalProperties);

        /// @brief Computes a curvilinear grid using the growing front method (spline2curvi).
        ///
        /// OrthogonalCurvilinearGrid algorithm.
        /// 1. Compute spline properties (the crossings)
        /// 2. Make all grid lines of the central spline
        /// 3. Adds artificial splines
        /// 4. Compute properties with artificial splines added
        /// 5. Compute the edge velocities
        /// 6. Grow layers
        /// @param curvilinearGrid The computed curvilinear grid
        void Compute(CurvilinearGrid& curvilinearGrid);

        /// @brief Initialize the OrthogonalCurvilinearGrid algorithm.
        void Initialize();

        /// @brief Performs one iteration for generating another layer on the advancing fronts
        /// @param layer The index of the layer to be generated
        void Iterate(int layer);

        /// @brief Get the curvilinear grid
        /// @param curvilinearGrid
        void ComputeCurvilinearGrid(CurvilinearGrid& curvilinearGrid);

        /// @brief For the central spline, computes the spline subdivisions along the spline (make_wholegridline)
        void MakeAllGridLines();

        std::vector<Point> m_gridLine;                        // coordinates of the first gridline (xg1, yg1)
        std::vector<double> m_gridLineDimensionalCoordinates; // center spline coordinates of the first gridline (sg1)
        std::vector<double> m_maximumGridHeights;             // maximum transversal grid height ()
        size_t m_numM = 0;                                    // Number of columns in the curvilinear grid

        std::shared_ptr<Splines> m_splines; // A pointer to spline

    private:
        /// @brief From the layer index gets the next grid layer and the transversal sublayer index (get_isub)
        /// @param[in] layer The current layer
        /// @param[out] gridLayer The next grid layer
        /// @param[out] subLayerIndex The transversal sub-layer index
        void GetSubIntervalAndGridLayer(int layer,
                                        int& gridLayer,
                                        int& subLayerIndex);

        /// @brief Grow layer at layer index
        /// @param layerIndex The layer index to grow
        void GrowLayer(int layerIndex);

        /// @brief Compute the maximum allowable grid layer growth time self crossings (comp_tmax_self)
        /// @param coordinates The coordinates to grow
        /// @param velocities The velocities
        /// @param maximumGridLayerGrowTime The maximum grow layer time
        void ComputeMaximumGridLayerGrowTime(const std::vector<Point>& coordinates,
                                             const std::vector<Point>& velocities,
                                             std::vector<double>& maximumGridLayerGrowTime) const;

        /// <summary>
        /// Copy growth velocities to the advancing front, add points at front corners corners (copy_vel_to_front)
        /// </summary>
        /// <param name="layerIndex"></param>
        /// <param name="previousVelocities"></param>
        /// <param name="numFrontPoints"></param>
        /// <param name="gridPointsIndices"></param>
        /// <param name="frontGridPoints"></param>
        /// <param name="velocities"></param>
        /// <returns></returns>
        void CopyVelocitiesToFront(const int layerIndex,
                                   const std::vector<Point>& previousVelocities,
                                   int& numFrontPoints,
                                   std::vector<std::vector<int>>& gridPointsIndices,
                                   std::vector<Point>& frontGridPoints,
                                   std::vector<Point>& velocities);

        /// <summary>
        /// Computes the points at front, which have to be moved.
        /// </summary>
        /// <param name="gridPointsIndices"></param>
        /// <param name="frontGridPoints"></param>
        /// <param name="numFrontPoints"></param>
        /// <returns></returns>
        void FindFront(std::vector<std::vector<int>>& gridPointsIndices,
                       std::vector<Point>& frontGridPoints,
                       int& numFrontPoints);

        /// @brief Compute growth velocity vectors at grid points (comp_vel)
        /// @param layerIndex
        /// @param velocityVector
        void ComputeVelocitiesAtGridPoints(int layerIndex, std::vector<Point>& velocityVector);

        /// <summary>
        /// Get left and right points at given layer for a given index (get_LR)
        /// </summary>
        /// <param name="gridPoints">The layer</param>
        /// <param name="index"></param>
        /// <param name="currentLeftIndex"></param>
        /// <param name="currentRightIndex"></param>
        /// <returns></returns>
        void GetNeighbours(const std::vector<Point>& gridPoints,
                           int index,
                           int& currentLeftIndex,
                           int& currentRightIndex) const;

        /// <summary>
        /// Compute the edge grow velocities (comp_edgevel)
        /// TODO: can this be split in compute heights and computeGrowFactors
        /// </summary>
        /// <param name="edgeVelocities"></param>
        /// <param name="growFactorOnSubintervalAndEdge"></param>
        /// <param name="numPerpendicularFacesOnSubintervalAndEdge"></param>
        /// <returns></returns>
        void ComputeEdgeVelocities(std::vector<double>& edgeVelocities,
                                   std::vector<std::vector<double>>& growFactorOnSubintervalAndEdge,
                                   std::vector<std::vector<int>>& numPerpendicularFacesOnSubintervalAndEdge);

        /// @brief Compute the grid grow factor for a given total grid height, first grid layer height and number of grid layers (comp_dgrow)
        /// @param[in] totalGridHeight
        /// @param[in] firstGridLayerHeight
        /// @param[in] numberOfGridLayers
        /// @param[out] result
        void ComputeGrowFactor(double totalGridHeight,
                               double firstGridLayerHeight,
                               int numberOfGridLayers,
                               double& result) const;

        /// <summary>
        /// Computes the exponential grid height
        /// </summary>
        /// <param name="aspectRatioGrowFactor"></param>
        /// <param name="firstGridLayerHeights"></param>
        /// <param name="numberOfGridLayers"></param>
        /// <returns></returns>
        [[nodiscard]] double ComputeTotalExponentialHeight(double aspectRatioGrowFactor,
                                                           double firstGridLayerHeights,
                                                           int numberOfGridLayers) const;

        /// <summary>
        /// Compute the number of grid layers for a given grow factor, first grid layer height and total grid height (comp_nfac)
        /// </summary>
        /// <param name="hhMaxRatio"></param>
        /// <returns></returns>
        [[nodiscard]] int ComputeNumberExponentialIntervals(double hhMaxRatio) const;

        /// @brief Computes the sub-interval velocities (left and right)
        /// @param[in] s
        /// @param[in] startGridLineIndex
        /// @param[in] endGridLineIndex
        /// @param[in] numHeights
        /// @param[in] numOtherSideHeights
        /// @param[in] firstHeight
        /// @param[in] gridLineIndex
        /// @param[in] otherGridLineIndex
        /// @param[out] numPerpendicularFacesOnSubintervalAndEdge
        /// @param[out] edgeVelocities
        /// @param[out] hh0MaxRatio
        void ComputeVelocitiesSubIntervals(int s,
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

        /// @brief Compute the grid heights at grid edges on the center spline (comp_gridheights)
        void ComputeGridHeights();

        /// @brief Find the nearest cross spline
        /// @param[in] s
        /// @param[in] j
        /// @param[in] numHeightsLeft
        /// @param[in] edgesCenterPoints
        /// @param[in] crossSplineLeftHeights
        /// @param[in] localValidSplineIndexes
        /// @param[out] localSplineDerivatives
        /// @param[out] crossingSplinesDimensionalCoordinates
        /// @param[out] heights
        void FindNearestCrossSplines(int s,
                                     int j,
                                     const std::vector<int>& numHeightsLeft,
                                     const std::vector<double>& edgesCenterPoints,
                                     const std::vector<std::vector<double>>& crossSplineLeftHeights,
                                     std::vector<int>& localValidSplineIndexes,
                                     std::vector<double>& localSplineDerivatives,
                                     std::vector<double>& crossingSplinesDimensionalCoordinates,
                                     std::vector<std::vector<double>>& heights);

        /// @brief Gets the valid spline indices
        /// @param[in] numValues
        /// @param[in] v
        /// @param[out] validIndices
        /// @param[out] numValid
        void GetValidSplineIndices(size_t numValues,
                                   const std::vector<int>& v,
                                   std::vector<int>& validIndices,
                                   size_t& numValid) const;

        /// <summary>
        /// Computes the intersection of two splines, one must have only two nodes (get_crosssplines)
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        void GetSplineIntersections(int index);

        /// @brief Generate a gridline on a spline with a prescribed maximum mesh width (make_gridline)
        /// @param[in] splineIndex
        /// @param[in] startingIndex
        /// @param[out] gridLine
        /// @param[out] adimensionalCoordinates
        /// @param[out] numM
        void MakeGridLine(int splineIndex,
                          int startingIndex,
                          std::vector<Point>& gridLine,
                          std::vector<double>& adimensionalCoordinates,
                          int& numM);

        /// @brief Compute the grid heights using ComputeSubHeights and calculates the maximum sub height (get_heights)
        void ComputeHeights();

        /// @brief Compute the heights at left and right of the center spline index.
        ///
        /// Heights are determined by the lengths of the left and right parts of the crossing splines (comp_subheights)
        ///
        /// @param centerSplineIndex
        /// @param crossingSplineLocalIndex
        void ComputeSubHeights(int centerSplineIndex, int crossingSplineLocalIndex);

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

        /// @brief Remove skewed cells and cells whose aspect ratio exceeds a prescibed value (postgrid)
        void RemoveSkinnyTriangles();

        /// <summary>
        /// Delete a spline
        /// </summary>
        /// <param name="splineIndex">The spline index to delete</param>
        /// <returns></returns>
        bool DeleteSpline(int splineIndex);

        /// @brief Allocate spline properties arrays
        void AllocateSplinesProperties();

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
        meshkernelapi::CurvilinearParametersNative m_curvilinearParametersNative;
        meshkernelapi::SplinesToCurvilinearParametersNative m_splinesToCurvilinearParametersNative;

        const int m_maxNumCenterSplineHeights = 10; // Nsubmax, naz number of different heights a cross spline can have (is determined by how many crossing spline the user can input)
        const int m_maxNUniformPart = 5;            // maximum number of layers in the uniform part
        bool m_growGridOutside = true;              // grow the grid outside the prescribed grid height
        double m_onTopOfEachOtherSquaredTolerance;  // On top of each other squared tolerance
        size_t m_numOriginalSplines = 0;            // The original number of splines
        int m_allocationSize = 5;                   // allocation cache size

        // Spline properties (first index is the spline number)
        std::vector<SplineTypes> m_type;
        std::vector<int> m_centralSplineIndex;                  // for each spline the index to its central
        std::vector<int> m_numCrossingSplines;                  // ncs num of cross splines
        std::vector<std::vector<int>> m_crossingSplinesIndices; // ics for each cross spline, the indices of the center splines

        std::vector<std::vector<bool>> m_isLeftOriented;                         // isLeftOriented cross spline is left to right(.true.) or not (.false.) w.r.t.center spline
        std::vector<std::vector<double>> m_crossSplineCoordinates;               // t center spline coordinates of cross splines
        std::vector<std::vector<double>> m_cosCrossingAngle;                     // cosPhi cosine of crossing angle
        std::vector<std::vector<std::vector<double>>> m_crossSplineLeftHeights;  // hL left - hand side grid heights at cross spline locations for each grid layer subinterval, hL(1, :) being the height of the first subinterval, etc.
        std::vector<std::vector<std::vector<double>>> m_crossSplineRightHeights; // hR right - hand side grid heights at cross spline locations for each grid layer subinterval, hR(1, :) being the height of the first subinterval, etc.
        std::vector<std::vector<int>> m_numCrossSplineLeftHeights;               // NsubL number of subintervals of grid layers at cross spline locations at the left - hand side of the spline, each having their own exponential grow factor
        std::vector<std::vector<int>> m_numCrossSplineRightHeights;              // NsubR number of subintervals of grid layers at cross spline locations at the right - hand side of the spline, each having their own exponential grow factor
        std::vector<int> m_numMSplines;                                          // mfac number of grid intervals on the spline
        std::vector<std::vector<int>> m_nfacL;                                   // nfacL number of grid layers in each subinterval at the left - hand side of the spline * not used yet*
        std::vector<std::vector<int>> m_nfacR;                                   // nfacR number of grid layers in each subinterval at the right - hand side of the spline * not used yet*
        std::vector<int> m_leftGridLineIndex;                                    // iL index in the whole gridline array of the first grid point on the left - hand side of the spline
        std::vector<int> m_rightGridLineIndex;                                   // iR index in the whole gridline array of the first grid point on the right - hand side of the spline
        std::vector<std::vector<double>> m_gridHeights;                          // heights of all grid elements

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
} // namespace meshkernel
