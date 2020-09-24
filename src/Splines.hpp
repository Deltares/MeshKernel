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
#include "Polygons.hpp"
#include "CurvilinearParametersNative.hpp"
#include "SplinesToCurvilinearParametersNative.hpp"

namespace MeshKernel
{
    class CurvilinearGrid;
    
    class Splines
    {

    public:

        /// <summary>
        /// Ctor
        /// </summary>
        /// <returns></returns>
        Splines();

        /// <summary>
        /// Ctor, set projection
        /// </summary>
        /// <param name="projection">The map projection</param>
        /// <returns></returns>
        Splines(Projections projection);

        /// <summary>
        /// Adds a new spline to m_splineCornerPoints 
        /// </summary>
        /// <param name="splines">The spline corner points</param>
        /// <param name="start">The starting index in splines</param>
        /// <param name="size">The end index splines</param>
        /// <returns>If the method succeeded</returns>
        bool AddSpline(const std::vector<Point>& splines, size_t start, size_t size);

        /// <summary>
        /// Second order derivative at spline corner points
        /// </summary>
        /// <param name="splines">The spline corner points</param>
        /// <param name="numNodes">The number of corner points</param>
        /// <param name="coordinatesDerivatives">The second order derivative at corner points</param>
        /// <returns>If the method succeeded</returns>
        static bool SecondOrderDerivative(const std::vector<Point>& splines,
                                          int numNodes,
                                          std::vector<Point>& coordinatesDerivatives);

        /// <summary>
        /// Second order derivative at spline corner points (result in a flat array)
        /// </summary>
        /// <param name="coordinates">The spline corner point coordinate (x or y)</param>
        /// <param name="numNodes">The number of corner points</param>
        /// <param name="coordinatesDerivatives">The second order derivative at corner points (x derivative or y derivative)</param>
        /// <returns>If the method succeeded</returns>
        static bool SecondOrderDerivative(const std::vector<double>& coordinates,
                                          int numNodes,
                                          std::vector<double>& coordinatesDerivatives);

        /// <summary>
        /// Computes the intersection of two splines (sect3r)
        /// </summary>
        /// <param name="first">The index of the first spline</param>
        /// <param name="second">The index of the second spline</param>
        /// <param name="crossProductIntersection">The cross product of the intersection</param>
        /// <param name="intersectionPoint">The intersection point</param>
        /// <param name="firstSplineRatio">The ratio of the first spline length where the intersection occours</param>
        /// <param name="secondSplineRatio">The ratio of the second spline length where the intersection occours</param>
        /// <returns>If a valid intersection is found</returns>
        bool GetSplinesIntersection(int first,
                                    int second,
                                    double& crossProductIntersection,
                                    Point& intersectionPoint,
                                    double& firstSplineRatio,
                                    double& secondSplineRatio);

        /// <summary>
        /// Computes the spline length in s coordinates (GETDIS)
        /// </summary>
        /// <param name="index">The spline index</param>
        /// <param name="startIndex">Adimensional start spline</param>
        /// <param name="endIndex">Adimensional end spline</param>
        /// <param name="numSamples">How many intervals to use between the startIndex and endIndex</param>
        /// <param name="accountForCurvature">Accounting for curvature</param>
        /// <param name="height">When accounting for curvature, the height to use</param>
        /// <param name="assignedDelta">When larger than zero, the number of intervals the spline is devided when computing the length</param>
        /// <returns>The computed length</returns>
        double GetSplineLength(int index,
                               double startIndex,
                               double endIndex,
                               int numSamples = 100,
                               bool accountForCurvature = false,
                               double height = 1.0,
                               double assignedDelta = -1);

        /// <summary>
        /// Compute the points on a spline lying at certain distance
        /// </summary>
        /// <param name="index">The spline index</param>
        /// <param name="distances">The distances</param>
        /// <param name="points">The resulting point along the spline</param>
        /// <returns>If the method succeeded</returns>
        bool InterpolatePointsOnSpline(int index, 
            double maximumGridHeight, 
            bool isSpacingCurvatureAdapted, 
            const std::vector<double>& distances, 
            std::vector<Point>& points,
            std::vector<double>& adimensionalDistances);

        std::vector<std::vector<Point>> m_splineNodes;           // The spline corner points
        std::vector<std::vector<Point>> m_splineDerivatives;     // The spline derivatives at the corner points  
        std::vector<int> m_numSplineNodes;                       // Number of spline nodes in each spline
        std::vector<int> m_numAllocatedSplineNodes;              // Number of allocated node in each spline
        std::vector<double> m_splinesLength;                     // The length of each spline
        Projections m_projection;                                // The map projection  
        size_t m_numSplines = 0;                                    // Current number of splines

        int m_numAllocatedSplines = 0;                           // Total number of allocated splines
        int m_allocationSize = 5;                                // allocation cache size
        
    private:

        /// <summary>
        /// Adds a new corner point in an existing spline
        /// </summary>
        /// <param name="splineIndex">The spline index</param>
        /// <param name="point">The point to add</param>
        /// <returns>If the method succeeded</returns>
        bool AddPointInExistingSpline(int splineIndex, const Point& point);
        
        /// <summary>
        /// Computes curvature in a spline point (comp_curv)
        /// </summary>
        /// <param name="splineIndex">the spline index</param>
        /// <param name="adimensionalPointCoordinate">The adimensional coordinate of the point along the spline</param>
        /// <param name="curvatureFactor">The computed curvature factor</param>
        /// <param name="normalVector">The computed normal vector</param>
        /// <param name="tangentialVector">The computed tangential vector</param>
        /// <returns>If the method succeeded</returns>
        bool ComputeCurvatureOnSplinePoint(int splineIndex,
                                           double adimensionalPointCoordinate,
                                           double& curvatureFactor,
                                           Point& normalVector,
                                           Point& tangentialVector);

        /// <summary>
        /// Delete a spline
        /// </summary>
        /// <param name="splineIndex">The index of the spline to delete</param>
        /// <returns>If the method succeeded</returns>
        bool DeleteSpline(int splineIndex);

        /// <summary>
        /// Allocate spline properties vectors 
        /// </summary>
        /// <returns>If the method succeeded</returns>
        bool AllocateSplinesProperties();

    };

    struct FuncDimensionalToAdimensionalDistance
    {
        FuncDimensionalToAdimensionalDistance(Splines* splines,
            int splineIndex,
            bool isSpacingCurvatureAdapted,
            double h) :
            m_spline(splines),
            m_splineIndex(splineIndex),
            m_isSpacingCurvatureAdapted(isSpacingCurvatureAdapted),
            m_h(h)
        {
        };

        void SetDimensionalDistance(double distance)
        {
            m_DimensionalDistance = distance;
        }

        // this is the function we want to find the root
        double operator()(double adimensionalDistancereferencePoint)
        {
            double distanceFromReferencePoint = m_spline->GetSplineLength(m_splineIndex, 0, adimensionalDistancereferencePoint, m_numSamples, m_isSpacingCurvatureAdapted, m_h, 0.1);
            distanceFromReferencePoint = std::abs(distanceFromReferencePoint - m_DimensionalDistance);
            return distanceFromReferencePoint;
        }

        Splines* m_spline;
        int m_splineIndex;
        bool m_isSpacingCurvatureAdapted;
        double m_h;
        int m_numSamples = 10;
        double m_DimensionalDistance;
    };

}


