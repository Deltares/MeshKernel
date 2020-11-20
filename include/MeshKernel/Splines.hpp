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
#include <MeshKernel/Entities.hpp>

namespace meshkernel
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
        explicit Splines(Projections projection);

        /// @brief Adds a new spline to m_splineCornerPoints
        /// @param[in] splines The spline corner points
        /// @param[in] start The starting index in splines
        /// @param[in] size The end index splines
        void AddSpline(const std::vector<Point>& splines, size_t start, size_t size);

        /// @brief Second order derivative at spline corner points
        /// @param[in] splines The spline corner points
        /// @param[in] numNodes The number of corner points
        /// @param[out] coordinatesDerivatives The second order derivative at corner points
        static void SecondOrderDerivative(const std::vector<Point>& splines,
                                          size_t numNodes,
                                          std::vector<Point>& coordinatesDerivatives);

        /// @brief Second order derivative at spline corner point coordinates
        /// @param[in] coordinates The spline corner point coordinate (x or y)
        /// @param[in] numNodes The number of corner points
        /// @param[out] coordinatesDerivatives The second order derivative at corner points (x derivative or y derivative)
        static void SecondOrderDerivative(const std::vector<double>& coordinates,
                                          size_t numNodes,
                                          std::vector<double>& coordinatesDerivatives);

        /// @brief Computes the intersection of two splines (sect3r)
        /// @param[in] first The index of the first spline
        /// @param[in] second The index of the second spline
        /// @param[out] crossProductIntersection The cross product of the intersection
        /// @param[out] intersectionPoint The intersection point
        /// @param[out] firstSplineRatio The ratio of the first spline length where the intersection occurs
        /// @param[out] secondSplineRatio The ratio of the second spline length where the intersection occurs
        /// @returns If a valid intersection is found
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
        /// <param name="assignedDelta">When larger than zero, the number of intervals the spline is divided when computing the length</param>
        /// <returns>The computed length</returns>
        [[nodiscard]] double GetSplineLength(int index,
                                             double startIndex,
                                             double endIndex,
                                             int numSamples = 100,
                                             bool accountForCurvature = false,
                                             double height = 1.0,
                                             double assignedDelta = -1);

        /// @brief Compute the points on a spline lying at certain distance
        /// @param[in] index The spline index
        /// @param[in] maximumGridHeight Maximum grid height
        /// @param[in] isSpacingCurvatureAdapted Is spacing-curvature adapted
        /// @param[in] distances The distances
        /// @param[out] points The resulting point along the spline
        /// @param[out] adimensionalDistances Adimensional distances
        void InterpolatePointsOnSpline(int index,
                                       double maximumGridHeight,
                                       bool isSpacingCurvatureAdapted,
                                       const std::vector<double>& distances,
                                       std::vector<Point>& points,
                                       std::vector<double>& adimensionalDistances);

        /// @brief Get the number of splines
        /// @return the number of splines
        int GetNumSplines() const { return static_cast<int>(m_splineNodes.size()); }

        std::vector<std::vector<Point>> m_splineNodes;       // The spline corner points
        std::vector<std::vector<Point>> m_splineDerivatives; // The spline derivatives at the corner points
        std::vector<double> m_splinesLength;                 // The length of each spline
        Projections m_projection = Projections::cartesian;   // The map projection

    private:
        /// @brief Adds a new corner point in an existing spline
        /// @param[in] splineIndex The spline index
        /// @param[in] point The point to add
        void AddPointInExistingSpline(int splineIndex, const Point& point);

        /// @brief Computes curvature in a spline point (comp_curv)
        /// @param[in] splineIndex the spline index
        /// @param[in] adimensionalPointCoordinate The adimensional coordinate of the point along the spline
        /// @param[out] curvatureFactor The computed curvature factor
        /// @param[out] normalVector The computed normal vector
        /// @param[out] tangentialVector The computed tangential vector
        void ComputeCurvatureOnSplinePoint(int splineIndex,
                                           double adimensionalPointCoordinate,
                                           double& curvatureFactor,
                                           Point& normalVector,
                                           Point& tangentialVector);

        /// @brief Delete a spline
        /// @param[in] splineIndex The index of the spline to delete
        void DeleteSpline(int splineIndex);

        /// @brief Allocate spline properties vectors
        void AllocateSplinesProperties();
    };

    struct FuncDimensionalToAdimensionalDistance
    {
        FuncDimensionalToAdimensionalDistance(Splines* splines,
                                              int splineIndex,
                                              bool isSpacingCurvatureAdapted,
                                              double h) : m_spline(splines),
                                                          m_splineIndex(splineIndex),
                                                          m_isSpacingCurvatureAdapted(isSpacingCurvatureAdapted),
                                                          m_h(h){};

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
        double m_DimensionalDistance = 0.0;
    };

} // namespace meshkernel
