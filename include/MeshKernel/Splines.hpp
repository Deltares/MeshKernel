//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2021.
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

    /// @brief A class describing splines
    ///
    /// This class stores the corner points of each spline.
    /// Besides the corner points, the derivatives at the corner points are also stored.
    /// The coordinates of the points between the corner points are computed in the
    /// method Splines::InterpolatePointsOnSpline.
    class Splines
    {

    public:
        /// @brief Ctor
        /// @returns
        Splines() = default;

        /// @brief Ctor, set projection
        /// @brief projection The map projection
        /// @returns
        explicit Splines(Projection projection);

        /// @brief Adds a new spline to m_splineCornerPoints
        /// @param[in] splines The spline corner points
        /// @param[in] start The starting index in splines
        /// @param[in] size The end index splines
        void AddSpline(const std::vector<Point>& splines, size_t start, size_t size);

        /// @brief Second order derivative at spline corner points, from the start node to the end node of the spline (splint)
        /// @param[in] splines The spline corner points
        /// @param[in] startIndex The start spline node
        /// @param[in] endIndex The end spline node
        /// @returns coordinatesDerivatives The second order derivative at corner points
        [[nodiscard]] static std::vector<Point> SecondOrderDerivative(const std::vector<Point>& splines, size_t startIndex, size_t endIndex);

        /// @brief Second order derivative at spline corner point coordinates (splint)
        /// @param[in] coordinates The spline corner point coordinate (x or y)
        /// @param[in] startIndex The start spline node
        /// @param[in] endIndex The end spline node
        /// @returns coordinatesDerivatives The second order derivative at corner points (x derivative or y derivative)
        [[nodiscard]] static std::vector<double> SecondOrderDerivative(const std::vector<double>& coordinates, size_t startIndex, size_t endIndex);

        /// @brief Computes the intersection of two splines (sect3r)
        /// @param[in] first The index of the first spline
        /// @param[in] second The index of the second spline
        /// @param[out] crossProductIntersection The cross product of the intersection
        /// @param[out] intersectionPoint The intersection point
        /// @param[out] firstSplineRatio The ratio of the first spline length where the intersection occurs
        /// @param[out] secondSplineRatio The ratio of the second spline length where the intersection occurs
        /// @returns If a valid intersection is found
        bool GetSplinesIntersection(size_t first,
                                    size_t second,
                                    double& crossProductIntersection,
                                    Point& intersectionPoint,
                                    double& firstSplineRatio,
                                    double& secondSplineRatio);

        /// @brief Computes the spline length in s coordinates (GETDIS)
        /// @brief index The spline index
        /// @brief startIndex Adimensional start spline
        /// @brief endIndex Adimensional end spline
        /// @brief numSamples How many intervals to use between the startIndex and endIndex
        /// @brief accountForCurvature Accounting for curvature
        /// @brief height When accounting for curvature, the height to use
        /// @brief assignedDelta When larger than zero, the number of intervals the spline is divided when computing the length
        /// @returns The computed length
        [[nodiscard]] double GetSplineLength(size_t index,
                                             double startIndex,
                                             double endIndex,
                                             size_t numSamples = 100,
                                             bool accountForCurvature = false,
                                             double height = 1.0,
                                             double assignedDelta = -1.0);

        /// @brief Compute the points on a spline lying at certain distance
        /// @param[in] index The spline index
        /// @param[in] maximumGridHeight Maximum grid height
        /// @param[in] isSpacingCurvatureAdapted Is spacing-curvature adapted
        /// @param[in] distances The distances
        /// @param[out] points The resulting point along the spline
        /// @param[out] adimensionalDistances Adimensional distances
        void InterpolatePointsOnSpline(size_t index,
                                       double maximumGridHeight,
                                       bool isSpacingCurvatureAdapted,
                                       const std::vector<double>& distances,
                                       std::vector<Point>& points,
                                       std::vector<double>& adimensionalDistances);

        /// @brief Get the number of splines
        /// @return the number of splines
        auto GetNumSplines() const { return m_splineNodes.size(); }

        std::vector<std::vector<Point>> m_splineNodes;       ///< The spline corner points
        std::vector<std::vector<Point>> m_splineDerivatives; ///< The spline derivatives at the corner points
        std::vector<double> m_splinesLength;                 ///< The length of each spline
        Projection m_projection = Projection::cartesian;     ///< The map projection

    private:
        /// @brief Adds a new corner point in an existing spline
        /// @param[in] splineIndex The spline index
        /// @param[in] point The point to add
        void AddPointInExistingSpline(size_t splineIndex, const Point& point);

        /// @brief Computes curvature in a spline point (comp_curv)
        /// @param[in] splineIndex the spline index
        /// @param[in] adimensionalPointCoordinate The adimensional coordinate of the point along the spline
        /// @param[out] curvatureFactor The computed curvature factor
        /// @param[out] normalVector The computed normal vector
        /// @param[out] tangentialVector The computed tangential vector
        void ComputeCurvatureOnSplinePoint(size_t splineIndex,
                                           double adimensionalPointCoordinate,
                                           double& curvatureFactor,
                                           Point& normalVector,
                                           Point& tangentialVector);

        /// @brief Delete a spline
        /// @param[in] splineIndex The index of the spline to delete
        void DeleteSpline(size_t splineIndex);

        /// @brief Allocate spline properties vectors
        void AllocateSplinesProperties();
    };

    /// @brief This struct is used to create a function for converting an adimensional distance to a dimensional one
    struct FuncAdimensionalToDimensionalDistance
    {
        /// @brief Constructor
        /// @param[in] splines Pointer to splines
        /// @param[in] splineIndex Spline index
        /// @param[in] isSpacingCurvatureAdapted Is spacing curvature adapted
        /// @param[in] h When accounting for curvature, the height to use
        FuncAdimensionalToDimensionalDistance(Splines* splines,
                                              size_t splineIndex,
                                              bool isSpacingCurvatureAdapted,
                                              double h) : m_spline(splines),
                                                          m_splineIndex(splineIndex),
                                                          m_isSpacingCurvatureAdapted(isSpacingCurvatureAdapted),
                                                          m_h(h){};
        /// @brief Set dimensional distance
        /// @param[in] distance Distance
        void SetDimensionalDistance(double distance)
        {
            m_DimensionalDistance = distance;
        }

        /// @brief This is the function we want to find the root of
        double operator()(double adimensionalDistanceReferencePoint)
        {
            double distanceFromReferencePoint = m_spline->GetSplineLength(m_splineIndex, 0.0, adimensionalDistanceReferencePoint, m_numSamples, m_isSpacingCurvatureAdapted, m_h, 0.1);
            distanceFromReferencePoint = std::abs(distanceFromReferencePoint - m_DimensionalDistance);
            return distanceFromReferencePoint;
        }

        Splines* m_spline;                  ///< Pointer to splines
        size_t m_splineIndex;               ///< Spline index
        bool m_isSpacingCurvatureAdapted;   ///< Is spacing curvature adapted
        double m_h;                         ///< When accounting for curvature, the height to use
        size_t m_numSamples = 10;           ///< Number of samples
        double m_DimensionalDistance = 0.0; ///< Dimensional distance
    };

} // namespace meshkernel
