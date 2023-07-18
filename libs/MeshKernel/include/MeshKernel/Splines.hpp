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

#include "MeshKernel/Utilities/LinearAlgebra.hpp"
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/LandBoundary.hpp>
#include <MeshKernel/Operations.hpp>

namespace meshkernel
{
    class CurvilinearGrid;

    /// @brief A class describing splines
    ///
    /// This class stores the corner points of each spline.
    /// Besides the corner points, the derivatives at the corner points are also stored.
    /// The coordinates of the points between the corner points are computed in the
    /// method Splines::ComputePointOnSplineFromAdimensionalDistance.
    class Splines
    {

    public:
        /// @brief Default constructor
        Splines() = default;

        /// @brief Ctor, set projection
        /// @brief[in] projection   The map projection
        explicit Splines(Projection projection);

        /// @brief Ctor from grids, each gridline is converted to spline, first the first m_n horizontal lines then the m_m vertical lines
        /// @brief[in] grid         The curvilinear grid
        explicit Splines(CurvilinearGrid const& grid);

        /// @brief Adds a new spline to m_splineCornerPoints
        /// @param[in] splines The spline corner points
        /// @param[in] start The starting index in splines
        /// @param[in] size The end index splines
        void AddSpline(const std::vector<Point>& splines, UInt start, UInt size);

        /// @brief Second order derivative at spline corner points, from the start node to the end node of the spline (splint)
        /// @param[in] splines The spline corner points
        /// @param[in] startIndex The start spline node
        /// @param[in] endIndex The end spline node
        /// @returns coordinatesDerivatives The second order derivative at corner points
        [[nodiscard]] static std::vector<Point> SecondOrderDerivative(const std::vector<Point>& splines, UInt startIndex, UInt endIndex);

        /// @brief Second order derivative at spline corner point coordinates (splint)
        /// @param[in] coordinates The spline corner point coordinate (x or y)
        /// @param[in] startIndex The start spline node
        /// @param[in] endIndex The end spline node
        /// @returns coordinatesDerivatives The second order derivative at corner points (x derivative or y derivative)
        [[nodiscard]] static std::vector<double> SecondOrderDerivative(const std::vector<double>& coordinates, UInt startIndex, UInt endIndex);

        static lin_alg::ColumnVector<double> ComputeSplineWeights(const std::vector<Point>& splinePoints,
                                                                  const UInt totalNumberOfPoints,
                                                                  const Projection projection);

        /// @brief Evaluate a spline function (splint)
        static Point evaluate(const std::vector<Point>& coordinates, const std::vector<Point>& secondDerivative, const double evaluationPoint);

        // sample_spline
        static void sampleSpline(const std::vector<Point>& controlPoints,
                                 const size_t intermediatePointCount,
                                 std::vector<Point>& samplePoints);

        // comp_afinespline
        static void compAfinespline(const UInt n, const UInt numRef, UInt& count, lin_alg::MatrixColMajor<double>& mat);

        static lin_alg::MatrixColMajor<double> ComputeInterpolationMatrix(const lin_alg::MatrixColMajor<double>& splineCoefficients,
                                                                          const lin_alg::ColumnVector<double>& weights);

        // snap_spline
        void snapSpline(const LandBoundary& landBoundary,
                        const size_t splineIndex);

        /// @brief Computes the intersection of two splines (sect3r)
        /// @param[in] first The index of the first spline
        /// @param[in] second The index of the second spline
        /// @param[out] crossProductIntersection The cross product of the intersection
        /// @param[out] intersectionPoint The intersection point
        /// @param[out] firstSplineRatio The ratio of the first spline length where the intersection occurs
        /// @param[out] secondSplineRatio The ratio of the second spline length where the intersection occurs
        /// @returns If a valid intersection is found
        bool GetSplinesIntersection(UInt first,
                                    UInt second,
                                    double& crossProductIntersection,
                                    Point& intersectionPoint,
                                    double& firstSplineRatio,
                                    double& secondSplineRatio);

        /// @brief Computes the spline length in s coordinates (GETDIS)
        /// @param[in] index The spline index
        /// @param[in] startAdimensionalCoordinate Adimensional start spline
        /// @param[in] endAdimensionalCoordinate Adimensional end spline
        /// @param[in] numSamples How many intervals to use between the startAdimensionalCoordinate and endAdimensionalCoordinate
        /// @param[in] accountForCurvature Accounting for curvature
        /// @param[in] height When accounting for curvature, the height to use
        /// @param[in] assignedDelta When larger than zero, the number of intervals the spline is divided when computing the length
        /// @returns The computed length
        [[nodiscard]] double ComputeSplineLength(UInt index,
                                                 double startAdimensionalCoordinate,
                                                 double endAdimensionalCoordinate,
                                                 UInt numSamples = 100,
                                                 bool accountForCurvature = false,
                                                 double height = 1.0,
                                                 double assignedDelta = -1.0) const;

        /// @brief Compute the points on a spline lying at certain distance
        /// @param[in] index The spline index
        /// @param[in] maximumGridHeight Maximum grid height
        /// @param[in] isSpacingCurvatureAdapted Is spacing-curvature adapted
        /// @param[in] distances The dimensional distances of each point
        /// @returns The points along the splineand the adimensional distances of each point
        std::tuple<std::vector<Point>, std::vector<double>> ComputePointOnSplineFromAdimensionalDistance(UInt index,
                                                                                                         double maximumGridHeight,
                                                                                                         bool isSpacingCurvatureAdapted,
                                                                                                         const std::vector<double>& distances);

        /// @brief Computes the point on a spline segment which is the closest to another point
        /// @param[in] index The spline index
        /// @param[in] startSplineSegment The begin of the spline segment to consider
        /// @param[in] endSplineSegment The end of the spline segment to consider
        /// @param[in] point The point to account for in the calculation
        /// @returns The point on a spline segment which is the closest to the input point
        Point ComputeClosestPointOnSplineSegment(UInt index, double startSplineSegment, double endSplineSegment, Point point);

        /// @brief Get the number of splines
        /// @return the number of splines
        auto GetNumSplines() const { return static_cast<UInt>(m_splineNodes.size()); }

        std::vector<std::vector<Point>> m_splineNodes;       ///< The spline corner points
        std::vector<std::vector<Point>> m_splineDerivatives; ///< The spline derivatives at the corner points
        std::vector<double> m_splinesLength;                 ///< The length of each spline
        Projection m_projection = Projection::cartesian;     ///< The map projection

    private:
        /// @brief Compute the spline sample points.
        static std::tuple<lin_alg::ColumnVector<double>, lin_alg::ColumnVector<double>> ComputeSamplePoints (const std::vector<Point>& splinePoints,
                                                                                                             const lin_alg::MatrixColMajor<double>& aMatrix);


        /// @brief Adds a new corner point in an existing spline
        /// @param[in] splineIndex The spline index
        /// @param[in] point The point to add
        void AddPointInExistingSpline(UInt splineIndex, const Point& point);

        /// @brief Computes curvature in a spline point (comp_curv)
        /// @param[in] splineIndex the spline index
        /// @param[in] adimensionalPointCoordinate The adimensional coordinate of the point along the spline
        /// @returns The computed curvatureFactor, normal vector and tangential vector
        std::tuple<Point, Point, double> ComputeCurvatureOnSplinePoint(UInt splineIndex,
                                                                       double adimensionalPointCoordinate) const;

        /// @brief Delete a spline
        /// @param[in] splineIndex The index of the spline to delete
        void DeleteSpline(UInt splineIndex);

        /// @brief Allocate spline properties vectors
        void AllocateSplinesProperties();
    };

    /// @brief This structure is used to create a function for converting an adimensional distance on a spline to a dimensional one
    struct FuncAdimensionalToDimensionalDistanceOnSpline
    {
        /// @brief Constructor
        /// @param[in] splines A pointer to splines
        /// @param[in] splineIndex The index of the current spline
        /// @param[in] isSpacingCurvatureAdapted Is spacing curvature adapted
        /// @param[in] h When accounting for curvature, the height to use
        FuncAdimensionalToDimensionalDistanceOnSpline(Splines* splines,
                                                      UInt splineIndex,
                                                      bool isSpacingCurvatureAdapted,
                                                      double h) : m_spline(splines),
                                                                  m_splineIndex(splineIndex),
                                                                  m_isSpacingCurvatureAdapted(isSpacingCurvatureAdapted),
                                                                  m_h(h)
        {
            if (m_spline == nullptr)
            {
                throw std::invalid_argument("FuncAdimensionalToDimensionalDistanceOnSpline::m_spline is nullptr");
            }
        };

        /// @brief Set dimensional distance
        /// @param[in] distance distance
        void SetDimensionalDistance(double distance)
        {
            m_DimensionalDistance = distance;
        }

        /// @brief This is the function we want to find the root of
        double operator()(double adimensionalDistanceReferencePoint) const
        {
            auto distanceFromReferencePoint = m_spline->ComputeSplineLength(m_splineIndex, 0.0, adimensionalDistanceReferencePoint, m_numSamples, m_isSpacingCurvatureAdapted, m_h, 0.1);
            distanceFromReferencePoint = std::abs(distanceFromReferencePoint - m_DimensionalDistance);
            return distanceFromReferencePoint;
        }

        Splines* m_spline = nullptr;        ///< Pointer to splines
        UInt m_splineIndex;                 ///< Spline index
        bool m_isSpacingCurvatureAdapted;   ///< Is spacing curvature adapted
        double m_h;                         ///< When accounting for curvature, the height to use
        UInt m_numSamples = 10;             ///< Number of samples
        double m_DimensionalDistance = 0.0; ///< Dimensional distance
    };

    /// @brief This structure is used to compute the point on a spline closest to another point
    struct FuncDistanceFromAPoint
    {
        /// @brief Constructor
        /// @param[in] splines A pointer to splines
        /// @param[in] splineIndex The index of the current spline
        /// @param[in] point The point from where the distance is calculated
        FuncDistanceFromAPoint(Splines* splines,
                               UInt splineIndex,
                               Point point) : m_spline(splines),
                                              m_splineIndex(splineIndex),
                                              m_point(point)
        {
            if (m_spline == nullptr)
            {
                throw std::invalid_argument("FuncDistanceFromAPoint::m_spline is nullptr");
            }
        };

        /// @brief This is the function we want to find the root of
        double operator()(double adimensionalDistanceReferencePoint) const
        {
            const auto pointOnSpline = ComputePointOnSplineAtAdimensionalDistance(m_spline->m_splineNodes[m_splineIndex], m_spline->m_splineDerivatives[m_splineIndex], adimensionalDistanceReferencePoint);
            return ComputeDistance(m_point, pointOnSpline, Projection::cartesian);
        }

        Splines* m_spline;                  ///< Pointer to splines
        UInt m_splineIndex;                 ///< Spline index
        Point m_point;                      ///< The point from where the distance is calculated
        double m_DimensionalDistance = 0.0; ///< Dimensional distance
    };

} // namespace meshkernel
