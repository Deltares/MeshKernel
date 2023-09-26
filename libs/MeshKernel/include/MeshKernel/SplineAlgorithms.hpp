//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2023.
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

#include <memory>

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/LandBoundary.hpp"
#include "MeshKernel/Utilities/LinearAlgebra.hpp"

namespace meshkernel
{

    /// @brief Provide algorithms operating on splines.
    class SplineAlgorithms
    {
    public:
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

        /// @brief Computes curvature in a spline point (comp_curv)
        /// @param[in] splinePoints the spline points
        /// @param[in] splineDerivative the spline derivative points
        /// @param[in] adimensionalPointCoordinate The adimensional coordinate of the point along the spline
        /// @param[in] projection The coordinate system projection
        /// @returns The computed curvatureFactor, normal vector and tangential vector
        static std::tuple<Point, Point, double>
        ComputeCurvatureOnSplinePoint(const std::vector<Point>& splinePoints,
                                      const std::vector<Point>& splineDerivative,
                                      double adimensionalPointCoordinate,
                                      const Projection::Type projection);

        /// @brief Evaluate a spline function (splint)
        ///
        /// @param [in] coordinates      The spline points
        /// @param [in] secondDerivative Second derivative of the spline
        /// @param [in] evaluationPoint  The point at which the spline is to be evaluated
        /// @returns the result of evaluating the spline at the point.
        static Point Evaluate(const std::vector<Point>& coordinates, const std::vector<Point>& secondDerivative, const double evaluationPoint);

        /// @brief Snap the spline to the land boundary (snap_spline)
        ///
        /// Uses an iterated Lagrange multiplier scheme to snap the points.
        /// @param[in, out] splinePoints The spline points, these will be updated
        /// @param[in] splineDerivative Derivative of the spline at the spline points
        /// @param[in] landBoundary The land boundary to which the spline will be snapped
        /// @param[in] projection The coordinate system projection
        /// @param[in] numberOfIterations The maximum number of iterations to be performed
        static void SnapSplineToBoundary(std::vector<Point>& splinePoints,
                                         const std::vector<Point>& splineDerivative,
                                         const LandBoundary& landBoundary,
                                         const Projection::Type projection,
                                         const int numberOfIterations = constants::numeric::defaultSnappingIterations);

    private:
        /// @brief Compute the spline weights (from snap_spline)
        ///
        /// @param [in] xf The x-coordinates
        /// @param [in] yf y-coordinates
        /// @param [in] projection
        /// @returns vector of the spline weights
        static lin_alg::ColVector<double> ComputeSplineWeights(const lin_alg::ColVector<double>& xf,
                                                               const lin_alg::ColVector<double>& yf,
                                                               const Projection::Type projection);
        /// @brief Compute the spline sample points.
        ///
        /// @param [in] splinePoints The spline points
        /// @param [in] aMatrix Matrix expressing the approximation of the spline evaluated at the sample points.
        /// @returns the matrix product with the spline points.
        static std::tuple<lin_alg::ColVector<double>, lin_alg::ColVector<double>> ComputeSamplePoints(const std::vector<Point>& splinePoints,
                                                                                                      const lin_alg::Matrix<double, Eigen::ColMajor>& aMatrix);

        /// @brief Sample the spline. (sample_spline)
        ///
        /// Sample points can be at a high resolution that the control points.
        /// @param [in] splinePoints Spline to be sampled
        /// @param [in] intermediatePointCount Number of sample points
        /// @param [out] samplePoints Result of sample the spline
        static void SampleSpline(const std::vector<Point>& splinePoints,
                                 const size_t intermediatePointCount,
                                 std::vector<Point>& samplePoints);

        /// @brief Compute the interpolation matrix (comp_afinespline)
        ///
        /// @param [in] numberOfSplinePoints
        /// @param [in] intervalRefinement
        /// @param [out] numberOfSamplePoints Number of sample points
        /// @param [out] interpolationMatrix Matrix containing the spline coefficients at the sample points.
        static void ComputeInterpolationMatrix(const Eigen::Index numberOfSplinePoints,
                                               const Eigen::Index intervalRefinement,
                                               Eigen::Index& numberOfSamplePoints,
                                               lin_alg::Matrix<double, Eigen::ColMajor>& interpolationMatrix);

        /// @brief Compute the inverse of the least squares matrix (from snap_spline)
        ///
        /// Computes the inverse of \f$a^t w a\f$, where \f$a\f$ is the matrix of spline coefficients, and \f$w\f$ the spline weights.
        /// The vector \f$w\f$, should be considered as a matrix with the weight values along the main diagonal.
        /// @param [in] splineCoefficients The spline coefficients
        /// @param [in] weights The slpine weights
        static lin_alg::Matrix<double, Eigen::ColMajor> ComputeLeastSquaresMatrixInverse(const lin_alg::Matrix<double, Eigen::ColMajor>& splineCoefficients,
                                                                                         const lin_alg::ColVector<double>& weights);
    };

} // namespace meshkernel
