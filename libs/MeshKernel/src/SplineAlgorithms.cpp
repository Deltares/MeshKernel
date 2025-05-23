//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2025.
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

#include "MeshKernel/SplineAlgorithms.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Operations.hpp"

#include <Eigen/Core>
#include <Eigen/LU>

using meshkernel::SplineAlgorithms;

std::vector<meshkernel::Point> SplineAlgorithms::SecondOrderDerivative(const std::vector<Point>& spline, size_t startIndex, size_t endIndex)
{
    const auto numNodes = endIndex - startIndex + 1;
    std::vector<Point> u(numNodes, {0.0, 0.0});
    std::vector<Point> coordinatesDerivative(numNodes, {0.0, 0.0});

    UInt index = 1;
    for (auto i = startIndex + 1; i < endIndex; i++)
    {
        const Point p = coordinatesDerivative[index - 1] * 0.5 + 2.0;
        coordinatesDerivative[index].x = -0.5 / p.x;
        coordinatesDerivative[index].y = -0.5 / p.y;

        const Point delta = spline[i + 1] - spline[i] - (spline[i] - spline[i - 1]);
        u[index] = (delta * 6.0 / 2.0 - u[index - 1] * 0.5) / p;
        index++;
    }

    coordinatesDerivative.back() = {0.0, 0.0};
    for (int i = static_cast<int>(coordinatesDerivative.size()) - 2; i >= 0; --i)
    {
        coordinatesDerivative[i] = coordinatesDerivative[i] * coordinatesDerivative[i + 1] + u[i];
    }

    return coordinatesDerivative;
}

std::vector<double> SplineAlgorithms::SecondOrderDerivative(const std::vector<double>& coordinates, size_t startIndex, size_t endIndex)
{
    const auto numNodes = endIndex - startIndex + 1;
    std::vector<double> u(numNodes, 0.0);
    std::vector<double> coordinatesDerivatives(numNodes, 0.0);

    UInt index = 1;
    for (auto i = startIndex + 1; i < endIndex; i++)
    {
        const double p = coordinatesDerivatives[index - 1] * 0.5 + 2.0;
        coordinatesDerivatives[index] = -0.5 / p;

        const double delta = coordinates[i + 1] - coordinates[i] - (coordinates[i] - coordinates[i - 1]);
        u[index] = (delta * 6.0 / 2.0 - u[index - 1] * 0.5) / p;
        index++;
    }

    coordinatesDerivatives.back() = 0.0;
    for (int i = static_cast<int>(coordinatesDerivatives.size()) - 2; i >= 0; --i)
    {
        coordinatesDerivatives[i] = coordinatesDerivatives[i] * coordinatesDerivatives[i + 1] + u[i];
    }

    return coordinatesDerivatives;
}

std::tuple<meshkernel::Point, meshkernel::Point, double>
SplineAlgorithms::ComputeCurvatureOnSplinePoint(const std::vector<Point>& splinePoints,
                                                const std::vector<Point>& splineDerivative,
                                                double adimensionalPointCoordinate,
                                                const Projection projection)
{
    if (splinePoints.empty())
    {
        return {};
    }

    const auto numNodesFirstSpline = static_cast<UInt>(splinePoints.size());
    auto const leftCornerPoint = std::max(std::min(static_cast<UInt>(std::floor(adimensionalPointCoordinate)), numNodesFirstSpline - 2), 0U);
    auto const rightCornerPoint = std::max(leftCornerPoint + 1, 0U);

    if (rightCornerPoint >= numNodesFirstSpline)
    {
        throw ConstraintError("Coordinate out of bounds, resulting in index out of bounds: coordinate = {}, end index = {}, size = {}",
                              adimensionalPointCoordinate,
                              rightCornerPoint,
                              numNodesFirstSpline);
    }

    const auto leftSegment = static_cast<double>(rightCornerPoint) - adimensionalPointCoordinate;
    const auto rightSegment = adimensionalPointCoordinate - static_cast<double>(leftCornerPoint);

    const auto pointCoordinate = ComputePointOnSplineAtAdimensionalDistance(splinePoints, splineDerivative, adimensionalPointCoordinate);
    if (!pointCoordinate.IsValid())
    {
        throw AlgorithmError("Could not interpolate spline points.");
    }

    Point p = splinePoints[rightCornerPoint] - splinePoints[leftCornerPoint] +
              (splineDerivative[leftCornerPoint] * (-3.0 * leftSegment * leftSegment + 1.0) +
               splineDerivative[rightCornerPoint] * (3.0 * rightSegment * rightSegment - 1.0)) /
                  6.0;

    Point pp = splineDerivative[leftCornerPoint] * leftSegment +
               splineDerivative[rightCornerPoint] * rightSegment;

    if (projection == Projection::spherical)
    {
        p.TransformSphericalToCartesian(pointCoordinate.y);
        pp.TransformSphericalToCartesian(pointCoordinate.y);
    }

    double curvatureFactor = std::abs(pp.x * p.y - pp.y * p.x) / std::pow((p.x * p.x + p.y * p.y + 1e-8), 1.5);

    const auto incrementedPointCoordinate = pointCoordinate + p * 1e-4;
    const auto normalVector = NormalVectorOutside(pointCoordinate, incrementedPointCoordinate, projection);

    const auto distance = ComputeDistance(pointCoordinate, incrementedPointCoordinate, projection);
    const auto dx = GetDx(pointCoordinate, incrementedPointCoordinate, projection);
    const auto dy = GetDy(pointCoordinate, incrementedPointCoordinate, projection);

    Point tangentialVector;
    tangentialVector.x = dx / distance;
    tangentialVector.y = dy / distance;

    return {normalVector, tangentialVector, curvatureFactor};
}

meshkernel::Point SplineAlgorithms::Evaluate(const std::vector<Point>& splinePoints, const std::vector<Point>& secondDerivative, const double evaluationPoint)
{
    // Constant used in dflowfm code: splint.f90
    constexpr double splfac = 1.0;
    constexpr double eps = 1.0e-5;

    // Start and end index of interval containing evaluationPoint
    size_t startIndex = static_cast<size_t>(evaluationPoint);

    Point result;

    if (std::abs(evaluationPoint - std::floor(evaluationPoint)) <= eps)
    {
        result = splinePoints[startIndex];
    }
    else
    {
        double a = static_cast<double>(startIndex + 1) - evaluationPoint;
        double b = evaluationPoint - static_cast<double>(startIndex);

        result = a * splinePoints[startIndex] + b * splinePoints[startIndex + 1];
        result += (splfac / 6.0) * ((a * a * a - a) * secondDerivative[startIndex] + (b * b * b - b) * secondDerivative[startIndex + 1]);
    }

    return result;
}

lin_alg::ColVector<double> SplineAlgorithms::ComputeSplineWeights(const lin_alg::ColVector<double>& xf,
                                                                  const lin_alg::ColVector<double>& yf,
                                                                  const Projection projection)
{

    const Eigen::Index totalNumberOfPoints = xf.size();

    lin_alg::ColVector<double> weights(totalNumberOfPoints);

    // Compute weights
    for (Eigen::Index i = 1; i <= totalNumberOfPoints; ++i)
    {
        Eigen::Index il = std::max<Eigen::Index>(1, i - 1);
        Eigen::Index ir = std::min<Eigen::Index>(totalNumberOfPoints, i + 1);
        Point p1{xf[il - 1], yf[il - 1]};
        Point p2{xf[ir - 1], yf[ir - 1]};
        weights(i - 1) = 1.0 / std::sqrt(ComputeDistance(p1, p2, projection) / static_cast<double>(ir - il));
    }

    return weights;
}

std::tuple<lin_alg::ColVector<double>, lin_alg::ColVector<double>> SplineAlgorithms::ComputeSamplePoints(const std::vector<Point>& splinePoints,
                                                                                                         const lin_alg::Matrix<double, Eigen::ColMajor>& aMatrix)
{

    lin_alg::ColVector<double> xf(aMatrix.rows());
    lin_alg::ColVector<double> yf(aMatrix.rows());

    for (int r = 0; r < aMatrix.rows(); ++r)
    {
        xf(r) = 0.0;
        yf(r) = 0.0;

        for (int c = 0; c < aMatrix.cols(); ++c)
        {
            xf(r) += aMatrix(r, c) * splinePoints[c].x;
            yf(r) += aMatrix(r, c) * splinePoints[c].y;
        }
    }

    return {xf, yf};
}

void SplineAlgorithms::SampleSpline(const std::vector<Point>& splinePoints,
                                    const size_t intermediatePointCount,
                                    std::vector<Point>& samplePoints)
{

    if (splinePoints.empty())
    {
        throw ConstraintError("Spline is empty");
    }

    const size_t sampleCount = splinePoints.size() + (splinePoints.size() - 1) * intermediatePointCount;

    samplePoints.resize(sampleCount);

    const auto secondDerivative = SecondOrderDerivative(splinePoints, 0, splinePoints.size() - 1);

    const double intermediatePointCountFloat = static_cast<double>(intermediatePointCount + 1);
    double evaluationPoint;

    size_t count = 0;

    for (size_t i = 0; i < splinePoints.size() - 1; ++i)
    {
        const double floatI = static_cast<double>(i);

        for (size_t j = 0; j <= intermediatePointCount; ++j)
        {
            evaluationPoint = floatI + static_cast<double>(j) / intermediatePointCountFloat;
            samplePoints[count] = Evaluate(splinePoints, secondDerivative, evaluationPoint);
            ++count;
        }
    }

    evaluationPoint = static_cast<double>(splinePoints.size() - 1);
    samplePoints[count] = Evaluate(splinePoints, secondDerivative, evaluationPoint);
}

void SplineAlgorithms::ComputeInterpolationMatrix(const Eigen::Index numberOfSplinePoints,
                                                  const Eigen::Index intervalRefinement,
                                                  Eigen::Index& numberOfSamplePoints,
                                                  lin_alg::Matrix<double, Eigen::ColMajor>& interpolationMatrix)
{

    if (numberOfSplinePoints < 1)
    {
        throw ConstraintError("Invalid spline point count: {}", numberOfSplinePoints);
    }

    numberOfSamplePoints = numberOfSplinePoints + (numberOfSplinePoints - 1) * intervalRefinement;
    interpolationMatrix.resize(numberOfSamplePoints, numberOfSplinePoints);

    Point p{0.0, 0.0};
    std::vector<Point> locations(static_cast<size_t>(numberOfSplinePoints), p);
    std::vector<Point> xf;

    for (size_t i = 0; i < locations.size(); ++i)
    {
        locations[i].x = 1.0;
        SampleSpline(locations, intervalRefinement, xf);

        // Assign values to column of matrix.
        for (UInt j = 0; j < xf.size(); ++j)
        {
            interpolationMatrix(j, i) = xf[j].x;
        }

        locations[i].x = 0.0;
    }
}

lin_alg::Matrix<double, Eigen::ColMajor> SplineAlgorithms::ComputeLeastSquaresMatrixInverse(const lin_alg::Matrix<double, Eigen::ColMajor>& splineCoefficients,
                                                                                            const lin_alg::ColVector<double>& weights)
{

    lin_alg::Matrix<double, Eigen::ColMajor> atwa(splineCoefficients.cols(), splineCoefficients.cols());

    // Compute matrix A^t W A
    // least squares
    for (int i = 0; i < splineCoefficients.cols(); ++i)
    {
        for (int j = 0; j < splineCoefficients.cols(); ++j)
        {
            atwa(i, j) = 0.0;

            for (int k = 0; k < splineCoefficients.rows(); ++k)
            {
                atwa(i, j) += splineCoefficients(k, i) * weights(k) * splineCoefficients(k, j);
            }
        }
    }

    return atwa.inverse();
}

void SplineAlgorithms::SnapSplineToBoundary(std::vector<Point>& splinePoints,
                                            const std::vector<Point>& splineDerivative,
                                            const LandBoundary& landBoundary,
                                            const Projection projection,
                                            const int numberOfIterations)
{

    if (splinePoints.empty())
    {
        throw ConstraintError("Empty spline");
    }

    if (splineDerivative.empty())
    {
        throw ConstraintError("Empty spline derivative");
    }

    if (splinePoints.size() != splineDerivative.size())
    {
        throw ConstraintError("Spline and derivative are not the same size: {} /= {}",
                              splinePoints.size(),
                              splineDerivative.size());
    }

    constexpr double tolerance = 1.0e-5;

    // Number of additional points between spline control points for sampled spline.
    constexpr Eigen::Index numberRefinements = 19;
    constexpr Eigen::Index numberOfConstraints = 2;

    // Save original spline end points
    Point startPoint = splinePoints.front();
    Point endPoint = splinePoints.back();
    Eigen::Index numberOfSamplePoints = 0;

    lin_alg::Matrix<double, Eigen::ColMajor> interpolationMatrix;

    // Compute spline coefficients.
    // Resizes the matrix to be numberOfSamplePoints by splinePoints.size()
    ComputeInterpolationMatrix(splinePoints.size(),
                               numberRefinements,
                               numberOfSamplePoints,
                               interpolationMatrix);

    // Returns two Eigen vectors, for the x- and y-sample points.
    auto [splineValuesX, splineValuesY] = ComputeSamplePoints(splinePoints, interpolationMatrix);

    lin_alg::ColVector<double> weights(ComputeSplineWeights(splineValuesX, splineValuesY, projection));

    // The tangent vector is unused
    auto [startNormal, startTangent, startCurvature] = ComputeCurvatureOnSplinePoint(splinePoints, splineDerivative, 0.0, projection);
    auto [endNormal, endTangent, endCurvature] = ComputeCurvatureOnSplinePoint(splinePoints, splineDerivative, static_cast<double>(splinePoints.size() - 1), projection);

    // (a^t w a)^-1
    lin_alg::Matrix<double, Eigen::ColMajor> leastSquaresMatrixInverse(ComputeLeastSquaresMatrixInverse(interpolationMatrix, weights));

    // Matrix and vectors used in the Lagrange multiplier method.
    lin_alg::Matrix<double, Eigen::ColMajor> bMatrix(numberOfConstraints, static_cast<Eigen::Index>(splinePoints.size()));
    lin_alg::Matrix<double, Eigen::ColMajor> cMatrix(numberOfConstraints, static_cast<Eigen::Index>(splinePoints.size()));
    lin_alg::ColVector<double> dVector(numberOfConstraints);
    // Lagrange multiplier values
    lin_alg::ColVector<double> constraintValues(numberOfConstraints);

    bMatrix.setZero();
    cMatrix.setZero();

    // Compute the matrix for the constraints.
    bMatrix(0, 0) = startNormal.y;
    cMatrix(0, 0) = -startNormal.x;
    dVector(0) = startNormal.y * startPoint.x - startNormal.x * startPoint.y;
    bMatrix(1, splinePoints.size() - 1) = endNormal.y;
    cMatrix(1, splinePoints.size() - 1) = -endNormal.x;
    dVector(1) = endNormal.y * endPoint.x - endNormal.x * endPoint.y;

    // eMatrix = bMatrix * leastSquaresMatrixInverse * bMatrix.transpose() + cMatrix * leastSquaresMatrixInverse * cMatrix.transpose()

    lin_alg::Matrix<double, Eigen::ColMajor> temp1 = leastSquaresMatrixInverse * bMatrix.transpose();
    lin_alg::Matrix<double, Eigen::ColMajor> temp2 = leastSquaresMatrixInverse * cMatrix.transpose();

    lin_alg::Matrix<double, Eigen::ColMajor> eMatrix = bMatrix * temp1 + cMatrix * temp2;

    constraintValues.setZero();
    // Inplace inversion of the e-matrix.
    eMatrix = eMatrix.inverse();

    std::vector<Point> nearestPoints(numberOfSamplePoints, Point());

    // Weighted evaluation of spline at sample points
    lin_alg::ColVector<double> atwxb(splinePoints.size());
    lin_alg::ColVector<double> atwyb(splinePoints.size());

    // Intermediate result vector when solving the
    lin_alg::ColVector<double> rhsx(splinePoints.size());
    lin_alg::ColVector<double> rhsy(splinePoints.size());

    lin_alg::ColVector<double> xVals(splinePoints.size());
    lin_alg::ColVector<double> yVals(splinePoints.size());

    lin_alg::ColVector<double> xValsOld(splinePoints.size());
    lin_alg::ColVector<double> yValsOld(splinePoints.size());

    int iterationCount = 1;

    for (size_t i = 0; i < splinePoints.size(); ++i)
    {
        xVals(i) = splinePoints[i].x;
        yVals(i) = splinePoints[i].y;
    }

    bool converged = false;

    while (!converged && iterationCount <= numberOfIterations)
    {

        for (int i = 0; i < numberOfSamplePoints; ++i)
        {
            Point point(splineValuesX(i), splineValuesY(i));

            double smallestDistance;
            double scaledDistance;
            UInt segmentIndex;

            landBoundary.FindNearestPoint(point, projection, nearestPoints[i], smallestDistance, segmentIndex, scaledDistance);
        }

        for (int i = 0; i < interpolationMatrix.cols(); ++i)
        {
            atwxb(i) = 0.0;
            atwyb(i) = 0.0;

            for (int j = 0; j < interpolationMatrix.rows(); ++j)
            {
                atwxb(i) += interpolationMatrix(j, i) * weights(j) * nearestPoints[j].x;
                atwyb(i) += interpolationMatrix(j, i) * weights(j) * nearestPoints[j].y;
            }
        }

        // Compute constraints.
        constraintValues = eMatrix * (bMatrix * (leastSquaresMatrixInverse * atwxb) + cMatrix * (leastSquaresMatrixInverse * atwyb) - dVector);

        rhsx = atwxb - bMatrix.transpose() * constraintValues;
        rhsy = atwyb - cMatrix.transpose() * constraintValues;

        xValsOld = xVals;
        yValsOld = yVals;

        xVals = leastSquaresMatrixInverse * rhsx;
        yVals = leastSquaresMatrixInverse * rhsy;

        splineValuesX = interpolationMatrix * xVals;
        splineValuesY = interpolationMatrix * yVals;

        // is there a better check for convergence?
        // AND what should the tolerance be?
        converged = (xVals - xValsOld).norm() + (yVals - yValsOld).norm() < tolerance;
        ++iterationCount;
    }

    // Copy vectors back to array of points.
    for (size_t i = 0; i < splinePoints.size(); ++i)
    {
        splinePoints[i].x = xVals(i);
        splinePoints[i].y = yVals(i);
    }
}
