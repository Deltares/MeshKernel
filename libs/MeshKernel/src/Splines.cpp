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

#include <Eigen/Core>
#include <Eigen/LU>

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/LandBoundaries.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Splines.hpp>

using meshkernel::Splines;

Splines::Splines(Projection projection) : m_projection(projection) {}

Splines::Splines(CurvilinearGrid const& grid)

{
    // first the m_n m_m-gridlines
    std::vector<std::vector<Point>> mGridLines(grid.m_numN, std::vector<Point>(grid.m_numM));
    for (UInt n = 0; n < grid.m_numN; ++n)
    {
        for (UInt m = 0; m < grid.m_numM; ++m)
        {
            mGridLines[n][m] = grid.m_gridNodes[m][n];
        }
        AddSpline(mGridLines[n], 0, static_cast<UInt>(mGridLines[n].size()));
    }

    // then the m_m m_n-gridlines
    std::vector<std::vector<Point>> nGridLines(grid.m_numM, std::vector<Point>(grid.m_numN));
    for (UInt m = 0; m < grid.m_numM; ++m)
    {
        AddSpline(grid.m_gridNodes[m], 0, static_cast<UInt>(grid.m_gridNodes[m].size()));
    }

    m_projection = grid.m_projection;
}

/// add a new spline, return the index
void Splines::AddSpline(const std::vector<Point>& splines, UInt start, UInt size)
{
    if (size == 0)
    {
        return;
    }

    // copy the spline nodes from start to start + size
    UInt count = 0;
    std::vector<Point> splinesNodes(size);
    for (auto i = start; i < start + size; ++i)
    {
        splinesNodes[count] = splines[i];
        count++;
    }
    m_splineNodes.emplace_back(splinesNodes);

    // compute second order derivatives
    std::vector<Point> splineDerivatives(splinesNodes.size());
    const auto indices = FindIndices(splinesNodes, 0, static_cast<UInt>(splinesNodes.size()), constants::missing::doubleValue);
    for (auto index : indices)
    {
        const auto& [startIndex, endIndex] = index;
        const auto derivatives = SecondOrderDerivative(splinesNodes, startIndex, endIndex);
        for (auto j = startIndex; j <= endIndex; ++j)
        {
            splineDerivatives[j] = derivatives[j - startIndex];
        }
    }
    m_splineDerivatives.emplace_back(splineDerivatives);

    m_splinesLength.emplace_back(ComputeSplineLength(GetNumSplines() - 1, 0.0, static_cast<double>(size - 1)));
}

void Splines::DeleteSpline(UInt splineIndex)
{
    m_splineNodes.erase(m_splineNodes.begin() + splineIndex);
    m_splineDerivatives.erase(m_splineDerivatives.begin() + splineIndex);
    m_splinesLength.erase(m_splinesLength.begin() + splineIndex);
}

void Splines::AddPointInExistingSpline(UInt splineIndex, const Point& point)
{
    if (splineIndex > GetNumSplines())
    {
        throw std::invalid_argument("Splines::AddPointInExistingSpline: Invalid spline index.");
    }
    m_splineNodes[splineIndex].emplace_back(point);
}

bool Splines::GetSplinesIntersection(UInt first,
                                     UInt second,
                                     double& crossProductIntersection,
                                     Point& intersectionPoint,
                                     double& firstSplineRatio,
                                     double& secondSplineRatio)
{
    double minimumCrossingDistance = std::numeric_limits<double>::max();
    double crossingDistance;
    UInt numCrossing = 0;
    double firstCrossingRatio = -1.0;
    double secondCrossingRatio = -1.0;
    UInt firstCrossingIndex = 0;
    UInt secondCrossingIndex = 0;
    Point closestIntersection;
    const auto numNodesFirstSpline = static_cast<UInt>(m_splineNodes[first].size());
    const auto numNodesSecondSpline = static_cast<UInt>(m_splineNodes[second].size());

    // First find a valid crossing, the closest to spline central point
    for (UInt n = 0; n < numNodesFirstSpline - 1; n++)
    {
        for (UInt nn = 0; nn < numNodesSecondSpline - 1; nn++)
        {
            Point intersection;
            double crossProduct;
            double firstRatio;
            double secondRatio;
            const bool areCrossing = AreSegmentsCrossing(m_splineNodes[first][n],
                                                         m_splineNodes[first][n + 1],
                                                         m_splineNodes[second][nn],
                                                         m_splineNodes[second][nn + 1],
                                                         false,
                                                         m_projection,
                                                         intersection,
                                                         crossProduct,
                                                         firstRatio,
                                                         secondRatio);

            if (areCrossing)
            {
                if (numNodesFirstSpline == 2)
                {
                    crossingDistance = std::min(minimumCrossingDistance, std::abs(firstRatio - 0.5));
                }
                else if (numNodesSecondSpline == 2)
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
                    firstCrossingIndex = n;            // TI0
                    secondCrossingIndex = nn;          // TJ0
                    firstCrossingRatio = firstRatio;   // SL
                    secondCrossingRatio = secondRatio; // SM
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

    double firstCrossing = IsEqual(firstCrossingRatio, -1.0) ? 0.0 : static_cast<double>(firstCrossingIndex) + firstCrossingRatio;
    double secondCrossing = IsEqual(secondCrossingRatio, -1.0) ? 0.0 : static_cast<double>(secondCrossingIndex) + secondCrossingRatio;

    // use bisection to find the intersection
    double squaredDistanceBetweenCrossings = std::numeric_limits<double>::max();
    const double maxSquaredDistanceBetweenCrossings = 1e-12;
    const double maxDistanceBetweenNodes = 0.0001;
    double firstRatioIterations = 1.0;
    double secondRatioIterations = 1.0;
    UInt numIterations = 0;
    while (squaredDistanceBetweenCrossings > maxSquaredDistanceBetweenCrossings && numIterations < 20)
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

        firstCrossing = std::max(0.0, std::min(firstCrossing, double(numNodesFirstSpline)));
        secondCrossing = std::max(0.0, std::min(secondCrossing, double(numNodesSecondSpline)));

        const double firstLeft = std::max(0.0, std::min(double(numNodesFirstSpline - 1), firstCrossing - firstRatioIterations / 2.0));
        const double firstRight = std::max(0.0, std::min(double(numNodesFirstSpline - 1), firstCrossing + firstRatioIterations / 2.0));

        const double secondLeft = std::max(0.0, std::min(double(numNodesSecondSpline - 1), secondCrossing - secondRatioIterations / 2.0));
        const double secondRight = std::max(0.0, std::min(double(numNodesSecondSpline - 1), secondCrossing + secondRatioIterations / 2.0));

        firstRatioIterations = firstRight - firstLeft;
        secondRatioIterations = secondRight - secondLeft;

        const auto firstLeftSplinePoint = ComputePointOnSplineAtAdimensionalDistance(m_splineNodes[first], m_splineDerivatives[first], firstLeft);
        if (!firstLeftSplinePoint.IsValid())
        {
            return false;
        }

        const auto firstRightSplinePoint = ComputePointOnSplineAtAdimensionalDistance(m_splineNodes[first], m_splineDerivatives[first], firstRight);
        if (!firstRightSplinePoint.IsValid())
        {
            return false;
        }

        const auto secondLeftSplinePoint = ComputePointOnSplineAtAdimensionalDistance(m_splineNodes[second], m_splineDerivatives[second], secondLeft);
        if (!secondLeftSplinePoint.IsValid())
        {
            return false;
        }

        const auto secondRightSplinePoint = ComputePointOnSplineAtAdimensionalDistance(m_splineNodes[second], m_splineDerivatives[second], secondRight);
        if (!secondRightSplinePoint.IsValid())
        {
            return false;
        }

        Point oldIntersection = closestIntersection;

        double crossProduct;
        double firstRatio = constants::missing::doubleValue;
        double secondRatio = constants::missing::doubleValue;
        const bool areCrossing = AreSegmentsCrossing(firstLeftSplinePoint,
                                                     firstRightSplinePoint,
                                                     secondLeftSplinePoint,
                                                     secondRightSplinePoint,
                                                     true,
                                                     m_projection,
                                                     closestIntersection,
                                                     crossProduct,
                                                     firstRatio,
                                                     secondRatio);

        // search close by
        if (firstRatio > -2.0 && firstRatio < 3.0 && secondRatio > -2.0 && secondRatio < 3.0)
        {
            const double previousFirstCrossing = firstCrossing;
            const double previousSecondCrossing = secondCrossing;

            firstCrossing = firstLeft + firstRatio * (firstRight - firstLeft);
            secondCrossing = secondLeft + secondRatio * (secondRight - secondLeft);

            firstCrossing = std::max(0.0, std::min(static_cast<double>(numNodesFirstSpline) - 1.0, firstCrossing));
            secondCrossing = std::max(0.0, std::min(static_cast<double>(numNodesSecondSpline) - 1.0, secondCrossing));

            if (areCrossing)
            {
                numCrossing = 1;
                crossProductIntersection = crossProduct;
            }

            if (std::abs(firstCrossing - previousFirstCrossing) > maxDistanceBetweenNodes ||
                std::abs(secondCrossing - previousSecondCrossing) > maxDistanceBetweenNodes)
            {
                squaredDistanceBetweenCrossings = ComputeSquaredDistance(oldIntersection, closestIntersection, m_projection);
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

    // not crossing
    return false;
}

double Splines::ComputeSplineLength(UInt index,
                                    double startAdimensionalCoordinate,
                                    double endAdimensionalCoordinate,
                                    UInt numSamples,
                                    bool accountForCurvature,
                                    double height,
                                    double assignedDelta) const
{
    if (m_splineNodes[index].empty())
    {
        return 0.0;
    }

    double delta = assignedDelta;
    UInt numPoints = static_cast<UInt>(endAdimensionalCoordinate / delta) + 1;
    if (delta < 0.0)
    {
        delta = 1.0 / static_cast<double>(numSamples);
        // TODO: Refactor or at least document the calculation of "numPoints"
        numPoints = static_cast<UInt>(std::max(std::floor(0.9999 + (endAdimensionalCoordinate - startAdimensionalCoordinate) / delta), 10.0));
        delta = (endAdimensionalCoordinate - startAdimensionalCoordinate) / static_cast<double>(numPoints);
    }

    // first point
    auto leftPoint = ComputePointOnSplineAtAdimensionalDistance(m_splineNodes[index], m_splineDerivatives[index], startAdimensionalCoordinate);
    double splineLength = 0.0;

    auto rightPointCoordinateOnSpline = startAdimensionalCoordinate;
    for (UInt p = 0; p < numPoints; ++p)
    {
        const double leftPointCoordinateOnSpline = rightPointCoordinateOnSpline;
        rightPointCoordinateOnSpline += delta;
        if (rightPointCoordinateOnSpline > endAdimensionalCoordinate)
        {
            rightPointCoordinateOnSpline = endAdimensionalCoordinate;
        }

        const auto rightPoint = ComputePointOnSplineAtAdimensionalDistance(m_splineNodes[index], m_splineDerivatives[index], rightPointCoordinateOnSpline);
        if (!rightPoint.IsValid())
        {
            continue;
        }

        double curvatureFactor = 0.0;
        if (accountForCurvature)
        {
            const auto [normalVector, tangentialVector, computedCurvatureFactor] = ComputeCurvatureOnSplinePoint(index, 0.5 * (rightPointCoordinateOnSpline + leftPointCoordinateOnSpline));
            curvatureFactor = computedCurvatureFactor;
        }
        splineLength = splineLength + ComputeDistance(leftPoint, rightPoint, m_projection) * (1.0 + curvatureFactor * height);
        leftPoint = rightPoint;
    }

    return splineLength;
}

std::tuple<meshkernel::Point, meshkernel::Point, double>
Splines::ComputeCurvatureOnSplinePoint(UInt splineIndex, double adimensionalPointCoordinate) const
{
    if (m_splineNodes[splineIndex].empty())
    {
        return {};
    }

    const auto numNodesFirstSpline = static_cast<UInt>(m_splineNodes[splineIndex].size());
    auto const leftCornerPoint = std::max(std::min(static_cast<UInt>(std::floor(adimensionalPointCoordinate)), numNodesFirstSpline - 1), 0U);
    auto const rightCornerPoint = std::max(leftCornerPoint + 1, 0U);

    const auto leftSegment = static_cast<double>(rightCornerPoint) - adimensionalPointCoordinate;
    const auto rightSegment = adimensionalPointCoordinate - static_cast<double>(leftCornerPoint);

    const auto pointCoordinate = ComputePointOnSplineAtAdimensionalDistance(m_splineNodes[splineIndex], m_splineDerivatives[splineIndex], adimensionalPointCoordinate);
    if (!pointCoordinate.IsValid())
    {
        throw AlgorithmError("Splines::ComputeCurvatureOnSplinePoint: Could not interpolate spline points.");
    }

    Point p = m_splineNodes[splineIndex][rightCornerPoint] - m_splineNodes[splineIndex][leftCornerPoint] +
              (m_splineDerivatives[splineIndex][leftCornerPoint] * (-3.0 * leftSegment * leftSegment + 1.0) +
               m_splineDerivatives[splineIndex][rightCornerPoint] * (3.0 * rightSegment * rightSegment - 1.0)) /
                  6.0;

    Point pp = m_splineDerivatives[splineIndex][leftCornerPoint] * leftSegment +
               m_splineDerivatives[splineIndex][rightCornerPoint] * rightSegment;

    if (m_projection == Projection::spherical)
    {
        p.TransformSphericalToCartesian(pointCoordinate.y);
        pp.TransformSphericalToCartesian(pointCoordinate.y);
    }

    double curvatureFactor = std::abs(pp.x * p.y - pp.y * p.x) / std::pow((p.x * p.x + p.y * p.y + 1e-8), 1.5);

    const auto incrementedPointCoordinate = pointCoordinate + p * 1e-4;
    const auto normalVector = NormalVectorOutside(pointCoordinate, incrementedPointCoordinate, m_projection);

    const auto distance = ComputeDistance(pointCoordinate, incrementedPointCoordinate, m_projection);
    const auto dx = GetDx(pointCoordinate, incrementedPointCoordinate, m_projection);
    const auto dy = GetDy(pointCoordinate, incrementedPointCoordinate, m_projection);

    Point tangentialVector;
    tangentialVector.x = dx / distance;
    tangentialVector.y = dy / distance;

    return {normalVector, tangentialVector, curvatureFactor};
}

std::vector<meshkernel::Point> Splines::SecondOrderDerivative(const std::vector<Point>& spline, UInt startIndex, UInt endIndex)
{
    const auto numNodes = endIndex - startIndex + 1;
    std::vector<Point> u(numNodes, {0.0, 0.0});
    std::vector<Point> coordinatesDerivative(numNodes, {0.0, 0.0});

    UInt index = 1;
    for (auto i = startIndex + 1; i < numNodes - 1; i++)
    {
        const Point p = coordinatesDerivative[index - 1] * 0.5 + 2.0;
        coordinatesDerivative[index].x = -0.5 / p.x;
        coordinatesDerivative[index].y = -0.5 / p.y;

        const Point delta = spline[i + 1] - spline[i] - (spline[i] - spline[i - 1]);
        u[index] = (delta * 6.0 / 2.0 - u[i - 1] * 0.5) / p;
        index++;
    }
    // TODO: C++ 20 for(auto& i :  views::reverse(vec))
    coordinatesDerivative.back() = {0.0, 0.0};
    for (auto i = numNodes - 2; i < numNodes; --i)
    {
        coordinatesDerivative[i] = coordinatesDerivative[i] * coordinatesDerivative[i + 1] + u[i];
    }

    return coordinatesDerivative;
}

std::vector<double> Splines::SecondOrderDerivative(const std::vector<double>& coordinates, UInt startIndex, UInt endIndex)
{
    const auto numNodes = endIndex - startIndex + 1;
    std::vector<double> u(numNodes, 0.0);
    std::vector<double> coordinatesDerivatives(numNodes, 0.0);

    UInt index = 1;
    for (auto i = startIndex + 1; i < numNodes - 1; i++)
    {
        const double p = coordinatesDerivatives[index - 1] * 0.5 + 2.0;
        coordinatesDerivatives[index] = -0.5 / p;

        const double delta = coordinates[i + 1] - coordinates[i] - (coordinates[i] - coordinates[i - 1]);
        u[index] = (delta * 6.0 / 2.0 - u[i - 1] * 0.5) / p;
        index++;
    }

    // TODO: C++ 20 for(auto& i :  views::reverse(vec))
    coordinatesDerivatives.back() = 0.0;
    for (auto i = numNodes - 2; i < numNodes; --i)
    {
        coordinatesDerivatives[i] = coordinatesDerivatives[i] * coordinatesDerivatives[i + 1] + u[i];
    }

    return coordinatesDerivatives;
}

meshkernel::Point Splines::evaluate(const std::vector<Point>& coordinates, const std::vector<Point>& secondDerivative, const double evaluationPoint)
{
    // Start and end index of interval containing evaluationPoint
    size_t startIndex = static_cast<size_t>(evaluationPoint);

    constexpr double splfac = 1.0;

    double a = static_cast<double>(startIndex + 1) - evaluationPoint;
    double b = evaluationPoint - static_cast<double>(startIndex);

    // TODO add check to see if point lies close to coordinate point

    Point result = a * coordinates[startIndex] + b * coordinates[startIndex + 1];
    result += (splfac / 6.0) * ((a * a * a - a) * secondDerivative[startIndex] + (b * b * b - b) * secondDerivative[startIndex + 1]);

    return result;
}

void Splines::sampleSpline(const std::vector<Point>& splinePoints,
                           const size_t intermediatePointCount,
                           std::vector<Point>& samplePoints)
{

    if (splinePoints.size() == 0)
    {
        throw meshkernel::ConstraintError("Spline is empty");
    }

    size_t sampleCount = splinePoints.size() + (splinePoints.size() - 1) * intermediatePointCount;

    samplePoints.resize(sampleCount);

    std::vector<Point> evaluatedSpline(splinePoints.size());
    // TODO Can I use the derivatives that are saved in the splines object?
    std::vector<Point> secondDerivative = SecondOrderDerivative(splinePoints, 0, splinePoints.size() - 1);

    double intermediatePointCountFloat = static_cast<double>(intermediatePointCount + 1);
    double evaluationPoint;

    size_t count = 0;

    for (size_t i = 0; i < splinePoints.size() - 1; ++i)
    {
        double floatI = static_cast<double>(i);

        for (size_t j = 0; j <= intermediatePointCount; ++j)
        {
            evaluationPoint = floatI + static_cast<double>(j) / intermediatePointCountFloat;
            samplePoints[count] = evaluate(splinePoints, secondDerivative, evaluationPoint);
            ++count;
        }
    }

    evaluationPoint = static_cast<double>(splinePoints.size() - 1);
    samplePoints[count] = evaluate(splinePoints, secondDerivative, evaluationPoint);
}

void Splines::compAfinespline(const UInt n, const UInt numRef, UInt& count, lin_alg::MatrixColMajor<double>& mat)
{

    if (n < 1)
    {
        throw meshkernel::ConstraintError(VariadicErrorMessage("Invalid spline point count: {}", n));
    }

    count = n + (n - 1) * numRef;

    // if ( Nr_in.lt.Nr ) then
    //    ierror = 2
    //    goto 1234
    // end if

    mat.resize(count, n);

    Point p;
    p.x = 0.0;
    p.y = 0.0;
    std::vector<Point> loc(n, p);
    std::vector<Point> xf;

    for (UInt i = 0; i < n; ++i)
    {
        loc[i].x = 1.0;
        sampleSpline(loc, numRef, xf);

        // Assign values to column of matrix.
        for (UInt j = 0; j < xf.size(); ++j)
        {
            mat(j, i) = xf[j].x;
        }

        loc[i].x = 0.0;
    }
}

lin_alg::MatrixColMajor<double> Splines::ComputeInterpolationMatrix(const lin_alg::MatrixColMajor<double>& splineCoefficients,
                                                                    const lin_alg::ColumnVector<double>& weights)
{

    lin_alg::MatrixColMajor<double> atwa(splineCoefficients.cols(), splineCoefficients.cols());

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

lin_alg::ColumnVector<double> Splines::ComputeSplineWeights(const std::vector<Point>& splinePoints,
                                                            const UInt totalNumberOfPoints,
                                                            const Projection projection)
{

    lin_alg::ColumnVector<double> weights(totalNumberOfPoints);

    // Compute weights
    for (size_t i = 1; i <= totalNumberOfPoints; ++i)
    {
        size_t il = std::max<size_t>(1, i - 1);
        size_t ir = std::min<size_t>(totalNumberOfPoints, i + 1);
        weights(i - 1) = 1.0 / std::sqrt(ComputeDistance(splinePoints[il - 1], splinePoints[ir - 1], projection) / static_cast<double>(ir - il));
    }

    return weights;
}

void Splines::snapSpline(const LandBoundaries& landBoundary,
                         const size_t splineIndex)
{
    if (splineIndex > GetNumSplines())
    {
        throw meshkernel::ConstraintError(VariadicErrorMessage("Invalid spline index: {}", splineIndex));
    }

    if (m_splineNodes[splineIndex].empty())
    {
        throw meshkernel::ConstraintError(VariadicErrorMessage("Empty spline at index: {}", splineIndex));
    }

    constexpr int MaxNumberOfIterations = 10;
    constexpr double tolerance = 1.0e-5;

    // Number of additional points between spline control points for sampled spline.
    constexpr UInt numberRefinements = 19;
    constexpr UInt numberOfConstraints = 2;

    std::vector<Point>& splinePoints(m_splineNodes[splineIndex]);
    std::vector<Point>& splineDerivative(m_splineDerivatives[splineIndex]);

    // Save original spline end points
    Point startPoint = splinePoints.front();
    Point endPoint = splinePoints.back();
    UInt totalNumberOfPoints = 0;

    lin_alg::MatrixColMajor<double> aMatrix;

    // Compute spline coefficients.
    // Resizes the matrix to be totalNumberOfPoints by splinePoints.size()
    compAfinespline(splinePoints.size(), numberRefinements, totalNumberOfPoints, aMatrix);

    // Vectors containing x- and y-coordinates
    // Simplifies linear algebra operations to have then in Eigen::vector format.
    // TODO Only problem, they are used only once, to compute a gemv. Is there an easier way?
    lin_alg::ColumnVector<double> vecX(splinePoints.size());
    lin_alg::ColumnVector<double> vecY(splinePoints.size());

    // extract spline points to vectors.
    for (size_t i = 0; i < splinePoints.size(); ++i)
    {
        vecX(i) = splinePoints[i].x;
        vecY(i) = splinePoints[i].y;
    }

    lin_alg::ColumnVector<double> xf(aMatrix * vecX);
    lin_alg::ColumnVector<double> yf(aMatrix * vecY);
    lin_alg::ColumnVector<double> xfOld(xf);
    lin_alg::ColumnVector<double> yfOld(yf);
    lin_alg::ColumnVector<double> weights(ComputeSplineWeights(splinePoints, totalNumberOfPoints, m_projection));

    auto [startCurvature, startNormal, startTangent] = ComputeCurvatureOnSplinePoint(splineIndex, 0.0);
    auto [endCurvature, endNormal, endTangent] = ComputeCurvatureOnSplinePoint(splineIndex, static_cast<double>(splinePoints.size() - 1));

    lin_alg::MatrixColMajor<double> atwaInverse(ComputeInterpolationMatrix(aMatrix, weights));
    lin_alg::MatrixColMajor<double> bMatrix(numberOfConstraints, splinePoints.size());
    lin_alg::MatrixColMajor<double> cMatrix(numberOfConstraints, splinePoints.size());
    lin_alg::ColumnVector<double> dVector(numberOfConstraints);
    lin_alg::ColumnVector<double> lambda(numberOfConstraints);

    // Compute the matrix for the constraints.
    bMatrix.setZero();
    cMatrix.setZero();

    bMatrix(0, 0) = startNormal.y;
    cMatrix(0, 0) = -startNormal.x;
    dVector(0) = startNormal.y * startPoint.x - startNormal.x * startPoint.y;
    bMatrix(1, splinePoints.size() - 1) = endNormal.y;
    cMatrix(1, splinePoints.size() - 1) = -endNormal.x;
    dVector(1) = endNormal.y * endPoint.x - endNormal.x * endPoint.y;

    lin_alg::MatrixColMajor<double> eMatrix = bMatrix * atwaInverse * bMatrix.transpose() + cMatrix * atwaInverse * cMatrix.transpose();

    lambda.setZero();
    // Inplace inversion of the e-matrix.
    eMatrix = eMatrix.inverse();

    bool converged = true;

    lin_alg::ColumnVector<double> xbVec(totalNumberOfPoints);
    lin_alg::ColumnVector<double> ybVec(totalNumberOfPoints);

    lin_alg::ColumnVector<double> atwxb(splinePoints.size());
    lin_alg::ColumnVector<double> atwyb(splinePoints.size());

    lin_alg::ColumnVector<double> rhsx(numberOfConstraints);
    lin_alg::ColumnVector<double> rhsy(numberOfConstraints);

    lin_alg::ColumnVector<double> xVals(numberOfConstraints);
    lin_alg::ColumnVector<double> yVals(numberOfConstraints);

    int iterationCount = 0;

    while (!converged && iterationCount <= MaxNumberOfIterations)
    {
        Point nearestPoint;

        // These variable are unused in the rest of the algorithm
        double smallestDistance;
        double scaledDistance;
        UInt segmentIndex;

        for (size_t i = 0; i < splinePoints.size(); ++i)
        {
            Point point(xf(i), yf(i));
            landBoundary.nearestPointOnLandBoundary(point, 2, nearestPoint, smallestDistance, segmentIndex, scaledDistance);
            xbVec(i) = nearestPoint.x;
            ybVec(i) = nearestPoint.y;
        }

        // TODO how to get a pointwise multiply here?
        // xbVec = weights.colwise() * xbVec.colwise();
        // Can then replace the two nested loops below with a gemv.

        // Could be replaced with gemv.
        for (size_t i = 0; i < splinePoints.size(); ++i)
        {
            atwxb(i) = 0.0;
            atwyb(i) = 0.0;

            for (size_t j = 0; j < splinePoints.size(); ++j)
            {
                atwxb(i) += aMatrix(j, i) * weights(j) * xbVec(j);
                atwyb(i) += aMatrix(j, i) * weights(j) * ybVec(j);
            }
        }

        // Compute constraints.
        lambda = eMatrix * (bMatrix * atwaInverse * atwxb + cMatrix * atwaInverse * atwyb - dVector);

        rhsx = atwxb - bMatrix.transpose() * lambda;
        rhsy = atwyb - cMatrix.transpose() * lambda;

        xVals = atwaInverse * rhsx;
        yVals = atwaInverse * rhsy;

        xf = aMatrix * xVals;
        yf = aMatrix * yVals;

        // TODO is there a better check for convergence?
        // AND what should the tolerance be?
        converged = (xf - xfOld).norm() + (yf - yfOld).norm() < tolerance;
        xfOld = xf;
        yfOld = yf;
        ++iterationCount;
    }

    if (converged)
    {
        // Copy vectors back to array of points.
        for (size_t i = 0; i < splinePoints.size(); ++i)
        {
            splinePoints[i].x = xf(i);
            splinePoints[i].y = yf(i);
        }

        splineDerivative = SecondOrderDerivative(splinePoints, 0, splinePoints.size() - 1);
    }
    else
    {
        throw AlgorithmError(VariadicErrorMessage("Failed to reach convergence for snapping spline {} to land boundary after {} iterations for {} tolerance",
                                                  splineIndex, MaxNumberOfIterations, tolerance));
    }
}

std::tuple<std::vector<meshkernel::Point>, std::vector<double>>
Splines::ComputePointOnSplineFromAdimensionalDistance(UInt index,
                                                      double maximumGridHeight,
                                                      bool isSpacingCurvatureAdapted,
                                                      const std::vector<double>& distances)
{

    std::vector<Point> points(distances.size());
    std::vector<double> adimensionalDistances(distances.size());

    FuncAdimensionalToDimensionalDistanceOnSpline func(this, index, isSpacingCurvatureAdapted, maximumGridHeight);
    const auto numNodes = static_cast<UInt>(m_splineNodes[index].size());
    for (UInt i = 0, size = static_cast<UInt>(distances.size()); i < size; ++i)
    {
        func.SetDimensionalDistance(distances[i]);
        adimensionalDistances[i] = FindFunctionRootWithGoldenSectionSearch(func, 0, static_cast<double>(numNodes) - 1.0);
        points[i] = ComputePointOnSplineAtAdimensionalDistance(m_splineNodes[index], m_splineDerivatives[index], adimensionalDistances[i]);
        if (!points[i].IsValid())
        {
            throw AlgorithmError("Splines::ComputePointOnSplineFromAdimensionalDistance: Could not interpolate spline points.");
        }
    }
    return {points, adimensionalDistances};
}

meshkernel::Point Splines::ComputeClosestPointOnSplineSegment(UInt index, double startSplineSegment, double endSplineSegment, Point point)
{
    FuncDistanceFromAPoint func(this, index, point);
    const auto adimensionalDistance = FindFunctionRootWithGoldenSectionSearch(func, startSplineSegment, endSplineSegment);
    return ComputePointOnSplineAtAdimensionalDistance(m_splineNodes[index], m_splineDerivatives[index], adimensionalDistance);
}
