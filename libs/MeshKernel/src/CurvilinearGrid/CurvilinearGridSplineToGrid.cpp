#include "MeshKernel/CurvilinearGrid/CurvilinearGridSplineToGrid.hpp"
#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Utilities/LinearAlgebra.hpp"

#include <algorithm>
#include <iomanip>

void meshkernel::CurvilinearGridSplineToGrid::Compute(const Splines& splines,
                                                      const CurvilinearParameters& curvilinearParameters,
                                                      CurvilinearGrid& grid) const
{

    if (splines.GetNumSplines() < 4)
    {
        throw ConstraintError("At least 4 splines are required to generate a grid, number of splines is {}", splines.GetNumSplines());
    }

    if (!CheckSplines(splines))
    {
        throw ConstraintError("One or more splines has one or fewer support points");
    }

    // Make copy of spline because the may be changed
    // i.e. the order may change and the number of spline points may increase.
    Splines splinesCopy(splines);
    VectorOfDoubleVectors splineIntersections(splines.GetNumSplines(), DoubleVector(splines.GetNumSplines(), 0.0));

    UInt mRefinement = curvilinearParameters.m_refinement;
    UInt nRefinement = curvilinearParameters.n_refinement;

    // mn12 in Fortran.
    VectorOfThreeInts splineInteraction(splines.GetNumSplines());
    UInt numMSplines = 0;

    ComputeSplineIntersections(splinesCopy, splineIntersections, numMSplines);
    OrderSplines(splinesCopy, numMSplines, splineIntersections);
    DetermineSplineOrientation(splinesCopy, numMSplines, splineIntersections, splineInteraction);
    GenerateGrid(splinesCopy, splineIntersections, splineInteraction, numMSplines, mRefinement, nRefinement, grid);
}

meshkernel::CurvilinearGrid meshkernel::CurvilinearGridSplineToGrid::Compute(const Splines& splines,
                                                                             const CurvilinearParameters& curvilinearParameters) const
{
    CurvilinearGrid grid(splines.m_projection);
    Compute(splines, curvilinearParameters, grid);
    return grid;
}

bool meshkernel::CurvilinearGridSplineToGrid::CheckSplines(const Splines& splines) const
{
    for (UInt splineIndex = 0; splineIndex < splines.GetNumSplines(); ++splineIndex)
    {
        if (splines.Size(splineIndex) < 2)
        {
            return false;
        }
    }

    return true;
}

void meshkernel::CurvilinearGridSplineToGrid::DetermineIntersection(Splines& splines,
                                                                    const UInt splineI,
                                                                    const UInt splineJ,
                                                                    UInt& numberTimesCrossing,
                                                                    double& crossProductOfIntersection,
                                                                    double& firstNormalisedIntersectionLength,
                                                                    double& secondNormalisedIntersectionLength) const
{
    Point intersectionPoint;

    // TODO Need way of getting the number of times the spline intersect.
    if (splines.GetSplinesIntersection(splineI, splineJ, crossProductOfIntersection, intersectionPoint, firstNormalisedIntersectionLength, secondNormalisedIntersectionLength))
    {
        numberTimesCrossing = 1;
    }
    else
    {
        secondNormalisedIntersectionLength = constants::missing::doubleValue;
        firstNormalisedIntersectionLength = constants::missing::doubleValue;
        crossProductOfIntersection = constants::missing::doubleValue;
        numberTimesCrossing = 0;
    }
}

void meshkernel::CurvilinearGridSplineToGrid::DoubleSplinePoints(Splines& splines) const
{
    std::vector<Point> doubledSplinePoints;
    doubledSplinePoints.reserve(2 * splines.Size(splines.MaxSizeIndex()) - 1);

    for (UInt splineIndex = 0; splineIndex < splines.GetNumSplines(); ++splineIndex)
    {
        doubledSplinePoints.resize(2 * splines.Size(splineIndex) - 1);

        // Copy original spline points to every second point
        for (UInt j = 0; j < splines.Size(splineIndex); ++j)
        {
            doubledSplinePoints[2 * j] = splines.m_splineNodes[splineIndex][j];
        }

        // Evaluate the spline at the intermediate points
        for (UInt j = 1; j < splines.Size(splineIndex); ++j)
        {
            double segmentMidPoint = static_cast<double>(j) - 0.5;
            doubledSplinePoints[2 * j - 1] = splines.Evaluate(splineIndex, segmentMidPoint);
        }

        // Replace the spline points with the new double size point set.
        splines.Replace(splineIndex, doubledSplinePoints);
    }
}

bool meshkernel::CurvilinearGridSplineToGrid::ComputeAndCheckIntersection(Splines& splines,
                                                                          const UInt splineI,
                                                                          const UInt splineJ,
                                                                          std::vector<int>& splineType,
                                                                          VectorOfDoubleVectors& splineIntersections) const
{
    UInt maxSplineSize = splines.Size(splines.MaxSizeIndex());

    double crossProduct;
    double normalisedIntersectionSplineI;
    double normalisedIntersectionSplineJ;
    UInt numberTimesCrossing = 0;

    // TODO need to know if two splines cross more than 1 time
    DetermineIntersection(splines, splineI, splineJ, numberTimesCrossing, crossProduct, normalisedIntersectionSplineI, normalisedIntersectionSplineJ);

    if (numberTimesCrossing == 1)
    {

        if (splineType[splineI] * splineType[splineJ] == 1)
        {
            if (maxSplineSize > MaximumNumberOfSplinePoints / 2)
            {
                throw AlgorithmError("splines both in m- and n-direction: spline {} and {}", splineI + 1, splineJ + 1);
            }
            else
            {
                return true;
            }
        }
        else if (splineType[splineI] == 0 && splineType[splineJ] == 0)
        {
            throw AlgorithmError("Both spline {} and {} are undefined", splineI + 1, splineJ + 1);
        }
        else if (splineType[splineJ] == 0)
        {
            splineType[splineJ] = -splineType[splineI];

            if (crossProduct * static_cast<double>(splineType[splineI]) < 0.0)
            {
                splines.Reverse(splineJ);
                // Reverse the normalised distance
                normalisedIntersectionSplineJ = static_cast<double>(splines.Size(splineJ)) - 1.0 - normalisedIntersectionSplineJ;
            }
        }
        else if (splineType[splineI] == 0)
        {
            splineType[splineI] = -splineType[splineJ];

            if (crossProduct * static_cast<double>(splineType[splineJ]) > 0.0)
            {
                splines.Reverse(splineI);
                // Reverse the normalised distance
                normalisedIntersectionSplineI = static_cast<double>(splines.Size(splineI)) - 1.0 - normalisedIntersectionSplineI;
            }
        }

        splineIntersections[splineI][splineJ] = normalisedIntersectionSplineI;
        splineIntersections[splineJ][splineI] = normalisedIntersectionSplineJ;
    }
    else if (numberTimesCrossing > 1)
    {
        if (maxSplineSize > MaximumNumberOfSplinePoints / 2)
        {
            throw AlgorithmError("Splines {} and {}, appear to intersect more that one time", splineI + 1, splineJ + 1);
        }
        else
        {
            return true;
        }
    }

    return false;
}

bool meshkernel::CurvilinearGridSplineToGrid::ComputeInteractions(Splines& splines,
                                                                  std::vector<int>& splineType,
                                                                  VectorOfDoubleVectors& splineIntersections) const
{

    std::ranges::fill(splineIntersections, DoubleVector(splines.GetNumSplines(), 0.0));

    for (UInt splineI = 0; splineI < splines.GetNumSplines(); ++splineI)
    {
        for (UInt splineJ = splineI + 1; splineJ < splines.GetNumSplines(); ++splineJ)
        {
            if (ComputeAndCheckIntersection(splines, splineI, splineJ, splineType, splineIntersections))
            {
                return true;
            }
        }
    }

    return false;
}

bool meshkernel::CurvilinearGridSplineToGrid::SplinesRemainUnlabeled(const std::vector<int>& splineType, UInt& unlabledSplineCount) const
{
    for (UInt splineIndex = 0; splineIndex < splineType.size(); ++splineIndex)
    {
        if (splineType[splineIndex] == 0)
        {
            ++unlabledSplineCount;

            if (unlabledSplineCount > MaximumCumulativeUnlabeledSplineCount)
            {
                throw AlgorithmError("Spline number {} cannot be attached to the grid", splineIndex + 1);
            }

            return true;
        }
    }

    return false;
}

void meshkernel::CurvilinearGridSplineToGrid::SortInteractionsOnSplineType(Splines& splines,
                                                                           std::vector<int>& splineType,
                                                                           VectorOfDoubleVectors& splineIntersections) const
{
    // Sort splines on type, first in list are the horizontal spines (n = constant)
    for (UInt i = 0; i < splines.GetNumSplines(); ++i)
    {
        if (splineType[i] == -1)
        {
            for (UInt j = i + 1; j < splines.GetNumSplines(); ++j)
            {
                if (splineType[j] == 1)
                {
                    splines.SwapSplines(i, j);

                    splineIntersections[i].swap(splineIntersections[j]);
                    SwapColumns(splineIntersections, i, j);

                    splineType[i] = 1;
                    splineType[j] = -1;
                    break;
                }
            }
        }
    }
}

bool meshkernel::CurvilinearGridSplineToGrid::SortSplines(Splines& splines,
                                                          const UInt outerStartIndex,
                                                          const UInt outerEndIndex,
                                                          const UInt innerStartIndex,
                                                          const UInt innerEndIndex,
                                                          VectorOfDoubleVectors& splineIntersections,
                                                          bool& splinesSwapped) const
{
    UInt count = 0;

    for (UInt i = outerStartIndex; i < outerEndIndex; ++i)
    {
        for (UInt j = innerStartIndex; j < innerEndIndex; ++j)
        {
            if (splineIntersections[i][j] == 0.0)
            {
                continue;
            }

            for (UInt k = j + 1; k < innerEndIndex; ++k)
            {
                if (splineIntersections[i][k] == 0.0)
                {
                    continue;
                }

                if (splineIntersections[i][j] > splineIntersections[i][k])
                {
                    splines.SwapSplines(j, k);

                    splineIntersections[j].swap(splineIntersections[k]);
                    SwapColumns(splineIntersections, j, k);

                    splinesSwapped = true;
                    ++count;

                    if (count > splines.GetNumSplines())
                    {
                        throw AlgorithmError("Problem in spline ordering, modify splines");
                    }

                    return false;
                }
            }
        }
    }

    return true;
}

void meshkernel::CurvilinearGridSplineToGrid::OrderSplines(Splines& splines,
                                                           const UInt numMSplines,
                                                           VectorOfDoubleVectors& splineIntersections) const
{
    bool repeatBoth = true;
    UInt iterations = 0;
    const UInt maxIterations = 10 * splines.GetNumSplines();

    while (repeatBoth && iterations <= maxIterations)
    {
        ++iterations;
        repeatBoth = false;
        bool splinesSorted = false;
        UInt sortingIterations = 0;
        const UInt maxSortingIterations = 100 * splines.GetNumSplines();

        while (!splinesSorted && sortingIterations <= maxSortingIterations)
        {
            ++sortingIterations;
            repeatBoth = false;
            splinesSorted = SortSplines(splines, 0, numMSplines, numMSplines, splines.GetNumSplines(), splineIntersections, repeatBoth);
        }

        sortingIterations = 0;
        splinesSorted = false;

        while (!splinesSorted && sortingIterations <= maxSortingIterations)
        {
            ++sortingIterations;
            splinesSorted = SortSplines(splines, numMSplines, splines.GetNumSplines(), 0, numMSplines, splineIntersections, repeatBoth);
        }
    }

    if (iterations > maxIterations)
    {
        throw AlgorithmError("Problem in spline ordering, maximum iterations exceeded, modify splines");
    }
}

void meshkernel::CurvilinearGridSplineToGrid::DetermineSplineOrientation(const Splines& splines,
                                                                         const UInt numMSplines,
                                                                         const VectorOfDoubleVectors& splineIntersections,
                                                                         VectorOfThreeInts& splineInteraction) const
{

    std::ranges::fill(splineInteraction, ArrayOfThree({0, 0, 0}));

    // Eerst alles ranken in N richting
    for (UInt i = 0; i < numMSplines; ++i)
    {
        int maxn = 0;

        for (UInt j = numMSplines; j < splines.GetNumSplines(); ++j)
        {
            maxn = 0;
            UInt jjlast = 0;

            for (UInt k = 0; k <= i; ++k)
            {
                if (splineIntersections[j][k] != 0.0)
                {
                    maxn = splineInteraction[jjlast][0] + 1;
                    jjlast = k;
                }
            }

            splineInteraction[j][1] = maxn;
        }

        maxn = 0;

        for (UInt j = numMSplines; j < splines.GetNumSplines(); ++j)
        {
            if (splineIntersections[j][i] != 0.0)
            {
                maxn = std::max(maxn, splineInteraction[j][1]);
            }
        }

        splineInteraction[i][0] = maxn;
    }

    // Dan alles ranken in M richting
    for (UInt i = numMSplines; i < splines.GetNumSplines(); ++i)
    {
        int maxm = 0;

        for (UInt j = 0; j < numMSplines; ++j)
        {
            maxm = 0;
            UInt iilast = numMSplines;

            for (UInt k = numMSplines; k <= i; ++k)
            {
                if (splineIntersections[j][k] != 0.0)
                {
                    maxm = splineInteraction[iilast][0] + 1;
                    iilast = k;
                }
            }

            splineInteraction[j][2] = maxm;
        }

        maxm = 0;

        for (UInt j = 0; j < numMSplines; ++j)
        {
            if (splineIntersections[j][i] != 0.0)
            {
                maxm = std::max(maxm, splineInteraction[j][2]);
            }
        }

        splineInteraction[i][0] = maxm;
    }

    for (UInt i = 0; i < splines.GetNumSplines(); ++i)
    {
        splineInteraction[i][1] = 0;
        splineInteraction[i][2] = 0;
    }

    // Daarna per spline begin- en eindpunt tellen, eerst N = constant

    for (UInt i = 0; i < numMSplines; ++i)
    {
        for (UInt j = numMSplines; j < splines.GetNumSplines(); ++j)
        {
            if (splineIntersections[i][j] != 0.0)
            {
                if (splineInteraction[i][1] == 0)
                {
                    splineInteraction[i][1] = splineInteraction[j][0];
                }

                splineInteraction[i][2] = splineInteraction[j][0];
            }
        }
    }

    // Dan M = constant

    for (UInt i = numMSplines; i < splines.GetNumSplines(); ++i)
    {
        for (UInt j = 0; j < numMSplines; ++j)
        {
            if (splineIntersections[i][j] != 0.0)
            {
                if (splineInteraction[i][1] == 0)
                {
                    splineInteraction[i][1] = splineInteraction[j][0];
                }

                splineInteraction[i][2] = splineInteraction[j][0];
            }
        }
    }
}

void meshkernel::CurvilinearGridSplineToGrid::ComputeSplineIntersections(Splines& splines,
                                                                         VectorOfDoubleVectors& splineIntersections,
                                                                         UInt& numMSplines) const
{
    std::vector<int> splineType(splines.GetNumSplines(), 0);

    UInt cumulativeUnlabelledSplineCount = 0;

    std::ranges::fill(splineType, 0);
    splineType[0] = 1;

    do
    {

        while (ComputeInteractions(splines, splineType, splineIntersections))
        {
            cumulativeUnlabelledSplineCount = 0;
            DoubleSplinePoints(splines);
            // Reset spline type info after doubling of the number of spline points
            std::ranges::fill(splineType, 0);
            splineType[0] = 1;
        }

    } while (SplinesRemainUnlabeled(splineType, cumulativeUnlabelledSplineCount));

    SortInteractionsOnSplineType(splines, splineType, splineIntersections);
    numMSplines = FindIndex(splineType, -1);
}

std::vector<double> meshkernel::CurvilinearGridSplineToGrid::CompressRow(const VectorOfDoubleVectors& splineIntersections,
                                                                         const UInt whichRow) const
{
    DoubleVector compressedRow;
    compressedRow.reserve(splineIntersections.size());

    for (UInt col = 0; col < splineIntersections.size(); ++col)
    {
        if (splineIntersections[whichRow][col] != 0.0)
        {
            compressedRow.push_back(splineIntersections[whichRow][col]);
        }
    }

    return compressedRow;
}

void meshkernel::CurvilinearGridSplineToGrid::ComputeSplineIntervalLength(const Splines& splines,
                                                                          const UInt whichSpline,
                                                                          double& normalisedDistance,
                                                                          double& intervalLength) const
{
    Point startPoint = splines.m_splineNodes[whichSpline][0];
    Point endPoint;
    double stepSize = 0.1;
    double intervalStart = 0.0;
    double intervalEnd;

    // in getdis.f90 this is: normalisedDistance = std::min(normalisedDistance, static_cast<double>(splines.m_splineNodes[whichSpline].size()))
    // without the -1. I think the fortran is incorrect, there should be the -1
    normalisedDistance = std::min(normalisedDistance, static_cast<double>(splines.Size(whichSpline) - 1));
    intervalLength = 0.0;

    do
    {
        intervalEnd = intervalStart + stepSize;

        if (intervalEnd < normalisedDistance)
        {
            endPoint = splines.Evaluate(whichSpline, intervalEnd);
        }
        else
        {
            endPoint = splines.Evaluate(whichSpline, normalisedDistance);
        }

        intervalLength += ComputeDistance(startPoint, endPoint, splines.m_projection);
        startPoint = endPoint;
        intervalStart = intervalEnd;

    } while (intervalEnd < normalisedDistance);
}

void meshkernel::CurvilinearGridSplineToGrid::ComputeExponentialDistances(const double factor,
                                                                          const double intervalStart,
                                                                          const double intervalEnd,
                                                                          DoubleVector& intervalDiscretisation) const
{
    double step = 1.0;
    intervalDiscretisation[0] = 0.0;

    for (UInt k = 0; k < intervalDiscretisation.size() - 1; ++k)
    {
        intervalDiscretisation[k + 1] = intervalDiscretisation[k] + step;
        step *= factor;
    }

    double stepSize = (intervalEnd - intervalStart) / intervalDiscretisation[intervalDiscretisation.size() - 1];

    for (double& distance : intervalDiscretisation)
    {
        distance = intervalStart + stepSize * distance;
    }
}

void meshkernel::CurvilinearGridSplineToGrid::ComputeDiscretisedDistances(const DoubleVector& splineIntervalLengths,
                                                                          const UInt refinementFactor,
                                                                          DoubleVector& discretisedDistances) const
{

    const double refinementFactorInv = 1.0 / static_cast<double>(refinementFactor);

    if (splineIntervalLengths.size() == 2)
    {
        for (UInt i = 0; i <= refinementFactor; ++i)
        {
            discretisedDistances[i] = splineIntervalLengths[0] + (splineIntervalLengths[1] - splineIntervalLengths[0]) * static_cast<double>(i) * refinementFactorInv;
        }
    }
    else if (splineIntervalLengths.size() >= 3)
    {
        DoubleVector intervalRatios(splineIntervalLengths.size());
        DoubleVector leftDiscretisations(refinementFactor + 1);
        DoubleVector rightDiscretisations(refinementFactor + 1);

        for (UInt i = 1; i < splineIntervalLengths.size() - 1; ++i)
        {
            intervalRatios[i] = (splineIntervalLengths[i + 1] - splineIntervalLengths[i]) / (splineIntervalLengths[i] - splineIntervalLengths[i - 1]);
        }

        intervalRatios[0] = intervalRatios[1];
        intervalRatios[splineIntervalLengths.size() - 1] = intervalRatios[splineIntervalLengths.size() - 2];

        for (UInt i = 0; i < splineIntervalLengths.size() - 1; ++i)
        {
            double rightFactor = std::pow(intervalRatios[i + 1], refinementFactorInv);
            ComputeExponentialDistances(rightFactor, splineIntervalLengths[i], splineIntervalLengths[i + 1], rightDiscretisations);

            double leftFactor = std::pow(intervalRatios[i], refinementFactorInv);
            ComputeExponentialDistances(leftFactor, splineIntervalLengths[i], splineIntervalLengths[i + 1], leftDiscretisations);

            for (UInt k = 0; k <= refinementFactor; ++k)
            {
                UInt kr = i * refinementFactor + k;
                double rightWeight = static_cast<double>(k) * refinementFactorInv;
                double leftWeight = 1.0 - rightWeight;
                double interp = rightWeight * rightDiscretisations[k] + leftWeight * leftDiscretisations[k];
                rightWeight = (interp - splineIntervalLengths[i]) / (splineIntervalLengths[i + 1] - splineIntervalLengths[i]);
                leftWeight = 1.0 - rightWeight;
                discretisedDistances[kr] = rightWeight * rightDiscretisations[k] + leftWeight * leftDiscretisations[k];
            }
        }
    }
}

meshkernel::Point meshkernel::CurvilinearGridSplineToGrid::ComputePoint(const Splines& splines,
                                                                        const UInt whichSpline,
                                                                        const DoubleVector& intersectionPoints,
                                                                        const double ssq) const
{
    Point result;

    double ax = intersectionPoints[0];
    double cx = intersectionPoints[intersectionPoints.size() - 1];

    FuncAdimensionalToDimensionalDistanceOnSpline func(splines, whichSpline, false, 0.0, ssq);
    double distance = FindFunctionRootWithGoldenSectionSearch(func, ax, cx);

    result = splines.Evaluate(whichSpline, distance);

    return result;
}

void meshkernel::CurvilinearGridSplineToGrid::PrepareNormalisedDistances(const UInt intervalStart,
                                                                         const UInt intervalEnd,
                                                                         DoubleVector& normalisedDistances) const
{
    bool valuesAreDecreasing = false;

    do
    {
        if (valuesAreDecreasing)
        {
            valuesAreDecreasing = false;

            for (UInt k = intervalStart + 1; k <= intervalEnd - 1; ++k)
            {
                normalisedDistances[k - 1] = 0.5 * (normalisedDistances[k - 2] + normalisedDistances[k]);
            }
        }

        for (UInt k = intervalStart; k < intervalEnd - 1; ++k)
        {
            if (normalisedDistances[k + 1] < normalisedDistances[k])
            {
                valuesAreDecreasing = true;
                break;
            }
        }

    } while (valuesAreDecreasing);
}

void meshkernel::CurvilinearGridSplineToGrid::GenerateGridPoints(const Splines& splines,
                                                                 const UInt whichSpline,
                                                                 const UInt refinementFactor,
                                                                 DoubleVector& intersectionPoints,
                                                                 std::vector<Point>& gridPoints) const
{
    // evaluate the spline at the intersection points
    bool curvatureAdapted = false;  // Parameter H is not used in getdis.f90
    double maximumGridHeight = 0.0; // constants::missing::doubleValue; // what values should this be?

    auto [splinePoints, distances] = splines.ComputePointOnSplineFromAdimensionalDistance(whichSpline, maximumGridHeight, curvatureAdapted, intersectionPoints);

    // The length of the spline at the intersection points
    DoubleVector splineIntervalLengths(intersectionPoints.size());

    for (UInt i = 0; i < intersectionPoints.size(); ++i)
    {
        ComputeSplineIntervalLength(splines, whichSpline, intersectionPoints[i], splineIntervalLengths[i]);
    }

    UInt numberOfPoints = (static_cast<UInt>(intersectionPoints.size()) - 1) * refinementFactor + 1;

    DoubleVector normalisedDistances(numberOfPoints, -999.0);

    if (intersectionPoints.size() >= 2)
    {
        ComputeDiscretisedDistances(splineIntervalLengths, refinementFactor, normalisedDistances);

        for (UInt l = 1; l <= intersectionPoints.size() - 1; ++l)
        {
            UInt intervalStart = refinementFactor * (l - 1);
            UInt intervalEnd = intervalStart + refinementFactor;

            PrepareNormalisedDistances(intervalStart, intervalEnd, normalisedDistances);
        }
    }
    else
    {
        normalisedDistances[0] = intersectionPoints[0];
    }

    gridPoints.resize(numberOfPoints);

    for (UInt k = 0; k < numberOfPoints; ++k)
    {
        Point p = ComputePoint(splines, whichSpline, intersectionPoints, normalisedDistances[k]);
        gridPoints[k] = p;
    }
}

void meshkernel::CurvilinearGridSplineToGrid::AssignBoundaryPoint(const UInt loopIndex,
                                                                  const UInt boundaryIndex,
                                                                  const UInt refinementFactor,
                                                                  std::vector<Point>& startBoundaryPoints,
                                                                  std::vector<Point>& endBoundaryPoints,
                                                                  const Point gridNode) const
{
    if (loopIndex == 0)
    {
        startBoundaryPoints[boundaryIndex] = gridNode;
    }
    else if (loopIndex == refinementFactor)
    {
        endBoundaryPoints[boundaryIndex] = gridNode;
    }
}

void meshkernel::CurvilinearGridSplineToGrid::AssignPatchGridPoints(const UInt i, const UInt j, const UInt nRefinement, const UInt mRefinement,
                                                                    const lin_alg::Matrix<Point>& patchNodes,
                                                                    lin_alg::Matrix<Point>& gridNodes) const
{

    for (UInt k2 = 0; k2 <= mRefinement; ++k2)
    {
        UInt ki = i * mRefinement + k2;

        for (UInt l2 = 1; l2 <= nRefinement; ++l2)
        {
            UInt lj = j * nRefinement + l2;

            if (!gridNodes(lj, ki).IsValid())
            {
                gridNodes(lj, ki) = patchNodes(k2, l2);
            }
        }
    }
}

void meshkernel::CurvilinearGridSplineToGrid::GenerateGridPointsAlongSpline(const Splines& splines,
                                                                            const VectorOfDoubleVectors& splineIntersections,
                                                                            const VectorOfThreeInts& splineInteraction,
                                                                            const UInt numMSplines,
                                                                            const UInt mRefinement,
                                                                            const UInt nRefinement,
                                                                            lin_alg::Matrix<Point>& gridNodes) const
{
    std::vector<Point> gridPoints;
    gridPoints.reserve(static_cast<UInt>(std::max(gridNodes.rows(), gridNodes.cols())));

    for (UInt splineIndex = 0; splineIndex < splines.GetNumSplines(); ++splineIndex)
    {
        DoubleVector intersectionPoints(CompressRow(splineIntersections, splineIndex));
        UInt position;
        UInt startIndex;
        UInt endIndex;
        UInt refinementFactor;

        if (splineIndex < numMSplines)
        {
            refinementFactor = mRefinement;
            position = (splineInteraction[splineIndex][0] - 1) * nRefinement;
            startIndex = (splineInteraction[splineIndex][1] - 1) * mRefinement;
            endIndex = (splineInteraction[splineIndex][2] - 1) * mRefinement + 1;
        }
        else
        {
            refinementFactor = nRefinement;
            position = (splineInteraction[splineIndex][0] - 1) * mRefinement;
            startIndex = (splineInteraction[splineIndex][1] - 1) * nRefinement;
            endIndex = (splineInteraction[splineIndex][2] - 1) * nRefinement + 1;
        }

        GenerateGridPoints(splines, splineIndex, refinementFactor, intersectionPoints, gridPoints);

        UInt index = 0;

        if (splineIndex < numMSplines)
        {
            for (UInt i = startIndex; i < endIndex; ++i)
            {
                gridNodes(position, i) = gridPoints[index];
                ++index;

                if (index == gridPoints.size())
                {
                    break;
                }
            }
        }
        else
        {
            for (UInt i = startIndex; i < endIndex; ++i)
            {
                gridNodes(i, position) = gridPoints[index];
                ++index;

                if (index == gridPoints.size())
                {
                    break;
                }
            }
        }
    }
}

void meshkernel::CurvilinearGridSplineToGrid::FillPatchesWithPoints(const Splines& splines,
                                                                    const UInt numMSplines,
                                                                    const UInt mRefinement,
                                                                    const UInt nRefinement,
                                                                    lin_alg::Matrix<Point>& gridNodes) const
{
    UInt numNSplines = splines.GetNumSplines() - numMSplines;
    UInt maximumPointCount = std::max((numMSplines - 1) * mRefinement, (numNSplines - 1) * nRefinement) + 1;

    std::vector<Point> eastBoundaryPoints(maximumPointCount);
    std::vector<Point> westBoundaryPoints(maximumPointCount);
    std::vector<Point> southBoundaryPoints(maximumPointCount);
    std::vector<Point> northBoundaryPoints(maximumPointCount);

    for (UInt i = 0; i < numNSplines - 1; ++i)
    {
        for (UInt j = 0; j < numMSplines - 1; ++j)
        {
            eastBoundaryPoints[1].SetInvalid();
            westBoundaryPoints[1].SetInvalid();
            southBoundaryPoints[1].SetInvalid();
            northBoundaryPoints[1].SetInvalid();

            for (UInt k = 0; k <= mRefinement; ++k)
            {
                UInt ki = i * mRefinement + k;

                for (UInt l = 0; l <= nRefinement; ++l)
                {
                    UInt lj = j * nRefinement + l;
                    Point gridNode = gridNodes(lj, ki);

                    if (!gridNode.IsValid())
                    {
                        continue;
                    }

                    AssignBoundaryPoint(k, l, mRefinement, eastBoundaryPoints, westBoundaryPoints, gridNode);
                    AssignBoundaryPoint(l, k, nRefinement, southBoundaryPoints, northBoundaryPoints, gridNode);
                }
            }

            bool isValid = eastBoundaryPoints[1].IsValid() &&
                           westBoundaryPoints[1].IsValid() &&
                           southBoundaryPoints[1].IsValid() &&
                           northBoundaryPoints[1].IsValid();

            if (!isValid)
            {
                continue;
            }

            auto interpolationResult = DiscretizeTransfinite(eastBoundaryPoints,
                                                             westBoundaryPoints,
                                                             southBoundaryPoints,
                                                             northBoundaryPoints,
                                                             splines.m_projection,
                                                             nRefinement,
                                                             mRefinement);

            AssignPatchGridPoints(i, j, nRefinement, mRefinement, interpolationResult, gridNodes);
        }
    }
}

void meshkernel::CurvilinearGridSplineToGrid::GenerateGrid(const Splines& splines,
                                                           const VectorOfDoubleVectors& splineIntersections,
                                                           const VectorOfThreeInts& splineInteraction,
                                                           const UInt numMSplines,
                                                           const UInt mRefinement,
                                                           const UInt nRefinement,
                                                           CurvilinearGrid& grid) const
{
    UInt cols = (splines.GetNumSplines() - numMSplines - 1) * mRefinement + 1;
    UInt rows = (numMSplines - 1) * nRefinement + 1;

    lin_alg::Matrix<Point> gridNodes(rows, cols);
    gridNodes.fill({constants::missing::doubleValue, constants::missing::doubleValue});

    GenerateGridPointsAlongSpline(splines, splineIntersections, splineInteraction, numMSplines, mRefinement, nRefinement, gridNodes);
    FillPatchesWithPoints(splines, numMSplines, mRefinement, nRefinement, gridNodes);

    grid.SetGridNodes(gridNodes);
}
