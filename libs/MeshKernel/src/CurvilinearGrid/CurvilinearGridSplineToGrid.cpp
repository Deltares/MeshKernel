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

    Splines splinesCopy(splines);
    EigenMatrix<double> splineIntersections(splines.GetNumSplines(), DoubleVector<double>(splines.GetNumSplines(), 0.0));
    // EigenMatrix<double> splineIntersections(splines.GetNumSplines(), splines.GetNumSplines());

    UInt mFac = curvilinearParameters.m_refinement;
    UInt nFac = curvilinearParameters.n_refinement;

    AnotherMatrix splineInteraction(splines.GetNumSplines()); // mn12

    UInt numMSplines = 0;

    ComputeSplineIntersections(splinesCopy, splineIntersections, numMSplines);
    OrderSplines(splinesCopy, numMSplines, splineIntersections);
    ComputeSplineInteractions(splinesCopy, numMSplines, splineIntersections, splineInteraction);
    splrgf(splinesCopy, splineIntersections, splineInteraction, grid, numMSplines, mFac, nFac);
}

meshkernel::UInt meshkernel::CurvilinearGridSplineToGrid::longestSplineLength(const Splines& splines) const
{
    UInt result = 0;

    for (UInt i = 0; i < splines.GetNumSplines(); ++i)
    {
        result = std::max(result, splines.Size(i));
    }

    return result;
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

void meshkernel::CurvilinearGridSplineToGrid::determineIntersection(Splines& splines,
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
    doubledSplinePoints.reserve(2 * splines.MaxSize() - 1);

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

bool meshkernel::CurvilinearGridSplineToGrid::ComputeInteractions(Splines& splines,
                                                                  std::vector<int>& splineType,
                                                                  EigenMatrix<double>& splineIntersections) const
{
    UInt maxSplineSize = splines.Size(splines.MaxSize());

    std::fill(splineIntersections.begin(), splineIntersections.end(), std::vector<double>(splines.GetNumSplines(), 0.0));

    for (UInt splineI = 0; splineI < splines.GetNumSplines(); ++splineI)
    {
        for (UInt splineJ = splineI + 1; splineJ < splines.GetNumSplines(); ++splineJ)
        {
            double crossProduct;
            double normalisedIntersectionSplineI;
            double normalisedIntersectionSplineJ;
            UInt numberTimesCrossing = 0;

            // TODO need to know if two splines cross more than 1 time
            determineIntersection(splines, splineI, splineJ, numberTimesCrossing, crossProduct, normalisedIntersectionSplineI, normalisedIntersectionSplineJ);

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

#ifdef USE_EIGEN
                splineIntersections(splineI, splineJ) = normalisedIntersectionSplineI;
                splineIntersections(splineJ, splineI) = normalisedIntersectionSplineJ;
#else
                splineIntersections[splineI][splineJ] = normalisedIntersectionSplineI;
                splineIntersections[splineJ][splineI] = normalisedIntersectionSplineJ;
#endif
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
                                                                           EigenMatrix<double>& splineIntersections) const
{
    // sorteren op type, eerst de horizontalen (N = CONSTANT)
    for (UInt i = 0; i < splines.GetNumSplines(); ++i)
    {
        if (splineType[i] == -1)
        {
            for (UInt j = i + 1; j < splines.GetNumSplines(); ++j)
            {
                if (splineType[j] == 1)
                {
                    splines.SwapSplines(i, j);

#ifdef USE_EIGEN
                    splineIntersections.row(i).swap(splineIntersections.row(j));
                    splineIntersections.col(i).swap(splineIntersections.col(j));
#else
                    splineIntersections[i].swap(splineIntersections[j]);
                    SwapColumns(splineIntersections, i, j);
#endif

                    splineType[i] = 1;
                    splineType[j] = -1;
                    break;
                }
            }
        }
    }
}

meshkernel::UInt meshkernel::CurvilinearGridSplineToGrid::GetNumberOfSplinesInDirectionM(const std::vector<int>& splineType) const
{
    UInt numMSplines = constants::missing::uintValue;

    for (UInt i = 0; i < splineType.size(); ++i)
    {
        if (splineType[i] == -1)
        {
            numMSplines = i;
            break;
        }
    }

    return numMSplines;
}

// bool meshkernel::CurvilinearGridSplineToGrid::SortSplinesOnM(Splines& splines,
//                                                              const UInt outerStartIndex,
//                                                              const UInt outerEndIndex,
//                                                              const UInt innerStartIndex,
//                                                              const UInt innerEndIndex,
//                                                              EigenMatrix<double>& splineIntersections,
//                                                              UInt& count) const
// {

//     for (UInt i = outerStartIndex; i < outerEndIndex; ++i)
//     {
//         for (UInt j = innerStartIndex; j < innerendIndex; ++j)
//         {
// #ifdef USE_EIGEN
//             if (splineIntersections(i, j) != 0.0)
// #else
//             if (splineIntersections[i][j] != 0.0)
// #endif
//             {
//                 for (UInt k = j + 1; k < splines.GetNumSplines(); ++k)
//                 {
// #ifdef USE_EIGEN
//                     if (splineIntersections(i, k) != 0.0)
// #else
//                     if (splineIntersections[i][k] != 0.0)
// #endif
//                     {

// #ifdef USE_EIGEN
//                         if (splineIntersections(i, j) > splineIntersections(i, k))
// #else
//                         if (splineIntersections[i][j] > splineIntersections[i][k])
// #endif
//                         {
//                             splines.SwapSplines(j, k);

// #ifdef USE_EIGEN
//                             splineIntersections.row(j).swap(splineIntersections.row(k));
//                             splineIntersections.col(j).swap(splineIntersections.col(k));
// #else
//                             splineIntersections[j].swap(splineIntersections[k]);
//                             SwapColumns(splineIntersections, j, k);
// #endif

//                             ++count;

//                             if (count > splines.GetNumSplines())
//                             {
//                                 // error message: PROBLEM IN SPLINE ORDERING, MODIFY SPLINES
//                             }

//                             return true;
//                         }
//                     }
//                 }
//             }
//         }
//     }

//     return false;
// }

bool meshkernel::CurvilinearGridSplineToGrid::SortSplines(Splines& splines,
                                                          const UInt outerStartIndex,
                                                          const UInt outerEndIndex,
                                                          const UInt innerStartIndex,
                                                          const UInt innerEndIndex,
                                                          EigenMatrix<double>& splineIntersections,
                                                          bool& jaChange) const
{
    UInt count = 0;

    for (UInt i = outerStartIndex; i < outerEndIndex; ++i)
    {
        for (UInt j = innerStartIndex; j < innerEndIndex; ++j)
        {

#ifdef USE_EIGEN
            if (splineIntersections(i, j) != 0.0)
#else
            if (splineIntersections[i][j] != 0.0)
#endif
            {
                for (UInt k = j + 1; k < innerEndIndex; ++k)
                {
#ifdef USE_EIGEN
                    if (splineIntersections(i, k) != 0.0)
#else
                    if (splineIntersections[i][k] != 0.0)
#endif
                    {

#ifdef USE_EIGEN
                        if (splineIntersections(i, j) > splineIntersections(i, k))
#else
                        if (splineIntersections[i][j] > splineIntersections[i][k])
#endif
                        {
                            splines.SwapSplines(j, k);

#ifdef USE_EIGEN
                            splineIntersections.row(j).swap(splineIntersections.row(k));
                            splineIntersections.col(j).swap(splineIntersections.col(k));
#else
                            splineIntersections[j].swap(splineIntersections[k]);
                            SwapColumns(splineIntersections, j, k);
#endif

                            jaChange = true;
                            ++count;

                            if (count > splines.GetNumSplines())
                            {
                                throw AlgorithmError("Problem in spline ordering, modify splines");
                            }

                            return true;
                        }
                    }
                }
            }
        }
    }

    return false;
}

void meshkernel::CurvilinearGridSplineToGrid::OrderSplines(Splines& splines,
                                                           const UInt numMSplines,
                                                           EigenMatrix<double>& splineIntersections) const
{
    bool repeatBoth = true;
    UInt iterations = 0;
    const UInt maxIterations = 10 * splines.GetNumSplines();

    while (repeatBoth && iterations <= maxIterations)
    {
        ++iterations;
        repeatBoth = false;
        bool repeatM = true;
        UInt sortingIterations = 0;
        const UInt maxSortingIterations = 100 * splines.GetNumSplines();

        while (repeatM && sortingIterations <= maxSortingIterations)
        {
            ++sortingIterations;
            repeatBoth = false;
            repeatM = SortSplines(splines, 0, numMSplines, numMSplines, splines.GetNumSplines(), splineIntersections, repeatBoth);
        }

        bool repeatN = true;
        sortingIterations = 0;

        while (repeatN && sortingIterations <= maxSortingIterations)
        {
            ++sortingIterations;
            repeatN = SortSplines(splines, numMSplines, splines.GetNumSplines(), 0, numMSplines, splineIntersections, repeatBoth);
        }
    }

    if (iterations > maxIterations)
    {
        throw AlgorithmError("Problem in spline ordering, maximum iterations exceeded, modify splines");
    }
}

void meshkernel::CurvilinearGridSplineToGrid::ComputeSplineInteractions(const Splines& splines,
                                                                        const UInt numMSplines,
                                                                        const EigenMatrix<double>& splineIntersections,
                                                                        AnotherMatrix& splineInteraction) const
{

    std::ranges::fill(splineInteraction, std::array<int, 3>({0, 0, 0}));

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
#ifdef USE_EIGEN
                if (splineIntersections(j, k) != 0.0)
#else
                if (splineIntersections[j][k] != 0.0)
#endif
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
#ifdef USE_EIGEN
            if (splineIntersections(j, i) != 0.0)
#else
            if (splineIntersections[j][i] != 0.0)
#endif
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
#ifdef USE_EIGEN
                if (splineIntersections(j, k) != 0.0)
#else
                if (splineIntersections[j][k] != 0.0)
#endif
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
#ifdef USE_EIGEN
            if (splineIntersections(j, i) != 0.0)
#else
            if (splineIntersections[j][i] != 0.0)
#endif
            {
                maxm = std::max(maxm, splineInteraction[j][2]);
            }
        }

        splineInteraction[i][0] = maxm;
    }

    // splineInteraction.col(1).fill(0);
    // splineInteraction.col(2).fill(0);

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
#ifdef USE_EIGEN
            if (splineIntersections(i, j) != 0.0)
#else
            if (splineIntersections[i][j] != 0.0)
#endif
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
#ifdef USE_EIGEN
            if (splineIntersections(i, j) != 0.0)
#else
            if (splineIntersections[i][j] != 0.0)
#endif
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
                                                                         EigenMatrix<double>& splineIntersections,
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
            // Reset spline type info after doubleing of the number of spline points
            std::ranges::fill(splineType, 0);
            splineType[0] = 1;
        }

    } while (SplinesRemainUnlabeled(splineType, cumulativeUnlabelledSplineCount));

    SortInteractionsOnSplineType(splines, splineType, splineIntersections);
    numMSplines = GetNumberOfSplinesInDirectionM(splineType);
}

std::vector<double> meshkernel::CurvilinearGridSplineToGrid::CompressRow(const EigenMatrix<double>& splineIntersections,
                                                                         const UInt whichRow) const
{
    std::vector<double> compressedRow;
    compressedRow.reserve(splineIntersections.size());

#ifdef USE_EIGEN
    for (UInt col = 0; col < splineIntersections.cols(); ++col)
    {
        if (splineIntersections(whichRow, col) != 0.0)
        {
            compressedRow.push_back(splineIntersections(whichRow, col));
        }
    }
#else
    for (UInt col = 0; col < splineIntersections.size(); ++col)
    {
        if (splineIntersections[whichRow][col] != 0.0)
        {
            compressedRow.push_back(splineIntersections[whichRow][col]);
        }
    }
#endif

    return compressedRow;
}

void meshkernel::CurvilinearGridSplineToGrid::getdis(const Splines& splines,
                                                     const UInt whichSpline,
                                                     double& tValue,
                                                     double& sValue) const
{
    double dt = 0.1;
    double t0 = 0.0;
    Point startPoint = splines.m_splineNodes[whichSpline][0];
    Point endPoint;
    double t1;

    // in getdis.f90 this is: tValue = std::min(tValue, static_cast<double>(splines.m_splineNodes[whichSpline].size()))
    // without the -1. I think the fortran is incorrect, there should be the -1
    tValue = std::min(tValue, static_cast<double>(splines.Size(whichSpline) - 1));
    sValue = 0.0;

    do
    {
        t1 = t0 + dt;

        if (t1 < tValue)
        {
            endPoint = splines.Evaluate(whichSpline, t1);
        }
        else
        {
            // TODO try to remove the 1e-5
            endPoint = splines.Evaluate(whichSpline, tValue - 1.0e-5);
        }

        sValue += ComputeDistance(startPoint, endPoint, splines.m_projection);
        startPoint = endPoint;
        t0 = t1;

    } while (t1 < tValue);
}

void meshkernel::CurvilinearGridSplineToGrid::makesr(const double ar,
                                                     const double s0,
                                                     const double s1,
                                                     std::vector<double>& sr) const
{
    double ds = 1.0;
    sr[0] = 0.0;

    for (UInt k = 0; k < sr.size() - 1; ++k)
    {
        sr[k + 1] = sr[k] + ds;
        ds *= ar;
    }

    double fac = (s1 - s0) / sr[sr.size() - 1];

    for (UInt k = 0; k < sr.size(); ++k)
    {
        sr[k] = s0 + fac * sr[k];
    }
}

void meshkernel::CurvilinearGridSplineToGrid::makessq(const std::vector<double>& s, // intersectionPoints
                                                      const UInt mnFac,
                                                      std::vector<double>& ssq) const
{
    if (s.size() == 2)
    {
        for (UInt i = 0; i <= mnFac; ++i)
        {
            ssq[i] = s[0] + (s[1] - s[0]) * static_cast<double>(i) / static_cast<double>(mnFac);
        }
    }
    else if (s.size() >= 3)
    {
        std::vector<double> a(s.size());
        std::vector<double> sl(mnFac + 1);
        std::vector<double> sr(mnFac + 1);

        for (UInt i = 1; i < s.size() - 1; ++i)
        {
            a[i] = (s[i + 1] - s[i]) / (s[i] - s[i - 1]);
        }

        a[0] = a[1];
        a[s.size() - 1] = a[s.size() - 2];

        for (UInt i = 0; i < s.size() - 1; ++i)
        {
            // TODO store value rather that re-compute each time.
            double ar = std::pow(a[i + 1], 1.0 / static_cast<double>(mnFac));
            makesr(ar, s[i], s[i + 1], sr);

            double al = std::pow(a[i], 1.0 / static_cast<double>(mnFac));
            makesr(al, s[i], s[i + 1], sl);

            // TODO try k = 0 and mnFac
            for (UInt k = 1; k <= mnFac + 1; ++k)
            {
                UInt kr = i * mnFac + k - 1;
                ar = static_cast<double>(k - 1) / static_cast<double>(mnFac);
                al = 1.0 - ar;
                // TODO store in temporary variable
                double interp = ar * sr[k - 1] + al * sl[k - 1];
                ar = (interp - s[i]) / (s[i + 1] - s[i]);
                al = 1.0 - ar;
                ssq[kr] = ar * sr[k - 1] + al * sl[k - 1];
            }
        }
    }
}

meshkernel::Point meshkernel::CurvilinearGridSplineToGrid::GetXy(const Splines& splines,
                                                                 const UInt whichSpline,
                                                                 const std::vector<double>& intersectionPoints,
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

void meshkernel::CurvilinearGridSplineToGrid::GenerateGridPoints(const Splines& splines,
                                                                 const UInt whichSpline,
                                                                 const UInt mnFac,
                                                                 std::vector<double>& intersectionPoints,
                                                                 std::vector<Point>& gridPoints) const
{
    // evaluate the spline at the intersection points
    bool curvatureAdapted = false;  // Parameter H is not used in getdis.f90
    double maximumGridHeight = 0.0; // constants::missing::doubleValue; // what values should this be?

    auto [splinePoints, distances] = splines.ComputePointOnSplineFromAdimensionalDistance(whichSpline, maximumGridHeight, curvatureAdapted, intersectionPoints);

    //--------------------------------
    // getdis

    std::vector<double> s(intersectionPoints.size());

    for (UInt i = 0; i < intersectionPoints.size(); ++i)
    {
        getdis(splines, whichSpline, intersectionPoints[i], s[i]);
    }

    UInt kmax = (static_cast<UInt>(intersectionPoints.size()) - 1) * mnFac + 1;

    std::vector<double> ssq(kmax, -999.0);

    if (intersectionPoints.size() >= 2)
    {
        makessq(s, mnFac, ssq);

        for (UInt l = 1; l <= intersectionPoints.size() - 1; ++l)
        {
            UInt k1 = mnFac * (l - 1);
            UInt k2 = k1 + mnFac;
            bool jaDip = false;

        label_23:

            if (jaDip)
            {
                // TODO check indices and loop bounds
                for (UInt k = k1 + 1; k <= k2 - 1; ++k)
                {
                    ssq[k - 1] = 0.5 * (ssq[k - 2] + ssq[k + 0]);
                }
            }

            for (UInt k = k1; k < k2 - 1; ++k)
            {
                if (ssq[k + 1] < ssq[k])
                {
                    jaDip = true;
                    goto label_23;
                }
            }
        }
    }
    else
    {
        ssq[0] = intersectionPoints[0];
    }

    gridPoints.resize(kmax);

    for (UInt k = 0; k < kmax; ++k)
    {
        Point p = GetXy(splines, whichSpline, intersectionPoints, ssq[k]);
        gridPoints[k] = p;
    }
}

void meshkernel::CurvilinearGridSplineToGrid::assignBoundaryPoint(const UInt loopIndex,
                                                                  const UInt boundaryIndex,
                                                                  const UInt mnFac,
                                                                  std::vector<Point>& startBoundaryPoints,
                                                                  std::vector<Point>& endBoundaryPoints,
                                                                  const Point gridNode) const
{
    if (loopIndex == 0)
    {
        startBoundaryPoints[boundaryIndex] = gridNode;
    }
    else if (loopIndex == mnFac)
    {
        endBoundaryPoints[boundaryIndex] = gridNode;
    }
}

void meshkernel::CurvilinearGridSplineToGrid::splrgf(Splines& splines,
                                                     const EigenMatrix<double>& splineIntersections,
                                                     const AnotherMatrix& splineInteraction,
                                                     CurvilinearGrid& grid [[maybe_unused]],
                                                     const UInt numMSplines,
                                                     const UInt mFac,
                                                     const UInt nFac) const
{
    std::vector<Point> gridPoints;
    UInt cols = (splines.GetNumSplines() - numMSplines - 1) * mFac + 1;
    UInt rows = (numMSplines - 1) * nFac + 1;

    UInt ns = 0;
    UInt ms = 0;

    lin_alg::Matrix<Point> gridNodes(rows, cols);
    gridNodes.fill({constants::missing::doubleValue, constants::missing::doubleValue});

    // TODO can this be changed to just 'i'
    for (UInt splineIndex = 0; splineIndex < splines.GetNumSplines(); ++splineIndex)
    {
        std::vector<double> intersectionPoints = CompressRow(splineIntersections, splineIndex);
        UInt position;
        UInt startIndex;
        UInt endIndex;
        UInt refinementFactor;

        if (splineIndex < numMSplines)
        {
            refinementFactor = mFac;
            position = (splineInteraction[splineIndex][0] - 1) * nFac;
            startIndex = (splineInteraction[splineIndex][1] - 1) * mFac;
            endIndex = (splineInteraction[splineIndex][2] - 1) * mFac + 1;
        }
        else
        {
            refinementFactor = nFac;
            position = (splineInteraction[splineIndex][0] - 1) * mFac;
            startIndex = (splineInteraction[splineIndex][1] - 1) * nFac;
            endIndex = (splineInteraction[splineIndex][2] - 1) * nFac + 1;
        }

        GenerateGridPoints(splines, splineIndex, refinementFactor, intersectionPoints, gridPoints);

        UInt k = 0;

        // TODO change loop variable name, too confusing ii, splineIndex, startIndex.
        for (UInt ii = startIndex; ii < endIndex; ++ii)
        {

            if (k < gridPoints.size())
            {
                if (splineIndex < numMSplines)
                {
                    gridNodes(position, ii) = gridPoints[k];
                }
                else
                {
                    gridNodes(ii, position) = gridPoints[k];
                }
            }

            ++k;
        }

        if (splineIndex < numMSplines)
        {
            ns = std::max(ns, static_cast<UInt>(splineInteraction[splineIndex][0]));
        }
        else
        {
            ms = std::max(ms, static_cast<UInt>(splineInteraction[splineIndex][0]));
        }
    }

    UInt nmax = std::max((ns - 1) * mFac, (ms - 1) * nFac) + 1;

    std::vector<Point> eastBoundaryPoints(nmax);
    std::vector<Point> westBoundaryPoints(nmax);
    std::vector<Point> southBoundaryPoints(nmax);
    std::vector<Point> northBoundaryPoints(nmax);

    for (UInt i = 0; i < ms - 1; ++i)
    {
        for (UInt j = 0; j < ns - 1; ++j)
        {
            eastBoundaryPoints[1].SetInvalid();
            westBoundaryPoints[1].SetInvalid();
            southBoundaryPoints[1].SetInvalid();
            northBoundaryPoints[1].SetInvalid();

            for (UInt k = 0; k <= mFac; ++k)
            {
                UInt ki = i * mFac + k;

                for (UInt l = 0; l <= nFac; ++l)
                {
                    UInt lj = j * nFac + l;
                    Point gridNode = gridNodes(lj, ki);

                    if (gridNode.IsValid())
                    {
                        assignBoundaryPoint(k, l, mFac, eastBoundaryPoints, westBoundaryPoints, gridNode);
                        assignBoundaryPoint(l, k, nFac, southBoundaryPoints, northBoundaryPoints, gridNode);
                    }
                }
            }

            bool isValid = eastBoundaryPoints[1].IsValid() &&
                           westBoundaryPoints[1].IsValid() &&
                           southBoundaryPoints[1].IsValid() &&
                           northBoundaryPoints[1].IsValid();

            if (isValid)
            {
                auto interpolationResult = DiscretizeTransfinite(eastBoundaryPoints,
                                                                 westBoundaryPoints,
                                                                 southBoundaryPoints,
                                                                 northBoundaryPoints,
                                                                 splines.m_projection,
                                                                 nFac,
                                                                 mFac);

                for (UInt k2 = 0; k2 <= mFac; ++k2)
                {
                    UInt ki = i * mFac + k2;

                    for (UInt l2 = 1; l2 <= nFac; ++l2)
                    {
                        UInt lj = j * nFac + l2;

                        if (!gridNodes(lj, ki).IsValid())
                        {
                            gridNodes(lj, ki) = interpolationResult(k2, l2);
                        }
                    }
                }
            }
        }
    }

    grid.SetGridNodes(gridNodes);
}
