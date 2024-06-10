#include "MeshKernel/CurvilinearGrid/CurvilinearGridSplineToGrid.hpp"
#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Utilities/LinearAlgebra.hpp"

#include <algorithm>
#include <iomanip>

void meshkernel::CurvilinearGridSplineToGrid::Compute(const Splines& splines, CurvilinearGrid& grid [[maybe_unused]]) const
{
    Splines splinesCopy(splines);
    EigenMatrix<double> splineIntersections(splines.GetNumSplines(), DoubleVector<double>(splines.GetNumSplines(), 0.0));
    // EigenMatrix<double> splineIntersections(splines.GetNumSplines(), splines.GetNumSplines());
    std::cout << "Compute:: number of splines: " << splines.GetNumSplines() << std::endl;

    AnotherMatrix mn12(splines.GetNumSplines());

    UInt numi;

    sectr(splinesCopy, splineIntersections, mn12, numi);
    splrgf(splinesCopy, splineIntersections, mn12, grid, numi);
    std::cout << std::endl;
    std::cout.precision(15);

    // for (UInt i = 0; i < splines.GetNumSplines(); ++i)
    // {
    //     for (UInt j = 0; j < splines.GetNumSplines(); ++j)
    //     {
    //         std::cout << std::setw(15) << i + 1 << "  " << std::setw(15) << j + 1 << "  " << splineIntersections(i, j) << "  ";
    //         std::cout << std::endl;
    //     }

    //     // std::cout << std::endl;
    // }
}

meshkernel::UInt meshkernel::CurvilinearGridSplineToGrid::longestSplineLength(const Splines& splines) const
{
    UInt result = 0;

    for (UInt i = 0; i < splines.GetNumSplines(); ++i)
    {
        result = std::max(result, static_cast<UInt>(splines.m_splineNodes[i].size()));
    }

    return result;
}

bool meshkernel::CurvilinearGridSplineToGrid::checkSplines(const Splines& splines [[maybe_unused]]) const
{
    return true;
}

void meshkernel::CurvilinearGridSplineToGrid::determineIntersection(Splines& splines,
                                                                    const UInt i,
                                                                    const UInt j,
                                                                    UInt& numberTimesCrossing,
                                                                    double& crossProductOfIntersection,
                                                                    double& firstNormalisedIntersectionLength,
                                                                    double& secondNormalisedIntersectionLength) const
{
    Point intersectionPoint;

    // TODO Need way of getting the number of times the spline intersect.
    if (splines.GetSplinesIntersection(i, j, crossProductOfIntersection, intersectionPoint, firstNormalisedIntersectionLength, secondNormalisedIntersectionLength))
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

void meshkernel::CurvilinearGridSplineToGrid::sectr(Splines& splines,
                                                    EigenMatrix<double>& splineIntersections,
                                                    AnotherMatrix& mn12,
                                                    UInt& numi) const
{
    // TODO need to know if two splines cross more than 1 time

    std::cout << "numebr of splines: " << splines.GetNumSplines() << std::endl;

    std::vector<int> ntyp(splines.GetNumSplines(), 0);
    // std::vector<int> ntyp(longestSplineLength(splines), 0);

    // EigenMatrix<double> splineIntersections(splines.GetNumSplines(), splines.GetNumSplines());

    if (!checkSplines(splines))
    {
        // Some error message
        // Or raise exception
        return;
    }

    bool doubleSupportPoints = false;

label_5:
    // Check the number of splines

    std::fill(splineIntersections.begin(), splineIntersections.end(), std::vector<double>(splines.GetNumSplines(), 0.0));
    // splineIntersections.fill(0.0);
    std::ranges::fill(ntyp, 0);
    ntyp[0] = 1;

    if (doubleSupportPoints)
    {
        // VERDUBBEL AANTAL STEUNPUNTEN ALS
    }

    UInt unNamed = 0;

label_6:

    for (UInt i = 0; i < splines.GetNumSplines(); ++i)
    {
        for (UInt j = i + 1; j < splines.GetNumSplines(); ++j)
        {
            double crp;
            double ti; // firstNormalisedIntersection
            double tj; // secondNormalisedIntersection
            UInt numberTimesCrossing = 0;

            determineIntersection(splines, i, j, numberTimesCrossing, crp, ti, tj); // firstNormalisedIntersection, secondNormalisedIntersection);

            // std::cout << " determine-intersection "
            //           << std::setw(12) << i + 1
            //           << std::setw(12) << j + 1
            //           << std::setw(12) << numberTimesCrossing << "  "
            //           << std::setw(15) << ti << "  "
            //           << std::setw(15) << tj << std::endl;

            if (numberTimesCrossing == 1)
            {
                // std::cout << "ntyp: " << ntyp[i] << "  " << ntyp[j] << std::endl;

                if (ntyp[i] * ntyp[j] == 1)
                {
                    // TODO What to do about this?
                    // if (numpx > nmax / 2)
                    // {
                    //     // QNERROR(' ', ' ', 'Spaghetty; spline both in m- and n-direction');
                    //     // MERR = MERR + 1;
                    //     return;
                    // }
                    // else
                    {
                        doubleSupportPoints = true;
                        std::cout << "first goto label_5" << std::endl;
                        goto label_5;
                    }
                }
                else if (ntyp[i] == 0 && ntyp[j] == 0)
                {
                    std::cout << "both undefined -- " << i << "  " << j << std::endl;
                    // error/warning
                    // both undefined yet.
                }
                else if (ntyp[j] == 0)
                {
                    ntyp[j] = -ntyp[i];

                    if (crp * ntyp[i] < 0)
                    {
                        // std::cout << "reversing spline A " << j << std::endl;
                        // switch
                        splines.Reverse(j);
                        // Reverse the normalised distance
                        tj = static_cast<double>(splines.m_splineNodes[j].size()) - 1.0 - tj;
                    }
                }
                else if (ntyp[i] == 0)
                {
                    ntyp[i] = -ntyp[j];

                    if (crp * ntyp[j] > 0)
                    {
                        std::cout << "reversing spline B " << i << std::endl;
                        // switch
                        splines.Reverse(i);
                        ti = static_cast<double>(splines.m_splineNodes[i].size()) - 1.0 - ti;
                    }
                }

#ifdef USE_EIGEN
                splineIntersections(i, j) = ti;
                splineIntersections(j, i) = tj;
#else
                splineIntersections[i][j] = ti;
                splineIntersections[j][i] = tj;
#endif
            }
            else if (numberTimesCrossing > 1)
            {

                doubleSupportPoints = true;
                std::cout << "second goto label_5" << std::endl;
                goto label_5;
                // // Error: 2 splines (i and j) intersect more than once
                // return;
            }
        }
    }

    // end of label_5 block

    // if (const auto iter = std::ranges::find (ntyp, 0);  iter != ntyp.end ())
    // {
    //     if ()
    // }

    for (UInt i = 0; i < splines.GetNumSplines(); ++i)
    {
        if (ntyp[i] == 0)
        {
            ++unNamed;

            if (unNamed > 1000)
            {
                // Error; one of the splines cannot be attached in the grid.
                return;
            }

            std::cout << "first goto label_6 -- " << unNamed << std::endl;
            goto label_6;
        }
    }

    // for (UInt i = 0; i < splines.GetNumSplines(); ++i)
    // {
    //     for (UInt j = 0; j < splines.GetNumSplines(); ++j)
    //     {
    //         std::cout << " before " << std::setw(15) << i + 1 << "  " << std::setw(15) << j + 1 << "  " << splineIntersections(i, j) << "  ";
    //         std::cout << std::endl;
    //     }

    //     // std::cout << std::endl;
    // }

    // end of label_6 block

    // sorteren op type, eerst de horizontalen (N = CONSTANT)
    for (UInt i = 0; i < splines.GetNumSplines(); ++i)
    {
        if (ntyp[i] == -1)
        {
            for (UInt j = i + 1; j < splines.GetNumSplines(); ++j)
            {
                if (ntyp[j] == 1)
                {
                    splines.SwapSplines(i, j);

                    // auto row = splineIntersections.row(i);
                    // splineIntersections.row(i) = splineIntersections.row(j);
                    // splineIntersections.row(j) = row;

                    // auto col = splineIntersections.col(i);
                    // splineIntersections.col(i) = splineIntersections.col(j);
                    // splineIntersections.col(j) = col;

#ifdef USE_EIGEN
                    splineIntersections.row(i).swap(splineIntersections.row(j));
                    splineIntersections.col(i).swap(splineIntersections.col(j));
#else
                    splineIntersections[i].swap(splineIntersections[j]);
                    SwapColumns(splineIntersections, i, j);
#endif

                    ntyp[i] = 1;
                    ntyp[j] = -1;
                    break;
                }
            }
        }
    }

    // std::cout << std::endl;

    // for (UInt i = 0; i < ntyp.size(); ++i)
    // {
    //     std::cout << "ntyp " << i << " = " << ntyp[i] << std::endl;
    // }

    // std::cout << std::endl;

    // Number of splines in m direction (or is it n)
    numi = constants::missing::uintValue;

    for (UInt i = 0; i < splines.GetNumSplines(); ++i)
    {
        if (ntyp[i] == -1)
        {
            numi = i;
            break;
        }
    }

    //--------------------------------

    std::cout << "numi " << numi << std::endl;

label_59:

    // sort de M
    bool jaChange = false;
    UInt count = 0;

    for (UInt i = 0; i < numi; ++i)
    {
        for (UInt j = numi; j < splines.GetNumSplines(); ++j)
        {
#ifdef USE_EIGEN
            if (splineIntersections(i, j) != 0.0)
#else
            if (splineIntersections[i][j] != 0.0)
#endif
            {
                for (UInt k = j + 1; k < splines.GetNumSplines(); ++k)
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

                            // std::cout << "swapping " << j + 1 << "  " << k + 1 << std::endl;

                            if (count > splines.GetNumSplines())
                            {
                                // error message: PROBLEM IN SPLINE ORDERING, MODIFY SPLINES
                            }

                            goto label_59;
                        }
                    }
                }
            }
        }
    }

    //--------------------------------

label_79:

    // sorteer de N
    count = 0;

    for (UInt i = numi; i < splines.GetNumSplines(); ++i)
    {
        for (UInt j = 0; j < numi; ++j)
        {

#ifdef USE_EIGEN
            if (splineIntersections(i, j) != 0.0)
#else
            if (splineIntersections[i][j] != 0.0)
#endif
            {
                for (UInt k = j + 1; k < numi; ++k)
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
                                // error message: PROBLEM IN SPLINE ORDERING, MODIFY SPLINES
                            }

                            goto label_79;
                        }
                    }
                }
            }
        }
    }

    if (jaChange)
    {
        goto label_59;
    }

    //--------------------------------
    // Initialiseer ranking, start en eind, 1,2,3

    // mn12.fill(std::{0, 0, 0});

    std::ranges::fill(mn12, std::array<int, 3>({0, 0, 0}));
    int maxn = 0;
    int maxm = 0;

    // Eerst alles ranken in N richting
    for (UInt i = 0; i < numi; ++i)
    {
        for (UInt j = numi; j < splines.GetNumSplines(); ++j)
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
                    maxn = mn12[jjlast][0] + 1;
                    jjlast = k;
                }
            }

            mn12[j][1] = maxn;
        }

        maxn = 0;

        for (UInt j = numi; j < splines.GetNumSplines(); ++j)
        {
#ifdef USE_EIGEN
            if (splineIntersections(j, i) != 0.0)
#else
            if (splineIntersections[j][i] != 0.0)
#endif
            {
                maxn = std::max(maxn, mn12[j][1]);
            }
        }

        mn12[i][0] = maxn;
    }

    // Dan alles ranken in M richting
    for (UInt i = numi; i < splines.GetNumSplines(); ++i)
    {
        for (UInt j = 0; j < numi; ++j)
        {
            maxm = 0;
            UInt iilast = numi;

            for (UInt k = numi; k <= i; ++k)
            {
#ifdef USE_EIGEN
                if (splineIntersections(j, k) != 0.0)
#else
                if (splineIntersections[j][k] != 0.0)
#endif
                {
                    maxm = mn12[iilast][0] + 1;
                    iilast = k;
                }
            }

            mn12[j][2] = maxm;
        }

        maxm = 0;

        for (UInt j = 0; j < numi; ++j)
        {
#ifdef USE_EIGEN
            if (splineIntersections(j, i) != 0.0)
#else
            if (splineIntersections[j][i] != 0.0)
#endif
            {
                maxm = std::max(maxm, mn12[j][2]);
            }
        }

        mn12[i][0] = maxm;
    }

    // mn12.col(1).fill(0);
    // mn12.col(2).fill(0);

    for (UInt i = 0; i < splines.GetNumSplines(); ++i)
    {
        mn12[i][1] = 0;
        mn12[i][2] = 0;
    }

    // Daarna per spline begin- en eindpunt tellen, eerst N = constant

    for (UInt i = 0; i < numi; ++i)
    {
        for (UInt j = numi; j < splines.GetNumSplines(); ++j)
        {
#ifdef USE_EIGEN
            if (splineIntersections(i, j) != 0.0)
#else
            if (splineIntersections[i][j] != 0.0)
#endif
            {
                if (mn12[i][1] == 0)
                {
                    mn12[i][1] = mn12[j][0];
                }

                mn12[i][2] = mn12[j][0];
            }
        }
    }

    // Dan M = constant

    for (UInt i = numi; i < splines.GetNumSplines(); ++i)
    {
        for (UInt j = 0; j < numi; ++j)
        {
#ifdef USE_EIGEN
            if (splineIntersections(i, j) != 0.0)
#else
            if (splineIntersections[i][j] != 0.0)
#endif
            {
                if (mn12[i][1] == 0)
                {
                    mn12[i][1] = mn12[j][0];
                }

                mn12[i][2] = mn12[j][0];
            }
        }
    }

    // for (UInt i = 0; i < splines.GetNumSplines(); ++i)
    // {
    //     std::cout << "mn12 " << std::setw(5) << mn12(i, 0) << "  " << std::setw(5) << mn12[i][ 1] << "  " << std::setw(5) << mn12[i][ 2] << "  " << std::endl;
    // }
}

std::vector<double> meshkernel::CurvilinearGridSplineToGrid::paktij(const EigenMatrix<double>& splineIntersections,
                                                                    const UInt whichRow) const
{
    std::vector<double> compressedRow;

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

    // if (tValue > static_cast<double>(splines.m_splineNodes[whichSpline].size() - 1))
    // {
    //     std::cout << "tValue = " << tValue << "  " << splines.m_splineNodes[whichSpline].size() << std::endl;
    // }

    // in getdis.f90 this is: tValue = std::min(tValue, static_cast<double>(splines.m_splineNodes[whichSpline].size()))
    // without the -1. I think the fortran is incorrect, there should be the -1
    tValue = std::min(tValue, static_cast<double>(splines.m_splineNodes[whichSpline].size() - 1));
    sValue = 0.0;

    do
    {
        t1 = t0 + dt;

        if (t1 < tValue)
        {
            endPoint = splines.Evaluate(whichSpline, t1);
            // std::cout << " end point A " << endPoint.x << "  " << endPoint.y << std::endl;
        }
        else
        {
            endPoint = splines.Evaluate(whichSpline, tValue - 1.0e-5);
            // std::cout << " end point B " << endPoint.x << "  " << endPoint.y << "  " << tValue << "  " << splines.m_splineNodes[whichSpline].size() << std::endl;
        }

        sValue += ComputeDistance(startPoint, endPoint, splines.m_projection);
        startPoint = endPoint;
        t0 = t1;

    } while (t1 < tValue);

    // std::cout << " sValue " << tValue << "  " << sValue << std::endl;
}

void meshkernel::CurvilinearGridSplineToGrid::makesr(const double ar,
                                                     const double s0,
                                                     const double s1,
                                                     std::vector<double>& sr) const
{
    double ds = 1.0;
    sr[0] = 0.0;

    // std::cout << " sr " << ar << " :";

    for (UInt k = 0; k < sr.size() - 1; ++k)
    {
        sr[k + 1] = sr[k] + ds;

        // std::cout << " {" << sr[k] << "  " << sr[k + 1] << "  " << ds << "}";

        ds *= ar;
    }

    // std::cout << std::endl;

    [[maybe_unused]] double fac = (s1 - s0) / sr[sr.size() - 1];

    for (UInt k = 0; k < sr.size(); ++k)
    {
        sr[k] = s0 + fac * sr[k];
    }
}

void meshkernel::CurvilinearGridSplineToGrid::makessq(const std::vector<double>& s, // intersectionPoints
                                                      const UInt mFac,
                                                      std::vector<double>& ssq) const
{
    if (s.size() == 2)
    {
        for (UInt i = 0; i <= mFac; ++i)
        {
            ssq[i] = s[0] + (s[1] - s[0]) * static_cast<double>(i) / static_cast<double>(mFac);
        }
    }
    else if (s.size() >= 3)
    {
        std::vector<double> a(s.size());
        std::vector<double> sl(mFac + 1);
        std::vector<double> sr(mFac + 1);

        // std::vector<double> glad (intersectionPoints.size ());

        // for (UInt i = 1; i < intersectionPoints.size () - 1; ++i)
        // {
        //     glad[i] = (s[i+1] - s[i]) / (s[i] - s[i - 1]);
        // }

        // std::cout << "a : ";

        for (UInt i = 1; i < s.size() - 1; ++i)
        {
            a[i] = (s[i + 1] - s[i]) / (s[i] - s[i - 1]);
            // std::cout << " {" << i << ", " << a[i] << "  " << s[i + 1] << "  " << s[i] << "  " << s[i - 1] << "}, ";
        }

        // std::cout << std::endl;

        a[0] = a[1];
        a[s.size() - 1] = a[s.size() - 2];

        // std::cout << " ssq = ";
        // UInt kr = 0;

        for (UInt i = 0; i < s.size() - 1; ++i)
        {
            // std::cout << a[i + 1] << " --- ";
            double ar = std::pow(a[i + 1], 1.0 / static_cast<double>(mFac));
            makesr(ar, s[i], s[i + 1], sr);

            double al = std::pow(a[i], 1.0 / static_cast<double>(mFac));
            makesr(al, s[i], s[i + 1], sl);

            for (UInt k = 1; k <= mFac + 1; ++k)
            {
                UInt kr = i * mFac + k - 1;
                // std::cout << kr << "  ";
                ar = static_cast<double>(k - 1) / static_cast<double>(mFac);
                al = 1.0 - ar;
                ssq[kr] = ar * sr[k - 1] + al * sl[k - 1];

                if (!std::isfinite(ssq[kr]))
                {
                    std::cout << "{ " << ar << "  " << al << "  " << k << "  " << ar << "  " << sr[k - 1] << "  " << al << "  " << sl[k - 1] << "} "; // << std::endl;
                }

                // std::cout << "{ "  << ssq[kr];
                ar = (ssq[kr] - s[i]) / (s[i + 1] - s[i]);
                al = 1.0 - ar;
                ssq[kr] = ar * sr[k - 1] + al * sl[k - 1];
                // ++kr;
            }
        }

        // std::cout << std::endl;
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
    // double bx = 0.5 * (ax + cx);
    // double tolerance = 0.00001;

    FuncAdimensionalToDimensionalDistanceOnSpline func(splines, whichSpline, false, 0.0);

    func.SetDimensionalDistance(ssq);

    double distance = FindFunctionRootWithGoldenSectionSearch(func, ax, cx);

    result = splines.Evaluate(whichSpline, distance);

    return result;
}

void meshkernel::CurvilinearGridSplineToGrid::makespl(const Splines& splines,
                                                      const UInt whichSpline,
                                                      const UInt mnFac,
                                                      std::vector<double>& intersectionPoints,
                                                      std::vector<Point>& gridPoints) const
{
    // evaluate the spline at the intersection points
    bool curvatureAdapted = false;                              // Parameter H is not used in getdis.f90
    double maximumGridHeight = constants::missing::doubleValue; // what values should this be?

    auto [splinePoints, distances] = splines.ComputePointOnSplineFromAdimensionalDistance(whichSpline, maximumGridHeight, curvatureAdapted, intersectionPoints);

    // std::cout << "makspl ";

    // for (UInt i = 0; i < distances.size(); ++i)
    // {
    //     std::cout << "{" << splinePoints[i].x << ", " << splinePoints[i].y << ", " << distances[i] << ", " << intersectionPoints[i] << "}, ";
    // }

    // std::cout << std::endl;

    //--------------------------------
    // getdis

    std::vector<double> s(intersectionPoints.size());

    for (UInt i = 0; i < intersectionPoints.size(); ++i)
    {
        getdis(splines, whichSpline, intersectionPoints[i], s[i]);
    }

    UInt kmax = (intersectionPoints.size() - 1) * mnFac + 1;

    std::vector<double> ssq(kmax, -999.0);

    if (intersectionPoints.size() >= 2)
    {
        makessq(s, mnFac, ssq);

        for (UInt l = 1; l <= intersectionPoints.size() - 1; ++l)
        {
            // TODO check the +1 here
            UInt k1 = mnFac * (l - 1); // + 1;
            UInt k2 = k1 + mnFac;
            bool jaDip = false;

        label_23:

            if (jaDip)
            {
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
        Point p = GetXy(splines, whichSpline, intersectionPoints, ssq[k]); // splines.ComputeClosestPointOnSplineSegment(whichSpline, 0.0, splineMax, )
        gridPoints[k] = p;
    }
}

void meshkernel::CurvilinearGridSplineToGrid::splrgf(Splines& splines,
                                                     const EigenMatrix<double>& splineIntersections,
                                                     const AnotherMatrix& mn12,
                                                     CurvilinearGrid& grid [[maybe_unused]],
                                                     const UInt numi) const
{
    UInt nFac = 5;
    UInt mFac = 5;
    // UInt nFac = 10;
    // UInt mFac = 10;
    std::vector<Point> gridPoints;
    // EigenMatrix<Point> gridNodes(612, 612);
    UInt cols = (numi - 1) * mFac;
    UInt rows = (splines.GetNumSplines() - numi) * nFac;

    std::cout << "estimated sizeA: " << rows << " " << cols << std::endl;
    std::cout << "estimated sizeB: " << (splines.GetNumSplines() - numi) * nFac << "  " << (numi - 1) * mFac << std::endl;

    UInt maxSize = 0;
    UInt ns = 0;
    UInt ms = 0;
    UInt ns2 = 0;
    UInt ms2 = 0;

    cols = 0;
    rows = 0;

    for (UInt i1 = 0; i1 < splines.GetNumSplines(); ++i1)
    {
        if (i1 < numi)
        {
            cols = std::max(cols, static_cast<UInt>(mn12[i1][0]));
        }
        else
        {
            rows = std::max(rows, static_cast<UInt>(mn12[i1][0]));
        }
    }

    std::cout << "estimated sizeC: " << rows * nFac << " " << cols * mFac << std::endl;
    rows *= nFac;
    cols *= mFac;

    // EigenMatrix<Point> gridNodes(612, 612);
    // lin_alg::Matrix<Point> gridNodes(11, 11);
    lin_alg::Matrix<Point> gridNodes(rows, cols);
    // lin_alg::Matrix<Point> gridNodes(320, 320);
    // lin_alg::Matrix<Point> gridNodes(611, 611);
    gridNodes.fill({constants::missing::doubleValue, constants::missing::doubleValue});

    // TOOD can this be changed to just 'i'
    for (UInt i1 = 0; i1 < splines.GetNumSplines(); ++i1)
    {
        std::vector<double> intersectionPoints = paktij(splineIntersections, i1);
        // std::vector<Point> pointsAlongSpline;
        // std::vector<double> normalisedDistances;
        [[maybe_unused]] UInt jj;
        [[maybe_unused]] UInt ii1;
        [[maybe_unused]] UInt ii2;

        // TODO Check that "+ 1" in the statements below, might be that it should be removed.
        if (i1 < numi)
        {
            // std::tie(pointsAlongSpline, normalisedDistances) = splines.ComputePointOnSplineFromAdimensionalDistance(i1, );
            makespl(splines, i1, mFac, intersectionPoints, gridPoints);
            jj = (mn12[i1][0] - 1) * nFac;
            ii1 = (mn12[i1][1] - 1) * mFac;
            ii2 = (mn12[i1][2] - 1) * mFac + 1;
        }
        else
        {
            // std::tie(pointsAlongSpline, normalisedDistances) = splines.ComputePointOnSplineFromAdimensionalDistance(i1, );
            makespl(splines, i1, nFac, intersectionPoints, gridPoints);
            jj = (mn12[i1][0] - 1) * mFac;
            ii1 = (mn12[i1][1] - 1) * nFac;
            ii2 = (mn12[i1][2] - 1) * nFac + 1;
        }

        maxSize = std::max(maxSize, jj);
        maxSize = std::max(maxSize, ii1);
        maxSize = std::max(maxSize, ii2);

        UInt k = 0;

        // TODO change loop variable name, too confusing ii, i1, ii1.
        for (UInt ii = ii1; ii < ii2; ++ii)
        {

            if (k < gridPoints.size())
            {
                if (i1 < numi)
                {
                    gridNodes(ii, jj) = gridPoints[k];

                    if (gridPoints[k].IsValid())
                    {
                        ns2 = std::max(ns2, ii);
                        ms2 = std::max(ms2, jj);
                    }
                }
                else
                {
                    gridNodes(jj, ii) = gridPoints[k];

                    if (gridPoints[k].IsValid())
                    {
                        ns2 = std::max(ns2, jj);
                        ms2 = std::max(ms2, ii);
                    }
                }
            }

            ++k;
        }

        std::cout << "i1: " << i1 << "  " << mn12[i1][0] << "  " << ns << "  " << ms << std::endl;

        if (i1 < numi)
        {
            ns = std::max(ns, static_cast<UInt>(mn12[i1][0]));
        }
        else
        {
            ms = std::max(ms, static_cast<UInt>(mn12[i1][0]));
        }
    }

    std::cout << "grid size: " << ms << " x " << ns << "  " << mFac * ms << " x " << nFac * ns << " --- " << ns2 << "  " << ms2 << std::endl;

    UInt ncr = (ns - 1) * nFac; // + 1;
    UInt mcr = (ms - 1) * mFac; // + 1;

    UInt nmax = std::max(ncr, mcr) + 1;

    std::vector<Point> x1(nmax);
    std::vector<Point> x2(nmax);
    std::vector<Point> x3(nmax);
    std::vector<Point> x4(nmax);

    for (UInt i = 0; i < ms - 1; ++i)
    {
        for (UInt j = 0; j < ns - 1; ++j)
        {
            x1[1].SetInvalid();
            x2[1].SetInvalid();
            x3[1].SetInvalid();
            x4[1].SetInvalid();

            // TOOD check k and l loop range
            for (UInt k = 0; k <= mFac; ++k)
            {
                for (UInt l = 0; l <= nFac; ++l)
                {
                    UInt ki = i * mFac + k;
                    UInt lj = j * mFac + l;
                    Point gridNode = gridNodes(ki, lj);

                    // if (!gridNode.IsValid())
                    // {
                    //     continue;
                    // }

                    if (gridNode.IsValid())
                    {
                        if (k == 0)
                        {
                            x1[l] = gridNode;
                        }
                        else if (k == mFac)
                        {
                            x2[l] = gridNode;
                        }
                        if (l == 0)
                        {
                            x3[k] = gridNode;
                        }
                        else if (l == nFac)
                        {
                            x4[k] = gridNode;
                        }
                    }
                }
            }

            // std::cout << "x1 -- ";

            // for (UInt kk = 0; kk < x1.size(); ++kk)
            // {
            //     std::cout << "{" << x1[kk].x << ", " << x1[kk].y << "} ";
            // }

            // std::cout << std::endl;
            // std::cout << "x2 -- ";

            // for (UInt kk = 0; kk < x2.size(); ++kk)
            // {
            //     std::cout << "{" << x2[kk].x << ", " << x2[kk].y << "} ";
            // }

            // std::cout << std::endl;
            // std::cout << "x3 -- ";

            // for (UInt kk = 0; kk < x3.size(); ++kk)
            // {
            //     std::cout << "{" << x3[kk].x << ", " << x3[kk].y << "} ";
            // }

            // std::cout << std::endl;
            // std::cout << "x4 -- ";

            // for (UInt kk = 0; kk < x4.size(); ++kk)
            // {
            //     std::cout << "{" << x4[kk].x << ", " << x4[kk].y << "} ";
            // }

            // std::cout << std::endl;

            UInt no = 0;

            if (!x1[1].IsValid())
            {
                no = 1;
            }
            if (!x2[1].IsValid())
            {
                no = 1;
            }
            if (!x3[1].IsValid())
            {
                no = 1;
            }
            if (!x4[1].IsValid())
            {
                no = 1;
            }

            if (no == 0)
            {
                auto interpolationResult = DiscretizeTransfinite(x1,
                                                                 x2,
                                                                 x3,
                                                                 x4,
                                                                 splines.m_projection,
                                                                 mFac,
                                                                 nFac);

                for (UInt k2 = 0; k2 <= mFac; ++k2)
                {
                    for (UInt l2 = 1; l2 <= nFac; ++l2)
                    {
                        UInt ki = i * mFac + k2;
                        UInt lj = j * nFac + l2;

                        if (!gridNodes(ki, lj).IsValid() && std::isfinite(interpolationResult(k2, l2).x) && std::isfinite(interpolationResult(k2, l2).y))
                        {
                            gridNodes(ki, lj) = interpolationResult(k2, l2);
                        }
                    }
                }
            }
        }
    }

    grid.SetGridNodes(gridNodes);

    std::cout << " maxSize " << maxSize << std::endl;
}
