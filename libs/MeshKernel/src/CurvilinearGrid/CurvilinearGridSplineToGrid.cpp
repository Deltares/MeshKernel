#include "MeshKernel/CurvilinearGrid/CurvilinearGridSplineToGrid.hpp"
#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Utilities/LinearAlgebra.hpp"

#include <algorithm>
#include <iomanip>

void meshkernel::CurvilinearGridSplineToGrid::Compute(const Splines& splines, CurvilinearGrid& grid [[maybe_unused]]) const
{
    Splines splinesCopy(splines);
    lin_alg::Matrix<double> splineIntersections(splines.GetNumSplines(), splines.GetNumSplines());
    std::cout << "Compute:: number of splines: " << splines.GetNumSplines() << std::endl;

    lin_alg::Matrix<int> mn12(splines.GetNumSplines(), 3);

    sectr(splinesCopy, splineIntersections, mn12);
    splrgf(splinesCopy, splineIntersections, mn12, grid);
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
                                                    lin_alg::Matrix<double>& splineIntersections,
                                                    lin_alg::Matrix<int>& mn12) const
{
    // TODO need to know if two splines cross more than 1 time

    std::cout << "numebr of splines: " << splines.GetNumSplines() << std::endl;

    std::vector<int> ntyp(longestSplineLength(splines), 0);

    // lin_alg::Matrix<double> splineIntersections(splines.GetNumSplines(), splines.GetNumSplines());

    if (!checkSplines(splines))
    {
        // Some error message
        // Or raise exception
        return;
    }

    bool doubleSupportPoints = false;

label_5:
    // Check the number of splines

    splineIntersections.fill(0.0);
    std::ranges::fill(ntyp, 0);
    ntyp[0] = 1;

    if (doubleSupportPoints)
    {
        // VERDUBBEL AANTAL STEUNPUNTEN ALS
    }

    UInt unNamed = 0;

label_6:
    std::cout.precision(15);

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

                splineIntersections(i, j) = ti;
                splineIntersections(j, i) = tj;
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

                    auto row = splineIntersections.row(i);
                    splineIntersections.row(i) = splineIntersections.row(j);
                    splineIntersections.row(j) = row;

                    auto col = splineIntersections.col(i);
                    splineIntersections.col(i) = splineIntersections.col(j);
                    splineIntersections.col(j) = col;

                    // splineIntersections.row(i).swap(splineIntersections.row(j));
                    // splineIntersections.col(i).swap(splineIntersections.col(j));
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

    UInt numi = constants::missing::uintValue;

    for (UInt i = 0; i < splines.GetNumSplines(); ++i)
    {
        if (ntyp[i] == 1)
        {
            numi = i;
        }
    }

    //--------------------------------
    std::cout << "numi " << numi << std::endl;

label_59:

    // sort de M
    bool jaChange = false;
    UInt count = 0;

    for (UInt i = 0; i <= numi; ++i)
    {

        for (UInt j = numi + 1; j < splines.GetNumSplines(); ++j)
        {
            if (splineIntersections(i, j) != 0.0)
            {
                for (UInt k = j + 1; k < splines.GetNumSplines(); ++k)
                {
                    if (splineIntersections(i, k) != 0.0)
                    {

                        if (splineIntersections(i, j) > splineIntersections(i, k))
                        {
                            splines.SwapSplines(j, k);
                            splineIntersections.row(j).swap(splineIntersections.row(k));
                            splineIntersections.col(j).swap(splineIntersections.col(k));
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

    for (UInt i = numi + 1; i < splines.GetNumSplines(); ++i)
    {

        for (UInt j = 0; j <= numi; ++j)
        {
            if (splineIntersections(i, j) != 0.0)
            {
                for (UInt k = j + 1; k <= numi; ++k)
                {
                    if (splineIntersections(i, k) != 0.0)
                    {

                        if (splineIntersections(i, j) > splineIntersections(i, k))
                        {
                            splines.SwapSplines(j, k);
                            splineIntersections.row(j).swap(splineIntersections.row(k));
                            splineIntersections.col(j).swap(splineIntersections.col(k));
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

    mn12.fill(0);
    int maxn = 0;
    int maxm = 0;

    // Eerst alles ranken in N richting
    for (UInt i = 0; i <= numi; ++i)
    {
        for (UInt j = numi + 1; j < splines.GetNumSplines(); ++j)
        {
            maxn = 0;
            UInt jjlast = 0;

            for (UInt k = 0; k <= i; ++k)
            {
                if (splineIntersections(j, k) != 0.0)
                {
                    maxn = mn12(jjlast, 0) + 1;
                    jjlast = k;
                }

                mn12(j, 1) = maxn;
            }
        }

        maxn = 0;

        for (UInt j = numi + 1; j < splines.GetNumSplines(); ++j)
        {
            if (splineIntersections(j, i) != 0.0)
            {
                maxn = std::max(maxn, mn12(j, 1));
            }
        }

        mn12(i, 0) = maxn;
    }

    // Dan alles ranken in M richting
    for (UInt i = numi + 1; i < splines.GetNumSplines(); ++i)
    {
        for (UInt j = 0; j <= numi; ++j)
        {
            maxm = 0;
            UInt iilast = 0;

            for (UInt k = numi + 1; k <= i; ++k)
            {
                if (splineIntersections(j, k) != 0.0)
                {
                    maxm = mn12(iilast, 0) + 1;
                    iilast = k;
                }

                mn12(j, 2) = maxm;
            }
        }

        maxm = 0;

        for (UInt j = 0; j <= numi; ++j)
        {
            if (splineIntersections(j, i) != 0.0)
            {
                maxm = std::max(maxm, mn12(j, 2));
            }
        }

        mn12(i, 0) = maxm;
    }

    // mn12.col(1).fill(0);
    // mn12.col(2).fill(0);

    for (UInt i = 0; i < splines.GetNumSplines(); ++i)
    {
        mn12(i, 1) = 0;
        mn12(i, 2) = 0;
    }

    // Daarna per spline begin- en eindpunt tellen, eerst N = constant

    for (UInt i = 0; i <= numi; ++i)
    {
        for (UInt j = numi + 1; j < splines.GetNumSplines(); ++j)
        {
            if (splineIntersections(i, j) != 0.0)
            {
                if (mn12(i, 1) == 0)
                {
                    mn12(i, 1) = mn12(j, 0);
                }

                mn12(i, 2) = mn12(j, 0);
            }
        }
    }

    // Dan M = constant

    for (UInt i = numi + 1; i < splines.GetNumSplines(); ++i)
    {
        for (UInt j = 0; j <= numi; ++j)
        {
            if (splineIntersections(i, j) != 0.0)
            {
                if (mn12(i, 1) == 0)
                {
                    mn12(i, 1) = mn12(j, 0);
                }

                mn12(i, 2) = mn12(j, 0);
            }
        }
    }

    // for (UInt i = 0; i < splines.GetNumSplines(); ++i)
    // {
    //     std::cout << "mn12 " << std::setw(5) << mn12(i, 0) << "  " << std::setw(5) << mn12(i, 1) << "  " << std::setw(5) << mn12(i, 2) << "  " << std::endl;
    // }
}

std::vector<double> meshkernel::CurvilinearGridSplineToGrid::paktij(const lin_alg::Matrix<double>& splineIntersections,
                                                                    const UInt whichRow) const
{
    std::vector<double> compressedRow;

    for (UInt col = 0; col < splineIntersections.cols(); ++col)
    {
        if (splineIntersections(whichRow, col) != 0.0)
        {
            compressedRow.push_back(splineIntersections(whichRow, col));
        }
    }

    return compressedRow;
}

void meshkernel::CurvilinearGridSplineToGrid::getdis(const Splines& splines,
                                                     const UInt whichSpline,
                                                     double& tValue,
                                                     double& sValue) const
{
    double ts = std::min(ts, stattic_cast<double>(splines.m_splineNodes[whichSpline].size()));
    double dt = 0.1;
    double xt0 = splines.m_splineNodes[whichSpline].x;
    double yt0 = splines.m_splineNodes[whichSpline].y;
    sValue = 0.0;
    tValue = 0.0;

    do
    {
        t1 = t0 + dt;

        if (t1 < ts)
        {
        }
        else
        {
        }

        ss += dbdistance();

    } while (t1 < ts);
}

void meshkernel::CurvilinearGridSplineToGrid::makespl(const Splines& splines,
                                                      const UInt whichSpline,
                                                      const std::vector<double>& intersectionPoints) const
{
    // evaluate the spline at the intersection points
    bool curvatureAdapted = false;                              // Parameter H is not used in getdis.f90
    double maximumGridHeight = constants::missing::doubleValue; // what values should this be?

    auto [splinePoints, distances] = splines.ComputePointOnSplineFromAdimensionalDistance(whichSpline, maximumGridHeight, curvatureAdapted, intersectionPoints);

    std::cout << "makspl ";

    for (UInt i = 0; i < distances.size(); ++i)
    {
        std::cout << "{" << splinePoints[i].x << ", " << splinePoints[i].y << ", " << distances[i] << "}, ";
    }

    std::cout << std::endl;

    //--------------------------------
    // getdis
    for (UInt i = 0; i < intersectionPoints.size(); ++i)
    {
        getdis(splines, whichSpline, splinePoints, t[i], s[i]);
    }
}

void meshkernel::CurvilinearGridSplineToGrid::splrgf(Splines& splines,
                                                     const lin_alg::Matrix<double>& splineIntersections,
                                                     const lin_alg::Matrix<int>& mn12,
                                                     CurvilinearGrid& grid [[maybe_unused]]) const
{
    UInt nfac = 10;
    UInt mfac = 10;

    // TOOD can this be changed to just 'i'
    for (UInt i1 = 0; i1 < splines.GetNumSplines(); ++i1)
    {
        std::vector<double> intersectionPoints = paktij(splineIntersections, i1);
        // std::vector<Point> pointsAlongSpline;
        // std::vector<double> normalisedDistances;
        [[maybe_unused]] UInt jj;
        [[maybe_unused]] UInt ii1;
        [[maybe_unused]] UInt ii2;

        // std::cout << "spline " << i1 << "  " << splines.m_splineNodes.size() << "  " << row.size() << " : -- ";

        // for (UInt i = 0; i < row.size(); ++i)
        // {
        //     std::cout << row[i] << "  ";
        // }

        // std::cout << std::endl;

        // TODO Check that "+ 1" in the statements below, might be that it should be removed.
        if (i1 + 1 <= intersectionPoints.size())
        {
            // std::tie(pointsAlongSpline, normalisedDistances) = splines.ComputePointOnSplineFromAdimensionalDistance(i1, );
            makespl(splines, i1, intersectionPoints);
            jj = (mn12(i1, 0) - 1) * nfac + 1;
            ii1 = (mn12(i1, 1) - 1) * mfac + 1;
            ii2 = (mn12(i1, 2) - 1) * mfac + 1;
        }
        else
        {
            // std::tie(pointsAlongSpline, normalisedDistances) = splines.ComputePointOnSplineFromAdimensionalDistance(i1, );
            makespl(splines, i1, intersectionPoints);
            jj = (mn12(i1, 0) - 1) * mfac + 1;
            ii1 = (mn12(i1, 1) - 1) * nfac + 1;
            ii2 = (mn12(i1, 2) - 1) * nfac + 1;
        }

        // UInt k = 0;

        // for (UInt ii = ii1; ii < ii2; ++ii)
        // {
        //     ++k;

        //     if (k <= l1max)
        //     {
        //         gridNodes(ii, jj) = Point(x1(k), y1(k));
        //     }
        //     else
        //     {
        //         gridNodes(jj, ii) = Point(x1(k), y1(k));
        //     }
        // }

        // if (i1 <= numi)
        // {
        //     ns = std::max(ns, mn12(i1, 0));
        // }
        // else
        // {
        //     ms = std::max(ms, mn12(i1, 0));
        // }
    }

    // UInt ncr = (ns - 1) * nfac + 1;
    // UInt mcr = (ms - 1) * mfac + 1;

    // std::vector<double> x1(nmax);
    // std::vector<double> x2(nmax);
    // std::vector<double> x3(nmax);
    // std::vector<double> x4(nmax);
    // std::vector<double> y1(nmax);
    // std::vector<double> y2(nmax);
    // std::vector<double> y3(nmax);
    // std::vector<double> y4(nmax);

    // for (UInt i = 0; i < ms - 1; ++i)
    // {
    //     for (UInt j = 0; j < ns - 1; ++j)
    //     {
    //         x1[1] = constants::missing::doubleValue;
    //         x2[1] = constants::missing::doubleValue;
    //         x3[1] = constants::missing::doubleValue;
    //         x4[1] = constants::missing::doubleValue;

    //         // TOOD check k and l loop range
    //         for (UInt k = 1; k < mfac + 1; ++k)
    //         {
    //             for (UInt l = 1; l < nfac + 1; ++l)
    //             {
    //                 ki = (i - 1) * mfac + k;
    //                 lj = (j - 1) * mfac + l;

    //                 if (xc(ki, lj) != constants::missing::doubleValue)
    //                 {

    //                     // TODO if loop k range changes this should change too
    //                     if (k == 1)
    //                     {
    //                         x1[l] = xc(ki, lj);
    //                         y1[l] = yc(ki, lj);
    //                     }
    //                     if (k == mfac + 1)
    //                     {
    //                         x2[l] = xc(ki, lj);
    //                         y2[l] = yc(ki, lj);
    //                     }
    //                     if (l == 1)
    //                     {
    //                         x3[k] = xc(ki, lj);
    //                         y3[k] = yc(ki, lj);
    //                     }
    //                     if (l == nfac + 1)
    //                     {
    //                         x4[k] = xc(ki, lj);
    //                         y4[k] = yc(ki, lj);
    //                     }
    //                 }
    //             }

    //         } // end k loop

    //         UInt no = 0;

    //         if (x1(1) == constants::missing::doubleValue)
    //         {
    //             no = 1;
    //         }

    //         if (x2(1) == constants::missing::doubleValue)
    //         {
    //             no = 1;
    //         }

    //         if (x3(1) == constants::missing::doubleValue)
    //         {
    //             no = 1;
    //         }

    //         if (x4(1) == constants::missing::doubleValue)
    //         {
    //             no = 1;
    //         }

    //         if (no == 0)
    //         {
    //             tranfn2();

    //             for (UInt k = 1; k <= mfac + 1; ++k)
    //             {
    //                 for (UInt l = 1; l <= nfac + 1; ++l)
    //                 {
    //                     ki = (i - 1) * mfac + k;
    //                     lj = (j - 1) * nfac + l;

    //                     if (xc(ki, lj) == constants::missing::doubleValue)
    //                     {
    //                         xs(ki, lj) = xh(k, l);
    //                         ys(ki, lj) = yh(k, l);
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }
}
