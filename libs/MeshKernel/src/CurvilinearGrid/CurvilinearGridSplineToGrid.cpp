#include "MeshKernel/CurvilinearGrid/CurvilinearGridSplineToGrid.hpp"
#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Utilities/LinearAlgebra.hpp"

#include <algorithm>

void meshkernel::CurvilinearGridSplineToGrid::Compute(const Splines& splines, CurvilinearGrid& grid [[maybe_unused]]) const
{
    Splines splinesCopy(splines);
    lin_alg::Matrix<double> splineIntersections(splines.GetNumSplines(), splines.GetNumSplines());

    sectr(splinesCopy, splineIntersections);

    for (UInt i = 0; i < splines.GetNumSplines(); ++i)
    {
        for (UInt j = 0; j < splines.GetNumSplines(); ++j)
        {
            std::cout << splineIntersections(i, j) << "  ";
        }

        std::cout << std::endl;
    }
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
        numberTimesCrossing = 0;
    }
}

void meshkernel::CurvilinearGridSplineToGrid::sectr(Splines& splines, lin_alg::Matrix<double>& splineIntersections) const
{
    // TODO need to know if two splines cross more than 1 time

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

    for (UInt i = 0; i < splines.GetNumSplines(); ++i)
    {
        for (UInt j = i + 1; j < splines.GetNumSplines(); ++j)
        {
            double crp;
            double ti; // firstNormalisedIntersection
            double tj; // secondNormalisedIntersection
            UInt numberTimesCrossing = 0;

            determineIntersection(splines, i, j, numberTimesCrossing, crp, ti, tj); // firstNormalisedIntersection, secondNormalisedIntersection);

            std::cout << " determineIntersection " << numberTimesCrossing << "  " << crp << "  " << ti << "  " << tj << std::endl;

            if (numberTimesCrossing == 1)
            {

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
                    // error/warning
                    // both undefined yet.
                }
                else if (ntyp[j] == 0)
                {
                    ntyp[j] = -ntyp[i];

                    if (crp * ntyp[i] < 0)
                    {
                        // switch
                        splines.Reverse(j);
                        tj = static_cast<double>(splines.m_splineNodes[j].size()) - 1.0 - tj;
                    }
                }
                else if (ntyp[j] == 0)
                {
                    ntyp[i] = -ntyp[j];

                    if (crp * ntyp[j] < 0)
                    {
                        // switch
                        splines.Reverse(i);
                        ti = static_cast<double>(splines.m_splineNodes[i].size()) - 1.0 - ti;
                    }
                }
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

            std::cout << "first goto label_6" << std::endl;
            goto label_6;
        }
    }

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
                    splineIntersections.row(i).swap(splineIntersections.row(j));
                    splineIntersections.col(i).swap(splineIntersections.col(j));
                    ntyp[i] = 1;
                    ntyp[j] = -1;
                }
            }
        }
    }

    UInt numi = constants::missing::uintValue;

    for (UInt i = 0; i < splines.GetNumSplines(); ++i)
    {
        if (ntyp[i] == 1)
        {
            numi = i;
        }
    }

    //--------------------------------

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

    // sort de N
    count = 0;

    for (UInt i = numi + 1; i < splines.GetNumSplines(); ++i)
    {

        for (UInt j = 0; j < numi; ++j)
        {
            if (splineIntersections(i, j) != 0.0)
            {
                for (UInt k = j + 1; k < numi; ++k)
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

    lin_alg::Matrix<int> mn12(splines.GetNumSplines(), 3);
    mn12.fill(0);
    int maxn = 0;
    int maxm = 0;

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
                maxm = std::max(maxm, mn12(j, 3));
            }
        }

        mn12(i, 0) = maxm;
    }

    mn12.col(1).fill(0);
    mn12.col(2).fill(0);

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

                mn12(i, 2) = mn12(j, 1);
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

                mn12(i, 2) = mn12(j, 1);
            }
        }
    }
}
