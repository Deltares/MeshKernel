// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Stichting Deltares
//
// C++20 translation of mshcopyboun.for and mshchkstapl.for from the SEPRAN library.
//
// Original authors:
//   mshcopyboun  – Guus Segal, 1999-2001
//   mshchkstapl  – Guus Segal, 2003

#include "SepranBoundary.hpp"
#include "SepranSort.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

namespace sepran
{

// ---------------------------------------------------------------------------
// copyBoundary  (mshcopyboun.for)
// ---------------------------------------------------------------------------

BoundaryCopyResult copyBoundary(int nbound,
                                std::span<const double> bcord,
                                std::span<double> coor,
                                std::span<int> kbndpt,
                                std::span<int> kbound,
                                bool fillKbound,
                                const SepranContext& /*ctx*/)
{
    // Initialize outputs
    int nbn    = 0;
    int jpnt   = 0;
    int inside = -1;

    if (nbound == 0)
        return {jpnt, inside, nbn};

    // bcord is stored interleaved: bcord[2*i]=x, bcord[2*i+1]=y (0-based node i)
    // Compute reference distance from first two points
    const double dx  = bcord[2] - bcord[0];
    const double dy  = bcord[3] - bcord[1];
    const double ref = 1.0e-5 * std::sqrt(dx * dx + dy * dy);

    // -------------------------------------------------------------------
    // Detect closed-loop boundaries and store the index of the last node
    // of each loop in ihelp.
    // -------------------------------------------------------------------
    std::vector<int> ihelp; // ihelp[j] = 0-based index of last node of loop j
    ihelp.reserve(1 + nbound / 2);

    int ifirst = 0;
    while (true)
    {
        const double xst = bcord[2 * ifirst];
        const double yst = bcord[2 * ifirst + 1];

        // Find the last node that coincides with the starting node
        int ilast = nbound - 1;
        for (int i = ifirst + 1; i < nbound; ++i)
        {
            const double dxi = bcord[2 * i] - xst;
            const double dyi = bcord[2 * i + 1] - yst;
            if (std::sqrt(dxi * dxi + dyi * dyi) <= ref)
                ilast = i;
        }

        ihelp.push_back(ilast);
        if (ilast == nbound - 1)
            break;
        ifirst = ilast + 1;
    }

    const int nparts = static_cast<int>(ihelp.size());

    // -------------------------------------------------------------------
    // Loop over all closed parts
    // -------------------------------------------------------------------
    int partStart = 0; // 0-based index of first node of current part
    for (int j = 0; j < nparts; ++j)
    {
        const int partEnd = ihelp[j]; // 0-based index of last node of current part

        // Start of this part
        const int jst = jpnt + 1; // 1-based first node number of this part

        for (int i = partStart; i <= partEnd; ++i)
        {
            const double xp = bcord[2 * i];
            const double yp = bcord[2 * i + 1];

            if (i == partStart)
            {
                // First node of a new closed loop
                if (inside < 1) ++inside;
                ++jpnt;
                coor[2 * (jpnt - 1)]     = xp; // 1-based jpnt → 0-based index
                coor[2 * (jpnt - 1) + 1] = yp;
                if (fillKbound)
                    kbndpt[i] = jpnt;
            }
            else
            {
                // Subsequent node in the loop
                if (fillKbound)
                {
                    ++nbn;
                    kbound[2 * (nbn - 1)]     = jpnt;
                    kbound[2 * (nbn - 1) + 1] = jpnt + 1;
                }

                if (i == partEnd)
                {
                    // Last node of closed loop — wraps back to first
                    if (fillKbound)
                    {
                        kbound[2 * (nbn - 1) + 1] = jst;
                        kbndpt[i]                  = jst;
                    }
                }
                else
                {
                    ++jpnt;
                    coor[2 * (jpnt - 1)]     = xp;
                    coor[2 * (jpnt - 1) + 1] = yp;
                    if (fillKbound)
                        kbndpt[i] = jpnt;
                }
            }
        }

        partStart = partEnd + 1;
    }

    return {jpnt, inside, nbn};
}

// ---------------------------------------------------------------------------
// checkStaple  (mshchkstapl.for)
// ---------------------------------------------------------------------------

void checkStaple(std::span<const int> kstapl, int lenkstapl, std::string_view text)
{
    if (lenkstapl == 0)
        return;

    // kstapl is stored flat: [n1_0, n2_0, n1_1, n2_1, ...]
    // Sort column 1 (first nodes) and column 2 (second nodes) separately,
    // then check that they contain the same set of node numbers.

    std::vector<int> iwork(lenkstapl);
    std::vector<int> ihelp(lenkstapl);
    std::vector<int> iseq1(lenkstapl);
    std::vector<int> iseq2(lenkstapl);

    // Column 1
    for (int i = 0; i < lenkstapl; ++i)
        iwork[i] = kstapl[2 * i];
    heapSort(iwork, ihelp);
    for (int i = 0; i < lenkstapl; ++i)
        iseq1[i] = iwork[ihelp[i]];

    // Column 2
    for (int i = 0; i < lenkstapl; ++i)
        iwork[i] = kstapl[2 * i + 1];
    heapSort(iwork, ihelp);
    for (int i = 0; i < lenkstapl; ++i)
        iseq2[i] = iwork[ihelp[i]];

    // Compare
    for (int i = 0; i < lenkstapl; ++i)
    {
        if (iseq1[i] != iseq2[i])
        {
            std::cerr << text << '\n'
                      << "Both columns of kstapl do not contain the same node numbers\n"
                      << "The " << i << "-th nodes are different\n";
            for (int k = 0; k < lenkstapl; ++k)
                std::cerr << k << ":  " << iseq1[k] << "  " << iseq2[k] << '\n';
            break; // Diagnostic only; do not throw (matches Fortran behaviour)
        }
    }
}

} // namespace sepran
