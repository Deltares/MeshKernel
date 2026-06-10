// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Stichting Deltares
//
// C++20 translation of mshcurvinters.for, mshcurvinters1.for, and mshcurvinters2.for
// from the SEPRAN library.
//
// Original author: Guus Segal, 2005.

#include "SepranCurveIntersection.hpp"
#include "SepranGeometry.hpp"

#include <array>
#include <stdexcept>
#include <string>

namespace sepran
{

// Helper: access the flat xbox 3D array stored in Fortran column-major order:
//   xbox(coord, curve, minmax)  →  xbox[coord + 2*curve + 2*ncurvs*minmax]
// coord  = 0 (x) or 1 (y)
// curve  = 0-based curve index
// minmax = 0 (minimum) or 1 (maximum)
static double xboxAt(std::span<const double> xbox, int ncurvs, int coord, int curve, int minmax)
{
    return xbox[coord + 2 * curve + 2 * ncurvs * minmax];
}

// ---------------------------------------------------------------------------
// curvIntersectionPair  (mshcurvinters1.for)
// ---------------------------------------------------------------------------

void curvIntersectionPair(int icurnr,
                          int jcurnr,
                          std::span<const int> icurvs,
                          std::span<const double> curves,
                          SepranContext& ctx)
{
    if (ctx.ierror != 0)
        return;

    const double eps = ctx.sqreps;

    // Start positions (0-based node index into flat curves array)
    const int istart = (icurnr == 0) ? 0 : icurvs[icurnr - 1];
    const int inodes = icurvs[icurnr] - istart;

    const int jstart = (jcurnr == 0) ? 0 : icurvs[jcurnr - 1];
    const int jnodes = icurvs[jcurnr] - jstart;

    // curves is stored interleaved: x at [2*k], y at [2*k+1] (0-based node index k)
    std::array<double, 2> x1{curves[2 * istart], curves[2 * istart + 1]};

    for (int j = 1; j < inodes; ++j)
    {
        const std::array<double, 2> x2{curves[2 * (istart + j)], curves[2 * (istart + j) + 1]};

        std::array<double, 2> x3{curves[2 * jstart], curves[2 * jstart + 1]};

        for (int k = 1; k < jnodes; ++k)
        {
            const std::array<double, 2> x4{curves[2 * (jstart + k)], curves[2 * (jstart + k) + 1]};

            double fact1, fact2;
            crossLine1(x1, x2, x3, x4, fact1, fact2, eps, ctx);

            // Interior intersection only (not at shared endpoints)
            if (fact1 > eps && fact1 < 1.0 - eps &&
                fact2 > eps && fact2 < 1.0 - eps)
            {
                ctx.ierror = 1;
                throw std::runtime_error(
                    "SEPRAN: curves " + std::to_string(icurnr) + " and " +
                    std::to_string(jcurnr) + " intersect (error 2787)");
            }

            x3 = x4;
        }
        x1 = x2;
    }
}

// ---------------------------------------------------------------------------
// curvIntersectionCheck  (mshcurvinters.for)
// ---------------------------------------------------------------------------

void curvIntersectionCheck(int ncurvs,
                           std::span<const int> iwork,
                           std::span<const double> xbox,
                           std::span<const int> icurvs,
                           std::span<const double> curves,
                           SepranContext& ctx)
{
    if (ctx.ierror != 0)
        return;

    for (int icurnr = 0; icurnr < ncurvs; ++icurnr)
    {
        if (iwork[icurnr] != 0)
            continue; // Only process single curves

        for (int jcurnr = icurnr + 1; jcurnr < ncurvs; ++jcurnr)
        {
            if (iwork[jcurnr] != 0)
                continue;

            // Quick bounding-box check (Fortran: go to 200 on no overlap)
            if (xboxAt(xbox, ncurvs, 0, icurnr, 0) > xboxAt(xbox, ncurvs, 0, jcurnr, 1))
                continue;
            if (xboxAt(xbox, ncurvs, 0, jcurnr, 0) > xboxAt(xbox, ncurvs, 0, icurnr, 1))
                continue;
            if (xboxAt(xbox, ncurvs, 1, icurnr, 0) > xboxAt(xbox, ncurvs, 1, jcurnr, 1))
                continue;
            if (xboxAt(xbox, ncurvs, 1, jcurnr, 0) > xboxAt(xbox, ncurvs, 1, icurnr, 1))
                continue;

            curvIntersectionPair(icurnr, jcurnr, icurvs, curves, ctx);
            if (ctx.ierror != 0)
                return;
        }
    }
}

// ---------------------------------------------------------------------------
// boundarySelfIntersectionCheck  (mshcurvinters2.for)
// ---------------------------------------------------------------------------

void boundarySelfIntersectionCheck(std::span<const int> kbound,
                                   int nbound,
                                   std::span<const double> coor,
                                   int isurnr,
                                   SepranContext& ctx)
{
    if (ctx.ierror != 0)
        return;

    constexpr double eps = 1.0e-3;

    for (int i = 0; i < nbound - 1; ++i)
    {
        // kbound stores 1-based node indices; convert to 0-based for coor access.
        const int n1i = kbound[2 * i]     - 1; // first node of edge i (0-based)
        const int n2i = kbound[2 * i + 1] - 1; // second node of edge i (0-based)

        const std::array<double, 2> x1{coor[2 * n1i], coor[2 * n1i + 1]};
        const std::array<double, 2> x2{coor[2 * n2i], coor[2 * n2i + 1]};

        for (int j = i + 1; j < nbound; ++j)
        {
            const int n1j = kbound[2 * j]     - 1;
            const int n2j = kbound[2 * j + 1] - 1;

            const std::array<double, 2> x3{coor[2 * n1j], coor[2 * n1j + 1]};
            const std::array<double, 2> x4{coor[2 * n2j], coor[2 * n2j + 1]};

            double fact1, fact2;
            crossLine1(x1, x2, x3, x4, fact1, fact2, eps, ctx);

            if (fact1 > eps && fact1 < 1.0 - eps &&
                fact2 > eps && fact2 < 1.0 - eps)
            {
                ctx.ierror = 1;
                throw std::runtime_error(
                    "SEPRAN: boundary of surface " + std::to_string(isurnr) +
                    " intersects itself at edges " + std::to_string(i) +
                    " and " + std::to_string(j) + " (error 2788)");
            }
        }
    }
}

} // namespace sepran
