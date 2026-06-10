// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Stichting Deltares
//
// C++20 translation of msho02, msho03, msho10, msho18, msho19, msho22
// from the SEPRAN library.
//
// Original author: Niek Praagman, 1989-2010.

#include "SepranTopology.hpp"

#include <algorithm>
#include <cmath>

namespace sepran
{

// ---------------------------------------------------------------------------
// insertNeighbour  (msho02.for)
// ---------------------------------------------------------------------------

void insertNeighbour(std::span<const int> istart,
                     std::span<int> ibuur,
                     int i,
                     int j)
{
    // Range of neighbour slots for node i (1-based i)
    const int ist = (i == 1) ? 0 : istart[i - 2];   // 0-based start
    const int ien = istart[i - 1];                   // exclusive end

    for (int kn = ist; kn < ien; ++kn)
    {
        const int ib = ibuur[kn];
        if (ib == 0)
        {
            ibuur[kn] = j;
            return;
        }
        if (ib == j)
            return; // Already present — double line, do nothing
    }
}

// ---------------------------------------------------------------------------
// nodeDistance  (msho03.for)
// ---------------------------------------------------------------------------

double nodeDistance(int i, int j, std::span<const double> coor)
{
    const double dx = coor[2 * (i - 1)]     - coor[2 * (j - 1)];
    const double dy = coor[2 * (i - 1) + 1] - coor[2 * (j - 1) + 1];
    return std::sqrt(dx * dx + dy * dy);
}

// ---------------------------------------------------------------------------
// addElement  (msho10.for)
// ---------------------------------------------------------------------------

void addElement(std::span<int> kelem, int& nelem,
                int i1, int i2, int jpn,
                std::span<int> kstapl, int& kstap,
                std::span<int> itri)
{
    // Append the new element
    kelem[3 * nelem]     = i1;
    kelem[3 * nelem + 1] = i2;
    kelem[3 * nelem + 2] = jpn;
    ++nelem;

    // Walk through kstapl (skip slot 0 — Fortran starts at index 2)
    // Remove edges i2-jpn and jpn-i1 if present; otherwise they are new.
    int idrie = 0; // Number of surviving edges
    int ieen  = 1; // 1 = edge jpn-i2 is new (needs adding)
    int itwe  = 1; // 1 = edge i1-jpn is new (needs adding)

    for (int ks = 1; ks < kstap; ++ks) // Fortran: do ks = 2, kstap
    {
        const int ia = kstapl[2 * ks];
        const int ib = kstapl[2 * ks + 1];
        if (ia == i2 && ib == jpn)
        {
            ieen = 0;
        }
        else if (ia == jpn && ib == i1)
        {
            itwe = 0;
        }
        else
        {
            kstapl[2 * idrie]     = ia;
            kstapl[2 * idrie + 1] = ib;
            ++idrie;
        }
    }

    if (ieen == 1)
    {
        kstapl[2 * idrie]     = jpn;
        kstapl[2 * idrie + 1] = i2;
        ++idrie;
        ++itri[i2  - 1];
        ++itri[jpn - 1];
    }
    else
    {
        --itri[i2  - 1];
        --itri[jpn - 1];
    }

    if (itwe == 1)
    {
        kstapl[2 * idrie]     = i1;
        kstapl[2 * idrie + 1] = jpn;
        ++idrie;
        ++itri[jpn - 1];
        ++itri[i1  - 1];
    }
    else
    {
        --itri[jpn - 1];
        --itri[i1  - 1];
    }

    kstap = idrie;
}

// ---------------------------------------------------------------------------
// triangleArea  (msho19.for)
// ---------------------------------------------------------------------------

double triangleArea(std::span<const double> coor, int i1, int i2, int i3)
{
    const double x1 = coor[2 * (i1 - 1)];
    const double y1 = coor[2 * (i1 - 1) + 1];
    const double x2 = coor[2 * (i2 - 1)];
    const double y2 = coor[2 * (i2 - 1) + 1];
    const double x3 = coor[2 * (i3 - 1)];
    const double y3 = coor[2 * (i3 - 1) + 1];
    return (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2.0;
}

// ---------------------------------------------------------------------------
// checkTriangleAngles  (msho22.for)
// ---------------------------------------------------------------------------

int checkTriangleAngles(double a, double b, double c)
{
    const double aa = a * a;
    const double bb = b * b;
    const double cc = c * c;

    constexpr double q = -0.01;

    const double angle1 = (cc + aa - bb) / (2.0 * c * a);
    const double angle2 = (aa + bb - cc) / (2.0 * a * b);
    const double angle3 = (bb + cc - aa) / (2.0 * b * c);

    if      (angle1 < q) return 1;
    else if (angle2 < q) return 2;
    else if (angle3 < q) return 3;
    return 0;
}

// ---------------------------------------------------------------------------
// laplacianSmoothing  (msho18.for)
// ---------------------------------------------------------------------------

void laplacianSmoothing(std::span<const int> kelem, int nelem,
                        std::span<double> coor, int npoint, int nipnt,
                        std::span<const int> iwork,
                        std::span<const int> ibuurp,
                        double& coars,
                        SepranContext& ctx)
{
    coars = 0.02 * coars;

    int jr = 1;
    bool finished = false;

    while (!finished)
    {
        ++jr;
        double afstnd = 0.0;

        for (int ikn = nipnt + 1; ikn <= npoint; ++ikn) // 1-based node index
        {
            // Compute total surface of all elements containing ikn
            double surf  = 0.0;
            double surf1 = 0.0;

            const double xh = coor[2 * (ikn - 1)];
            const double yh = coor[2 * (ikn - 1) + 1];

            for (int i = 0; i < nelem; ++i)
            {
                const int ia = kelem[3 * i];
                const int ib = kelem[3 * i + 1];
                const int ic = kelem[3 * i + 2];
                if (ia == ikn || ib == ikn || ic == ikn)
                {
                    surf1 = triangleArea(coor, ia, ib, ic);
                    surf += surf1;
                }
            }

            const double surfre = surf;

            // CSR start/end for node ikn
            const int jstart = (ikn == 1) ? 0 : iwork[ikn - 2];
            const int jeind  = iwork[ikn - 1];

            if (jstart >= jeind)
                continue;

            // Compute barycentre of neighbours
            double xz    = 0.0;
            double yz    = 0.0;
            int    nbuur = jeind - jstart;

            for (int i = jstart; i < jeind; ++i)
            {
                const int jkn = ibuurp[i];
                if (jkn > 0)
                {
                    xz += coor[2 * (jkn - 1)];
                    yz += coor[2 * (jkn - 1) + 1];
                }
                else
                {
                    --nbuur;
                }
            }
            xz /= nbuur;
            yz /= nbuur;

            // Try over-relaxed move
            const double xn = xh + 1.62 * (xz - xh);
            const double yn = yh + 1.62 * (yz - yh);

            coor[2 * (ikn - 1)]     = xn;
            coor[2 * (ikn - 1) + 1] = yn;

            // Recompute surface and check validity
            surf = 0.0;
            for (int i = 0; i < nelem; ++i)
            {
                const int ia = kelem[3 * i];
                const int ib = kelem[3 * i + 1];
                const int ic = kelem[3 * i + 2];
                if (ia == ikn || ib == ikn || ic == ikn)
                {
                    surf1 = triangleArea(coor, ia, ib, ic);
                    surf += std::abs(surf1);
                }
            }

            int iallow = 0;
            if (std::abs(surf - surfre) < 100.0 * ctx.epsmac)
                iallow = 2;

            if (iallow != 2)
            {
                // Over-relaxation not accepted — try barycentre
                coor[2 * (ikn - 1)]     = xz;
                coor[2 * (ikn - 1) + 1] = yz;

                surf = 0.0;
                for (int i = 0; i < nelem; ++i)
                {
                    const int ia = kelem[3 * i];
                    const int ib = kelem[3 * i + 1];
                    const int ic = kelem[3 * i + 2];
                    if (ia == ikn || ib == ikn || ic == ikn)
                    {
                        surf1 = triangleArea(coor, ia, ib, ic);
                        surf += surf1;
                    }
                }
                if (std::abs(surf - surfre) < 100.0 * ctx.epsmac)
                    iallow = 1;

                if (iallow == 0)
                {
                    // Neither move is valid — restore old position
                    coor[2 * (ikn - 1)]     = xh;
                    coor[2 * (ikn - 1) + 1] = yh;
                }
            }

            // Track maximum displacement
            const double cx  = coor[2 * (ikn - 1)];
            const double cy  = coor[2 * (ikn - 1) + 1];
            const double af  = (xh - cx) * (xh - cx) + (yh - cy) * (yh - cy);
            if (af > afstnd) afstnd = af;
        }

        if (jr > 2 && std::sqrt(afstnd) < coars) finished = true;
        if (jr == 20) finished = true;
    }
}

} // namespace sepran
