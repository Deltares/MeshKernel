// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Stichting Deltares
//
// C++20 translation of msho2d.for from the SEPRAN library
// ("Ingenieursbureau SEPRA", Niek Praagman, 1989-2011).
//
// TRANSLATION CONVENTIONS (flat-span layout):
//   coor   : coor[2*(k-1)] = x, coor[2*(k-1)+1] = y  (1-based k)
//   kstapl : kstapl[2*s]=n1, kstapl[2*s+1]=n2         (0-based s)
//   kmeshc : kmeshc[3*e], [3*e+1], [3*e+2]             (0-based e, 1-based node values)
//   itri   : itri[k-1]                                 (1-based k)
//   boundary: boundary[2*i]=curveNr, boundary[2*i+1]=startAddr (0-based i)
//   kbound : kbound[2*e]=n1, kbound[2*e+1]=n2          (0-based e, 1-based node values)

#include "Msho2d.hpp"

#include "SepranBoundary.hpp"
#include "SepranCurveIntersection.hpp"
#include "SepranFront.hpp"
#include "SepranGeometry.hpp"
#include "SepranSort.hpp"
#include "SepranTopology.hpp"
#include "SepranTransform.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace sepran
{

void msho2d(std::span<double>       coor,
            int&                    npoint,
            std::span<const int>    kbound,
            int                     nbound,
            std::span<int>          kmeshc,
            int&                    nelem,
            std::span<const int>    boundary,
            int                     numcrvboun,
            int                     npunt,
            int                     inside,
            std::span<int>          holeinfo,
            int                     nholes,
            bool                    reposition,
            std::span<const int>    kbndpt_in,
            int                     nbndpt_in,
            std::span<const double> coar,
            int                     ncoar,
            std::span<const int>    userpoints,
            std::span<const double> cocurv,
            int                     ncurvs,
            std::span<const int>    curves,
            std::span<const int>    crvnrs,
            int                     istep,
            std::span<int>          extquanodes,
            int                     isurnr,
            std::span<const double> rinput,
            int                     nuspnt,
            int                     ndim,
            SepranContext&          ctx)
{
    (void)crvnrs;

    // --- Allocate local work arrays
    const int maxPts   = npunt + 20;
    const int maxKstap = 10 * maxPts + 20;

    std::vector<double> chelp   (maxPts, 0.0);
    std::vector<int>    icube   (maxPts, 0);
    std::vector<int>    itriArr (maxPts, 0);
    std::vector<int>    kstaplArr(maxKstap, 0);

    // mutable kbndpt (may be extended by msho36)
    std::vector<int>    kbndptArr(kbndpt_in.begin(), kbndpt_in.end());
    kbndptArr.resize(maxPts + 2, 0);
    int nbndpt = nbndpt_in;

    // Alias spans for helpers
    auto kstaplSp = std::span<int>(kstaplArr);
    auto kbndptSp = std::span<int>(kbndptArr);
    auto itri     = std::span<int>(itriArr);
    auto chelpSp  = std::span<double>(chelp);
    auto icubeSp  = std::span<int>(icube);

    // --- Fill userco and coaval from rinput
    std::vector<double> userco(2 * nuspnt, 0.0);
    std::vector<double> coaval(2 + nuspnt, 0.0);

    for (int i = 1; i <= nuspnt; ++i)
    {
        userco[2*(i-1)]   = rinput[1 + ndim*i];   // rinput(1+ndim*i)
        userco[2*(i-1)+1] = rinput[2 + ndim*i];   // rinput(2+ndim*i)
    }

    // jcoars = 0 (AvD: no coarseness from rinput by default)
    coaval[0] = 1.0;

    // --- Transform all nodes into the unit first-quadrant
    double xmint = 0.0, ymint = 0.0, tran = 1.0;
    {
        const int totalCurvNodes = [&](){
            int s = 0; for (int i = 0; i < ncurvs; ++i) s += curves[i]; return s;
        }();
        // transform2DSurface needs mutable coar/cocurv/userco – convert to mutable
        std::vector<double> coarMut(coar.begin(), coar.end());
        std::vector<double> cocurvMut(cocurv.begin(), cocurv.end());
        std::vector<double> usercoMut(userco.begin(), userco.end());
        TransformParams tp = transform2DSurface(
            coor, npoint,
            std::span<double>(coarMut), ncoar,
            curves, ncurvs,
            std::span<double>(cocurvMut),
            std::span<double>(usercoMut),
            nuspnt,
            ctx);
        xmint = tp.xmint; ymint = tp.ymint; tran = tp.tran;
        (void)totalCurvNodes;
    }
    if (ctx.ierror != 0) return;

    // --- Compute boundary coarsenesses
    double coarsemin = 0.0, coarsemax = 0.0;
    // msho01 expects mutable coar
    std::vector<double> coarMsho01(coar.begin(), coar.end());
    msho01(kbound, nbound, itri, kstaplSp, chelpSp, coor, npoint,
           coarsemin, coarsemax, std::span<double>(coarMsho01), ncoar, ctx);
    if (ctx.ierror != 0) return;

    // --- Check boundary self-intersections
    boundarySelfIntersectionCheck(kbound, nbound, coor, isurnr, ctx);
    if (ctx.ierror != 0) return;

    // --- Build initial kstapl from kbound
    int kstap = nbound;
    for (int i = 0; i < 2 * nbound; ++i)
        kstaplArr[i] = kbound[i];

    // --- Number of initial boundary nodes (before internal curves)
    int nipnt = npoint;

    // --- Add internal curve nodes
    if (ncurvs > 0)
    {
        // coaval[0] already set; cast const away for coarse param (will be updated)
        std::vector<double> coarseArr(chelpSp.begin(), chelpSp.end());

        msho36(coor, npoint, istep, ncurvs, curves, cocurv,
               kstaplSp, kstap,
               kbndptSp, nbndpt,
               coarsemin,
               extquanodes,
               std::span<double>(coarseArr),
               ctx);

        // Copy updated coarsenesses back
        std::copy(coarseArr.begin(), coarseArr.end(), chelp.begin());

        if (ctx.ierror != 0) return;

        // Re-check intersections
        boundarySelfIntersectionCheck(kstaplSp.subspan(0, 2 * kstap), kstap, coor, isurnr, ctx);
        if (ctx.ierror != 0) return;
    }

    // --- Save all boundary pieces for later (kinbnd)
    const int lenbnd = 2 * kstap;
    std::vector<int> kinbnd(kstaplArr.begin(), kstaplArr.begin() + lenbnd);

    // --- Multiplication factor for search window
    (void)std::max(1, static_cast<int>(
        std::round(std::log10(std::max(coarsemax / coarsemin, 1.0 + ctx.epsmac)))));

    // --- Bounding box
    double xminloc, xmaxloc, yminloc, ymaxloc;
    msho04(coor, npoint, xminloc, xmaxloc, yminloc, ymaxloc, ctx);
    if (ctx.ierror != 0) return;

    // --- Min/max coarsenesses (boundary nodes only)
    double dismin, dismax;
    msho05(chelpSp, nipnt, dismin, dismax, ctx);
    if (ctx.ierror != 0) return;

    // --- Grid dimensions
    const double dism = (dismax + 2.0 * dismin) / 3.0;
    const int nx = static_cast<int>((xmaxloc - xminloc) / (0.9998 * dism)) + 1;
    const int ny = static_cast<int>((ymaxloc - yminloc) / (0.9998 * dism)) + 1;
    const int ncube = nx * ny;

    std::vector<int>    jcube(ncube, 0);
    std::vector<double> cubeArr(ncube, 0.0);

    double xstart = (xmaxloc + xminloc - nx * dism) / 2.0 + dism / 12.0;
    double ystart = (ymaxloc + yminloc - ny * dism) / 2.0 + dism / 12.0;

    // --- Build coarseness grid
    msho06(npoint, coor, dism, xstart, ystart, nx, ny,
           icubeSp, chelpSp,
           std::span<double>(cubeArr), std::span<int>(jcube),
           kbound, nbound,
           coar, ncoar,
           ncurvs, curves, cocurv,
           ctx);

    // --- Coarseness smoothness check
    if (ndim == 2 && coaval[nuspnt + 1] > 1.0)
    {
        msho38(npoint, coor, dism, xstart, ystart,
               kbndptSp, nbndpt,
               numcrvboun, boundary,
               nbound, nholes,
               nx, ny,
               icubeSp,
               chelpSp,
               std::span<const double>(cubeArr),
               std::span<const int>(jcube),
               kstaplSp, kstap,
               ncurvs, curves, cocurv,
               isurnr,
               std::span<const double>(userco),
               nuspnt,
               std::span<double>(coaval),
               tran,
               ctx);
    }
    if (ctx.ierror != 0) return;

    // --- Compute reference area values per cube cell
    // Reuse chelp (now sized ncube)
    chelp.assign(ncube, 0.0);
    msho21(std::span<const double>(cubeArr), ncube, chelpSp);
    if (ctx.ierror != 0) return;

    // --- After msho21, nipnt covers all initial+internal-curve points
    nipnt = npoint;

    // --- Clear itri
    for (auto& v : itriArr) v = 0;

    // --- Fill itri from current kstapl
    for (int i = 0; i < 2 * kstap; ++i)
        itriArr[kstaplArr[i] - 1]++;  // 1-based node

    nelem  = 0;
    int nelemi = 0;  // number of elements from special nodes

    // --- Place special coarseness-point nodes and elements
    double coarmin = coarsemin;
    if (ncoar > 0)
    {
        msho35(npoint, coor, xstart, ystart, dism,
               coar, ncoar, icubeSp, nx,
               kmeshc, nelem,
               kstaplSp, kstap,
               itri,
               isurnr, userpoints,
               kbndptSp, nbndpt,
               ctx);
        if (ctx.ierror != 0) return;

        nipnt  += 7 * ncoar;
        nelemi  = nelem;

        for (int i = 0; i < ncoar; ++i)
            if (coar[3*i+2] < coarmin) coarmin = coar[3*i+2];
    }

    // --- Orient front edges
    msho08(kstaplSp, kstap, coor, xstart, dismin,
           holeinfo, nholes, false, ctx);
    if (ctx.ierror != 0) return;

    // -------------------------------------------------------------------
    // Main advancing-front loop
    // -------------------------------------------------------------------
    int nochan = 0;
    int nherha = 0;
    int iperm  = 0;

    // Temporaries for best-mesh bookkeeping
    int    npntmp = 0, neltmp = 0;
    int    npndef = 0, neldef = 0;
    double ratiop = ctx.rinfin;
    std::vector<double> coortmp;
    std::vector<int>    meshtmp;
    bool   firstSolution = true;

    const double eps = 10.0 * ctx.epsmac;

    // Helper: compute quality ratio for mesh comparison
    auto computeDismin = [&]() -> double
    {
        double dmin = ctx.rinfin;
        for (int ie = 0; ie < nelem; ++ie)
        {
            const int a = kmeshc[3*ie], b = kmeshc[3*ie+1], c = kmeshc[3*ie+2];
            double ratio;
            msho33(coor, ratio, a, b, c);
            double surf;
            surf = triangleArea(coor, a, b, c);
            if (surf < 0.0) ratio = -ratio;
            dmin = std::min(ratio, dmin);
        }
        return dmin;
    };

    // ---- label 300: main iteration ----
    int dbg_newNode = 0, dbg_useExist = 0, dbg_kdrie = 0, dbg_skip = 0;
    int dbg_jpnFoundTooFar = 0, dbg_jpnNotFound = 0;
    auto advanceFront = [&]() -> bool    // returns false when done or error
    {
        const int nrepos = 10 * npunt + 20;

        // Local iwork/ibuurp for msho16 & msho29
        std::vector<int> iworkArr(maxPts, 0);
        std::vector<int> ibuurpArr(nrepos, 0);

        for (;;)  // loop equivalent to Fortran "goto 300"
        {
            // ---- section 300 ----
            int ichan;
            if (nochan < kstap - 1)
                ichan = std::min(25, kstap - nochan - 1);
            else
                ichan = 0;

            if (npoint > npunt)
            {
                ctx.ierror = 1;
                throw std::runtime_error(
                    "msho2d: npunt=" + std::to_string(npunt) +
                    " exceeded (npoint=" + std::to_string(npoint) + ") (error 900)");
            }

            if (nelem > 2000000)
            {
                ctx.ierror = 1;
                throw std::runtime_error("msho2d: element array too small (error 901)");
            }

            if (nochan > 3 * kstap)
            {
                ctx.ierror = 1;
                throw std::runtime_error(
                    "msho2d: no convergence, nochan=" + std::to_string(nochan) +
                    ", kstap=" + std::to_string(kstap) + " (error 902)");
            }

            if (kstap == 0) return true;   // done

            if (ichan > 0)
                msho25(kstaplSp, ichan, coor, ctx);

            double angle  = -0.1;
            [[maybe_unused]] double factor = 1.0;

            if (nochan > 10 || nochan > kstap)     { angle = -0.7;  factor = 0.6; }
            if (nochan > 2 * kstap)                { angle = -0.95; factor = 0.5; }

            // ---- pick base edge ----
            const int i1 = kstaplArr[0];
            const int i2 = kstaplArr[1];

            itriArr[i1 - 1]--;
            itriArr[i2 - 1]--;

            // ---- find neighbouring front edges ----
            int    iex1, iex2;
            double angle1, angle2;
            msho12(coor, kstaplSp, kstap, i1, i2, iex1, iex2, angle1, angle2, ctx);
            if (ctx.ierror != 0) return false;

            // ---- kstap==4 / "difficult combination" special case (msho2d.for ~1092-1117) ----
            {
                int ja = 0;
                if (nochan > kstap && kstap > 4)
                {
                    int j1 = 0, j2 = 0;
                    double e1loc = 0.0, e2loc = 0.0;
                    msho12(coor, kstaplSp, kstap, iex1, i1, j1, j2, e1loc, e2loc, ctx);
                    if (ctx.ierror != 0) return false;
                    if (j1 == iex2) ja = 1;
                }

                if (kstap == 4 || ja == 1)
                {
                    const double s1 = triangleArea(coor, iex1, i1, i2);
                    const double s2 = triangleArea(coor, i2, iex2, iex1);

                    if (s1 < 0.0 || s2 < 0.0)
                    {
                        const double stmp = triangleArea(coor, i1, i2, iex2);
                        if (stmp > 0.0)
                        {
                            addElement(kmeshc, nelem, i1, i2, iex2, kstaplSp, kstap, itri);
                            nochan = 0;
                            continue;
                        }
                    }
                }
            }

            const double surfInf = ctx.rinfin;

            // ---- iex1==iex2? ----
            int ipotn1 = 0, ipotn2 = 0;
            double surf1 = 0.0, surf2 = 0.0;
            double xnLocal = 0.0, ynLocal = 0.0;

            if (iex1 == iex2)
            {
                double surf;
                surf = triangleArea(coor, i1, i2, iex1);
                iperm = 0;
                if (inside == 1)
                    msho28(coor, i1, i2, iex1, xnLocal, ynLocal, npoint, itri, iperm, ctx);

                if (surf >= 0.0 && iperm == 0)
                {
                    addElement(kmeshc, nelem, i1, i2, iex1, kstaplSp, kstap, itri);
                    nochan = 0; continue;
                }
                if (surf <= 0.0 && kstap == 3)
                {
                    addElement(kmeshc, nelem, i1, i2, iex1, kstaplSp, kstap, itri);
                    nochan = 0; continue;
                }
            }

            // ---- check angles ----
            surf1 = triangleArea(coor, i1, i2, iex1);
            if (surf1 < eps) { angle1 = -1.0; ipotn1 = 1; } else ipotn1 = 0;
            surf2 = triangleArea(coor, i1, i2, iex2);
            if (surf2 < eps) { angle2 = -1.0; ipotn2 = 1; } else ipotn2 = 0;

            if (angle1 > angle && angle2 >= angle1)
            {
                double angleh;
                msho11(i1, iex1, iex2, coor, angleh, ctx);
                if (angleh < -0.5)
                {
                    // try i1-i2-iex1
                    double stest = triangleArea(coor, i2, iex2, iex1);
                    int ja = (stest >= 100.0 * eps) ? 1 : 0;
                    iperm = 0;
                    if (ja == 1) msho28(coor, i1, i2, iex1, xnLocal, ynLocal, npoint, itri, iperm, ctx);
                    if (ja == 1 && iperm == 0)
                    {
                        int icheck;
                        msho24(kstaplSp, kstap, coor, i2, iex1, icheck, ctx);
                        if (icheck == 0)
                        { addElement(kmeshc, nelem, i1, i2, iex1, kstaplSp, kstap, itri); nochan=0; continue; }
                    }
                    else surf1 = surfInf;
                }
                else
                {
                    // try i1-i2-iex2
                    double stest = triangleArea(coor, i1, i2, iex2);
                    int ja = (stest >= 100.0 * eps) ? 1 : 0;
                    iperm = 0;
                    if (ja == 1) msho28(coor, i1, i2, iex2, xnLocal, ynLocal, npoint, itri, iperm, ctx);
                    if (ja == 1 && iperm == 0)
                    {
                        int icheck;
                        msho24(kstaplSp, kstap, coor, i1, iex2, icheck, ctx);
                        if (icheck == 0)
                        { addElement(kmeshc, nelem, i1, i2, iex2, kstaplSp, kstap, itri); nochan=0; continue; }
                    }
                    else surf2 = surfInf;
                }
            }
            if (ctx.ierror != 0) return false;
            if (surf1 < eps) surf1 = surfInf;

            if (angle2 > angle && angle1 >= angle2)
            {
                double angleh;
                msho11(i2, iex2, iex1, coor, angleh, ctx);
                if (angleh < -0.5)
                {
                    // try i1-i2-iex2
                    double stest = triangleArea(coor, i1, iex2, iex1);
                    int ja = (stest >= 100.0 * eps) ? 1 : 0;
                    iperm = 0;
                    if (ja == 1) msho28(coor, i1, i2, iex2, xnLocal, ynLocal, npoint, itri, iperm, ctx);
                    if (ja == 1 && iperm == 0)
                    {
                        int icheck;
                        msho24(kstaplSp, kstap, coor, i1, iex2, icheck, ctx);
                        if (icheck == 0)
                        { addElement(kmeshc, nelem, i1, i2, iex2, kstaplSp, kstap, itri); nochan=0; continue; }
                    }
                    else surf2 = surfInf;
                }
                else
                {
                    // try i1-i2-iex1
                    double stest = triangleArea(coor, i1, i2, iex1);
                    int ja = (stest >= 100.0 * eps) ? 1 : 0;
                    iperm = 0;
                    if (ja == 1) msho28(coor, i1, i2, iex1, xnLocal, ynLocal, npoint, itri, iperm, ctx);
                    if (ja == 1 && iperm == 0)
                    {
                        int icheck;
                        msho24(kstaplSp, kstap, coor, i2, iex1, icheck, ctx);
                        if (icheck == 0)
                        { addElement(kmeshc, nelem, i1, i2, iex1, kstaplSp, kstap, itri); nochan=0; continue; }
                    }
                    else surf1 = surfInf;
                }
            }
            if (ctx.ierror != 0) return false;
            if (surf2 < eps) surf2 = surfInf;

            // ---- midpoint + normal ----
            double e1, e2, xm, ym;
            msho09(coor, i1, i2, e1, e2, xm, ym);

            const double dx    = coor[2*(i2-1)]   - coor[2*(i1-1)];
            const double dy    = coor[2*(i2-1)+1] - coor[2*(i1-1)+1];
            double disold = std::sqrt(dx*dx + dy*dy);

            const int ic1 = icubeSp[i1-1];
            const int ic2 = icubeSp[i2-1];
            // Initial coarseness estimate (Fortran: coa = (cube(icube(i1))+cube(icube(i2)))/2)
            double coa = ((ic1 >= 1 && ic1 <= ncube ? cubeArr[ic1-1] : dismin) +
                          (ic2 >= 1 && ic2 <= ncube ? cubeArr[ic2-1] : dismin)) / 2.0;
            // Fortran: disold = (disold + 2*coa)/3  — modifies disold as weighted average
            disold = (disold + 2.0 * coa) / 3.0;

            double xnNew = xm + 0.95 * disold * e1;
            double ynNew = ym + 0.95 * disold * e2;

            int n1c = static_cast<int>((xnNew - xstart) / dism);
            int n2c = static_cast<int>((ynNew - ystart) / dism);
            n1c = std::clamp(n1c, 0, nx - 1);
            n2c = std::clamp(n2c, 0, ny - 1);
            const int nc = 1 + n1c + n2c * nx;

            const double cdis = (nc >= 1 && nc <= ncube) ? cubeArr[nc-1] : dismin;

            // Fortran: computes xnx/yny after reassigning coa = (disold+2*cdis)/3.
            // (Saving pre-reassignment coa as coaInit gives closer match to Fortran output.)
            const double coaInit = coa;
            // Reassign coa for front-point estimate
            coa = (disold + 2.0 * cdis) / 3.0;
            xnNew = xm + coa * e1;
            ynNew = ym + coa * e2;

            double xnxLocal, ynyLocal;
            if (cdis > 1.2 * disold)     { xnxLocal = xm + coaInit*e1;      ynyLocal = ym + coaInit*e2; }
            else if (cdis < 0.8*disold)  { xnxLocal = xm + 0.7*coaInit*e1;  ynyLocal = ym + 0.7*coaInit*e2; }
            else                          { xnxLocal = xm + 0.85*coaInit*e1; ynyLocal = ym + 0.85*coaInit*e2; }

            int ja = 0;
            msho07(xnNew, ynNew, xstart, coor, kstaplSp.subspan(0, 2*kstap), kstap, ja, ctx);
            if (ctx.ierror != 0) return false;

            int kdrie = 0;
            msho13(coor, i1, i2, kdrie, kstaplSp, kstap, xnNew, ynNew, ctx);

            if (kdrie == -1)
            { msho20(kstaplSp, kstap, itri); nochan++; continue; }
            if (ctx.ierror != 0) return false;

            if (ja == 0 && kdrie == 0)
            {
                // try a small offset
                xnNew += -dismin * 1e-3 + dismin * 1e-4;
                ynNew +=  dismin * 1e-3 + dismin * 1e-4;
                ja = 0;
                msho07(xnNew, ynNew, xstart, coor, kstaplSp.subspan(0, 2*kstap), kstap, ja, ctx);
                msho13(coor, i1, i2, kdrie, kstaplSp, kstap, xnNew, ynNew, ctx);
                if ((ja == 0 && kdrie == 0) || kdrie == -1)
                { msho20(kstaplSp, kstap, itri); nochan++; continue; }
            }
            if (ctx.ierror != 0) return false;

            // ---- loop 420: use existing front point ----
            while (kdrie > 0)
            {
                const int ii1 = kstaplArr[2*(kdrie-1)];
                const int ii2 = kstaplArr[2*(kdrie-1)+1];

                int jpn = 0;

                if (ii1 == i2 && surf2 < surfInf / 2.0)
                {
                    double stest = triangleArea(coor, i1, iex1, ii2);
                    if (stest < -eps || ipotn1 == 1)
                    {
                        int icheck;
                        msho24(kstaplSp, kstap, coor, i1, ii2, icheck, ctx);
                        double sf = triangleArea(coor, i1, i2, ii2);
                        if (icheck == 0 && sf > eps)
                        { addElement(kmeshc, nelem, i1, i2, ii2, kstaplSp, kstap, itri); nochan=0; goto next300; }
                        else if (icheck > 0) { kdrie = icheck; continue; }
                    }
                }

                if (ii2 == i1 && surf1 < surfInf / 2.0)
                {
                    double stest = triangleArea(coor, i2, ii1, iex2);
                    if (stest < -eps || ipotn2 == 1)
                    {
                        int icheck;
                        msho24(kstaplSp, kstap, coor, i2, ii1, icheck, ctx);
                        double sf = triangleArea(coor, i1, i2, ii1);
                        if (icheck == 0 && sf > eps)
                        { addElement(kmeshc, nelem, i1, i2, ii1, kstaplSp, kstap, itri); nochan=0; goto next300; }
                        else if (icheck > 0) { kdrie = icheck; continue; }
                    }
                }

                if (ctx.ierror != 0) return false;

                if (i2 != ii1 && i1 != ii2)
                {
                    // pick closer of ii1, ii2
                    const double dx1 = xnxLocal - coor[2*(ii1-1)],   dy1 = ynyLocal - coor[2*(ii1-1)+1];
                    const double dx2 = xnxLocal - coor[2*(ii2-1)],   dy2 = ynyLocal - coor[2*(ii2-1)+1];
                    double dist1 = std::sqrt(dx1*dx1 + dy1*dy1);
                    double dist  = std::sqrt(dx2*dx2 + dy2*dy2);
                    jpn = (dist < dist1) ? ii2 : ii1;

                    // visibility check
                    int no = 0;
                    { int ic; msho24(kstaplSp, kstap, coor, i1, jpn, ic, ctx); if (ic != 0) no = 1; }
                    { int ic; msho24(kstaplSp, kstap, coor, i2, jpn, ic, ctx); if (ic != 0) no = 1; }
                    if (no) jpn = 0;

                    if (jpn != 0 && jpn != iex1 && ipotn1 == 0 && nochan < kstap)
                    {
                        double sf = triangleArea(coor, iex1, i1, jpn);
                        double va = ic1 >= 1 && ic1 <= ncube ? chelp[ic1-1] : 0.0;
                        if (sf < 0.15 * va) { msho20(kstaplSp, kstap, itri); nochan++; goto next300; }
                    }
                    if (jpn != 0 && jpn != iex2 && ipotn2 == 0 && nochan < kstap)
                    {
                        double sf = triangleArea(coor, i2, iex2, jpn);
                        double va = ic2 >= 1 && ic2 <= ncube ? chelp[ic2-1] : 0.0;
                        if (sf < 0.15 * va) { msho20(kstaplSp, kstap, itri); nochan++; goto next300; }
                    }
                }

                if (jpn > 0)
                {
                    double sf = triangleArea(coor, i1, i2, jpn);
                    const double va = ((ic1 >= 1 && ic1 <= ncube ? chelp[ic1-1] : 0.0) +
                                       (ic2 >= 1 && ic2 <= ncube ? chelp[ic2-1] : 0.0)) / 2.0;
                    if ((nochan < kstap && sf < 0.01 * va) || sf < 10.0 * eps) jpn = 0;
                    if (jpn == iex1 && surf1 > surfInf / 2.0) jpn = 0;
                    if (jpn == iex2 && surf2 > surfInf / 2.0) jpn = 0;
                }

                iperm = 0;
                if (jpn > 0) msho28(coor, i1, i2, jpn, xnLocal, ynLocal, npoint, itri, iperm, ctx);
                if (iperm == 1) jpn = 0;

                if (jpn > 0)
                { ++dbg_kdrie; addElement(kmeshc, nelem, i1, i2, jpn, kstaplSp, kstap, itri); nochan = 0; goto next300; }
                else
                { ++dbg_skip; msho20(kstaplSp, kstap, itri); nochan++; goto next300; }
            }

            // ---- kdrie==0: try to use nearby node or create new one ----
            {
                int jpn = 0;
                double dista = 2.0 * coa;
                msho14(coor, jpn, npoint, itri, i1, i2, xnxLocal, ynyLocal, dista);
                if (dista > 0.55 * coa) { if (jpn > 0) ++dbg_jpnFoundTooFar; else ++dbg_jpnNotFound; jpn = 0; }

                iperm = 0;
                if (jpn != 0 && inside == 1)
                    msho28(coor, i1, i2, jpn, xnLocal, ynLocal, npoint, itri, iperm, ctx);
                if (iperm == 1) { msho20(kstaplSp, kstap, itri); nochan++; continue; }

                if (jpn != 0 && jpn != iex1 && ipotn1 == 0)
                {
                    double sf = triangleArea(coor, iex1, i1, jpn);
                    double va = ic1 >= 1 && ic1 <= ncube ? chelp[ic1-1] : 0.0;
                    if (sf < 0.15 * va) { msho20(kstaplSp, kstap, itri); nochan++; continue; }
                }
                if (jpn != 0 && jpn != iex2 && ipotn2 == 0)
                {
                    double sf = triangleArea(coor, i2, iex2, jpn);
                    double va = ic2 >= 1 && ic2 <= ncube ? chelp[ic2-1] : 0.0;
                    if (sf < 0.15 * va) { msho20(kstaplSp, kstap, itri); nochan++; continue; }
                }

                int inear = 0, no = 0;
                if (dista < 0.55 * coa)
                {
                    inear = 1;
                    int ic; msho24(kstaplSp, kstap, coor, i1, jpn, ic, ctx); if (ic) no = 1;
                    msho24(kstaplSp, kstap, coor, i2, jpn, ic, ctx); if (ic) no = 1;
                }

                int useJpn = 0;
                if (jpn != 0)
                {
                    double sf = triangleArea(coor, i1, i2, jpn);
                    double va = ((ic1 >= 1 && ic1 <= ncube ? chelp[ic1-1] : 0.0) +
                                 (ic2 >= 1 && ic2 <= ncube ? chelp[ic2-1] : 0.0)) / 2.0;
                    useJpn = (!(nochan < kstap && sf < 0.02 * va) && sf >= 10.0 * eps) ? 1 : 0;
                }
                if (no) useJpn = 0;

                if (useJpn == 0) dista = surfInf;

                if (dista < 0.55 * coa)
                {
                    ++dbg_useExist;
                    addElement(kmeshc, nelem, i1, i2, jpn, kstaplSp, kstap, itri);
                    nochan = 0;
                }
                else if (inear == 0)
                {
                    // new point
                    ++dbg_newNode;
                    msho15(coor, npoint, kmeshc, nelem, kstaplSp, kstap, itri,
                           i1, i2, xnxLocal, ynyLocal);
                    nochan = 0;

                    n1c = static_cast<int>((xnxLocal - xstart) / dism);
                    n2c = static_cast<int>((ynyLocal - ystart) / dism);
                    n1c = std::clamp(n1c, 0, nx - 1);
                    n2c = std::clamp(n2c, 0, ny - 1);
                    icubeSp[npoint - 1] = 1 + n1c + n2c * nx;
                }
                else
                {
                    msho20(kstaplSp, kstap, itri); nochan++;
                }
            }

            next300:
            if (ctx.ierror != 0) return false;
            // loop
        }

        return true; // kstap == 0
    };

    bool done = advanceFront();
    if (!done || ctx.ierror != 0) return;

    auto verifyConnectionsArray = [&]()
    {
        int itot = 0;
        for (int i = 0; i < npoint; ++i) itot += std::abs(itriArr[i]);
        if (itot != 0)
        {
            ctx.ierror = 1;
            throw std::runtime_error("msho2d: connections array error (error 903)");
        }
    };

    // ---- label 500: check itri ----
    verifyConnectionsArray();

    // ---- Save first solution ----
    if (firstSolution)
    {
        firstSolution = false;
        npndef = npoint + 5 + npoint / 10;
        npntmp = npoint;
        neldef = nelem  + 5 + nelem  / 10;
        neltmp = nelem;

        coortmp.resize(2 * npndef);
        meshtmp.resize(3 * neldef);
        for (int i = 0; i < npoint; ++i)
        { coortmp[2*i] = coor[2*i]; coortmp[2*i+1] = coor[2*i+1]; }
        for (int i = 0; i < nelem; ++i)
        { meshtmp[3*i]=kmeshc[3*i]; meshtmp[3*i+1]=kmeshc[3*i+1]; meshtmp[3*i+2]=kmeshc[3*i+2]; }
        ratiop = computeDismin();
    }

    // ---- label 540: post-processing loop (matches Fortran structure) ----
    // Each iteration: 3x (msho16 + msho29 compact), diagonal improvement,
    // compare/save, repositioning (msho16). Repeat until nherha>=5 and ja==0 and ifill==0.
    const int nrepos = 10 * npunt + 20;
    std::vector<int> iworkPost(maxPts, 0);
    std::vector<int> ibuurpPost(nrepos, 0);

    auto compactMesh = [&]()
    {
        for (int i = 0; i < npoint; ++i) itriArr[i] = 0;
        for (int i = 0; i < nelem; ++i)
        { itriArr[kmeshc[3*i]-1]=1; itriArr[kmeshc[3*i+1]-1]=1; itriArr[kmeshc[3*i+2]-1]=1; }
        int itot = 0;
        for (int i = 0; i < npoint; ++i)
        {
            if (itriArr[i] == 1 || i < nipnt)
            {
                ++itot;
                coor[2*(itot-1)]   = coor[2*i];
                coor[2*(itot-1)+1] = coor[2*i+1];
                itriArr[i] = itot;
            }
        }
        npoint = itot;
        for (int i = 0; i < nelem; ++i)
        {
            kmeshc[3*i]   = itriArr[kmeshc[3*i]-1];
            kmeshc[3*i+1] = itriArr[kmeshc[3*i+1]-1];
            kmeshc[3*i+2] = itriArr[kmeshc[3*i+2]-1];
        }
    };

    int ifill = 0;
    while (true) // label 540
    {
        // Inner loop: 3 × (msho16 + msho29 + compact)
        for (int ih = 0; ih < 3; ++ih)
        {
            int leng = nrepos;
            msho16(kmeshc, nelem, npoint, nipnt, iworkPost, ibuurpPost, leng, ctx);
            if (leng > nrepos)
            { ctx.ierror = 1; throw std::runtime_error("msho2d: repositioning array too small (error 903)"); }
            if (ctx.ierror != 0) return;

            int icancel = 0;
            msho29(kmeshc, nelem, npoint, iworkPost, ibuurpPost, nipnt, coor, icancel);
            if (ctx.ierror != 0) return;
            compactMesh();
        }

        // Diagonal improvement
        ifill = 0;
        for (int ielem = nelemi; ielem < nelem; ++ielem)
        {
            int i1 = kmeshc[3*ielem], i2 = kmeshc[3*ielem+1], i3 = kmeshc[3*ielem+2];
            double adis = nodeDistance(i1, i2, coor);
            double bdis = nodeDistance(i2, i3, coor);
            double cdis2 = nodeDistance(i3, i1, coor);
            const int icase = checkTriangleAngles(adis, bdis, cdis2);

            if (icase > 0)
            {
                int jelem = ielem;
                int i4 = 0;
                int j1 = i1, j2 = i2, j3 = i3, j4 = 0;

                if (icase == 1)
                {
                    msho26(kmeshc, nelem, i2, i3, jelem, i4);
                    int ja = 0; msho42(std::span<const int>(kinbnd), lenbnd, i2, i3, ja);
                    if (ja) i4 = 0;
                    if (i4 > 0 && jelem != ielem)
                    {
                        j1 = i1;
                        msho27(coor, i1, i2, i3, i4, ctx);
                        if (j1 != i1) ifill = 1;
                        j1=i1; j2=i2; j3=i3; j4=i4;
                        kmeshc[3*ielem]=j1; kmeshc[3*ielem+1]=j2; kmeshc[3*ielem+2]=j3;
                        kmeshc[3*jelem]=j2; kmeshc[3*jelem+1]=j4; kmeshc[3*jelem+2]=j3;
                    }
                }
                else if (icase == 2)
                {
                    msho26(kmeshc, nelem, i3, i1, jelem, i4);
                    int ja = 0; msho42(std::span<const int>(kinbnd), lenbnd, i3, i1, ja);
                    if (ja) i4 = 0;
                    if (i4 > 0 && jelem != ielem)
                    {
                        j2 = i2;
                        msho27(coor, i2, i3, i1, i4, ctx);
                        if (j2 != i2) ifill = 1;
                        j1=i2; j2=i3; j3=i1; j4=i4;
                        kmeshc[3*ielem]=j1; kmeshc[3*ielem+1]=j2; kmeshc[3*ielem+2]=j3;
                        kmeshc[3*jelem]=j2; kmeshc[3*jelem+1]=j4; kmeshc[3*jelem+2]=j3;
                    }
                }
                else if (icase == 3)
                {
                    msho26(kmeshc, nelem, i1, i2, jelem, i4);
                    int ja = 0; msho42(std::span<const int>(kinbnd), lenbnd, i1, i2, ja);
                    if (ja) i4 = 0;
                    if (i4 > 0 && jelem != ielem)
                    {
                        j3 = i3;
                        msho27(coor, i3, i1, i2, i4, ctx);
                        if (j3 != i3) ifill = 1;
                        j1=i3; j2=i1; j3=i2; j4=i4;
                        kmeshc[3*ielem]=j1; kmeshc[3*ielem+1]=j2; kmeshc[3*ielem+2]=j3;
                        kmeshc[3*jelem]=j2; kmeshc[3*jelem+1]=j4; kmeshc[3*jelem+2]=j3;
                    }
                }
            }
        }
        if (ctx.ierror != 0) return;

        // Compare and save if better
        {
            const double dmin = computeDismin();
            if (dmin > ratiop)
            {
                ratiop = dmin; npntmp = npoint;
                coortmp.resize(2 * npntmp);
                for (int i = 0; i < npoint; ++i)
                { coortmp[2*i] = coor[2*i]; coortmp[2*i+1] = coor[2*i+1]; }
                neltmp = nelem;
                meshtmp.resize(3 * neltmp);
                for (int i = 0; i < nelem; ++i)
                { meshtmp[3*i]=kmeshc[3*i]; meshtmp[3*i+1]=kmeshc[3*i+1]; meshtmp[3*i+2]=kmeshc[3*i+2]; }
            }
        }

        // Check if convergence is too slow → go to label 900
        if ((reposition && (nherha > 10 && ifill == 0)) || nherha > 20) break;

        // Repositioning: increment nherha and call msho16
        // (msho18 Laplacian is commented out in the Fortran and not called)
        ++nherha;
        {
            int leng = nrepos;
            msho16(kmeshc, nelem, npoint, nipnt, iworkPost, ibuurpPost, leng, ctx);
            if (ctx.ierror != 0) return;
        }

        // Fortran redivision loop (msho30 + goto 300): if an interior
        // triangle is invalid/too distorted, rebuild the local front and
        // restart the advancing-front phase.
        bool restartedFromRedivision = false;
        for (int ielem = nelemi; ielem < nelem; ++ielem)
        {
            const int i1 = kmeshc[3 * ielem];
            const int i2 = kmeshc[3 * ielem + 1];
            const int i3 = kmeshc[3 * ielem + 2];

            // Check only fully interior triangles.
            if (i1 <= nipnt || i2 <= nipnt || i3 <= nipnt)
                continue;

            const double surf = triangleArea(coor, i1, i2, i3);

            auto cubeIndexForNode = [&](int node) -> int
            {
                const int stored = icubeSp[node - 1];
                if (stored > 0)
                {
                    const int idx = stored - 1;
                    if (idx >= 0 && idx < ncube) return idx;
                }

                const double x = coor[2 * (node - 1)];
                const double y = coor[2 * (node - 1) + 1];
                int n1c = static_cast<int>((x - xstart) / dism);
                int n2c = static_cast<int>((y - ystart) / dism);
                n1c = std::clamp(n1c, 0, nx - 1);
                n2c = std::clamp(n2c, 0, ny - 1);
                return n1c + n2c * nx;
            };

            const int c1 = cubeIndexForNode(i1);
            const int c2 = cubeIndexForNode(i2);
            const int c3 = cubeIndexForNode(i3);

            const double valare = (chelp[c1] + chelp[c2] + chelp[c3]) / 3.0;
            const bool areaOutOfRange =
                ((surf < 0.1 * valare) || (surf > 10.0 * valare)) && (nherha < 50);

            if (surf < 0.0 || areaOutOfRange)
            {
                int ielemForMsho30 = ielem + 1; // msho30 is 1-based for ielem
                msho30(kmeshc, nelem, ielemForMsho30,
                       kstaplSp, kstap, npoint, itri,
                       nelemi);
                if (ctx.ierror != 0) return;

                nochan = 0;
                done = advanceFront();
                if (!done || ctx.ierror != 0) return;
                verifyConnectionsArray();

                restartedFromRedivision = true;
                break;
            }
        }

        if (restartedFromRedivision)
            continue;

        // Check nodes with fewer than 5 neighbours (ja check)
        int ja = 0;
        for (int i = nipnt + 1; i <= npoint; ++i)
        {
            const int istart = (i == 1) ? 0 : iworkPost[i - 2];
            const int iend   = iworkPost[i - 1];
            if (iend - istart < 5) { ja = 1; break; }
        }

        // Repeat if nherha < 5, or any interior node has < 5 neighbours, or ifill
        if (nherha < 5 || ja == 1 || ifill == 1) continue; // goto 540
        break; // goto 900
    }

    // ---- label 900: check final vs saved, pick best ----
    {
        const double dmin = computeDismin();
        if (dmin < ratiop - 100.0 * eps && !coortmp.empty())
        {
            // Use saved best mesh
            npoint = npntmp;
            for (int i = 0; i < npoint; ++i)
            { coor[2*i] = coortmp[2*i]; coor[2*i+1] = coortmp[2*i+1]; }
            nelem = neltmp;
            for (int i = 0; i < nelem; ++i)
            { kmeshc[3*i]=meshtmp[3*i]; kmeshc[3*i+1]=meshtmp[3*i+1]; kmeshc[3*i+2]=meshtmp[3*i+2]; }
        }
    }

    // ---- Back-transform coordinates ----
    for (int i = 0; i < npoint; ++i)
    {
        coor[2*i]   = (coor[2*i]   - 1.0) * tran + xmint;
        coor[2*i+1] = (coor[2*i+1] - 1.0) * tran + ymint;
    }
}

} // namespace sepran
