//----- AGPL --------------------------------------------------------------------
//
//  Copyright (C)  Stichting Deltares, 2017-2026.
//
//  This programme is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Affero General Public License as
//  published by the Free Software Foundation version 3.
//
//  This programme is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Affero General Public License for more details.
//
//  You should have received a copy of the GNU Affero General Public License
//  along with This programme. If not, see <http://www.gnu.org/licenses/>.
//
//  contact: delft3d.support@deltares.nl
//  Stichting Deltares
//  P.O. Box 177
//  2600 MH Delft, The Netherlands
//
//  All indications and logos of, and references to, "Delft3D",
//  "D-Flow Flexible Mesh" and "Deltares" are registered trademarks of Stichting
//  Deltares, and remain the property of Stichting Deltares. All rights reserved.
//
//-------------------------------------------------------------------------------
// SepranFront.cpp - all translated functions
//
// C++20 translation of msho01, msho09, msho11, msho12, msho13, msho14,
// msho17, msho20, msho21, msho24, msho25, msho26, msho27, msho28, msho29,
// msho30, msho33, msho34, msho39, msho41, msho42, msho75
// from the SEPRAN library ("Ingenieursbureau SEPRA", Niek Praagman, 1989-2010).
//
// TRANSLATION CONVENTIONS (flat-span layout):
//   coor   : interleaved flat; coor[2*(k-1)] = x, coor[2*(k-1)+1] = y (1-based k)
//   kstapl : flat pairs;       kstapl[2*s]=n1, kstapl[2*s+1]=n2  (0-based edge s)
//   kelem  : flat triples;     kelem[3*e], [3*e+1], [3*e+2]       (0-based elem e)
//   itri   : flat;             itri[k-1] for 1-based node k
//   istart : CSR ptr;          istart[i-1] = cumulative neighbour count through node i

#include "SepranFront.hpp"
#include "SepranContext.hpp"
#include "SepranGeometry.hpp"
#include "SepranTopology.hpp"
#include "SepranBoundary.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace sepran
{

// ---------------------------------------------------------------------------
// Forward declaration of internal helper (msho75)
// ---------------------------------------------------------------------------
static void checkSegmentIntersection(double xi, double yi, double xj, double yj,
                                     double x1, double y1, double x2, double y2,
                                     int& ih, const SepranContext& ctx);

// ---------------------------------------------------------------------------
// msho01 – Check boundary elements; fill coarse array for each point
// ---------------------------------------------------------------------------
//
// Parameters (C++ 0-based spans):
//   kbound     i    flat boundary element pairs (kbound[2*e]=n1, [2*e+1]=n2)
//   nbound     i    number of boundary line segments
//   istart    i/o   CSR row-pointer array (length npoint)
//   ibuur     i/o   CSR adjacency list (pre-allocated)
//   coarse     o    coarseness value per node (length npoint)
//   coor       i    node coordinate array (interleaved flat)
//   npoint     i    number of nodes
//   coarsemin  o    minimum boundary coarseness
//   coarsemax  o    maximum boundary coarseness
//   coar      i/o   special-point array coar[3*(i-1)+2] = coarseness (1-based i)
//   ncoar      i    number of special points in coar
// ---------------------------------------------------------------------------
void msho01(std::span<const int>    kbound,
            int                     nbound,
            std::span<int>          istart,
            std::span<int>          ibuur,
            std::span<double>       coarse,
            std::span<const double> coor,
            int                     npoint,
            double&                 coarsemin,
            double&                 coarsemax,
            std::span<double>       coar,
            int                     ncoar,
            SepranContext&          ctx)
{
    coarsemin =  ctx.rinfin;
    coarsemax = -ctx.rinfin;

    // --- Count neighbours: istart[i] = degree of node i+1 (temporarily)
    for (int i = 0; i < npoint; ++i) istart[i] = 0;

    for (int e = 0; e < nbound; ++e)
    {
        const int i1 = kbound[2 * e];
        const int i2 = kbound[2 * e + 1];
        istart[i1 - 1]++;
        istart[i2 - 1]++;
    }

    // --- Convert counts to cumulative CSR pointers
    int itotal = 0;
    for (int i = 0; i < npoint; ++i)
    {
        itotal   += istart[i];
        istart[i] = itotal;
    }

    // --- Clear ibuur then fill via insertNeighbour
    for (int k = 0; k < itotal; ++k) ibuur[k] = 0;

    for (int e = 0; e < nbound; ++e)
    {
        const int i1 = kbound[2 * e];
        const int i2 = kbound[2 * e + 1];
        insertNeighbour(istart, ibuur, i1, i2);
        insertNeighbour(istart, ibuur, i2, i1);
    }

    // --- Check neighbour counts and compute per-node coarseness
    int no = 0;
    for (int i = 1; i <= npoint; ++i)
    {
        const int ni   = istart[i - 1];
        const int ntal = ni - no;

        if (ntal != 0 && ntal != 2)
        {
            ctx.ierror = 1;
            throw std::runtime_error(
                "msho01: boundary point has wrong number of neighbours (error 905)");
        }

        double afstan = 0.0;
        if (ntal > 0)
        {
            for (int k = no; k < no + ntal; ++k)
            {
                const int jknoop = ibuur[k];
                if (jknoop > 0)
                    afstan += nodeDistance(i, jknoop, coor);
            }
            coarse[i - 1] = afstan / ntal;
            coarsemin = std::min(coarsemin, coarse[i - 1]);
            coarsemax = std::max(coarsemax, coarse[i - 1]);
        }
        else
        {
            coarse[i - 1] = 0.0;
        }
        no = ni;
    }

    // --- Back-fill missing coarsenesses in coar
    const double coarsemean = 0.5 * (coarsemin + coarsemax);
    for (int i = 0; i < ncoar; ++i)
    {
        // coar(3,i+1) in Fortran (1-based) = coar[3*i + 2] in C++
        if (coar[3 * i + 2] <= 1.0e-6) coar[3 * i + 2] = coarsemean;
    }
}

// ---------------------------------------------------------------------------
// msho09 – Perpendicular unit vector and midpoint of segment i→j
// ---------------------------------------------------------------------------
//
//   coor  i   node coordinate array (interleaved flat)
//   i     i   first  node (1-based)
//   j     i   second node (1-based)
//   e1    o   x-component of unit perpendicular vector
//   e2    o   y-component of unit perpendicular vector
//   xm    o   x-coordinate of midpoint (slightly offset toward i)
//   ym    o   y-coordinate of midpoint
// ---------------------------------------------------------------------------
void msho09(std::span<const double> coor,
            int                     i,
            int                     j,
            double&                 e1,
            double&                 e2,
            double&                 xm,
            double&                 ym)
{
    const double xi = coor[2 * (i - 1)];
    const double yi = coor[2 * (i - 1) + 1];
    const double xj = coor[2 * (j - 1)];
    const double yj = coor[2 * (j - 1) + 1];

    e1 = -(yj - yi);
    e2 =   xj - xi;

    const double eleng = std::sqrt(e1 * e1 + e2 * e2);
    e1 /= eleng;
    e2 /= eleng;

    // Slightly asymmetric midpoint (intentional — matches Fortran exactly)
    xm = 0.5000001234 * xi + 0.4999998766 * xj;
    ym = 0.5000001234 * yi + 0.4999998766 * yj;
}

// ---------------------------------------------------------------------------
// msho11 – Cosine of angle between line (i1→i2) and line (i2→i3)
// ---------------------------------------------------------------------------
//
//   i1, i2, i3  i   node indices (1-based)
//   coor        i   node coordinate array (interleaved flat)
//   angle       o   dot product (cosine); set to -1.5 if points coincide
// ---------------------------------------------------------------------------
void msho11(int                     i1,
            int                     i2,
            int                     i3,
            std::span<const double> coor,
            double&                 angle,
            const SepranContext&    ctx)
{
    const double x1 = coor[2 * (i1 - 1)],     y1 = coor[2 * (i1 - 1) + 1];
    const double x2 = coor[2 * (i2 - 1)],     y2 = coor[2 * (i2 - 1) + 1];
    const double x3 = coor[2 * (i3 - 1)],     y3 = coor[2 * (i3 - 1) + 1];

    // Direction i2→i3
    double e1 = x3 - x2, e2 = y3 - y2;
    double eleng  = std::sqrt(e1 * e1 + e2 * e2);
    double h1     = (x2 + x3) * 0.5;
    double h2     = (y2 + y3) * 0.5;
    double hnorm  = std::sqrt(h1 * h1 + h2 * h2);

    if (eleng < 1.0e2 * ctx.epsmac * hnorm)
    {
        angle = -1.5;
        return;
    }
    e1 /= eleng;  e2 /= eleng;

    // Direction i2→i1
    double e3 = x1 - x2, e4 = y1 - y2;
    eleng = std::sqrt(e3 * e3 + e4 * e4);
    h1    = (x1 + x2) * 0.5;
    h2    = (y1 + y2) * 0.5;
    hnorm = std::sqrt(h1 * h1 + h2 * h2);

    if (eleng < 1.0e2 * ctx.epsmac * hnorm)
    {
        angle = -1.5;
        return;
    }
    e3 /= eleng;  e4 /= eleng;

    angle = e1 * e3 + e2 * e4;
}

// ---------------------------------------------------------------------------
// msho12 – Adjacent-line search and angle computation at endpoints of i1–i2
// ---------------------------------------------------------------------------
//
//   coor    i    node coordinate array (interleaved flat)
//   kstapl  i    boundary edge list (flat pairs, kstap edges)
//   kstap   i    number of edges in kstapl
//   i1      i    first  node of base edge (1-based)
//   i2      i    second node of base edge (1-based)
//   iex1    o    node of best adjacent edge at i1
//   iex2    o    node of best adjacent edge at i2
//   angle1  o    best angle at i1
//   angle2  o    best angle at i2
// ---------------------------------------------------------------------------
void msho12(std::span<const double> coor,
            std::span<const int>    kstapl,
            int                     kstap,
            int                     i1,
            int                     i2,
            int&                    iex1,
            int&                    iex2,
            double&                 angle1,
            double&                 angle2,
            SepranContext&          ctx)
{
    iex1 = 0;  iex2 = 0;
    angle1 = -1.2;  angle2 = -1.2;

    for (int s = 0; s < kstap; ++s)
    {
        const int ii1 = kstapl[2 * s];
        const int ii2 = kstapl[2 * s + 1];

        // Check edge ii1→ii2: does it adjoin i1 (i.e. ii2 == i1)?
        if (i1 == ii2)
        {
            const double surf = triangleArea(coor, ii1, i1, i2);
            double angle;
            if (surf <= 0.0)
                angle = -1.1;
            else
                msho11(ii1, i1, i2, coor, angle, ctx);

            if (angle > angle1)
            {
                iex1   = ii1;
                angle1 = angle;
            }
        }

        // Check edge ii1→ii2: does it adjoin i2 (i.e. ii1 == i2)?
        if (i2 == ii1)
        {
            const double surf = triangleArea(coor, i1, i2, ii2);
            double angle;
            if (surf <= 0.0)
                angle = -1.1;
            else
                msho11(i1, i2, ii2, coor, angle, ctx);

            if (angle > angle2)
            {
                iex2   = ii2;
                angle2 = angle;
            }
        }
    }

    if (iex1 == 0 || iex2 == 0)
    {
        ctx.ierror = 1;
        throw std::runtime_error("msho12: internal error in triangle (error 2433)");
    }
}

// ---------------------------------------------------------------------------
// msho13 – Check whether lines from i1/i2 to new point (xn,yn) cross
//          any existing boundary edge.
//
//   coor    i    node coordinate array (2-component interleaved flat)
//   i1, i2  i    base edge nodes (1-based)
//   kdrie   o    index (1-based) of first crossing edge; 0 = none; -1 = ambiguous
//   kstapl  i    boundary edge list (flat pairs)
//   kstap   i    number of edges in kstapl
//   xn, yn  i    proposed new point coordinates
// ---------------------------------------------------------------------------
void msho13(std::span<const double> coor,
            int                     i1,
            int                     i2,
            int&                    kdrie,
            std::span<const int>    kstapl,
            int                     kstap,
            double                  xn,
            double                  yn,
            const SepranContext&    ctx)
{
    const double eps = 1.0e4 * ctx.epsmac;

    const double xi1 = coor[2 * (i1 - 1)],     yi1 = coor[2 * (i1 - 1) + 1];
    const double xi2 = coor[2 * (i2 - 1)],     yi2 = coor[2 * (i2 - 1) + 1];

    const double xmin = std::min({xi1, xi2, xn});
    const double xmax = std::max({xi1, xi2, xn});
    const double ymin = std::min({yi1, yi2, yn});
    const double ymax = std::max({yi1, yi2, yn});

    kdrie = 0;

    // Helper: update kdrie when a new crossing is found at edge s (0-based)
    // Returns true if an ambiguous double-crossing (-1) was detected.
    auto updateKdrie = [&](int s, double xMid, double yMid) -> bool
    {
        const int oneBasedS = s + 1;
        if (kdrie == 0)
        {
            kdrie = oneBasedS;
            return false;
        }
        if (kdrie > 0 && kdrie != oneBasedS)
        {
            // Double crossing — pick the one whose midpoint is closer to edge i1-i2 midpoint
            const int ipn1_k = kstapl[2 * (kdrie - 1)];
            const int ipn2_k = kstapl[2 * (kdrie - 1) + 1];
            const double xpn1_k = (coor[2*(ipn1_k-1)] + coor[2*(ipn2_k-1)]) * 0.5;
            const double ypn1_k = (coor[2*(ipn1_k-1)+1] + coor[2*(ipn2_k-1)+1]) * 0.5;
            const double xpn2   = (xi1 + xi2) * 0.5;
            const double ypn2   = (yi1 + yi2) * 0.5;
            const double dis1   = (xpn1_k - xpn2)*(xpn1_k - xpn2) + (ypn1_k - ypn2)*(ypn1_k - ypn2);

            const int ipn1_s = kstapl[2 * s];
            const int ipn2_s = kstapl[2 * s + 1];
            const double xpn1_s = (coor[2*(ipn1_s-1)] + coor[2*(ipn2_s-1)]) * 0.5;
            const double ypn1_s = (coor[2*(ipn1_s-1)+1] + coor[2*(ipn2_s-1)+1]) * 0.5;
            const double dis2   = (xpn1_s - xpn2)*(xpn1_s - xpn2) + (ypn1_s - ypn2)*(ypn1_s - ypn2);

            if (dis2 < (1.0 - eps) * dis1)
                kdrie = oneBasedS;
            else if (dis2 < (1.0 + eps) * dis1)
            {
                kdrie = -1;
                return true;  // ambiguous
            }
            // else keep old kdrie
        }
        return false;
        (void)xMid; (void)yMid;
    };

    // Loop over boundary edges, skipping index 0 (the current base edge)
    for (int s = 1; s < kstap; ++s)
    {
        const int ii1 = kstapl[2 * s];
        const int ii2 = kstapl[2 * s + 1];

        const double x1 = coor[2 * (ii1 - 1)],     y1 = coor[2 * (ii1 - 1) + 1];
        const double x2 = coor[2 * (ii2 - 1)],     y2 = coor[2 * (ii2 - 1) + 1];

        const double xmi = std::min(x1, x2), xma = std::max(x1, x2);
        const double ymi = std::min(y1, y2), yma = std::max(y1, y2);

        if (xmi > xmax || xma < xmin || ymi > ymax || yma < ymin)
            continue;

        // --- Check line i1 → new point
        if (i1 != ii1 && i1 != ii2)
        {
            const double dis1 = std::abs(xi1 - x1) + std::abs(yi1 - y1);
            const double dis2 = std::abs(xi1 - x2) + std::abs(yi1 - y2);
            const double h1   = std::abs(xi1) + std::abs(x1) + std::abs(yi1) + std::abs(y1);
            const double h2   = std::abs(xi1) + std::abs(x2) + std::abs(yi1) + std::abs(y2);

            if (dis1 > eps * h1 && dis2 > eps * h2)
            {
                int ih;
                checkSegmentIntersection(xi1, yi1, xn, yn, x1, y1, x2, y2, ih, ctx);
                if (ih == 0)
                    if (updateKdrie(s, 0, 0)) return;
            }
        }

        // --- Check line i2 → new point
        if (i2 != ii1 && i2 != ii2)
        {
            const double dis1 = std::abs(xi2 - x1) + std::abs(yi2 - y1);
            const double dis2 = std::abs(xi2 - x2) + std::abs(yi2 - y2);
            const double h1   = std::abs(xi2) + std::abs(x1) + std::abs(yi2) + std::abs(y1);
            const double h2   = std::abs(xi2) + std::abs(x2) + std::abs(yi2) + std::abs(y2);

            if (dis1 > eps * h1 && dis2 > eps * h2)
            {
                int ih;
                checkSegmentIntersection(xi2, yi2, xn, yn, x1, y1, x2, y2, ih, ctx);
                if (ih == 0)
                    if (updateKdrie(s, 0, 0)) return;
            }
        }

        // --- Check line midpoint(i1,i2) → new point
        double xh = (xi1 + xi2) * 0.5;
        double yh = (yi1 + yi2) * 0.5;
        xh = eps * xn + (1.0 - eps) * xh;
        yh = eps * yn + (1.0 - eps) * yh;

        {
            int ih;
            checkSegmentIntersection(xh, yh, xn, yn, x1, y1, x2, y2, ih, ctx);
            if (ih == 0)
                if (updateKdrie(s, 0, 0)) return;
        }
    }
}

// ---------------------------------------------------------------------------
// msho14 – Find the nearest boundary point to proposed new point (xn, yn)
//
//   coor    i   node coordinate array (interleaved flat)
//   jpn     o   node number of nearest boundary point (1-based; 0 if none)
//   npoint  i   total number of nodes
//   itri    i   per-node boundary membership counter
//   i1, i2  i   base edge nodes (excluded from search)
//   xn, yn  i   proposed new point
//   dista  i/o  current best distance (updated if closer point found)
// ---------------------------------------------------------------------------
void msho14(std::span<const double> coor,
            int&                    jpn,
            int                     npoint,
            std::span<const int>    itri,
            int                     i1,
            int                     i2,
            double                  xn,
            double                  yn,
            double&                 dista)
{
    jpn = 0;
    for (int ih = 1; ih <= npoint; ++ih)
    {
        if (ih == i1 || ih == i2 || itri[ih - 1] == 0) continue;

        const double xi   = coor[2 * (ih - 1)];
        const double yi   = coor[2 * (ih - 1) + 1];
        const double dx   = xn - xi;
        const double dy   = yn - yi;
        const double dist = std::sqrt(dx * dx + dy * dy);

        if (dist < dista)
        {
            jpn   = ih;
            dista = dist;
        }
    }
}

// ---------------------------------------------------------------------------
// msho17 – Insert neighbour ij into the ibuurp array for node ih
//          (uses a different pointer structure than msho02/insertNeighbour)
//
//   ibuurp i/o  neighbour array
//   iwork   i   row-pointer array; iwork[i-1] = end of neighbours for node i
//   ih      i   node whose neighbour list is extended (1-based)
//   ij      i   neighbour node to insert (1-based)
// ---------------------------------------------------------------------------
void msho17(std::span<int>       ibuurp,
            std::span<const int> iwork,
            int                  ih,
            int                  ij)
{
    // jstart will be 1-based into ibuurp
    int jstart = (ih != 1) ? iwork[ih - 2] : 0;

    for (;;)
    {
        ++jstart;
        const int val = ibuurp[jstart - 1];
        if (val == ij) return;   // already present
        if (val != 0)  continue; // slot occupied, advance
        ibuurp[jstart - 1] = ij; // insert
        return;
    }
}

// ---------------------------------------------------------------------------
// msho20 – Cyclic rotation of kstapl: move first edge to last position
//          and update itri accordingly.
//
//   kstapl i/o  boundary edge list (flat pairs)
//   kstap   i   number of edges
//   itri   i/o  per-node boundary membership counter
// ---------------------------------------------------------------------------
void msho20(std::span<int> kstapl,
            int            kstap,
            std::span<int> itri)
{
    const int i1 = kstapl[0];
    const int i2 = kstapl[1];

    // Shift elements left by one position
    for (int ikl = 1; ikl < kstap; ++ikl)
    {
        kstapl[2 * (ikl - 1)]     = kstapl[2 * ikl];
        kstapl[2 * (ikl - 1) + 1] = kstapl[2 * ikl + 1];
    }

    kstapl[2 * (kstap - 1)]     = i1;
    kstapl[2 * (kstap - 1) + 1] = i2;

    itri[i1 - 1]++;
    itri[i2 - 1]++;
}

// ---------------------------------------------------------------------------
// msho21 – Compute reference triangle-area values for each cube cell
//
//   cube    i   coarseness per cell (length ncube)
//   ncube   i   number of cells
//   refvol  o   reference surface value = 0.5 * coarseness^2
// ---------------------------------------------------------------------------
void msho21(std::span<const double> cube,
            int                     ncube,
            std::span<double>       refvol)
{
    for (int i = 0; i < ncube; ++i)
        refvol[i] = 0.5 * cube[i] * cube[i];
}

// ---------------------------------------------------------------------------
// msho24 – Check whether line i1–i2 shares a point with any boundary edge
//
//   kstapl  i   boundary edge list (flat pairs, kstap edges)
//   kstap   i   number of boundary edges
//   coor    i   node coordinate array (interleaved flat)
//   i1, i2  i   line nodes (1-based)
//   icheck  o   -1: only Plaxis endpoints coincide;
//                0: no common point found;
//               >0: 1-based index of intersecting edge
// ---------------------------------------------------------------------------
void msho24(std::span<const int>    kstapl,
            int                     kstap,
            std::span<const double> coor,
            int                     i1,
            int                     i2,
            int&                    icheck,
            const SepranContext&    ctx)
{
    const double eps = 10.0 * ctx.epsmac;

    icheck = -1;

    const double x1 = coor[2 * (i1 - 1)],     y1 = coor[2 * (i1 - 1) + 1];
    const double x2 = coor[2 * (i2 - 1)],     y2 = coor[2 * (i2 - 1) + 1];

    const double xmin = std::min(x1, x2), xmax = std::max(x1, x2);
    const double ymin = std::min(y1, y2), ymax = std::max(y1, y2);

    for (int il = 0; il < kstap; ++il)
    {
        const int ia = kstapl[2 * il];
        const int ib = kstapl[2 * il + 1];

        // Skip if boundary edge shares a node with i1–i2
        if (ia == i1 || ia == i2) continue;
        if (ib == i1 || ib == i2) continue;

        const double xa = coor[2 * (ia - 1)],     ya = coor[2 * (ia - 1) + 1];
        const double xb = coor[2 * (ib - 1)],     yb = coor[2 * (ib - 1) + 1];

        // Coincidence of ia with i1
        double dis = (x1 - xa)*(x1 - xa) + (y1 - ya)*(y1 - ya);
        if (dis < eps) continue;

        // Coincidence of ia with i2 (Plaxis point)
        dis = (x2 - xa)*(x2 - xa) + (y2 - ya)*(y2 - ya);
        if (dis < eps) { icheck = il + 1; return; }

        // Coincidence of ib with i1
        dis = (x1 - xb)*(x1 - xb) + (y1 - yb)*(y1 - yb);
        if (dis < eps) continue;

        // Coincidence of ib with i2 (Plaxis point)
        dis = (x2 - xb)*(x2 - xb) + (y2 - yb)*(y2 - yb);
        if (dis < eps) { icheck = il + 1; return; }

        // Bounding-box pre-check
        if (std::min(xa, xb) > xmax) continue;
        if (std::max(xa, xb) < xmin) continue;
        if (std::min(ya, yb) > ymax) continue;
        if (std::max(ya, yb) < ymin) continue;

        int ih;
        checkSegmentIntersection(xa, ya, xb, yb, x1, y1, x2, y2, ih, ctx);
        if (ih == 0) { icheck = il + 1; return; }
    }

    icheck = 0;
}

// ---------------------------------------------------------------------------
// msho25 – Move the shortest edge in kstapl to position 0
//
//   kstapl i/o  boundary edge list (flat pairs)
//   kstap   i   number of edges
//   coor    i   node coordinate array (interleaved flat)
// ---------------------------------------------------------------------------
void msho25(std::span<int>          kstapl,
            int                     kstap,
            std::span<const double> coor,
            const SepranContext&    ctx)
{
    int    ielem = -1;
    double afref = ctx.rinfin;

    for (int il = 0; il < kstap; ++il)
    {
        const int    ia = kstapl[2 * il];
        const int    ib = kstapl[2 * il + 1];
        const double xa = coor[2 * (ia - 1)],     ya = coor[2 * (ia - 1) + 1];
        const double xb = coor[2 * (ib - 1)],     yb = coor[2 * (ib - 1) + 1];
        const double dx = xb - xa, dy = yb - ya;
        const double afst = dx * dx + dy * dy;

        if (afst < 0.98 * afref)
        {
            afref = afst;
            ielem = il;
        }
    }

    // Swap element ielem with element 0 (only if ielem is not already 0)
    if (ielem > 0)
    {
        const int i1 = kstapl[2 * ielem];
        const int i2 = kstapl[2 * ielem + 1];
        kstapl[2 * ielem]     = kstapl[0];
        kstapl[2 * ielem + 1] = kstapl[1];
        kstapl[0] = i1;
        kstapl[1] = i2;
    }
}

// ---------------------------------------------------------------------------
// msho26 – Find the triangle (other than ielem) that has edge i2–i3 as a side
//
//   kelem  i   element array (flat triples of 1-based node indices)
//   nelem  i   number of elements
//   i2     i   first  node of shared edge (1-based)
//   i3     i   second node of shared edge (1-based)
//   jelem i/o  i: known element containing i2–i3; o: neighbouring element (1-based; 0 if none)
//   i4     o   third node of the neighbouring triangle (1-based; 0 if not found)
// ---------------------------------------------------------------------------
void msho26(std::span<const int> kelem,
            int                  nelem,
            int                  i2,
            int                  i3,
            int&                 jelem,
            int&                 i4)
{
    const int ielem = jelem;
    jelem = 0;
    i4    = 0;

    for (int j = 0; j < nelem; ++j)
    {
        if (j == ielem) continue; // skip known element (0-based)

        const int j1 = kelem[3 * j];
        const int j2 = kelem[3 * j + 1];
        const int j3 = kelem[3 * j + 2];

        if      (j1 == i3 && j2 == i2) { jelem = j; i4 = j3; return; }
        else if (j2 == i3 && j3 == i2) { jelem = j; i4 = j1; return; }
        else if (j3 == i3 && j1 == i2) { jelem = j; i4 = j2; return; }
    }
}

// ---------------------------------------------------------------------------
// msho27 – Select best diagonal in quadrilateral (i1,i2,i3,i4)
//          Current diagonal is i2–i3; may be replaced by i1–i4.
//
//   coor         i/o  node coordinate array (interleaved flat)
//   i1,i2,i3,i4 i/o  quadrilateral nodes; may be reordered on output
// ---------------------------------------------------------------------------
void msho27(std::span<const double> coor,
            int&                    i1,
            int&                    i2,
            int&                    i3,
            int&                    i4,
            const SepranContext&    ctx)
{
    const int j1 = i1, j2 = i2, j3 = i3, j4 = i4;

    const double xa = coor[2 * (i1 - 1)],     ya = coor[2 * (i1 - 1) + 1];
    const double xb = coor[2 * (i2 - 1)],     yb = coor[2 * (i2 - 1) + 1];
    const double xc = coor[2 * (i3 - 1)],     yc = coor[2 * (i3 - 1) + 1];
    const double xd = coor[2 * (i4 - 1)],     yd = coor[2 * (i4 - 1) + 1];

    const double det = (xa - xb) * (ya - yc) - (xa - xc) * (ya - yb);

    if (std::abs(det) < 10.0 * ctx.epsmac)
    {
        // Triangle 1-2-3 degenerate — use diagonal 1–4
        i1 = j2; i2 = j4; i3 = j1; i4 = j3;
    }
    else
    {
        // Circumcentre of triangle i1-i2-i3
        const double xm = ((xa*xa + ya*ya) * (yb - yc) +
                           (xb*xb + yb*yb) * (yc - ya) +
                           (xc*xc + yc*yc) * (ya - yb)) / (2.0 * det);
        const double ym = ((xa*xa + ya*ya) * (xb - xc) +
                           (xb*xb + yb*yb) * (xc - xa) +
                           (xc*xc + yc*yc) * (xa - xb)) / (-2.0 * det);

        const double af1 = (xa - xm)*(xa - xm) + (ya - ym)*(ya - ym);
        const double af2 = (xd - xm)*(xd - xm) + (yd - ym)*(yd - ym);

        if (af2 < af1)
        {
            // d is inside circumcircle — swap diagonal
            i1 = j2; i2 = j4; i3 = j1; i4 = j3;
        }
    }
}

// ---------------------------------------------------------------------------
// msho28 – Check whether any boundary point lies strictly inside triangle
//          (i1, i2, i3) or triangle (i1, i2, (xn,yn)) when i3==0.
//
//   coor    i   node coordinate array (interleaved flat)
//   i1,i2   i   first two nodes (1-based)
//   i3      i   third node (1-based); if 0 use (xn,yn) instead
//   xn, yn  i   coordinates of virtual third vertex (used when i3==0)
//   npoint  i   total node count
//   itri    i   per-node boundary membership counter
//   iperm   o   1 = at least one interior point found; 0 = none
// ---------------------------------------------------------------------------
void msho28(std::span<const double> coor,
            int                     i1,
            int                     i2,
            int                     i3,
            double                  xn,
            double                  yn,
            int                     npoint,
            std::span<const int>    itri,
            int&                    iperm,
            const SepranContext&    ctx)
{
    iperm = 1;
    const double eps = 1.0e-2 * ctx.epsmac;

    const double x1 = coor[2 * (i1 - 1)],     y1 = coor[2 * (i1 - 1) + 1];
    const double x2 = coor[2 * (i2 - 1)],     y2 = coor[2 * (i2 - 1) + 1];
    double x3, y3;
    if (i3 == 0) { x3 = xn; y3 = yn; }
    else         { x3 = coor[2*(i3-1)]; y3 = coor[2*(i3-1)+1]; }

    const double norm = std::abs(x1)+std::abs(x2)+std::abs(x3)+
                        std::abs(y1)+std::abs(y2)+std::abs(y3);
    const double det  = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2);

    for (int i = 1; i <= npoint; ++i)
    {
        if (i == i1 || i == i2 || i == i3 || itri[i - 1] == 0) continue;

        const double xm = coor[2 * (i - 1)];
        const double ym = coor[2 * (i - 1) + 1];

        double opp = (x1*(y2 - ym) + x2*(ym - y1) + xm*(y1 - y2)) / det;
        if (opp < -eps || opp > 1.0 + eps) continue;

        opp = (x2*(y3 - ym) + x3*(ym - y2) + xm*(y2 - y3)) / det;
        if (opp < -eps || opp > 1.0 + eps) continue;

        opp = (x3*(y1 - ym) + x1*(ym - y3) + xm*(y3 - y1)) / det;
        if (opp < -eps || opp > 1.0 + eps) continue;

        // Check for "double" point coinciding with corner i1
        double dist = std::abs(x1 - xm) + std::abs(y1 - ym);
        if (dist < eps * norm) continue;

        dist = std::abs(x2 - xm) + std::abs(y2 - ym);
        if (dist < eps * norm) continue;

        // Point i is strictly inside and not a double — found
        return; // iperm == 1
    }

    iperm = 0; // no interior point found
}

// ---------------------------------------------------------------------------
// msho29 – Remove internal nodes that have exactly 3 or 4 neighbours
//
//   kelem   i/o  element array (flat triples)
//   nelem   i/o  number of elements
//   npoint   i   total node count
//   iwork    i   CSR row-pointer array for ibuurp (length npoint)
//   ibuurp   i   CSR adjacency list (ibuurp[iwork[i-2]..iwork[i-1]-1] for node i)
//   nipnt    i   number of boundary nodes (indices 1..nipnt are fixed)
//   coor     i   node coordinate array (interleaved flat)
//   icancel i/o  counter of cancelled points
// ---------------------------------------------------------------------------
void msho29(std::span<int>          kelem,
            int&                    nelem,
            int                     npoint,
            std::span<const int>    iwork,
            std::span<const int>    ibuurp,
            int                     nipnt,
            std::span<const double> coor,
            int&                    icancel)
{
    icancel = 0;

    for (int i = nipnt + 1; i <= npoint; ++i)
    {
        const int prevPtr = (i > 1) ? iwork[i - 2] : 0;
        const int np      = iwork[i - 1] - prevPtr;

        if (np == 3)
        {
            const int ia = ibuurp[prevPtr];
            const int ib = ibuurp[prevPtr + 1];
            const int ic = ibuurp[prevPtr + 2];

            icancel++;

            // Remove all triangles containing node i; keep the rest
            int nel = 0;
            for (int j = 0; j < nelem; ++j)
            {
                const int e1 = kelem[3*j], e2 = kelem[3*j+1], e3 = kelem[3*j+2];
                if (e1 != i && e2 != i && e3 != i)
                {
                    kelem[3*nel]   = e1;
                    kelem[3*nel+1] = e2;
                    kelem[3*nel+2] = e3;
                    nel++;
                }
            }

            // Add one new triangle ia-ib-ic (with correct orientation)
            const double opp = triangleArea(coor, ia, ib, ic);
            if (opp > 0.0)
            {
                kelem[3*nel]   = ia;
                kelem[3*nel+1] = ib;
                kelem[3*nel+2] = ic;
            }
            else
            {
                kelem[3*nel]   = ia;
                kelem[3*nel+1] = ic;
                kelem[3*nel+2] = ib;
            }
            nelem = nel + 1;
        }
        else if (np == 4)
        {
            // Read the 4 neighbours for the neighbour-check only
            const int ia_nb = ibuurp[prevPtr];
            const int ib_nb = ibuurp[prevPtr + 1];
            const int ic_nb = ibuurp[prevPtr + 2];
            const int id_nb = ibuurp[prevPtr + 3];

            // Skip if any neighbour itself has 3 neighbours, or a lower-index
            // neighbour with 4 (would already have been removed)
            auto neighbourCheck = [&](int nx) -> bool
            {
                if (nx <= nipnt) return false;
                const int prevNx = (nx > 1) ? iwork[nx - 2] : 0;
                const int iex    = iwork[nx - 1] - prevNx;
                return (iex == 3 || (nx < i && iex == 4));
            };

            if (neighbourCheck(ia_nb) || neighbourCheck(ib_nb) ||
                neighbourCheck(ic_nb) || neighbourCheck(id_nb))
                continue;

            // Snapshot before compacting so we can roll back if ring is degenerate.
            const int nelemBefore = nelem;
            std::vector<int> kelemBackup(3 * nelemBefore);
            std::copy_n(kelem.begin(), 3 * nelemBefore, kelemBackup.begin());

            icancel += 1000;

            // Direct Fortran translation: compact kelem removing triangles that
            // contain node i, and simultaneously find ia,ib (first triangle) and
            // ic,id ("opposite" triangle – the one sharing neither ia nor ib).
            int ia = 0, ib = 0, ic = 0, id = 0;
            int iex = 0;
            int nel = 0;

            for (int j = 0; j < nelem; ++j)
            {
                const int e1 = kelem[3*j], e2 = kelem[3*j+1], e3 = kelem[3*j+2];

                if (e1 != i && e2 != i && e3 != i)
                {
                    kelem[3*nel]   = e1;
                    kelem[3*nel+1] = e2;
                    kelem[3*nel+2] = e3;
                    nel++;
                }
                else
                {
                    // The two ring nodes in this triangle (CCW relative order)
                    int n1, n2;
                    if      (e1 == i) { n1 = e2; n2 = e3; }
                    else if (e2 == i) { n1 = e3; n2 = e1; }
                    else              { n1 = e1; n2 = e2; }

                    if (iex == 0)
                    {
                        ia = n1; ib = n2;
                        iex = 1;
                    }
                    else if (iex == 1)
                    {
                        // Accept the "opposite" triangle: neither ring node may
                        // equal ia or ib (direct translation of Fortran condition).
                        if (n1 != ia && n1 != ib && n2 != ia && n2 != ib)
                        {
                            ic = n1; id = n2;
                            iex = 2;
                        }
                    }
                }
            }

            // Degenerate ring: restore and skip (Fortran logs error but we skip).
            if (ia == 0 || ib == 0 || ic == 0 || id == 0)
            {
                std::copy_n(kelemBackup.begin(), 3 * nelemBefore, kelem.begin());
                nelem = nelemBefore;
                icancel -= 1000;
                continue;
            }

            // Choose best diagonal and write two new triangles (Fortran §msho27)
            msho27(coor, ia, ib, id, ic, SepranContext{});

            kelem[3*nel]   = ia;
            kelem[3*nel+1] = ib;
            kelem[3*nel+2] = id;
            nel++;

            kelem[3*nel]   = ib;
            kelem[3*nel+1] = ic;
            kelem[3*nel+2] = id;
            nel++;

            nelem = nel;
        }
    }
}

// ---------------------------------------------------------------------------
// msho30 – Refill kstapl and itri after removing triangle ielem and its
//           neighbours.  Elements nelmfix+1..nelem that share a node with
//           ielem are removed; the remaining ones are kept.
//
//   kelem   i/o  element array (flat triples)
//   nelem   i/o  number of elements
//   ielem    i   1-based index of central triangle to remove
//   kstapl   o   rebuilt boundary edge list (flat pairs)
//   kstap    o   number of edges in kstapl
//   npoint   i   total node count
//   itri     o   per-node boundary membership counter (rebuilt)
//   nelmfix  i   number of elements that are fixed (indices 1..nelmfix)
// ---------------------------------------------------------------------------
void msho30(std::span<int>  kelem,
            int&            nelem,
            int             ielem,
            std::span<int>  kstapl,
            int&            kstap,
            int             npoint,
            std::span<int>  itri,
            int             nelmfix)
{
    const int ia = kelem[3 * (ielem - 1)];
    const int ib = kelem[3 * (ielem - 1) + 1];
    const int ic = kelem[3 * (ielem - 1) + 2];

    kstap   = 0;
    int nextra = nelmfix;

    // Helper: add or cancel an edge in kstapl
    // Adds edge (n1→n2); if reverse (n2→n1) already in kstapl it removes it.
    auto addOrCancelEdge = [&](int n1, int n2)
    {
        // Handle Plaxis negative nodes: add directly
        if (n1 <= -1 && n2 <= -1)
        {
            kstapl[2 * kstap]     = n1;
            kstapl[2 * kstap + 1] = n2;
            kstap++;
            return;
        }

        // Check for reverse edge
        int jtal = 0;
        bool cancelled = false;
        for (int is = 0; is < kstap; ++is)
        {
            const int j1 = kstapl[2 * is];
            const int j2 = kstapl[2 * is + 1];
            if (j1 == n2 && j2 == n1)
            {
                cancelled = true; // remove this edge (don't copy)
            }
            else
            {
                kstapl[2 * jtal]     = j1;
                kstapl[2 * jtal + 1] = j2;
                jtal++;
            }
        }
        kstap = jtal;
        if (!cancelled)
        {
            kstapl[2 * kstap]     = n1;
            kstapl[2 * kstap + 1] = n2;
            kstap++;
        }
    };

    for (int i = nelmfix; i < nelem; ++i) // 0-based; Fortran: nelmfix+1..nelem
    {
        const int i1 = kelem[3 * i];
        const int i2 = kelem[3 * i + 1];
        const int i3 = kelem[3 * i + 2];

        const bool shared = (i1 == ia || i1 == ib || i1 == ic ||
                             i2 == ia || i2 == ib || i2 == ic ||
                             i3 == ia || i3 == ib || i3 == ic);

        if (shared)
        {
            addOrCancelEdge(i1, i2);
            addOrCancelEdge(i2, i3);
            addOrCancelEdge(i3, i1);
        }
        else
        {
            // Keep this element in the nextra zone
            kelem[3 * nextra]     = i1;
            kelem[3 * nextra + 1] = i2;
            kelem[3 * nextra + 2] = i3;
            nextra++;
        }
    }
    nelem = nextra;

    // Rebuild itri from kstapl
    for (int i = 0; i < npoint; ++i) itri[i] = 0;
    for (int s = 0; s < kstap; ++s)
    {
        itri[kstapl[2*s]     - 1]++;
        itri[kstapl[2*s + 1] - 1]++;
    }
}

// ---------------------------------------------------------------------------
// msho33 – Triangle quality ratio: 2 * r_in / r_out
//
//   coor         i   node coordinate array (interleaved flat)
//   i1, i2, i3   i   triangle nodes (1-based)
//   ratio        o   quality measure in [0,1]; equilateral → 1
// ---------------------------------------------------------------------------
void msho33(std::span<const double> coor,
            double&                 ratio,
            int                     i1,
            int                     i2,
            int                     i3)
{
    const double x1 = coor[2*(i1-1)],     y1 = coor[2*(i1-1)+1];
    const double x2 = coor[2*(i2-1)],     y2 = coor[2*(i2-1)+1];
    const double x3 = coor[2*(i3-1)],     y3 = coor[2*(i3-1)+1];

    const double s1 = std::sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
    const double s2 = std::sqrt((x3-x2)*(x3-x2) + (y3-y2)*(y3-y2));
    const double s3 = std::sqrt((x1-x3)*(x1-x3) + (y1-y3)*(y1-y3));

    const double s   = 0.5 * (s1 + s2 + s3);
    const double opp = std::sqrt(s * (s - s1) * (s - s2) * (s - s3));

    const double ro = (s1 * s2 * s3) / (4.0 * opp);
    const double ri = opp / s;

    ratio = 2.0 * ri / ro;
}

// ---------------------------------------------------------------------------
// msho34 – Find the user-point indices ius1, ius2 of the start and end of
//           the boundary curve on which node kpoint lies.
//
//   kpoint       i   node with overly large coarseness (1-based)
//   coor         i   node coordinate array (interleaved flat)
//   ncurvs       i   number of internal curves
//   curves       i   node counts per internal curve
//   kbndpt       i   ordered boundary node array (all curves concatenated)
//   nbndpt       i   length of kbndpt
//   numcurvboun  i   number of boundary curves
//   boundary     i   flat (2×(numcurvboun+...)) start/end indices into kbndpt
//                    boundary[2*(i-1)] = curve number,
//                    boundary[2*(i-1)+1] = start position in kbndpt
//   nbound       i   number of nodes in boundary curves - 1
//   nholes       i   number of holes
//   userco       i   user-specified coordinates (2-component interleaved)
//   nuspnt       i   number of user points
//   coaval       i   coarseness values (length nuspnt+1)
//   ius1         o   start user-point index (1-based)
//   ius2         o   end   user-point index (1-based)
// ---------------------------------------------------------------------------
void msho34(int                     kpoint,
            std::span<const double> coor,
            int                     ncurvs,
            std::span<const int>    curves,
            std::span<const int>    kbndpt,
            int                     nbndpt,
            int                     numcurvboun,
            std::span<const int>    boundary,
            int                     nbound,
            int                     nholes,
            std::span<const double> userco,
            int                     nuspnt,
            std::span<const double> coaval,
            int&                    ius1,
            int&                    ius2)
{
    (void)nbndpt; (void)coaval; // unused in this path

    const double xk = coor[2*(kpoint-1)];
    const double yk = coor[2*(kpoint-1)+1];
    (void)xk; (void)yk;

    bool   found    = false;
    int    bndStart = 0;  // renamed from local 'istart' to avoid confusion
    double xst = 9999.0, yst = 9999.0;
    double xen = 9999.0, yen = 9999.0;

    int n2_prev = 0; // tracks n2 across loops for the internal-curves continuation

    // --- Search through boundary curves
    for (int i = 1; i <= numcurvboun && !found; ++i)
    {
        int n1, n2;
        if (i < numcurvboun)
        {
            // boundary(2,i) and boundary(2,i+1) in Fortran
            n1 = boundary[2*(i-1) + 1];
            n2 = boundary[2*i     + 1];
        }
        else
        {
            n1 = boundary[2*(i-1) + 1];
            n2 = nbound + 1 + nholes;
        }

        if (bndStart == 0) bndStart = kbndpt[n1 - 1]; // kbndpt(n1) 1-based

        if (n2 > n1 + 1 && kbndpt[n2 - 2] == bndStart) // kbndpt(n2-1)
        {
            n2 = n2 - 1;
            bndStart = 0;
        }

        for (int j = n1; j <= n2; ++j)
        {
            if (kpoint == kbndpt[j - 1]) // kbndpt(j) 1-based
            {
                xst = coor[2*(kbndpt[n1-1]-1)];
                yst = coor[2*(kbndpt[n1-1]-1)+1];
                xen = coor[2*(kbndpt[n2-1]-1)];
                yen = coor[2*(kbndpt[n2-1]-1)+1];
                found = true;
                break;
            }
        }
        n2_prev = n2;
    }

    // --- If not found on boundary curves, search internal curves
    for (int i = 1; i <= ncurvs && !found; ++i)
    {
        const int n1 = n2_prev + 1;
        const int n2 = n1 + curves[i - 1] - 1;

        for (int j = n1; j <= n2; ++j)
        {
            if (kpoint == kbndpt[j - 1])
            {
                xst = coor[2*(kbndpt[n1-1]-1)];
                yst = coor[2*(kbndpt[n1-1]-1)+1];
                xen = coor[2*(kbndpt[n2-1]-1)];
                yen = coor[2*(kbndpt[n2-1]-1)+1];
                found = true;
                break;
            }
        }
        n2_prev = n2;
    }

    if (!found)
        throw std::runtime_error("msho34: node kpoint not found on any curve");

    // --- Map start/end coordinates to user-point indices
    ius1 = 0; ius2 = 0;
    for (int i = 1; i <= nuspnt; ++i)
    {
        const double xt = userco[2*(i-1)];
        const double yt = userco[2*(i-1)+1];

        if (std::sqrt((xst-xt)*(xst-xt)+(yst-yt)*(yst-yt)) < 1.0e-8) ius1 = i;
        if (std::sqrt((xen-xt)*(xen-xt)+(yen-yt)*(yen-yt)) < 1.0e-8) ius2 = i;
    }
}

// ---------------------------------------------------------------------------
// msho39 – Check whether coarsenesses cmin and cmax at two endpoints that
//           are dist apart are geometrically compatible.
//
//   cmax      i/o  coarseness at far end; may be reduced if not compatible
//   cmin       i   coarseness at near end
//   dist       i   Euclidean distance between the two points
//   maxratio   i   maximum allowed ratio between successive element lengths
//   iallow     o   1 = compatible; 0 = cmax was adjusted
// ---------------------------------------------------------------------------
void msho39(double&            cmax,
            double             cmin,
            double             dist,
            double             maxratio,
            int&               iallow)
{
    iallow = 1;

    if (dist < cmin)
    {
        iallow = 0;
        cmax   = dist;
    }
    else if (dist > cmin)
    {
        double afst = 0.6  * cmin;
        double som  = afst;
        while (som < dist)
        {
            afst *= maxratio;
            som  += afst;
        }
        if (afst < 0.5 * cmax)
        {
            iallow = 0;
            cmax   = afst;
        }
    }
}

// ---------------------------------------------------------------------------
// msho41 – Compute the maximum allowed coarseness cmax in point B given
//           coarseness cmin in point A at Euclidean distance dist.
//
//   cmin      i    coarseness at A
//   cmax     i/o   i: user-given coarseness at B; o: adjusted max allowed
//   dist      i    Euclidean distance A→B
//   maxratio  i    maximum ratio between consecutive element lengths
// ---------------------------------------------------------------------------
void msho41(double  cmin,
            double& cmax,
            double  dist,
            double  maxratio)
{
    if (dist < cmin)
    {
        cmax = cmin;
    }
    else if (dist > cmin)
    {
        double afst = 0.65 * cmin;
        double som  = afst;
        while (som < dist)
        {
            afst *= maxratio;
            som  += afst;
        }
        if (afst < cmax) cmax = afst;
    }
}

// ---------------------------------------------------------------------------
// msho42 – Check whether the line segment i→j is part of an internal boundary
//
//   kbndpt  i   flat array of internal boundary pairs (node-pair sequence)
//   lenbnd  i   length of kbndpt
//   i       i   first  node (1-based)
//   j       i   second node (1-based)
//   ja      o   1 = i–j is an internal boundary segment; 0 = not
// ---------------------------------------------------------------------------
void msho42(std::span<const int> kbndpt,
            int                  lenbnd,
            int                  i,
            int                  j,
            int&                 ja)
{
    ja = 0;
    for (int k = 0; k + 1 < lenbnd; k += 2)
    {
        const int ia = kbndpt[k];
        const int ib = kbndpt[k + 1];
        if ((ia == i && ib == j) || (ia == j && ib == i))
        {
            ja = 1;
            return;
        }
    }
}

// ---------------------------------------------------------------------------
// msho75 (internal) – Check whether line segment (xi,yi)→(xj,yj) and
//                     segment (x1,y1)→(x2,y2) share an interior point.
//
//   ih  o   0 = common point exists; 1 = no common point
// ---------------------------------------------------------------------------
static void checkSegmentIntersection(double              xi,
                                     double              yi,
                                     double              xj,
                                     double              yj,
                                     double              x1,
                                     double              y1,
                                     double              x2,
                                     double              y2,
                                     int&                ih,
                                     const SepranContext& ctx)
{
    const double eps = 10.0 * ctx.epsmac;
    ih = 1; // default: no intersection

    // Bounding-box pre-check
    const double xmin = std::min(xi, xj), xmax = std::max(xi, xj);
    const double ymin = std::min(yi, yj), ymax = std::max(yi, yj);
    const double xmi  = std::min(x1, x2), xma  = std::max(x1, x2);
    const double ymi  = std::min(y1, y2), yma  = std::max(y1, y2);

    if (xmi > xmax || xma < xmin || ymi > ymax || yma < ymin) return;

    const bool xiNeqXj = std::abs(xi - xj) > eps * (std::abs(xi) + std::abs(xj));
    const bool x1NeqX2 = std::abs(x1 - x2) > eps * (std::abs(x1) + std::abs(x2));

    if (xiNeqXj)
    {
        if (x1NeqX2)
        {
            // General case: neither line is vertical
            const double r1 = (y1*x2 - y2*x1) / (x2 - x1);
            const double r2 = (yi*xj - yj*xi) / (xj - xi);
            const double r3 = (yj - yi) / (xj - xi);
            const double r4 = (y2 - y1) / (x2 - x1);

            if (std::abs(r3 - r4) > eps)
            {
                const double xs = (r1 - r2) / (r3 - r4);
                const bool inIJ = (xi < xs && xj > xs) || (xj < xs && xi > xs) ||
                                  std::abs(xi - xs) < eps || std::abs(xj - xs) < eps;
                const bool in12 = (x1 < xs && x2 > xs) || (x2 < xs && x1 > xs) ||
                                  std::abs(x1 - xs) < eps || std::abs(x2 - xs) < eps;
                if (inIJ && in12) ih = 0;
            }
        }
        else
        {
            // x1 == x2: line 1-2 is vertical, line i-j is not
            if (std::abs(yi - yj) < eps)
            {
                // i-j is horizontal, 1-2 is vertical
                const double xs = x1, ys = yi;
                const bool ok = !((xi < xs && xj < xs) || (xi > xs && xj > xs) ||
                                  (y1 < ys && y2 < ys) || (y1 > ys && y2 > ys));
                if (ok) ih = 0;
            }
            else
            {
                // General vertical 1-2 case
                const double ys = ((yj - yi) * x1 + yi * xj - yj * xi) / (xj - xi);
                const bool in12 = (y1 < ys && y2 > ys) || (y2 < ys && y1 > ys) ||
                                  std::abs(y1 - ys) < eps || std::abs(y2 - ys) < eps;
                const bool inIJ = (yi < ys && yj > ys) || (yj < ys && yi > ys) ||
                                  std::abs(yi - ys) < eps || std::abs(yj - ys) < eps;
                if (in12 && inIJ) ih = 0;
            }
        }
    }
    else
    {
        // xi == xj: line i-j is vertical
        if (x1NeqX2)
        {
            if (std::abs(y1 - y2) < eps)
            {
                // i-j vertical, 1-2 horizontal
                const double xs = xi, ys = y1;
                const bool ok = !((yi < ys && yj < ys) || (yi > ys && yj > ys) ||
                                  (x1 < xs && x2 < xs) || (x1 > xs && x2 > xs));
                if (ok) ih = 0;
            }
            else
            {
                const double ys = ((y2 - y1) * xi + y1*x2 - y2*x1) / (x2 - x1);
                const bool inIJ = (yi < ys && yj > ys) || (yj < ys && yi > ys) ||
                                  std::abs(yi - ys) < eps || std::abs(yj - ys) < eps;
                const bool in12 = (y1 < ys && y2 > ys) || (y2 < ys && y1 > ys) ||
                                  std::abs(y1 - ys) < eps || std::abs(y2 - ys) < eps;
                if (inIJ && in12) ih = 0;
            }
        }
        else
        {
            // Both lines are vertical
            if (std::abs(x1 - xi) < eps * (std::abs(x1) + std::abs(xi))) ih = 0;
        }
    }
}

// ---------------------------------------------------------------------------
// msho04 – Bounding box of all npoint nodes
// ---------------------------------------------------------------------------
void msho04(std::span<const double> coor,
            int npoint,
            double& xmin, double& xmax,
            double& ymin, double& ymax,
            const SepranContext& ctx)
{
    xmax = -ctx.rinfin; xmin = ctx.rinfin;
    ymax = -ctx.rinfin; ymin = ctx.rinfin;

    for (int i = 0; i < npoint; ++i)
    {
        const double x = coor[2 * i];
        const double y = coor[2 * i + 1];
        xmin = std::min(xmin, x); xmax = std::max(xmax, x);
        ymin = std::min(ymin, y); ymax = std::max(ymax, y);
    }
}

// ---------------------------------------------------------------------------
// msho05 – Extreme coarsenesses in dist array
// ---------------------------------------------------------------------------
void msho05(std::span<const double> dist,
            int npoint,
            double& dismin, double& dismax,
            const SepranContext& ctx)
{
    const double eps = 10.0 * ctx.epsmac;
    dismin =  ctx.rinfin;
    dismax = -ctx.rinfin;

    for (int i = 0; i < npoint; ++i)
    {
        const double d = dist[i];
        if (d > eps) { dismin = std::min(dismin, d); dismax = std::max(dismax, d); }
    }
}

// ---------------------------------------------------------------------------
// msho06 – Fill coarseness grid (cube/jcube) from boundary and internal data
// ---------------------------------------------------------------------------
void msho06(int npoint,
            std::span<const double> coor,
            double dist, double xstart, double ystart,
            int nx, int ny,
            std::span<int> icube,
            std::span<const double> chelp,
            std::span<double> cube,
            std::span<int> jcube,
            std::span<const int> kbound,
            int nbound,
            std::span<const double> coar,
            int ncoar,
            int ncurvs,
            std::span<const int> curves,
            std::span<const double> cocurvs,
            const SepranContext& ctx)
{
    const double eps = 10.0 * ctx.epsmac;
    double coa = ctx.rinfin;

    // --- Assign each boundary node to a cube
    for (int i = 0; i < npoint; ++i)
    {
        const double xp = coor[2 * i];
        const double yp = coor[2 * i + 1];
        const int n1 = static_cast<int>((xp - xstart) / dist);
        const int n2 = static_cast<int>((yp - ystart) / dist);
        const int nc = 1 + n1 + n2 * nx;
        icube[i] = nc;
        if (chelp[i] > eps && chelp[i] < coa) coa = chelp[i];
    }

    // --- Initialise cube coarsenesses: average chelp of nodes in each cube
    for (int i = 0; i < nx * ny; ++i) { jcube[i] = 0; cube[i] = 0.0; }

    for (int i = 0; i < nx * ny; ++i)
    {
        int ntal = 0;
        double coarse = 0.0;
        for (int ik = 0; ik < npoint; ++ik)
        {
            if (chelp[ik] > eps && icube[ik] == i + 1)
            {
                ++ntal;
                coarse += chelp[ik];
            }
        }
        if (ntal > 0) { cube[i] = coarse / ntal; jcube[i] = 1; }
    }

    // --- Fill boundary-edge cubes along path of each edge
    for (int e = 0; e < nbound; ++e)
    {
        const int i1 = kbound[2 * e];
        const int i2 = kbound[2 * e + 1];
        const double cgem = (chelp[i1 - 1] + chelp[i2 - 1]) * 0.5;

        const int n1i1 = static_cast<int>((coor[2*(i1-1)]   - xstart) / dist);
        const int n2i1 = static_cast<int>((coor[2*(i1-1)+1] - ystart) / dist);
        const int n1i2 = static_cast<int>((coor[2*(i2-1)]   - xstart) / dist);
        const int n2i2 = static_cast<int>((coor[2*(i2-1)+1] - ystart) / dist);

        for (int n1 = std::min(n1i1,n1i2); n1 <= std::max(n1i1,n1i2); ++n1)
        for (int n2 = std::min(n2i1,n2i2); n2 <= std::max(n2i1,n2i2); ++n2)
        {
            const int nc = 1 + n1 + n2 * nx;
            if (cube[nc-1] < eps) { cube[nc-1] = cgem; jcube[nc-1] = 1; }
        }
    }

    // --- Internal curves
    if (ncurvs > 0)
    {
        int nnodes = 0;
        for (int i = 0; i < ncurvs; ++i)
        {
            for (int j = nnodes; j < nnodes + curves[i] - 1; ++j)
            {
                const double dx = cocurvs[2*(j+1)] - cocurvs[2*j];
                const double dy = cocurvs[2*(j+1)+1] - cocurvs[2*j+1];
                const double afst = std::sqrt(dx*dx + dy*dy);
                const double xp = (cocurvs[2*(j+1)] + cocurvs[2*j]) * 0.5;
                const double yp = (cocurvs[2*(j+1)+1] + cocurvs[2*j+1]) * 0.5;
                const int n1 = static_cast<int>((xp - xstart) / dist);
                const int n2 = static_cast<int>((yp - ystart) / dist);
                const int nc = 1 + n1 + n2 * nx;
                if (jcube[nc-1] == 0)      { cube[nc-1] = afst;                    jcube[nc-1] = 1; }
                else if (jcube[nc-1] > 0)  { cube[nc-1] = (cube[nc-1] + afst) * 0.5; }
            }
            nnodes += curves[i];
        }
    }

    // --- Extra coarseness points
    for (int i = 0; i < ncoar; ++i)
    {
        const double xp = coar[3*i];
        const double yp = coar[3*i+1];
        const int n1 = static_cast<int>((xp - xstart) / dist);
        const int n2 = static_cast<int>((yp - ystart) / dist);
        const int nc = 1 + n1 + n2 * nx;
        cube[nc-1]  = coar[3*i+2];
        jcube[nc-1] = 1;
    }

    // --- Interpolate empty interior cubes
    for (int i = 1; i <= nx; ++i)
    for (int j = 1; j <= ny; ++j)
    {
        const int nrkube = i + (j-1)*nx;
        if (jcube[nrkube-1] != 0 || cube[nrkube-1] >= eps) continue;

        double coxmin = coa, coxmax = coa;
        int ikxmin = 0,   ikxmax = 0;

        for (int ik = (j-1)*nx + 1; ik < nrkube; ++ik)
            if (jcube[ik-1]==1) { ikxmin = ik; coxmin = cube[ik-1]; }
        for (int ik = j*nx; ik > nrkube; --ik)
            if (jcube[ik-1]==1) { ikxmax = ik; coxmax = cube[ik-1]; break; }

        double coymin = coa, coymax = coa;
        int ikymin = 0,   ikymax = 0;

        for (int ik = i; ik < nrkube; ik += nx)
            if (jcube[ik-1]==1) { ikymin = ik; coymin = cube[ik-1]; }
        for (int ik = i + (ny-1)*nx; ik > nrkube; ik -= nx)
            if (jcube[ik-1]==1) { ikymax = ik; coymax = cube[ik-1]; break; }

        double deel = 0.0, cox = 0.0, coy = 0.0;
        const int jblokx = ikxmax * ikxmin;
        if (jblokx != 0)
        {
            if (jblokx > 0)
            {
                const double alpha = std::pow(coxmax/coxmin, 1.0/(ikxmax-ikxmin));
                cox = coxmin * std::pow(alpha, static_cast<double>(nrkube-ikxmin));
            }
            else
            {
                const double alpha = (nrkube-ikxmin) * (1.0/(ikxmax-ikxmin));
                cox = coxmin + alpha * (coxmax-coxmin);
            }
            ++deel;
            cube[nrkube-1] = cox;
        }

        const int jbloky = ikymax * ikymin;
        if (jbloky != 0)
        {
            if (jbloky > 0)
            {
                const double alpha = std::pow(coymax/coymin,
                                              1.0/(static_cast<double>(ikymax-ikymin)/nx));
                coy = coymin * std::pow(alpha, static_cast<double>(nrkube-ikymin)/nx);
            }
            else
            {
                const double alpha = (nrkube-ikymin) * (1.0/(ikymax-ikymin));
                coy = coymin + alpha * (coymax-coymin);
            }
            ++deel;
            cube[nrkube-1] = (cube[nrkube-1] + coy) / deel;
        }

        if (deel > 0.5) jcube[nrkube-1] = 2;
    }

    // --- Fill empty cubes with a background value
    double valmin = ctx.rinfin, valmax = 0.0;
    for (int i = 0; i < nx*ny; ++i)
    {
        if (cube[i] > eps) { valmin = std::min(valmin,cube[i]); valmax = std::max(valmax,cube[i]); }
    }
    coa = (valmax + 3.0*valmin) / 4.0;
    for (int i = 0; i < nx*ny; ++i)
        if (cube[i] < valmin) cube[i] = coa;

    // --- Clamp: no cube may exceed 2.75× its minimum neighbour
    bool changed = true;
    while (changed && nx > 2 && ny > 2)
    {
        changed = false;
        for (int i = nx-2; i >= 1; --i)
        for (int j = ny-2; j >= 1; --j)
        {
            const int nrkube = i + j*nx; // 0-based
            if (jcube[nrkube] != 2) continue;
            const double vmin = std::min({cube[nrkube-nx], cube[nrkube-1],
                                          cube[nrkube+1],  cube[nrkube+nx]});
            if (cube[nrkube] > 2.75*vmin) { cube[nrkube] = 2.75*vmin; changed = true; }
        }
    }
}

// ---------------------------------------------------------------------------
// msho07 – Point-in-polygon test via ray casting
// ---------------------------------------------------------------------------
void msho07(double xcub, double ycub,
            double xmini,
            std::span<const double> coor,
            std::span<const int> kbound,
            int nbound,
            int& ja,
            const SepranContext& ctx)
{
    const double xleft = xmini - 10.0;
    const double yleft = ycub;
    int isnij = 0;
    ja = 0;

    for (int i = 0; i < nbound; ++i)
    {
        const int i1 = kbound[2 * i];
        const int i2 = kbound[2 * i + 1];

        const double x1 = coor[2*(i1-1)],     y1 = coor[2*(i1-1)+1];
        const double x2 = coor[2*(i2-1)],     y2 = coor[2*(i2-1)+1];

        const double xmin = std::min(x1,x2);
        const double ymin = std::min(y1,y2);
        const double ymax = std::max(y1,y2);

        if (xcub < xmin) continue;
        if (ycub < ymin) continue;
        if (ycub > ymax) continue;

        if (!segmentsIntersect(x1, y1, x2, y2, xleft, yleft, xcub, ycub, ctx))
            continue;

        ++isnij;
    }

    if ((isnij % 2) != 0) ja = 1;
}

// ---------------------------------------------------------------------------
// msho08 – Check and fix orientation of front edges
// ---------------------------------------------------------------------------
void msho08(std::span<int> kstapl,
            int kstap,
            std::span<const double> coor,
            double xstart, double dismin,
            std::span<int> holeinfo,
            int nholes,
            bool check,
            SepranContext& ctx)
{
    constexpr double eps = 1.0e-5;
    int icount = 0;

    for (int i = 0; i < kstap; ++i)
    {
        const int i1 = kstapl[2 * i];
        const int i2 = kstapl[2 * i + 1];

        kstapl[2 * kstap + i] = 0;

        double e1, e2, xm, ym;
        msho09(coor, i1, i2, e1, e2, xm, ym);

        double xn = xm + eps * e1 * dismin;
        double yn = ym + eps * e2 * dismin;

        int ja = 0;
        msho07(xn, yn, xstart, coor, kstapl.subspan(0, 2 * kstap), kstap, ja, ctx);

        if (ja == 0)
        {
            if (check)
            {
                ctx.ierror = 1;
                throw std::runtime_error(
                    "msho08: boundary edge (" + std::to_string(i1) + "," +
                    std::to_string(i2) + ") has wrong orientation (error 2434)");
            }
            std::swap(kstapl[2 * i], kstapl[2 * i + 1]);
        }

        if (icount < nholes && holeinfo[2 * icount + 1] == i1)  // holeinfo(2,icount) 1-based
        {
            if (ja == 0) holeinfo[2 * icount] = -holeinfo[2 * icount]; // negate
            if (icount < nholes) ++icount;
        }

        if (ja == 1)
        {
            xn = xm - eps * e1 * dismin;
            yn = ym - eps * e2 * dismin;
            ja = 0;
            msho07(xn, yn, xstart, coor, kstapl.subspan(0, 2 * kstap), kstap, ja, ctx);
            if (ja == 1) kstapl[2 * kstap + i] = 1; // double edge
        }
    }

    // --- Repair connectivity for double (internal) edges
    int iloop   = 0;
    int ichange = 1;

    while (ichange == 1)
    {
        ichange = 0;

        for (int i = 0; i < kstap - 1; ++i)
        {
            if (kstapl[2 * kstap + i] != 0) continue;

            const int i1 = kstapl[2 * i];
            const int i2 = kstapl[2 * i + 1];

            int ih1 = 1, ih2 = 1;

            for (int j = 0; j < kstap; ++j)
            {
                if (kstapl[2 * kstap + j] != 0 || j == i) continue;
                const int j1 = kstapl[2 * j];
                const int j2 = kstapl[2 * j + 1];
                if (j1 == i2) ih2 = 0;
                if (j2 == i1) ih1 = 0;
            }

            if (ih1 + ih2 > 0)
            {
                for (int j = 0; j < kstap; ++j)
                {
                    if (kstapl[2 * kstap + j] != 1) continue;

                    const int j1 = kstapl[2 * j];
                    const int j2 = kstapl[2 * j + 1];

                    if ((j1 == i2 && ih2 > 0) || (j2 == i1 && ih1 > 0))
                    {
                        kstapl[2 * kstap + j] = 0;
                        ih1 = 0; ih2 = 0;
                    }
                    else if ((j1 == i1 && ih1 > 0) || (j2 == i2 && ih2 > 0))
                    {
                        kstapl[2 * kstap + j] = 0;
                        std::swap(kstapl[2 * j], kstapl[2 * j + 1]);
                        ichange = 1;
                        ih1 = 0; ih2 = 0;
                    }
                }
            }
        }

        ++iloop;
        if (iloop > 1000)
        {
            ctx.ierror = 1;
            throw std::runtime_error("msho08: internal error - failed to repair front connectivity (error 1274)");
        }
    }

    checkStaple(kstapl.subspan(0, 2 * kstap), kstap, "Final check in msho08");
}

// ---------------------------------------------------------------------------
// msho15 – Commit new node and triangle; update kstapl and itri
// ---------------------------------------------------------------------------
void msho15(std::span<double> coor,
            int& npoint,
            std::span<int> kelem,
            int& nelem,
            std::span<int> kstapl,
            int& kstap,
            std::span<int> itri,
            int i1, int i2,
            double xnx, double yny)
{
    ++npoint;
    coor[2 * (npoint - 1)]     = xnx;
    coor[2 * (npoint - 1) + 1] = yny;

    ++nelem;
    kelem[3 * (nelem - 1)]     = i1;
    kelem[3 * (nelem - 1) + 1] = i2;
    kelem[3 * (nelem - 1) + 2] = npoint;

    // Remove first edge (left-shift by one pair)
    const int twoKstap = 2 * kstap;
    for (int ib = 0; ib < twoKstap - 2; ++ib)
        kstapl[ib] = kstapl[ib + 2];

    const int base = twoKstap - 2;
    ++kstap;
    kstapl[base]     = i1;
    kstapl[base + 1] = npoint;
    kstapl[base + 2] = npoint;
    kstapl[base + 3] = i2;

    itri[npoint - 1] = 2;
    ++itri[i1 - 1];
    ++itri[i2 - 1];
}

// ---------------------------------------------------------------------------
// msho16 – Build full CSR neighbour structure for all mesh nodes
// ---------------------------------------------------------------------------
void msho16(std::span<const int> kelem,
            int nelem, int npoint, int nipnt,
            std::span<int> iwork,
            std::span<int> ibuurp,
            int& leng,
            SepranContext& ctx)
{
    // Count appearances per node in element array
    for (int i = 0; i < npoint; ++i) iwork[i] = 0;

    for (int i = 0; i < nelem * 3; ++i)
        iwork[kelem[i] - 1]++;  // 1-based node index

    // Add 1 for every boundary node
    for (int i = 0; i < nipnt; ++i)
        iwork[i]++;

    // Convert to cumulative CSR pointers
    int isum = 0;
    for (int i = 0; i < npoint; ++i)
    {
        iwork[i] += isum;
        isum      = iwork[i];
    }

    if (isum > leng)
    {
        ctx.ierror = 1;
        throw std::runtime_error("msho16: declared length of ibuurp is too small (error 903)");
    }

    for (int i = 0; i < isum; ++i) ibuurp[i] = 0;
    leng = isum;

    // Fill ibuurp with unique neighbours
    for (int i = 0; i < nelem; ++i)
    {
        const int is = 3 * i;
        const int i1 = kelem[is];
        const int i2 = kelem[is + 1];
        const int i3 = kelem[is + 2];
        msho17(ibuurp, iwork, i1, i2);
        msho17(ibuurp, iwork, i1, i3);
        msho17(ibuurp, iwork, i2, i1);
        msho17(ibuurp, iwork, i2, i3);
        msho17(ibuurp, iwork, i3, i1);
        msho17(ibuurp, iwork, i3, i2);
    }
}

// ---------------------------------------------------------------------------
// msho35 – Place fixed hexagonal clusters for user coarseness points.
//
// Translated from Fortran msho35 (SEPRAN, Niek Praagman, 2005-2008).
// For each special coarseness point in coar[], place the centre node plus
// 6 ring nodes in coor, then add 6 triangles and 6 front edges.
// ---------------------------------------------------------------------------
void msho35(int                     npoint,
            std::span<double>       coor,
            double                  xstart,
            double                  ystart,
            double                  dismax,
            std::span<const double> coar,
            int                     ncoar,
            std::span<int>          icube,
            int                     nx,
            std::span<int>          kelem,
            int&                    nelem,
            std::span<int>          kstapl,
            int&                    kstap,
            std::span<int>          itri,
            int                     isurnr,
            std::span<const int>    userpoints,
            std::span<int>          kbndpt,
            int&                    nbndpt,
            SepranContext&          ctx)
{
    const int npoin = npoint;
    // npoint will be incremented as we add nodes; track via local counters
    int np = npoint + ncoar; // centre nodes are placed at npoin+1 .. npoin+ncoar

    for (int i = 1; i <= ncoar; ++i)
    {
        if (ctx.ierror != 0) return;

        // --- Index of 6 ring nodes for coar-point i (0-based i here)
        const int npoine = npoin + ncoar + (i - 1) * 6;

        const double dist = coar[3 * (i - 1) + 2];
        const double xp   = coar[3 * (i - 1)];
        const double yp   = coar[3 * (i - 1) + 1];
        const double eps  = 1.0e-10 * dist;

        // --- Check point inside domain
        int ja = 0;
        msho07(xp + eps, yp + eps, xstart, coor, kstapl.subspan(0, 2 * kstap), kstap, ja, ctx);

        if (ja == 0)
        {
            ctx.ierror = 1;
            throw std::runtime_error(
                "msho35: user point " + std::to_string(userpoints[i - 1]) +
                " (surface " + std::to_string(isurnr) + ") not inside region (error 1344)");
        }

        // --- Store centre node
        nbndpt++;
        kbndpt[nbndpt - 1] = npoin + i;

        coor[2 * (npoin + i - 1)]     = xp;
        coor[2 * (npoin + i - 1) + 1] = yp;

        const int n1 = static_cast<int>((xp - xstart) / dismax);
        const int n2 = static_cast<int>((yp - ystart) / dismax);
        icube[npoin + i - 1] = 1 + n1 + n2 * nx;

        // --- Six ring points: clockwise hexagon
        const std::array<double, 12> offsets = {
            -0.40, +0.67, +0.40, +0.67, -0.80, 0.0,
            +0.80, 0.0,  -0.40, -0.67, +0.40, -0.67
        };

        for (int k = 0; k < 6; ++k)
        {
            const double xpoin = xp + offsets[2 * k]     * dist;
            const double ypoin = yp + offsets[2 * k + 1] * dist;

            ja = 0;
            msho07(xpoin + eps, ypoin + eps, xstart, coor,
                   kstapl.subspan(0, 2 * kstap), kstap, ja, ctx);

            if (ja == 0)
            {
                ctx.ierror = 1;
                throw std::runtime_error(
                    "msho35: user point " + std::to_string(userpoints[i - 1]) +
                    " (surface " + std::to_string(isurnr) + ") ring node " +
                    std::to_string(k + 1) + " outside region (error 1346)");
            }

            ++np;
            coor[2 * (np - 1)]     = xpoin;
            coor[2 * (np - 1) + 1] = ypoin;
            icube[np - 1]          = 1 + n1 + n2 * nx;
        }

        // --- Add 6 triangles (Fortran kelem is 1-based column-major, same node numbering)
        // Ring order (1-based relative to npoine): 1,3,5,6,4,2
        const std::array<int, 6> ring = { npoine + 1, npoine + 3, npoine + 5,
                                           npoine + 6, npoine + 4, npoine + 2 };
        const int centre = npoin + i;

        for (int k = 0; k < 6; ++k)
        {
            const int r1 = ring[k];
            const int r2 = ring[(k + 1) % 6];
            ++nelem;
            kelem[3 * (nelem - 1)]     = centre;
            kelem[3 * (nelem - 1) + 1] = r1;
            kelem[3 * (nelem - 1) + 2] = r2;
        }

        // --- Adjust boundary array (6 new front edges = ring perimeter)
        for (int k = 0; k < 6; ++k)
        {
            const int r1 = ring[k];
            const int r2 = ring[(k + 1) % 6];
            ++kstap;
            kstapl[2 * (kstap - 1)]     = r1;
            kstapl[2 * (kstap - 1) + 1] = r2;
        }

        for (int k = 0; k < 6; ++k)
            itri[npoine + k] = 2;
    }
}

// ---------------------------------------------------------------------------
// msho36 – Register internal curve nodes into coor, kstapl, kbndpt.
//
// Translated from Fortran msho36 (SEPRAN, Niek Praagman, 2008-2009).
// ---------------------------------------------------------------------------
void msho36(std::span<double>       coor,
            int&                    npoint,
            int                     istep,
            int                     ncurvs,
            std::span<const int>    curves,
            std::span<const double> cocurvs,
            std::span<int>          kstapl,
            int&                    kstap,
            std::span<int>          kbndpt,
            int&                    nbndpt,
            double                  coarsemin,
            std::span<int>          extquanodes,
            std::span<double>       coarse,
            SepranContext&          ctx)
{
    const double eps = 1.0e-6 * coarsemin;
    int nnodes = 0;
    int ibound = 1; // 1-based index into extquanodes (for quadratic)

    for (int i = 0; i < ncurvs; ++i)
    {
        if (ctx.ierror != 0) return;

        // Loop over edges of curve i (node indices are 0-based into cocurvs)
        for (int j = nnodes; j < nnodes + curves[i] - 1; j += istep)
        {
            ++npoint;
            ++nbndpt;
            kbndpt[nbndpt - 1] = npoint;

            coor[2 * (npoint - 1)]     = cocurvs[2 * j];
            coor[2 * (npoint - 1) + 1] = cocurvs[2 * j + 1];

            const double x = cocurvs[2 * j];
            const double y = cocurvs[2 * j + 1];

            const double dx = cocurvs[2 * (j + 1)]     - x;
            const double dy = cocurvs[2 * (j + 1) + 1] - y;
            coarse[npoint - 1] = std::sqrt(dx * dx + dy * dy);

            // Check for coincidence with earlier node
            int npn = 0;
            for (int iext = 0; iext < npoint - 1; ++iext)
            {
                const double ddx = x - coor[2 * iext];
                const double ddy = y - coor[2 * iext + 1];
                if (std::sqrt(ddx * ddx + ddy * ddy) < eps && npn == 0)
                {
                    npn    = iext + 1; // 1-based
                    --npoint;
                    kbndpt[nbndpt - 1] = npn;
                }
            }

            // For quadratic: extra midpoint
            if (istep == 2)
            {
                ++npoint;
                ++nbndpt;
                kbndpt[nbndpt - 1]         = npoint;
                coor[2 * (npoint - 1)]     = cocurvs[2 * (j + 1)];
                coor[2 * (npoint - 1) + 1] = cocurvs[2 * (j + 1) + 1];
            }

            // Augment kstapl (front) in both directions
            ++kstap;
            if (npn > 0)
            {
                kstapl[2 * (kstap - 1)]     = npn;
                kstapl[2 * (kstap - 1) + 1] = npoint + 1;
                ++kstap;
                kstapl[2 * (kstap - 1)]     = npoint + 1;
                kstapl[2 * (kstap - 1) + 1] = npn;

                if (istep == 2)
                {
                    extquanodes[ibound + 1] = npn;
                    extquanodes[ibound + 2] = npoint;
                    extquanodes[ibound + 3] = npoint + 1;
                    ibound += 3;
                }
            }
            else
            {
                kstapl[2 * (kstap - 1)]     = npoint + 1 - istep;
                kstapl[2 * (kstap - 1) + 1] = npoint + 1;
                ++kstap;
                kstapl[2 * (kstap - 1)]     = npoint + 1;
                kstapl[2 * (kstap - 1) + 1] = npoint + 1 - istep;

                if (istep == 2)
                {
                    extquanodes[ibound + 1] = npoint - 1;
                    extquanodes[ibound + 2] = npoint;
                    extquanodes[ibound + 3] = npoint + 1;
                    ibound += 3;
                }
            }
        }

        // --- Last node of curve i
        const int jlast = nnodes + curves[i] - 1; // 0-based
        const double xlast = cocurvs[2 * jlast];
        const double ylast = cocurvs[2 * jlast + 1];

        int npn = 0;
        for (int iext = 0; iext < npoint - 1; ++iext)
        {
            const double ddx = xlast - coor[2 * iext];
            const double ddy = ylast - coor[2 * iext + 1];
            if (std::sqrt(ddx * ddx + ddy * ddy) < eps && npn == 0)
                npn = iext + 1; // 1-based
        }

        if (npn == 0)
        {
            ++npoint;
            ++nbndpt;
            kbndpt[nbndpt - 1]         = npoint;
            coor[2 * (npoint - 1)]     = xlast;
            coor[2 * (npoint - 1) + 1] = ylast;

            // Coarseness = distance from previous to last
            const int jprev = jlast - 1;
            const double dx = xlast - cocurvs[2 * jprev];
            const double dy = ylast - cocurvs[2 * jprev + 1];
            coarse[npoint - 1] = std::sqrt(dx * dx + dy * dy);
        }
        else
        {
            ++nbndpt;
            kbndpt[nbndpt - 1] = npn;

            // Adjust last pair in kstapl: endpoint was written as npoint+1
            // but the real node is npn
            kstapl[2 * (kstap - 1)]     = npn;
            kstapl[2 * (kstap - 2) + 1] = npn;
        }

        nnodes += curves[i];
    }

    if (istep == 2 && !extquanodes.empty())
        extquanodes[0] = ibound;
}

// ---------------------------------------------------------------------------
// msho38 – Check/adjust coarseness smoothness for all boundary and internal
//           curves and user-specified coarseness points.
//
// Translated from Fortran msho38 (SEPRAN, Niek Praagman, 2009).
// ---------------------------------------------------------------------------
void msho38(int                     npoint,
            std::span<const double> coor,
            double                  dist,
            double                  xstart,
            double                  ystart,
            std::span<const int>    kbndpt,
            int                     nbndpt,
            int                     numcurvboun,
            std::span<const int>    boundary,
            int                     nbound,
            int                     nholes,
            int                     nx,
            int                     ny,
            std::span<const int>    icube,
            std::span<const double> coarse,
            std::span<const double> cube,
            std::span<const int>    jcube,
            std::span<const int>    kstapl,
            int                     kstap,
            int                     ncurvs,
            std::span<const int>    curves,
            std::span<const double> cocurv,
            int                     isurnr,
            std::span<const double> userco,
            int                     nuspnt,
            std::span<double>       coaval,
            double                  tran,
            const SepranContext&    ctx)
{
    (void)dist; (void)xstart; (void)ystart;
    (void)nx;   (void)ny;
    (void)icube; (void)cube; (void)jcube;
    (void)ncurvs; (void)cocurv;
    (void)isurnr;

    const double eps = 1.0e-9 * coaval[0];

    // maxratio stored in coaval[nuspnt+1] (0-based, matching coaval(nuspnt+2) Fortran)
    const double maxratio = coaval[nuspnt + 1];
    const double csmall   = [&]()
    {
        double v = ctx.rinfin;
        for (int i = 0; i < nx * ny; ++i)
            if (cube[i] < v && jcube[i] > 0) v = cube[i];
        return v;
    }();

    int iperm = 1;
    (void)iperm;

    for (int i = 0; i < npoint - 1; ++i)
    for (int j = i + 1; j < npoint; ++j)
    {
        if (coarse[i] <= eps || coarse[j] <= eps) continue;

        // Euclidean distance
        const double afst = nodeDistance(i + 1, j + 1, coor); // 1-based

        double cmax, cmin;
        int kpoint;
        if (coarse[i] > coarse[j])
        { cmax = coarse[i]; cmin = coarse[j]; kpoint = i + 1; }
        else
        { cmax = coarse[j]; cmin = coarse[i]; kpoint = j + 1; }

        const double coa = std::abs((cmax - cmin) / (cmax + cmin));

        if (afst < 0.2 * csmall)
        {
            int icheck = 0;
            msho24(kstapl, kstap, coor, i + 1, j + 1, icheck, ctx);
            if (icheck == 0)
            {
                int ius1 = 0, ius2 = 0;
                msho34(kpoint, coor, ncurvs, curves, kbndpt, nbndpt,
                       numcurvboun, boundary, nbound, nholes,
                       userco, nuspnt, coaval, ius1, ius2);
                // warning, no coaval adjustment for too-small distance
            }
        }
        else if (coa > 0.1)
        {
            int iallow = 1;
            msho39(cmax, cmin, afst, maxratio, iallow);

            int icheck = 0;
            if (iallow == 0)
                msho24(kstapl, kstap, coor, i + 1, j + 1, icheck, ctx);

            if (iallow == 0 && icheck == 0)
            {
                iperm = 0;

                msho41(cmin, cmax, afst, maxratio);

                int ius1 = 0, ius2 = 0;
                msho34(kpoint, coor, ncurvs, curves, kbndpt, nbndpt,
                       numcurvboun, boundary, nbound, nholes,
                       userco, nuspnt, coaval, ius1, ius2);

                if (ius1 == 0 || ius2 == 0) continue;

                const double val1 = coaval[0] * coaval[ius1] / tran;
                const double val2 = coaval[0] * coaval[ius2] / tran;

                if (cmax < val1 && cmax < val2)
                {
                    coaval[ius1] = cmax * tran / coaval[0];
                    coaval[ius2] = cmax * tran / coaval[0];
                }
                else
                {
                    const double x1 = userco[2 * (ius1 - 1)];
                    const double y1 = userco[2 * (ius1 - 1) + 1];
                    const double x2 = coor[2 * (kpoint - 1)];
                    const double y2 = coor[2 * (kpoint - 1) + 1];
                    const double x3 = userco[2 * (ius2 - 1)];
                    const double y3 = userco[2 * (ius2 - 1) + 1];

                    const double af1 = std::sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
                    const double af2 = std::sqrt((x2-x3)*(x2-x3)+(y2-y3)*(y2-y3));
                    double dif       = val1 + af1 * (val2 - val1) / (af1 + af2) - cmax;

                    double v1 = val1, v2 = val2;
                    if (v1 - dif < 0.5 * cmax) v1 = dif + 0.5 * cmax;
                    if (v2 - dif < 0.5 * cmax) v2 = dif + 0.5 * cmax;

                    coaval[ius1] = tran / coaval[0] * (v1 - dif);
                    coaval[ius2] = tran / coaval[0] * (v2 - dif);
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// msho40 – Full mesh self-intersection check.
//
// Translated from Fortran msho40 (SEPRAN, Niek Praagman, 2009-2010).
// Returns iperm = 1 if intersections are found, 0 if clean.
// ---------------------------------------------------------------------------
void msho40(std::span<const double> coor,
            std::span<const int>    kelem,
            int npoint, int nelem,
            int& iperm)
{
    constexpr double eps = 1.0e-5;
    iperm = 0;

    // --- Barycentre check: no barycentre of one triangle may lie inside another
    for (int ie = 0; ie < nelem; ++ie)
    {
        const int i1 = kelem[3 * ie];
        const int i2 = kelem[3 * ie + 1];
        const int i3 = kelem[3 * ie + 2];

        const double x1 = coor[2*(i1-1)], y1 = coor[2*(i1-1)+1];
        const double x2 = coor[2*(i2-1)], y2 = coor[2*(i2-1)+1];
        const double x3 = coor[2*(i3-1)], y3 = coor[2*(i3-1)+1];

        const double xm = (x1 + x2 + x3) / 3.0;
        const double ym = (y1 + y2 + y3) / 3.0;

        for (int je = 0; je < nelem; ++je)
        {
            if (je == ie) continue;

            const int j1 = kelem[3 * je];
            const int j2 = kelem[3 * je + 1];
            const int j3 = kelem[3 * je + 2];

            const double a1 = coor[2*(j1-1)], b1 = coor[2*(j1-1)+1];
            const double a2 = coor[2*(j2-1)], b2 = coor[2*(j2-1)+1];
            const double a3 = coor[2*(j3-1)], b3 = coor[2*(j3-1)+1];

            const double det = a1*(b2-b3) + a2*(b3-b1) + a3*(b1-b2);
            if (std::abs(det) < eps * eps) continue;

            const double o1 = (a1*(b2-ym) + a2*(ym-b1) + xm*(b1-b2)) / det;
            if (o1 < -eps || o1 > 1.0 + eps) continue;
            const double o2 = (a2*(b3-ym) + a3*(ym-b2) + xm*(b2-b3)) / det;
            if (o2 < -eps || o2 > 1.0 + eps) continue;
            const double o3 = (a3*(b1-ym) + a1*(ym-b3) + xm*(b3-b1)) / det;
            if (o3 < -eps || o3 > 1.0 + eps) continue;

            // Barycentre is inside triangle je
            iperm = 1;
        }
    }

    if (iperm == 1) return;

    // --- Node check: no node may lie strictly inside a triangle
    for (int i = 0; i < npoint; ++i)
    {
        const double xm = coor[2 * i];
        const double ym = coor[2 * i + 1];

        for (int ie = 0; ie < nelem; ++ie)
        {
            const int i1 = kelem[3 * ie];
            const int i2 = kelem[3 * ie + 1];
            const int i3 = kelem[3 * ie + 2];

            // Node is a vertex of this element – skip
            if (i + 1 == i1 || i + 1 == i2 || i + 1 == i3) continue;

            const double x1 = coor[2*(i1-1)], y1 = coor[2*(i1-1)+1];
            const double x2 = coor[2*(i2-1)], y2 = coor[2*(i2-1)+1];
            const double x3 = coor[2*(i3-1)], y3 = coor[2*(i3-1)+1];

            // Skip if node is numerically identical to a vertex (Plaxis double pts)
            const double d1 = (x1-xm)*(x1-xm)+(y1-ym)*(y1-ym);
            if (d1 < eps*(std::abs(x1)+std::abs(xm)+std::abs(y1)+std::abs(ym))) continue;
            const double d2 = (x2-xm)*(x2-xm)+(y2-ym)*(y2-ym);
            if (d2 < eps*(std::abs(x2)+std::abs(xm)+std::abs(y2)+std::abs(ym))) continue;
            const double d3 = (x3-xm)*(x3-xm)+(y3-ym)*(y3-ym);
            if (d3 < eps*(std::abs(x3)+std::abs(xm)+std::abs(y3)+std::abs(ym))) continue;

            const double det = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2);
            if (std::abs(det) < eps * eps) continue;

            const double o1 = (x1*(y2-ym) + x2*(ym-y1) + xm*(y1-y2)) / det;
            if (o1 < -eps || o1 > 1.0 + eps) continue;
            const double o2 = (x2*(y3-ym) + x3*(ym-y2) + xm*(y2-y3)) / det;
            if (o2 < -eps || o2 > 1.0 + eps) continue;
            const double o3 = (x3*(y1-ym) + x1*(ym-y3) + xm*(y3-y1)) / det;
            if (o3 < -eps || o3 > 1.0 + eps) continue;

            iperm = 1;
        }

        if (iperm == 1) return;
    }
}

} // namespace sepran
