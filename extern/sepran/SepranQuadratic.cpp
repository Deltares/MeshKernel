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
//
// C++20 translation of msho31.for, msho32.for, msh401-403/406/416 from SEPRAN.
//
// TRANSLATION CONVENTIONS:
//   coor    : coor[2*(k-1)] = x, coor[2*(k-1)+1] = y   (1-based k)
//   kelem   : flat, kelem[3*e], [3*e+1], [3*e+2]         (0-based e, 1-based nodes)
//   kbound  : kbound[2*e], kbound[2*e+1]                  (0-based e, 1-based nodes)

#include "SepranQuadratic.hpp"
#include "SepranContext.hpp"

#include <algorithm>
#include <cassert>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace sepran
{

// ---------------------------------------------------------------------------
// Internal helpers (translated from msh402, msh403, msh406, msh416)
// ---------------------------------------------------------------------------

/// msh416 – Diagonal endpoints of quadrilateral i1-i2-i3-i4.
/// Returns (ie,it) where ie=min, it=max of the chosen diagonal.
static void msh416(int i1, int i2, int i3, int i4, int& ie, int& it)
{
    const int im = std::min({i1, i2, i3, i4});
    if (im == i1 || im == i3) { ie = std::min(i1,i3); it = std::max(i1,i3); }
    else                       { ie = std::min(i2,i4); it = std::max(i2,i4); }
}

/// msh402 – Store the pair (ilow,ihigh) in the neighbour arrays.
///
/// istart[ihigh-1] is a count of neighbours stored so far when ja==1
/// (the "fill" phase); when ja==0 (the "query" phase after compaction),
/// istart holds CSR start addresses.
///
/// kelemh[0] is the count of overflow entries; overflow pairs are stored
/// from kelemh[2] onward as (ihigh, ilow) pairs.
static void msh402(int npoint,
                   int ibrp,
                   std::span<int> ibrpnt,
                   std::span<int> istart,
                   int i1, int i2,
                   std::span<int> kelemh,
                   SepranContext& ctx)
{
    int ja;
    if (i1 < 0) { ja = 0; i1 = -i1; }
    else          { ja = 1; }

    // Determine low and high node numbers
    int ilow, ihigh;
    if (i1 > i2)
    {
        if (i1 > npoint) { ilow = i1; ihigh = i2; }
        else              { ilow = i2; ihigh = i1; }
    }
    else
    {
        if (i2 > npoint) { ilow = i2; ihigh = i1; }
        else              { ilow = i1; ihigh = i2; }
    }

    int jstart, jstend;
    if (ja == 1)
    {
        jstart = ibrp * (ihigh - 1);      // 0-based index
        jstend = jstart + ibrp - 1;
    }
    else
    {
        jstart = istart[ihigh - 1] - 1;   // 0-based
        jstend = istart[ihigh]     - 2;   // 0-based
    }

    for (;;)
    {
        const int ih = ibrpnt[jstart];

        if (ih == 0)
        {
            // empty slot – store here
            ibrpnt[jstart] = ilow;
            if (ja == 1) istart[ihigh - 1]++;
            return;
        }

        if (ih == ilow) return;  // already stored

        ++jstart;

        if (jstart > jstend)
        {
            if (ja == 0)
            {
                ctx.ierror = 1;
                throw std::runtime_error("msh402: no place left for neighbour (error 1274)");
            }

            // Overflow: use kelemh
            int j = kelemh[0];
            // Check if already in kelemh
            for (int il = 0; il < j; ++il)
            {
                if (std::abs(kelemh[2*(il+1)]   - ihigh) == 0 &&
                    std::abs(kelemh[2*(il+1)+1] - ilow ) == 0)
                    return;
            }
            kelemh[2 * (j + 1)]     = ihigh;
            kelemh[2 * (j + 1) + 1] = ilow;
            kelemh[0]                = j + 1;
            istart[ihigh - 1]++;
            return;
        }
    }
}

/// msh403 – Find the mid-side node between i1 and i3 in ibrpnt.
static void msh403(int npoint, int i1, int i3,
                   std::span<const int> ibrpnt,
                   std::span<const int> istart,
                   int& i2,
                   SepranContext& ctx)
{
    i2 = 0;

    int ilow, ihigh;
    if (i1 > i3)
    {
        if (i1 > npoint) { ilow = i1; ihigh = i3; }
        else              { ilow = i3; ihigh = i1; }
    }
    else
    {
        if (i3 > npoint) { ilow = i3; ihigh = i1; }
        else              { ilow = i1; ihigh = i3; }
    }

    int jstart = istart[ihigh - 1] - 1;   // 0-based
    const int jstend = istart[ihigh] - 2; // 0-based

    while (jstart <= jstend)
    {
        if (ibrpnt[jstart] == ilow)
        {
            i2 = ibrpnt[jstart + 1];
            return;
        }
        jstart += 2;
    }

    // Not found
    ctx.ierror = 1;
    throw std::runtime_error("msh403: intermediate neighbour not found (error 949)");
}

/// msh406 – Write 6 node numbers into a quadratic-triangle element.
/// kelem layout: flat, 6 entries per element (0-based element j).
static void msh406(std::span<int> kelem, int j,
                   int i1, int i2, int i3, int i4, int i5, int i6)
{
    const int base = 6 * j;
    kelem[base + 0] = i1;
    kelem[base + 1] = i2;
    kelem[base + 2] = i3;
    kelem[base + 3] = i4;
    kelem[base + 4] = i5;
    kelem[base + 5] = i6;
}

/// msh401 – Fill neighbour arrays for all linear triangles (ishape=3).
static void msh401(int npoint,
                   std::span<const int> kmeshc,
                   int nelem,
                   std::span<int> istart,
                   int ibrp,
                   std::span<int> ibrpnt,
                   std::span<int> kelemh,
                   SepranContext& ctx)
{
    for (int ielem = 0; ielem < nelem; ++ielem)
    {
        const int base = 3 * ielem;
        const int n1 = kmeshc[base];
        const int n2 = kmeshc[base + 1];
        const int n3 = kmeshc[base + 2];

        msh402(npoint, ibrp, ibrpnt, istart, n1, n2, kelemh, ctx);
        msh402(npoint, ibrp, ibrpnt, istart, n2, n3, kelemh, ctx);
        msh402(npoint, ibrp, ibrpnt, istart, n3, n1, kelemh, ctx);
        if (ctx.ierror != 0) return;
    }
}

// ---------------------------------------------------------------------------
// msho31 – Linear to quadratic triangles
// ---------------------------------------------------------------------------
void msho31(std::span<double>    coor,
            int&                 npoint,
            std::span<int>       kelem,
            int                  nelem,
            std::span<const int> kbound,
            int                  nbn,
            std::span<const int> extquanodes,
            SepranContext&       ctx)
{
    const int ibrp   = 10;
    const int ibrlen = ibrp * (npoint + 1);

    std::vector<int> istart (npoint + 2, 0);
    std::vector<int> ibrpnt (ibrlen,     0);
    std::vector<int> kelemh (npoint + 2, 0);

    // --- Fill neighbour graph
    msh401(npoint, kelem, nelem, istart, ibrp, ibrpnt, kelemh, ctx);
    if (ctx.ierror != 0) return;

    // --- Compact: collect all non-zero neighbours into a contiguous prefix
    int nummer = 0;
    for (int i = 0; i < npoint; ++i)
    {
        const int jstart = ibrp * i;
        const int jstend = jstart + ibrp;
        for (int ii = jstart; ii < jstend; ++ii)
        {
            if (ibrpnt[ii] != 0)
            {
                ibrpnt[nummer] = ibrpnt[ii];
                ++nummer;
            }
        }
    }
    const int itotalCompact = nummer;

    if (ibrlen < 2 * itotalCompact)
    {
        ctx.ierror = 1;
        throw std::runtime_error(
            "msho31: not enough space for neighbours (error 1275). Need " +
            std::to_string(2 * itotalCompact) + ", have " + std::to_string(ibrlen));
    }

    // --- Shift to upper half of ibrpnt and rebuild istart as CSR pointers
    int npt    = ibrlen - 1;  // 0-based end
    nummer    = itotalCompact - 1;
    int iextra = 0;

    for (int i = npoint - 1; i >= 0; --i)
    {
        int nbr = istart[i];

        if (nbr > ibrp)
        {
            iextra = 1;
            for (int ib = 0; ib < nbr - ibrp; ++ib)
            {
                ibrpnt[npt] = 0;
                --npt;
            }
            nbr = ibrp;
        }

        for (int ib = 0; ib < nbr; ++ib)
        {
            ibrpnt[npt] = ibrpnt[nummer];
            --npt;
            --nummer;
        }

        if (npt < nummer)
        {
            ctx.ierror = 1;
            throw std::runtime_error("msho31: too few positions in ibrpnt (error 1274)");
        }
    }

    // Account for extra neighbours in kelemh
    int itotal = itotalCompact;
    if (kelemh[0] != 0)
    {
        if (iextra == 0)
        {
            ctx.ierror = 1;
            throw std::runtime_error("msho31: extra points not recognised (error 1274)");
        }
        itotal = kelemh[0] + itotalCompact;
        for (int i = 0; i < kelemh[0]; ++i)
        {
            const int ia = -kelemh[2*(i+1)];     // negative flag for msh402
            const int ib =  kelemh[2*(i+1) + 1];
            msh402(npoint, ibrp, ibrpnt, istart, ia, ib, kelemh, ctx);
            if (ctx.ierror != 0) return;
        }
    }

    if (ibrlen < 2 * itotal)
    {
        ctx.ierror = 1;
        throw std::runtime_error("msho31: not enough space for neighbours (error 1275)");
    }

    // --- Convert istart to proper CSR pointers
    itotal = istart[npoint] - 1; // total edges
    int jstend = 0;
    for (int i = 0; i < npoint; ++i)
    {
        const int jstart = jstend + 1;
        jstend += istart[i];
        istart[i] = jstart;
    }
    istart[npoint] = jstend + 1;

    // --- Expand: shift pairs to make room for new-node slot
    //     ibrpnt[2*k-2] = neighbour, ibrpnt[2*k-1] = midnode (initially 0)
    for (int i = itotal; i >= 1; --i)
    {
        ibrpnt[2 * i - 1] = 0;
        ibrpnt[2 * i - 2] = ibrpnt[i - 1 + npt + 1]; // reindex from compacted region
    }

    // Adjust istart to 2-based offsets
    for (int i = 1; i <= npoint; ++i)
        istart[i] = 2 * istart[i] - 1;  // 1-based Fortran style; adjust for 0-based below

    // Actually, let's do this cleanly: offset all istart so [i-1] points correctly
    // istart is still 1-based internally; convert to 0-based access:
    // We'll use a lambda to abstract the conversion.

    // --- Reuse existing mid-side nodes on boundary
    for (int i = 0; i < nbn; ++i)
    {
        int ip1 = kbound[2 * i];
        int ip2 = kbound[2 * i + 1];
        int ip3 = ip1 + 1;  // pre-existing midpoint assumed to be ip1+1

        // larger node is "ihigh"
        if (ip1 < ip2) { std::swap(ip1, ip2); }

        const int js  = istart[ip1 - 1] - 1;   // 0-based
        const int jse = istart[ip1]     - 2;    // 0-based
        for (int ii = js; ii <= jse; ii += 2)
        {
            if (ibrpnt[ii] == ip2)
            {
                ibrpnt[ii + 1] = ip3;
                break;
            }
        }
    }

    // --- Reuse existing mid-side nodes on internal quadratic lines
    if (!extquanodes.empty() && extquanodes[0] > 1)
    {
        const int ileng = (extquanodes[0] - 1) / 3;
        for (int i = 0; i < ileng; ++i)
        {
            int ip1 = extquanodes[3 * i + 1];   // start node
            int ip2 = extquanodes[3 * i + 3];   // end node
            int ip3 = extquanodes[3 * i + 2];   // mid node

            if (ip1 < ip2) { std::swap(ip1, ip2); }

            const int js  = istart[ip1 - 1] - 1;
            const int jse = istart[ip1]     - 2;
            for (int ii = js; ii <= jse; ii += 2)
            {
                if (ibrpnt[ii] == ip2)
                {
                    ibrpnt[ii + 1] = ip3;
                    break;
                }
            }
        }
    }

    // --- Create new mid-side nodes for all internal edges
    const int npnold = npoint;
    for (int i = 0; i < npnold; ++i)
    {
        const int jstart_i = istart[i] - 1;   // 0-based
        const int jstend_i = istart[i + 1] - 2;
        for (int ii = jstart_i; ii <= jstend_i; ii += 2)
        {
            const int ip2 = ibrpnt[ii];
            const int ip3 = ibrpnt[ii + 1];

            if (ip2 > 0 && ip3 == 0)
            {
                ++npoint;
                ibrpnt[ii + 1] = npoint;

                // Midpoint coordinates (1-based node indices)
                coor[2 * (npoint - 1)]     = (coor[2 * i]           + coor[2 * (ip2 - 1)])     / 2.0;
                coor[2 * (npoint - 1) + 1] = (coor[2 * i + 1]       + coor[2 * (ip2 - 1) + 1]) / 2.0;
            }
        }
    }

    // --- Build quadratic elements (backward so we don't clobber linear data)
    // The input kelem has 3 entries per element; we need 6.
    // Process backward: fill kelem[6*(i-1)..6*i-1] from kelem[3*(i-1)..3*i-1]
    for (int i = nelem - 1; i >= 0; --i)
    {
        int ip[7] = {};
        ip[1] = kelem[3 * i];
        ip[3] = kelem[3 * i + 1];
        ip[5] = kelem[3 * i + 2];

        msh403(npnold, ip[1], ip[3], ibrpnt, istart, ip[2], ctx);
        msh403(npnold, ip[3], ip[5], ibrpnt, istart, ip[4], ctx);
        msh403(npnold, ip[5], ip[1], ibrpnt, istart, ip[6], ctx);
        if (ctx.ierror != 0) return;

        msh406(kelem, i, ip[1], ip[2], ip[3], ip[4], ip[5], ip[6]);
    }
}

// ---------------------------------------------------------------------------
// msho32 – Recover linear triangles from quadratic element array
// ---------------------------------------------------------------------------
void msho32(std::span<int> kelem,
            int nelem,
            int inpelm,
            int& npunt)
{
    for (int i = 0; i < nelem; ++i)
    {
        const int lin = 3 * i;
        const int qua = inpelm * i;
        kelem[lin]     = kelem[qua];
        kelem[lin + 1] = kelem[qua + 2];
        kelem[lin + 2] = kelem[qua + 4];
    }

    npunt = 0;
    for (int i = 0; i < 3 * nelem; ++i)
        if (kelem[i] > npunt) npunt = kelem[i];
}

} // namespace sepran
