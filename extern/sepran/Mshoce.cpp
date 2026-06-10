// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Stichting Deltares
//
// C++20 translation of mshoce.for from the SEPRAN library
// ("Ingenieursbureau SEPRA", Niek Praagman, 1989-2010).

#include "Mshoce.hpp"

#include "Msho2d.hpp"
#include "SepranBoundary.hpp"
#include "SepranContext.hpp"
#include "SepranQuadratic.hpp"

#include <numeric>
#include <vector>

namespace sepran
{

void mshoce(bool                    jnew,
            std::span<double>       coor,
            std::span<int>          kmeshc,
            int                     inpelm,
            int                     nbound,
            std::span<const double> bcord,
            std::span<int>          kbndpt,
            std::span<const int>    boundary,
            int                     numcurvboun,
            int&                    npoint,
            int&                    nelem,
            std::span<int>          holeinfo,
            int                     nholes,
            int                     ncoar,
            std::span<const double> coar,
            std::span<const int>    userpoints,
            int                     isurnr,
            int                     numextcurves,
            std::span<const int>    numnodextcurvs,
            std::span<const int>    curvenumbers,
            std::span<const double> rinput,
            int                     nuspnt,
            int                     ndim,
            SepranContext&          ctx)
{
    // --- Estimated total node count (= capacity of coor array)
    const int npunt = npoint;

    // --- Allocate local boundary-edge array
    std::vector<int> kbound(2 * npunt, 0);

    // --- Copy boundary coordinates into coor (output node array)
    BoundaryCopyResult bcr = copyBoundary(
        nbound,
        bcord,
        coor,
        kbndpt,
        std::span<int>(kbound),
        true,
        ctx);

    if (ctx.ierror != 0) return;

    const int jpnt   = bcr.jpnt;
    const int inside = bcr.inside;
    int       nbn    = bcr.nbn;

    if (jnew)
    {
        // ---- Generating a new mesh ----------------------------------------
        npoint = jpnt;
        int istep = 1;

        // --- Allocate extquanodes help array
        std::vector<int> extquanodes;

        if (inpelm == 6 || inpelm == 7)
        {
            // Quadratic: double-up boundary connectivity (every other edge)
            istep = 2;
            nbn   = nbn / 2;

            for (int i = 0; i < nbn; ++i)
            {
                kbound[2 * i]     = kbound[4 * i];
                kbound[2 * i + 1] = kbound[4 * i + 3];
            }

            // Compute extra nodes on internal curves
            int nnodes = 0;
            for (int i = 0; i < numextcurves; ++i)
                nnodes += numnodextcurvs[i];
            if (nnodes == 0) nnodes = 1;

            extquanodes.assign(2 * nnodes, 0);
            extquanodes[0] = 0;
        }
        else
        {
            extquanodes.assign(1, 0);
        }

        // --- Build internal-curve coordinate span
        // bcord layout: [2*nbound outer nodes][2*totalInternalNodes internal nodes]
        const std::span<const double> cocurv = bcord.subspan(2 * nbound);

        int nbndpt_loc = nbound;

        // --- Main triangulation
        msho2d(coor, npoint,
               std::span<const int>(kbound.data(), 2 * nbn), nbn,
               kmeshc, nelem,
               boundary, numcurvboun,
               npunt,
               inside,
               holeinfo, nholes,
               true,   // reposition = .true.
               kbndpt, nbndpt_loc,
               coar, ncoar,
               userpoints,
               cocurv,
               numextcurves,
               numnodextcurvs,
               curvenumbers,
               istep,
               std::span<int>(extquanodes),
               isurnr,
               rinput,
               nuspnt, ndim,
               ctx);

        if (ctx.ierror != 0) return;

        // --- Upgrade to quadratic elements if needed
        if (inpelm == 6 || inpelm == 7)
        {
            msho31(coor, npoint, kmeshc, nelem,
                   std::span<const int>(kbound.data(), 2 * nbn), nbn,
                   std::span<const int>(extquanodes),
                   ctx);

            if (ctx.ierror != 0) return;
            // (inpelm==7 centroid placement omitted – not needed for MeshKernel)
        }
    }
    else
    {
        // ---- Repositioning only (jnew == false) -----------------------------
        // mshrep is not translated; only handle quadratic renaming below.

        int nbndpt_loc2 = jpnt;
        (void)nbndpt_loc2;

        if (inpelm == 6 || inpelm == 7)
        {
            int npunt_out = 0;
            msho32(kmeshc, nelem, inpelm, npunt_out);
            if (ctx.ierror != 0) return;

            npoint = (npunt_out > jpnt) ? npunt_out : jpnt;
            nbn   = nbn / 2;

            for (int i = 0; i < nbn; ++i)
            {
                kbound[2 * i]     = kbound[4 * i];
                kbound[2 * i + 1] = kbound[4 * i + 3];
            }

            // Compute extra nodes
            int nnodes = 0;
            for (int i = 0; i < numextcurves; ++i)
                nnodes += numnodextcurvs[i];
            if (nnodes == 0) nnodes = 1;

            std::vector<int> extquanodes2(2 * nnodes, 0);
            extquanodes2[0] = 0;

            msho31(coor, npoint, kmeshc, nelem,
                   std::span<const int>(kbound.data(), 2 * nbn), nbn,
                   std::span<const int>(extquanodes2),
                   ctx);
            if (ctx.ierror != 0) return;
        }
    }

    // --- For quadratic elements: re-copy the boundary coords to get exact positions
    if (ctx.igobs == 0 && inpelm > 4)
    {
        std::vector<int> kboundTmp(2 * npunt, 0);
        copyBoundary(
            nbound, bcord, coor,
            kbndpt,
            std::span<int>(kboundTmp),
            false,
            ctx);
    }
}

} // namespace sepran
