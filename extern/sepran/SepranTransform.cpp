// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Stichting Deltares
//
// C++20 translation of mshtrans2dsur.for from the SEPRAN library.
//
// Original author: Guus Segal, 2009.

#include "SepranTransform.hpp"

#include <algorithm>
#include <limits>

namespace sepran
{

TransformParams transform2DSurface(std::span<double> coor,
                                   int npoint,
                                   std::span<double> coar,
                                   int ncoar,
                                   std::span<const int> curves,
                                   int ncurvs,
                                   std::span<double> cocurvs,
                                   std::span<double> userco,
                                   int nuspnt,
                                   SepranContext& ctx)
{
    if (ctx.ierror != 0)
        return {};

    // --- Find bounding box of all boundary nodes ---
    double xmint =  ctx.rinfin;
    double xmax  = -ctx.rinfin;
    double ymint =  ctx.rinfin;
    double ymax  = -ctx.rinfin;

    for (int i = 0; i < npoint; ++i)
    {
        xmint = std::min(xmint, coor[2 * i]);
        xmax  = std::max(xmax,  coor[2 * i]);
        ymint = std::min(ymint, coor[2 * i + 1]);
        ymax  = std::max(ymax,  coor[2 * i + 1]);
    }

    // --- Choose uniform scale factor ---
    const double tran = (xmax - xmint > ymax - ymint) ? (xmax - xmint) : (ymax - ymint);

    // --- Scale boundary node coordinates ---
    for (int i = 0; i < npoint; ++i)
    {
        coor[2 * i]     = 1.0 + (coor[2 * i]     - xmint) / tran;
        coor[2 * i + 1] = 1.0 + (coor[2 * i + 1] - ymint) / tran;
    }

    // --- Scale user points ---
    for (int i = 0; i < nuspnt; ++i)
    {
        userco[2 * i]     = 1.0 + (userco[2 * i]     - xmint) / tran;
        userco[2 * i + 1] = 1.0 + (userco[2 * i + 1] - ymint) / tran;
    }

    // --- Scale internal coarseness points ---
    // coar layout: [x0,y0,c0, x1,y1,c1, ...] (0-based)
    for (int i = 0; i < ncoar; ++i)
    {
        coar[3 * i]     = 1.0 + (coar[3 * i]     - xmint) / tran;
        coar[3 * i + 1] = 1.0 + (coar[3 * i + 1] - ymint) / tran;
        coar[3 * i + 2] /= tran; // scale the coarseness value itself
    }

    // --- Scale internal curve nodes ---
    // curves[i] = node count of curve i (0-based).
    // cocurvs stores all curve nodes concatenated, interleaved x,y.
    if (ncurvs > 0)
    {
        int offset = 0;
        for (int i = 0; i < ncurvs; ++i)
        {
            const int nnodes = curves[i];
            for (int k = 0; k < nnodes; ++k)
            {
                cocurvs[2 * (offset + k)]     = 1.0 + (cocurvs[2 * (offset + k)]     - xmint) / tran;
                cocurvs[2 * (offset + k) + 1] = 1.0 + (cocurvs[2 * (offset + k) + 1] - ymint) / tran;
            }
            offset += nnodes;
        }
    }

    return {xmint, ymint, tran};
}

void reverseTransform2DSurface(std::span<double> coor, int npoint, const TransformParams& tp)
{
    for (int i = 0; i < npoint; ++i)
    {
        coor[2 * i]     = tp.xmint + (coor[2 * i]     - 1.0) * tp.tran;
        coor[2 * i + 1] = tp.ymint + (coor[2 * i + 1] - 1.0) * tp.tran;
    }
}

} // namespace sepran
