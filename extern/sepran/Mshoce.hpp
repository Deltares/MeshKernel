// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Stichting Deltares
//
// C++20 translation of mshoce.for from the SEPRAN library.

#pragma once

#include "SepranContext.hpp"

#include <span>

namespace sepran
{

/// \brief Main entry point for SEPRAN surface triangulation.
///
/// Translated from Fortran \c mshoce (SEPRAN, Niek Praagman, 1989-2010).
///
/// Orchestrates boundary copying, advancing-front triangulation (msho2d),
/// and optional upgrade to quadratic triangles (msho31).
///
/// \param jnew         i   true = generate new mesh; false = reposition only.
/// \param coor         o   Flat interleaved output node coords (length >= 2*npunt).
/// \param kmeshc       o   Flat element array. Linear: 3 nodes/elem; quadratic: 6.
/// \param inpelm       i   Nodes per element: 3 (linear), 6 or 7 (quadratic).
/// \param nbound       i   Number of nodes described in bcord.
/// \param bcord        i   Flat 2*nbound boundary coord array: [x0,y0, x1,y1, ...].
///                         Internal-curve coords follow immediately at bcord[2*nbound+].
/// \param kbndpt       o   Boundary-point node-number array (length >= nbound).
/// \param boundary     i   Flat 2*numcurvboun array: [curveNr, kbndptStart] per curve.
/// \param numcurvboun  i   Number of curves in the boundary description.
/// \param npoint       o   Number of nodes in the output mesh.
/// \param nelem        o   Number of elements in the output mesh.
/// \param holeinfo     i   Flat 2*(nholes+2) hole information array.
/// \param nholes       i   Number of holes.
/// \param ncoar        i   Number of internal coarseness-control points.
/// \param coar         i   Flat 3*ncoar array: [x, y, coarseness] per point.
/// \param userpoints   i   User sequence numbers for the ncoar special points.
/// \param isurnr       i   Surface sequence number (used for diagnostics).
/// \param numextcurves i   Number of extra internal curves.
/// \param numnodextcurvs i Node counts per internal curve (length numextcurves).
/// \param curvenumbers i   User curve numbers (length numextcurves).
/// \param rinput       i   Real input array with coarseness/user-data.
/// \param nuspnt       i   Number of user-prescribed nodes.
/// \param ndim         i   Spatial dimension of the original problem.
/// \param ctx          i/o SEPRAN error/constant context.
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
            SepranContext&          ctx);

} // namespace sepran
