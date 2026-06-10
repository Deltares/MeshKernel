// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Stichting Deltares
//
// C++20 translation of mshcurvinters.for, mshcurvinters1.for, and mshcurvinters2.for
// from the SEPRAN library.

#pragma once

#include "SepranContext.hpp"

#include <span>

namespace sepran
{

/// \brief Check that no two curves in a set mutually intersect.
///
/// Translated from Fortran \c mshcurvinters (SEPRAN, Guus Segal, 2005).
///
/// \param ncurvs   Number of curves.
/// \param iwork    Work array of length ncurvs.
///                 iwork[i] = -1: empty, 0: single, 1: compound.
/// \param xbox     Flat 3D bounding-box array, column-major layout:
///                 xbox[ 0 + 2*(j) + 2*ncurvs*0 ] = min-x of curve j (0-based j)
///                 xbox[ 1 + 2*(j) + 2*ncurvs*0 ] = min-y of curve j
///                 xbox[ 0 + 2*(j) + 2*ncurvs*1 ] = max-x of curve j
///                 xbox[ 1 + 2*(j) + 2*ncurvs*1 ] = max-y of curve j
/// \param icurvs   Accumulated node counts per curve (length ncurvs).
///                 icurvs[j] is the total number of nodes in curves 0..j.
/// \param curves   Flat coordinate array for all curve nodes, interleaved x,y,
///                 1-based node indices (length 2 * total_nodes).
/// \param ctx      SEPRAN context (reads ctx.ierror; sets ctx.ierror on error).
void curvIntersectionCheck(int ncurvs,
                           std::span<const int> iwork,
                           std::span<const double> xbox,
                           std::span<const int> icurvs,
                           std::span<const double> curves,
                           SepranContext& ctx);

/// \brief Check that two specific curves do not intersect edge-by-edge.
///
/// Translated from Fortran \c mshcurvinters1 (SEPRAN, Guus Segal, 2005).
///
/// Throws std::runtime_error if an intersection is detected.
///
/// \param icurnr   0-based index of the first curve.
/// \param jcurnr   0-based index of the second curve.
/// \param icurvs   Accumulated node counts per curve (length >= max(icurnr,jcurnr)+1).
/// \param curves   Flat coord array (interleaved x,y), 1-based node access.
/// \param ctx      SEPRAN context.
void curvIntersectionPair(int icurnr,
                          int jcurnr,
                          std::span<const int> icurvs,
                          std::span<const double> curves,
                          SepranContext& ctx);

/// \brief Check that the boundary edges of a surface do not mutually intersect.
///
/// Translated from Fortran \c mshcurvinters2 (SEPRAN, Guus Segal, 2005).
///
/// Sets ctx.ierror and throws std::runtime_error if self-intersection is found.
///
/// \param kbound   Flat edge-node array (1-based node indices), length 2*nbound.
///                 kbound[2*i]   = first node of edge i (0-based i).
///                 kbound[2*i+1] = second node of edge i.
/// \param nbound   Number of boundary edges.
/// \param coor     Flat node coordinate array, interleaved x,y (1-based node indices).
///                 coor[2*(k-1)]   = x of node k (1-based k).
///                 coor[2*(k-1)+1] = y of node k.
/// \param isurnr   Surface sequence number (used in error message).
/// \param ctx      SEPRAN context.
void boundarySelfIntersectionCheck(std::span<const int> kbound,
                                   int nbound,
                                   std::span<const double> coor,
                                   int isurnr,
                                   SepranContext& ctx);

} // namespace sepran
