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
// C++20 translation of mshcopyboun.for and mshchkstapl.for from the SEPRAN library.

#pragma once

#include "SepranContext.hpp"

#include <span>
#include <string_view>

namespace sepran
{

/// \brief Result of copyBoundary.
struct BoundaryCopyResult
{
    int jpnt;   ///< Last node number used (1-based).
    int inside; ///< Count of internal regions: -1 for outer only, 0+ with holes.
    int nbn;    ///< Number of boundary edges written into kbound.
};

/// \brief Copy boundary coordinates into the node coordinate array.
///
/// Translated from Fortran \c mshcopyboun (SEPRAN, Guus Segal, 1999-2001).
///
/// Fills \p coor (first nbn boundary nodes), \p kbndpt (boundary node numbers),
/// and \p kbound (boundary edge connectivity).
///
/// Array index convention (same as Fortran):
/// - \p bcord and \p coor are stored interleaved: [x0,y0, x1,y1, ...] (0-based node access).
/// - Node number values stored in \p kbndpt and \p kbound are **1-based** (matching the
///   convention used throughout the SEPRAN internal arrays).
///
/// \param nbound       Number of nodes in bcord (boundary description).
/// \param bcord        Flat boundary coordinate array (length 2*nbound), interleaved x,y.
/// \param coor         Output node coordinate array (pre-allocated, length >= 2*nbound).
/// \param kbndpt       Output: node number for each boundary description point (length nbound).
/// \param kbound       Output: flat edge-node pairs, 1-based (length >= 2*nbound).
/// \param fillKbound   If true, fill kbndpt and kbound; otherwise only fill coor.
/// \param ctx          SEPRAN context.
/// \return BoundaryCopyResult with jpnt, inside, nbn.
BoundaryCopyResult copyBoundary(int nbound,
                                std::span<const double> bcord,
                                std::span<double> coor,
                                std::span<int> kbndpt,
                                std::span<int> kbound,
                                bool fillKbound,
                                const SepranContext& ctx);

/// \brief Consistency check on the kstapl (staple) array.
///
/// Translated from Fortran \c mshchkstapl (SEPRAN, Guus Segal, 2003).
///
/// Verifies that every node on the left-hand side of a kstapl pair also appears
/// on the right-hand side.  Prints a diagnostic message if not (debug aid; does
/// not throw).
///
/// \param kstapl       Flat edge pairs, stored as [n1_0, n2_0, n1_1, n2_1, ...].
///                     Node indices are 1-based. Length 2*lenkstapl.
/// \param lenkstapl    Number of pairs in kstapl.
/// \param text         Descriptive label printed if an inconsistency is found.
void checkStaple(std::span<const int> kstapl,
                 int lenkstapl,
                 std::string_view text);

} // namespace sepran
