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
// C++20 translation of msho31.for, msho32.for, msh401.for, msh402.for,
// msh403.for, msh406.for, msh416.for from the SEPRAN library.

#pragma once

#include "SepranContext.hpp"

#include <span>

namespace sepran
{

/// \brief Convert linear triangle mesh to quadratic (6-node) triangles.
///
/// Translated from Fortran \c msho31 (SEPRAN, Niek Praagman, 1990-2010).
///
/// New mid-side nodes are added to \p coor; \p kelem is expanded from
/// 3 to 6 nodes per element.  Already-existing mid-side nodes on the
/// boundary (given in \p kbound) and on internal quadratic lines
/// (\p extquanodes) are reused.
///
/// On entry \p kelem holds 3 node-indices per element (flat, 0-based element,
/// 1-based node values).  On exit it holds 6 node-indices per element
/// (must be pre-allocated to at least 6*nelem entries).
///
/// \param coor         i/o  Interleaved coords (extended in place).
/// \param npoint       i/o  Node count (extended).
/// \param kelem        i/o  Element connectivity array.
/// \param nelem        i    Number of linear triangles.
/// \param kbound       i    Flat boundary-edge pairs (1-based), length 2*nbn.
/// \param nbn          i    Number of boundary edges.
/// \param extquanodes  i    Quadratic-internal-line data (from msho36).
/// \param ctx          i/o  SEPRAN context.
void msho31(std::span<double>    coor,
            int&                 npoint,
            std::span<int>       kelem,
            int                  nelem,
            std::span<const int> kbound,
            int                  nbn,
            std::span<const int> extquanodes,
            SepranContext&       ctx);

/// \brief Recover linear triangles from a quadratic element array.
///
/// Translated from Fortran \c msho32 (SEPRAN, Niek Praagman, 1991-2005).
///
/// Squashes \p kelem from \p inpelm nodes per element to 3, retaining only
/// vertices 1, 3, 5 (0-based positions 0, 2, 4) per element.
///
/// \param kelem    i/o  Element array (inpelm*nelem entries; overwritten with 3*nelem).
/// \param nelem    i    Number of elements.
/// \param inpelm   i    Nodes per quadratic element (typically 6).
/// \param npunt    o    Highest node number in the resulting linear mesh.
void msho32(std::span<int> kelem,
            int nelem,
            int inpelm,
            int& npunt);

} // namespace sepran
