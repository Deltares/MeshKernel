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
// C++20 translation of msho2d.for from the SEPRAN library
// ("Ingenieursbureau SEPRA", Niek Praagman, 1989-2011).

#pragma once

#include "SepranContext.hpp"

#include <span>

namespace sepran
{

/// \brief Main advancing-front triangulation driver for a single 2D surface.
///
/// Translated from Fortran \c msho2d (SEPRAN, Niek Praagman).
///
/// \param coor          i/o  Interleaved flat coord array: coor[2*(k-1)]=x, coor[2*(k-1)+1]=y.
/// \param npoint        i/o  Number of nodes (in: boundary nodes; out: all nodes).
/// \param kbound        i    Flat boundary-edge pairs (1-based node indices), length 2*nbound.
/// \param nbound        i    Number of boundary edges.
/// \param kmeshc        o    Output triangle array: 3 node indices per element (1-based), flat.
/// \param nelem         o    Number of triangles generated.
/// \param boundary      i    Flat 2*numcrvboun array: [curveNr, startAddrInKbndpt] per curve.
/// \param numcrvboun    i    Number of curves in the boundary description.
/// \param npunt         i    Estimated total number of nodes in the final mesh.
/// \param inside        i    1 = there are interior areas (holes/constraints).
/// \param holeinfo      i    Flat 2*(nholes+2) array with hole boundary info.
/// \param nholes        i    Number of holes.
/// \param reposition    i    Whether Laplacian repositioning is allowed.
/// \param kbndpt        i    Node-number sequence for boundary curves.
/// \param nbndpt        i    Length of kbndpt.
/// \param coar          i    Flat 3*ncoar array: [x, y, coarseness] per special point.
/// \param ncoar         i    Number of special (user coarseness) points.
/// \param userpoints    i    User-given sequence numbers for the ncoar special points.
/// \param cocurv        i    Flat 2*totalCurvNodes array of internal-curve node coordinates.
/// \param ncurvs        i    Number of internal curves.
/// \param curves        i    Node counts per internal curve (length ncurvs).
/// \param crvnrs        i    User-given curve numbers (length ncurvs).
/// \param istep         i    1 = linear elements, 2 = quadratic.
/// \param extquanodes   i    Help array for quadratic internal line data.
/// \param isurnr        i    Surface sequence number (for diagnostics).
/// \param rinput        i    Real input array: [1..ndim*nuspnt+2+nuspnt] coarseness/user-data.
/// \param nuspnt        i    Number of user-prescribed nodes.
/// \param ndim          i    Spatial dimension of the original problem.
/// \param ctx           i/o  SEPRAN error/constant context.
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
            std::span<const int>    kbndpt,
            int                     nbndpt,
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
            SepranContext&          ctx);

} // namespace sepran
