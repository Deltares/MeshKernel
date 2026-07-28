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
// C++20 translation of mshtrans2dsur.for from the SEPRAN library.

#pragma once

#include "SepranContext.hpp"

#include <span>

namespace sepran
{

/// \brief Parameters computed by the 2D coordinate normalisation step.
struct TransformParams
{
    double xmint; ///< X translation (minimum x before scaling).
    double ymint; ///< Y translation (minimum y before scaling).
    double tran;  ///< Scale factor (maximum extent in either dimension).
};

/// \brief Normalise coordinates so the domain fits in the first quadrant with
///        unit extent.
///
/// Translated from Fortran \c mshtrans2dsur (SEPRAN, Guus Segal, 2009).
///
/// All coordinate arrays are modified in-place.  After the call, every
/// coordinate c becomes 1 + (c - cmin) / tran.
///
/// Array conventions (shared with the rest of sepran):
/// - \p coor: interleaved flat [x0,y0, x1,y1, ...], 0-based, length 2*npoint.
/// - \p coar: flat [x0,y0,coarseness0, x1,y1,coarseness1, ...], length 3*ncoar.
/// - \p cocurvs: interleaved flat, length 2 * total_curve_nodes.
/// - \p userco: interleaved flat, length 2*nuspnt.
/// - \p curves: node counts per internal curve, length ncurvs.
///
/// \param coor     Node coordinates to normalise (in-out).
/// \param npoint   Number of nodes in coor.
/// \param coar     Internal-point coordinate + coarseness array (in-out).
/// \param ncoar    Number of entries in coar.
/// \param curves   Node count per internal curve (length ncurvs).
/// \param ncurvs   Number of internal curves.
/// \param cocurvs  Coordinates of nodes on internal curves (in-out).
/// \param userco   User-prescribed node coordinates (in-out), length 2*nuspnt.
/// \param nuspnt   Number of user points.
/// \param ctx      SEPRAN context.
/// \return TransformParams with xmint, ymint and tran.
TransformParams transform2DSurface(std::span<double> coor,
                                   int npoint,
                                   std::span<double> coar,
                                   int ncoar,
                                   std::span<const int> curves,
                                   int ncurvs,
                                   std::span<double> cocurvs,
                                   std::span<double> userco,
                                   int nuspnt,
                                   SepranContext& ctx);

/// \brief Reverse the transformation applied by transform2DSurface.
///
/// Applies the inverse: c_orig = xmint + (c_scaled - 1) * tran.
///
/// \param coor    Node coordinate array to reverse-transform (in-out).
/// \param npoint  Number of nodes.
/// \param tp      Parameters returned by transform2DSurface.
void reverseTransform2DSurface(std::span<double> coor,
                               int npoint,
                               const TransformParams& tp);

} // namespace sepran
