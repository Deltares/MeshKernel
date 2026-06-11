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
// C++20 translation of mshcrossline.for, mshcrossline1.for, and msho75.for
// from the SEPRAN library.

#pragma once

#include "SepranContext.hpp"

#include <array>

namespace sepran
{

/// \brief Find the parametric intersection of the infinite lines x1–x2 and x3–x4.
///
/// Translated from Fortran \c mshcrossline (SEPRAN, Guus Segal, 2000-2008).
///
/// The intersection point on the first line is defined by x1 + fact1*(x2-x1).
/// The intersection point on the second line is defined by x3 + fact2*(x4-x3).
/// If the lines are (nearly) parallel, both facts are set to -ctx.rinfin.
///
/// \param x1, x2  Endpoints of first line segment (2D coordinates).
/// \param x3, x4  Endpoints of second line segment (2D coordinates).
/// \param fact1   Output: parametric factor for x1–x2.
/// \param fact2   Output: parametric factor for x3–x4.
/// \param eps     Local accuracy tolerance.
/// \param ctx     SEPRAN context (reads ctx.epsmac, ctx.rinfin).
void crossLine(const std::array<double, 2>& x1,
               const std::array<double, 2>& x2,
               const std::array<double, 2>& x3,
               const std::array<double, 2>& x4,
               double& fact1,
               double& fact2,
               double eps,
               const SepranContext& ctx);

/// \brief Like crossLine, but first tests whether the bounding rectangles overlap.
///
/// Translated from Fortran \c mshcrossline1 (SEPRAN, Guus Segal, 2003).
///
/// If the bounding boxes of x1–x2 and x3–x4 do not overlap, both facts are
/// returned as -ctx.rinfin without calling the full intersection computation.
///
/// \param x1, x2  Endpoints of first line segment.
/// \param x3, x4  Endpoints of second line segment.
/// \param fact1   Output: parametric factor for x1–x2.
/// \param fact2   Output: parametric factor for x3–x4.
/// \param eps     Local accuracy tolerance.
/// \param ctx     SEPRAN context.
void crossLine1(const std::array<double, 2>& x1,
                const std::array<double, 2>& x2,
                const std::array<double, 2>& x3,
                const std::array<double, 2>& x4,
                double& fact1,
                double& fact2,
                double eps,
                const SepranContext& ctx);

/// \brief Check whether line segment i–j and line segment 1–2 share a common point.
///
/// Translated from Fortran \c msho75 (SEPRAN, Niek Praagman, 1996-2008).
///
/// \return true if the two segments intersect (ih==0 in Fortran), false otherwise.
bool segmentsIntersect(double xi, double yi,
                       double xj, double yj,
                       double x1, double y1,
                       double x2, double y2,
                       const SepranContext& ctx);

} // namespace sepran
