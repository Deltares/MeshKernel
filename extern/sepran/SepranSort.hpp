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
// C++20 translation of chsort.for from the SEPRAN library.

#pragma once

#include <span>

namespace sepran
{

/// \brief Compute a permutation index that sorts keysrt in ascending order.
///
/// Translated from Fortran \c chsort (SEPRAN library, Numerical-Recipes heap sort).
///
/// \param keysrt  Integer sort keys (0-based, length n).
/// \param kgrad   Output permutation: kgrad[i] is the 0-based index of the
///                i-th element in sorted order.  Must be the same size as keysrt.
void heapSort(std::span<const int> keysrt, std::span<int> kgrad);

} // namespace sepran
