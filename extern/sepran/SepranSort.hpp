// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Stichting Deltares
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
