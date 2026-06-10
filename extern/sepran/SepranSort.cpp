// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Stichting Deltares
//
// C++20 translation of chsort.for from the SEPRAN library.
//
// Original Fortran: chsort (programmer Jos van Kan / Onno Hoitinga, 1987-1998)
// Numerical-Recipes heap sort algorithm (Cambridge University Press, 1989).

#include "SepranSort.hpp"

#include <algorithm>
#include <cassert>
#include <numeric>

namespace sepran
{

void heapSort(std::span<const int> keysrt, std::span<int> kgrad)
{
    const int n = static_cast<int>(keysrt.size());
    assert(kgrad.size() == keysrt.size());

    // Initialise permutation to identity (Fortran: do j=1,npoint kgrad(j)=j)
    std::iota(kgrad.begin(), kgrad.end(), 0);

    if (n <= 1)
        return;

    // Heap sort — direct translation of the Fortran, converted to 0-based indices.
    // The Fortran algorithm works on a 1-based kgrad array; here every index is
    // decremented by 1.
    int l  = n / 2;   // Fortran: l = npoint/2+1, then l=l-1 on first entry → same
    int ir = n - 1;   // Fortran: ir = npoint  (0-based: n-1)

    while (true)
    {
        int kgradt;
        int q;

        if (l > 0)
        {
            --l;
            kgradt = kgrad[l];
            q      = keysrt[kgradt];
        }
        else
        {
            kgradt    = kgrad[ir];
            q         = keysrt[kgradt];
            kgrad[ir] = kgrad[0];
            --ir;
            if (ir == 0)
            {
                kgrad[0] = kgradt;
                return;
            }
        }

        int i = l;
        int j = l + l + 1; // Fortran j = l+l, but kgrad is 0-based → j = 2*l+1

        while (j <= ir)
        {
            if (j < ir)
            {
                if (keysrt[kgrad[j]] < keysrt[kgrad[j + 1]])
                    ++j;
            }
            if (q < keysrt[kgrad[j]])
            {
                kgrad[i] = kgrad[j];
                i        = j;
                j        = j + j + 1;
            }
            else
            {
                break;
            }
        }
        kgrad[i] = kgradt;
    }
}

} // namespace sepran
