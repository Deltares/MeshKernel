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
