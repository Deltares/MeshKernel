//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2023.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation version 3.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// contact: delft3d.support@deltares.nl
// Stichting Deltares
// P.O. Box 177
// 2600 MH Delft, The Netherlands
//
// All indications and logos of, and references to, "Delft3D" and "Deltares"
// are registered trademarks of Stichting Deltares, and remain the property of
// Stichting Deltares. All rights reserved.
//
//------------------------------------------------------------------------------

#pragma once

#include <cmath>
#include <concepts>
#include <limits>

namespace meshkernel
{

    /// @brief Generic function for determining if two floating point values are equal
    /// @param[value] The value to compare
    /// @param[ref_value] The reference value to compare to
    /// @param[relative_tol] Relative tolerance to which the values are compared.
    /// @return Boolean indicating whether the value and reference value are equal to a relative tolerance.
    template <std::floating_point T>
    static bool IsEqual(const T value, T ref_value, T relative_tol = 10.0 * std::numeric_limits<T>::epsilon())
    {
        if (value == ref_value)
        {
            return true;
        }

        const T abs_diff = std::abs(value - ref_value);
        const T abs_value = std::abs(value);
        const T abs_ref_value = std::abs(ref_value);

        return abs_diff < relative_tol * std::min(abs_value, abs_ref_value);
    }

} // namespace meshkernel
