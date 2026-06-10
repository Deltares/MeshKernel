// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Stichting Deltares
//
// C++20 translation of the Fortran module mshconstants.f90 and msherror module
// from the SEPRAN library ("Ingenieursbureau SEPRA").
//
// This header is intentionally header-only (no .cpp needed).

#pragma once

namespace sepran
{

/// \brief Runtime context carrying numerical constants and error state.
///
/// In the original Fortran this was split across two Fortran modules:
///   - mshconstants: EPSMAC, SQREPS, RINFIN, IREFWR, ITIME, IGOBS, JTIMES
///   - msherror:     IERROR
///
/// In C++ we combine both into a single context object that is passed (by
/// reference) through the call stack instead of being global module state.
struct SepranContext
{
    double epsmac  = 1.0e-15; ///< Machine epsilon (EPSMAC in Fortran)
    double sqreps  = 1.0e-15; ///< Square root of machine epsilon (SQREPS)
    double rinfin  = 1.0e77;  ///< Large sentinel value (RINFIN)

    int itime      = 1;       ///< Timing flag (ITIME)
    int igobs      = 1;       ///< Observer flag (IGOBS)
    int jtimes     = 1;       ///< Timing multiplier (JTIMES)

    int ierror     = 0;       ///< Error flag; non-zero indicates an error (IERROR)

    /// \brief Return a default-initialized context with the standard SEPRAN constants.
    static SepranContext defaults() noexcept
    {
        return SepranContext{};
    }
};

} // namespace sepran
