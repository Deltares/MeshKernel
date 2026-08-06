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
