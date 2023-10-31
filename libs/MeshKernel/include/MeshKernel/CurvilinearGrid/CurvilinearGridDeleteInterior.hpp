//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2021.
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

#include "CurvilinearGridAlgorithm.hpp"

namespace meshkernel
{

    /// @brief A class implementing the curvilinear grid delete interior algorithm.
    class CurvilinearGridDeleteInterior : public CurvilinearGridAlgorithm
    {
    public:
        /// @brief Class constructor
        /// @param[in] grid The input curvilinear grid
        CurvilinearGridDeleteInterior(CurvilinearGrid& grid);

        /// @brief Set all the nodes interior to the block to be invalid.
        /// @note If any index is the null value or is out of range a ConstraintError will be thrown
        void Compute() override;
    };
} // namespace meshkernel
