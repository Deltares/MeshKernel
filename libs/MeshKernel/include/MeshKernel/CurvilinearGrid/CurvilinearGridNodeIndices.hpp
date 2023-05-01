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
#include <stddef.h>

#include <MeshKernel/Constants.hpp>

namespace meshkernel
{

    /// @brief A struct describing a node in the curvilinear grid in terms of node indices
    struct CurvilinearGridNodeIndices
    {
        /// @brief Default constructor sets the indices to invalid
        CurvilinearGridNodeIndices() : m_m(constants::missing::sizetValue), m_n(constants::missing::sizetValue){};

        /// @brief Constructor sets indices from values
        /// @param[in] m The m index
        /// @param[in] n The n index
        CurvilinearGridNodeIndices(size_t m, size_t n) : m_m(m), m_n(n){};

        /// @brief Determines if one of the indices  equals to \p missingValue
        [[nodiscard]] bool IsValid(const size_t missingValue = constants::missing::sizetValue) const
        {
            return m_m != missingValue && m_n != missingValue;
        }

        /// @brief Overloads equality with another CurvilinearGridNodeIndices
        bool operator==(const CurvilinearGridNodeIndices& rhs) const = default;

        /// @brief Inquires if another node is on the same grid line of the current node
        /// @param[in] rhs The node to inquire
        /// @return True if on the same grid line, false otherwise
        [[nodiscard]] bool IsOnTheSameGridLine(const CurvilinearGridNodeIndices& rhs) const
        {
            return m_m == rhs.m_m || m_n == rhs.m_n;
        }

        size_t m_m; ///< Columns
        size_t m_n; ///< Rows
    };
} // namespace meshkernel
