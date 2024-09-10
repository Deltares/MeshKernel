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
#include <MeshKernel/Exceptions.hpp>

namespace meshkernel
{

    /// @brief A struct describing a node in the curvilinear grid in terms of node indices
    struct CurvilinearGridNodeIndices
    {
        /// @brief Default constructor sets the indices to invalid
        CurvilinearGridNodeIndices() : m_n(constants::missing::uintValue), m_m(constants::missing::uintValue) {}

        /// @brief Constructor setting n and m indices
        /// @param[in] m The m index
        /// @param[in] n The n index
        CurvilinearGridNodeIndices(UInt n, UInt m) : m_n(n), m_m(m) {}

        /// @brief Determines if one of the indices  equals to \p missingValue
        [[nodiscard]] bool IsValid(const UInt missingValue = constants::missing::uintValue) const
        {
            return m_m != missingValue && m_n != missingValue;
        }

        CurvilinearGridNodeIndices& operator+=(const CurvilinearGridNodeIndices& val)
        {
            if (!IsValid())
            {
                throw ConstraintError("Invalid node index");
            }

            if (!val.IsValid())
            {
                throw ConstraintError("Invalid node index increment");
            }

            m_n += val.m_n;
            m_m += val.m_m;
            return *this;
        }

        /// @brief Overloads <=> operator
        auto operator<=>(const CurvilinearGridNodeIndices& rhs) const
        {
            if (m_n == rhs.m_n)
            {
                return m_m <=> rhs.m_m; // Return result of m_m comparison if not equal
            }
            return m_n <=> rhs.m_n; // Compare m_n if m_m values are equal
        }

        /// @brief Overloads equality with another CurvilinearGridNodeIndices
        bool operator==(const CurvilinearGridNodeIndices& rhs) const = default;

        CurvilinearGridNodeIndices& operator-=(const CurvilinearGridNodeIndices& val)
        {
            if (!IsValid())
            {
                throw ConstraintError("Invalid node index");
            }

            if (!val.IsValid())
            {
                throw ConstraintError("Invalid node index increment");
            }

            m_n -= val.m_n;
            m_m -= val.m_m;
            return *this;
        }

        /// @brief Inquires if another node is on the same grid line of the current node
        /// @param[in] rhs The node to inquire
        /// @return True if on the same grid line, false otherwise
        [[nodiscard]] bool IsOnTheSameGridLine(const CurvilinearGridNodeIndices& rhs) const
        {
            return m_n == rhs.m_n || m_m == rhs.m_m;
        }

        UInt m_n; ///< Rows
        UInt m_m; ///< Columns
    };
} // namespace meshkernel
