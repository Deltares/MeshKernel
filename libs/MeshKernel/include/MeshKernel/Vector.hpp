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

namespace meshkernel
{
    /// @brief A class defining a vector
    class Vector
    {
    public:
        /// @brief Class constructor
        ///
        /// @param[in] x The x coordinate of the vector
        /// @param[in] y The y coordinate of the vector
        Vector(const double x, const double y)
            : m_x(x),
              m_y(y)
        {
        }

        /// @brief Gets the x coordinate of the vector
        /// @returns x The x coordinate of the vector
        [[nodiscard]] double x() const
        {
            return m_x;
        }

        /// @brief Gets the y coordinate of the vector
        /// @returns x The y coordinate of the vector
        [[nodiscard]] double y() const
        {
            return m_y;
        }

    private:
        double m_x; ///< The x coordinate of the vector
        double m_y  ///< The y coordinate of the vector
    };
} // namespace meshkernel
