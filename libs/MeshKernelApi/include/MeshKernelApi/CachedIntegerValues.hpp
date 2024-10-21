//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2024.
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

#include <cstring>
#include <vector>

namespace meshkernelapi
{
    /// @brief Caches x- and y-coordinate values for various algorithms
    class CachedIntegerValues
    {
    public:
        /// @brief Default constructor
        CachedIntegerValues() {}

        /// @brief Construct with point values
        CachedIntegerValues(const std::vector<int>& values);

        /// @brief Destructor
        virtual ~CachedIntegerValues() = default;

        /// @brief Number of points saved
        int Size() const
        {
            return static_cast<int>(m_values.size());
        }

        /// @brief Copy cached points to geometry
        void Copy(int* buffer) const;

    protected:
        /// @brief Reset the saved integer values.
        void Reset(std::vector<int>&& values);

    private:
        std::vector<int> m_values; ///< integer values
    };
} // namespace meshkernelapi
