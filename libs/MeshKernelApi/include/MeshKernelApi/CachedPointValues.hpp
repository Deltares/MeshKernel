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

#include "MeshKernel/Point.hpp"

#include "MeshKernelApi/GeometryList.hpp"

namespace meshkernelapi
{
    /// @brief Caches x- and y-coordinate values for various algorithms
    class CachedPointValues
    {
    public:
        /// @brief Default constructor
        CachedPointValues() {}

        /// @brief Construct with point values
        CachedPointValues(const std::vector<meshkernel::Point>& coordinates);

        /// @brief Destructor
        virtual ~CachedPointValues() = default;

        /// @brief Number of points saved
        int Size() const
        {
            return static_cast<int>(m_coordsX.size());
        }

        /// @brief Copy cached points to geometry
        void Copy(const GeometryList& geometry) const;

    protected:
        /// @brief Reset the saved coordinate values.
        void Reset(std::vector<double>&& xValues,
                   std::vector<double>&& yValues);

    private:
        std::vector<double> m_coordsX; ///< x-coordinate values
        std::vector<double> m_coordsY; ///< y-coordinate values
    };
} // namespace meshkernelapi
