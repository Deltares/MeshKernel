//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2020.
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

#include <vector>
#include <MeshKernel/Entities.hpp>

namespace meshkernel
{
    class TriangulationInterpolation
    {

    public:
        /// @brief Constructor
        /// @param[in] locations interpolation points (where the values should be computed)
        /// @param[in] samples  Values to use for the interpolation
        /// @param[in] projection Projection to use (\ref Projections)
        TriangulationInterpolation(const std::vector<Point>& locations,
                                   const std::vector<Sample>& samples,
                                   Projections projection);

        /// @brief Compute results on the interpolation points
        void Compute();

        /// @brief Get the results
        /// @return
        [[nodiscard]] const auto& GetResults() const
        {
            return m_results;
        }

    private:
        const std::vector<Point>& m_locations;
        const std::vector<Sample>& m_samples;
        Projections m_projection;
        std::vector<double> m_results;
    };

}; // namespace meshkernel
