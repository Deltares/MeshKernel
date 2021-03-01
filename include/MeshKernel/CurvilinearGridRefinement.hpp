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

#include <MeshKernel/Entities.hpp>
#include <memory>

namespace meshkernel
{
    class CurvilinearGrid;

    /// @brief A class implementing the curvilinear grid refinement algorithm
    class CurvilinearGridRefinement
    {
    public:
        /// @brief Class constructor
        ///
        /// \p firstPoint and \p secondPoint must lie on the same gridline
        /// @param[in] grid The input curvilinear grid
        /// @param[in] firstPoint  The first node of the segment defining the refinement zone
        /// @param[in] secondPoint The second node of the segment defining the refinement zone
        /// @param[in] mRefinement The selected number of refinement lines along the m-direction
        /// @param[in] nRefinement The selected number of refinement lines along the n-direction
        CurvilinearGridRefinement(const std::shared_ptr<CurvilinearGrid>& grid,
                                  const Point& firstPoint,
                                  const Point& secondPoint,
                                  size_t mRefinement,
                                  size_t nRefinement);

        /// @brief Refine the curvilinear grid
        void Compute();

    private:
        /// @brief Computes the m and n grid lines and spline derivatives in separate vectors
        /// @return The m grid lines and the spline derivatives, the n grid lines and the splines derivatives
        [[nodiscard]] std::tuple<std::vector<std::vector<Point>>,
                                 std::vector<std::vector<Point>>,
                                 std::vector<std::vector<Point>>,
                                 std::vector<std::vector<Point>>>
        ComputeGridLinesAndSplinesDerivatives() const;

        /// @brief Compute spline derivatives along a gridline, also accounting for missing values
        /// @param[in] gridLine The input gridline
        /// @returns The spline derivatives
        std::vector<Point> ComputeSplineDerivatesAlongGridLine(const std::vector<Point>& gridLine) const;

        std::shared_ptr<CurvilinearGrid> m_grid; ///< A pointer to the curvilinear grid to modify
        Point m_firstPoint;                      ///< The first node of the segment defining the refinement zone
        Point m_secondPoint;                     ///< The second node of the segment defining the refinement zone
        size_t m_refinement;                     ///< The selected number of refinement lines along the m-direction
        size_t n_refinement;                     ///< The selected number of refinement lines along the n-direction
    };
} // namespace meshkernel
