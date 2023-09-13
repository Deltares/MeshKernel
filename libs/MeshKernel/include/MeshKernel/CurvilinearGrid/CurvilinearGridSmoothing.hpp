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

#include <memory>

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridAlgorithm.hpp>
#include <MeshKernel/Entities.hpp>

namespace meshkernel
{

    /// @brief A class implementing the curvilinear grid smoothing algorithm
    class CurvilinearGridSmoothing : public CurvilinearGridAlgorithm
    {
    public:
        /// @brief Class constructor
        /// @param[in] grid                        The input curvilinear grid
        /// @param[in] smoothingIterations         The number of smoothing iterations to perform
        CurvilinearGridSmoothing(std::shared_ptr<CurvilinearGrid> grid, UInt smoothingIterations);

        /// @brief Compute curvilinear grid block smoothing (modifies the m_grid nodal values)
        /// @return The smoothed grid
        CurvilinearGrid Compute() override;

        /// @brief Compute curvilinear grid line smoothing. The algorithm smooths the grid along the direction specified by the line.
        /// The line must be an m or n grid line of the curvilinear grid.
        /// @return The smoothed grid
        CurvilinearGrid ComputeDirectional();

    private:
        /// @brief Solve one iteration of block smoothing
        void Solve();

        /// @brief Solve one iteration of directional smoothing
        void SolveDirectional();

        /// @brief Projects a point on the closest grid boundary
        /// @param[in] point The point to project
        /// @param[in] m The current m coordinate on the boundary of the curvilinear grid
        /// @param[in] n The current n coordinate on the boundary of the curvilinear grid
        void ProjectPointOnClosestGridBoundary(Point const& point, UInt m, UInt n);

        UInt m_smoothingIterations;              ///< The orthogonalization parameters
        lin_alg::Matrix<Point> m_gridNodesCache; ///< A cache for storing current iteration node positions
    };
} // namespace meshkernel
