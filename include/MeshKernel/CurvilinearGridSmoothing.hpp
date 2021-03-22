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

#include <MeshKernel/Entities.hpp>

namespace meshkernel
{
    class CurvilinearGrid;

    /// @brief A class implementing the curvilinear grid smoothing algorithm
    class CurvilinearGridSmoothing
    {
    public:
        /// @brief Class constructor
        /// @param[in] grid                        The input curvilinear grid
        /// @param[in] smoothingIterations         The number of smoothing iterations to perform
        /// @param[in] firstCornerPoint            The first point defining the smoothing bounding box
        /// @param[in] secondCornerPoint           The second point defining the smoothing bounding box
        CurvilinearGridSmoothing(std::shared_ptr<CurvilinearGrid> grid,
                                 size_t smoothingIterations,
                                 const Point& firstCornerPoint,
                                 const Point& secondCornerPoint);

        /// @brief Compute curvilinear grid smoothing (modifies the m_grid nodal values)
        void Compute();

    private:
        /// @brief Solve one smoothing iteration
        void Solve();

        std::shared_ptr<CurvilinearGrid> m_grid; ///< A pointer to the curvilinear grid to modify
        size_t m_smoothingIterations;            ///< The orthogonalization parameters
        Point m_firstCornerPoint;                ///< The first point defining the orthogonalization bounding box
        Point m_secondCornerPoint;               ///< The second point defining the orthogonalization bounding box

        size_t m_minM; ///< The minimum m grid index of the orthogonalization bounding box
        size_t m_minN; ///< The minimum n grid index of the orthogonalization bounding box
        size_t m_maxM; ///< The maximum m grid index of the orthogonalization bounding box
        size_t m_maxN; ///< The maximum n grid index of the orthogonalization bounding box

        std::vector<std::vector<Point>> m_gridNodesCache; ///< A cache for storing current iteration node positions
    };
} // namespace meshkernel