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

#include <MeshKernel/CurvilinearGrid.hpp>
#include <MeshKernel/Entities.hpp>

namespace meshkernel
{

    /// @brief A class implementing the curvilinear grid smoothing algorithm
    class CurvilinearGridSmoothing
    {
    public:
        /// @brief Class constructor
        /// @param[in] grid                        The input curvilinear grid
        /// @param[in] smoothingIterations         The number of smoothing iterations to perform
        CurvilinearGridSmoothing(std::shared_ptr<CurvilinearGrid> grid,
                                 size_t smoothingIterations);

        /// @brief Compute curvilinear grid block smoothing (modifies the m_grid nodal values)
        void Compute();

        /// @brief Compute curvilinear grid line smoothing. The algorithm smooths the grid along the direction specified by the line.
        /// The line must be an m or n grid line of the curvilinear grid. The grid is smoothed in the region (influence zone) specified by two corner points.
        /// @param[in] firstLinePoint The first point of the line
        /// @param[in] secondLinePoint The first point of the line
        /// @param[in] leftPointRegion The left point of the smoothing region
        /// @param[in] rightPointRegion The right point of the smoothing region
        void ComputeLineSmoothOnRegion(Point const& firstLinePoint, Point const& secondLinePoint, Point const& leftPointRegion, Point const& rightPointRegion);

        /// @brief Sets the orthogonalization block (TODO: Create base class for curvi orthogonalization and smoothing)
        /// @param[in] firstCornerPoint            The first point defining the orthogonalization bounding box
        /// @param[in] secondCornerPoint           The second point defining the orthogonalization bounding box
        void SetBlock(Point const& firstCornerPoint, Point const& secondCornerPoint);

    private:
        /// @brief Solve one smoothing iteration
        void Solve();

        /// @brief Projects a point on the closest grid boundary
        /// @param point The point to project
        /// @param m The current m coordinate on the boundary of the curvilinear grid
        /// @param n The current n coordinate on the boundary of the curvilinear grid
        void ProjectPointOnClosestGridBoundary(Point const& point, size_t m, size_t n);

        /// @brief Function for computing the smoothing factor at the current location given a line and an area of influence (SMEERFUNCTIE)
        /// @param currentPoint
        /// @param firstLinePointIndices
        /// @param secondLinePointIndices
        /// @param leftPointInfluenceZone
        /// @param rightPointInfluenceZone
        /// @return
        std::tuple<double, double, double> ComputeSmoothingFactors(CurvilinearGrid::NodeIndices const& currentPoint,
                                                                   const CurvilinearGrid::NodeIndices& firstLinePointIndices,
                                                                   const CurvilinearGrid::NodeIndices& secondLinePointIndices,
                                                                   const CurvilinearGrid::NodeIndices& leftPointInfluenceZone,
                                                                   const CurvilinearGrid::NodeIndices& rightPointInfluenceZone) const;

        std::shared_ptr<CurvilinearGrid> m_grid; ///< A pointer to the curvilinear grid to modify
        size_t m_smoothingIterations;            ///< The orthogonalization parameters

        size_t m_minM; ///< The minimum m grid index of the smoothing block
        size_t m_minN; ///< The minimum n grid index of the smoothing block
        size_t m_maxM; ///< The maximum m grid index of the smoothing block
        size_t m_maxN; ///< The maximum n grid index of the smoothing block

        std::vector<std::vector<Point>> m_gridNodesCache; ///< A cache for storing current iteration node positions
    };
} // namespace meshkernel