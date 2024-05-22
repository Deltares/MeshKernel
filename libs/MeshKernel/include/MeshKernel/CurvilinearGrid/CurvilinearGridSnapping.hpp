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

#include <memory>
#include <tuple>
#include <vector>

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridNodeIndices.hpp>
#include <MeshKernel/CurvilinearGrid/MeshSmoothingCalculator.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/UndoActions/UndoAction.hpp>

namespace meshkernel
{

    /// @brief Smoothly snap the grid to a land boundary or spline.
    class CurvilinearGridSnapping
    {
    public:
        /// @brief constructor
        /// @param [in] grid         The input curvilinear grid
        /// @param [in] landBoundary The land boundary to which the grid is to be snapped.
        /// @param [in] points       The points used to control the snapping and smoothing.
        CurvilinearGridSnapping(CurvilinearGrid& grid,
                                const std::vector<Point>& points);

        /// @brief Default destructor
        virtual ~CurvilinearGridSnapping() = default;

        /// @brief Executes the snapping and smoothing algorithm
        [[nodiscard]] UndoActionPtr Compute();

    private:
        /// @brief Tolerance to determine if point is on (close to) boundary
        static constexpr double epsilon = 1.0e-5;

        /// @brief Smoothing region scaling.
        ///
        /// Used when there are exactly 2 domain boundary points selected, to include a predefined smoothing region.
        /// Value for nump from modgr4.f90
        static constexpr UInt predefinedSmoothingRegionFactor = 80;

        /// @brief Smoothing region scaling.
        ///
        /// Used when there are more than 2 domain boundary points selected, to include a user defined smoothing region.
        /// Value from modgr4.f90
        static constexpr UInt userDefinedSmoothingRegionFactor = 10000;

        /// @brief Allocate the grid smoothing calculator
        std::unique_ptr<MeshSmoothingCalculator> AllocateGridSmoothingCalculator() const;

        /// @brief Allocate the grid smoothing calculator
        ///
        /// Abstract factory like method allocating the required smoothing calculator.
        virtual std::unique_ptr<MeshSmoothingCalculator> AllocateGridSmoothingCalculator(const CurvilinearGrid& originalGrid,
                                                                                         const CurvilinearGrid& snappedGrid) const = 0;

        /// @brief Find the closest point to the current point given.
        virtual Point FindNearestPoint(const Point& currentPoint) const = 0;

        /// @brief Compute the loop bounds for the smoothing.
        /// @param [in] snappedNodeIndex Index of the grid point snapped to the land boundary or spline
        std::tuple<CurvilinearGridNodeIndices, CurvilinearGridNodeIndices> ComputeLoopBounds(const CurvilinearGridNodeIndices& snappedNodeIndex) const;

        /// @brief Apply the smoothing to the grid
        /// @param [in] snappedNodeIndex Index of the grid point snapped to the land boundary or spline
        /// @param [in] smoothingFactor  Calculate how much smoothing translation is required for this index
        void ApplySmoothingToGrid(const CurvilinearGridNodeIndices& snappedNodeIndex,
                                  const MeshSmoothingCalculator& smoothingFactor);

        /// @brief Initialise member variables
        void Initialise();

        /// @brief The grid to be smoothed
        CurvilinearGrid& m_grid;

        /// @brief Copy of the grid to be smoothed
        const CurvilinearGrid m_originalGrid;

        /// @brief The control points for the snapping
        const std::vector<Point> m_controlPoints;

        /// @brief The start point index of the grid line to be snapped
        CurvilinearGridNodeIndices m_lineStartIndex;

        /// @brief The end point index of the grid line to be snapped
        CurvilinearGridNodeIndices m_lineEndIndex;

        /// @brief The start points index of the grid line to be snapped
        CurvilinearGridNodeIndices m_smoothingRegionIndicator;

        /// @brief Index of lower left point for smoothing region
        CurvilinearGridNodeIndices m_indexBoxLowerLeft;

        /// @brief Index of upper right point for smoothing region
        CurvilinearGridNodeIndices m_indexBoxUpperRight;
    };

} // namespace meshkernel
