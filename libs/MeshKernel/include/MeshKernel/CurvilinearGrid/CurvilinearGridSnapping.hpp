//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2023.
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
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLine.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridNodeIndices.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/LandBoundary.hpp>
#include <MeshKernel/Splines.hpp>

namespace meshkernel
{

    /// @brief Computes smoothing factor for a point in the grid.
    class MeshSmoothingCalculator
    {
    public:
        /// @brief Default destructor
        virtual ~MeshSmoothingCalculator() = default;

        /// @brief Compute the mesh smoothing factor.
        ///
        /// @param [in] snappedNodeIndex Index of the snapped grid point
        /// @param [in] gridLinePointIndex Current point on the smoothing grid line.
        virtual double compute(const CurvilinearGridNodeIndices& snappedNodeIndex,
                               const CurvilinearGridNodeIndices& gridLinePointIndex) const = 0;
    };

    /// @brief Computes the directional smoothing factor for a point in the grid.
    ///
    /// The size of the smoothing region is determine by the user selected points.
    class DirectionalSmoothingCalculator : public MeshSmoothingCalculator
    {
    public:
        /// @brief Constructor
        /// @param [in] grid The starting (before smoothing) grid
        /// @param [in] lowerLeft       Index of lower left point of smoothing region
        /// @param [in] upperRight      Index of upper right point of smoothing region
        /// @param [in] regionIndicator
        DirectionalSmoothingCalculator(const CurvilinearGrid& grid,
                                       const CurvilinearGridNodeIndices& lowerLeft,
                                       const CurvilinearGridNodeIndices& upperRight,
                                       const CurvilinearGridNodeIndices& regionIndicator);

        /// @brief Compute the directional smoothing factor.
        /// @param [in] snappedNodeIndex Index of the snapped grid point
        double compute(const CurvilinearGridNodeIndices& snappedNodeIndex,
                       const CurvilinearGridNodeIndices& gridLinePointIndex) const override;

    private:
        /// @brief The original grid before smoothing
        const CurvilinearGrid& m_grid;

        /// @brief Index of lower left point of smoothing region
        CurvilinearGridNodeIndices m_indexBoxLowerLeft;

        /// @brief Index of upper right point of smoothing region
        CurvilinearGridNodeIndices m_indexBoxUpperRight;

        /// @brief Indicator for the smoothing direction.
        CurvilinearGridNodeIndices m_smoothingRegionIndicator;
    };

    /// @brief Computes the non-directional smoothing factor for a point in the grid.
    ///
    /// The size of the smoothing region is pre-determined.
    class NonDirectionalSmoothingCalculator : public MeshSmoothingCalculator
    {
    public:
        /// @brief Constructor
        /// @param [in] The starting (before smoothing) grid
        /// @param [in] The grid to which smoothing is to be applied
        /// @param [in] The landboundary
        NonDirectionalSmoothingCalculator(const CurvilinearGrid& originalGrid,
                                          const CurvilinearGrid& snappedGrid,
                                          const LandBoundary& landBoundary);

        /// @brief Compute the non-direcitonal smoothing factor.
        /// @param [in] snappedNodeIndex Index of the snapped grid point
        double compute(const CurvilinearGridNodeIndices& snappedNodeIndex,
                       const CurvilinearGridNodeIndices& gridLinePointIndex) const override;

    private:
        /// @brief Viewing windows aspect ratio, value is used for the snapping.
        ///
        /// Value from editgridlineblok.f90
        static constexpr double aspectRatio = 990.0 / 1600.0;

        /// @brief How much to enlarge the size of the smoothing region bounding box dimensions.
        static constexpr double smoothingRegionEnlargementFactor = 1.2;

        /// @brief Compute the minimum smoothing region radius.
        ///
        /// Used to set the m_smoothingRegionMinimum member.
        static double CalculateSmoothingRegion(const CurvilinearGrid& grid,
                                               const LandBoundary& landBoundary);

        /// @brief The original grid before smoothing
        const CurvilinearGrid& m_originalGrid;

        /// @brief The grid to which the smoothing is to be applied.
        const CurvilinearGrid& m_snappedGrid;

        /// @brief The minimum smoothing region radius
        double m_smoothingRegionMinimum{0.0};
    };

    /// @brief Smoothly snap the grid to a land boundary or spline.
    class CurvilinearGridSnapping : public CurvilinearGridAlgorithm
    {
    public:
        /// @brief constructor
        /// @param [in] grid         The input curvilinear grid
        /// @param [in] landBoundary The land boundary to which the grid is to be snapped.
        /// @param [in] points       The points used to control the snapping and smoothing.
        CurvilinearGridSnapping(std::shared_ptr<CurvilinearGrid> grid,
                                const LandBoundary& landBoundary,
                                const std::vector<Point>& points);

        /// @brief Executes the snapping and smoothing algorithm
        CurvilinearGrid Compute() override;

    private:
        /// @brief Tolerance to determine if point is on (close to) boundary
        static constexpr double epsilon = 1.0e-5;

        /// @brief Smoothing region scaling.
        ///
        /// Used when there are exactly 2 domain boundary points selected, to include a predefined smoothing region.
        /// Value for nump from modgr4.f90
        static constexpr UInt predefinedSmootingRegionFactor = 80;

        /// @brief Smoothing region scaling.
        ///
        /// Used when there are more than 2 domain boundary points selected, to include a user defined smoothing region.
        /// Value from modgr4.f90
        static constexpr UInt userDefinedSmootingRegionFactor = 10000;

        /// @brief Compute the loop bounds for the smoothing.
        /// @param [in] snappedNodeIndex Index of the grid point snapped to the land boundary or spline
        std::tuple<CurvilinearGridNodeIndices, CurvilinearGridNodeIndices> ComputeLoopBounds(const CurvilinearGridNodeIndices& snappedNodeIndex) const;

        /// @brief Apply the smoothing to the grid
        /// @param [in] snappedNodeIndex Index of the grid point snapped to the land boundary or spline
        /// @param [in] smoothingFactor   Calculator to compute the amount of smoothing to be applied at a particular point in the grid
        void ApplySmoothingToGrid(const CurvilinearGridNodeIndices& snappedNodeIndex,
                                  const MeshSmoothingCalculator& smoothingFactor);

        /// @brief Initialise member variables
        void Initialise();

        /// @brief The grid to be smoothed
        const CurvilinearGrid& m_originalGrid;

        /// @brief The land boundary to which the grid is to be snapped.
        const LandBoundary& m_landBoundary;

        /// @brief The control points for the snapping.
        const std::vector<Point> m_points;

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
