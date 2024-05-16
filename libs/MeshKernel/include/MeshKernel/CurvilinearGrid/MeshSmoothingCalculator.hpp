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

#include "MeshKernel/BoundingBox.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridNodeIndices.hpp"
#include "MeshKernel/LandBoundary.hpp"
#include "MeshKernel/Splines.hpp"

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
        /// @param [in] lowerLeft       Index of lower left point of smoothing region
        /// @param [in] upperRight      Index of upper right point of smoothing region
        /// @param [in] regionIndicator
        DirectionalSmoothingCalculator(const CurvilinearGridNodeIndices& lowerLeft,
                                       const CurvilinearGridNodeIndices& upperRight,
                                       const CurvilinearGridNodeIndices& regionIndicator);

        /// @brief Compute the directional smoothing factor.
        /// @param [in] snappedNodeIndex Index of the snapped grid point
        double compute(const CurvilinearGridNodeIndices& snappedNodeIndex,
                       const CurvilinearGridNodeIndices& gridLinePointIndex) const override;

    private:
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

        NonDirectionalSmoothingCalculator(const CurvilinearGrid& originalGrid,
                                          const CurvilinearGrid& snappedGrid,
                                          const Splines& spline);

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
        static double CalculateSmoothingRegion(const BoundingBox& gridBoundingBox,
                                               const BoundingBox& entityBoundingBox);

        /// @brief The original grid before smoothing
        const CurvilinearGrid& m_originalGrid;

        /// @brief The grid to which the smoothing is to be applied.
        const CurvilinearGrid& m_snappedGrid;

        /// @brief The minimum smoothing region radius
        double m_smoothingRegionMinimum{0.0};
    };

} // namespace meshkernel
