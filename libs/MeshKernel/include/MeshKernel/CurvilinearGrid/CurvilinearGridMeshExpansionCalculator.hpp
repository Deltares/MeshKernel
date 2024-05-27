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

    /// @brief Computes expansion factor for a point in the grid.
    class CurvilinearGridMeshExpansionCalculator
    {
    public:
        /// @brief Default destructor
        virtual ~CurvilinearGridMeshExpansionCalculator() = default;

        /// @brief Compute the mesh expansion factor.
        ///
        /// @param [in] snappedNodeIndex Index of the snapped grid point
        /// @param [in] gridLinePointIndex Current point on the expansion grid line.
        virtual double compute(const CurvilinearGridNodeIndices& snappedNodeIndex,
                               const CurvilinearGridNodeIndices& gridLinePointIndex) const = 0;
    };

    /// @brief Computes the directional expansion factor for a point in the grid.
    ///
    /// The size of the expansion region is determine by the user selected points.
    class UserDefinedRegionExpasionCalculator : public CurvilinearGridMeshExpansionCalculator
    {
    public:
        /// @brief Constructor
        /// @param [in] lowerLeft       Index of lower left point of expansion region
        /// @param [in] upperRight      Index of upper right point of expansion region
        /// @param [in] regionIndicator
        UserDefinedRegionExpasionCalculator(const CurvilinearGridNodeIndices& lowerLeft,
                                            const CurvilinearGridNodeIndices& upperRight,
                                            const CurvilinearGridNodeIndices& regionIndicator);

        /// @brief Compute the directional expansion factor.
        /// @param [in] snappedNodeIndex Index of the snapped grid point
        double compute(const CurvilinearGridNodeIndices& snappedNodeIndex,
                       const CurvilinearGridNodeIndices& gridLinePointIndex) const override;

    private:
        /// @brief Index of lower left point of expansion region
        CurvilinearGridNodeIndices m_indexBoxLowerLeft;

        /// @brief Index of upper right point of expansion region
        CurvilinearGridNodeIndices m_indexBoxUpperRight;

        /// @brief Indicator for the expansion direction.
        CurvilinearGridNodeIndices m_expansionRegionIndicator;
    };

    /// @brief Computes the non-directional expansion factor for a point in the grid.
    ///
    /// The size of the expansion region is pre-determined.
    class DefaultRegionExpasionCalculator : public CurvilinearGridMeshExpansionCalculator
    {
    public:
        /// @brief Constructor
        /// @param [in] The starting (before expansion) grid
        /// @param [in] The grid to which expansion is to be applied
        /// @param [in] The landboundary
        DefaultRegionExpasionCalculator(const CurvilinearGrid& originalGrid,
                                        const CurvilinearGrid& snappedGrid,
                                        const LandBoundary& landBoundary);

        DefaultRegionExpasionCalculator(const CurvilinearGrid& originalGrid,
                                        const CurvilinearGrid& snappedGrid,
                                        const Splines& spline);

        /// @brief Compute the non-direcitonal expansion factor.
        /// @param [in] snappedNodeIndex Index of the snapped grid point
        double compute(const CurvilinearGridNodeIndices& snappedNodeIndex,
                       const CurvilinearGridNodeIndices& gridLinePointIndex) const override;

    private:
        /// @brief Viewing windows aspect ratio, value is used for the snapping.
        ///
        /// Value from editgridlineblok.f90
        static constexpr double aspectRatio = 990.0 / 1600.0;

        /// @brief How much to enlarge the size of the expansion region bounding box dimensions.
        static constexpr double expansionRegionEnlargementFactor = 1.2;

        /// @brief Compute the minimum expansion region radius.
        ///
        /// Used to set the m_expansionRegionMinimum member.
        static double CalculateExpansionRegion(const BoundingBox& gridBoundingBox,
                                               const BoundingBox& entityBoundingBox);

        /// @brief The original grid before expansion
        const CurvilinearGrid& m_originalGrid;

        /// @brief The grid to which the expansion is to be applied.
        const CurvilinearGrid& m_snappedGrid;

        /// @brief The minimum expansion region radius
        double m_expansionRegionMinimum{0.0};
    };

} // namespace meshkernel
