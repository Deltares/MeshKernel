//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2022.
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

    class MeshSnappingCalculator
    {
    public:
        virtual ~MeshSnappingCalculator() = default;

        virtual double compute(const CurvilinearGridNodeIndices& currentPointIndex,
                               const CurvilinearGridNodeIndices& gridLineIndex) const = 0;
    };

    class SmeerFunctie : public MeshSnappingCalculator
    {
    public:
        SmeerFunctie(const CurvilinearGrid& grid,
                     const CurvilinearGridNodeIndices& lowerLeft,
                     const CurvilinearGridNodeIndices& upperRight,
                     const CurvilinearGridNodeIndices& regionSize);

        double compute(const CurvilinearGridNodeIndices& currentPointIndex,
                       const CurvilinearGridNodeIndices& gridLineIndex) const override;

    private:
        const CurvilinearGrid& m_grid;
        CurvilinearGridNodeIndices m_indexBoxLowerLeft;
        CurvilinearGridNodeIndices m_indexBoxUpperRight;
        CurvilinearGridNodeIndices m_snappedRegionSize;
    };

    class NotSmeerFunctie : public MeshSnappingCalculator
    {
    public:
        NotSmeerFunctie(const CurvilinearGrid& originalGrid,
                        const CurvilinearGrid& snappedGrid,
                        const LandBoundary& landBoundary);

        double compute(const CurvilinearGridNodeIndices& currentPointIndex,
                       const CurvilinearGridNodeIndices& gridLineIndex) const override;

    private:
        const CurvilinearGrid& m_originalGrid;
        const CurvilinearGrid& m_snappedGrid;
        double m_dsix;
    };

    class CurvilinearGridSnapping : public CurvilinearGridAlgorithm
    {
    public:
        /// @brief Class constructor
        /// @param[in] grid                        The input curvilinear grid
        CurvilinearGridSnapping(std::shared_ptr<CurvilinearGrid> grid,
                                const LandBoundary& lb,
                                const std::vector<Point>& points);

        /// @brief Executes the snapping algorithm
        CurvilinearGrid Compute() override;

        // TODO How to add this?
        // 1. two different classes, constructed with land bondary or spline
        // 2. two constructors (or snap to functions) that create an adaptor.

        // void SnapTo(const LandBoundary& landboundary);

        // void SnapTo(const Splines& splines,
        //             const UInt whichSpline);

    protected:
        /// @brief Snap the point to the object.
        // TODO what is best to do here?
        // 1. virtual function
        // 2. Adaptor for LandBoundary and Splines
        // virtual Point SnapPoint(const Point& pnt) const = 0;

    private:
        std::tuple<CurvilinearGridNodeIndices, CurvilinearGridNodeIndices> ComputeLoopBounds(const CurvilinearGridNodeIndices& currentNodeIndex,
                                                                                             const CurvilinearGridNodeIndices& snappedRegionSize) const;

        void ModifyField(const CurvilinearGridNodeIndices& currentNodeIndex,
                         const CurvilinearGridNodeIndices& snappedRegionSize,
                         const MeshSnappingCalculator& updateFactor);

        void Initialise();

        // SHould this be a shared or raw poitner or a reference.
        const CurvilinearGrid& m_originalGrid;
        const LandBoundary& m_landBoundary;
        const std::vector<Point> m_points;

        // TODO change names, especially m_snappedRegionSize.
        CurvilinearGridNodeIndices m_lineStartIndex;
        CurvilinearGridNodeIndices m_lineEndIndex;
        CurvilinearGridNodeIndices m_snappedRegionSize;
        // is this the snapping region?
        CurvilinearGridNodeIndices m_indexBoxLowerLeft;
        CurvilinearGridNodeIndices m_indexBoxUpperRight;
    };

} // namespace meshkernel
