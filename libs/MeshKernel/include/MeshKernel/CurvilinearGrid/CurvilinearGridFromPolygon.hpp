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

namespace meshkernel
{
    // Forward declarations
    class Mesh2D;
    class Polygon;
    class CurvilinearGrid;

    /// @brief A class used to compute a CurvilinearGrid from a polygon
    class CurvilinearGridFromPolygon
    {
    public:
        /// @param polygon The input polygon
        explicit CurvilinearGridFromPolygon(const Polygon& polygon);

        /// @brief Compute curvilinear in a polygon (pol2curvi)
        /// @returns The computed curvilinear grid
        std::unique_ptr<CurvilinearGrid> Compute(UInt firstNode, UInt secondNode, UInt thirdNode, bool useFourthSide) const;

        /// @brief Compute curvilinear in a triangle (pol2curvi_tri)
        /// @returns The computed curvilinear grid
        std::unique_ptr<CurvilinearGrid> Compute(UInt firstNode, UInt secondNode, UInt thirdNode) const;

    private:
        /// &brief Fill boundary coordinates
        void AssignPolygonPointsToSegment(UInt nodeIndex,
                                          UInt numPointsSide,
                                          int direction,
                                          std::vector<Point>& sideToFill) const;

        void ComputeNumberOfMNodes(const UInt firstNode,
                                   const UInt secondNode,
                                   const UInt numPolygonNodes,
                                   int& direction,
                                   UInt& numMNodes) const;

        void ComputeNumberOfNNodes(const UInt secondNode,
                                   const UInt thirdNode,
                                   const UInt numPolygonNodes,
                                   const int direction,
                                   UInt& numNNodes) const;

        const Polygon& m_polygon; /// Reference to a polygon
    };

} // namespace meshkernel
