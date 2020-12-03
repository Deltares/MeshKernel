//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2020.
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

#include <vector>

#include <MeshKernel/Entities.hpp>

namespace meshkernel
{
    /// @brief A class describing polygons
    class Polygons
    {
    public:
        /// @brief Default constructor
        Polygons() = default;

        /// @brief Constructor
        /// @param[in] polygon The polygon nodes
        /// @param[in] projection The projection to use
        Polygons(const std::vector<Point>& polygon,
                 Projection projection);

        /// @brief Creates points inside the polygon using triangulation (the edges size determines how many points will be generated)
        /// @returns The generated points
        std::vector<std::vector<Point>> ComputePointsInPolygons() const;

        /// @brief Refines the polygon edges with additional nodes, from the start to the end index (refinepolygonpart)
        /// @param[in] startIndex The start index
        /// @param[in] endIndex The end index
        /// @param[in] refinementDistance The chosen refinement distance
        /// @return refinedPolygon The computed polygon
        std::vector<Point> RefineFirstPolygon(int startIndex, int endIndex, double refinementDistance) const;

        /// @brief Makes a new polygon from an existing one, by offsetting it by a distance (copypol)
        /// @param[in] distance The offset distance
        /// @param[in] Inner Offset inwards or outward
        /// @return The new offset polygon
        Polygons OffsetCopy(double distance, bool Inner) const;

        /// @brief Checks if a point is included in a given polygon
        /// @param[in] point The point to check
        /// @param[in] polygonIndex The index of the polygon to account for
        /// @return True if it is included, false otherwise
        bool IsPointInPolygon(Point point, int polygonIndex) const;

        /// @brief Checks if a point is included in any of the polygons (dbpinpol_optinside_perpol)
        /// @param[in] point The point to check
        /// @return True if it is included, false otherwise
        bool IsPointInPolygons(Point point) const;

        /// @brief Checks if the polygon is empty
        /// @return True if it is empty, false otherwise
        bool IsEmpty() const;

        /// @brief Gets the number of polygon nodes
        /// @return the number of polygon nodes
        [[nodiscard]] auto GetNumNodes() const { return m_nodes.size(); }

        std::vector<Point> m_nodes;                 ///< The polygon nodes
        Projection m_projection;                    ///< The current projection
        std::vector<std::vector<size_t>> m_indices; ///< Start-end indices of each polygon in m_nodes

    private:
        /// @brief Computes the perimeter of a closed polygon
        /// @param[in] polygonNodes The polygon nodes to use in the computation
        /// @return perimeter The computed polygon perimeter
        double PerimeterClosedPolygon(const std::vector<Point>& polygonNodes) const;

        /// @brief Computes the lengths of the polygon edges
        /// @param[in] polygonNodes The polygon nodes to use in the computation
        /// @return edgeLengths The length of each polygon edge
        std::vector<double> PolygonEdgeLengths(const std::vector<Point>& polygonNodes) const;

        /// @brief Computes the maximum edge length
        /// @param[in] polygonNodes The polygon to use in the computation
        /// @return maximumEdgeLength The maximum edge length of the polygon
        double MaximumEdgeLength(const std::vector<Point>& polygonNodes) const;
    };
} // namespace meshkernel
