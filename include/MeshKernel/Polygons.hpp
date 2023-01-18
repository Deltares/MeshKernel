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

#include <MeshKernel/Entities.hpp>
#include <unordered_map>

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
        [[nodiscard]] std::vector<std::vector<Point>> ComputePointsInPolygons() const;

        /// @brief Refines the polygon edges with additional nodes, from the start to the end index (refinepolygonpart)
        /// @param[in] startIndex The start index
        /// @param[in] endIndex The end index
        /// @param[in] refinementDistance The chosen refinement distance
        /// @return refinedPolygon The computed polygon
        [[nodiscard]] std::vector<Point> RefineFirstPolygon(size_t startIndex, size_t endIndex, double refinementDistance) const;

        /// @brief Makes a new polygon from an existing one, by offsetting it by a distance (copypol)
        /// @param[in] distance The offset distance
        /// @param[in] innerAndOuter Offset inwards or outward
        /// @return The new offset polygon
        [[nodiscard]] Polygons OffsetCopy(double distance, bool innerAndOuter) const;

        /// @brief Checks if a point is included in a given polygon.
        /// When the polygon is empty, the point is always included by default
        /// @param[in] point The point to check
        /// @param[in] polygonIndex The index of the polygon to account for
        /// @return True if it is included, false otherwise
        [[nodiscard]] bool IsPointInPolygon(Point const& point, size_t polygonIndex) const;

        /// @brief Checks if a point is included in any of the polygons (dbpinpol_optinside_perpol)
        /// @param[in] point The point to check
        /// @return The index of a polygon where the point is included or if none has been found, constants::missing::sizetValue
        [[nodiscard]] std::tuple<bool, size_t> IsPointInPolygons(Point point) const;

        /// @brief For each point, compute the index of the polygon including it
        /// @param[in] point The vector of points
        /// @return A vector of booleans to indicate if the point is in polygon
        [[nodiscard]] std::vector<bool> PointsInPolygons(const std::vector<Point>& point) const;

        /// @brief Checks if the polygon is empty
        /// @return True if it is empty, false otherwise
        bool IsEmpty() const;

        /// @brief Gives the number of polygons
        /// @return Number of polygons
        [[nodiscard]] size_t GetNumPolygons() const;

        /// @brief Gets the number of polygon nodes
        /// @return The number of polygon nodes
        [[nodiscard]] auto GetNumNodes() const { return m_nodes.size(); }

        /// @brief Gets the projection
        /// @return The projection
        [[nodiscard]] Projection GetProjection() const { return m_projection; }

        /// @brief Gets the start-end indices of each outer polygon
        /// @param[in] i Outer polygon index
        /// @return Pair of start and end indices
        [[nodiscard]] std::pair<size_t, size_t> const& OuterIndices(size_t i) const
        {
            return m_outer_polygons_indices[i];
        }

        /// @brief Gets the nodes of the polygon
        /// @return Vector of nodes of the polygon
        [[nodiscard]] std::vector<Point> const& Nodes() const { return m_nodes; }

        /// @brief Gets the coordinates of a node by index
        /// @param[in] i Node index
        /// @return Node coordinates
        [[nodiscard]] Point const& Node(size_t i) const { return m_nodes[i]; }

    private:
        std::vector<Point> m_nodes;                                                                  ///< The polygon nodes
        Projection m_projection;                                                                     ///< The current projection
        std::vector<std::pair<size_t, size_t>> m_outer_polygons_indices;                             ///< Start-end indices of each outer polygon in m_nodes
        std::unordered_map<size_t, std::vector<std::pair<size_t, size_t>>> m_inner_polygons_indices; ///< For each outer polygon, the indices of each inner polygon

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
