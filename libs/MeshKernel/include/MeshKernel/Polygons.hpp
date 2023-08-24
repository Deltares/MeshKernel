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

#include <unordered_map>
#include <vector>

#include <MeshKernel/BoundingBox.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/Polygon.hpp>
#include <MeshKernel/PolygonalEnclosure.hpp>

namespace meshkernel
{

    /// @brief Forward declaration of the LandBoundary
    class LandBoundary;

    /// @brief A class containing a list of polygonaly enclosed regions
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

        /// @brief Get index of the polygon, and map the start- and end-index to the start- and end-index of the local polygon nodes
        // This is bollocks
        // what about nodes for the interior polygons
        std::tuple<UInt, UInt, UInt> PolygonIndex(UInt startIndex, UInt endIndex) const;

        /// @brief Get the polygonal enclosure at the index
        /// @note Will throw ConstraintError exception if polygons.IsEmpty = true or index is out of range
        const PolygonalEnclosure& Enclosure(const UInt index) const;

        /// @brief Creates points inside the polygon using triangulation (the edges size determines how many points will be generated)
        /// @returns The generated points
        [[nodiscard]] std::vector<std::vector<Point>> ComputePointsInPolygons() const;

        /// @brief Refines the polygon edges with additional nodes, from the start to the end index (refinepolygonpart)
        /// @param[in] startIndex The start index
        /// @param[in] endIndex The end index
        /// @param[in] refinementDistance The chosen refinement distance
        /// @return refinedPolygon The computed polygon
        [[nodiscard]] std::vector<Point> RefineFirstPolygon(UInt startIndex, UInt endIndex, double refinementDistance) const;

        /// @brief Refines the polygon edges with additional nodes, from the start to the end index (refinepolygonpart)
        /// @param[in] polygonIndex The polygon index
        /// @param[in] startIndex The start index for the node array for the polygon
        /// @param[in] endIndex The end index for the node array for the polygon
        /// @param[in] refinementDistance The chosen refinement distance
        /// @return refinedPolygon The computed polygon
        [[nodiscard]] std::vector<Point> RefinePolygon(UInt polygonIndex, UInt startIndex, UInt endIndex, double refinementDistance) const;

        /// @brief Makes a new polygon from an existing one, by offsetting it by a distance (copypol)
        /// @param[in] distance The offset distance
        /// @param[in] innerAndOuter Offset inwards or outward
        /// @return The new offset polygon
        [[nodiscard]] Polygons OffsetCopy(double distance, bool innerAndOuter) const;

        /// @brief Snap the polygon to the land boundary
        ///
        /// The polygon points are snapped to the closest point on the land boundary.
        /// @param[in] landBoundary The land boundary to which the polygon should be snapped.
        /// @param[in] startIndex The start index
        /// @param[in] endIndex The end index
        void SnapToLandBoundary(const LandBoundary& landBoundary, UInt startIndex, UInt endIndex);

        /// @brief Checks if a point is included in a given polygon.
        /// When the polygon is empty, the point is always included by default
        /// @param[in] point The point to check
        /// @param[in] polygonIndex The index of the polygon to account for
        /// @return True if it is included, false otherwise
        [[nodiscard]] bool IsPointInPolygon(Point const& point, UInt polygonIndex) const;

        // TODO can reduce the result of this function to only a UInt (valid value => found, invalid => not found)
        /// @brief Checks if a point is included in any of the polygons (dbpinpol_optinside_perpol)
        /// @param[in] point The point to check
        /// @return The index of a polygon where the point is included or if none has been found, constants::missing::sizetValue
        [[nodiscard]] std::tuple<bool, UInt> IsPointInPolygons(const Point& point) const;

        /// @brief For each point, compute the index of the polygon including it
        /// @param[in] point The vector of points
        /// @return A vector of booleans to indicate if the point is in polygon
        [[nodiscard]] std::vector<bool> PointsInPolygons(const std::vector<Point>& point) const;

        /// @brief Checks if the polygon is empty
        /// @return True if it is empty, false otherwise
        [[nodiscard]] bool IsEmpty() const;

        /// @brief Gives the number of polygons
        /// @return Number of polygons
        [[nodiscard]] UInt GetNumPolygons() const;

        /// @brief Gets the number of polygon nodes
        /// @return The number of polygon nodes
        [[nodiscard]] size_t GetNumNodes() const;

        /// @brief Gets the projection
        /// @return The projection
        [[nodiscard]] Projection GetProjection() const { return m_projection; }

        /// @brief Gets the nodes of the polygon
        /// @return Vector of nodes of the polygon
        [[nodiscard]] std::vector<Point> const& Nodes() const { return m_nodes; }

        /// @brief Gets the coordinates of a node by index
        /// @param[in] i Node index
        /// @return Node coordinates
        [[nodiscard]] Point const& Node(UInt i) const { return m_nodes[i]; }

        /// @brief Gets the bounding box for the polygon index i
        /// @param[in] polygonIndex Outer polygon index
        /// @return The bounding box
        [[nodiscard]] BoundingBox GetBoundingBox(UInt polygonIndex) const;

    private:
        std::vector<PolygonalEnclosure> m_enclosures;                ///< List of polygons
        std::vector<Point> m_nodes;                                  ///< The polygon nodes
        Projection m_projection;                                     ///< The current projection
        std::vector<std::pair<UInt, UInt>> m_outer_polygons_indices; ///< Start-end indices of each outer polygon in m_nodes

        /// @brief Computes the perimeter of a closed polygon
        /// @param[in] polygonNodes The polygon nodes to use in the computation
        /// @return perimeter The computed polygon perimeter
        double PerimeterClosedPolygon(const std::vector<Point>& polygonNodes) const;

        /// @brief Computes the lengths of the polygon edges
        /// @param[in] polygonNodes The polygon nodes to use in the computation
        /// @return edgeLengths The length of each polygon edge
        std::vector<double> PolygonEdgeLengths(const std::vector<Point>& polygonNodes) const;
    };
} // namespace meshkernel

inline const meshkernel::PolygonalEnclosure& meshkernel::Polygons::Enclosure(const UInt index) const
{
    if (IsEmpty())
    {
        throw ConstraintError(VariadicErrorMessage("Enclosures list is empty."));
    }

    if (index >= m_enclosures.size())
    {
        throw ConstraintError(VariadicErrorMessage("Invalid enclosure index: {}, maximum index: {}",
                                                   index, m_enclosures.size() - 1));
    }

    return m_enclosures[index];
}
