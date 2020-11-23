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
    class Mesh;

    class Polygons
    {
    public:
        /// @brief Default constructor
        Polygons();

        /// @brief Constructor
        /// @param[in] polygon The polygon nodes
        /// @param[in] projection The projection to use
        Polygons(const std::vector<Point>& polygon,
                 Projections projection);

        /// @brief Convert all mesh boundaries to a polygon, including holes (copynetboundstopol)
        /// @param[in] mesh The input mesh
        /// @returns meshBoundaryPolygon The resulting polygon mesh boundary
        std::vector<Point> MeshBoundaryToPolygon(Mesh& mesh) const;

        /// @brief Creates points inside the polygon using triangulation
        /// @param[out] generatedPoints the generated points
        std::vector<std::vector<Point>> ComputePointsInPolygons() const;

        /// @brief Refines the polygon edges with additional nodes, from the start to the end index (refinepolygonpart)
        /// @param[in] startIndex The start index
        /// @param[in] endIndex The end index
        /// @param[in] refinementDistance The chosen refinement distance
        /// @param[out] refinedPolygon The computed polygon
        void RefinePolygonPart(int startIndex, int endIndex, double refinementDistance, std::vector<Point>& refinedPolygon);

        /// @brief Makes a new polygon from an existing one, by offsetting it by a distance (copypol)
        /// @param[in] distance The offset distance
        /// @param[in] Inner Inner or outer polygon offset
        /// @param[out] newPolygon the new polygon
        void OffsetCopy(double distance, bool Inner, Polygons& newPolygon);

        /// @brief Checks if a point is included in a given polygon
        /// @param[in] point The point to check
        /// @param[in] polygonIndex The index of the polygon to account for
        /// @return True if it is included, false otherwise
        bool IsPointInPolygon(const Point& point, int polygonIndex) const;

        /// @brief Checks if a point is included in any of the polygons (dbpinpol_optinside_perpol)
        /// @param[in] point The point to check
        /// @return True if it is included, false otherwise
        bool IsPointInPolygons(const Point& point) const;

        /// @brief Checks if the polygon is empty
        /// @return True if it is empty, false otherwise
        bool Polygons::IsEmpty() const;

        /// @brief Gets the number of polygon nodes
        /// @return the number of polygon nodes
        [[nodiscard]] auto GetNumNodes() const { return m_nodes.size(); }

        std::vector<Point> m_nodes;                 // The polygon nodes
        Projections m_projection;                   // The current projection
        std::vector<std::vector<size_t>> m_indices; // Start-end of polygon nodes in m_nodes

    private:
        /// @brief Computes the perimeter of a closed polygon
        /// @param[in] polygonNodes The polygon nodes to use in the computation
        /// @return perimeter The computed perimeter
        double PerimeterClosedPolygon(const std::vector<Point>& polygonNodes) const;

        /// @brief Computes the lengths of the polygon edges
        /// @param[in] polygonNodes The polygon nodes to use in the computation
        /// @return edgeLengths The length of each edge
        std::vector<double> PolygonEdgeLengths(const std::vector<Point>& polygonNodes) const;

        /// @brief Computes the maximum edge length
        /// @param[in] localPolygon The polygon to use in the computation
        /// @param[in] numPoints The number of polygon points
        /// @param[out] maximumEdgeLength The maximum edge length
        void MaximumEdgeLength(const std::vector<Point>& localPolygon, size_t numPoints, double& maximumEdgeLength) const;

        /// @brief Constructs a polygon from the meshboundary, by walking through the mesh
        /// @param[in] mesh The input mesh
        /// @param[in] isVisited the visited mesh nodes
        /// @param[in] nodeIndex the node where to initialize the algorithm
        /// @param[in] currentNode the current node
        /// @param[out] meshBoundaryPolygon The resulting polygon points
        void WalkBoundaryFromNode(const Mesh& mesh,
                                  std::vector<bool>& isVisited,
                                  int& currentNode,
                                  std::vector<Point>& meshBoundaryPolygon) const;
    };
} // namespace meshkernel
