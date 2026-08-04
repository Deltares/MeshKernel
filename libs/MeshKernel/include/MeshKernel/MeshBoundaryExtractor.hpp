//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2026.
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

#include <tuple>
#include <vector>

#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/Point.hpp"

namespace meshkernel
{

    /// @brief Extract the boundary polygon from the mesh
    class MeshBoundaryExtractor
    {
    public:
        /// @brief Extract all boundaries as a single sequence of points, separated by an invalid point
        static std::vector<Point> ExtractConcatenated(const Mesh2D& mesh, BoundarySelection boundaryType = BoundarySelection::All);

        /// @brief Extract all boundaries keeping them separated and
        ///
        /// The result consists of an array of each of the boundary polygons
        /// Additionally, an array indicating if the boundary polygon is a exterior boundary or not.
        /// True => is-exterior, False => otherwise
        static std::tuple<std::vector<std::vector<Point>>, std::vector<bool>> Extract(const Mesh2D& mesh);

    private:
        /// @brief Temporary struct, used when computing the boundaries
        struct BoundaryEdge
        {
            UInt edgeId;        ///< Id of edge
            UInt neighbourNode; ///< Id of node at opposite end of edge
            UInt leftFace;      ///< Store face mapping on the edge for easy retrieval during polygon trace
            double angle;       ///< Angle of the edge pointing away from the pivot node
        };

        /// @brief Ensure the angle lies between 0 .. 2 pi.
        static double NormalizeAngle(double angle);

        /// @brief Find the edge that has the smallest angle [0 .. 2pi), to the incident edge
        static UInt FindEdgeWithMinumumAngle(const std::vector<BoundaryEdge>& boundaryEdges,
                                             const std::vector<bool>& edgeVisited,
                                             const double incomingAngle);

        /// @brief Construct mapping from node-id to all impinging edges for boundary all edges
        static void FindAllBoundarEdges(const std::vector<Point>& nodes,
                                        const std::vector<Edge>& edges,
                                        const std::vector<std::array<UInt, 2>>& edgesFaces,
                                        std::unordered_map<UInt, std::vector<BoundaryEdge>>& boundaryAdjacency);

        /// @brief Append the boundary polygon to the set of all boundary polygon
        ///
        /// The boundary polygons may be reversed if found to be in clockwise direction
        static void Append(const Point& centre,
                           const Projection projection,
                           std::vector<Point>& boundaryPolygon,
                           std::vector<bool>& isExterior,
                           std::vector<std::vector<Point>>& separatedBoundaryPolygons);

        /// @brief Find boundary loops
        ///
        /// Any boundary loops found may need to be processed further as they may themselves contain sub-loops
        static void FindBoundaryPolygons(const std::vector<Point>& nodes,
                                         const std::vector<Edge>& edges,
                                         const std::vector<std::array<UInt, 2>>& edgesFaces,
                                         std::vector<std::vector<Point>>& allPolygons,
                                         std::vector<std::vector<UInt>>& allTouchedFaces);

        /// @brief Separate polygons that contains multiple sub-polygons and determine externality
        ///
        /// allBoundaryPolygons is not const because it may be updated.
        /// @note allBoundaryPolygons should not be accessed after calling this function
        static std::tuple<std::vector<std::vector<meshkernel::Point>>, std::vector<bool>>
        SeparateAndDetermineExternality(const Mesh2D& mesh,
                                        std::vector<std::vector<Point>>& allBoundaryPolygons,
                                        const std::vector<std::vector<UInt>>& allTouchedFaces);
    };

} // namespace meshkernel
