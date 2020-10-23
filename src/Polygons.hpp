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
#include "Entities.hpp"

namespace meshkernel
{
    class Mesh;

    class Polygons
    {
    public:
        Polygons();

        Polygons(const std::vector<Point>& polygon,
                 Projections projection);

        /// <summary>
        /// (copynetboundstopol)
        /// </summary>
        bool MeshBoundaryToPolygon(Mesh& mesh,
                                   std::vector<Point>& meshBoundaryPolygon,
                                   int& numNodesBoundaryPolygons);

        /// <summary>
        /// create a set of points in a polygon
        /// </summary>
        /// <param name="generatedPoints"></param>
        /// <returns></returns>
        bool CreatePointsInPolygons(std::vector<std::vector<Point>>& generatedPoints);

        std::vector<Point> m_nodes; // Polygon nodes
        Projections m_projection;
        int m_numNodes = 0;                      // NPL
        int m_numAllocatedNodes = 0;             // MAXPOL
        std::vector<std::vector<int>> m_indices; // start-end of polygon nodes in m_nodes
        int m_allocationSize = 100;

        /// perimeter closed polygon
        bool PerimeterClosedPolygon(const std::vector<Point>& localPolygon, int numPoints, double& perimeter);

        /// refinepolygonpart
        bool RefinePolygonPart(int startIndex, int endIndex, double refinementDistance, std::vector<Point>& refinedPolygon);

        /// refinepolygonpart
        bool PolygonEdgeLengths(const std::vector<Point>& localPolygon, std::vector<double>& edgeLengths) const;

        ///copypol, copy and move a polygon orthogonally
        bool OffsetCopy(double distance, bool Inner, Polygons& newPolygon);

        int GetNumNodes() const { return m_numNodes; }

        bool IsPointInPolygon(const Point& point, int polygonIndex) const;

        // dbpinpol_optinside_perpol
        bool IsPointInPolygons(const Point& point) const;

    private:
        /// maximum edge length of a given polygon
        bool MaximumEdgeLength(const std::vector<Point>& localPolygon, int numPoints, double& maximumEdgeLength);

        bool WalkBoundaryFromNode(const Mesh& mesh,
                                  std::vector<bool>& isVisited,
                                  int& nodeIndex,
                                  int& currentNode,
                                  std::vector<Point>& meshBoundaryPolygon) const;
    };

} // namespace meshkernel
