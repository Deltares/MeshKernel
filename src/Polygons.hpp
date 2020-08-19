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
#include <algorithm>

#include "Entities.hpp"
#include "Operations.cpp"

extern "C"
{
    void Triangulation(int *jatri, double* xs, double* ys, int* ns, int* indx, int* numtri, int* edgeidx, int* numedge, int* triedge, double* xs3, double* ys3, int* ns3, double* trisize);
}

namespace GridGeom
{
    class Mesh;

    class Polygons
    {
    public:

        Polygons();

        bool Set(const std::vector<Point>& polygon, Projections projection);

        /// copynetboundstopol
        bool MeshBoundaryToPolygon(Mesh& mesh,
            int counterClockWise,
            std::vector<Point>& meshBoundaryPolygon,
            int& numNodesBoundaryPolygons);

        /// create a set of points in a polygon 
        bool CreatePointsInPolygons(std::vector<std::vector<Point>>& generatedPoints);

        std::vector<Point> m_nodes;                // Polygon nodes
        int m_numNodes = 0;                        // NPL
        int m_numAllocatedNodes = 0;               // MAXPOL
        std::vector<std::vector<int>> m_indexses;  // start-end of polygon nodes in m_nodes
        int m_allocationSize = 100;
        Projections m_projection;

        /// perimeter closed polygon
        bool PerimeterClosedPolygon(const std::vector<Point>& localPolygon, int numPoints, double& perimeter);

        /// refinepolygonpart
        bool RefinePart(int startIndex, int endIndex, double refinementDistance, std::vector<Point>& refinedPolygon);

        /// refinepolygonpart
        bool EdgeLengths(const std::vector<Point>& localPolygon, std::vector<double>& edgeLengths );

        ///copypol, copy and move a polygon orthogonally
        bool OffsetCopy(int nodeIndex, double distance, bool Inner, Polygons& newPolygon);

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
            std::vector<Point>& meshBoundaryPolygon);

    };

}