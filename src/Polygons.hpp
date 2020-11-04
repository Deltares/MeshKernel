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

extern "C"
{
    void Triangulation(int* jatri, double* xs, double* ys, int* ns, int* indx, int* numtri, int* edgeidx, int* numedge, int* triedge, double* xs3, double* ys3, int* ns3, double* trisize);
}

namespace meshkernel
{
    class Mesh;

    class Polygons
    {
    public:
        Polygons();

        Polygons(const std::vector<Point>& polygon,
                 Projections projection);

        /// @brief (copynetboundstopol)
        /// @param[out] meshBoundaryPolygon
        /// @param[out] numNodesBoundaryPolygons
        void MeshBoundaryToPolygon(Mesh& mesh,
                                   std::vector<Point>& meshBoundaryPolygon,
                                   int& numNodesBoundaryPolygons);

        /// @brief Create a set of points in a polygon
        /// @param[out] generatedPoints
        void CreatePointsInPolygons(std::vector<std::vector<Point>>& generatedPoints);

        std::vector<Point> m_nodes; // Polygon nodes
        Projections m_projection;
        int m_numNodes = 0;                      // NPL
        int m_numAllocatedNodes = 0;             // MAXPOL
        std::vector<std::vector<int>> m_indices; // start-end of polygon nodes in m_nodes
        int m_allocationSize = 100;

        /// @brief perimeter closed polygon
        void PerimeterClosedPolygon(const std::vector<Point>& localPolygon, int numPoints, double& perimeter);

        /// @brief refinepolygonpart
        void RefinePolygonPart(int startIndex, int endIndex, double refinementDistance, std::vector<Point>& refinedPolygon);

        /// @brief refinepolygonpart
        void PolygonEdgeLengths(const std::vector<Point>& localPolygon, std::vector<double>& edgeLengths) const;

        ///copypol, copy and move a polygon orthogonally
        void OffsetCopy(double distance, bool Inner, Polygons& newPolygon);

        [[nodiscard]] int GetNumNodes() const { return m_numNodes; }

        bool IsPointInPolygon(const Point& point, int polygonIndex) const;

        // dbpinpol_optinside_perpol
        bool IsPointInPolygons(const Point& point) const;

        /// @brief Checks if the polygon is empty
        /// @return Boolean whether the polygon is empty or not
        bool Polygons::IsEmpty() const;

    private:
        /// @brief maximum edge length of a given polygon
        void MaximumEdgeLength(const std::vector<Point>& localPolygon, int numPoints, double& maximumEdgeLength);

        // TODO: Document
        void WalkBoundaryFromNode(const Mesh& mesh,
                                  std::vector<bool>& isVisited,
                                  int& nodeIndex,
                                  int& currentNode,
                                  std::vector<Point>& meshBoundaryPolygon) const;
    };

} // namespace meshkernel
