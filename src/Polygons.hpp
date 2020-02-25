#pragma once

#include <vector>
#include "Mesh.hpp"

namespace GridGeom
{
    class Polygons
    {
    public:

        Polygons();

        bool Set(const std::vector<Point>& polygon);

        /// copynetboundstopol
        bool MeshBoundaryToPolygon(const Mesh& mesh,
            int counterClockWise,
            int setMeshState,
            std::vector<Point>& meshBoundaryPolygon,
            int& numNodesBoundaryPolygons);

        /// create a set of points in a polygon 
        bool CreatePointsInPolygons(std::vector<std::vector<Point>>);

        /// perimeter closed polygon
        bool PerimeterClosedPolygon(const std::vector<Point>& localPolygon,  int numPoints, double& perimeter);

        /// maximum edge length of a given polygon
        bool MaximumEdgeLength(const std::vector<Point>& localPolygon, int numPoints, double& maximumEdgeLength);

        std::vector<Point> m_nodes;             // Polygon nodes
        int m_numNodes;                         // NPL
        int m_numAllocatedNodes;                // MAXPOL
        int m_allocationSize = 100;
        Projections m_projection;

    private:

        bool WalkBoundary(const Mesh& mesh,
            std::vector<bool>& isVisited,
            int& numNodesBoundaryPolygon,
            int& currentNode,
            int meshBoundaryPolygonSize,
            std::vector<Point>& meshBoundaryPolygon);

    };

}