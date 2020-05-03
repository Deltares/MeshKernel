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
        int m_numNodes;                            // NPL
        int m_numAllocatedNodes;                   // MAXPOL
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