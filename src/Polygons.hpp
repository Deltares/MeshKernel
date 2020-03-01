#pragma once
#include "GeometryListNative.hpp"

#include <vector>
extern "C"
{
    void TRICALL(int *jatri, double* xs, double* ys, int* ns, int* indx, int* numtri, int* edgeidx, int* numedge, int* triedge, double* xs3, double* ys3, int* ns3, double* trisize);
}

namespace GridGeom
{
    class Mesh;

    class Polygons
    {
    public:

        Polygons();

        bool Set(const std::vector<Point>& polygon, Projections projection);

        bool Set(const GridGeomApi::GeometryListNative& geometryListNative, Projections projection);

        /// copynetboundstopol
        bool MeshBoundaryToPolygon(const Mesh& mesh,
            int counterClockWise,
            int setMeshState,
            std::vector<Point>& meshBoundaryPolygon,
            int& numNodesBoundaryPolygons);

        /// create a set of points in a polygon 
        bool CreatePointsInPolygons(std::vector<std::vector<Point>>& generatedPoints);

        std::vector<Point> m_nodes;             // Polygon nodes
        int m_numNodes;                         // NPL
        int m_numAllocatedNodes;                // MAXPOL
        int m_allocationSize = 100;
        Projections m_projection;

        /// perimeter closed polygon
        bool PerimeterClosedPolygon(const std::vector<Point>& localPolygon, int numPoints, double& perimeter);

    private:

        /// maximum edge length of a given polygon
        bool MaximumEdgeLength(const std::vector<Point>& localPolygon, int numPoints, double& maximumEdgeLength);

        bool WalkBoundary(const Mesh& mesh,
            std::vector<bool>& isVisited,
            int& numNodesBoundaryPolygon,
            int& currentNode,
            int meshBoundaryPolygonSize,
            std::vector<Point>& meshBoundaryPolygon);

    };

}