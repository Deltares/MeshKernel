#pragma once

#include <vector>
#include "Entities.hpp"

namespace GridGeom
{
    class Polygons;
    class Mesh;

    class LandBoundaries
    {

    public:

        LandBoundaries();

        /// admin_landboundary_segments
        bool Set(const std::vector<Point>& landBoundary);

        /// admin_landboundary_segments
        /// The land boundary will be split into segments that are within the polygon, and either close or not to the mesh boundary
        /// TODO: ? Why splitting in two segments is required?
        bool Administrate(Mesh& mesh, Polygons& polygons);

        /// find_nearest_meshline
        /// Find the mesh boundary line closest to the land boundary
        bool FindNearestMeshBoundary(const Mesh& mesh, const Polygons& polygon, int snapping);

        /// snap_to_landboundary
        /// snap netnodes to land boundary segment
        bool SnapMeshToLandBoundaries(Mesh& mesh);

        std::vector<int> m_meshNodesLandBoundarySegments; // lanseg_map mesh nodes to land boundary

    private:

        /// connect_boundary_paths, build an additional boundary for not assigned nodes  
        bool AssignSegmentsToAllMeshNodes(const Mesh& mesh, int edgeIndex, bool initialize, std::vector<int>& nodes, int numNodes);

        /// add_land, add new land boundary segment that connects two others
        bool AddLandBoundary(const Mesh& mesh, const std::vector<int>& nodesLoc, int numNodesLoc, int nodeIndex);

        /// make_path
        /// Assigns to each mesh node a land boundary a segment index (m_nodeLandBoundarySegments) 
        bool MakePath(const Mesh& mesh,
            const Polygons& polygons,
            int landBoundarySegment,
            bool meshBoundOnly,
            int& numNodesInPath,
            int& numRejectedNodesInPath);

        /// masknodes
        /// mask the mesh nodes to be considered in the shortest path algorithm for the current segmentIndex
        /// is setting leftIndex, rightIndex, leftEdgeRatio, rightEdgeRatio 
        bool ComputeMask(const Mesh& mesh,
            const Polygons& polygon,
            int segmentIndex,
            bool meshBoundOnly,
            int& leftIndex,
            int& rightIndex,
            double& leftEdgeRatio,
            double& rightEdgeRatio,
            int& startLandBoundaryIndex,
            int& endLandBoundaryIndex);

        /// maskcells
        /// mask the faces that are intersected by the land boundary
        bool MaskFaces(const Mesh& mesh,
            const bool& meshBoundOnly,
            std::vector<int>& landBoundaryFaces,
            int startNodeLandBoundaryIndex,
            int endNodeLandBoundaryindex,
            int& leftIndex,
            int& rightIndex,
            double& leftEdgeRatio,
            double& rightEdgeRatio);

        /// linkcrossedbyland
        /// check if a mesh edge is close to a land boundary segment
        bool IsMeshEdgeCloseToLandBoundaries(const Mesh& mesh,
            int edgeIndex,
            int startNodeLandBoundaryIndex,
            int endNodeLandBoundaryIndex,
            bool meshBoundOnly,
            int& leftIndex,
            int& rightIndex,
            double& leftEdgeRatio,
            double& rightEdgeRatio,
            int& landBoundaryNode);

        /// get_kstartend2
        /// Finds the start and end mesh node. These are the nodes that are
        /// on a edge close to the land boundary segment
        bool FindStartEndMeshNodes(const Mesh& mesh,
            const Polygons& polygons,
            int startLandBoundaryIndex,
            int endLandBoundaryIndex,
            int leftIndex,
            int rightIndex,
            double leftEdgeRatio,
            double rightEdgeRatio,
            int& startMeshNode,
            int& endMeshNode);

        /// Shortest_path
        /// connect mesh nodes starting from startMeshNode, using Dijkstra's shortest path algorithm
        /// the distance of each edge is the edge length multiplied by the distance from the land boundary
        bool ShortestPath(const Mesh& mesh, const Polygons& polygons, int landBoundarySegment,
            int startLandBoundaryIndex, int endLandBoundaryIndex, int startMeshNode, bool meshBoundOnly, std::vector<int>& connectedNodes);

        /// toland, compute the nearest point on the land boundary
        bool NearestLandBoundaryNode(const Projections& projection,
            const Point& node,
            int startLandBoundaryIndex,
            int endLandBoundaryIndex,
            double& minimumDistance,
            Point& pointOnLandBoundary,
            int& nearestLandBoundaryNodeIndex,
            double& edgeRatio);

        /// cellcrossedbyland
        /// TODO: it could be moved to generic operations
        bool IsFaceCrossedByLandBoundaries(const Mesh& mesh, int face, int startLandBoundaryIndex, int endLandBoundaryIndex);

        std::vector<Point> m_nodes;                       // XLAN, YLAN, ZLAN
        int m_numAllocatedNodes;                          // MAXLAN
        int m_numNode;                                    // actual MXLAN
        int m_numNodesLoc;                                // MXLAN_loc

        int m_numSegments = 0;                            // Nlanseg, number of land boundary segments 
        std::vector<std::vector<int>> m_segmentIndices;   // lanseg_startend
        std::vector<std::vector<double>> m_nodesLand;     // !node to land boundary segment mapping

        std::vector<int> m_nodeMask;                      // nodemask, masking the net nodes
        std::vector<int> m_faceMask;                      // masking faces
        std::vector<int> m_edgeMask;                      // masking edges

        bool m_landMask = true;
        bool m_addLandboundaries = true;
        int m_numFacesMasked = 0;
        int m_maskDepth = 0;

        // caches
        std::vector<Point> m_polygonNodesCache;          // array of points (e.g. points of a face)
        std::vector<double> m_nodesMinDistances;
        const int m_allocationSize = 10000;               // allocation size for allocateVector

        // Parameters
        const double m_closeToLandBoundaryFactor = 5.0;   // close - to - landboundary tolerance, measured in number of meshwidths
        const double m_closeWholeMeshFactor = 1.0;        // close - to - landboundary tolerance, measured in number of meshwidths
        const double m_minDistanceFromLandFactor = 2.0;
    };

}