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

namespace GridGeom
{
    class Polygons;
    class Mesh;

    class LandBoundaries
    {

    public:

        /// <summary>
        /// Default Ctor
        /// </summary>
        /// <returns></returns>
        LandBoundaries();

        /// <summary>
        /// Set the landboundaries parameters (admin_landboundary_segments)
        /// </summary>
        /// <param name="landBoundary"></param>
        /// <param name="mesh"></param>
        /// <param name="polygons"></param>
        /// <returns></returns>
        bool Set(const std::vector<Point>& landBoundary, 
                 std::shared_ptr<Mesh> mesh, 
                 std::shared_ptr<Polygons> polygons);

        /// <summary>
        /// The land boundary will be split into segments that are within the polygon, and either close or not to the mesh boundary (admin_landboundary_segments)
        /// TODO: ? Why splitting in two segments is required?
        /// </summary>
        /// <param name="polygons"></param>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool Administrate();

        /// <summary>
        /// Find the mesh boundary line closest to the land boundary (find_nearest_meshline)
        /// </summary>
        /// <param name="mesh"></param>
        /// <param name="polygon"></param>
        /// <param name="snapping"></param>
        /// <returns></returns>
        bool FindNearestMeshBoundary( int snapping);

        /// <summary>
        /// Snap mesh nodes to land boundaies (snap_to_landboundary)
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool SnapMeshToLandBoundaries();

        std::vector<int> m_meshNodesLandBoundarySegments; // lanseg_map, mesh nodes to land boundary mapping

    private:

        /// <summary>
        /// Build an additional boundary for not assigned nodes (connect_boundary_paths) 
        /// </summary>
        /// <param name="mesh"></param>
        /// <param name="edgeIndex"></param>
        /// <param name="initialize"></param>
        /// <param name="nodes"></param>
        /// <param name="numNodes"></param>
        /// <returns></returns>
        bool AssignSegmentsToAllMeshNodes(int edgeIndex, 
                                          bool initialize, 
                                          std::vector<int>& nodes, 
                                          int numNodes);
        /// <summary>
        /// add new land boundary segment that connects two others (add_land)
        /// </summary>
        /// <param name="mesh"></param>
        /// <param name="nodesLoc"></param>
        /// <param name="numNodesLoc"></param>
        /// <param name="nodeIndex"></param>
        /// <returns></returns>
        bool AddLandBoundary(const std::vector<int>& nodesLoc, 
                             int numNodesLoc, 
                             int nodeIndex);

        /// <summary>
        /// Assigns to each mesh node a land boundary a segment index (m_nodeLandBoundarySegments) 
        /// </summary>
        /// <param name="mesh"></param>
        /// <param name="polygons"></param>
        /// <param name="landBoundarySegment"></param>
        /// <param name="meshBoundOnly"></param>
        /// <param name="numNodesInPath"></param>
        /// <param name="numRejectedNodesInPath"></param>
        /// <returns></returns>
        bool MakePath(int landBoundarySegment,
                      bool meshBoundOnly,
                      int& numNodesInPath,
                      int& numRejectedNodesInPath);

        /// <summary>
        /// Mask the mesh nodes to be considered in the shortest path algorithm for the current segmentIndex.
        /// Is setting leftIndex, rightIndex, leftEdgeRatio, rightEdgeRatio.
        /// </summary>
        /// <param name="mesh"></param>
        /// <param name="polygon"></param>
        /// <param name="segmentIndex"></param>
        /// <param name="meshBoundOnly"></param>
        /// <param name="leftIndex"></param>
        /// <param name="rightIndex"></param>
        /// <param name="leftEdgeRatio"></param>
        /// <param name="rightEdgeRatio"></param>
        /// <param name="startLandBoundaryIndex"></param>
        /// <param name="endLandBoundaryIndex"></param>
        /// <returns></returns>
        bool ComputeMask(int segmentIndex,
                         bool meshBoundOnly,
                         int& leftIndex,
                         int& rightIndex,
                         double& leftEdgeRatio,
                         double& rightEdgeRatio,
                         int& startLandBoundaryIndex,
                         int& endLandBoundaryIndex);

        /// <summary>
        /// Mask the faces that are intersected by the land boundary (maskcells)
        /// </summary>
        /// <param name="mesh"></param>
        /// <param name="meshBoundOnly"></param>
        /// <param name="landBoundaryFaces"></param>
        /// <param name="startNodeLandBoundaryIndex"></param>
        /// <param name="endNodeLandBoundaryindex"></param>
        /// <param name="leftIndex"></param>
        /// <param name="rightIndex"></param>
        /// <param name="leftEdgeRatio"></param>
        /// <param name="rightEdgeRatio"></param>
        /// <returns></returns>
        bool MaskFaces(const bool& meshBoundOnly,
            std::vector<int>& landBoundaryFaces,
            int startNodeLandBoundaryIndex,
            int endNodeLandBoundaryindex,
            int& leftIndex,
            int& rightIndex,
            double& leftEdgeRatio,
            double& rightEdgeRatio);

        /// <summary>
        /// Check if a mesh edge is close to a land boundary segment (linkcrossedbyland)
        /// </summary>
        /// <param name="mesh"></param>
        /// <param name="edgeIndex"></param>
        /// <param name="startNodeLandBoundaryIndex"></param>
        /// <param name="endNodeLandBoundaryIndex"></param>
        /// <param name="meshBoundOnly"></param>
        /// <param name="leftIndex"></param>
        /// <param name="rightIndex"></param>
        /// <param name="leftEdgeRatio"></param>
        /// <param name="rightEdgeRatio"></param>
        /// <param name="landBoundaryNode"></param>
        /// <returns></returns>
        bool IsMeshEdgeCloseToLandBoundaries(int edgeIndex,
                                             int startNodeLandBoundaryIndex,
                                             int endNodeLandBoundaryIndex,
                                             bool meshBoundOnly,
                                             int& leftIndex,
                                             int& rightIndex,
                                             double& leftEdgeRatio,
                                             double& rightEdgeRatio,
                                             int& landBoundaryNode);

        /// <summary>
        /// Finds the start and end mesh node. These are the nodes that are on a edge close to the land boundary segment (get_kstartend2)
        /// </summary>
        /// <param name="mesh"></param>
        /// <param name="polygons"></param>
        /// <param name="startLandBoundaryIndex"></param>
        /// <param name="endLandBoundaryIndex"></param>
        /// <param name="leftIndex"></param>
        /// <param name="rightIndex"></param>
        /// <param name="leftEdgeRatio"></param>
        /// <param name="rightEdgeRatio"></param>
        /// <param name="startMeshNode"></param>
        /// <param name="endMeshNode"></param>
        /// <returns></returns>
        bool FindStartEndMeshNodes(int startLandBoundaryIndex,
                                   int endLandBoundaryIndex,
                                   int leftIndex,
                                   int rightIndex,
                                   double leftEdgeRatio,
                                   double rightEdgeRatio,
                                   int& startMeshNode,
                                   int& endMeshNode);

        /// <summary>
        /// Connect mesh nodes starting from startMeshNode, using Dijkstra's shortest path algorithm. 
        /// The distance of each edge is the edge length multiplied by the distance from the land boundary
        /// </summary>
        /// <param name="mesh"></param>
        /// <param name="polygons"></param>
        /// <param name="landBoundarySegment"></param>
        /// <param name="startLandBoundaryIndex"></param>
        /// <param name="endLandBoundaryIndex"></param>
        /// <param name="startMeshNode"></param>
        /// <param name="meshBoundOnly"></param>
        /// <param name="connectedNodes"></param>
        /// <returns></returns>
        bool ShortestPath(int landBoundarySegment,
                          int startLandBoundaryIndex,
                          int endLandBoundaryIndex,
                          int startMeshNode,
                          bool meshBoundOnly,
                          std::vector<int>& connectedNodes);

        /// <summary>
        /// Compute the nearest point on the land boundary (toland)
        /// </summary>
        /// <param name="projection"></param>
        /// <param name="node"></param>
        /// <param name="startLandBoundaryIndex"></param>
        /// <param name="endLandBoundaryIndex"></param>
        /// <param name="minimumDistance"></param>
        /// <param name="pointOnLandBoundary"></param>
        /// <param name="nearestLandBoundaryNodeIndex"></param>
        /// <param name="edgeRatio"></param>
        /// <returns></returns>
        bool NearestLandBoundaryNode(const Projections& projection,
                                     const Point& node,
                                     int startLandBoundaryIndex,
                                     int endLandBoundaryIndex,
                                     double& minimumDistance,
                                     Point& pointOnLandBoundary,
                                     int& nearestLandBoundaryNodeIndex,
                                     double& edgeRatio);

        /// <summary>
        /// TODO: It could be moved to generic operations (cellcrossedbyland)
        /// </summary>
        /// <param name="mesh"></param>
        /// <param name="face"></param>
        /// <param name="startLandBoundaryIndex"></param>
        /// <param name="endLandBoundaryIndex"></param>
        /// <returns></returns>
        bool IsFaceCrossedByLandBoundaries(int face, 
                                           int startLandBoundaryIndex, 
                                           int endLandBoundaryIndex);

        std::shared_ptr<Mesh>     m_mesh;                 // A pointer to mesh 
        std::shared_ptr<Polygons> m_polygons;             // A pointer to polygons

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
        double m_closeFactor = 5.0;
    };

}