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

#include <memory>
#include <vector>

#include <MeshKernel/Entities.hpp>

namespace meshkernel
{
    class Polygons;
    class Mesh;

    /// @brief A class describing land boundaries, which are used to visualise the land-water interface
    class LandBoundaries
    {

    public:
        /// Enumerator describing the options how to project to the land boundary
        enum class ProjectToLandBoundaryOption
        {
            DoNotProjectToLandBoundary = 0,
            ToOriginalNetBoundary = 1,
            OuterMeshBoundaryToLandBoundary = 2,
            InnerAndOuterMeshBoundaryToLandBoundary = 3,
            WholeMesh = 4
        };

        /// @brief Default Ctor
        /// @brief landBoundary
        /// @brief mesh
        /// @brief polygons
        /// @returns
        LandBoundaries(const std::vector<Point>& landBoundary,
                       std::shared_ptr<Mesh> mesh,
                       std::shared_ptr<Polygons> polygons);

        /// @brief The land boundary will be split into segments that are within the polygon, and either close or not to the mesh boundary (admin_landboundary_segments)
        void Administrate();

        /// @brief Find the mesh boundary line closest to the land boundary (find_nearest_meshline)
        /// @param[in] projectToLandBoundaryOption the option to use to project to land boundary
        void FindNearestMeshBoundary(ProjectToLandBoundaryOption projectToLandBoundaryOption);

        /// @brief Snap mesh nodes to land boundaries (snap_to_landboundary)
        void SnapMeshToLandBoundaries();

        /// @brief Gets the number of nodes
        /// @return the number of nodes
        int GetNumNodes() const { return static_cast<int>(m_nodes.size()); };

        std::vector<int> m_meshNodesLandBoundarySegments; ///< lanseg_map, mesh nodes to land boundary mapping

    private:
        /// @brief Build an additional boundary for not assigned nodes (connect_boundary_paths)
        /// @param[in] edgeIndex
        /// @param[in] initialize
        /// @param[in] nodes
        /// @param[in] numNodes
        void AssignSegmentsToMeshNodes(int edgeIndex,
                                       bool initialize,
                                       std::vector<int>& nodes,
                                       int numNodes);

        /// @brief Add new land boundary segment that connects two others (add_land)
        /// @param[in] nodesLoc
        /// @param[in] numNodesLoc
        /// @param[in] nodeIndex
        void AddLandBoundary(const std::vector<int>& nodesLoc,
                             int numNodesLoc,
                             int nodeIndex);

        /// @brief Assigns to each mesh node a land boundary segment index ()
        /// @param[in] landBoundarySegment
        /// @param[in] meshBoundOnly
        /// @param[out] numNodesInPath
        /// @param[out] numRejectedNodesInPath
        void MakePath(int landBoundarySegment,
                      bool meshBoundOnly,
                      int& numNodesInPath,
                      int& numRejectedNodesInPath);

        /// @brief Mask the mesh nodes to be considered in the shortest path algorithm for the current segmentIndex.
        /// Is setting leftIndex, rightIndex, leftEdgeRatio, rightEdgeRatio (masknodes).
        /// @param[in] segmentIndex
        /// @param[in] meshBoundOnly
        /// @param[in] startLandBoundaryIndex
        /// @param[in] endLandBoundaryIndex
        /// @param[out] leftIndex
        /// @param[out] rightIndex
        /// @param[out] leftEdgeRatio
        /// @param[out] rightEdgeRatio
        void ComputeMask(int segmentIndex,
                         bool meshBoundOnly,
                         int startLandBoundaryIndex,
                         int endLandBoundaryIndex,
                         int& leftIndex,
                         int& rightIndex,
                         double& leftEdgeRatio,
                         double& rightEdgeRatio);

        /// @brief Mask the faces that are intersected by the land boundary (maskcells)
        /// @param[in] meshBoundOnly
        /// @param[in] landBoundaryFaces
        /// @param[in] startNodeLandBoundaryIndex
        /// @param[in] endNodeLandBoundaryindex
        /// @param[out] leftIndex
        /// @param[out] rightIndex
        /// @param[out] leftEdgeRatio
        /// @param[out] rightEdgeRatio
        void MaskFaces(bool meshBoundOnly,
                       std::vector<int>& landBoundaryFaces,
                       int startNodeLandBoundaryIndex,
                       int endNodeLandBoundaryindex,
                       int& leftIndex,
                       int& rightIndex,
                       double& leftEdgeRatio,
                       double& rightEdgeRatio);

        /// @brief Check if a mesh edge is close to a land boundary segment (linkcrossedbyland)
        /// @param[in] edgeIndex
        /// @param[in] startNodeLandBoundaryIndex
        /// @param[in] endNodeLandBoundaryIndex
        /// @param[in] meshBoundOnly
        /// @param[out] leftIndex
        /// @param[out] rightIndex
        /// @param[out] leftEdgeRatio
        /// @param[out] rightEdgeRatio
        /// @param[out] landBoundaryNode
        [[nodiscard]] bool IsMeshEdgeCloseToLandBoundaries(int edgeIndex,
                                                           int startNodeLandBoundaryIndex,
                                                           int endNodeLandBoundaryIndex,
                                                           bool meshBoundOnly,
                                                           int& leftIndex,
                                                           int& rightIndex,
                                                           double& leftEdgeRatio,
                                                           double& rightEdgeRatio,
                                                           int& landBoundaryNode);

        /// @brief Finds the start and end mesh node.
        /// These are the nodes that are on a edge close to the land boundary segment (get_kstartend2)
        /// @param[in] endLandBoundaryIndex
        /// @param[in] leftIndex
        /// @param[in] rightIndex
        /// @param[in] leftEdgeRatio
        /// @param[in] rightEdgeRatio
        /// @param[out] startMeshNode
        /// @param[out] endMeshNode
        void FindStartEndMeshNodes(int endLandBoundaryIndex,
                                   int leftIndex,
                                   int rightIndex,
                                   double leftEdgeRatio,
                                   double rightEdgeRatio,
                                   int& startMeshNode,
                                   int& endMeshNode);

        /// @brief Finds the start and end mesh node from given edges.
        /// @param[in] startEdge
        /// @param[in] endEdge
        /// @param[in] startPoint
        /// @param[in] endPoint
        /// @param[out] startMeshNode
        /// @param[out] endMeshNode
        void FindStartEndMeshNodesFromEdges(int startEdge,
                                            int endEdge,
                                            Point startPoint,
                                            Point endPoint,
                                            int& startMeshNode,
                                            int& endMeshNode) const;

        /// @brief Connect mesh nodes starting from startMeshNode, using Dijkstra's shortest path algorithm.
        /// The distance of each edge is the edge length multiplied by the distance from the land boundary
        /// @brief mesh
        /// @brief polygons
        /// @brief landBoundarySegment
        /// @brief startLandBoundaryIndex
        /// @brief endLandBoundaryIndex
        /// @brief startMeshNode
        /// @brief meshBoundOnly
        /// @brief connectedNodes
        /// @returns
        void ShortestPath(int landBoundarySegment,
                          int startLandBoundaryIndex,
                          int endLandBoundaryIndex,
                          int startMeshNode,
                          bool meshBoundOnly,
                          std::vector<int>& connectedNodes);

        /// @brief Compute the nearest node on the land boundary (toland)
        /// @param projection
        /// @param node
        /// @param startLandBoundaryIndex
        /// @param endLandBoundaryIndex
        /// @param minimumDistance
        /// @param pointOnLandBoundary
        /// @param nearestLandBoundaryNodeIndex
        /// @param edgeRatio
        void NearestLandBoundaryNode(const Projection& projection,
                                     const Point& node,
                                     int startLandBoundaryIndex,
                                     int endLandBoundaryIndex,
                                     double& minimumDistance,
                                     Point& pointOnLandBoundary,
                                     int& nearestLandBoundaryNodeIndex,
                                     double& edgeRatio);

        /// @brief (cellcrossedbyland)
        /// @param face
        /// @param startLandBoundaryIndex
        /// @param endLandBoundaryIndex
        [[nodiscard]] bool IsFaceCrossedByLandBoundaries(size_t face,
                                                         size_t startLandBoundaryIndex,
                                                         size_t endLandBoundaryIndex);

        std::shared_ptr<Mesh> m_mesh;                      // A pointer to mesh
        std::shared_ptr<Polygons> m_polygons;              // A pointer to polygons
        std::vector<Point> m_nodes;                        // XLAN, YLAN, ZLAN
        std::vector<Point> m_polygonNodesCache;            // array of points (e.g. points of a face)
        std::vector<std::vector<size_t>> m_segmentIndices; // lanseg_startend
        std::vector<std::vector<double>> m_nodesLand;      // node to land boundary segment mapping

        std::vector<int> m_nodeMask; // nodemask, masking the net nodes
        std::vector<int> m_faceMask; // masking faces
        std::vector<int> m_edgeMask; // masking edges

        bool m_landMask = true;
        bool m_addLandboundaries = true;
        int m_numFacesMasked = 0;
        int m_maskDepth = 0;

        // caches
        std::vector<double> m_nodesMinDistances;
        const size_t m_allocationSize = 10000; // allocation size for allocateVector

        // Parameters
        const double m_closeToLandBoundaryFactor = 5.0; // close - to - landboundary tolerance, measured in number of meshwidths
        const double m_closeWholeMeshFactor = 1.0;      // close - to - landboundary tolerance, measured in number of meshwidths
        const double m_minDistanceFromLandFactor = 2.0;
        double m_closeFactor = 5.0;
    };

} // namespace meshkernel
