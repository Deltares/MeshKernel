//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2021.
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
    class Mesh2D;

    /// @brief A class describing land boundaries.
    /// These are used to visualise the land-water interface.
    ///
    /// The main responsibility of this class is to store the land boundary polygons,
    /// categorize them based on their proximity to a mesh
    /// and provide the functionality to assign each mesh node to the appropriate land boundary polyline.
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
        /// @param[in] landBoundary
        /// @param[in] mesh
        /// @param[in] polygons
        /// @returns
        LandBoundaries(const std::vector<Point>& landBoundary,
                       std::shared_ptr<Mesh2D> mesh,
                       std::shared_ptr<Polygons> polygons);

        /// @brief The land boundary will be split into segments that are within the polygon,
        /// and either close or not to the mesh boundary (admin_landboundary_segments)
        ///
        /// This method uses a Point vector member variable and identifies
        /// the start-end points of each land boundary polyline with the requirement
        /// that all polyline nodes are close enough to the mesh boundary and is inside the polygon.
        /// \image html LandBoundarySegmentation_step1.jpg  "Land boundary segmentation"
        void Administrate();

        /// @brief Find the mesh boundary line closest to the land boundary (find_nearest_meshline)
        /// @param[in] projectToLandBoundaryOption the option to use to project to land boundary
        void FindNearestMeshBoundary(ProjectToLandBoundaryOption projectToLandBoundaryOption);

        /// @brief Snap mesh nodes to land boundaries (snap_to_landboundary)
        void SnapMeshToLandBoundaries();

        /// @brief Gets the number of nodes
        /// @return the number of nodes
        auto GetNumNodes() const { return m_nodes.size(); }

        std::vector<size_t> m_meshNodesLandBoundarySegments; ///< lanseg_map, mesh nodes to land boundary mapping

    private:
        /// @brief Build an additional boundary for not assigned nodes (connect_boundary_paths)
        /// @param[in] edgeIndex
        /// @param[in] initialize
        /// @param[in] nodes
        /// @param[in] numNodes
        void AssignSegmentsToMeshNodes(size_t edgeIndex,
                                       bool initialize,
                                       std::vector<size_t>& nodes,
                                       size_t numNodes);

        /// @brief Add new land boundary segment that connects two others (add_land)
        /// @param[in] nodesLoc
        /// @param[in] numNodesLoc
        /// @param[in] nodeIndex
        void AddLandBoundary(const std::vector<size_t>& nodesLoc,
                             size_t numNodesLoc,
                             size_t nodeIndex);

        /// @brief Assigns to each mesh node a land boundary segment index ()
        /// @param[in] landBoundarySegment
        /// @param[in] meshBoundOnly
        /// @param[out] numNodesInPath
        /// @param[out] numRejectedNodesInPath
        void MakePath(size_t landBoundarySegment,
                      bool meshBoundOnly,
                      size_t& numNodesInPath,
                      size_t& numRejectedNodesInPath);

        /// @brief Mask the mesh nodes to be considered in the shortest path algorithm for the current segmentIndex.
        /// It is setting leftIndex, rightIndex, leftEdgeRatio, rightEdgeRatio (masknodes).
        /// @param[in] segmentIndex
        /// @param[in] meshBoundOnly
        /// @param[in] startLandBoundaryIndex
        /// @param[in] endLandBoundaryIndex
        /// @param[out] leftIndex
        /// @param[out] rightIndex
        /// @param[out] leftEdgeRatio
        /// @param[out] rightEdgeRatio
        //// \image html LandBoundaryNodeFlagging_step2.jpg  "Flag the mesh node close to the land boundary"
        void ComputeMask(size_t segmentIndex,
                         bool meshBoundOnly,
                         size_t startLandBoundaryIndex,
                         size_t endLandBoundaryIndex,
                         size_t& leftIndex,
                         size_t& rightIndex,
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
                       std::vector<size_t>& landBoundaryFaces,
                       size_t startNodeLandBoundaryIndex,
                       size_t endNodeLandBoundaryindex,
                       size_t& leftIndex,
                       size_t& rightIndex,
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
        [[nodiscard]] bool IsMeshEdgeCloseToLandBoundaries(size_t edgeIndex,
                                                           size_t startNodeLandBoundaryIndex,
                                                           size_t endNodeLandBoundaryIndex,
                                                           bool meshBoundOnly,
                                                           size_t& leftIndex,
                                                           size_t& rightIndex,
                                                           double& leftEdgeRatio,
                                                           double& rightEdgeRatio,
                                                           size_t& landBoundaryNode);

        /// @brief Finds the start and end mesh node.
        /// These are the nodes that are on a edge close to the land boundary segment (get_kstartend2)
        /// @param[in] endLandBoundaryIndex
        /// @param[in] leftIndex
        /// @param[in] rightIndex
        /// @param[in] leftEdgeRatio
        /// @param[in] rightEdgeRatio
        /// @param[out] startMeshNode
        /// @param[out] endMeshNode
        //// \image html LandBoundaryDijkstra_step4.jpg  "Compute the land boundary representation on the mesh using the Djikstra shortest path algorithm."
        void FindStartEndMeshNodes(size_t endLandBoundaryIndex,
                                   size_t leftIndex,
                                   size_t rightIndex,
                                   double leftEdgeRatio,
                                   double rightEdgeRatio,
                                   size_t& startMeshNode,
                                   size_t& endMeshNode);

        /// @brief Finds the start and end mesh node from given edges.
        /// @param[in] startEdge
        /// @param[in] endEdge
        /// @param[in] startPoint
        /// @param[in] endPoint
        /// @param[out] startMeshNode
        /// @param[out] endMeshNode
        /// \image html LandBoundaryStartEndNodes_step3.jpg  "Find the start and end mesh nodes of the land boundary on the mesh."
        void FindStartEndMeshNodesFromEdges(size_t startEdge,
                                            size_t endEdge,
                                            Point startPoint,
                                            Point endPoint,
                                            size_t& startMeshNode,
                                            size_t& endMeshNode) const;

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
        void ShortestPath(size_t landBoundarySegment,
                          size_t startLandBoundaryIndex,
                          size_t endLandBoundaryIndex,
                          size_t startMeshNode,
                          bool meshBoundOnly,
                          std::vector<size_t>& connectedNodes);

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
                                     size_t startLandBoundaryIndex,
                                     size_t endLandBoundaryIndex,
                                     double& minimumDistance,
                                     Point& pointOnLandBoundary,
                                     size_t& nearestLandBoundaryNodeIndex,
                                     double& edgeRatio);

        /// @brief (cellcrossedbyland)
        /// @param face
        /// @param startLandBoundaryIndex
        /// @param endLandBoundaryIndex
        [[nodiscard]] bool IsFaceCrossedByLandBoundaries(size_t face,
                                                         size_t startLandBoundaryIndex,
                                                         size_t endLandBoundaryIndex);

        std::shared_ptr<Mesh2D> m_mesh;                    ///< A pointer to mesh
        std::shared_ptr<Polygons> m_polygons;              ///< A pointer to polygons
        std::vector<Point> m_nodes;                        ///< XLAN, YLAN, ZLAN
        std::vector<Point> m_polygonNodesCache;            ///< array of points (e.g. points of a face)
        std::vector<std::vector<size_t>> m_segmentIndices; ///< lanseg_startend
        std::vector<std::vector<double>> m_nodesLand;      ///< node to land boundary segment mapping

        std::vector<size_t> m_nodeMask; ///< nodemask, masking the net nodes
        std::vector<size_t> m_faceMask; ///< masking faces
        std::vector<size_t> m_edgeMask; ///< masking edges

        bool m_landMask = true;          ///< Whether land masks were given
        bool m_addLandboundaries = true; ///< Whether to add land boundaries
        size_t m_numFacesMasked = 0;     ///< Number of masked faces
        size_t m_maskDepth = 0;          ///< Mask depth

        // caches
        std::vector<double> m_nodesMinDistances; ///< Min distances of nodes

        // Parameters
        const double m_closeToLandBoundaryFactor = 5.0; ///< close - to - landboundary tolerance, measured in number of meshwidths
        const double m_closeWholeMeshFactor = 1.0;      ///< close - to - landboundary tolerance, measured in number of meshwidths
        const double m_minDistanceFromLandFactor = 2.0; ///< Minimum distance from land factor
        double m_closeFactor = 5.0;                     ///< Factor to determine minimal distance from mesh nodes
    };

} // namespace meshkernel
