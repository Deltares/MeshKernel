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

        /// @brief Default constructor
        /// @brief landBoundary The points describing the landboundaries
        /// @brief mesh The 2d mesh
        /// @brief polygons Account for land boundaries within the polygon
        LandBoundaries(const std::vector<Point>& landBoundary,
                       std::shared_ptr<Mesh2D> mesh,
                       std::shared_ptr<Polygons> polygons);

        /// @brief The land boundary segments close enough to the mesh boundary are flagged (admin_landboundary_segments)
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
        void AssignLandBoundaryPolylineToMeshNodes(size_t edgeIndex,
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

        /// @brief Assigns to each mesh node a landBoundaryPolyline
        /// @param[in] landBoundaryPolyline The current landboundary segment
        /// @param[in] meshBoundOnly Account only for mesh boundaries
        /// @param[out] numNodesInPath The number of mesh nodes for this path
        /// @param[out] numRejectedNodesInPath The number of rejected nodes in path
        void MakePath(size_t landBoundaryPolyline,
                      size_t& numNodesInPath,
                      size_t& numRejectedNodesInPath);

        /// @brief Mask the mesh nodes to be considered in the shortest path algorithm for the current landBoundaryPolyline.
        /// Is setting leftIndex, rightIndex, leftEdgeRatio, rightEdgeRatio (masknodes).
        /// @param[in] landBoundaryPolyline
        /// @param[in] meshBoundOnly
        /// @param[out] leftIndex
        /// @param[out] rightIndex
        /// @param[out] leftEdgeRatio
        /// @param[out] rightEdgeRatio
        void ComputeMeshNodeMask(size_t landBoundaryPolyline);

        /// @brief Mask all face close to a land boundary, starting from a seed of others and growing from there (maskcells)
        /// @param[in] meshBoundOnly Account only for mesh boundary
        /// @param[in] initialFaces The initial face seed
        /// @param[in] startNodeLandBoundaryIndex The start index of the current land boundary node to account for
        /// @param[in] endNodeLandBoundaryindex The end index of the land boundary node to account for
        /// @param[out] leftIndex
        /// @param[out] rightIndex
        /// @param[out] leftEdgeRatio
        /// @param[out] rightEdgeRatio
        void MaskMeshFaceMask(std::vector<size_t>& initialFaces,
                              size_t landBoundaryPolyline);

        /// @brief Check if a mesh edge is close to a land boundary segment (linkcrossedbyland)
        /// @param[in] meshEdgeIndex
        /// @param[in] startNodeLandBoundaryIndex
        /// @param[in] endNodeLandBoundaryIndex
        /// @param[in] meshBoundOnly
        /// @param[out] leftIndex
        /// @param[out] rightIndex
        /// @param[out] leftEdgeRatio
        /// @param[out] rightEdgeRatio
        /// @return the closest land boundary node
        [[nodiscard]] size_t IsMeshEdgeCloseToLandBoundaries(size_t meshEdgeIndex,
                                                             size_t landBoundaryPolyline);

        /// @brief Finds the start and end mesh node.
        /// These are the nodes that are on a edge close to the land boundary segment (get_kstartend2)
        /// @param[in] landBoundaryPolyline
        /// @returns the start and end index on a mesh
        std::tuple<size_t, size_t> FindStartEndMeshNodes(size_t landBoundaryPolyline);

        /// @brief Finds the start and end mesh node from given edges.
        /// @param[in] startEdge
        /// @param[in] endEdge
        /// @param[in] startPoint
        /// @param[in] endPoint
        /// @param[out] startMeshNode
        /// @param[out] endMeshNode
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
        /// @brief landBoundaryPolyline
        /// @brief startLandBoundaryIndex
        /// @brief endLandBoundaryIndex
        /// @brief startMeshNode
        /// @brief meshBoundOnly
        /// @brief connectedNodes
        /// @returns
        std::vector<size_t> ShortestPath(size_t landBoundaryPolyline, size_t startMeshNode);

        /// @brief Compute the nearest land boundary segment (toland)
        /// @param projection
        /// @param node
        /// @param startLandBoundaryIndex
        /// @param endLandBoundaryIndex
        /// @param minimumDistance
        /// @param pointOnLandBoundary
        /// @param nearestLandBoundaryNodeIndex
        /// @param edgeRatio
        void NearestLandBoundarySegment(const Point& node,
                                        size_t landBoundaryPolyline,
                                        double& minimumDistance,
                                        Point& pointOnLandBoundary,
                                        size_t& nearestLandBoundaryNodeIndex,
                                        double& edgeRatio);

        std::shared_ptr<Mesh2D> m_mesh;                         // A pointer to mesh
        std::shared_ptr<Polygons> m_polygons;                   // A pointer to polygons
        std::vector<Point> m_nodes;                             // XLAN, YLAN, ZLAN
        std::vector<Point> m_polygonNodesCache;                 // array of points (e.g. points of a face)
        std::vector<std::vector<size_t>> m_validLandBoundaries; // lanseg_startend
        std::vector<std::vector<double>> m_nodesLand;           // node to land boundary segment mapping
        std::vector<size_t> m_nodeFaceIndices;                  // For each node, the indices of the faces including them

        std::vector<size_t> m_nodeMask; // nodemask, masking the net nodes
        std::vector<bool> m_faceMask;   // masking faces
        std::vector<size_t> m_edgeMask; // masking edges

        bool m_landMask = true;
        bool m_addLandboundaries = true;

        // caches
        std::vector<double> m_nodesMinDistances;

        // Parameters
        const double m_closeToLandBoundaryFactor = 5.0; // close - to - landboundary tolerance, measured in number of meshwidths
        const double m_closeWholeMeshFactor = 1.0;      // close - to - landboundary tolerance, measured in number of meshwidths
        const double m_minDistanceFromLandFactor = 2.0;
        double m_closeFactor = 5.0;

        //findOnlyOuterMeshBoundary
        bool m_findOnlyOuterMeshBoundary = false;
    };

} // namespace meshkernel
