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

        /// @brief Default constructor
        LandBoundaries() = default;

        /// @brief Default constructor
        /// @param[in] landBoundary A vector of points defining the land boundary
        /// @param[in] mesh The current 2d mesh
        /// @param[in] polygons A polygon for selecting part of the land boundaries
        LandBoundaries(const std::vector<Point>& landBoundary,
                       std::shared_ptr<Mesh2D> mesh,
                       std::shared_ptr<Polygons> polygons);

        /// @brief The portion of the boundary segments close enough to the mesh boundary are flagged (admin_landboundary_segments)
        ///
        /// This method uses a Point vector member variable and identifies
        /// the start-end points of each land boundary polyline with the requirement
        /// that all polyline nodes are close enough to the mesh boundary and is inside the polygon.
        /// \image html LandBoundarySegmentation_step1.jpg  "Land boundary segmentation"
        void Administrate();

        /// @brief Find the mesh boundary line closest to the land boundary (find_nearest_meshline)
        /// @param[in] projectToLandBoundaryOption the option to use to project to land boundary
        void FindNearestMeshBoundary(ProjectToLandBoundaryOption projectToLandBoundaryOption);

        /// @brief Snap the mesh nodes to land boundaries (snap_to_landboundary)
        void SnapMeshToLandBoundaries();

        /// @brief Gets the number of land boundary nodes
        /// @return the number of land boundary nodes
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

        /// @brief Assigns to each mesh node a landboundaryIndex
        /// @param[in] landboundaryIndex The current landboundary segment
        /// @param[out] numNodesInPath The number of mesh nodes for this path
        /// @param[out] numRejectedNodesInPath The number of rejected nodes in path
        void MakePath(size_t landboundaryIndex,
                      size_t& numNodesInPath,
                      size_t& numRejectedNodesInPath);

        /// @brief Mask the mesh nodes to be considered in the shortest path algorithm for the current land boundary polyline (masknodes).
        /// @param[in] landboundaryIndex The land boundary polyline  index
        void ComputeMeshNodeMask(size_t landboundaryIndex);

        /// @brief Mask the mesh nodes to be considered in the shortest path algorithm for the current segmentIndex.
        /// It is setting leftIndex, rightIndex, leftEdgeRatio, rightEdgeRatio (masknodes).
        //// \image html LandBoundaryNodeFlagging_step2.jpg  "Flag the mesh node close to the land boundary"
        /// @param[in] segmentIndex
        /// @param[in] meshBoundOnly
        /// @param[in] startLandBoundaryIndex
        /// @param[in] endLandBoundaryIndex
        /// @param[out] leftIndex
        /// @param[out] rightIndex
        /// @param[out] leftEdgeRatio
        /// @param[out] rightEdgeRatio
        void ComputeMask(size_t segmentIndex,
                         bool meshBoundOnly,
                         size_t startLandBoundaryIndex,
                         size_t endLandBoundaryIndex,
                         size_t& leftIndex,
                         size_t& rightIndex,
                         double& leftEdgeRatio,
                         double& rightEdgeRatio);

        /// @brief Mask all face close to a land boundary, starting from a seed of others and growing from there (maskcells)
        /// @param[in] landboundaryIndex The land boundary polyline index
        /// @param[in] initialFaces The initial face seeds
        void MaskMeshFaceMask(size_t landboundaryIndex, std::vector<size_t>& initialFaces);

        /// @brief Check if a mesh edge is close to a land boundary segment (linkcrossedbyland)
        /// @param[in] landboundaryIndex The land boundary polyline index
        /// @param[in] edge The mesh edge to inquire
        /// @return the closest land boundary node
        [[nodiscard]] size_t IsMeshEdgeCloseToLandBoundaries(size_t landboundaryIndex, size_t edge);

        /// @brief Finds the start and the end mesh node indices which correspond to a landboundary polyline.
        /// These are the nodes that are on a edge close to the land boundary segment (get_kstartend2)
        ///
        /// \image html LandBoundaryDijkstra_step4.jpg  "Compute the land boundary representation on the mesh using the Djikstra shortest path algorithm."
        /// @param[in] landboundaryIndex The land boundary polyline index
        /// @returns the start and the end mesh nodes indices
        std::tuple<size_t, size_t> FindStartEndMeshNodesDijkstraAlgorithm(size_t landboundaryIndex);

        /// @brief Finds the edge nodes closest to a point
        ///
        /// \image html LandBoundaryStartEndNodes_step3.jpg  "Find the start and end mesh nodes of the land boundary on the mesh."
        /// @param[in] edge The edge index
        /// @param[in] point The point to inquire
        size_t FindStartEndMeshNodesFromEdges(size_t edge, Point point) const;

        /// @brief Connect mesh nodes close to the landboundaryIndex using Dijkstra's algorithm
        /// @param[in] landboundaryIndex The index of a valid landboundary
        /// @param[in] startMeshNode the starting point
        /// @returns A vector of connected edge indices for each node
        std::vector<size_t> ShortestPath(size_t landboundaryIndex, size_t startMeshNode);

        /// @brief Compute the nearest land boundary segment (toland)
        /// @param[in] landBoundaryIndex The land boundary index
        /// @param[in] node The node
        /// @returns A tuple containing the distance of the node from the land boundary, the projected node on the land boundary, the closest land boundary node,
        /// The length of the segment from the starting point to the projected point expressed as an edge ratio
        std::tuple<double, Point, size_t, double> NearestLandBoundarySegment(int landBoundaryIndex, const Point& node);

        std::shared_ptr<Mesh2D> m_mesh;                         ///< A pointer to mesh
        std::shared_ptr<Polygons> m_polygons;                   ///< A pointer to polygons
        std::vector<Point> m_nodes;                             ///< XLAN, YLAN, ZLAN
        std::vector<Point> m_polygonNodesCache;                 ///< array of points (e.g. points of a face)
        std::vector<std::vector<size_t>> m_validLandBoundaries; ///< lanseg_startend
        std::vector<std::vector<double>> m_nodesLand;           ///< node to land boundary segment mapping
        std::vector<size_t> m_nodeFaceIndices;                  ///< For each node, the indices of the faces including them

        std::vector<size_t> m_nodeMask; ///< nodemask, masking the net nodes
        std::vector<bool> m_faceMask;   ///< masking faces
        std::vector<size_t> m_edgeMask; ///< masking edges

        bool m_landMask = true;          ///< Land mask
        bool m_addLandboundaries = true; ///< Whether to add land boundaries

        // caches
        std::vector<double> m_nodesMinDistances; ///< the minimum distances to land boundaries

        // Parameters
        const double m_closeToLandBoundaryFactor = 5.0; ///< close - to - landboundary tolerance, measured in number of meshwidths
        const double m_closeWholeMeshFactor = 1.0;      ///< close - to - landboundary tolerance, measured in number of meshwidths
        const double m_minDistanceFromLandFactor = 2.0; ///< Minimal distance from land factor
        double m_closeFactor = 5.0;                     ///< Close factor

        //findOnlyOuterMeshBoundary
        bool m_findOnlyOuterMeshBoundary = false; ///< Whether to find only outer mesh boundary
    };

} // namespace meshkernel
