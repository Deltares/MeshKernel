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

#include <MeshKernel/AveragingInterpolation.hpp>
#include <MeshKernel/Parameters.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/Utilities/RTreeBase.hpp>

namespace meshkernel
{
    // Forward declarations
    class Mesh2D;

    /// @brief A class used to refine a Mesh2D instance
    ///
    /// Mesh refinement operates on Mesh2D and is based on
    /// iteratively splitting the edges until the desired level of refinement
    /// or the maximum number of iterations is reached.
    /// Refinement can be based on samples or based on a polygon.
    /// The refinement based on samples uses
    /// the averaging interpolation algorithm to compute the level of refinement
    /// from the samples to the centers of the edges.
    /// At a high level, the mesh refinement is performed as follow:
    ///
    /// -   Flag the nodes inside the refinement polygon.
    ///
    /// -   Flag all face nodes of the faces not fully included in the polygon.
    ///
    /// -   Execute the refinement iterations
    ///
    ///     1.  For each edge store the index of its neighboring edge sharing a
    ///         hanging node (the so-called brother edge). This is required for
    ///         the following steps because edges with hanging nodes will not be
    ///         divided further.
    ///
    ///     2.  Compute edge and face refinement masks from the samples.
    ///
    ///     3.  Compute if a face should be divided based on the computed
    ///         refinement value.
    ///
    ///     4.  Split the face by dividing the edges.
    ///
    /// -   Connect the hanging nodes if required, thus forming triangular faces
    ///     in the transition area.
    ///
    /// As with OrthogonalizationAndSmoothing, MeshRefinement modifies an
    /// existing Mesh2D instance.
    class MeshRefinement
    {
        /// @brief Enumerator describing the face location types
        enum class FaceLocation
        {
            Land = 1,
            Water = 2,
            LandWater = 3
        };

        /// @brief Enumerator describing the different refinement types
        enum class RefinementType
        {
            WaveCourant = 1,
            RefinementLevels = 2,
            RidgeDetection = 3
        };

    public:
        /// @brief The constructor for refining based on samples
        /// @param[in] mesh The mesh to be refined
        /// @param[in] interpolant The averaging interpolation to use
        /// @param[in] meshRefinementParameters The mesh refinement parameters
        MeshRefinement(std::shared_ptr<Mesh2D> mesh,
                       std::shared_ptr<MeshInterpolation> interpolant,
                       const MeshRefinementParameters& meshRefinementParameters);

        /// @brief The constructor for refining based on samples
        /// @param[in] mesh The mesh to be refined
        /// @param[in] interpolant The averaging interpolation to use
        /// @param[in] meshRefinementParameters The mesh refinement parameters
        /// @param[in] useNodalRefinement Use nodal refinement
        MeshRefinement(std::shared_ptr<Mesh2D> mesh,
                       std::shared_ptr<MeshInterpolation> interpolant,
                       const MeshRefinementParameters& meshRefinementParameters,
                       bool useNodalRefinement);

        /// @brief The constructor for refining based on polygons
        /// @param[in] mesh The mesh to be refined
        /// @param[in] polygon The polygon where to refine
        /// @param[in] meshRefinementParameters The mesh refinement parameters
        MeshRefinement(std::shared_ptr<Mesh2D> mesh,
                       const Polygons& polygon,
                       const MeshRefinementParameters& meshRefinementParameters);

        /// @brief Compute mesh refinement (refinecellsandfaces2).
        ///
        /// Steps:
        /// 1. Masks the node to be refined (those inside a polygon)
        /// 2. Find the brother edges, the edge sharing a hanging node, FindBrotherEdges
        /// 3. Mask nodes at the polygon perimeter, ComputeNodeMaskAtPolygonPerimeter
        /// 4. Do refinement iterations
        ///    -# Find the brother edges, FindBrotherEdges
        ///    -# Compute the edge refinement mask based on samples, ComputeRefinementMasksFromSamples
        ///    -# Compute the edge refinement mask based on polygon, ComputeEdgesRefinementMask
        ///    -# Compute if a face should be split, ComputeIfFaceShouldBeSplit
        ///    -# Compute face by splitting edges, RefineFacesBySplittingEdges
        /// 5. Connect hanging nodes if requested, DeleteIsolatedHangingnodes, connect_hanging_nodes
        void Compute();

    private:
        /// @brief Finds if two edges are brothers, sharing an hanging node. Can be moved to Mesh2D
        void FindBrotherEdges();

        /// @brief Modifies m_nodeMask, all nodes of the faces intersecting the polygon perimeter will get value of -2 (set_initial_mask)
        ///        The mask value of the other nodes will not be modified.
        void ComputeNodeMaskAtPolygonPerimeter();

        /// @brief Computes the edge and face refinement mask from sample values (compute_jarefine_poly)
        void ComputeRefinementMasksFromSamples();

        /// @brief Computes the edge and face refinement mask from the minimum edge size
        void ComputeRefinementMaskFromEdgeSize();

        /// @brief Computes the refinement mask for refinement based on levels
        void ComputeRefinementMasksForRefinementLevels(UInt face,
                                                       size_t& numberOfEdgesToRefine,
                                                       std::vector<UInt>& edgeToRefine) const;

        /// @brief Computes the refinement mask for refinement based on wave courant criteria
        void ComputeRefinementMasksForWaveCourant(UInt face,
                                                  size_t& numberOfEdgesToRefine,
                                                  std::vector<UInt>& edgeToRefine);

        /// @brief Computes the refinement mask for refinement based on ridge detection
        void ComputeRefinementMasksForRidgeDetection(UInt face,
                                                     size_t& numberOfEdgesToRefine,
                                                     std::vector<UInt>& edgeToRefine) const;

        /// @brief Computes refinement masks (compute_jarefine_poly)
        ///        Face nodes, edge and edge lengths are stored in local caches. See Mesh2D.FaceClosedPolygon method
        /// @param face The face index
        void ComputeRefinementMasksFromSamples(UInt face);

        /// Computes the edge refinement mask (comp_jalink)
        void ComputeEdgesRefinementMask();

        /// @brief Finds the hanging nodes in a face (find_hangingnodes)
        /// @param[in] face The current face index
        void FindHangingNodes(UInt face);

        /// @brief Get the number of hanging nodes
        /// @returns The number of hanging nodes
        UInt CountHangingNodes() const;

        /// @brief Get the number of hanging nodes
        /// @returns The number of hanging nodes
        UInt CountHangingEdges() const;

        /// @brief Get the number of hanging nodes
        /// @param[in] face The current face index
        /// @returns The number of hanging nodes
        UInt CountEdgesToRefine(UInt face) const;

        /// Deletes isolated hanging nodes(remove_isolated_hanging_nodes)
        /// @returns Number of deleted isolated hanging nodes
        [[nodiscard]] UInt DeleteIsolatedHangingnodes();

        /// @brief Connect the hanging nodes with triangles (connect_hanging_nodes)
        void ConnectHangingNodes();

        /// @brief Smooth the face and edge refinement masks (smooth_jarefine)
        void SmoothRefinementMasks();

        /// @brief Computes m_faceMask, if a face must be split later on (split_cells)
        void ComputeIfFaceShouldBeSplit();

        /// @brief The refinement operation by splitting the face (refine_cells)
        void RefineFacesBySplittingEdges();

        /// @brief Compute if an edge must be refined based on the face location type and the Courant criteria
        /// @param edge The index of the edge to be refined
        /// @param faceLocationType The face location type
        /// @returns If the edge should be refined
        bool IsEdgeToBeRefinedBasedFaceLocationTypeAndCourantCriteria(UInt edge, FaceLocation faceLocationType) const;

        /// @brief Compute if an edge must be refined based on the face location type
        /// @param edge The index of the edge to be refined
        /// @param depthValues The depth value
        /// @returns If the edge should be refined based on Courant criteria
        bool IsRefineNeededBasedOnCourantCriteria(UInt edge, double depthValues) const;

        /// @brief Compute the face location type based on the depths values on the nodes
        void ComputeFaceLocationTypes();

        /// @brief Compute which edge will be below the minimum size after refinement
        void ComputeEdgeBelowMinSizeAfterRefinement();

        std::shared_ptr<RTreeBase> m_samplesRTree; ///< The sample node RTree

        std::vector<int> m_faceMask;                           ///< Compute face without hanging nodes (1), refine face with hanging nodes (2), do not refine cell at all (0) or refine face outside polygon (-2)
        std::vector<int> m_edgeMask;                           ///< If 0, edge is not split
        std::vector<bool> m_isEdgeBelowMinSizeAfterRefinement; ///< If an edge is below the minimum size after refinement
        std::vector<int> m_nodeMask;                           ///< The node mask used in the refinement process
        std::vector<UInt> m_brotherEdges;                      ///< The index of the brother edge for each edge

        /// Local caches
        std::vector<bool> m_isHangingNodeCache;       ///< Cache for maintaining if node is hanging
        std::vector<bool> m_isHangingEdgeCache;       ///< Cache for maintaining if edge is hanging
        std::vector<Point> m_polygonNodesCache;       ///< Cache for maintaining polygon nodes
        std::vector<UInt> m_localNodeIndicesCache;    ///< Cache for maintaining local node indices
        std::vector<UInt> m_globalEdgeIndicesCache;   ///< Cache for maintaining edge indices
        std::vector<UInt> m_refineEdgeCache;          ///< Cache for the edges to be refined
        std::vector<FaceLocation> m_faceLocationType; ///< Cache for the face location types

        RefinementType m_refinementType = RefinementType::WaveCourant; ///< The type of refinement to use
        bool m_directionalRefinement = false;                          ///< Whether there is directional refinement

        std::shared_ptr<Mesh2D> m_mesh;                             ///< Pointer to the mesh
        std::shared_ptr<MeshInterpolation> m_interpolant;           ///< Pointer to the AveragingInterpolation instance
        Polygons m_polygons;                                        ///< Polygons
        MeshRefinementParameters m_meshRefinementParameters;        ///< The mesh refinement parameters
        bool m_useNodalRefinement = false;                          ///< Use refinement based on interpolated values at nodes
        const double m_mergingDistance = 0.001;                     ///< The distance for merging two edges
    };
} // namespace meshkernel
