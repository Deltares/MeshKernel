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

#include <vector>

#include <MeshKernel/AveragingInterpolation.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/RTree.hpp>
#include <MeshKernelApi/InterpolationParameters.hpp>
#include <MeshKernelApi/SampleRefineParameters.hpp>

namespace meshkernel
{
    // Forward declarations
    class Mesh2D;

    /// @brief A class used to refine a Mesh2D
    class MeshRefinement
    {
        /// @brief Enumerator describing the different refinement types
        enum class RefinementType
        {
            RidgeRefinement = 1,
            WaveCourant = 2,
            RefinementLevels = 3
        };

    public:
        /// @brief The constructor for refining based on samples
        /// @param[in] mesh The mesh to be refined
        /// @param[in] averaging The averaging interpolation to use
        /// @param[in] sampleRefineParameters Refinement based on samples parameters
        /// @param[in] interpolationParameters Interpolation parameters
        explicit MeshRefinement(std::shared_ptr<Mesh2D> mesh,
                                std::shared_ptr<AveragingInterpolation> averaging,
                                const meshkernelapi::SampleRefineParameters& sampleRefineParameters,
                                const meshkernelapi::InterpolationParameters& interpolationParameters);

        /// @brief The constructor for refining based on polygons
        /// @param[in] mesh The mesh to be refined
        /// @param[in] polygon The polygon where to refine
        /// @param[in] interpolationParameters Interpolation parameters
        explicit MeshRefinement(std::shared_ptr<Mesh2D> mesh,
                                const Polygons& polygon,
                                const meshkernelapi::InterpolationParameters& interpolationParameters);

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
        /// 5. Connect hanging nodes if requested, DeleteIsolatedHangingnodes, ConnectHangingNodes
        void Compute();

    private:
        /// @brief Finds if two edges are brothers, sharing an hanging node. Can be moved to Mesh2D
        void FindBrotherEdges();

        /// @brief Modifies m_mesh.m_nodeMask, all nodes of the faces intersecting the polygon perimeter will get value of -2 (set_initial_mask)
        ///        The mask value of the other nodes will not be modified.
        void ComputeNodeMaskAtPolygonPerimeter();

        /// @brief Computes the edge and face refinement mask from sample values (compute_jarefine_poly)
        void ComputeRefinementMasksFromSamples();

        /// @brief Computes the number of edges that should be refined in a face (compute_jarefine_poly)
        ///        Face nodes, edge and edge lengths are stored in local caches. See Mesh2D.FaceClosedPolygon method
        /// @param face The number of face nodes
        /// @param refineEdgeCache 1 if the edge should be refined, 0 otherwise
        /// @param numEdgesToBeRefined The computed number of edges to refined
        void ComputeEdgesRefinementMaskFromSamples(size_t face,
                                                   std::vector<size_t>& refineEdgeCache,
                                                   size_t& numEdgesToBeRefined);

        /// Computes the edge refinement mask (comp_jalink)
        void ComputeEdgesRefinementMask();

        /// @brief Finds the hanging nodes in a face (find_hangingnodes)
        /// @param[in] face
        /// @param[out] numHangingEdges
        /// @param[out] numHangingNodes
        /// @param[out] numEdgesToRefine
        void FindHangingNodes(size_t face,
                              size_t& numHangingEdges,
                              size_t& numHangingNodes,
                              size_t& numEdgesToRefine);

        /// Deletes isolated hanging nodes(remove_isolated_hanging_nodes)
        /// @returns Number of deleted isolated hanging nodes
        [[nodiscard]] size_t DeleteIsolatedHangingnodes();

        /// @brief Connect the hanging nodes with triangles (connect_hanging_nodes)
        void ConnectHangingNodes();

        /// @brief Smooth the face refinement mask (smooth_jarefine)
        void SmoothEdgeRefinementMask() const;

        /// @brief Computes m_faceMask, if a face must be split later on (split_cells)
        void ComputeIfFaceShouldBeSplit();

        /// @brief The refinement operation by splitting the face (refine_cells)
        /// @param[in] numEdgesBeforeRefinement Number of edges before the refinement
        void RefineFacesBySplittingEdges(size_t numEdgesBeforeRefinement);

        RTree m_samplesRTree; ///< The sample node RTree

        std::vector<int> m_faceMask;        ///< Compute face without hanging nodes (1), refine face with hanging nodes (2), do not refine cell at all (0) or refine face outside polygon (-2)
        std::vector<int> m_edgeMask;        ///< If 0, edge is not split
        std::vector<size_t> m_brotherEdges; ///< The index of the brother edge for each edge

        /// Local caches
        std::vector<bool> m_isHangingNodeCache;       ///< Cache for maintaining if node is hanging
        std::vector<bool> m_isHangingEdgeCache;       ///< Cache for maintaining if edge is hanging
        std::vector<Point> m_polygonNodesCache;       ///< Cache for maintaining polygon nodes
        std::vector<size_t> m_localNodeIndicesCache;  ///< Cache for maintaining local node indices
        std::vector<size_t> m_globalEdgeIndicesCache; ///< Cache for maintaining edge indices

        RefinementType m_refinementType = RefinementType::WaveCourant; ///< The type of refinement to use
        bool m_directionalRefinement = false;                          ///< Whether there is directional refinement
        bool m_useMassCenters = false;                                 ///< Split cells on the mass centers

        std::shared_ptr<Mesh2D> m_mesh;                                   ///< Pointer to the mesh
        std::shared_ptr<AveragingInterpolation> m_averaging = nullptr;    ///< Pointer to the AveragingInterpolation instance
        Polygons m_polygons;                                              ///< Polygons
        meshkernelapi::SampleRefineParameters m_sampleRefineParameters;   ///< The sample parameters
        meshkernelapi::InterpolationParameters m_interpolationParameters; ///< The interpolation parameters
    };
} // namespace meshkernel
