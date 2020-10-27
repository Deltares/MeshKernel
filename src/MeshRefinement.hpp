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
#include "SampleRefineParametersNative.hpp"
#include "InterpolationParametersNative.hpp"
#include "Entities.hpp"
#include "SpatialTrees.hpp"

namespace meshkernel
{
    // Forward declarations
    class Mesh;
    class Polygons;
    class AveragingInterpolation;

    class MeshRefinement
    {
        enum class RefinementType
        {
            RidgeRefinement,
            WaveCourant,
            RefinementLevels
        };

    public:
        /// <summary>
        /// Constructor, pass a mesh reference
        /// </summary>
        /// <param name="mesh">The mesh to be refined</param>
        /// <returns></returns>
        explicit MeshRefinement(std::shared_ptr<Mesh> mesh, std::shared_ptr<AveragingInterpolation> averaging);

        explicit MeshRefinement(std::shared_ptr<Mesh> mesh);

        /// <summary>
        /// Refine a mesh (refinecellsandfaces2). Steps:
        /// 1. Masks the node to be refined (those inside a polygon)
        /// 2. Find the brother edges, the edge sharing a hanging node, FindBrotherEdges
        /// 3. Mask nodes at the polygon perimeter, ComputeNodeMaskAtPolygonPerimeter
        /// 4. Do refinement iterations
        ///    4.1 Find the brother edges, FindBrotherEdges
        ///    4.2 Compute the edge refinement mask based on samples, ComputeRefinementMasksFromSamples
        ///    4.3 Compute the edge refinement mask based on polygon, ComputeEdgesRefinementMask
        ///    4.3 Compute if a face should be splitted, ComputeIfFaceShouldBeSplitted
        ///    4.4 Refine face by splitting edges, RefineFacesBySplittingEdges
        /// 5. Connect hanging nodes if requested, RemoveIsolatedHangingnodes, ConnectHangingNodes
        /// </summary>
        /// <param name="sample">The samples with refinement levels (option 1, refine based on sample)</param>
        /// <param name="polygon">The polygon where to perform refinement (option 2, refine in polygon)</param>
        /// <param name="sampleRefineParametersNative">Refinement based on samples parameters</param>
        /// <param name="interpolationParametersNative">Interpolation parameters</param>
        /// <returns>If the method succeeded</returns>
        bool Refine(const Polygons& polygon,
                    const meshkernelapi::SampleRefineParametersNative& sampleRefineParametersNative,
                    const meshkernelapi::InterpolationParametersNative& interpolationParametersNative);

    private:
        /// <summary>
        /// Finds if two edges are brothers, sharing an hanging node. Can be moved to Mesh
        /// </summary>
        /// <returns>If the method succeeded</returns>
        bool FindBrotherEdges();

        /// <summary>
        /// Modifies m_mesh.m_nodeMask, all nodes of the faces intersecting the polygon perimeter will get value of -2 (set_initial_mask)
        /// The mask value of the other nodes will not be modified.
        /// </summary>
        /// <returns>If the method succeeded</returns>
        bool ComputeNodeMaskAtPolygonPerimeter();

        /// <summary>
        /// Computes the edge and face refinement mask from sample values (compute_jarefine_poly)
        /// </summary>
        /// <param name="samples"> The sample to use for computing masking</param>
        /// <returns>If the method succeeded</returns>
        bool ComputeRefinementMasksFromSamples();

        /// <summary>
        /// Computes the number of edges that should be refined in a face (compute_jarefine_poly)
        /// Face nodes, edge and edge lenghts are stored in local caches. See Mesh.FaceClosedPolygon method
        /// </summary>
        /// <param name="numPolygonNodes">The number of face nodes</param>
        /// <param name="samples"> The samples to use for refinement</param>
        /// <param name="refineEdgeCache"> 1 if the edge should be refined, 0 otherwise</param>
        /// <param name="numEdgesToBeRefined"> The computed number of edges to refined</param>
        /// <returns>If the method succeeded</returns>
        bool ComputeEdgesRefinementMaskFromSamples(int face,
                                                   std::vector<int>& refineEdgeCache,
                                                   int& numEdgesToBeRefined);

        /// <summary>
        /// Computes the edge refinement mask (comp_jalink)
        /// </summary>
        /// <returns>If the method succeeded</returns>
        bool ComputeEdgesRefinementMask();

        /// <summary>
        /// Finds the hanging nodes in a face (find_hangingnodes)
        /// </summary>
        /// <param name="face"></param>
        /// <param name="numHangingEdges"></param>
        /// <param name="numHangingNodes"></param>
        /// <param name="numEdgesToRefine"></param>
        /// <returns>If the method succeeded</returns>
        bool FindHangingNodes(int face,
                              int& numHangingEdges,
                              int& numHangingNodes,
                              int& numEdgesToRefine);

        /// <summary>
        /// Remove isolated hanging nodes(remove_isolated_hanging_nodes)
        /// </summary>
        /// <param name="numRemovedIsolatedHangingNodes"></param>
        /// <returns>If the method succeeded</returns>
        bool RemoveIsolatedHangingnodes(int& numRemovedIsolatedHangingNodes);

        /// <summary>
        /// Connect the hanging nodes with triangles (connect_hanging_nodes)
        /// </summary>
        /// <returns>If the method succeeded</returns>
        bool ConnectHangingNodes();

        /// <summary>
        /// Smooth the face refinement mask (smooth_jarefine)
        /// </summary>
        /// <returns>If the method succeeded</returns>
        bool SmoothEdgeRefinementMask() const;

        /// <summary>
        /// Computes m_faceMask, if a face must be splitted later on (split_cells)
        /// </summary>
        /// <returns>If the method succeeded</returns>
        bool ComputeIfFaceShouldBeSplitted();

        /// <summary>
        /// The refinement operation by splitting the face (refine_cells)
        /// </summary>
        /// <param name="numEdgesBeforeRefinemet">Numer of edges before the refinement</param>
        /// <returns>If the method succeeded</returns>
        bool RefineFacesBySplittingEdges(int numEdgesBeforeRefinemet);

        /// <summary>
        /// Compute the refinement value at the face center of mass
        /// </summary>
        /// <param name="numPolygonNodes">The number of polygon nodes</param>
        /// <param name="samples">The number of samples</param>
        /// <param name="averagingMethod">The averaging method to used</param>
        /// <param name="centerOfMass">Tha face center of mass</param>
        /// <returns>The refinement value at the face center of mass</returns>
        double ComputeFaceRefinementFromSamples(int numPolygonNodes,
                                                const std::vector<Sample>& samples,
                                                AveragingMethod averagingMethod,
                                                Point centerOfMass);
        /// The sample node RTree
        SpatialTrees::RTree m_samplesRTree;

        std::vector<int> m_faceMask;     /// Refine face without hanging nodes (1), refine face with hanging nodes (2), do not refine cell at all (0) or refine face outside polygon (-2)
        std::vector<int> m_edgeMask;     /// If 0, edge is not splitted
        std::vector<int> m_brotherEdges; /// The index of the brother edge for each edge

        /// Local caches
        std::vector<int> m_refineEdgeCache;
        std::vector<bool> m_isHangingNodeCache;
        std::vector<bool> m_isHangingEdgeCache;
        std::vector<Point> m_polygonNodesCache;
        std::vector<int> m_localNodeIndicesCache;
        std::vector<int> m_edgeIndicesCache;

        std::vector<bool> m_subtractedSample; /// Is the sample value subtracted (e.g. in refinement by levels)

        double m_deltaTimeMaxCourant = 0.0;
        double m_minimumFaceSize = 5e4;
        bool m_directionalRefinement = false;
        bool m_refineOutsideFace = false;
        bool m_connectHangingNodes = true;
        bool m_refineIntersectedFaces = false;
        int m_maxNumberOfRefinementIterations = 10;

        RefinementType m_refinementType; /// The type of refinement to use

        std::shared_ptr<Mesh> m_mesh;
        std::shared_ptr<AveragingInterpolation> m_averaging;
    };
} // namespace meshkernel
