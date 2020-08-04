#pragma once

#include <vector>
#include "SampleRefineParametersNative.hpp"
#include "InterpolationParametersNative.hpp"
#include "Entities.hpp"
#include "SpatialTrees.hpp"

namespace GridGeom 
{

    class Mesh;
    class Polygons;

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
        /// Constructor, store a mesh reference
        /// </summary>
        /// <param name="mesh">Mesh to be refined</param>
        /// <returns></returns>
        MeshRefinement(Mesh& mesh);

        /// <summary>
        /// Refine a mesh (refinecellsandfaces2). Steps:
        /// 1. Masks the node to be refined (those inside a polygon)
        /// 2. Find the brother edges (FindBrotherEdges)
        /// 3. Mask nodes at the polygon perimeter (ComputeNodeMaskAtPolygonPerimeter)
        /// 4. Do mesh refinement iterations
        ///    4.1 Find the brother edges (FindBrotherEdges)
        ///    4.2 Compute edge refinement mask based on samples (ComputeRefinementMasksFromSamples)
        ///    4.3 Compute edge refinement mask based on polygon (ComputeEdgesRefinementMask)
        ///    4.3 Compute if a face should be splitted (ComputeIfFaceShouldBeSplitted)
        ///    4.4 Refine face by splitting edges (RefineFacesBySplittingEdges)
        /// 5. Connect hanging nodes if requested (RemoveIsolatedHangingnodes, ConnectHangingNodes)
        /// </summary>
        /// <param name="sample">The samples with values used for refinement (option 1, refine based on sample)</param>
        /// <param name="polygon">The samples with values used for refinement (option 2, refine in polygon)</param>
        /// <param name="sampleRefineParametersNative">Sample's related parameters</param>
        /// <param name="interpolationParametersNative">Interpolation parameters</param>
        /// <returns>If the operation succeeded</returns>
        bool Refine(std::vector<Sample>& sample,
                    const Polygons& polygon,
                    GridGeomApi::SampleRefineParametersNative& sampleRefineParametersNative,
                    GridGeomApi::InterpolationParametersNative& interpolationParametersNative);

    private:

        /// <summary>
        /// Finds if two edges are brothers, for example sharing an hanging node.
        /// </summary>
        /// <returns>If the operation succeeded</returns>
        bool FindBrotherEdges();

        /// <summary>
        /// Modifies the initial m_mesh.m_nodeMask, all mesh nodes of faces at the polygon perimeter included in the polygon will get a node mask value of -2 (set_initial_mask)
        /// The mask value of the other nodes will not be modified.
        /// </summary>
        /// <returns>If the operation succeeded</returns>
        bool ComputeNodeMaskAtPolygonPerimeter();

        /// <summary>
        /// Computes the edge and face refinement mask from sample values (compute_jarefine_poly)
        /// </summary>
        /// <param name="samples"> the sample to use for computing masking</param>
        /// <returns>If the operation succeeded</returns>
        bool ComputeRefinementMasksFromSamples(std::vector<Sample>& samples);

        /// <summary>
        /// Computes the number of edges that should be refined in a face (compute_jarefine_poly)
        /// Face nodes, edge and edge lenghts are stored in local caches. See Mesh.FaceClosedPolygon method
        /// </summary>
        /// <param name="numPolygonNodes">The number of face nodes</param>
        /// <param name="samples"> The samples to use for refinement</param>
        /// <param name="numEdgesToBeRefined"> The computed numebr of edges to be refined</param>
        /// <returns>If the operation succeeded</returns>
        bool ComputeEdgesRefinementMaskFromSamples(int numPolygonNodes,
                                                   std::vector<Sample>& samples,
                                                   int& numEdgesToBeRefined);

        /// <summary>
        /// Computes the edge refinement mask (comp_jalink)
        /// </summary>
        /// <returns>If the operation succeeded</returns>
        bool ComputeEdgesRefinementMask();

        /// <summary>
        /// Finds the hanging nodes in a face. Should this be a mesh class responsability? (find_hangingnodes) 
        /// </summary>
        /// <param name="faceIndex"></param>
        /// <param name="numHangingEdges"></param>
        /// <param name="numHangingNodes"></param>
        /// <param name="numEdgesToRefine"></param>
        /// <returns>If the operation succeeded</returns>
        bool FindHangingNodes(int faceIndex,
                              int& numHangingEdges,
                              int& numHangingNodes,
                              int& numEdgesToRefine);

        /// <summary>
        /// Remove isolated hanging nodes(remove_isolated_hanging_nodes)
        /// </summary>
        /// <param name="numRemovedIsolatedHangingNodes"></param>
        /// <returns>If the operation succeeded</returns>
        bool RemoveIsolatedHangingnodes(int& numRemovedIsolatedHangingNodes);

        /// <summary>
        /// Connect the hanging nodes with triangles (connect_hanging_nodes)
        /// </summary>
        /// <returns>If the operation succeeded</returns>
        bool ConnectHangingNodes();

        /// <summary>
        /// Smooth the face refinement mask (smooth_jarefine)
        /// </summary>
        /// <returns>If the operation succeeded</returns>
        bool SmoothEdgeRefinementMask();

        /// <summary>
        /// Computes m_faceMask, if a face must be splitted later on (split_cells)
        /// </summary>
        /// <returns>If the operation succeeded</returns>
        bool ComputeIfFaceShouldBeSplitted();

        /// <summary>
        /// Actual refinement operation by splitting the face (refine_cells)
        /// </summary>
        /// <param name="numEdgesBeforeRefinemet">Numer of edges before the refinement</param>
        /// <returns>If the operation succeeded</returns>
        bool RefineFacesBySplittingEdges(int numEdgesBeforeRefinemet);

        /// <summary>
        /// Compute the refinement value at the face center of mass
        /// </summary>
        /// <param name="numPolygonNodes">The number of polygon nodes</param>
        /// <param name="samples">The number of samples</param>
        /// <param name="averagingMethod">the averaging method to be used</param>
        /// <param name="centerOfMass">Tha face center of mass</param>
        /// <returns>The refinement value at the face center of mass</returns>
        double ComputeFaceRefinementFromSamples(int numPolygonNodes, 
                                                const std::vector<Sample>& samples, 
                                                AveragingMethod averagingMethod, 
                                                Point centerOfMass);
        // samples RTree
        SpatialTrees::RTree m_samplesRTree;              

        // refine cell without hanging nodes (1), refine cell with hanging nodes (2), do not refine cell at all (0) or refine cell outside polygon (-2)
        std::vector<int>    m_faceMask;  

        std::vector<int>    m_edgeMask;
        std::vector<int>    m_brotherEdges;
        std::vector<int>    m_refineEdgeCache;
        std::vector<bool>   m_isHangingNodeCache;
        std::vector<bool>   m_isHangingEdgeCache;
        std::vector<Point>  m_polygonNodesCache;
        std::vector<int>    m_localNodeIndexsesCache;
        std::vector<int>    m_edgeIndexsesCache;
        std::vector<double> m_polygonEdgesLengthsCache;
        std::vector<bool>   m_subtractedSample;
       
        double              m_deltaTimeMaxCourant = 0.0;
        double              m_minimumFaceSize = 5e4;
        bool                m_directionalRefinement = false;
        bool                m_refineOutsideFace = false;
        bool                m_connectHangingNodes = true;
        bool                m_refineIntersectedFaces = false;
        int                 m_maxNumberOfRefinementIterations = 10;
        RefinementType      m_refinementType;

        Mesh& m_mesh;
    };
}
