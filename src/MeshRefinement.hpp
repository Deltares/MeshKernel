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
        /// Constructor, store a reference of mesh
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        MeshRefinement(Mesh& mesh);

        /// <summary>
        /// Refine a mesh (refinecellsandfaces2)
        /// </summary>
        /// <param name="sample">The samples with values used for refinement (option 1, refine based on sample)</param>
        /// <param name="polygon">The samples with values used for refinement (option 2, refine in polygon)</param>
        /// <param name="sampleRefineParametersNative">Sample's related parameters</param>
        /// <param name="interpolationParametersNative">Interpolation parameters</param>
        /// <returns></returns>
        bool Refine(std::vector<Sample>& sample,
                    const Polygons& polygon,
                    GridGeomApi::SampleRefineParametersNative& sampleRefineParametersNative,
                    GridGeomApi::InterpolationParametersNative& interpolationParametersNative);

    private:

        /// <summary>
        /// Finds where the current edges originates from (find_linkbrothers)
        /// </summary>
        /// <returns></returns>
        bool FindParentEdges();

        /// <summary>
        /// Modifies the m_mesh.m_nodeMask, the mask where to perform the refinement (set_initial_mask)
        /// </summary>
        /// <returns></returns>
        bool ComputeNodeMask();

        /// <summary>
        /// Computes the edge and face refinement mask from samples (compute_jarefine_poly)
        /// </summary>
        /// <param name="samples"> the sample to use for computing masking</param>
        /// <returns></returns>
        bool ComputeMaskFromSamples(std::vector<Sample>& samples);

        /// <summary>
        /// Computes the edge and face refinement mask from samples for a single face (compute_jarefine_poly)
        /// Face nodes, edge and edge lenghts are stored in local caches. See Mesh.FaceClosedPolygon function
        /// </summary>
        /// <param name="numPolygonNodes"></param>
        /// <param name="samples"></param>
        /// <param name="numEdgesToBeRefined"></param>
        /// <returns></returns>
        bool ComputeEdgesRefinementFromSamplesSingleFace(int numPolygonNodes,
            std::vector<Sample>& samples,
            int& numEdgesToBeRefined);

        /// <summary>
        /// Computes the edge refinement mask (comp_jalink)
        /// </summary>
        /// <returns></returns>
        bool ComputeEdgesRefinementMask();

        /// <summary>
        /// Finds the hanging nodes in a face. Should this be a mesh class responsability? (find_hangingnodes) 
        /// </summary>
        /// <param name="faceIndex"></param>
        /// <param name="numHangingEdges"></param>
        /// <param name="numHangingNodes"></param>
        /// <param name="numEdgesToRefine"></param>
        /// <returns></returns>
        bool FindHangingNodes(int faceIndex,
            int& numHangingEdges,
            int& numHangingNodes,
            int& numEdgesToRefine);

        ///remove_isolated_hanging_nodes
        bool RemoveIsolatedHangingnodes(int& numRemovedIsolatedHangingNodes);

        ///connect_hanging_nodes
        bool ConnectHangingNodes();

        ///TODO: smooth_jarefine
        bool SmoothEdgeRefinementMask();

        ///split_cells
        bool SplitFaces();

        ///refine_cells
        bool RefineFacesBySplittingEdges(int numEdgesBeforeRefinemet);

        // compute the refinament value at the center of mass
        double ComputeFaceRefinementFromSamples(int numPolygonNodes, const std::vector<Sample>& samples, AveragingMethod averagingMethod, Point centerOfMass);

        // mesh R-Tree
        Mesh& m_mesh;

        // samples R-Tree
        SpatialTrees::RTree m_rtree;

        std::vector<int> m_faceMask; //refine cell without hanging nodes (1), refine cell with hanging nodes (2), do not refine cell at all (0) or refine cell outside polygon (-2)
        std::vector<int> m_edgeMask;
        std::vector<int> m_brotherEdges;
        std::vector<int> m_refineEdgeCache;
        std::vector<bool> m_isHangingNodeCache;
        std::vector<bool> m_isHangingEdgeCache;
        std::vector<Point> m_polygonNodesCache;
        std::vector<int> m_localNodeIndexsesCache;
        std::vector<int> m_edgeIndexsesCache;
        std::vector<double> m_polygonEdgesLengthsCache;
        std::vector<bool> m_subtractedSample;
       
        double m_deltaTimeMaxCourant = 0.0;
        double m_minimumFaceSize = 5e4;
        bool m_directionalRefinement = false;
        bool m_refineOutsideFace = false;
        bool m_connectHangingNodes = true;
        bool m_refineIntersectedFaces = false;
        int m_maxNumberOfRefinementIterations = 10;
        RefinementType m_refinementType;

    };
}
