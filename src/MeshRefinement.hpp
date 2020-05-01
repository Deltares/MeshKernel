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

    public:

        MeshRefinement(Mesh& mesh);

        ///refinecellsandfaces2
        bool Refine(std::vector<Sample>& sample,
            const Polygons& polygon,
            GridGeomApi::SampleRefineParametersNative& sampleRefineParametersNative,
            GridGeomApi::InterpolationParametersNative& interpolationParametersNative);

    private:

        ///find_linkbrothers
        bool FindBrotherEdges();

        ///set_initial_mask
        // do not refine, based on the criterion :
        //    cells with hanging nodes
        //    cells that are crossed by the selecting polygon
        //    ensure that no crossed cells have hanging nodes
        bool ComputeInitialRefinementMask();

        ///compute_jarefine_poly
        bool ComputeEdgeAndFaceRefinementMaskFromSamples(std::vector<Sample>& polygon);

        ///compute_jarefine_poly
        bool ComputeLocalEdgeRefinementFromSamples(int faceindex, 
            int numPolygonNodes,
            const std::vector<Sample>& samples,
            int refineType,
            int& numEdgesToBeRefined);

        ///comp_jalink
        bool ComputeEdgesRefinementMask();

        ///find_hangingnodes
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

        bool RefineFaces(int numEdgesBeforeRefinemet);

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
       
        double m_deltaTimeMaxCourant = 0.0;
        double m_minimumFaceSize = 5e4;
        bool m_directionalRefinement = false;
        bool m_refineOutsideFace = false;
        bool m_connectHangingNodes = true;
        int m_maxNumberOfRefinementIterations = 10;

        enum RefinementType 
        {
            RidgeRefinement = 1,
            WaveCourant = 2,
            MeshWidth = 3
        };


    };
}
