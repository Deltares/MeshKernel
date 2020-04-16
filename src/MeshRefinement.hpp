#pragma once

#include <vector>
#include "MakeGridParametersNative.hpp"
#include "GeometryListNative.hpp"
#include "SampleRefineParametersNative.hpp"
#include "InterpolationParametersNative.hpp"
#include "Entities.hpp"
#include "Mesh.hpp"
#include "SpatialTrees.hpp"

namespace GridGeom 
{

    class Mesh;
    class Polygons;
    class Sample;

    class MeshRefinement
    {

    public:

        MeshRefinement(Mesh& mesh): 
            m_mesh(mesh)
        {
            // all gets refined
            m_faceMask.resize(m_mesh.GetNumFaces(), 1);
            m_edgeMask.resize(m_mesh.GetNumEdges(), -1);
            m_refineEdgeCache.resize(maximumNumberOfEdgesPerFace);
            m_isHangingNodeCache.resize(maximumNumberOfNodesPerFace, false);
            m_isHangingEdgeCache.resize(maximumNumberOfEdgesPerFace, false);
            m_polygonNodesCache.resize(maximumNumberOfNodesPerFace); 
            m_localNodeIndexsesCache.resize(maximumNumberOfNodesPerFace,intMissingValue);
            m_edgeIndexsesCache.resize(maximumNumberOfEdgesPerFace, intMissingValue);
        };

        ///refinecellsandfaces2
        bool RefineMeshBasedOnPoints(std::vector<Sample>& sample,
            const Polygons& polygon,
            GridGeomApi::SampleRefineParametersNative& sampleRefineParametersNative,
            GridGeomApi::InterpolationParametersNative& interpolationParametersNative);

    private:


        bool FindBrotherEdges(std::vector<int>& brotherEdges);

        ///set_initial_mask
        // do not refine, based on the criterion :
        //    cells with hanging nodes
        //    cells that are crossed by the selecting polygon
        //    ensure that no crossed cells have hanging nodes
        bool ComputeInitialRefinementMask();

        ///compute_jarefine_poly
        bool ComputeRefinementFromSamples(std::vector<Sample>& polygon);

        ///compute_jarefine_poly
        bool ComputeRefinementInPolygon(int numPolygonNodes,
            const std::vector<Sample>& samples,
            int refineType,
            bool& performRefinement);

        ///comp_jalink
        bool ComputeEdgesRefinementMask();

        ///find_hangingnodes
        bool FindHangingNodes(int faceIndex,
            int& numHangingEdges,
            int& numHangingNodes,
            int& numEdgesToRefine);

        ///TODO: smooth_jarefine

        ///split_cells
        bool SplitFaces();

        bool RefineFaces();

        Mesh& m_mesh;
        double m_deltaTimeMaxCourant;

        std::vector<int> m_faceMask;
        std::vector<int> m_edgeMask;
        std::vector<int> m_refineEdgeCache;
        std::vector<bool> m_isHangingNodeCache;
        std::vector<bool> m_isHangingEdgeCache;
        std::vector<Point> m_polygonNodesCache;
        std::vector<int> m_localNodeIndexsesCache;
        std::vector<int> m_edgeIndexsesCache;
        

        GridGeom::SpatialTrees::RTree m_rtree;

        double m_minimumFaceSize = 5e4;
        bool m_directionalRefinement = false;
        int m_maxNumEdges = 6;

        enum RefinementType 
        {
            RidgeRefinement = 1,
            WaveCourant = 2,
            MeshWidth = 3
        };


    };
}
