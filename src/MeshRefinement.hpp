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

        MeshRefinement(Mesh& mesh, Polygons& polygons): m_mesh(mesh)
        {
            // all gets refined
            m_faceMask.resize(m_mesh.m_numFaces, 1);
            m_nodeMask.resize(m_mesh.m_nodes.size(), 1);

            //polygons.

        };

        ///refinecellsandfaces2
        bool RefineMeshBasedOnPoints(std::vector<Sample>& sample,
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

        ///comp_jalink
        bool ComputeEdgesRefinementMask(std::vector<Point>& polygon);

        ///compute_jarefine_poly
        bool ComputeRefinementInPolygon(int numPolygonNodes,
            const std::vector<Sample>& samples,
            std::vector<Point>& polygon,
            const SpatialTrees::RTree& rtree,
            double deltaCourant,
            int refineType);

        bool ComputeEdgeLengths();

        Mesh& m_mesh;
        std::vector<int> m_faceMask;
        std::vector<int> m_nodeMask;
        std::vector<double> edgeLengths;

        enum RefinementType 
        {
            RidgeRefinement = 1,
            WaveCourant = 2,
            MeshWidth = 3
        };


    };
}
