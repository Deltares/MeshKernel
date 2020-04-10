#pragma once

#include <vector>
#include "MakeGridParametersNative.hpp"
#include "GeometryListNative.hpp"
#include "SampleRefineParametersNative.hpp"
#include "InterpolationParametersNative.hpp"
#include "Entities.hpp"
#include "Mesh.hpp"

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
        bool RefineMeshBasedOnPoints(std::vector<Point>& points,
            GridGeomApi::SampleRefineParametersNative& sampleRefineParametersNative,
            GridGeomApi::InterpolationParametersNative& interpolationParametersNative);

    private:


        bool FindBrotherEdges(std::vector<int>& brotherEdges);

        ///set_initial_mask
        // do not refine, based on the criterion :
        //    cells with hanging nodes
        //    cells that are crossed by the selecting polygon
        //    ensure that no crossed cells have hanging nodes
        bool ComputeInitialRefinementMask(std::vector<int>& brotherEdges);

        ///compute_jarefine_poly
        bool ComputeRefinementFromSamples(std::vector<Sample>& polygon);

        ///comp_jalink
        bool ComputeEdgesRefinementMask(std::vector<Point>& polygon);

        ///split_cells
        bool SplitFaces();

        ///refine_cells
        bool RefineFaces();

        /// remove isolated hanging nodes
        bool RemoveHangingNodes();

        /// connect hanging nodes
        bool ConnectHangingNodes();

        Mesh& m_mesh;
        std::vector<int> m_faceMask;
        std::vector<int> m_nodeMask;


    };
}
