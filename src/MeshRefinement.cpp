#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include "MeshRefinement.hpp"
#include "Mesh.hpp"
#include "Entities.hpp"
#include "SpatialTrees.hpp"
#include "Operations.cpp"


bool GridGeom::MeshRefinement::RefineMeshBasedOnPoints(std::vector<Point>& points, 
    GridGeomApi::SampleRefineParametersNative& sampleRefineParametersNative,
    GridGeomApi::InterpolationParametersNative& interpolationParametersNative)
{

    bool courantNetwork = true;
    bool waveCourantRefinementType = true;
    // triangular interpolation


    m_mesh.FindFaces();

    // build the R-Tree
    GridGeom::SpatialTrees::RTree rtree;
    rtree.BuildTree(points, m_mesh.m_projection);

    //get mesh bounds

    //find_link_broders
    std::vector<int> brotherEdges;
    FindBrotherEdges(brotherEdges);

    const int level = 10;
    for (int l = 0; l < level; l++)
    {


    }

    return true;
}

bool GridGeom::MeshRefinement::FindBrotherEdges(std::vector<int>& brotherEdges)
{

    //resize brotherEdges
    brotherEdges.resize(m_mesh.m_edges.size(), -1);

    for (int n = 0; n < m_mesh.m_nodes.size(); n++)
    {
        int ee = 0;
        for (int e = 0; e < m_mesh.m_nodesNumEdges[n]; e++)
        {
            ee = e + 1;
            if (ee >= m_mesh.m_nodesNumEdges[n]) 
            {
                ee = ee - m_mesh.m_nodesNumEdges[n];
            }

            int firstEdgeIndex = m_mesh.m_nodesEdges[n][e];
            if (m_mesh.m_edgesNumFaces[firstEdgeIndex]<1)
            {
                continue;
            }

            int secondEdgeIndex = m_mesh.m_nodesEdges[n][ee];
            if (m_mesh.m_edgesNumFaces[secondEdgeIndex]<1)
            {
                continue;
            }

            auto firstEdgeLeftFace = m_mesh.m_edgesFaces[firstEdgeIndex][0];
            auto firstEdgeRighFace = m_mesh.m_edgesFaces[firstEdgeIndex].size() == 1 ? firstEdgeLeftFace : m_mesh.m_edgesFaces[firstEdgeIndex][1];

            auto secondEdgeLeftFace = m_mesh.m_edgesFaces[secondEdgeIndex][0];
            auto secondEdgeRighFace = m_mesh.m_edgesFaces[secondEdgeIndex].size() == 1 ? secondEdgeLeftFace : m_mesh.m_edgesFaces[secondEdgeIndex][1];
         
            if (firstEdgeLeftFace!=secondEdgeLeftFace && 
                firstEdgeLeftFace!=secondEdgeRighFace &&
                firstEdgeRighFace!=secondEdgeLeftFace &&
                firstEdgeRighFace!=secondEdgeRighFace)
            {
                continue;
            }

            //check if node k is in the middle
            auto firstEdgeOtherNode = m_mesh.m_edges[firstEdgeIndex].first + m_mesh.m_edges[firstEdgeIndex].second - n;
            auto secondEdgeOtherNode = m_mesh.m_edges[secondEdgeIndex].first + m_mesh.m_edges[secondEdgeIndex].second - n;
            auto centre = (m_mesh.m_nodes[firstEdgeOtherNode] + m_mesh.m_nodes[secondEdgeOtherNode])*0.5;

            if (m_mesh.m_projection == Projections::spherical) 
            {
                double middleLatitude;
                bool successful = ComputeMiddleLatitude(m_mesh.m_nodes[firstEdgeOtherNode].y, m_mesh.m_nodes[secondEdgeOtherNode].y, middleLatitude);
                if (!successful) 
                {
                    return false;
                }

                //TODO: FINISH FOR SPHERICAL
            }

            //compute tolerance
            auto firstEdgeLength = Distance(m_mesh.m_nodes[firstEdgeOtherNode], m_mesh.m_nodes[n], m_mesh.m_projection);
            auto secondEdgeLength = Distance(m_mesh.m_nodes[secondEdgeOtherNode], m_mesh.m_nodes[n], m_mesh.m_projection);
            auto tolerance = 0.0001 * std::max(firstEdgeLength, secondEdgeLength);

            auto distanceFromCentre = Distance(centre, m_mesh.m_nodes[n], m_mesh.m_projection);
            if (distanceFromCentre<tolerance)
            {
                brotherEdges[firstEdgeIndex] = secondEdgeIndex;
                brotherEdges[secondEdgeIndex] = firstEdgeIndex;
            }
        }
    }

    return true;
}

bool GridGeom::MeshRefinement::ComputeInitialRefinementMask(std::vector<int>& brotherEdges) 
{
    bool repeat = true;

    while (repeat)
    {
        repeat = false;
        bool hanging = false;
        bool crossing = false;
        for (int f = 0; f < m_mesh.m_numFaces; f++)
        {
            int numnodes = m_mesh.m_facesNodes[f].size();
            for (int n = 0; n < numnodes; n++)
            {
                int nodeIndex = m_mesh.m_facesNodes[f][n];
                int nn = n + 1;
                if (nn > numnodes)
                {
                    nn = nn - numnodes;
                }

                int firstEdgeIndex = m_mesh.m_facesEdges[f][n];
                int secondEdgeIndex = m_mesh.m_facesEdges[f][nn];

                if (brotherEdges[firstEdgeIndex] == secondEdgeIndex || brotherEdges[secondEdgeIndex] == firstEdgeIndex)
                {
                    hanging = true;
                }

                if (m_nodeMask[nodeIndex] == 0)
                {
                    crossing = true;
                }
            }

            if (crossing)
            {
                m_faceMask[f] = 0;
                for (int n = 0; n < numnodes; n++)
                {
                    int nodeIndex = m_mesh.m_facesNodes[f][n];
                    if (m_nodeMask[nodeIndex] == 1)
                    {
                        m_nodeMask[nodeIndex] == -2;
                        repeat = true;;
                    }

                }
            }
        }
    }

    return true;
}

///compute_jarefine
bool GridGeom::MeshRefinement::ComputeRefinementFromSamples(std::vector<Sample>& samples)
{



    return true;
};

///comp_jalink
bool GridGeom::MeshRefinement::ComputeEdgesRefinementMask(std::vector<Point>& polygon)
{
    return true;
};

///split_cells
bool GridGeom::MeshRefinement::SplitFaces()
{
    return true;
};

///refine_cells
bool GridGeom::MeshRefinement::RefineFaces()
{
    return true;
};

/// remove isolated hanging nodes
bool GridGeom::MeshRefinement::RemoveHangingNodes()
{
    return true;
};

/// connect hanging nodes
bool GridGeom::MeshRefinement::ConnectHangingNodes()
{
    return true;
};

