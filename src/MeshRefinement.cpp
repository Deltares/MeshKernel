#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include "MeshRefinement.hpp"
#include "Mesh.hpp"
#include "Entities.hpp"
#include "SpatialTrees.hpp"
#include "Operations.cpp"


bool GridGeom::MeshRefinement::RefineMeshBasedOnPoints(std::vector<Sample>& sample,
    GridGeomApi::SampleRefineParametersNative& sampleRefineParametersNative,
    GridGeomApi::InterpolationParametersNative& interpolationParametersNative)
{

    bool courantNetwork = true;
    bool waveCourantRefinementType = true;
    // triangular interpolation
    m_mesh.FindFaces();

    // compute all edges lengths
    ComputeEdgeLengths();

    // build the R-Tree
    GridGeom::SpatialTrees::RTree rtree;
    std::vector<Point> points;
    rtree.BuildTree(points, m_mesh.m_projection);

    //get mesh bounds

    //find_link_broders
    m_mesh.FindBrotherEdges();

    const int level = 10;
    for (int l = 0; l < level; l++)
    {
        ComputeRefinementFromSamples(sample);
    }

    return true;
}


bool GridGeom::MeshRefinement::ComputeInitialRefinementMask()
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

                if (m_mesh.m_brotherEdges[firstEdgeIndex] == secondEdgeIndex || m_mesh.m_brotherEdges[secondEdgeIndex] == firstEdgeIndex)
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
    std::vector<Point> polygonNodesCache;
    for (int f = 0; f < m_mesh.m_numFaces; f++)
    {
        m_mesh.FacePolygon(f, polygonNodesCache);

    
    }
    return true;
};


bool GridGeom::MeshRefinement::ComputeRefinementInPolygon(int numPolygonNodes, 
    const std::vector<Sample>& samples, 
    std::vector<Point>& polygon,
    const SpatialTrees::RTree& rtree,
    double deltaCourant, 
    int refineType)
{
    Point centerOfMass;
    double area;
    bool successful = faceAreaAndCenterOfMass(polygon, numPolygonNodes, area, centerOfMass, m_mesh.m_projection);

    std::vector<double> edgesLengths(numPolygonNodes);
    for (int i = 0; i < numPolygonNodes; i++)
    {
        int nextNode = i + 1; if (nextNode >= numPolygonNodes) nextNode -= numPolygonNodes;
        double distance = Distance(polygon[i], polygon[nextNode], m_mesh.m_projection);
        edgesLengths[i] = distance;
    }
    auto minmax = std::minmax_element(edgesLengths.begin(), edgesLengths.end());
    double minLenght = *minmax.first;
    double maxLenght = *minmax.second;

    if (deltaCourant > 0 || refineType == RefinementType::WaveCourant)
    {
        double result;
        bool success = Averaging(samples, numPolygonNodes, polygon, centerOfMass, m_mesh.m_projection, rtree, AveragingMethod::KdTree, result);

        if (result == doubleMissingValue || success == false) 
        {
            return true;
        }

    }


    return true;
}

bool GridGeom::MeshRefinement::ComputeEdgeLengths()
{
    auto numEdges = m_mesh.m_edges.size();
    edgeLengths.resize(numEdges, doubleMissingValue);
    for (int e = 0; e < numEdges; e++)
    {
        int first = m_mesh.m_edges[e].first;
        int second = m_mesh.m_edges[e].second;
        edgeLengths[e] = Distance(m_mesh.m_nodes[first], m_mesh.m_nodes[second], m_mesh.m_projection);
    }
    return true;
}
