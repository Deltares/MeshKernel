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
    const Polygons& polygon,
    GridGeomApi::SampleRefineParametersNative& sampleRefineParametersNative,
    GridGeomApi::InterpolationParametersNative& interpolationParametersNative)
{

    bool courantNetwork = true;
    bool waveCourantRefinementType = true;
    m_deltaTimeMaxCourant = sampleRefineParametersNative.MaximumTimeStepInCourantGrid;

    // triangular interpolation
    m_mesh.FindFaces();

    // select nodes inside polygon
    m_mesh.SelectNodesInPolygon(polygon, true);

    // build the R-Tree
    std::vector<Point> points(sample.size());
    for (int i = 0; i < sample.size(); i++)
    {
        points[i] = { sample[i].x, sample[i].y };
    }
    m_rtree.BuildTree(points, m_mesh.m_projection);

    //find_link_broders
    m_mesh.FindBrotherEdges();

    //set_initial_mask
    ComputeInitialRefinementMask();

    const int numLevels = 10;
    int numFacesAfterRefinement = m_mesh.GetNumFaces();
    for (int level = 0; level < numLevels; level++)
    {
        bool successful = ComputeRefinementFromSamples(sample);
        if (!successful)
        {
            return false;
        }
        for (int i = 0; i < m_edgeMask.size(); i++)
        {
            m_edgeMask[i] = -m_edgeMask[i];
        }

        //disable direct refinement of cells with one or more inactive nodes
        if (level > 0)
        {
            //first level: disable cells with all nodes inactive only
            for (int f = 0; f < m_mesh.GetNumFaces(); f++)
            {
                for (int n = 0; n < m_mesh.GetNumFaces(); n++)
                {
                    auto nodeIndex = m_mesh.m_facesNodes[f][n];
                    if (m_mesh.m_nodeMask[nodeIndex] != 1)
                    {
                        m_faceMask[f] = 0;
                        break;
                    }
                }
            } 
        }
        else 
        {
            for (int f = 0; f < m_mesh.GetNumFaces(); f++)
            {
                bool activeNodeFound = false;
                for (int n = 0; n < m_mesh.GetNumFaces(); n++)
                {
                    auto nodeIndex = m_mesh.m_facesNodes[f][n];
                    if (m_mesh.m_nodeMask[nodeIndex] != 0 && m_mesh.m_nodeMask[nodeIndex] != -2)
                    {
                        //active node found : discard this cell and continue with next
                        activeNodeFound = true;
                        break;
                    }
                }
                if (!activeNodeFound) 
                {
                    m_faceMask[f] = 0;
                }
            }
        }

        successful = ComputeEdgesRefinementMask();
        if (!successful) 
        {
            return false;
        }

        successful = SplitFaces();
        if (!successful)
        {
            return false;
        }
        int numFacesToRefine = 0;
        for (auto f = 0; f < m_mesh.GetNumFaces(); f++)
        {
            if (m_faceMask[f] != 0) 
            {
                numFacesToRefine++;
            }
        }
        if (numFacesToRefine == 0) 
        {
            return true;
        }
        numFacesAfterRefinement = numFacesAfterRefinement * 4;
        
        // do refinement
    }

    return true;
}


bool GridGeom::MeshRefinement::RefineFaces() 
{

    //Add new nodes where required
    std::vector<int> localFaceNodes(m_maxNumEdges, intMissingValue);
    std::vector<bool> ishanging(m_maxNumEdges, false);
    std::vector<Point> facePolygonWithoutHangingNodes(m_maxNumEdges);
    const int maxNumEdgesAfterRefinement = 10;
    for (int e = 0; e < m_mesh.GetNumEdges(); e++)
    {
        //Compute the center of the edge
        auto firstNodeIndex = m_mesh.m_edges[e].first;
        auto secondNodeIndex = m_mesh.m_edges[e].second;
        auto firstNode = m_mesh.m_nodes[firstNodeIndex];
        auto secondNode = m_mesh.m_nodes[secondNodeIndex];

        Point middle{ (firstNode.x + secondNode.x) / 2.0, (firstNode.y + secondNode.y) / 2.0 };
        if (m_mesh.m_projection == Projections::spherical) 
        {
            double middlelatitude;
            bool successful = ComputeMiddleLatitude(firstNode.y, secondNode.y, middle.y);
            if (!successful) 
            {
                return false;
            }

            if (std::abs(firstNode.x - secondNode.x) > 180.0) 
            {
                middle.x += 180.0;
            }

            // fix at poles
            bool firstNodeAtPole = std::abs(std::abs(firstNode.y) - 90.0) < dtol_pole;
            bool secondNodeAtPole = std::abs(std::abs(secondNode.y) - 90.0) < dtol_pole;
            if (firstNodeAtPole && !secondNodeAtPole) 
            {
                middle.x = secondNode.x;
            }
            else if (!firstNodeAtPole && secondNodeAtPole)
            {
                middle.x = firstNode.x;
            }
        }


        int newNodeIndex;
        m_mesh.InsertNode(middle, newNodeIndex);
        m_edgeMask[e] = newNodeIndex;

        // set mask on the new node
        m_mesh.m_nodeMask[newNodeIndex] = 1;
        if (m_mesh.m_nodeMask[firstNodeIndex] == 0 && m_mesh.m_nodeMask[secondNodeIndex] == 0) 
        {
            m_mesh.m_nodeMask[newNodeIndex] = 0;
        }
        else if (m_mesh.m_nodeMask[firstNodeIndex] != 1 || m_mesh.m_nodeMask[secondNodeIndex] != 1)
        {
            m_mesh.m_nodeMask[newNodeIndex] = -1;
        }
    }

    for (int f = 0; f < m_mesh.GetNumFaces(); f++)
    {
        if (m_faceMask[f]==0)
        {
            continue;
        }
        
        bool isFullFaceNotInPolygon = m_mesh.IsFullFaceNotInPolygon(f);

        auto numNodes = m_mesh.GetNumFaceNodes(f);

        bool successful = m_mesh.FacePolygon(f, m_polygonNodesCache, m_localNodeIndexsesCache, m_edgeIndexsesCache);
        if (!successful) 
        {
            return false;
        }


        int numBrotherEdges = 0;
        int numNonHangingNodes = 0;
        std::fill(localFaceNodes.begin(), localFaceNodes.end(), intMissingValue);
        for (int n = 0; n < numNodes; n++)
        {
            auto firstEdge = NextCircularBackwardIndex(n, numNodes);
            auto secondEdge = NextCircularForwardIndex(n, numNodes);

            int mappedNode = m_localNodeIndexsesCache[n];
            auto edgeIndex = m_mesh.m_facesEdges[f][mappedNode];

            mappedNode = m_localNodeIndexsesCache[firstEdge];
            auto firstEdgeIndex = m_mesh.m_facesEdges[f][mappedNode];

            mappedNode = m_localNodeIndexsesCache[secondEdge];
            auto secondEdgeIndex = m_mesh.m_facesEdges[f][mappedNode];

            if (edgeIndex < 0) 
            {
                continue;
            }

            if (m_mesh.m_brotherEdges[edgeIndex] != secondEdgeIndex && secondEdgeIndex > 0)
            {
                numBrotherEdges++;
                int newNode;
                m_mesh.FindCommonNode(edgeIndex, m_mesh.m_brotherEdges[edgeIndex], newNode);
                numNonHangingNodes++;
                localFaceNodes[numNonHangingNodes] = newNode;
            }
            else if (m_mesh.m_brotherEdges[edgeIndex] != firstEdgeIndex || m_mesh.m_brotherEdges[edgeIndex] == intMissingValue)
            {
                if (m_edgeMask[edgeIndex] != 0) 
                {
                    numNonHangingNodes++;
                    localFaceNodes[numNonHangingNodes] = m_edgeMask[edgeIndex];
                }
            }

            if (numNonHangingNodes >= maxNumEdgesAfterRefinement) 
            {
                return true;
            }

            // check if start of this link is hanging
            if (m_mesh.m_brotherEdges[edgeIndex]== firstEdgeIndex && firstEdgeIndex>=0)
            {
                ishanging[n] = true;
            }
        }


        //compute new center node : circumcenter without hanging nodes for quads, c / g otherwise
        int numvalidNodes = 0;
        for (int n = 0; n < numNodes; n++)
        {    
            if (!ishanging[n])
            {
                facePolygonWithoutHangingNodes[numvalidNodes] = m_polygonNodesCache[n];
                numvalidNodes++;

                int mappedNode = m_localNodeIndexsesCache[n];
                if (mappedNode > 0) 
                {

                
                }
                else 
                {
                
                }

            }

        }






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
        for (int f = 0; f < m_mesh.GetNumFaces(); f++)
        {
            int numnodes = m_mesh.GetNumFaceNodes(f);
            for (int n = 0; n < numnodes; n++)
            {
                int nodeIndex = m_mesh.m_facesNodes[f][n];
                auto nn = NextCircularForwardIndex(n, numnodes);


                int firstEdgeIndex = m_mesh.m_facesEdges[f][n];
                int secondEdgeIndex = m_mesh.m_facesEdges[f][nn];

                if (m_mesh.m_brotherEdges[firstEdgeIndex] == secondEdgeIndex || m_mesh.m_brotherEdges[secondEdgeIndex] == firstEdgeIndex)
                {
                    hanging = true;
                }

                if (m_mesh.m_nodeMask[nodeIndex] == 0)
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
                    if (m_mesh.m_nodeMask[nodeIndex] == 1)
                    {
                        m_mesh.m_nodeMask[nodeIndex] == -2;
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
    for (int f = 0; f < m_mesh.GetNumFaces(); f++)
    {
        m_mesh.FacePolygon(f, m_polygonNodesCache, m_localNodeIndexsesCache, m_edgeIndexsesCache);

        int numHangingEdges;
        int numHangingNodes;
        int numEdgesToRefine;
        bool successful = FindHangingNodes(f, numHangingEdges, numHangingNodes, numEdgesToRefine);
        if (successful) 
        {
            return false;
        }

        int numPolygonNodes = m_mesh.GetNumFaces();

        bool performRefinement;
        ComputeRefinementInPolygon(numPolygonNodes, samples, RefinementType::WaveCourant, performRefinement);
        if (performRefinement)
        {
            for (int n = 0; n < m_mesh.GetNumNodes(); n++)
            {
                int node = m_localNodeIndexsesCache[n];
                if (node < 0) 
                {
                    continue;
                }
                int edgeIndex = m_mesh.m_facesEdges[f][node];
                if (edgeIndex >= 0)
                {
                    m_edgeMask[edgeIndex] = 1;
                }
            }
        }
    }
    return true;
};

bool GridGeom::MeshRefinement::FindHangingNodes(int faceIndex,
    int& numHangingEdges,
    int& numHangingNodes,
    int& numEdgesToRefine)
{

    numEdgesToRefine = 0;
    numHangingEdges = 0;
    numHangingNodes = 0;
    auto numFaceNodes = m_mesh.GetNumFaceNodes(faceIndex);
    std::fill(m_isHangingNodeCache.begin(), m_isHangingNodeCache.end(), false);
    std::fill(m_isHangingEdgeCache.begin(), m_isHangingEdgeCache.end(), false);

    int kknod = 0;

    for (int n = 0; n < numFaceNodes; n++)
    {
        auto edgeIndex = m_mesh.m_facesEdges[faceIndex][n];
        if (m_edgeMask[edgeIndex] != 0)
        {
            numEdgesToRefine += 1;
        }

        // check if the brother link is in the cell
        if (m_mesh.m_brotherEdges[edgeIndex] != intMissingValue)
        {
            auto e = NextCircularBackwardIndex(n, numFaceNodes);
            auto ee = NextCircularForwardIndex(n, numFaceNodes);

            int commonNode = intMissingValue;
            if (m_mesh.m_brotherEdges[edgeIndex] == m_mesh.m_facesEdges[faceIndex][e])
            {

            }
            else if (m_mesh.m_brotherEdges[edgeIndex] == m_mesh.m_facesEdges[faceIndex][ee])
            {

            }

            if (commonNode != intMissingValue)
            {
                m_isHangingEdgeCache[n] = true;
                numHangingEdges++;
                for (int nn = 0; nn < numFaceNodes; nn++)
                {
                    kknod = kknod + 1;
                    if (kknod >= numFaceNodes)
                    {
                        kknod = kknod - numFaceNodes;
                    }

                    if (m_mesh.m_facesNodes[faceIndex][n] == commonNode && m_isHangingNodeCache[kknod] == 0)
                    {
                        numHangingNodes++;
                        m_isHangingNodeCache[kknod] = true;
                    }
                }
            }
        }
    }

    return true;
}


bool GridGeom::MeshRefinement::ComputeRefinementInPolygon(int numPolygonNodes,
    const std::vector<Sample>& samples,
    int refineType,
    bool& performRefinement)
{
    performRefinement = false;
    std::fill(m_refineEdgeCache.begin(), m_refineEdgeCache.end(), 0);

    //find center of mass
    Point centerOfMass;
    double area;
    bool successful = faceAreaAndCenterOfMass(m_polygonNodesCache, numPolygonNodes, area, centerOfMass, m_mesh.m_projection);

    std::vector<double> edgesLengths(numPolygonNodes);
    for (int i = 0; i < numPolygonNodes; i++)
    {
        int nextNode = i + 1; if (nextNode >= numPolygonNodes) nextNode -= numPolygonNodes;
        double distance = Distance(m_polygonNodesCache[i], m_polygonNodesCache[nextNode], m_mesh.m_projection);
        edgesLengths[i] = distance;
    }
    auto minmax = std::minmax_element(edgesLengths.begin(), edgesLengths.end());
    double minLenght = *minmax.first;
    double maxLenght = *minmax.second;

    double refinementValue;
    if (m_deltaTimeMaxCourant > 0 || refineType == RefinementType::WaveCourant)
    {

        bool success = Averaging(samples, numPolygonNodes, m_polygonNodesCache, centerOfMass, m_mesh.m_projection, m_rtree, AveragingMethod::KdTree, refinementValue);
        if (!success)
        {
            return true;
        }
    }
    else
    {
        refinementValue = 0.0;
    }


    if (refinementValue == doubleMissingValue)
    {
        return true;
    }

    // always wave courant
    int numEdgesToBeRefined = 0;
    for (int i = 0; i < numPolygonNodes; i++)
    {
        if (edgesLengths[i] < mergingDistance)
        {
            numEdgesToBeRefined++;
            continue;
        }

        double c = std::sqrt(gravity*std::abs(refinementValue));
        double waveCourant = c * m_deltaTimeMaxCourant / edgesLengths[i];
        double newEdgeLength = 0.5 * edgesLengths[i];

        if (waveCourant < 1 && std::abs(newEdgeLength - m_minimumFaceSize) < std::abs(edgesLengths[i] - m_minimumFaceSize))
        {
            numEdgesToBeRefined++;
            m_refineEdgeCache[i] = 1;
        }
        else
        {
            m_refineEdgeCache[i] = 0;
        }
    }

    if (numEdgesToBeRefined > 0)
    {
        numEdgesToBeRefined = 0;
        for (int i = 0; i < numPolygonNodes; i++)
        {
            if (m_refineEdgeCache[i] == 1 || m_isHangingNodeCache[i])
            {
                numEdgesToBeRefined++;
            }
        }
    }

    if (!m_directionalRefinement)
    {
        if (numEdgesToBeRefined == numPolygonNodes)
        {
            for (int i = 0; i < numPolygonNodes; i++)
            {
                if (!m_isHangingNodeCache[i])
                {
                    m_refineEdgeCache[i] = 1;
                }
            }
        }
        else
        {
            numEdgesToBeRefined = 0;
            std::fill(m_refineEdgeCache.begin(), m_refineEdgeCache.end(), 0);
        }
    }

    if (numEdgesToBeRefined == 0)
    {
        performRefinement = false;
        return true;
    }
    if (numEdgesToBeRefined == 1)
    {
        performRefinement = false;
        std::fill(m_refineEdgeCache.begin(), m_refineEdgeCache.end(), 0);
        return true;
    }

    performRefinement = true;
    return true;

}

bool  GridGeom::MeshRefinement::ComputeEdgesRefinementMask() 
{
    bool repeat = true;
    int iter = 0;
    const int numMaxIterations = 6;
    std::vector<int> isQuadEdge;
    std::vector<int> numOfEdges(maximumNumberOfEdgesPerFace);

    while (repeat && iter < numMaxIterations)
    {
        repeat = false;
        iter++;

        for (int f = 0; f < m_mesh.GetNumFaces(); f++)
        {
            if (m_faceMask[f] == 0)
            {
                continue;
            }

            int numHangingEdges;
            int numHangingNodes;
            int numEdgesToRefine;
            bool successful = FindHangingNodes(f, numHangingEdges, numHangingNodes, numEdgesToRefine);
            if (!successful)
            {
                return false;
            }

            int numFaceNodes = m_mesh.GetNumFaceNodes(f);
            int numNodesEffective = numFaceNodes - numHangingNodes;

            if (numNodesEffective != 4)
            {
                //non-quads
                for (int n = 0; n < numFaceNodes; n++)
                {
                    auto e = NextCircularBackwardIndex(n, numFaceNodes);
                    auto ee = NextCircularForwardIndex(n, numFaceNodes);

                    int edgeIndex = m_mesh.m_facesEdges[f][n];

                    int firstEdgeIndex = m_mesh.m_facesEdges[f][e];
                    int secondEdgeIndex = m_mesh.m_facesEdges[f][ee];

                    if (m_mesh.m_brotherEdges[edgeIndex] != firstEdgeIndex && m_mesh.m_brotherEdges[edgeIndex] != secondEdgeIndex)
                    {
                        m_edgeMask[edgeIndex] = 1;
                    }
                }
            }
            else
            {

                // number the links in the cell, links that share a hanging node will have the same number

                int num = 0;
                //non-quads
                for (int n = 0; n < numFaceNodes; n++)
                {
                    int edgeIndex = m_mesh.m_facesEdges[f][n];
                    numOfEdges[n] = num;

                    if (m_edgeMask[edgeIndex] != 0)
                    {
                        isQuadEdge[num] = m_edgeMask[edgeIndex];
                    }

                    auto ee = NextCircularForwardIndex(n, numFaceNodes);

                    int secondEdgeIndex = m_mesh.m_facesEdges[f][ee];

                    if (n != numFaceNodes - 1 && m_mesh.m_brotherEdges[edgeIndex] != secondEdgeIndex)
                    {
                        num++;
                    }
                    else if (m_mesh.m_brotherEdges[edgeIndex] != secondEdgeIndex)
                    {
                        isQuadEdge[num] = 1;
                    }
                }

                if (num != 4)
                {
                    return false;
                }

                int numEdgesToRefine = 0;
                int firstEdgeIndex = 0;
                int secondEdgeIndex = 0;
                for (int i = 0; i < 4; i++)
                {
                    if (isQuadEdge[i] != 0)
                    {
                        numEdgesToRefine++;
                        if (firstEdgeIndex == 0)
                        {
                            firstEdgeIndex = i;
                        }
                        else if (secondEdgeIndex == 0)
                        {
                            secondEdgeIndex = i;
                        }
                    }
                }

                int numEdgesFirstSecond = secondEdgeIndex - firstEdgeIndex;
                bool refineAllEdges = false;
                if (refineAllEdges == 2 && (numEdgesFirstSecond == 1 || numEdgesFirstSecond == 3))
                {
                    repeat = true;
                    refineAllEdges = true;
                }

                for (int n = 0; n < numFaceNodes; n++)
                {

                    int edgeIndex = m_mesh.m_facesEdges[f][n];
                    if (refineAllEdges != true && m_edgeMask[edgeIndex] != -1)
                    {
                        continue;

                    }

                    auto e = NextCircularBackwardIndex(n, numFaceNodes);
                    auto ee = NextCircularForwardIndex(n, numFaceNodes);

                    int secondEdgeIndex = m_mesh.m_facesEdges[f][ee];

                    if (isQuadEdge[n] != isQuadEdge[e] && isQuadEdge[n] != isQuadEdge[ee])
                    {
                        m_edgeMask[edgeIndex] = 1;
                    }

                }
            }
        }
    }

    if (repeat) 
    {
        // solution did not converge
        return false;
    }

    //only keep jalink=1, set other values to 0
    for (int i = 0; i < m_edgeMask.size(); i++)
    {
        if (m_edgeMask[i] != 1) 
        {
            m_edgeMask[i] = 0;
        }
    }

    return true;
}

bool  GridGeom::MeshRefinement::SplitFaces()
{
    const int maxiter = 1000;
    int num = 0;
    int iter = 0;
    while (num != 0)
    {
        iter++;
        if (iter > maxiter)
        {
            break;
        }

        num = 0;
        for (int f = 0; f < m_mesh.GetNumFaces(); f++)
        {
            if (m_faceMask[f] != 0 && m_faceMask[f] != -1)
            {
                continue;
            }
            int numHangingEdges;
            int numHangingNodes;
            int numEdgesToRefine;
            bool successful = FindHangingNodes(f, numHangingEdges, numHangingNodes, numEdgesToRefine);

            bool isSplittingRequired = false;

            // check if the edge has a brother edge and needs to be refined
            auto numFaceNodes = m_mesh.GetNumFaceNodes(f);
            for (int n = 0; n < numFaceNodes; n++)
            {
                int edgeIndex = m_mesh.m_facesEdges[f][n];
                if (m_isHangingEdgeCache[edgeIndex] && m_edgeMask[edgeIndex] > 0)
                {
                    isSplittingRequired = true;
                }
            }

            //compute the effective face type
            if (numHangingEdges % 2 == 1)
            {
                //uneven number of brotherlinks   
                // TODO: ADD DOT
            }

            int numNodesEffective = numFaceNodes - numHangingNodes / 2;
            if (numFaceNodes + numEdgesToRefine > m_maxNumEdges ||          // would result in unsupported cells after refinement
                numFaceNodes - numHangingNodes - numEdgesToRefine < 1 ||  // cells with only one unrefined edge
                numNodesEffective == numEdgesToRefine)                    // refine all edges
            {
                isSplittingRequired = true;
            }

            if (isSplittingRequired)
            {
                if (m_faceMask[f] != -1)
                {
                    m_faceMask[f] = 2;
                }
                else
                {
                    m_faceMask[f] = -2;
                }

                for (int n = 0; n < numFaceNodes; n++)
                {
                    int edgeIndex = m_mesh.m_facesEdges[f][n];
                    if (!m_isHangingEdgeCache[edgeIndex] && m_edgeMask[edgeIndex] == 0)
                    {
                        m_edgeMask[edgeIndex] = 1;
                        num++;
                    }
                    if (iter == maxiter)
                    {
                        //TODO: ADD DOT
                    }
                }

            }
        }
    }

    return true;
}