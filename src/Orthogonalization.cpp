#include <vector>
#include <algorithm>
#include <numeric>
#include <string>

#include "Operations.cpp"
#include "Smoother.hpp"
#include "Orthogonalization.hpp"
#include "Entities.hpp"
#include "Mesh.hpp"

bool GridGeom::Orthogonalization::Set(Mesh& mesh,
    int& isTriangulationRequired,
    int& isAccountingForLandBoundariesRequired,
    int& projectToLandBoundaryOption,
    GridGeomApi::OrthogonalizationParametersNative& orthogonalizationParametersNative,
    const Polygons& polygon,
    std::vector<Point>& landBoundaries)
{
    m_maxNumNeighbours = *(std::max_element(mesh.m_nodesNumEdges.begin(), mesh.m_nodesNumEdges.end()));
    m_maxNumNeighbours += 1;
    m_nodesNodes.resize(mesh.GetNumNodes() , std::vector<int>(m_maxNumNeighbours, intMissingValue));
    m_wOrthogonalizer.resize(mesh.GetNumNodes() , std::vector<double>(m_maxNumNeighbours, 0.0));
    m_rhsOrthogonalizer.resize(mesh.GetNumNodes() , std::vector<double>(2, 0.0));
    m_aspectRatios.resize(mesh.GetNumEdges(), 0.0);
    m_polygons = polygon;
    m_smoother = Smoother();

    // Sets the node mask
    mesh.MaskNodesInPolygons(m_polygons, true);
    // Flag nodes outside the polygon as corner points
    for (auto n = 0; n < mesh.GetNumNodes(); n++)
    {
        if (mesh.m_nodeMask[n] == 0) 
        {
            mesh.m_nodesTypes[n] = 3;
        }
    }

    //for each node, determine the neighbouring nodes
    for (auto n = 0; n < mesh.GetNumNodes() ; n++)
    {
        for (auto nn = 0; nn < mesh.m_nodesNumEdges[n]; nn++)
        {
            Edge edge = mesh.m_edges[mesh.m_nodesEdges[n][nn]];
            m_nodesNodes[n][nn] = edge.first + edge.second - n;
        }
    }

    // TODO: calculate volume weights for areal smoother
    m_mumax = (1.0 - m_smoothorarea) * 0.5;
    m_mu = std::min(1e-2, m_mumax);
    m_orthogonalCoordinates.resize(mesh.GetNumNodes() );

    // in this case the nearest point is the point itself
    m_nearestPoints.resize(mesh.GetNumNodes() );
    std::iota(m_nearestPoints.begin(), m_nearestPoints.end(), 0);

    // back-up original nodes, for projection on original mesh boundary
    m_originalNodes = mesh.m_nodes;
    m_orthogonalCoordinates = mesh.m_nodes;

    // algorithm settings
    m_orthogonalizationToSmoothingFactor = orthogonalizationParametersNative.OrthogonalizationToSmoothingFactor;
    m_orthogonalizationToSmoothingFactorBoundary = orthogonalizationParametersNative.OrthogonalizationToSmoothingFactorBoundary;
    m_smoothorarea = orthogonalizationParametersNative.Smoothorarea;
    m_orthogonalizationOuterIterations = orthogonalizationParametersNative.OuterIterations;
    m_orthogonalizationBoundaryIterations = orthogonalizationParametersNative.BoundaryIterations;
    m_orthogonalizationInnerIterations = orthogonalizationParametersNative.InnerIterations;

    m_landBoundaries.Set(landBoundaries);

    m_isTriangulationRequired = isTriangulationRequired;

    m_isAccountingForLandBoundariesRequired = isAccountingForLandBoundariesRequired;

    m_projectToLandBoundaryOption = projectToLandBoundaryOption;

    // project on land boundary
    if (m_projectToLandBoundaryOption >= 1)
    {
        // account for enclosing polygon
        m_landBoundaries.Administrate(mesh, m_polygons);
        m_landBoundaries.FindNearestMeshBoundary(mesh, m_polygons, m_projectToLandBoundaryOption);
    }

    // for spherical accurate computations we need to call orthonet_comp_ops 
    if(mesh.m_projection == Projections::sphericalAccurate)
    {
        if(m_orthogonalizationToSmoothingFactor<1.0)
        {
            bool successful = PrapareOuterIteration(mesh);
            if (!successful)
            {
                return false;
            }
        }

        m_localCoordinatesIndexes.resize(mesh.GetNumNodes() + 1);
        m_localCoordinatesIndexes[0] = 1;
        for (int n = 0; n < mesh.GetNumNodes(); ++n)
        {
            m_localCoordinatesIndexes[n + 1] = m_localCoordinatesIndexes[n] + std::max(mesh.m_nodesNumEdges[n] + 1, m_smoother.GetNumConnectedNodes(n));
        }

        m_localCoordinates.resize(m_localCoordinatesIndexes.back() - 1, { doubleMissingValue, doubleMissingValue });
    }

    return true;
}

bool GridGeom::Orthogonalization::Compute(Mesh& mesh)
{
    bool successful = true;

    for (auto outerIter = 0; outerIter < m_orthogonalizationOuterIterations; outerIter++)
    {
        if (successful)
        {
            successful = PrapareOuterIteration(mesh);
        }
        for (auto boundaryIter = 0; boundaryIter < m_orthogonalizationBoundaryIterations; boundaryIter++)
        {
            for (auto innerIter = 0; innerIter < m_orthogonalizationInnerIterations; innerIter++)
            {
                if (successful)
                {
                    successful = InnerIteration(mesh);
                }        
            } // inner iteration
        } // boundary iter

        //update mu
        if (successful)
        {
            successful = FinalizeOuterIteration(mesh);
        }    
    }// outer iter

    DeallocateLinearSystem();

    return true;
}


bool GridGeom::Orthogonalization::PrapareOuterIteration(const Mesh& mesh) 
{

    bool successful = true;

    // compute aspect ratios
    if (successful)
    {
        successful = AspectRatio(mesh);
    }

    // compute weights and rhs of orthogonalizer
    if (successful)
    {
        successful = ComputeWeightsAndRhsOrthogonalizer(mesh);
    }

    // computes the smoother weights
    if (successful)
    {
        successful = m_smoother.Compute(mesh);
    }
     
    // allocate linear system for smoother and orthogonalizer
    if (successful)
    {
        successful = AllocateLinearSystem(mesh);
    }

    // compute linear system terms for smoother and orthogonalizer
    if (successful)
    {
        successful = ComputeLinearSystemTerms(mesh);
    }
    return successful;
}

bool GridGeom::Orthogonalization::AllocateLinearSystem(const Mesh& mesh) 
{
    bool successful = true;
    // reallocate caches
    if (successful && m_nodeCacheSize == 0)
    {
        m_compressedRhs.resize(mesh.GetNumNodes() * 2);
        std::fill(m_compressedRhs.begin(), m_compressedRhs.end(), 0.0);

        m_compressedEndNodeIndex.resize(mesh.GetNumNodes());
        std::fill(m_compressedEndNodeIndex.begin(), m_compressedEndNodeIndex.end(), 0.0);

        m_compressedStartNodeIndex.resize(mesh.GetNumNodes() );
        std::fill(m_compressedStartNodeIndex.begin(), m_compressedStartNodeIndex.end(), 0.0);

        for (int n = 0; n < mesh.GetNumNodes() ; n++)
        {
            m_compressedEndNodeIndex[n] = m_nodeCacheSize;
            m_nodeCacheSize += std::max(mesh.m_nodesNumEdges[n] + 1, m_smoother.GetNumConnectedNodes(n));
            m_compressedStartNodeIndex[n] = m_nodeCacheSize;
        }

        m_compressedNodesNodes.resize(m_nodeCacheSize);
        m_compressedWeightX.resize(m_nodeCacheSize);
        m_compressedWeightY.resize(m_nodeCacheSize);
    }
    return successful;
}


bool GridGeom::Orthogonalization::DeallocateLinearSystem() 
{
    m_compressedRhs.resize(0);
    m_compressedEndNodeIndex.resize(0);
    m_compressedStartNodeIndex.resize(0);
    m_compressedNodesNodes.resize(0);
    m_compressedWeightX.resize(0);
    m_compressedWeightY.resize(0);
    m_nodeCacheSize = 0;

    return true;
}

bool GridGeom::Orthogonalization::FinalizeOuterIteration(Mesh& mesh)
{
    m_mu = std::min(2.0 * m_mu, m_mumax);

    //compute new faces circumcenters
    if (!m_keepCircumcentersAndMassCenters)
    {
        mesh.ComputeFaceCircumcentersMassCentersAndAreas();
    }

    return true;
}

bool GridGeom::Orthogonalization::ComputeLinearSystemTerms(const Mesh& mesh)
{
    double max_aptf = std::max(m_orthogonalizationToSmoothingFactorBoundary, m_orthogonalizationToSmoothingFactor);
#pragma omp parallel for
	for (int n = 0; n < mesh.GetNumNodes() ; n++)
    {    
        if ((mesh.m_nodesTypes[n] != 1 && mesh.m_nodesTypes[n] != 2) || mesh.m_nodesNumEdges[n] < 2)
        {
            continue;
        }
        if (m_keepCircumcentersAndMassCenters != false && (mesh.m_nodesNumEdges[n] != 3 || mesh.m_nodesNumEdges[n] != 1))
        {
            continue;
        }

        const double atpfLoc  = mesh.m_nodesTypes[n] == 2 ? max_aptf : m_orthogonalizationToSmoothingFactor;
        const double atpf1Loc = 1.0 - atpfLoc;
        double mumat    = m_mu;
        int maxnn = m_compressedStartNodeIndex[n] - m_compressedEndNodeIndex[n];
        for (int nn = 1, cacheIndex = m_compressedEndNodeIndex[n]; nn < maxnn; nn++, cacheIndex++)
        {
            double wwx = 0.0;
            double wwy = 0.0;

            // Smoother
            if (atpf1Loc > 0.0 && mesh.m_nodesTypes[n] == 1)
            {
                wwx = atpf1Loc * m_smoother.GetWeight(n, nn);
                wwy = atpf1Loc * m_smoother.GetWeight(n, nn);
            }
            
            // Orthogonalizer
            if (nn < mesh.m_nodesNumEdges[n] + 1)
            {
                wwx += atpfLoc * m_wOrthogonalizer[n][nn - 1];
                wwy += atpfLoc * m_wOrthogonalizer[n][nn - 1];
                m_compressedNodesNodes[cacheIndex] = m_nodesNodes[n][nn - 1];
            }
            else
            {
                m_compressedNodesNodes[cacheIndex] = m_smoother.GetCoonectedNodeIndex(n,nn);
            }

            m_compressedWeightX[cacheIndex] = wwx;
            m_compressedWeightY[cacheIndex] = wwy;
        }
        int firstCacheIndex = n * 2;
        m_compressedRhs[firstCacheIndex] = atpfLoc * m_rhsOrthogonalizer[n][0];
        m_compressedRhs[firstCacheIndex+1] = atpfLoc * m_rhsOrthogonalizer[n][1];
	}

	return true;
}


bool GridGeom::Orthogonalization::InnerIteration(Mesh& mesh)
{
#pragma omp parallel for
	for (int n = 0; n < mesh.GetNumNodes() ; n++)
    {
	    UpdateNodeCoordinates(n, mesh);   
    }
	
    // update mesh node coordinates
    mesh.m_nodes = m_orthogonalCoordinates;

    // project on the original net boundary
    ProjectOnOriginalMeshBoundary(mesh);

    // compute local coordinates
    ComputeCoordinates(mesh);

    // project on land boundary
    if (m_projectToLandBoundaryOption >= 1)
    {
        m_landBoundaries.SnapMeshToLandBoundaries(mesh);
    }

    return true;
}

bool GridGeom::Orthogonalization::ProjectOnOriginalMeshBoundary(Mesh& mesh)
{
    Point firstPoint;
    Point secondPoint;
    Point thirdPoint;
    Point normalSecondPoint;
    Point normalThirdPoint;

    for (auto n = 0; n < mesh.GetNumNodes() ; n++)
    {
        int nearestPointIndex = m_nearestPoints[n];
        if (mesh.m_nodesTypes[n] == 2 && mesh.m_nodesNumEdges[n] > 0 && mesh.m_nodesNumEdges[nearestPointIndex] > 0)
        {
            firstPoint = mesh.m_nodes[n];
            int numEdges = mesh.m_nodesNumEdges[nearestPointIndex];
            int numNodes = 0;
            int leftNode;
            int rightNode;
            for (auto nn = 0; nn < numEdges; nn++)
            {
                auto edgeIndex = mesh.m_nodesEdges[nearestPointIndex][nn];
                if (mesh.m_edgesNumFaces[edgeIndex] == 1)
                {
                    numNodes++;
                    if (numNodes == 1)
                    {
                        leftNode = m_nodesNodes[n][nn];
                        if (leftNode == intMissingValue)
                        {
                            return false;
                        }
                        secondPoint = m_originalNodes[leftNode];
                    }
                    else if (numNodes == 2)
                    {
                        rightNode = m_nodesNodes[n][nn];
                        if (rightNode == intMissingValue)
                        {
                            return false;
                        }
                        thirdPoint = m_originalNodes[rightNode];
                    }
                }
            }

            //Project the moved boundary point back onto the closest original edge (either between 0 and 2 or 0 and 3)
            double rl2;
            double dis2 = DistanceFromLine(firstPoint, m_originalNodes[nearestPointIndex], secondPoint, normalSecondPoint, rl2, mesh.m_projection);

            double rl3;
            double dis3 = DistanceFromLine(firstPoint, m_originalNodes[nearestPointIndex], thirdPoint, normalThirdPoint, rl3, mesh.m_projection);

            if (dis2 < dis3)
            {
                mesh.m_nodes[n] = normalSecondPoint;
                if (rl2 > 0.5 && mesh.m_nodesTypes[n] != 3)
                {
                    m_nearestPoints[n] = leftNode;
                }
            }
            else
            {
                mesh.m_nodes[n] = normalThirdPoint;
                if (rl3 > 0.5 && mesh.m_nodesTypes[n] != 3)
                {
                    m_nearestPoints[n] = rightNode;
                }
            }
        }
    }
    return true;
}


bool GridGeom::Orthogonalization::ComputeCoordinates(const Mesh& mesh)
{
    if(mesh.m_projection == Projections::sphericalAccurate)
    {
        //TODO : missing implementation
        return true;
    }

    return true;
}

bool GridGeom::Orthogonalization::AspectRatio(const Mesh& mesh)
{
    std::vector<std::vector<double>> averageEdgesLength(mesh.GetNumEdges(), std::vector<double>(2, doubleMissingValue));
    std::vector<double> averageFlowEdgesLength(mesh.GetNumEdges(), doubleMissingValue);
    std::vector<bool> curvilinearGridIndicator(mesh.GetNumNodes() , true);
    std::vector<double> edgesLength(mesh.GetNumEdges(), 0.0);

    for (auto e = 0; e < mesh.GetNumEdges(); e++)
    {
        auto first = mesh.m_edges[e].first;
        auto second = mesh.m_edges[e].second;

        if (first == second) continue;
        double edgeLength = Distance(mesh.m_nodes[first], mesh.m_nodes[second], mesh.m_projection);
        edgesLength[e] = edgeLength;

        Point leftCenter;
        Point rightCenter;
        if (mesh.m_edgesNumFaces[e] > 0)
        {
            leftCenter = mesh.m_facesCircumcenters[mesh.m_edgesFaces[e][0]];
        }
        else
        {
            leftCenter = mesh.m_nodes[first];
        }

        //find right cell center, if it exists
        if (mesh.m_edgesNumFaces[e] == 2)
        {
            rightCenter = mesh.m_facesCircumcenters[mesh.m_edgesFaces[e][1]];
        }
        else
        {
            //otherwise, make ghost node by imposing boundary condition
            double dinry = InnerProductTwoSegments(mesh.m_nodes[first], mesh.m_nodes[second], mesh.m_nodes[first], leftCenter, mesh.m_projection);
            dinry = dinry / std::max(edgeLength * edgeLength, minimumEdgeLength);

            double x0_bc = (1.0 - dinry) * mesh.m_nodes[first].x + dinry * mesh.m_nodes[second].x;
            double y0_bc = (1.0 - dinry) * mesh.m_nodes[first].y + dinry * mesh.m_nodes[second].y;
            rightCenter.x = 2.0 * x0_bc - leftCenter.x;
            rightCenter.y = 2.0 * y0_bc - leftCenter.y;
        }

        averageFlowEdgesLength[e] = Distance(leftCenter, rightCenter, mesh.m_projection);
    }

    // Compute normal length
    for (int f = 0; f < mesh.GetNumFaces(); f++)
    {
        auto numberOfFaceNodes = mesh.GetNumFaceEdges(f);
        if (numberOfFaceNodes < 3) continue;

        for (int n = 0; n < numberOfFaceNodes; n++)
        {
            if (numberOfFaceNodes != 4) curvilinearGridIndicator[mesh.m_facesNodes[f][n]] = false;
            auto edgeIndex = mesh.m_facesEdges[f][n];

            if (mesh.m_edgesNumFaces[edgeIndex] < 1) continue;

            //get the other links in the right numbering
            //TODO: ask why only 3 are requested, why an average length stored in averageEdgesLength is needed?
            //int kkm1 = n - 1; if (kkm1 < 0) kkm1 = kkm1 + numberOfFaceNodes;
            //int kkp1 = n + 1; if (kkp1 >= numberOfFaceNodes) kkp1 = kkp1 - numberOfFaceNodes;
            //
            //std::size_t klinkm1 = mesh.m_facesEdges[f][kkm1];
            //std::size_t klinkp1 = mesh.m_facesEdges[f][kkp1];
            //

            double edgeLength = edgesLength[edgeIndex];
            if (edgeLength != 0.0)
            {
                m_aspectRatios[edgeIndex] = averageFlowEdgesLength[edgeIndex] / edgeLength;
            }

            //quads
            if (numberOfFaceNodes == 4)
            {
                int kkp2 = n + 2; if (kkp2 >= numberOfFaceNodes) kkp2 = kkp2 - numberOfFaceNodes;
                auto klinkp2 = mesh.m_facesEdges[f][kkp2];
                edgeLength = 0.5 * (edgesLength[edgeIndex] + edgesLength[klinkp2]);
            }

            if (averageEdgesLength[edgeIndex][0] == doubleMissingValue)
            {
                averageEdgesLength[edgeIndex][0] = edgeLength;
            }
            else
            {
                averageEdgesLength[edgeIndex][1] = edgeLength;
            }
        }
    }

    if (curvilinearToOrthogonalRatio == 1.0)
        return true;

    for (auto e = 0; e < mesh.GetNumEdges(); e++)
    {
        auto first = mesh.m_edges[e].first;
        auto second = mesh.m_edges[e].second;

        if (first < 0 || second < 0) continue;
        if (mesh.m_edgesNumFaces[e] < 1) continue;
        // Consider only quads
        if (!curvilinearGridIndicator[first] || !curvilinearGridIndicator[second]) continue;

        if (mesh.m_edgesNumFaces[e] == 1)
        {
            if (averageEdgesLength[e][0] != 0.0 && averageEdgesLength[e][0] != doubleMissingValue)
            {
                m_aspectRatios[e] = averageFlowEdgesLength[e] / averageEdgesLength[e][0];
            }
        }
        else
        {
            if (averageEdgesLength[e][0] != 0.0 && averageEdgesLength[e][1] != 0.0 &&
                averageEdgesLength[e][0] != doubleMissingValue && averageEdgesLength[e][1] != doubleMissingValue)
            {
                m_aspectRatios[e] = curvilinearToOrthogonalRatio * m_aspectRatios[e] +
                    (1.0 - curvilinearToOrthogonalRatio) * averageFlowEdgesLength[e] / (0.5 * (averageEdgesLength[e][0] + averageEdgesLength[e][1]));
            }
        }
    }
    return true;
}

bool GridGeom::Orthogonalization::ComputeWeightsAndRhsOrthogonalizer(const Mesh& mesh)
{
    std::fill(m_rhsOrthogonalizer.begin(), m_rhsOrthogonalizer.end(), std::vector<double>(2, 0.0));
    for (auto n = 0; n < mesh.GetNumNodes() ; n++)
    {
        if (mesh.m_nodesTypes[n] != 1 && mesh.m_nodesTypes[n] != 2)
        {
            continue;
        }

        for (auto nn = 0; nn < mesh.m_nodesNumEdges[n]; nn++)
        {

            const auto edgeIndex = mesh.m_nodesEdges[n][nn];
            double aspectRatio = m_aspectRatios[edgeIndex];
            m_wOrthogonalizer[n][nn] = 0.0;

            if (aspectRatio != doubleMissingValue)
            {
                // internal nodes
                m_wOrthogonalizer[n][nn] = aspectRatio;

                if (mesh.m_edgesNumFaces[edgeIndex] == 1)
                {
                    // boundary nodes
                    m_wOrthogonalizer[n][nn] = 0.5 * aspectRatio;

                    // compute the edge length
                    Point neighbouringNode = mesh.m_nodes[m_nodesNodes[n][nn]];
                    double neighbouringNodeDistance = Distance(neighbouringNode, mesh.m_nodes[n], mesh.m_projection);
                    double aspectRatioByNodeDistance = aspectRatio * neighbouringNodeDistance;

                    auto leftFace = mesh.m_edgesFaces[edgeIndex][0];
                    bool flippedNormal;
                    Point normal;
                    NormalVectorInside(mesh.m_nodes[n], neighbouringNode, mesh.m_facesMassCenters[leftFace], normal, flippedNormal, mesh.m_projection);
                    
                    if(mesh.m_projection==Projections::spherical && mesh.m_projection != Projections::sphericalAccurate)
                    {
                        normal.x = normal.x * std::cos(degrad_hp * 0.5 * (mesh.m_nodes[n].y + neighbouringNode.y));
                    }

                    m_rhsOrthogonalizer[n][0] +=  neighbouringNodeDistance * normal.x * 0.5;
                    m_rhsOrthogonalizer[n][1] +=  neighbouringNodeDistance * normal.y * 0.5;
                }

            }
        }

        // normalize
        double factor = std::accumulate(m_wOrthogonalizer[n].begin(), m_wOrthogonalizer[n].end(), 0.0);
        if (std::abs(factor) > 1e-14)
        {
            factor = 1.0 / factor;
            for (auto& w : m_wOrthogonalizer[n]) w = w * factor;
            m_rhsOrthogonalizer[n][0] = factor * m_rhsOrthogonalizer[n][0];
            m_rhsOrthogonalizer[n][1] = factor * m_rhsOrthogonalizer[n][1];
        }

    }
    return true;
}

bool GridGeom::Orthogonalization::GetOrthogonality(const Mesh& mesh, double* orthogonality)
{
    for(int e=0; e < mesh.GetNumEdges() ; e++)
    {
        orthogonality[e] = doubleMissingValue;
        int firstVertex = mesh.m_edges[e].first;
        int secondVertex = mesh.m_edges[e].second;

        if (firstVertex!=0 && secondVertex !=0)
        {
            if (e < mesh.GetNumEdges() && mesh.m_edgesNumFaces[e]==2 )
            {
                orthogonality[e] = NormalizedInnerProductTwoSegments(mesh.m_nodes[firstVertex], mesh.m_nodes[secondVertex],
                    mesh.m_facesCircumcenters[mesh.m_edgesFaces[e][0]], mesh.m_facesCircumcenters[mesh.m_edgesFaces[e][1]], mesh.m_projection);
                if (orthogonality[e] != doubleMissingValue)
                {
                    orthogonality[e] = std::abs(orthogonality[e]);

                }
            }
        }
    }
    return true;
}

bool GridGeom::Orthogonalization::GetSmoothness(const Mesh& mesh, double* smoothness)
{
    for (int e = 0; e < mesh.GetNumEdges(); e++)
    {
        smoothness[e] = doubleMissingValue;
        int firstVertex = mesh.m_edges[e].first;
        int secondVertex = mesh.m_edges[e].second;

        if (firstVertex != 0 && secondVertex != 0)
        {
            if (e < mesh.GetNumEdges() && mesh.m_edgesNumFaces[e] == 2)
            {
                const auto leftFace = mesh.m_edgesFaces[e][0];
                const auto rightFace = mesh.m_edgesFaces[e][1];
                const auto leftFaceArea = mesh.m_faceArea[leftFace];
                const auto rightFaceArea = mesh.m_faceArea[rightFace];

                if (leftFaceArea< minimumCellArea || rightFaceArea< minimumCellArea)
                {
                    smoothness[e] = rightFaceArea / leftFaceArea;
                }
                if (smoothness[e] < 1.0) 
                {
                    smoothness[e] = 1.0 / smoothness[e];
                }
            }
        }
    }
    return true;
}

bool GridGeom::Orthogonalization::UpdateNodeCoordinates(int nodeIndex, const Mesh& mesh)
{
    int numConnectedNodes = m_compressedStartNodeIndex[nodeIndex] - m_compressedEndNodeIndex[nodeIndex];    
    double dx0 = 0.0;
    double dy0 = 0.0;
    double increments[2]{ 0.0, 0.0 };
    for (int nn = 1, cacheIndex = m_compressedEndNodeIndex[nodeIndex]; nn < numConnectedNodes; nn++, cacheIndex++)
    {
        ComputeLocalIncrements(m_compressedWeightX[cacheIndex], m_compressedWeightY[cacheIndex], m_compressedNodesNodes[cacheIndex], nodeIndex, mesh, dx0, dy0, increments);
    }

    if (increments[0] <= 1e-8 || increments[1] <= 1e-8)
    {
        return true;
    }

    int firstCacheIndex = nodeIndex * 2;
    dx0 = (dx0 + m_compressedRhs[firstCacheIndex]) / increments[0];
    dy0 = (dy0 + m_compressedRhs[firstCacheIndex + 1]) / increments[1];
    constexpr double relaxationFactor = 0.75;
    if (mesh.m_projection == Projections::cartesian || mesh.m_projection == Projections::spherical)
    {
        double x0 = mesh.m_nodes[nodeIndex].x + dx0;
        double y0 = mesh.m_nodes[nodeIndex].y + dy0;
        static constexpr double relaxationFactorCoordinates = 1.0 - relaxationFactor;

        m_orthogonalCoordinates[nodeIndex].x = relaxationFactor * x0 + relaxationFactorCoordinates * mesh.m_nodes[nodeIndex].x;
        m_orthogonalCoordinates[nodeIndex].y = relaxationFactor * y0 + relaxationFactorCoordinates * mesh.m_nodes[nodeIndex].y;
    }
    if (mesh.m_projection == Projections::sphericalAccurate)
    {
        Point localPoint{ relaxationFactor * dx0, relaxationFactor * dy0 };

        double exxp[3];
        double eyyp[3];
        double ezzp[3];
        ComputeThreeBaseComponents(mesh.m_nodes[nodeIndex], exxp, eyyp, ezzp);

        //get 3D-coordinates in rotated frame
        Cartesian3DPoint cartesianLocalPoint;
        SphericalToCartesian(localPoint, cartesianLocalPoint);

        //project to fixed frame
        Cartesian3DPoint transformedCartesianLocalPoint;
        transformedCartesianLocalPoint.x = exxp[0] * cartesianLocalPoint.x + eyyp[0] * cartesianLocalPoint.y + ezzp[0] * cartesianLocalPoint.z;
        transformedCartesianLocalPoint.y = exxp[1] * cartesianLocalPoint.x + eyyp[1] * cartesianLocalPoint.y + ezzp[1] * cartesianLocalPoint.z;
        transformedCartesianLocalPoint.z = exxp[2] * cartesianLocalPoint.x + eyyp[2] * cartesianLocalPoint.y + ezzp[2] * cartesianLocalPoint.z;

        //tranform to spherical coordinates
        CartesianToSpherical(transformedCartesianLocalPoint, mesh.m_nodes[nodeIndex].x, m_orthogonalCoordinates[nodeIndex]);
    }
    return true;
}

bool GridGeom::Orthogonalization::ComputeLocalIncrements(double wwx, double wwy, int currentNode, int n, const Mesh& mesh, double& dx0, double& dy0, double* increments)
{

    double wwxTransformed;
    double wwyTransformed;
    if (mesh.m_projection == Projections::cartesian)
    {
        wwxTransformed = wwx;
        wwyTransformed = wwy;
        dx0 = dx0 + wwxTransformed * (mesh.m_nodes[currentNode].x - mesh.m_nodes[n].x);
        dy0 = dy0 + wwyTransformed * (mesh.m_nodes[currentNode].y - mesh.m_nodes[n].y);
    }

    if (mesh.m_projection == Projections::spherical)
    {
        wwxTransformed = wwx * earth_radius * degrad_hp *
            std::cos(0.5 * (mesh.m_nodes[n].y + mesh.m_nodes[currentNode].y) * degrad_hp);
        wwyTransformed = wwy * earth_radius * degrad_hp;

        dx0 = dx0 + wwxTransformed * (mesh.m_nodes[currentNode].x - mesh.m_nodes[n].x);
        dy0 = dy0 + wwyTransformed * (mesh.m_nodes[currentNode].y - mesh.m_nodes[n].y);

    }
    if (mesh.m_projection == Projections::sphericalAccurate)
    {
        wwxTransformed = wwx * earth_radius * degrad_hp;
        wwyTransformed = wwy * earth_radius * degrad_hp;

        dx0 = dx0 + wwxTransformed * m_localCoordinates[m_localCoordinatesIndexes[n] + currentNode - 1].x;
        dy0 = dy0 + wwyTransformed * m_localCoordinates[m_localCoordinatesIndexes[n] + currentNode - 1].y;
    }

    increments[0] += wwxTransformed;
    increments[1] += wwyTransformed;

    return true;
}
