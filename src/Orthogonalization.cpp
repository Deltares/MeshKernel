#include <vector>
#include <algorithm>
#include <numeric>
#include <string>

#include "Operations.cpp"
#include "Orthogonalization.hpp"

bool GridGeom::Orthogonalization::Initialize(const Mesh& mesh,
    int& isTriangulationRequired,
    int& isAccountingForLandBoundariesRequired,
    int& projectToLandBoundaryOption,
    GridGeomApi::OrthogonalizationParametersNative& orthogonalizationParametersNative,
    GridGeomApi::GeometryListNative& geometryListNativePolygon,
    GridGeomApi::GeometryListNative& geometryListNativeLandBoundaries)
{
    m_maxNumNeighbours = *(std::max_element(mesh.m_nodesNumEdges.begin(), mesh.m_nodesNumEdges.end()));
    m_maxNumNeighbours += 1;
    m_nodesNodes.resize(mesh.m_nodes.size(), std::vector<int>(m_maxNumNeighbours, intMissingValue));
    m_weights.resize(mesh.m_nodes.size(), std::vector<double>(m_maxNumNeighbours, 0.0));
    m_rightHandSide.resize(mesh.m_nodes.size(), std::vector<double>(2, 0.0));
    m_aspectRatios.resize(mesh.m_edges.size(), 0.0);

    //for each node, determine the neighbouring nodes
    for (std::size_t n = 0; n < mesh.m_nodes.size(); n++)
    {
        for (std::size_t nn = 0; nn < mesh.m_nodesNumEdges[n]; nn++)
        {
            Edge edge = mesh.m_edges[mesh.m_nodesEdges[n][nn]];
            std::size_t neighbour = edge.first == n ? edge.second : edge.first;
            m_nodesNodes[n][nn] = neighbour;
        }
    }

    // computes the number of nodes for each face
    computeFacesNumEdges(mesh);

    // before iteration
    if (m_smoothorarea != 1)
    {
        m_ww2x.resize(mesh.m_nodes.size(), std::vector<double>(m_maximumNumConnectedNodes, 0.0));
        m_ww2y.resize(mesh.m_nodes.size(), std::vector<double>(m_maximumNumConnectedNodes, 0.0));
        // TODO: calculate volume weights
    }

    m_mumax = (1.0 - m_smoothorarea) * 0.5;
    m_mu = std::min(1e-2, m_mumax);
    m_orthogonalCoordinates.resize(mesh.m_nodes.size());

    // in this case the nearest point is the point itself
    m_nearestPoints.resize(mesh.m_nodes.size());
    std::iota(m_nearestPoints.begin(), m_nearestPoints.end(), 0);

    // back-up original nodes, for projection on original mesh boundary
    m_originalNodes = mesh.m_nodes;
    m_orthogonalCoordinates = mesh.m_nodes;

    // algorithm settings
    m_orthogonalizationToSmoothingFactor = orthogonalizationParametersNative.OrthogonalizationToSmoothingFactor;
    m_orthogonalizationToSmoothingFactorBoundary = orthogonalizationParametersNative.OrthogonalizationToSmoothingFactorBoundary;
    m_smoothorarea = orthogonalizationParametersNative.Smoothorarea;

    return true;
}

bool GridGeom::Orthogonalization::Iterate(Mesh& mesh)
{
    bool state = true;

    for (std::size_t outerIter = 0; outerIter < m_orthogonalizationOuterIterations; outerIter++)
    {
        if (state)
        {
            state = PrapareOuterIteration(mesh);
        }
        for (std::size_t boundaryIter = 0; boundaryIter < m_orthogonalizationBoundaryIterations; boundaryIter++)
        {
            for (std::size_t innerIter = 0; innerIter < m_orthogonalizationInnerIterations; innerIter++)
            {
                if (state)
                {
                    state = InnerIteration(mesh);
                }
                    
            } // inner iter, inner iteration
        } // boundary iter

        //update mu
        if (state)
        {
            state = FinalizeOuterIteration(mesh);
        }
            
    }// outer iter

    deallocateCaches();

    return true;
}


bool GridGeom::Orthogonalization::PrapareOuterIteration(const Mesh& mesh) 
{

    bool state = true;

    //compute aspect ratios
    if (state)
    {
        state = aspectRatio(mesh);
    }

    //compute weights orthogonalizer
    if (state)
    {
        state = computeWeightsOrthogonalizer(mesh);
    }

    //compute weights operators for smoother
    if (state)
    {
        state = computeSmootherOperators(mesh);
    }

    //compute weights smoother
    if (state)
    {
        state = computeWeightsSmoother(mesh);
    }
    
    //allocate orthogonalization caches
    if (state)
    {
        state = allocateCaches(mesh);
    }

    if (state)
    {
        state = computeIncrements(mesh);
    }
    return state;
}


bool GridGeom::Orthogonalization::allocateCaches(const Mesh& mesh) 
{
    bool state = true;
    // reallocate caches
    if (state && m_cacheSize == 0)
    {
        m_increments.resize(mesh.m_nodes.size() * 2);
        m_rightHandSideCache.resize(mesh.m_nodes.size() * 2);
        m_startCacheIndex.resize(mesh.m_nodes.size());
        m_endCacheIndex.resize(mesh.m_nodes.size());

        for (int n = 0; n < mesh.m_nodes.size(); n++)
        {
            m_startCacheIndex[n] = m_cacheSize;
            m_cacheSize += std::max(mesh.m_nodesNumEdges[n] + 1, m_numConnectedNodes[n]);
            m_endCacheIndex[n] = m_cacheSize;
        }
        m_k1.resize(m_cacheSize);
        m_wwx.resize(m_cacheSize);
        m_wwy.resize(m_cacheSize);
    }
    return state;
}


bool GridGeom::Orthogonalization::deallocateCaches() 
{
    m_increments.resize(0);
    m_rightHandSideCache.resize(0);
    m_startCacheIndex.resize(0);
    m_endCacheIndex.resize(0);
    m_k1.resize(0);
    m_wwx.resize(0);
    m_wwy.resize(0);
    m_cacheSize = 0;

    return true;
}

bool GridGeom::Orthogonalization::FinalizeOuterIteration(Mesh& mesh)
{
    m_mu = std::min(2.0 * m_mu, m_mumax);

    //compute new faces circumcenters
    if (m_keepCircumcentersAndMassCenters != true)
    {
        mesh.FaceCircumcenters(1.0);
        mesh.FacesAreasAndMassCenters();
    }

    return true;
}

bool GridGeom::Orthogonalization::computeIncrements(const Mesh& mesh)
{
    double max_aptf = std::max(m_orthogonalizationToSmoothingFactorBoundary, m_orthogonalizationToSmoothingFactor);
    
#pragma omp parallel for
	for (int n = 0; n < mesh.m_nodes.size(); n++)
    {
        int firstCacheIndex = n * 2;
        double increments[2]{ 0.0, 0.0 };
        if ((mesh.m_nodesTypes[n] != 1 && mesh.m_nodesTypes[n] != 2) || mesh.m_nodesNumEdges[n] < 2)
        {
            continue;
        }
        if (m_keepCircumcentersAndMassCenters != false && (mesh.m_nodesNumEdges[n] != 3 || mesh.m_nodesNumEdges[n] != 1))
        {
            continue;
        }

        double atpfLoc = mesh.m_nodesTypes[n] == 2 ? max_aptf : m_orthogonalizationToSmoothingFactor;
        double atpf1Loc = 1.0 - atpfLoc;

        double mumat = m_mu;
        if (!m_ww2x.empty() && m_ww2x[n][0] != 0.0 && m_ww2y[n][0] != 0.0)
        {
            mumat = m_mu * m_ww2Global[n][0] / std::max(m_ww2x[n][0], m_ww2y[n][0]);
        }

        increments[0] = 0.0;
        increments[1] = 0.0;

        int maxnn = m_endCacheIndex[n] - m_startCacheIndex[n];
        for (int nn = 1, cacheIndex = m_startCacheIndex[n]; nn < maxnn; nn++, cacheIndex++)
        {
            double wwx = 0.0;
            double wwy = 0.0;
            // Smoother
            if (atpf1Loc > 0.0 && mesh.m_nodesTypes[n] == 1)
            {
                if(!m_ww2x.empty())
                {
                    wwx = atpf1Loc * (mumat * m_ww2x[n][nn] + m_ww2Global[n][nn]);
                    wwy = atpf1Loc * (mumat * m_ww2y[n][nn] + m_ww2Global[n][nn]);
                }
                else
                {
                    wwx = atpf1Loc * m_ww2Global[n][nn];
                    wwy = atpf1Loc * m_ww2Global[n][nn];
                }
            }
            
            // Orthogonalizer
            if (nn < mesh.m_nodesNumEdges[n] + 1)
            {
                wwx += atpfLoc * m_weights[n][nn - 1];
                wwy += atpfLoc * m_weights[n][nn - 1];
                m_k1[cacheIndex] = m_nodesNodes[n][nn - 1];
            }
            else
            {
                m_k1[cacheIndex] = m_connectedNodes[n][nn];
            }
            m_wwx[cacheIndex] = wwx;
            m_wwy[cacheIndex] = wwy;
            orthogonalizationComputeIncrements(wwx, wwy, increments, mesh.m_projection);
        }

        m_increments[firstCacheIndex] = increments[0];
        m_increments[firstCacheIndex + 1] = increments[1];

        m_rightHandSideCache[firstCacheIndex] = atpfLoc * m_rightHandSide[n][0];
        m_rightHandSideCache[firstCacheIndex + 1] = atpfLoc * m_rightHandSide[n][1];
	}
	
	return true;
}


bool GridGeom::Orthogonalization::InnerIteration(Mesh& mesh)
{
#pragma omp parallel for
	for (int n = 0; n < mesh.m_nodes.size(); n++)
    {
        double dx0 = 0.0;
        double dy0 = 0.0;
        int firstCacheIndex = n * 2;
        int maxnn = m_endCacheIndex[n] - m_startCacheIndex[n];
        for (int nn = 1, cacheIndex = m_startCacheIndex[n]; nn < maxnn; nn++, cacheIndex++)
        {
            orthogonalizationComputeDeltasDxDy(m_k1[cacheIndex], n, m_wwx[cacheIndex], m_wwy[cacheIndex], mesh.m_nodes, dx0, dy0, mesh.m_projection);
        }
        
        if (std::abs(m_increments[firstCacheIndex]) > 1e-8 && std::abs(m_increments[firstCacheIndex +1]) > 1e-8)
        {
            dx0 = (dx0 + m_rightHandSideCache[firstCacheIndex]) / m_increments[firstCacheIndex];
            dy0 = (dy0 + m_rightHandSideCache[firstCacheIndex +1]) / m_increments[firstCacheIndex +1];
        }

        orthogonalizationComputeCoordinates(dx0, dy0, mesh.m_nodes[n], m_orthogonalCoordinates[n], mesh.m_projection);
    }
	
    // update mesh node coordinates
    mesh.m_nodes = m_orthogonalCoordinates;

    // project on the original net boundary
    projectOnBoundary(mesh);

    // project on land boundary
    //m_landBoundaries->SnapMeshToLandBoundaries(mesh);

    return true;
}

bool GridGeom::Orthogonalization::SnapToLandBoundary(Mesh& mesh, const LandBoundaries& landBoundaries)
{
    bool successful = false;

    return successful;
}


bool GridGeom::Orthogonalization::projectOnBoundary(Mesh& mesh)
{
    Point firstPoint;
    Point secondPoint;
    Point thirdPoint;
    Point normalSecondPoint;
    Point normalThirdPoint;
    int leftNode;
    int rightNode;

    for (std::size_t n = 0; n < mesh.m_nodes.size(); n++)
    {
        int nearestPointIndex = m_nearestPoints[n];
        if (mesh.m_nodesTypes[n] == 2 && mesh.m_nodesNumEdges[n] > 0 && mesh.m_nodesNumEdges[nearestPointIndex] > 0)
        {
            firstPoint = mesh.m_nodes[n];
            int numEdges = mesh.m_nodesNumEdges[nearestPointIndex];
            int numNodes = 0;
            for (std::size_t nn = 0; nn < numEdges; nn++)
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

            //Project the moved boundary point back onto the closest ORIGINAL edge(netlink) (either between 0 and 2 or 0 and 3)
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


bool GridGeom::Orthogonalization::computeWeightsSmoother(const Mesh& mesh)
{
    std::vector<std::vector<double>> J(mesh.m_nodes.size(), std::vector<double>(4, 0)); //Jacobian
    std::vector<std::vector<double>> Ginv(mesh.m_nodes.size(), std::vector<double>(4, 0)); //mesh monitor matrices

    for (std::size_t n = 0; n < mesh.m_nodes.size(); n++)
    {
        if (mesh.m_nodesTypes[n] != 1 && mesh.m_nodesTypes[n] != 2 && mesh.m_nodesTypes[n] != 4) continue;
        int currentTopology = m_nodeTopologyMapping[n];
        orthogonalizationComputeJacobian(n, m_Jxi[currentTopology], m_Jeta[currentTopology], m_topologyConnectedNodes[currentTopology], m_numTopologyNodes[currentTopology], mesh.m_nodes, J[n], mesh.m_projection);
    }

    // TODO: Account for samples: call orthonet_comp_Ginv(u, ops, J, Ginv)
    for (std::size_t n = 0; n < mesh.m_nodes.size(); n++)
    {
        Ginv[n][0] = 1.0;
        Ginv[n][1] = 0.0;
        Ginv[n][2] = 0.0;
        Ginv[n][3] = 1.0;
    }

    m_ww2Global.resize(mesh.m_nodes.size(), std::vector<double>(m_maximumNumConnectedNodes, 0));
    std::vector<double> a1(2);
    std::vector<double> a2(2);

    // matrices for dicretization
    std::vector<double> DGinvDxi(4, 0.0);
    std::vector<double> DGinvDeta(4, 0.0);
    std::vector<double> currentGinv(4, 0.0);
    std::vector<double> GxiByDivxi(m_maximumNumConnectedNodes, 0.0);
    std::vector<double> GxiByDiveta(m_maximumNumConnectedNodes, 0.0);
    std::vector<double> GetaByDivxi(m_maximumNumConnectedNodes, 0.0);
    std::vector<double> GetaByDiveta(m_maximumNumConnectedNodes, 0.0);
    for (std::size_t n = 0; n < mesh.m_nodes.size(); n++)
    {

        if (mesh.m_nodesNumEdges[n] < 2) continue;

        // Internal nodes and boundary nodes
        if (mesh.m_nodesTypes[n] == 1 || mesh.m_nodesTypes[n] == 2)
        {
            int currentTopology = m_nodeTopologyMapping[n];

            orthogonalizationComputeJacobian(n, m_Jxi[currentTopology], m_Jeta[currentTopology], m_topologyConnectedNodes[currentTopology], m_numTopologyNodes[currentTopology], mesh.m_nodes, J[n], mesh.m_projection);

            //compute the contravariant base vectors
            double determinant = J[n][0] * J[n][3] - J[n][3] * J[n][1];
            if (determinant == 0.0)
            {
                continue;
            }

            a1[0] = J[n][3] / determinant;
            a1[1] = -J[n][2] / determinant;
            a2[0] = -J[n][1] / determinant;
            a2[1] = J[n][0] / determinant;

            //m << J[n][0], J[n][3], J[n][2], J[n][4];
            //Eigen::JacobiSVD<Eigen::MatrixXf> svd(m, Eigen::ComputeThinU | Eigen::ComputeThinV);

            std::fill(DGinvDxi.begin(), DGinvDxi.end(), 0.0);
            std::fill(DGinvDeta.begin(), DGinvDeta.end(), 0.0);
            for (int i = 0; i < m_numTopologyNodes[currentTopology]; i++)
            {
                DGinvDxi[0] += Ginv[m_topologyConnectedNodes[currentTopology][i]][0] * m_Jxi[currentTopology][i];
                DGinvDxi[1] += Ginv[m_topologyConnectedNodes[currentTopology][i]][1] * m_Jxi[currentTopology][i];
                DGinvDxi[2] += Ginv[m_topologyConnectedNodes[currentTopology][i]][2] * m_Jxi[currentTopology][i];
                DGinvDxi[3] += Ginv[m_topologyConnectedNodes[currentTopology][i]][3] * m_Jxi[currentTopology][i];

                DGinvDeta[0] += Ginv[m_topologyConnectedNodes[currentTopology][i]][0] * m_Jeta[currentTopology][i];
                DGinvDeta[1] += Ginv[m_topologyConnectedNodes[currentTopology][i]][1] * m_Jeta[currentTopology][i];
                DGinvDeta[2] += Ginv[m_topologyConnectedNodes[currentTopology][i]][2] * m_Jeta[currentTopology][i];
                DGinvDeta[3] += Ginv[m_topologyConnectedNodes[currentTopology][i]][3] * m_Jeta[currentTopology][i];
            }

            // compute current Ginv
            currentGinv[0] = Ginv[n][0];
            currentGinv[1] = Ginv[n][1];
            currentGinv[2] = Ginv[n][2];
            currentGinv[3] = Ginv[n][3];
            double currentGinvFactor = matrixNorm(a1, a2, currentGinv);

            // compute small matrix operations
            std::fill(GxiByDivxi.begin(), GxiByDivxi.end(), 0.0);
            std::fill(GxiByDiveta.begin(), GxiByDiveta.end(), 0.0);
            std::fill(GetaByDivxi.begin(), GetaByDivxi.end(), 0.0);
            std::fill(GetaByDiveta.begin(), GetaByDiveta.end(), 0.0);
            for (int i = 0; i < m_numTopologyNodes[currentTopology]; i++)
            {
                for (int j = 0; j < m_Divxi[currentTopology].size(); j++)
                {
                    GxiByDivxi[i] += m_Gxi[currentTopology][j][i] * m_Divxi[currentTopology][j];
                    GxiByDiveta[i] += m_Gxi[currentTopology][j][i] * m_Diveta[currentTopology][j];
                    GetaByDivxi[i] += m_Geta[currentTopology][j][i] * m_Divxi[currentTopology][j];
                    GetaByDiveta[i] += m_Geta[currentTopology][j][i] * m_Diveta[currentTopology][j];
                }
            }
            for (int i = 0; i < m_numTopologyNodes[currentTopology]; i++)
            {
                m_ww2Global[n][i] -= matrixNorm(a1, a1, DGinvDxi) * m_Jxi[currentTopology][i] +
                    matrixNorm(a1, a2, DGinvDeta) * m_Jxi[currentTopology][i] +
                    matrixNorm(a2, a1, DGinvDxi) * m_Jeta[currentTopology][i] +
                    matrixNorm(a2, a2, DGinvDeta) * m_Jeta[currentTopology][i];
                m_ww2Global[n][i] += (matrixNorm(a1, a1, currentGinv) * GxiByDivxi[i] +
                    matrixNorm(a1, a2, currentGinv) * GxiByDiveta[i] +
                    matrixNorm(a2, a1, currentGinv) * GetaByDivxi[i] +
                    matrixNorm(a2, a2, currentGinv) * GetaByDiveta[i]);
            }

            double alpha = 0.0;
            for (int i = 1; i < m_numTopologyNodes[currentTopology]; i++)
            {
                alpha = std::max(alpha, -m_ww2Global[n][i]) / std::max(1.0, m_ww2[currentTopology][i]);
            }

            double sumValues = 0.0;
            for (int i = 1; i < m_numTopologyNodes[currentTopology]; i++)
            {
                m_ww2Global[n][i] = m_ww2Global[n][i] + alpha * std::max(1.0, m_ww2[currentTopology][i]);
                sumValues += m_ww2Global[n][i];
            }
            m_ww2Global[n][0] = -sumValues;
            for (int i = 0; i < m_numTopologyNodes[currentTopology]; i++)
            {
                m_ww2Global[n][i] = -m_ww2Global[n][i] / (-sumValues + 1e-8);
            }
        }

    }
    return true;
}


bool GridGeom::Orthogonalization::computeSmootherOperators(const Mesh& mesh)
{
    //allocate small administration arrays only once
    std::vector<int> sharedFaces(maximumNumberOfEdgesPerNode, -1); //icell
    std::vector<std::size_t> connectedNodes(maximumNumberOfConnectedNodes, 0); //adm%kk2
    std::vector<std::vector<std::size_t>> faceNodeMapping(maximumNumberOfConnectedNodes, std::vector<std::size_t>(maximumNumberOfNodesPerFace, 0));//kkc
    std::vector<double> xi(maximumNumberOfConnectedNodes, 0.0);
    std::vector<double> eta(maximumNumberOfConnectedNodes, 0.0);

    m_numConnectedNodes.resize(mesh.m_nodes.size(), 0.0);
    m_connectedNodes.resize(mesh.m_nodes.size(), std::vector<std::size_t>(maximumNumberOfConnectedNodes));
    bool state = initializeTopologies(mesh);
    for (std::size_t n = 0; n < mesh.m_nodes.size(); n++)
    {
        int numSharedFaces = 0;
        int numConnectedNodes = 0;
        if (state)
        {
            state = orthogonalizationAdministration(mesh, n, sharedFaces, numSharedFaces, connectedNodes, numConnectedNodes, faceNodeMapping);
        }

        if (state)
        {
            state = computeXiEta(mesh, n, sharedFaces, numSharedFaces, connectedNodes, numConnectedNodes, faceNodeMapping, xi, eta);
        }
            
        if (state)
        {
            state = saveTopology(n, sharedFaces, numSharedFaces, connectedNodes, numConnectedNodes, faceNodeMapping, xi, eta);
        }
        
        if (state)
        {
            m_maximumNumConnectedNodes = std::max(m_maximumNumConnectedNodes, numConnectedNodes);
            m_maximumNumSharedFaces = std::max(m_maximumNumSharedFaces, numSharedFaces);
        }
            
    }

    if (!state)
    {
        return false;
    }

    // allocate local operators for unique topologies
    m_Az.resize(m_numTopologies);
    m_Gxi.resize(m_numTopologies);
    m_Geta.resize(m_numTopologies);
    m_Divxi.resize(m_numTopologies);
    m_Diveta.resize(m_numTopologies);
    m_Jxi.resize(m_numTopologies);
    m_Jeta.resize(m_numTopologies);
    m_ww2.resize(m_numTopologies);

    // allocate caches
    m_boundaryEdges.resize(2, -1);
    m_leftXFaceCenter.resize(maximumNumberOfEdgesPerNode, 0.0);
    m_leftYFaceCenter.resize(maximumNumberOfEdgesPerNode, 0.0);
    m_rightXFaceCenter.resize(maximumNumberOfEdgesPerNode, 0.0);
    m_rightYFaceCenter.resize(maximumNumberOfEdgesPerNode, 0.0);
    m_xis.resize(maximumNumberOfEdgesPerNode, 0.0);
    m_etas.resize(maximumNumberOfEdgesPerNode, 0.0);

    std::vector<bool> isNewTopology(m_numTopologies, true);
    for (std::size_t n = 0; n < mesh.m_nodes.size(); n++)
    {
        int currentTopology = m_nodeTopologyMapping[n];

        if (isNewTopology[currentTopology])
        {
            isNewTopology[currentTopology] = false;
            // Compute node operators
            if (state)
            {
                state = allocateNodeOperators(currentTopology);
            }
            if (state)
            {
                state = computeOperatorsNode(mesh, n,
                    m_numTopologyNodes[currentTopology], m_topologyConnectedNodes[currentTopology],
                    m_numTopologyFaces[currentTopology], m_topologySharedFaces[currentTopology],
                    m_topologyXi[currentTopology], m_topologyEta[currentTopology],
                    m_topologyFaceNodeMapping[currentTopology]);
            } 
        }
    }

    //have side effects only for spherical coordinates
    if (state)
    {
        state = orthogonalizationComputeLocalCoordinates(mesh.m_nodesNumEdges, m_numConnectedNodes, m_localCoordinates, mesh.m_projection);
    }

    return state;
}


bool GridGeom::Orthogonalization::computeOperatorsNode(const Mesh& mesh, const int currentNode,
    const std::size_t& numConnectedNodes, const std::vector<std::size_t>& connectedNodes,
    const std::size_t& numSharedFaces, const std::vector<int>& sharedFaces,
    const std::vector<double>& xi, const std::vector<double>& eta,
    const std::vector<std::vector<std::size_t>>& faceNodeMapping
)
{
    // the current topology index
    int topologyIndex = m_nodeTopologyMapping[currentNode];

    for (int f = 0; f < numSharedFaces; f++)
    {
        if (sharedFaces[f] < 0 || mesh.m_nodesTypes[currentNode] == 3) continue;

        int edgeLeft = f + 1;
        int edgeRight = edgeLeft + 1;
        if (edgeRight > numSharedFaces)
        {
            edgeRight -= numSharedFaces;
        }

        double edgeLeftSquaredDistance = std::sqrt(xi[edgeLeft] * xi[edgeLeft] + eta[edgeLeft] * eta[edgeLeft] + 1e-16);
        double edgeRightSquaredDistance = std::sqrt(xi[edgeRight] * xi[edgeRight] + eta[edgeRight] * eta[edgeRight] + 1e-16);
        double cDPhi = (xi[edgeLeft] * xi[edgeRight] + eta[edgeLeft] * eta[edgeRight]) / (edgeLeftSquaredDistance * edgeRightSquaredDistance);

        int numFaceNodes = m_faceNumNodes[sharedFaces[f]];


        if (numFaceNodes == 3)
        {
            // determine the index of the current stencil node
            int nodeIndex = FindIndex(mesh.m_facesNodes[sharedFaces[f]], std::size_t(currentNode));

            int nodeLeft = nodeIndex - 1; if (nodeLeft < 0)nodeLeft += numFaceNodes;
            int nodeRight = nodeIndex + 1; if (nodeRight >= numFaceNodes) nodeRight -= numFaceNodes;
            double alpha = 1.0 / (1.0 - cDPhi * cDPhi + 1e-8);
            double alphaLeft = 0.5 * (1.0 - edgeLeftSquaredDistance / edgeRightSquaredDistance * cDPhi) * alpha;
            double alphaRight = 0.5 * (1.0 - edgeRightSquaredDistance / edgeLeftSquaredDistance * cDPhi) * alpha;

            m_Az[topologyIndex][f][faceNodeMapping[f][nodeIndex]] = 1.0 - (alphaLeft + alphaRight);
            m_Az[topologyIndex][f][faceNodeMapping[f][nodeLeft]] = alphaLeft;
            m_Az[topologyIndex][f][faceNodeMapping[f][nodeRight]] = alphaRight;
        }
        else
        {
            for (int i = 0; i < faceNodeMapping[f].size(); i++)
            {
                m_Az[topologyIndex][f][faceNodeMapping[f][i]] = 1.0 / double(numFaceNodes);
            }
        }
    }

    // initialize caches
    std::fill(m_boundaryEdges.begin(), m_boundaryEdges.end(), -1);
    std::fill(m_leftXFaceCenter.begin(), m_leftXFaceCenter.end(), 0.0);
    std::fill(m_leftYFaceCenter.begin(), m_leftYFaceCenter.end(), 0.0);
    std::fill(m_rightXFaceCenter.begin(), m_rightXFaceCenter.end(), 0.0);
    std::fill(m_rightYFaceCenter.begin(), m_rightYFaceCenter.end(), 0.0);
    std::fill(m_xis.begin(), m_xis.end(), 0.0);
    std::fill(m_etas.begin(), m_etas.end(), 0.0);

    int faceRightIndex = 0;
    int faceLeftIndex = 0;
    double xiBoundary = 0.0;
    double etaBoundary = 0.0;

    for (int f = 0; f < numSharedFaces; f++)
    {
        std::size_t edgeIndex = mesh.m_nodesEdges[currentNode][f];
        int otherNode = mesh.m_edges[edgeIndex].first + mesh.m_edges[edgeIndex].second - currentNode;
        int leftFace = mesh.m_edgesFaces[edgeIndex][0];
        faceLeftIndex = FindIndex(sharedFaces, leftFace);

        // face not found, this happens when the cell is outside of the polygon
        if (sharedFaces[faceLeftIndex] != leftFace)
        {
            return false;
        }

        //by construction
        double xiOne = xi[f + 1];
        double etaOne = eta[f + 1];

        double leftRightSwap = 1.0;
        double leftXi = 0.0;
        double leftEta = 0.0;
        double rightXi = 0.0;
        double rightEta = 0.0;
        double alpha_x = 0.0;
        if (mesh.m_edgesNumFaces[edgeIndex] == 1)
        {
            if (m_boundaryEdges[0] < 0)
            {
                m_boundaryEdges[0] = f;
            }
            else
            {
                m_boundaryEdges[1] = f;
            }

            // find the boundary cell in the icell array assume boundary at the right
            // swap Left and Right if the boundary is at the left with I_SWAP_LR
            if (f != faceLeftIndex) leftRightSwap = -1.0;

            for (int i = 0; i < numConnectedNodes; i++)
            {
                leftXi += xi[i] * m_Az[topologyIndex][faceLeftIndex][i];
                leftEta += eta[i] * m_Az[topologyIndex][faceLeftIndex][i];
                m_leftXFaceCenter[f] += mesh.m_nodes[connectedNodes[i]].x * m_Az[topologyIndex][faceLeftIndex][i];
                m_leftYFaceCenter[f] += mesh.m_nodes[connectedNodes[i]].y * m_Az[topologyIndex][faceLeftIndex][i];
            }


            double alpha = leftXi * xiOne + leftEta * etaOne;
            alpha = alpha / (xiOne * xiOne + etaOne * etaOne);

            alpha_x = alpha;
            xiBoundary = alpha * xiOne;
            etaBoundary = alpha * etaOne;

            rightXi = 2.0 * xiBoundary - leftXi;
            rightEta = 2.0 * etaBoundary - leftEta;

            double xBc = (1.0 - alpha) * mesh.m_nodes[currentNode].x + alpha * mesh.m_nodes[otherNode].x;
            double yBc = (1.0 - alpha) * mesh.m_nodes[currentNode].y + alpha * mesh.m_nodes[otherNode].y;
            m_leftYFaceCenter[f] = 2.0 * xBc - m_leftXFaceCenter[f];
            m_rightYFaceCenter[f] = 2.0 * yBc - m_leftYFaceCenter[f];
        }
        else
        {
            faceLeftIndex = f;
            faceRightIndex = faceLeftIndex - 1; if (faceRightIndex < 0)faceRightIndex += numSharedFaces;

            if (faceRightIndex < 0) continue;

            int faceLeft = sharedFaces[faceLeftIndex];
            int faceRight = sharedFaces[faceRightIndex];

            if ((faceLeft != mesh.m_edgesFaces[edgeIndex][0] && faceLeft != mesh.m_edgesFaces[edgeIndex][1]) ||
                (faceRight != mesh.m_edgesFaces[edgeIndex][0] && faceRight != mesh.m_edgesFaces[edgeIndex][1]))
            {
                return false;
            }

            int numNodesLeftFace = m_faceNumNodes[faceLeft];
            int numNodesRightFace = m_faceNumNodes[faceRight];

            for (int i = 0; i < numConnectedNodes; i++)
            {
                leftXi += xi[i] * m_Az[topologyIndex][faceLeftIndex][i];
                leftEta += eta[i] * m_Az[topologyIndex][faceLeftIndex][i];
                rightXi += xi[i] * m_Az[topologyIndex][faceRightIndex][i];
                rightEta += eta[i] * m_Az[topologyIndex][faceRightIndex][i];

                m_leftXFaceCenter[f] += mesh.m_nodes[connectedNodes[i]].x * m_Az[topologyIndex][faceLeftIndex][i];
                m_leftYFaceCenter[f] += mesh.m_nodes[connectedNodes[i]].y * m_Az[topologyIndex][faceLeftIndex][i];
                m_leftYFaceCenter[f] += mesh.m_nodes[connectedNodes[i]].x * m_Az[topologyIndex][faceRightIndex][i];
                m_rightYFaceCenter[f] += mesh.m_nodes[connectedNodes[i]].y * m_Az[topologyIndex][faceRightIndex][i];
            }
        }

        m_xis[f] = 0.5 * (leftXi + rightXi);
        m_etas[f] = 0.5 * (leftEta + rightEta);

        double exiLR = rightXi - leftXi;
        double eetaLR = rightEta - leftEta;
        double exi01 = xiOne;
        double eeta01 = etaOne;

        double fac = 1.0 / std::abs(exi01 * eetaLR - eeta01 * exiLR + 1e-16);
        double facxi1 = -eetaLR * fac * leftRightSwap;

        double facxi0 = -facxi1;
        double faceta1 = exiLR * fac * leftRightSwap;
        double faceta0 = -faceta1;
        double facxiR = eeta01 * fac * leftRightSwap;
        double facxiL = -facxiR;
        double facetaR = -exi01 * fac * leftRightSwap;
        double facetaL = -facetaR;

        //boundary link
        if (mesh.m_edgesNumFaces[edgeIndex] == 1)
        {
            facxi1 = facxi1 - facxiL * 2.0 * alpha_x;
            facxi0 = facxi0 - facxiL * 2.0 * (1.0 - alpha_x);
            facxiL = facxiL + facxiL;
            //note that facxiR does not exist
            faceta1 = faceta1 - facetaL * 2.0 * alpha_x;
            faceta0 = faceta0 - facetaL * 2.0 * (1.0 - alpha_x);
            facetaL = facetaL + facetaL;
        }

        int node1 = f + 1;
        int node0 = 0;
        for (int i = 0; i < numConnectedNodes; i++)
        {
            m_Gxi[topologyIndex][f][i] = facxiL * m_Az[topologyIndex][faceLeftIndex][i];
            m_Geta[topologyIndex][f][i] = facetaL * m_Az[topologyIndex][faceLeftIndex][i];
            if (mesh.m_edgesNumFaces[edgeIndex] == 2)
            {
                m_Gxi[topologyIndex][f][i] = m_Gxi[topologyIndex][f][i] + facxiR * m_Az[topologyIndex][faceRightIndex][i];
                m_Geta[topologyIndex][f][i] = m_Geta[topologyIndex][f][i] + facetaR * m_Az[topologyIndex][faceRightIndex][i];
            }
        }


        m_Gxi[topologyIndex][f][node1] = m_Gxi[topologyIndex][f][node1] + facxi1;
        m_Geta[topologyIndex][f][node1] = m_Geta[topologyIndex][f][node1] + faceta1;

        m_Gxi[topologyIndex][f][node0] = m_Gxi[topologyIndex][f][node0] + facxi0;
        m_Geta[topologyIndex][f][node0] = m_Geta[topologyIndex][f][node0] + faceta0;

        //fill the node - based gradient matrix
        m_Divxi[topologyIndex][f] = -eetaLR * leftRightSwap;
        m_Diveta[topologyIndex][f] = exiLR * leftRightSwap;

        // boundary link
        if (mesh.m_edgesNumFaces[edgeIndex] == 1)
        {
            m_Divxi[topologyIndex][f] = 0.5 * m_Divxi[topologyIndex][f] + etaBoundary * leftRightSwap;
            m_Diveta[topologyIndex][f] = 0.5 * m_Diveta[topologyIndex][f] - xiBoundary * leftRightSwap;
        }


    }

    double volxi = 0.0;
    for (int i = 0; i < mesh.m_nodesNumEdges[currentNode]; i++)
    {
        volxi += 0.5 * (m_Divxi[topologyIndex][i] * m_xis[i] + m_Diveta[topologyIndex][i] * m_etas[i]);
    }
    if (volxi == 0.0)volxi = 1.0;

    for (int i = 0; i < mesh.m_nodesNumEdges[currentNode]; i++)
    {
        m_Divxi[topologyIndex][i] = m_Divxi[topologyIndex][i] / volxi;
        m_Diveta[topologyIndex][i] = m_Diveta[topologyIndex][i] / volxi;
    }

    //compute the node-to-node gradients
    for (int f = 0; f < numSharedFaces; f++)
    {
        // internal edge
        if (mesh.m_edgesNumFaces[mesh.m_nodesEdges[currentNode][f]] == 2)
        {
            int rightNode = f - 1;
            if (rightNode < 0)
            {
                rightNode += mesh.m_nodesNumEdges[currentNode];
            }
            for (int i = 0; i < numConnectedNodes; i++)
            {
                m_Jxi[topologyIndex][i] += m_Divxi[topologyIndex][f] * 0.5 * (m_Az[topologyIndex][f][i] + m_Az[topologyIndex][rightNode][i]);
                m_Jeta[topologyIndex][i] += m_Diveta[topologyIndex][f] * 0.5 * (m_Az[topologyIndex][f][i] + m_Az[topologyIndex][rightNode][i]);
            }
        }
        else
        {
            m_Jxi[topologyIndex][0] = m_Jxi[topologyIndex][0] + m_Divxi[topologyIndex][f] * 0.5;
            m_Jxi[topologyIndex][f + 1] = m_Jxi[topologyIndex][f + 1] + m_Divxi[topologyIndex][f] * 0.5;
            m_Jeta[topologyIndex][0] = m_Jeta[topologyIndex][0] + m_Diveta[topologyIndex][f] * 0.5;
            m_Jeta[topologyIndex][f + 1] = m_Jeta[topologyIndex][f + 1] + m_Diveta[topologyIndex][f] * 0.5;
        }
    }

    //compute the weights in the Laplacian smoother
    std::fill(m_ww2[topologyIndex].begin(), m_ww2[topologyIndex].end(), 0.0);
    for (int n = 0; n < mesh.m_nodesNumEdges[currentNode]; n++)
    {
        for (int i = 0; i < numConnectedNodes; i++)
        {
            m_ww2[topologyIndex][i] += m_Divxi[topologyIndex][n] * m_Gxi[topologyIndex][n][i] + m_Diveta[topologyIndex][n] * m_Geta[topologyIndex][n][i];
        }
    }

    return true;
}


bool GridGeom::Orthogonalization::computeXiEta(const Mesh& mesh,
    int currentNode,
    const std::vector<int>& sharedFaces,
    const int& numSharedFaces,
    const std::vector<std::size_t>& connectedNodes,
    const std::size_t& numConnectedNodes,
    const std::vector<std::vector<std::size_t>>& faceNodeMapping,
    std::vector<double>& xi,
    std::vector<double>& eta)
{
    std::fill(xi.begin(), xi.end(), 0.0);
    std::fill(eta.begin(), eta.end(), 0.0);
    // the angles for the squared nodes connected to the stencil nodes, first the ones directly connected, then the others
    std::vector<double> thetaSquare(numConnectedNodes, doubleMissingValue);
    // for each shared face, a bollean indicating if it is squared or not
    std::vector<bool> isSquareFace(numSharedFaces, false);

    int numNonStencilQuad = 0;

    //loop over the connected edges
    for (int f = 0; f < numSharedFaces; f++)
    {
        std::size_t edgeIndex = mesh.m_nodesEdges[currentNode][f];
        std::size_t nextNode = connectedNodes[f + 1]; // the first entry is always the stencil node 
        int faceLeft = mesh.m_edgesFaces[edgeIndex][0];
        int faceRigth = faceLeft;

        if (mesh.m_edgesNumFaces[edgeIndex] == 2)
            faceRigth = mesh.m_edgesFaces[edgeIndex][1];

        //check if it is a rectangular node (not currentNode itself) 
        bool isSquare = true;
        for (int e = 0; e < mesh.m_nodesNumEdges[nextNode]; e++)
        {
            std::size_t edge = mesh.m_nodesEdges[nextNode][e];
            for (int ff = 0; ff < mesh.m_edgesNumFaces[edge]; ff++)
            {
                std::size_t face = mesh.m_edgesFaces[edge][ff];
                if (face != faceLeft && face != faceRigth)
                {
                    isSquare = isSquare && mesh.m_facesNodes[face].size() == 4;
                }
            }
            if (!isSquare)
            {
                break;
            }
        }

        //Compute optimal angle thetaSquare based on the node type
        int leftFaceIndex = f - 1; if (leftFaceIndex < 0) leftFaceIndex = leftFaceIndex + numSharedFaces;
        if (isSquare)
        {
            if (mesh.m_nodesTypes[nextNode] == 1 || mesh.m_nodesTypes[nextNode] == 4)
            {
                // Inner node
                numNonStencilQuad = mesh.m_nodesNumEdges[nextNode] - 2;
                thetaSquare[f + 1] = (2.0 - double(numNonStencilQuad) * 0.5) * M_PI;
            }
            if (mesh.m_nodesTypes[nextNode] == 2)
            {
                // boundary node
                numNonStencilQuad = mesh.m_nodesNumEdges[nextNode] - 1 - mesh.m_edgesNumFaces[edgeIndex];
                thetaSquare[f + 1] = (1.0 - double(numNonStencilQuad) * 0.5) * M_PI;
            }
            if (mesh.m_nodesTypes[nextNode] == 3)
            {
                //corner node
                thetaSquare[f + 1] = 0.5 * M_PI;
            }

            if (sharedFaces[f] > 1)
            {
                if (m_faceNumNodes[sharedFaces[f]] == 4) numNonStencilQuad += 1;
            }
            if (sharedFaces[leftFaceIndex] > 1)
            {
                if (m_faceNumNodes[sharedFaces[leftFaceIndex]] == 4) numNonStencilQuad += 1;
            }
            if (numNonStencilQuad > 3)
            {
                isSquare = false;
            }
        }

        isSquareFace[f] = isSquareFace[f] || isSquare;
        isSquareFace[leftFaceIndex] = isSquareFace[leftFaceIndex] || isSquare;
    }

    for (int f = 0; f < numSharedFaces; f++)
    {
        // boundary face
        if (sharedFaces[f] < 0) continue;

        // non boundary face 
        if (m_faceNumNodes[sharedFaces[f]] == 4)
        {
            for (int n = 0; n < m_faceNumNodes[sharedFaces[f]]; n++)
            {
                if (faceNodeMapping[f][n] <= numSharedFaces) continue;
                thetaSquare[faceNodeMapping[f][n]] = 0.5 * M_PI;
            }
        }
    }

    // Compute internal angle
    int numSquaredTriangles = 0.0;
    int numTriangles = 0;
    double phiSquaredTriangles = 0.0;
    double phiQuads = 0.0;
    double phiTriangles = 0.0;
    double phiTot = 0.0;
    numNonStencilQuad = 0;
    for (int f = 0; f < numSharedFaces; f++)
    {
        // boundary face
        if (sharedFaces[f] < 0) continue;

        int numFaceNodes = m_faceNumNodes[sharedFaces[f]];
        double phi = optimalEdgeAngle(numFaceNodes);

        if (isSquareFace[f] || numFaceNodes == 4)
        {
            int nextNode = f + 2;
            if (nextNode > numSharedFaces)
            {
                nextNode = nextNode - numSharedFaces;
            }
            bool isBoundaryEdge = mesh.m_edgesNumFaces[mesh.m_nodesEdges[currentNode][f]] == 1;
            phi = optimalEdgeAngle(numFaceNodes, thetaSquare[f + 1], thetaSquare[nextNode], isBoundaryEdge);
            if (numFaceNodes == 3)
            {
                numSquaredTriangles += 1;
                phiSquaredTriangles += phi;
            }
            else if (numFaceNodes == 4)
            {
                numNonStencilQuad += 1;
                phiQuads += phi;
            }
        }
        else
        {
            numTriangles += 1;
            phiTriangles += phi;
        }
        phiTot += phi;
    }


    double factor = 1.0;
    if (mesh.m_nodesTypes[currentNode] == 2) factor = 0.5;
    if (mesh.m_nodesTypes[currentNode] == 3) factor = 0.25;
    double mu = 1.0;
    double muSquaredTriangles = 1.0;
    double muTriangles = 1.0;
    double minPhi = 15.0 / 180.0 * M_PI;
    if (numTriangles > 0)
    {
        muTriangles = (factor * 2.0 * M_PI - (phiTot - phiTriangles)) / phiTriangles;
        muTriangles = std::max(muTriangles, double(numTriangles) * minPhi / phiTriangles);
    }
    else if (numSquaredTriangles > 0)
    {
        muSquaredTriangles = std::max(factor * 2.0 * M_PI - (phiTot - phiSquaredTriangles), double(numSquaredTriangles) * minPhi) / phiSquaredTriangles;
    }

    if (phiTot > 1e-18)
    {
        mu = factor * 2.0 * M_PI / (phiTot - (1.0 - muTriangles) * phiTriangles - (1.0 - muSquaredTriangles) * phiSquaredTriangles);
    }
    else if (numSharedFaces > 0)
    {
        //TODO: add logger and cirr(xk(k0), yk(k0), ncolhl)
        std::string message{ "fatal error in computeXiEta: phiTot=0'" };
        m_nodeXErrors.push_back(mesh.m_nodes[currentNode].x);
        m_nodeXErrors.push_back(mesh.m_nodes[currentNode].y);
        return false;
    }

    double phi0 = 0.0;
    double dPhi0 = 0.0;
    double dPhi = 0.0;
    double dTheta = 0.0;
    for (int f = 0; f < numSharedFaces; f++)
    {
        phi0 = phi0 + 0.5 * dPhi;
        if (sharedFaces[f] < 0)
        {
            if (mesh.m_nodesTypes[currentNode] == 2)
            {
                dPhi = M_PI;
            }
            else if (mesh.m_nodesTypes[currentNode] == 3)
            {
                dPhi = 1.5 * M_PI;
            }
            else
            {
                //TODO: add logger and cirr(xk(k0), yk(k0), ncolhl)
                std::string message{ "fatal error in computeXiEta: inappropriate fictitious boundary cell" };
                m_nodeXErrors.push_back(mesh.m_nodes[currentNode].x);
                m_nodeXErrors.push_back(mesh.m_nodes[currentNode].y);
                return false;
            }
            phi0 = phi0 + 0.5 * dPhi;
            continue;
        }

        int numFaceNodes = m_faceNumNodes[sharedFaces[f]];
        if (numFaceNodes > maximumNumberOfEdgesPerNode)
        {
            //TODO: error
            return false;
        }

        dPhi0 = optimalEdgeAngle(numFaceNodes);
        if (isSquareFace[f])
        {
            int nextNode = f + 2; if (nextNode > numSharedFaces) nextNode = nextNode - numSharedFaces;
            bool isBoundaryEdge = mesh.m_edgesNumFaces[mesh.m_nodesEdges[currentNode][f]] == 1;
            dPhi0 = optimalEdgeAngle(numFaceNodes, thetaSquare[f + 1], thetaSquare[nextNode], isBoundaryEdge);
            if (numFaceNodes == 3)
            {
                dPhi0 = muSquaredTriangles * dPhi0;
            }
        }
        else if (numFaceNodes == 3)
        {
            dPhi0 = muTriangles * dPhi0;
        }

        dPhi = mu * dPhi0;
        phi0 = phi0 + 0.5 * dPhi;

        // determine the index of the current stencil node
        int nodeIndex = FindIndex(mesh.m_facesNodes[sharedFaces[f]], std::size_t(currentNode));

        // optimal angle
        dTheta = 2.0 * M_PI / double(numFaceNodes);

        // orientation of the face (necessary for folded cells)
        int previousNode = nodeIndex + 1; if (previousNode >= numFaceNodes) previousNode -= numFaceNodes;
        int nextNode = nodeIndex - 1; if (nextNode < 0) nextNode += numFaceNodes;
        if ((faceNodeMapping[f][nextNode] - faceNodeMapping[f][previousNode]) == -1 ||
            (faceNodeMapping[f][nextNode] - faceNodeMapping[f][previousNode]) == mesh.m_nodesNumEdges[currentNode])
        {
            dTheta = -dTheta;
        }

        double aspectRatio = (1.0 - std::cos(dTheta)) / std::sin(std::abs(dTheta)) * std::tan(0.5 * dPhi);
        double radius = std::cos(0.5 * dPhi) / (1.0 - cos(dTheta));

        for (int n = 0; n < numFaceNodes; n++)
        {
            double theta = dTheta * (n - nodeIndex);
            double xip = radius - radius * std::cos(theta);
            double ethap = -radius * std::sin(theta);

            xi[faceNodeMapping[f][n]] = xip * std::cos(phi0) - aspectRatio * ethap * std::sin(phi0);
            eta[faceNodeMapping[f][n]] = xip * std::sin(phi0) + aspectRatio * ethap * std::cos(phi0);
        }
    }

    return true;
}


bool GridGeom::Orthogonalization::computeFacesNumEdges(const Mesh& mesh)
{
    // Cache the result for later calls.
    m_faceNumNodes.resize(mesh.m_numFaces, false);
    for (int f = 0; f < mesh.m_numFaces; f++)
    {
        m_faceNumNodes[f] = mesh.m_facesNodes[f].size();
    }
    return true;
}


bool GridGeom::Orthogonalization::orthogonalizationAdministration(const Mesh& mesh, const int currentNode, std::vector<int>& sharedFaces, int& numSharedFaces, std::vector<std::size_t>& connectedNodes, int& numConnectedNodes, std::vector<std::vector<std::size_t>>& faceNodeMapping)
{
    for (auto& f : sharedFaces)  f = -1;

    if (mesh.m_nodesNumEdges[currentNode] < 2) return true;

    // 1. For the currentNode, find the shared faces
    int newFaceIndex = -999;
    numSharedFaces = 0;
    for (int e = 0; e < mesh.m_nodesNumEdges[currentNode]; e++)
    {
        std::size_t firstEdge = mesh.m_nodesEdges[currentNode][e];

        int secondEdgeIndex = e + 1;

        if (secondEdgeIndex >= mesh.m_nodesNumEdges[currentNode])
            secondEdgeIndex = 0;

        std::size_t secondEdge = mesh.m_nodesEdges[currentNode][secondEdgeIndex];

        if (mesh.m_edgesNumFaces[firstEdge] < 1 || mesh.m_edgesNumFaces[secondEdge] < 1) continue;

        // find the face that the first and the second edge share
        int firstFaceIndex = std::max(std::min(mesh.m_edgesNumFaces[firstEdge], std::size_t(2)), std::size_t(1)) - 1;
        int secondFaceIndex = std::max(std::min(mesh.m_edgesNumFaces[secondEdge], std::size_t(2)), std::size_t(1)) - 1;

        if (mesh.m_edgesFaces[firstEdge][0] != newFaceIndex &&
            (mesh.m_edgesFaces[firstEdge][0] == mesh.m_edgesFaces[secondEdge][0] || mesh.m_edgesFaces[firstEdge][0] == mesh.m_edgesFaces[secondEdge][secondFaceIndex]))
        {
            newFaceIndex = mesh.m_edgesFaces[firstEdge][0];
        }
        else if (mesh.m_edgesFaces[firstEdge][firstFaceIndex] != newFaceIndex &&
            (mesh.m_edgesFaces[firstEdge][firstFaceIndex] == mesh.m_edgesFaces[secondEdge][0] || mesh.m_edgesFaces[firstEdge][firstFaceIndex] == mesh.m_edgesFaces[secondEdge][secondFaceIndex]))
        {
            newFaceIndex = mesh.m_edgesFaces[firstEdge][firstFaceIndex];
        }
        else
        {
            newFaceIndex = -999;
        }

        //corner face (already found in the first iteration)
        if (mesh.m_nodesNumEdges[currentNode] == 2 && e == 1 && mesh.m_nodesTypes[currentNode] == 3)
        {
            if (sharedFaces[0] == newFaceIndex) newFaceIndex = -999;
        }
        sharedFaces[numSharedFaces] = newFaceIndex;
        numSharedFaces += 1;
    }

    if (numSharedFaces < 1) return true;

    for (auto& v : connectedNodes) v = 0;
    int connectedNodesIndex = 0;
    connectedNodes[connectedNodesIndex] = currentNode;

    // edge connected nodes
    for (int e = 0; e < mesh.m_nodesNumEdges[currentNode]; e++)
    {
        std::size_t edgeIndex = mesh.m_nodesEdges[currentNode][e];
        std::size_t node = mesh.m_edges[edgeIndex].first + mesh.m_edges[edgeIndex].second - currentNode;
        connectedNodesIndex++;
        connectedNodes[connectedNodesIndex] = node;
        if (m_connectedNodes[currentNode].size() < connectedNodesIndex + 1)
        {
            m_connectedNodes[currentNode].resize(connectedNodesIndex);
        }
        m_connectedNodes[currentNode][connectedNodesIndex] = node;
    }

    // for each face store the positions of the its nodes in the connectedNodes (compressed array)
    if (faceNodeMapping.size() < numSharedFaces)
    {
        faceNodeMapping.resize(numSharedFaces, std::vector<std::size_t>(maximumNumberOfNodesPerFace, 0));
    }


    for (int f = 0; f < numSharedFaces; f++)
    {
        int faceIndex = sharedFaces[f];
        if (faceIndex < 0) continue;

        // find the stencil node position  in the current face
        std::size_t faceNodeIndex = 0;
        std::size_t numFaceNodes = m_faceNumNodes[faceIndex];
        for (int n = 0; n < numFaceNodes; n++)
        {
            if (mesh.m_facesNodes[faceIndex][n] == currentNode)
            {
                faceNodeIndex = n;
                break;
            }
        }

        for (int n = 0; n < numFaceNodes; n++)
        {

            if (faceNodeIndex >= numFaceNodes)
            {
                faceNodeIndex -= numFaceNodes;
            }


            int node = mesh.m_facesNodes[faceIndex][faceNodeIndex];

            bool isNewNode = true;
            for (int n = 0; n < connectedNodesIndex + 1; n++)
            {
                if (node == connectedNodes[n])
                {
                    isNewNode = false;
                    faceNodeMapping[f][faceNodeIndex] = n;
                    break;
                }
            }

            if (isNewNode)
            {
                connectedNodesIndex++;
                connectedNodes[connectedNodesIndex] = node;
                faceNodeMapping[f][faceNodeIndex] = connectedNodesIndex;
                if (m_connectedNodes[currentNode].size() < connectedNodesIndex + 1)
                {
                    m_connectedNodes[currentNode].resize(connectedNodesIndex);
                }
                m_connectedNodes[currentNode][connectedNodesIndex] = node;
            }

            //update node index
            faceNodeIndex += 1;
        }
    }

    // compute the number of connected nodes
    numConnectedNodes = connectedNodesIndex + 1;

    //update connected nodes kkc
    m_numConnectedNodes[currentNode] = numConnectedNodes;

    return true;
}


double GridGeom::Orthogonalization::optimalEdgeAngle(int numFaceNodes, double theta1, double theta2, bool isBoundaryEdge)
{
    double angle = M_PI * (1 - 2.0 / double(numFaceNodes));

    if (theta1 != -1 && theta2 != -1 && numFaceNodes == 3)
    {
        angle = 0.25 * M_PI;
        if (theta1 + theta2 == M_PI && !isBoundaryEdge)
        {
            angle = 0.5 * M_PI;
        }
    }
    return angle;
}


bool GridGeom::Orthogonalization::aspectRatio(const Mesh& mesh)
{
    std::vector<std::vector<double>> averageEdgesLength(mesh.m_edges.size(), std::vector<double>(2, doubleMissingValue));
    std::vector<double> averageFlowEdgesLength(mesh.m_edges.size(), doubleMissingValue);
    std::vector<bool> curvilinearGridIndicator(mesh.m_nodes.size(), true);
    std::vector<double> edgesLength(mesh.m_edges.size(), 0.0);

    for (std::size_t e = 0; e < mesh.m_edges.size(); e++)
    {
        std::size_t first = mesh.m_edges[e].first;
        std::size_t second = mesh.m_edges[e].second;

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
            double dinry = innerProductTwoSegments(mesh.m_nodes[first], mesh.m_nodes[second], mesh.m_nodes[first], leftCenter, mesh.m_projection);
            dinry = dinry / std::max(edgeLength * edgeLength, minimumEdgeLength);

            double x0_bc = (1.0 - dinry) * mesh.m_nodes[first].x + dinry * mesh.m_nodes[second].x;
            double y0_bc = (1.0 - dinry) * mesh.m_nodes[first].y + dinry * mesh.m_nodes[second].y;
            rightCenter.x = 2.0 * x0_bc - leftCenter.x;
            rightCenter.y = 2.0 * y0_bc - leftCenter.y;
        }

        averageFlowEdgesLength[e] = Distance(leftCenter, rightCenter, mesh.m_projection);
    }

    // Compute normal length
    for (int f = 0; f < mesh.m_facesNodes.size(); f++)
    {
        std::size_t numberOfFaceNodes = mesh.m_facesNodes[f].size();
        if (numberOfFaceNodes < 3) continue;

        for (int n = 0; n < numberOfFaceNodes; n++)
        {
            if (numberOfFaceNodes != 4) curvilinearGridIndicator[mesh.m_facesNodes[f][n]] = false;
            std::size_t edgeIndex = mesh.m_facesEdges[f][n];

            if (mesh.m_edgesNumFaces[edgeIndex] < 1) continue;

            //get the other links in the right numbering
            //TODO: ask why only 3 are requested, why an average lenght stored in averageEdgesLength is needed?
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
                std::size_t klinkp2 = mesh.m_facesEdges[f][kkp2];
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

    for (std::size_t e = 0; e < mesh.m_edges.size(); e++)
    {
        std::size_t first = mesh.m_edges[e].first;
        std::size_t second = mesh.m_edges[e].second;

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

bool GridGeom::Orthogonalization::computeWeightsOrthogonalizer(const Mesh& mesh)
{
    double localOrthogonalizationToSmoothingFactor = 1.0;
    double localOrthogonalizationToSmoothingFactorSymmetric = 1.0 - localOrthogonalizationToSmoothingFactor;
    constexpr double mu = 1.0;
    Point normal;

    std::fill(m_rightHandSide.begin(), m_rightHandSide.end(), std::vector<double>(2, 0.0));
    for (std::size_t n = 0; n < mesh.m_nodes.size(); n++)
    {
        if (mesh.m_nodesTypes[n] != 1 && mesh.m_nodesTypes[n] != 2)
        {
            continue;
        }

        for (std::size_t nn = 0; nn < mesh.m_nodesNumEdges[n]; nn++)
        {
            std::size_t edgeIndex = mesh.m_nodesEdges[n][nn];
            double aspectRatio = m_aspectRatios[edgeIndex];

            if (aspectRatio != doubleMissingValue)
            {
                m_weights[n][nn] = localOrthogonalizationToSmoothingFactor * aspectRatio + localOrthogonalizationToSmoothingFactorSymmetric * mu;

                if (mesh.m_edgesNumFaces[edgeIndex] == 1)
                {
                    //boundary nodes
                    Point neighbouringNode = mesh.m_nodes[m_nodesNodes[n][nn]];
                    double neighbouringNodeDistance = Distance(neighbouringNode, mesh.m_nodes[n], mesh.m_projection);
                    double aspectRatioByNodeDistance = aspectRatio * neighbouringNodeDistance;

                    std::size_t leftFace = mesh.m_edgesFaces[edgeIndex][0];
                    bool flippedNormal;
                    normalVectorInside(mesh.m_nodes[n], neighbouringNode, mesh.m_facesMassCenters[leftFace], normal, flippedNormal, mesh.m_projection);

                    m_rightHandSide[n][0] += localOrthogonalizationToSmoothingFactor * neighbouringNodeDistance * normal.x / 2.0 +
                        localOrthogonalizationToSmoothingFactorSymmetric * aspectRatioByNodeDistance * normal.x * 0.5 / mu;

                    m_rightHandSide[n][1] += localOrthogonalizationToSmoothingFactor * neighbouringNodeDistance * normal.y / 2.0 +
                        localOrthogonalizationToSmoothingFactorSymmetric * aspectRatioByNodeDistance * normal.y * 0.5 / mu;

                    m_weights[n][nn] = localOrthogonalizationToSmoothingFactor * 0.5 * aspectRatio +
                        localOrthogonalizationToSmoothingFactorSymmetric * 0.5 * mu;
                }

            }
            else
            {
                m_weights[n][nn] = 0.0;
            }
        }

        // normalize
        double factor = std::accumulate(m_weights[n].begin(), m_weights[n].end(), 0.0);
        if (std::abs(factor) > 1e-14)
        {
            factor = 1.0 / factor;
            for (auto& w : m_weights[n]) w = w * factor;
            m_rightHandSide[n][0] = factor * m_rightHandSide[n][0];
            m_rightHandSide[n][1] = factor * m_rightHandSide[n][1];
        }

    }
    return true;
}


double GridGeom::Orthogonalization::matrixNorm(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& matCoefficents)
{
    double norm = (matCoefficents[0] * x[0] + matCoefficents[1] * x[1]) * y[0] + (matCoefficents[2] * x[0] + matCoefficents[3] * x[1]) * y[1];
    return norm;
}


bool GridGeom::Orthogonalization::initializeTopologies(const Mesh& mesh)
{
    // topology 
    m_numTopologies = 0;
    m_nodeTopologyMapping.resize(mesh.m_nodes.size(), -1);
    m_numTopologyNodes.resize(m_topologyInitialSize, -1);
    m_numTopologyFaces.resize(m_topologyInitialSize, -1);
    m_topologyXi.resize(m_topologyInitialSize, std::vector<double>(maximumNumberOfConnectedNodes, 0));
    m_topologyEta.resize(m_topologyInitialSize, std::vector<double>(maximumNumberOfConnectedNodes, 0));
    m_topologySharedFaces.resize(m_topologyInitialSize, std::vector<int>(maximumNumberOfConnectedNodes, -1));
    m_topologyConnectedNodes.resize(m_topologyInitialSize, std::vector<std::size_t>(maximumNumberOfConnectedNodes, -1));
    m_topologyFaceNodeMapping.resize(m_topologyInitialSize, std::vector<std::vector<std::size_t>>(maximumNumberOfConnectedNodes, std::vector<std::size_t>(maximumNumberOfConnectedNodes, -1)));

    return true;
}


bool GridGeom::Orthogonalization::allocateNodeOperators(const int topologyIndex)
{
    int numSharedFaces = m_numTopologyFaces[topologyIndex];
    int numConnectedNodes = m_numTopologyNodes[topologyIndex];
    // will reallocate only if necessary
    m_Az[topologyIndex].resize(numSharedFaces, std::vector<double>(numConnectedNodes, 0.0));
    m_Gxi[topologyIndex].resize(numSharedFaces, std::vector<double>(numConnectedNodes, 0.0));
    m_Geta[topologyIndex].resize(numSharedFaces, std::vector<double>(numConnectedNodes, 0.0));
    m_Divxi[topologyIndex].resize(numSharedFaces);
    m_Diveta[topologyIndex].resize(numSharedFaces);
    m_Jxi[topologyIndex].resize(numConnectedNodes);
    m_Jeta[topologyIndex].resize(numConnectedNodes);
    m_ww2[topologyIndex].resize(numConnectedNodes);

    return true;
}


bool GridGeom::Orthogonalization::saveTopology(int currentNode,
    const std::vector<int>& sharedFaces,
    int numSharedFaces,
    const std::vector<std::size_t>& connectedNodes,
    int numConnectedNodes,
    const std::vector<std::vector<std::size_t>>& faceNodeMapping,
    const std::vector<double>& xi,
    const std::vector<double>& eta)
{
    bool isNewTopology = true;
    for (int topo = 0; topo < m_numTopologies; topo++)
    {
        if (numSharedFaces != m_numTopologyFaces[topo] || numConnectedNodes != m_numTopologyNodes[topo])
        {
            continue;
        }

        isNewTopology = false;
        for (int n = 1; n < numConnectedNodes; n++)
        {
            double thetaLoc = std::atan2(eta[n], xi[n]);
            double thetaTopology = std::atan2(m_topologyEta[topo][n], m_topologyXi[topo][n]);
            if (std::abs(thetaLoc - thetaTopology) > m_thetaTolerance)
            {
                isNewTopology = true;
                break;
            }
        }

        if (!isNewTopology)
        {
            m_nodeTopologyMapping[currentNode] = topo;
            break;
        }
    }

    if (isNewTopology)
    {
        m_numTopologies += 1;

        if (m_numTopologies > m_numTopologyNodes.size())
        {
            m_numTopologyNodes.resize(int(m_numTopologies * 1.5), 0);
            m_numTopologyFaces.resize(int(m_numTopologies * 1.5), 0);
            m_topologyXi.resize(int(m_numTopologies * 1.5), std::vector<double>(maximumNumberOfConnectedNodes, 0));
            m_topologyEta.resize(int(m_numTopologies * 1.5), std::vector<double>(maximumNumberOfConnectedNodes, 0));

            m_topologySharedFaces.resize(int(m_numTopologies * 1.5), std::vector<int>(maximumNumberOfEdgesPerNode, -1));
            m_topologyConnectedNodes.resize(int(m_numTopologies * 1.5), std::vector<std::size_t>(maximumNumberOfConnectedNodes, -1));
            m_topologyFaceNodeMapping.resize(int(m_numTopologies * 1.5), std::vector<std::vector<std::size_t>>(maximumNumberOfConnectedNodes, std::vector<std::size_t>(maximumNumberOfConnectedNodes, -1)));
        }

        int topologyIndex = m_numTopologies - 1;
        m_numTopologyNodes[topologyIndex] = numConnectedNodes;
        m_topologyConnectedNodes[topologyIndex] = connectedNodes;
        m_numTopologyFaces[topologyIndex] = numSharedFaces;
        m_topologySharedFaces[topologyIndex] = sharedFaces;
        m_topologyXi[topologyIndex] = xi;
        m_topologyEta[topologyIndex] = eta;
        m_topologyFaceNodeMapping[topologyIndex] = faceNodeMapping;
        m_nodeTopologyMapping[currentNode] = topologyIndex;
    }

    return true;
}

bool GridGeom::Orthogonalization::GetOrthogonality(const Mesh& mesh, double* orthogonality)
{
    for(int e=0; e < mesh.m_edges.size() ; e++)
    {
        orthogonality[e] = doubleMissingValue;
        int firstVertex = mesh.m_edges[e].first;
        int secondVertex = mesh.m_edges[e].second;

        if (firstVertex!=0 && secondVertex !=0)
        {
            if (e < mesh.m_edgesNumFaces.size() && mesh.m_edgesNumFaces[e]==2 )
            {
                orthogonality[e] = normalizedInnerProductTwoSegments(mesh.m_nodes[firstVertex], mesh.m_nodes[secondVertex],
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
    for (int e = 0; e < mesh.m_edges.size(); e++)
    {
        smoothness[e] = doubleMissingValue;
        int firstVertex = mesh.m_edges[e].first;
        int secondVertex = mesh.m_edges[e].second;

        if (firstVertex != 0 && secondVertex != 0)
        {
            if (e < mesh.m_edgesNumFaces.size() && mesh.m_edgesNumFaces[e] == 2)
            {
                int leftFace = mesh.m_edgesFaces[e][0];
                int rightFace = mesh.m_edgesFaces[e][1];
                double leftFaceArea = mesh.m_faceArea[leftFace];
                double rightFaceArea = mesh.m_faceArea[rightFace];

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