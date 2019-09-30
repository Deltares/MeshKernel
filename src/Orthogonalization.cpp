#pragma once

#define _USE_MATH_DEFINES
#include <vector>
#include <algorithm>
#include <numeric>
#include "Operations.cpp"
#include "Mesh.hpp"
#include "Mesh.cpp"

using namespace GridGeom;

template<typename Mesh>
class Orthogonalization
{
public:

    typedef typename Mesh::Point Point;
    typedef Operations<Point> Operations;

    enum class NodeTypes
    {
        internalNode,
        onRing,
        cornerNode,
        other
    };

    std::vector<std::vector<std::vector<double>>> m_Az;
    std::vector<std::vector<std::vector<double>>> m_Gxi;
    std::vector<std::vector<std::vector<double>>> m_Geta;
    std::vector<std::vector<double>> m_Divxi;
    std::vector<std::vector<double>> m_Diveta;
    std::vector<std::vector<double>> m_Jxi;
    std::vector<std::vector<double>> m_Jeta;
    std::vector<std::vector<double>> m_ww2;

    int m_numTopologies = 0;
    std::vector<int> m_nodeTopologyMapping;
    std::vector<int> m_numTopologyNodes;
    std::vector<int> m_numTopologyFaces;
    std::vector<std::vector<double>> m_topologyXi;
    std::vector<std::vector<double>> m_topologyEta;
    std::vector<std::vector<int>> m_topologySharedFaces;
    std::vector<std::vector<std::vector<size_t>>> m_topologyFaceNodeMapping;
    std::vector < std::vector<size_t>>  m_topologyConnectedNodes;

    std::vector<double> m_aspectRatios;
    std::vector<std::vector<double>> m_ww2Global;
    std::vector<size_t> m_numConnectedNodes;             // nmk2, determined from local node administration
    std::vector<std::vector<size_t>> m_connectedNodes;   // kk2, determined from local node administration
    std::vector<int> m_localCoordinates;                 // iloc

    // run-time options
    bool m_keepCircumcentersAndMassCenters = false;
    double m_atpf = 0.975;                               // Factor(0. <= ATPF <= 1.) between grid smoothing and grid ortho resp.
    double m_atpf_boundary = 1.0;                        // minimum ATPF on the boundary
    double m_smoothorarea = 1.0;                         // Factor between smoother(1.0) and area - homogenizer(0.0)

    bool initialize(const Mesh& mesh)
    {
        m_maxNumNeighbours = *(std::max_element(mesh.m_nodesNumEdges.begin(), mesh.m_nodesNumEdges.end()));
        m_maxNumNeighbours += 1;
        m_nodesNodes.resize(mesh.m_nodes.size(), std::vector<int>(m_maxNumNeighbours, intMissingValue));
        m_weights.resize(mesh.m_nodes.size(), std::vector<double>(m_maxNumNeighbours, 0.0));
        m_rightHandSide.resize(mesh.m_nodes.size(), std::vector<double>(2, 0.0));
        m_aspectRatios.resize(mesh.m_edges.size(), 0.0);

        //for each node, determine the neighbouring nodes
        for (size_t n = 0; n < mesh.m_nodes.size(); n++)
        {
            for (size_t nn = 0; nn < mesh.m_nodesNumEdges[n]; nn++)
            {
                Edge edge = mesh.m_edges[mesh.m_nodesEdges[n][nn]];
                size_t neighbour = edge.first == n ? edge.second : edge.first;
                m_nodesNodes[n][nn] = neighbour;
            }
        }

        // classify the nodes
        classifyNodes(mesh);

        // computes the number of nodes for each face
        computeFacesNumEdges(mesh);

        return true;
    } 

    bool iterate(Mesh& mesh)
    {
        std::vector< std::vector<double>> ww2x;
        std::vector< std::vector<double>> ww2y;
        if (m_smoothorarea != 1)
        {
            ww2x.resize(mesh.m_nodes.size(), std::vector<double>(m_maximumNumConnectedNodes, 0.0));
            ww2y.resize(mesh.m_nodes.size(), std::vector<double>(m_maximumNumConnectedNodes, 0.0));
            // TODO: calculate volume weights
        }

        double mumax = (1.0 - m_smoothorarea) * 0.5;
        double mu = std::min(1e-2, mumax);
        std::vector<double> rightHandSide(2);
        std::vector<double> increments(2);
        std::vector<Point> orthogonalCoordinates(mesh.m_nodes.size());

        // in this case the nearest point is the point itself
        std::vector<int> nearestPoints(mesh.m_nodes.size());
        std::iota(nearestPoints.begin(), nearestPoints.end(), 0);

        // back-up original nodes, for projection on original mesh boundary
        std::vector<Point> originalNodes{mesh.m_nodes};

        for (size_t outerIter = 0; outerIter < orthogonalizationOuterIterations; outerIter++)
        {
            //compute aspect ratios
            aspectRatio(mesh);

            //compute weights orthogonalizer
            computeWeightsOrthogonalizer(mesh);

            //compute weights operators for smoother
            computeSmootherOperators(mesh);

            //compute weights smoother
            computeWeightsSmoother(mesh);

            for (size_t boundaryIter = 0; boundaryIter < orthogonalizationBoundaryIterations; boundaryIter++)
            {
                for (size_t innerIter = 0; innerIter < orthogonalizationInnerIterations; innerIter++)
                {
                    for (int n = 0; n < mesh.m_nodes.size(); n++)
                    {
                        if ((m_nodesTypes[n] != 1 && m_nodesTypes[n] != 2) || mesh.m_nodesNumEdges[n] < 2)
                        {
                            continue;
                        }
                        if (m_keepCircumcentersAndMassCenters != false && (mesh.m_nodesNumEdges[n] != 3 || mesh.m_nodesNumEdges[n] != 1))
                        {
                            continue;
                        }

                        double atpfLoc = m_nodesTypes[n] == 2 ? std::max(m_atpf_boundary, m_atpf) : m_atpf;
                        double atpf1Loc = 1.0 - atpfLoc;

                        double mumat = mu;
                        if ((!ww2x.empty() && ww2x[n][0] != 0.0) && (!ww2y.empty() && ww2y[n][0] != 0.0))
                        {
                            mumat = mu * m_ww2Global[n][0] / std::max(ww2x[n][0], ww2y[n][0]);
                        }

                        int maxnn = std::max(mesh.m_nodesNumEdges[n] + 1, m_numConnectedNodes[n]);
                        double dx0 = 0.0;
                        double dy0 = 0.0;
                        increments[0] = 0.0;
                        increments[1] = 0.0;
                        for (int nn = 1; nn < maxnn; nn++)
                        {
                            double wwx = 0.0;
                            double wwy = 0.0;

                            // Smoother
                            if (atpf1Loc > 0.0 && m_nodesTypes[n] == 1 && !ww2x.empty() && !ww2y.empty())
                            {
                                wwx = atpf1Loc * (mumat * ww2x[n][nn] + m_ww2Global[n][nn]);
                                wwy = atpf1Loc * (mumat * ww2y[n][nn] + m_ww2Global[n][nn]);
                            }
                            else if (atpf1Loc > 0.0 && m_nodesTypes[n] == 1 && ww2x.empty() && ww2y.empty())
                            {
                                wwx = atpf1Loc * m_ww2Global[n][nn];
                                wwy = atpf1Loc * m_ww2Global[n][nn];
                            }

                            // Orthogonalizer
                            int k1;
                            if (nn < mesh.m_nodesNumEdges[n] + 1)
                            {
                                wwx += atpfLoc * m_weights[n][nn - 1];
                                wwy += atpfLoc * m_weights[n][nn - 1];
                                k1 = m_nodesNodes[n][nn - 1];
                            }
                            else
                            {
                                k1 = m_connectedNodes[n][nn];
                            }

                            Operations::orthogonalizationComputeDeltas(k1, n, wwx, wwy, mesh.m_nodes, dx0, dy0, increments);
                        }



                        // combine rhs
                        rightHandSide[0] = atpfLoc * m_rightHandSide[n][0];
                        rightHandSide[1] = atpfLoc * m_rightHandSide[n][1];

                        if (std::abs(increments[0]) > 1e-8 && std::abs(increments[1]) > 1e-8)
                        {
                            dx0 = (dx0 + rightHandSide[0]) / increments[0];
                            dy0 = (dy0 + rightHandSide[1]) / increments[1];
                        }
                        Operations::orthogonalizationComputeCoordinates(dx0, dy0, mesh.m_nodes[n], orthogonalCoordinates[n]);
                    } // n, iteration over nodes

                    // update mesh node coordinates
                    mesh.m_nodes = orthogonalCoordinates;

                    // project on the original net boundary
                    projectOnBoundary(mesh, nearestPoints, originalNodes);

                } // inner iter, inner iteration

            } // boundary iter

            //update mu
            mu = std::min(2.0 * mu, mumax);

            //compute new faces circumcenters
            if (m_keepCircumcentersAndMassCenters != true)
            {
                mesh.faceCircumcenters(1.0);
                mesh.facesAreasAndMassCenters();
            }
        }// outer iter

        return true;
    }

private:

    static constexpr int m_topologyInitialSize = 10;
    static constexpr double m_thetaTolerance = 1e-4;

    std::vector<std::vector<double>>  m_weights;
    std::vector<std::vector<double>>  m_rightHandSide;
    std::vector<int> m_nodesTypes;                             //types of nodes,  1=internal, 2=on ring, 3=corner point, 0/-1=other (e.g. 1d)
    std::vector<int> m_faceNumNodes;                           //number of face nodes

    int m_maximumNumConnectedNodes = 0;
    int m_maximumNumSharedFaces = 0;
    size_t m_maxNumNeighbours;
    std::vector< std::vector<int>> m_nodesNodes;            //node neighbours 

    //local caches (avoid re-allocation)
    std::vector<int> m_boundaryEdges;
    std::vector<double> m_leftXFaceCenter;
    std::vector<double> m_leftYFaceCenter;
    std::vector<double> m_rightXFaceCenter;
    std::vector<double> m_rightYFaceCenter;
    std::vector<double> m_xis;
    std::vector<double> m_etas;

    //orthonet_project_on_boundary: project boundary-nodes back to the boundary of an original net
    bool projectOnBoundary(Mesh& mesh, std::vector<int>& nearestPoints, const std::vector<Point>& originalNodes)
    {
        Point firstPoint;
        Point secondPoint;
        Point thirdPoint;
        Point normalSecondPoint;
        Point normalThirdPoint;
        int leftNode;
        int rightNode;

        for (size_t n = 0; n < mesh.m_nodes.size(); n++)
        {
            int nearestPointIndex = nearestPoints[n];
            if (m_nodesTypes[n] == 2 && mesh.m_nodesNumEdges[n] > 0 && mesh.m_nodesNumEdges[nearestPointIndex] > 0)
            {
                firstPoint = mesh.m_nodes[n];
                int numEdges = mesh.m_nodesNumEdges[nearestPointIndex];
                int numNodes = 0;
                for (size_t nn = 0; nn < numEdges; nn++)
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
                            secondPoint = originalNodes[leftNode];
                        }
                        else if (numNodes == 2)
                        {
                            rightNode = m_nodesNodes[n][nn];
                            if (rightNode == intMissingValue)
                            {
                                return false;
                            }
                            thirdPoint = originalNodes[rightNode];
                        }
                    }
                }

                //Project the moved boundary point back onto the closest ORIGINAL edge(netlink) (either between 0 and 2 or 0 and 3)
                double rl2;
                double dis2 = Operations::distanceFromLine(firstPoint, originalNodes[nearestPointIndex], secondPoint, normalSecondPoint, rl2);

                double rl3;
                double dis3 = Operations::distanceFromLine(firstPoint, originalNodes[nearestPointIndex], thirdPoint, normalThirdPoint, rl3);

                if (dis2 < dis3) 
                {
                    mesh.m_nodes[n] = normalSecondPoint;
                    if (rl2 > 0.5 && m_nodesTypes[n] != 3) 
                    {
                        nearestPoints[n] = leftNode;
                    }
                }
                else 
                {
                    mesh.m_nodes[n] = normalThirdPoint;
                    if (rl3 > 0.5 && m_nodesTypes[n] != 3)
                    {
                        nearestPoints[n] = rightNode;
                    }
                }
            }
        }
        return true;
    }


    // orthonet_compweights_smooth
    // inverse - mapping elliptic smoother
    // computes weight ww in :
    // sum_kk ww(kk, k0)* x1(kk2(kk, k0)) = 0
    // sum_kk ww(kk, k0) * y1(kk2(kk, k0)) = 0
    bool computeWeightsSmoother(const Mesh& mesh)
    {
        std::vector<std::vector<double>> J(mesh.m_nodes.size(), std::vector<double>(4, 0)); //Jacobian
        std::vector<std::vector<double>> Ginv(mesh.m_nodes.size(), std::vector<double>(4, 0)); //mesh monitor matrices

        for (size_t n = 0; n < mesh.m_nodes.size(); n++)
        {
            if (m_nodesTypes[n] != 1 && m_nodesTypes[n] != 2 && m_nodesTypes[n] != 4) continue;
            int currentTopology = m_nodeTopologyMapping[n];
            Operations::orthogonalizationComputeJacobian(n, m_Jxi[currentTopology], m_Jeta[currentTopology], m_topologyConnectedNodes[currentTopology], m_numTopologyNodes[currentTopology], mesh.m_nodes, J[n]);
        }

        // TODO: Account for samples: call orthonet_comp_Ginv(u, ops, J, Ginv)
        for (size_t n = 0; n < mesh.m_nodes.size(); n++)
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
        for (size_t n = 0; n < mesh.m_nodes.size(); n++)
        {

            if (mesh.m_nodesNumEdges[n] < 2) continue;

            // Internal nodes and boundary nodes
            if (m_nodesTypes[n] == 1 || m_nodesTypes[n] == 2)
            {
                int currentTopology = m_nodeTopologyMapping[n];

                Operations::orthogonalizationComputeJacobian(n, m_Jxi[currentTopology], m_Jeta[currentTopology], m_topologyConnectedNodes[currentTopology], m_numTopologyNodes[currentTopology], mesh.m_nodes, J[n]);

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
    
    bool computeSmootherOperators(const Mesh& mesh)
    {
        //allocate small administration arrays only once
        std::vector<int> sharedFaces(maximumNumberOfEdgesPerNode, -1); //icell
        std::vector<size_t> connectedNodes(maximumNumberOfConnectedNodes, 0); //adm%kk2
        std::vector<std::vector<size_t>> faceNodeMapping(maximumNumberOfConnectedNodes, std::vector<size_t>(maximumNumberOfNodesPerFace, 0));//kkc
        std::vector<double> xi(maximumNumberOfConnectedNodes, 0.0);
        std::vector<double> eta(maximumNumberOfConnectedNodes, 0.0);

        m_numConnectedNodes.resize(mesh.m_nodes.size(), 0.0);
        m_connectedNodes.resize(mesh.m_nodes.size(), std::vector<size_t>(maximumNumberOfConnectedNodes));
        initializeTopologies(mesh);
        for (size_t n = 0; n < mesh.m_nodes.size(); n++)
        {
            int numSharedFaces = 0;
            int numConnectedNodes = 0;
            orthogonalizationAdministration(mesh, n, sharedFaces, numSharedFaces, connectedNodes, numConnectedNodes, faceNodeMapping); 
            computeXiEta(mesh, n, sharedFaces, numSharedFaces, connectedNodes, numConnectedNodes, faceNodeMapping, xi, eta);
            saveTopology(n, sharedFaces, numSharedFaces, connectedNodes, numConnectedNodes, faceNodeMapping, xi, eta);
            m_maximumNumConnectedNodes = std::max(m_maximumNumConnectedNodes, numConnectedNodes);
            m_maximumNumSharedFaces = std::max(m_maximumNumSharedFaces, numSharedFaces);
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
        for (size_t n = 0; n < mesh.m_nodes.size(); n++)
        {
            int currentTopology = m_nodeTopologyMapping[n];

            if (isNewTopology[currentTopology])
            {
                isNewTopology[currentTopology] = false;
                // Compute node operators
                allocateNodeOperators(currentTopology);
                computeOperatorsNode(mesh, n,
                    m_numTopologyNodes[currentTopology], m_topologyConnectedNodes[currentTopology],
                    m_numTopologyFaces[currentTopology], m_topologySharedFaces[currentTopology],
                    m_topologyXi[currentTopology], m_topologyEta[currentTopology],
                    m_topologyFaceNodeMapping[currentTopology]);
            }
        }

        //have side effects only for spherical coordinates
        Operations::orthogonalizationComputeLocalCoordinates(mesh.m_nodesNumEdges, m_numConnectedNodes, m_localCoordinates);

        return true;
    }

    //orthonet_comp_operators
    //compute coefficient matrix G of gradient at link
    //(d Phi / d xi)_l = sum_{ k = 1 } ^ nmk2 Gxi_k, l  Phi_k
    //(d Phi / d eta)_l = sum_{ k = 1 } ^ nmk2 Geta_k, l Phi_k
    //compute coefficientmatrix Div of gradient in node
    //d Phi / d xi = sum_{ l = 1 } ^ nmk Divxi_l Phi_l
    //d Phi / d eta = sum_{ l = 1 } ^ nmk Diveta_l Phi_l
    //compute coefficientmatrix Az of cell - center in cell
    //Phi_c = sum_{ l - 1 } ^ nmk Az_l Phi_l
    //Gxi, Geta, Divxi, Divetaand Az are stored in(type tops) op
    bool computeOperatorsNode(const Mesh& mesh, const int currentNode, 
        const size_t& numConnectedNodes, const std::vector<size_t>& connectedNodes,
        const size_t& numSharedFaces, const std::vector<int>& sharedFaces,
        const std::vector<double>& xi, const std::vector<double>& eta,
        const std::vector<std::vector<size_t>>& faceNodeMapping
    )
    {
        // the current topology index
        int topologyIndex = m_nodeTopologyMapping[currentNode];

        // TODO: AVOID DEEP COPIES
        auto Az = m_Az[topologyIndex];
        auto Gxi = m_Gxi[topologyIndex];
        auto Geta = m_Geta[topologyIndex];
        auto Divxi = m_Divxi[topologyIndex];
        auto Diveta = m_Diveta[topologyIndex];
        auto Jxi = m_Jxi[topologyIndex];
        auto Jeta = m_Jeta[topologyIndex];
        auto ww2 = m_ww2[topologyIndex];

        for (int f = 0; f < numSharedFaces; f++)
        {
            if (sharedFaces[f] < 0 || m_nodesTypes[currentNode] == 3) continue;

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
                int nodeIndex = findIndex(mesh.m_facesNodes[sharedFaces[f]], size_t(currentNode));

                int nodeLeft = nodeIndex - 1; if (nodeLeft < 0)nodeLeft += numFaceNodes;
                int nodeRight = nodeIndex + 1; if (nodeRight >= numFaceNodes) nodeRight -= numFaceNodes;
                double alpha = 1.0 / (1.0 - cDPhi * cDPhi + 1e-8);
                double alphaLeft = 0.5 * (1.0 - edgeLeftSquaredDistance / edgeRightSquaredDistance * cDPhi) * alpha;
                double alphaRight = 0.5 * (1.0 - edgeRightSquaredDistance / edgeLeftSquaredDistance * cDPhi) * alpha;

                Az[f][faceNodeMapping[f][nodeIndex]] = 1.0 - (alphaLeft + alphaRight);
                Az[f][faceNodeMapping[f][nodeLeft]] = alphaLeft;
                Az[f][faceNodeMapping[f][nodeRight]] = alphaRight;
            }
            else
            {
                for (int i = 0; i < faceNodeMapping[f].size(); i++)
                {
                    Az[f][faceNodeMapping[f][i]] = 1.0 / double(numFaceNodes);
                }
            }
        }

        // initialize caches
        std::fill(m_boundaryEdges.begin(), m_boundaryEdges.end(), -1);
        std::fill(m_leftXFaceCenter.begin(), m_leftXFaceCenter.end(), 0.0);
        std::fill(m_leftYFaceCenter.begin(), m_leftYFaceCenter.end(), 0.0);
        std::fill(m_rightXFaceCenter.begin(), m_rightXFaceCenter.end(), 0.0);
        std::fill(m_rightYFaceCenter.begin(), m_rightYFaceCenter.end(),0.0);
        std::fill(m_xis.begin(), m_xis.end(), 0.0);
        std::fill(m_etas.begin(), m_etas.end(), 0.0);

        int faceRightIndex = 0;
        int faceLeftIndex = 0;
        double xiBoundary = 0.0;
        double etaBoundary = 0.0;

        for (int f = 0; f < numSharedFaces; f++)
        {
            size_t edgeIndex = mesh.m_nodesEdges[currentNode][f];
            int otherNode = mesh.m_edges[edgeIndex].first + mesh.m_edges[edgeIndex].second - currentNode;
            int leftFace = mesh.m_edgesFaces[edgeIndex][0];
            faceLeftIndex = findIndex(sharedFaces, leftFace);

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
                    leftXi += xi[i] * Az[faceLeftIndex][i];
                    leftEta += eta[i] * Az[faceLeftIndex][i];
                    m_leftXFaceCenter[f] += mesh.m_nodes[connectedNodes[i]].x * Az[faceLeftIndex][i];
                    m_leftYFaceCenter[f] += mesh.m_nodes[connectedNodes[i]].y * Az[faceLeftIndex][i];
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
                    leftXi += xi[i] * Az[faceLeftIndex][i];
                    leftEta += eta[i] * Az[faceLeftIndex][i];
                    rightXi += xi[i] * Az[faceRightIndex][i];
                    rightEta += eta[i] * Az[faceRightIndex][i];

                    m_leftXFaceCenter[f] += mesh.m_nodes[connectedNodes[i]].x * Az[faceLeftIndex][i];
                    m_leftYFaceCenter[f] += mesh.m_nodes[connectedNodes[i]].y * Az[faceLeftIndex][i];
                    m_leftYFaceCenter[f] += mesh.m_nodes[connectedNodes[i]].x * Az[faceRightIndex][i];
                    m_rightYFaceCenter[f] += mesh.m_nodes[connectedNodes[i]].y * Az[faceRightIndex][i];
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
                Gxi[f][i] = facxiL * Az[faceLeftIndex][i];
                Geta[f][i] = facetaL * Az[faceLeftIndex][i];
                if (mesh.m_edgesNumFaces[edgeIndex] == 2)
                {
                    Gxi[f][i] = Gxi[f][i] + facxiR * Az[faceRightIndex][i];
                    Geta[f][i] = Geta[f][i] + facetaR * Az[faceRightIndex][i];
                }
            }


            Gxi[f][node1] = Gxi[f][node1] + facxi1;
            Geta[f][node1] = Geta[f][node1] + faceta1;

            Gxi[f][node0] = Gxi[f][node0] + facxi0;
            Geta[f][node0] = Geta[f][node0] + faceta0;

            //fill the node - based gradient matrix
            Divxi[f] = -eetaLR * leftRightSwap;
            Diveta[f] = exiLR * leftRightSwap;

            // boundary link
            if (mesh.m_edgesNumFaces[edgeIndex] == 1)
            {
                Divxi[f] = 0.5 * Divxi[f] + etaBoundary * leftRightSwap;
                Diveta[f] = 0.5 * Diveta[f] - xiBoundary * leftRightSwap;
            }


        }

        double volxi = 0.0;
        for (int i = 0; i < mesh.m_nodesNumEdges[currentNode]; i++)
        {
            volxi += 0.5 * (Divxi[i] * m_xis[i] + Diveta[i] * m_etas[i]);
        }
        if (volxi == 0.0)volxi = 1.0;

        for (int i = 0; i < mesh.m_nodesNumEdges[currentNode]; i++)
        {
            Divxi[i] = Divxi[i] / volxi;
            Diveta[i] = Diveta[i] / volxi;
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
                    Jxi[i] += Divxi[f] * 0.5 * (Az[f][i] + Az[rightNode][i]);
                    Jeta[i] += Diveta[f] * 0.5 * (Az[f][i] + Az[rightNode][i]);
                }
            }
            else
            {
                Jxi[0] = Jxi[0] + Divxi[f] * 0.5;
                Jxi[f + 1] = Jxi[f + 1] + Divxi[f] * 0.5;
                Jeta[0] = Jeta[0] + Diveta[f] * 0.5;
                Jeta[f + 1] = Jeta[f + 1] + Diveta[f] * 0.5;
            }
        }

        //compute the weights in the Laplacian smoother
        std::fill(ww2.begin(), ww2.end(), 0.0);
        for (int n = 0; n < mesh.m_nodesNumEdges[currentNode]; n++)
        {
            for (int i = 0; i < numConnectedNodes; i++)
            {
                ww2[i] += Divxi[n] * Gxi[n][i] + Diveta[n] * Geta[n][i];
            }
        }

        m_Az[topologyIndex] = Az;
        m_Gxi[topologyIndex] = Gxi;
        m_Geta[topologyIndex] = Geta;
        m_Divxi[topologyIndex] = Divxi;
        m_Diveta[topologyIndex] = Diveta;
        m_Jxi[topologyIndex] = Jxi;
        m_Jeta[topologyIndex] = Jeta;
        m_ww2[topologyIndex] = ww2;

        return true;
    }

    //orthonet_assign_xieta
    // assign xiand eta to all nodes in the stencil
    bool computeXiEta(const Mesh& mesh,
        const int currentNode,
        const std::vector<int>& sharedFaces,
        const int& numSharedFaces,
        const std::vector<size_t>& connectedNodes,
        const size_t& numConnectedNodes,
        const std::vector<std::vector<size_t>>& faceNodeMapping,
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

        if (currentNode == 7)
        {
            std::cout << "Debug" << std::endl;
        }


        //loop over the connected edges
        for (int f = 0; f < numSharedFaces; f++)
        {
            size_t edgeIndex = mesh.m_nodesEdges[currentNode][f];
            size_t nextNode = connectedNodes[f + 1]; // the first entry is always the stencil node 
            int faceLeft = mesh.m_edgesFaces[edgeIndex][0];
            int faceRigth = faceLeft;

            if (mesh.m_edgesNumFaces[edgeIndex] == 2)
                faceRigth = mesh.m_edgesFaces[edgeIndex][1];

            //check if it is a rectangular node (not currentNode itself) 
            bool isSquare = true;
            for (int e = 0; e < mesh.m_nodesNumEdges[nextNode]; e++)
            {
                size_t edge = mesh.m_nodesEdges[nextNode][e];
                for (int ff = 0; ff < mesh.m_edgesNumFaces[edge]; ff++)
                {
                    size_t face = mesh.m_edgesFaces[edge][ff];
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
                if (m_nodesTypes[nextNode] == 1 || m_nodesTypes[nextNode] == 4)
                {
                    // Inner node
                    numNonStencilQuad = mesh.m_nodesNumEdges[nextNode] - 2;
                    thetaSquare[f + 1] = (2.0 - double(numNonStencilQuad) * 0.5) * M_PI;
                }
                if (m_nodesTypes[nextNode] == 2)
                {
                    // boundary node
                    numNonStencilQuad = mesh.m_nodesNumEdges[nextNode] - 1 - mesh.m_edgesNumFaces[edgeIndex];
                    thetaSquare[f + 1] = (1.0 - double(numNonStencilQuad) * 0.5) * M_PI;
                }
                if (m_nodesTypes[nextNode] == 3)
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
        if (m_nodesTypes[currentNode] == 2) factor = 0.5;
        if (m_nodesTypes[currentNode] == 3) factor = 0.25;
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
            //TODO: error
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
                if (m_nodesTypes[currentNode] == 2)
                {
                    dPhi = M_PI;
                }
                else if (m_nodesTypes[currentNode] == 3)
                {
                    dPhi = 1.5 * M_PI;
                }
                else
                {
                    //TODO: error
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
            int nodeIndex = findIndex(mesh.m_facesNodes[sharedFaces[f]], size_t(currentNode));

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

    bool computeFacesNumEdges(const Mesh& mesh)
    {
        // Cache the result for later calls.
        m_faceNumNodes.resize(mesh.m_numFaces, false);
        for (int f = 0; f < mesh.m_numFaces; f++)
        {
            m_faceNumNodes[f] = mesh.m_facesNodes[f].size();
        }
        return true;
    }

    // computes the shared faces and the connected nodes of a stencil node and the faceNodeMapping in the connectedNodes array for each shared face.
    bool orthogonalizationAdministration(const Mesh& mesh, const int currentNode, std::vector<int>& sharedFaces, int& numSharedFaces, std::vector<size_t>& connectedNodes, int& numConnectedNodes, std::vector<std::vector<size_t>>& faceNodeMapping)
    {
        for (auto& f : sharedFaces)  f = -1;

        if (mesh.m_nodesNumEdges[currentNode] < 2) return true;

        // 1. For the currentNode, find the shared faces
        int newFaceIndex = -999;
        numSharedFaces = 0;
        for (int e = 0; e < mesh.m_nodesNumEdges[currentNode]; e++)
        {
            size_t firstEdge = mesh.m_nodesEdges[currentNode][e];

            int secondEdgeIndex = e + 1;

            if (secondEdgeIndex >= mesh.m_nodesNumEdges[currentNode])
                secondEdgeIndex = 0;

            size_t secondEdge = mesh.m_nodesEdges[currentNode][secondEdgeIndex];

            if (mesh.m_edgesNumFaces[firstEdge] < 1 || mesh.m_edgesNumFaces[secondEdge] < 1) continue;

            // find the face that the first and the second edge share
            int firstFaceIndex = std::max(std::min(mesh.m_edgesNumFaces[firstEdge], size_t(2)), size_t(1)) - 1;
            int secondFaceIndex = std::max(std::min(mesh.m_edgesNumFaces[secondEdge], size_t(2)), size_t(1)) - 1;

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
            if (mesh.m_nodesNumEdges[currentNode] == 2 && e == 1 && m_nodesTypes[currentNode] == 3)
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
            size_t edgeIndex = mesh.m_nodesEdges[currentNode][e];
            size_t node = mesh.m_edges[edgeIndex].first + mesh.m_edges[edgeIndex].second - currentNode;
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
            faceNodeMapping.resize(numSharedFaces, std::vector<size_t>(maximumNumberOfNodesPerFace, 0));
        }


        for (int f = 0; f < numSharedFaces; f++)
        {
            int faceIndex = sharedFaces[f];
            if (faceIndex < 0) continue;

            // find the stencil node position  in the current face
            size_t faceNodeIndex = 0;
            size_t numFaceNodes = m_faceNumNodes[faceIndex];
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

    double optimalEdgeAngle(const int numFaceNodes, const double theta1 = -1, const double theta2 = -1, bool isBoundaryEdge = false)
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

    // compute link - based aspect ratios
    // orthonet_compute_aspect
    bool aspectRatio(const Mesh& mesh)
    {
        std::vector<std::vector<double>> averageEdgesLength(mesh.m_edges.size(), std::vector<double>(2, doubleMissingValue));
        std::vector<double> averageFlowEdgesLength(mesh.m_edges.size(), doubleMissingValue);
        std::vector<bool> curvilinearGridIndicator(mesh.m_nodes.size(), true);
        std::vector<double> edgesLength(mesh.m_edges.size(), 0.0);

        for (size_t e = 0; e < mesh.m_edges.size(); e++)
        {
            size_t first = mesh.m_edges[e].first;
            size_t second = mesh.m_edges[e].second;

            if (first == second) continue;
            double edgeLength = Operations::distance(mesh.m_nodes[first], mesh.m_nodes[second]);
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
                double dinry = Operations::innerProductTwoSegments(mesh.m_nodes[first], mesh.m_nodes[second], mesh.m_nodes[first], leftCenter);
                dinry = dinry / std::max(edgeLength * edgeLength, minimumEdgeLength);

                double x0_bc = (1.0 - dinry) * mesh.m_nodes[first].x + dinry * mesh.m_nodes[second].x;
                double y0_bc = (1.0 - dinry) * mesh.m_nodes[first].y + dinry * mesh.m_nodes[second].y;
                rightCenter.x = 2.0 * x0_bc - leftCenter.x;
                rightCenter.y = 2.0 * y0_bc - leftCenter.y;
            }

            averageFlowEdgesLength[e] = Operations::distance(leftCenter, rightCenter);
        }

        // Compute normal length
        for (int f = 0; f < mesh.m_facesNodes.size(); f++)
        {
            size_t numberOfFaceNodes = mesh.m_facesNodes[f].size();
            if (numberOfFaceNodes < 3) continue;

            for (int n = 0; n < numberOfFaceNodes; n++)
            {
                if (numberOfFaceNodes != 4) curvilinearGridIndicator[mesh.m_facesNodes[f][n]] = false;
                size_t edgeIndex = mesh.m_facesEdges[f][n];

                if (mesh.m_edgesNumFaces[edgeIndex] < 1) continue;

                //get the other links in the right numbering
                //TODO: ask why only 3 are requested, why an average lenght stored in averageEdgesLength is needed?
                //int kkm1 = n - 1; if (kkm1 < 0) kkm1 = kkm1 + numberOfFaceNodes;
                //int kkp1 = n + 1; if (kkp1 >= numberOfFaceNodes) kkp1 = kkp1 - numberOfFaceNodes;
                //
                //size_t klinkm1 = mesh.m_facesEdges[f][kkm1];
                //size_t klinkp1 = mesh.m_facesEdges[f][kkp1];
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
                    size_t klinkp2 = mesh.m_facesEdges[f][kkp2];
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

        for (size_t e = 0; e < mesh.m_edges.size(); e++)
        {
            size_t first = mesh.m_edges[e].first;
            size_t second = mesh.m_edges[e].second;

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

    //MAKENETNODESCODING
    bool classifyNodes(const Mesh& mesh)
    {
        m_nodesTypes.resize(mesh.m_nodes.size(), 0);

        for (int e = 0; e < mesh.m_edges.size(); e++)
        {
            size_t first = mesh.m_edges[e].first;
            size_t second = mesh.m_edges[e].second;

            if (mesh.m_edgesNumFaces[e] == 0)
            {
                m_nodesTypes[first] = -1;
                m_nodesTypes[second] = -1;
            }
            else if (mesh.m_edgesNumFaces[e] == 1)
            {
                m_nodesTypes[first] += 1;
                m_nodesTypes[second] += 1;
            }
        }

        for (int n = 0; n < mesh.m_nodes.size(); n++)
        {
            if (m_nodesTypes[n] == 1 || m_nodesTypes[n] == 2)
            {
                if (mesh.m_nodesNumEdges[n] == 2)
                {
                    //corner point
                    m_nodesTypes[n] = 3;
                }
                else {}
            }
            else if (m_nodesTypes[n] > 2)
            {
                // corner point
                m_nodesTypes[n] = 3;
            }
            else if (m_nodesTypes[n] != -1)
            {
                //internal node
                m_nodesTypes[n] = 1;
            }

            if (mesh.m_nodesNumEdges[n] < 2)
            {
                //hanging node
                m_nodesTypes[n] = -1;
            }
        }
        return true;
    }

    template<typename T>
    inline T findIndex(const std::vector<T>& vec, const T& el)
    {
        T index = 0;
        for (int n = 0; n < vec.size(); n++)
        {
            if (vec[n] == el)
            {
                index = n;
                break;
            }
        }
        return index;
    }

    //orthonet_compweights
    // compute weights wwand right - hand side rhs in orthogonizer :
    // sum_kk ww(kk, k0)* (x1(kk1(kk, k0)) - x1(k0)) = rhs(1, k0)
    // sum_kk ww(kk, k0) * (y1(kk1(kk, k0)) - y1(k0)) = rhs(2, k0)
    bool computeWeightsOrthogonalizer(const Mesh& mesh)
    {
        double localOrthogonalizationToSmoothingFactor = 1.0;
        double localOrthogonalizationToSmoothingFactorSymmetric = 1.0 - localOrthogonalizationToSmoothingFactor;
        constexpr double mu = 1.0;
        Point normal;

        std::fill(m_rightHandSide.begin(), m_rightHandSide.end(), std::vector<double>(2, 0.0));
        for (size_t n = 0; n < mesh.m_nodes.size(); n++)
        {
            if (m_nodesTypes[n] != 1 && m_nodesTypes[n] != 2)
            {
                continue; 
            }

            for (size_t nn = 0; nn < mesh.m_nodesNumEdges[n]; nn++)
            {
                size_t edgeIndex = mesh.m_nodesEdges[n][nn];
                double aspectRatio = m_aspectRatios[edgeIndex];

                if (aspectRatio != doubleMissingValue)
                {
                    m_weights[n][nn] = localOrthogonalizationToSmoothingFactor * aspectRatio + localOrthogonalizationToSmoothingFactorSymmetric * mu;

                    if (mesh.m_edgesNumFaces[edgeIndex] == 1)
                    {
                        //boundary nodes
                        Point neighbouringNode = mesh.m_nodes[m_nodesNodes[n][nn]];
                        double neighbouringNodeDistance = Operations::distance(neighbouringNode, mesh.m_nodes[n]);
                        double aspectRatioByNodeDistance = aspectRatio * neighbouringNodeDistance;

                        size_t leftFace = mesh.m_edgesFaces[edgeIndex][0];
                        bool flippedNormal;
                        Operations::normalVectorInside(mesh.m_nodes[n], neighbouringNode, mesh.m_facesMassCenters[leftFace], normal, flippedNormal);

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

    double matrixNorm(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& matCoefficents)
    {
        double norm = (matCoefficents[0] * x[0] + matCoefficents[1] * x[1]) * y[0] + (matCoefficents[2] * x[0] + matCoefficents[3] * x[1]) * y[1];
        return norm;
    }

    bool initializeTopologies(const Mesh& mesh)
    {
        // topology 
        m_numTopologies = 0;
        m_nodeTopologyMapping.resize(mesh.m_nodes.size(), -1);
        m_numTopologyNodes.resize(m_topologyInitialSize, -1);
        m_numTopologyFaces.resize(m_topologyInitialSize, -1);
        m_topologyXi.resize(m_topologyInitialSize, std::vector<double>(maximumNumberOfConnectedNodes, 0));
        m_topologyEta.resize(m_topologyInitialSize, std::vector<double>(maximumNumberOfConnectedNodes, 0));
        m_topologySharedFaces.resize(m_topologyInitialSize, std::vector<int>(maximumNumberOfConnectedNodes, -1));
        m_topologyConnectedNodes.resize(m_topologyInitialSize, std::vector<size_t>(maximumNumberOfConnectedNodes, -1));
        m_topologyFaceNodeMapping.resize(m_topologyInitialSize, std::vector<std::vector<size_t>>(maximumNumberOfConnectedNodes, std::vector<size_t>(maximumNumberOfConnectedNodes, -1)));

        return true;
    }

    bool allocateNodeOperators(const int topologyIndex)
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

    // save only the unique topologies
    bool saveTopology(const int currentNode,
        const std::vector<int>& sharedFaces,
        const int numSharedFaces,
        const std::vector<size_t>& connectedNodes,
        const int numConnectedNodes,
        const std::vector<std::vector<size_t>>& faceNodeMapping,
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
                m_topologyConnectedNodes.resize(int(m_numTopologies * 1.5), std::vector<size_t>(maximumNumberOfConnectedNodes, -1));
                m_topologyFaceNodeMapping.resize(int(m_numTopologies * 1.5), std::vector<std::vector<size_t>>(maximumNumberOfConnectedNodes, std::vector<size_t>(maximumNumberOfConnectedNodes, -1)));
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
};