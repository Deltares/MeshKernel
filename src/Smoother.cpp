//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2020.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation version 3.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// contact: delft3d.support@deltares.nl
// Stichting Deltares
// P.O. Box 177
// 2600 MH Delft, The Netherlands
//
// All indications and logos of, and references to, "Delft3D" and "Deltares"
// are registered trademarks of Stichting Deltares, and remain the property of
// Stichting Deltares. All rights reserved.
//
//------------------------------------------------------------------------------

#pragma once

#include <vector>
#include <algorithm>
#include <stdexcept>

#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/Operations.cpp>
#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Mesh.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Smoother.hpp>

meshkernel::Smoother::Smoother(std::shared_ptr<Mesh> mesh) : m_mesh(mesh)
{
}

void meshkernel::Smoother::Compute()
{
    // compute smoother topologies
    ComputeTopologies();

    // compute smoother operators
    ComputeOperators();

    // compute weights
    ComputeWeights();
}

void meshkernel::Smoother::ComputeTopologies()
{
    Initialize();

    for (auto n = 0; n < m_mesh->GetNumNodes(); n++)
    {
        int numSharedFaces = 0;
        int numConnectedNodes = 0;

        std::fill(m_sharedFacesCache.begin(), m_sharedFacesCache.end(), -1);
        std::fill(m_connectedNodesCache.begin(), m_connectedNodesCache.end(), 0);
        NodeAdministration(n, numSharedFaces, numConnectedNodes);

        std::fill(m_xiCache.begin(), m_xiCache.end(), 0.0);
        std::fill(m_etaCache.begin(), m_etaCache.end(), 0.0);
        ComputeNodeXiEta(n, numSharedFaces, numConnectedNodes);

        SaveNodeTopologyIfNeeded(n, numSharedFaces, numConnectedNodes);

        m_maximumNumConnectedNodes = std::max(m_maximumNumConnectedNodes, numConnectedNodes);
        m_maximumNumSharedFaces = std::max(m_maximumNumSharedFaces, numSharedFaces);
    }
}

void meshkernel::Smoother::ComputeOperators()
{
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
    m_boundaryEdgesCache.resize(2, -1);
    m_leftXFaceCenterCache.resize(maximumNumberOfEdgesPerNode, 0.0);
    m_leftYFaceCenterCache.resize(maximumNumberOfEdgesPerNode, 0.0);
    m_rightXFaceCenterCache.resize(maximumNumberOfEdgesPerNode, 0.0);
    m_rightYFaceCenterCache.resize(maximumNumberOfEdgesPerNode, 0.0);
    m_xisCache.resize(maximumNumberOfEdgesPerNode, 0.0);
    m_etasCache.resize(maximumNumberOfEdgesPerNode, 0.0);

    std::vector<bool> isNewTopology(m_numTopologies, true);

    for (auto n = 0; n < m_mesh->GetNumNodes(); n++)
    {
        if (m_mesh->m_nodesTypes[n] != 1 && m_mesh->m_nodesTypes[n] != 2 && m_mesh->m_nodesTypes[n] != 3 && m_mesh->m_nodesTypes[n] != 4)
        {
            continue;
        }

        // for each node, the associated topology
        int currentTopology = m_nodeTopologyMapping[n];

        if (isNewTopology[currentTopology])
        {
            isNewTopology[currentTopology] = false;

            // Compute node operators

            AllocateNodeOperators(currentTopology);
            ComputeOperatorsNode(n);
        }
    }
}

void meshkernel::Smoother::ComputeWeights()
{
    std::vector<std::vector<double>> J(m_mesh->GetNumNodes(), std::vector<double>(4, 0));    // Jacobian
    std::vector<std::vector<double>> Ginv(m_mesh->GetNumNodes(), std::vector<double>(4, 0)); // Mesh monitor matrices

    for (auto n = 0; n < m_mesh->GetNumNodes(); n++)
    {
        if (m_mesh->m_nodesTypes[n] != 1 && m_mesh->m_nodesTypes[n] != 2 && m_mesh->m_nodesTypes[n] != 4)
        {
            continue;
        }
        ComputeJacobian(n, J[n]);
    }

    // TODO: Account for samples: call orthonet_comp_Ginv(u, ops, J, Ginv)
    for (auto n = 0; n < m_mesh->GetNumNodes(); n++)
    {
        Ginv[n][0] = 1.0;
        Ginv[n][1] = 0.0;
        Ginv[n][2] = 0.0;
        Ginv[n][3] = 1.0;
    }

    m_weights.resize(m_mesh->GetNumNodes(), std::vector<double>(m_maximumNumConnectedNodes, 0.0));
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
    for (auto n = 0; n < m_mesh->GetNumNodes(); n++)
    {

        if (m_mesh->m_nodesNumEdges[n] < 2)
            continue;

        // Internal nodes and boundary nodes
        if (m_mesh->m_nodesTypes[n] == 1 || m_mesh->m_nodesTypes[n] == 2)
        {
            const auto currentTopology = m_nodeTopologyMapping[n];

            ComputeJacobian(n, J[n]);

            //compute the contravariant base vectors
            const double determinant = J[n][0] * J[n][3] - J[n][3] * J[n][1];
            if (determinant == 0.0)
            {
                continue;
            }

            a1[0] = J[n][3] / determinant;
            a1[1] = -J[n][2] / determinant;
            a2[0] = -J[n][1] / determinant;
            a2[1] = J[n][0] / determinant;

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
                m_weights[n][i] -= MatrixNorm(a1, a1, DGinvDxi) * m_Jxi[currentTopology][i] +
                                   MatrixNorm(a1, a2, DGinvDeta) * m_Jxi[currentTopology][i] +
                                   MatrixNorm(a2, a1, DGinvDxi) * m_Jeta[currentTopology][i] +
                                   MatrixNorm(a2, a2, DGinvDeta) * m_Jeta[currentTopology][i];
                m_weights[n][i] += MatrixNorm(a1, a1, currentGinv) * GxiByDivxi[i] +
                                   MatrixNorm(a1, a2, currentGinv) * GxiByDiveta[i] +
                                   MatrixNorm(a2, a1, currentGinv) * GetaByDivxi[i] +
                                   MatrixNorm(a2, a2, currentGinv) * GetaByDiveta[i];
            }

            double alpha = 0.0;
            for (int i = 1; i < m_numTopologyNodes[currentTopology]; i++)
            {
                alpha = std::max(alpha, -m_weights[n][i]) / std::max(1.0, m_ww2[currentTopology][i]);
            }

            double sumValues = 0.0;
            for (int i = 1; i < m_numTopologyNodes[currentTopology]; i++)
            {
                m_weights[n][i] = m_weights[n][i] + alpha * std::max(1.0, m_ww2[currentTopology][i]);
                sumValues += m_weights[n][i];
            }
            m_weights[n][0] = -sumValues;
            for (int i = 0; i < m_numTopologyNodes[currentTopology]; i++)
            {
                m_weights[n][i] = -m_weights[n][i] / (-sumValues + 1e-8);
            }
        }
    }
}

void meshkernel::Smoother::ComputeOperatorsNode(int currentNode)
{
    // the current topology index
    const int currentTopology = m_nodeTopologyMapping[currentNode];

    for (int f = 0; f < m_numTopologyFaces[currentTopology]; f++)
    {
        if (m_topologySharedFaces[currentTopology][f] < 0 || m_mesh->m_nodesTypes[currentNode] == 3)
        {
            continue;
        }

        int edgeLeft = f + 1;
        int edgeRight = edgeLeft + 1;
        if (edgeRight > m_numTopologyFaces[currentTopology])
        {
            edgeRight -= m_numTopologyFaces[currentTopology];
        }

        const auto xiLeft = m_topologyXi[currentTopology][edgeLeft];
        const auto xiRight = m_topologyXi[currentTopology][edgeRight];
        const auto etaLeft = m_topologyEta[currentTopology][edgeLeft];
        const auto etaRight = m_topologyEta[currentTopology][edgeRight];

        const double edgeLeftSquaredDistance = std::sqrt(xiLeft * xiLeft + etaLeft * etaLeft + 1e-16);
        const double edgeRightSquaredDistance = std::sqrt(xiRight * xiRight + etaRight * etaRight + 1e-16);
        const double cPhi = (xiLeft * xiRight + etaLeft * etaRight) / (edgeLeftSquaredDistance * edgeRightSquaredDistance);
        const auto numFaceNodes = m_mesh->GetNumFaceEdges(m_topologySharedFaces[currentTopology][f]);

        // the value of xi and eta needs to be estimated at the circumcenters, calculated the contributions of each node
        if (numFaceNodes == 3)
        {
            // for triangular faces
            int nodeIndex = FindIndex(m_mesh->m_facesNodes[m_topologySharedFaces[currentTopology][f]], currentNode);
            const auto nodeLeft = NextCircularBackwardIndex(nodeIndex, numFaceNodes);
            const auto nodeRight = NextCircularForwardIndex(nodeIndex, numFaceNodes);

            double alpha = 1.0 / (1.0 - cPhi * cPhi + 1e-8);
            double alphaLeft = 0.5 * (1.0 - edgeLeftSquaredDistance / edgeRightSquaredDistance * cPhi) * alpha;
            double alphaRight = 0.5 * (1.0 - edgeRightSquaredDistance / edgeLeftSquaredDistance * cPhi) * alpha;

            m_Az[currentTopology][f][m_topologyFaceNodeMapping[currentTopology][f][nodeIndex]] = 1.0 - (alphaLeft + alphaRight);
            m_Az[currentTopology][f][m_topologyFaceNodeMapping[currentTopology][f][nodeLeft]] = alphaLeft;
            m_Az[currentTopology][f][m_topologyFaceNodeMapping[currentTopology][f][nodeRight]] = alphaRight;
        }
        else
        {
            // for non-triangular faces
            for (const auto& element : m_topologyFaceNodeMapping[currentTopology][f])
            {
                m_Az[currentTopology][f][element] = 1.0 / static_cast<double>(numFaceNodes);
            }
        }
    }

    // Initialize caches
    std::fill(m_boundaryEdgesCache.begin(), m_boundaryEdgesCache.end(), -1);
    std::fill(m_leftXFaceCenterCache.begin(), m_leftXFaceCenterCache.end(), 0.0);
    std::fill(m_leftYFaceCenterCache.begin(), m_leftYFaceCenterCache.end(), 0.0);
    std::fill(m_rightXFaceCenterCache.begin(), m_rightXFaceCenterCache.end(), 0.0);
    std::fill(m_rightYFaceCenterCache.begin(), m_rightYFaceCenterCache.end(), 0.0);
    std::fill(m_xisCache.begin(), m_xisCache.end(), 0.0);
    std::fill(m_etasCache.begin(), m_etasCache.end(), 0.0);

    int faceRightIndex = 0;
    int faceLeftIndex = 0;
    double xiBoundary = 0.0;
    double etaBoundary = 0.0;

    for (int f = 0; f < m_numTopologyFaces[currentTopology]; f++)
    {
        auto edgeIndex = m_mesh->m_nodesEdges[currentNode][f];
        int otherNode = m_mesh->m_edges[edgeIndex].first + m_mesh->m_edges[edgeIndex].second - currentNode;
        int leftFace = m_mesh->m_edgesFaces[edgeIndex][0];
        faceLeftIndex = FindIndex(m_topologySharedFaces[currentTopology], leftFace);

        if (m_topologySharedFaces[currentTopology][faceLeftIndex] != leftFace)
        {
            throw std::invalid_argument("Smoother::ComputeOperatorsNode: Face could not be found, this happens when the cell is outside of the polygon.");
        }

        //by construction
        double xiOne = m_topologyXi[currentTopology][f + 1];
        double etaOne = m_topologyEta[currentTopology][f + 1];

        double leftRightSwap = 1.0;
        double leftXi = 0.0;
        double leftEta = 0.0;
        double rightXi = 0.0;
        double rightEta = 0.0;
        double alpha_x = 0.0;
        if (m_mesh->IsEdgeOnBoundary(edgeIndex))
        {
            // Boundary face
            if (m_boundaryEdgesCache[0] < 0)
            {
                m_boundaryEdgesCache[0] = f;
            }
            else
            {
                m_boundaryEdgesCache[1] = f;
            }

            // Swap left and right if the boundary is at the left
            if (f != faceLeftIndex)
            {
                leftRightSwap = -1.0;
            }

            // Compute the face circumcenter
            for (int i = 0; i < m_numTopologyNodes[currentTopology]; i++)
            {
                leftXi += m_topologyXi[currentTopology][i] * m_Az[currentTopology][faceLeftIndex][i];
                leftEta += m_topologyEta[currentTopology][i] * m_Az[currentTopology][faceLeftIndex][i];
                m_leftXFaceCenterCache[f] += m_mesh->m_nodes[m_topologyConnectedNodes[currentTopology][i]].x * m_Az[currentTopology][faceLeftIndex][i];
                m_leftYFaceCenterCache[f] += m_mesh->m_nodes[m_topologyConnectedNodes[currentTopology][i]].y * m_Az[currentTopology][faceLeftIndex][i];
            }

            double alpha = leftXi * xiOne + leftEta * etaOne;
            alpha = alpha / (xiOne * xiOne + etaOne * etaOne);

            alpha_x = alpha;
            xiBoundary = alpha * xiOne;
            etaBoundary = alpha * etaOne;

            rightXi = 2.0 * xiBoundary - leftXi;
            rightEta = 2.0 * etaBoundary - leftEta;

            const double xBc = (1.0 - alpha) * m_mesh->m_nodes[currentNode].x + alpha * m_mesh->m_nodes[otherNode].x;
            const double yBc = (1.0 - alpha) * m_mesh->m_nodes[currentNode].y + alpha * m_mesh->m_nodes[otherNode].y;
            m_leftYFaceCenterCache[f] = 2.0 * xBc - m_leftXFaceCenterCache[f];
            m_rightYFaceCenterCache[f] = 2.0 * yBc - m_leftYFaceCenterCache[f];
        }
        else
        {
            faceLeftIndex = f;
            faceRightIndex = NextCircularBackwardIndex(faceLeftIndex, m_numTopologyFaces[currentTopology]);

            if (faceRightIndex < 0)
                continue;

            auto faceLeft = m_topologySharedFaces[currentTopology][faceLeftIndex];
            auto faceRight = m_topologySharedFaces[currentTopology][faceRightIndex];

            if ((faceLeft != m_mesh->m_edgesFaces[edgeIndex][0] && faceLeft != m_mesh->m_edgesFaces[edgeIndex][1]) ||
                (faceRight != m_mesh->m_edgesFaces[edgeIndex][0] && faceRight != m_mesh->m_edgesFaces[edgeIndex][1]))
            {
                throw std::invalid_argument("Smoother::ComputeOperatorsNode: Invalid argument.");
            }

            for (int i = 0; i < m_numTopologyNodes[currentTopology]; i++)
            {
                leftXi += m_topologyXi[currentTopology][i] * m_Az[currentTopology][faceLeftIndex][i];
                leftEta += m_topologyEta[currentTopology][i] * m_Az[currentTopology][faceLeftIndex][i];
                rightXi += m_topologyXi[currentTopology][i] * m_Az[currentTopology][faceRightIndex][i];
                rightEta += m_topologyEta[currentTopology][i] * m_Az[currentTopology][faceRightIndex][i];

                m_leftXFaceCenterCache[f] += m_mesh->m_nodes[m_topologyConnectedNodes[currentTopology][i]].x * m_Az[currentTopology][faceLeftIndex][i];
                m_leftYFaceCenterCache[f] += m_mesh->m_nodes[m_topologyConnectedNodes[currentTopology][i]].y * m_Az[currentTopology][faceLeftIndex][i];
                m_rightXFaceCenterCache[f] += m_mesh->m_nodes[m_topologyConnectedNodes[currentTopology][i]].x * m_Az[currentTopology][faceRightIndex][i];
                m_rightYFaceCenterCache[f] += m_mesh->m_nodes[m_topologyConnectedNodes[currentTopology][i]].y * m_Az[currentTopology][faceRightIndex][i];
            }
        }

        m_xisCache[f] = 0.5 * (leftXi + rightXi);
        m_etasCache[f] = 0.5 * (leftEta + rightEta);

        const double exiLR = rightXi - leftXi;
        const double eetaLR = rightEta - leftEta;
        const double exi01 = xiOne;
        const double eeta01 = etaOne;

        const double fac = 1.0 / std::abs(exi01 * eetaLR - eeta01 * exiLR + 1e-16);
        double facxi1 = -eetaLR * fac * leftRightSwap;

        double facxi0 = -facxi1;
        double faceta1 = exiLR * fac * leftRightSwap;
        double faceta0 = -faceta1;
        const double facxiR = eeta01 * fac * leftRightSwap;
        double facxiL = -facxiR;
        const double facetaR = -exi01 * fac * leftRightSwap;
        double facetaL = -facetaR;

        //boundary link
        if (m_mesh->IsEdgeOnBoundary(edgeIndex))
        {
            facxi1 += -facxiL * 2.0 * alpha_x;
            facxi0 += -facxiL * 2.0 * (1.0 - alpha_x);
            facxiL += facxiL;
            //note that facxiR does not exist
            faceta1 += -facetaL * 2.0 * alpha_x;
            faceta0 += -facetaL * 2.0 * (1.0 - alpha_x);
            facetaL = 2.0 * facetaL;
        }

        int node1 = f + 1;
        int node0 = 0;
        for (int i = 0; i < m_numTopologyNodes[currentTopology]; i++)
        {
            m_Gxi[currentTopology][f][i] = facxiL * m_Az[currentTopology][faceLeftIndex][i];
            m_Geta[currentTopology][f][i] = facetaL * m_Az[currentTopology][faceLeftIndex][i];
            if (!m_mesh->IsEdgeOnBoundary(edgeIndex))
            {
                m_Gxi[currentTopology][f][i] += facxiR * m_Az[currentTopology][faceRightIndex][i];
                m_Geta[currentTopology][f][i] += facetaR * m_Az[currentTopology][faceRightIndex][i];
            }
        }

        m_Gxi[currentTopology][f][node1] += facxi1;
        m_Geta[currentTopology][f][node1] += faceta1;

        m_Gxi[currentTopology][f][node0] += facxi0;
        m_Geta[currentTopology][f][node0] += faceta0;

        //fill the node-based gradient matrix
        m_Divxi[currentTopology][f] = -eetaLR * leftRightSwap;
        m_Diveta[currentTopology][f] = exiLR * leftRightSwap;

        // boundary link
        if (m_mesh->IsEdgeOnBoundary(edgeIndex))
        {
            m_Divxi[currentTopology][f] = 0.5 * m_Divxi[currentTopology][f] + etaBoundary * leftRightSwap;
            m_Diveta[currentTopology][f] = 0.5 * m_Diveta[currentTopology][f] - xiBoundary * leftRightSwap;
        }
    }

    double volxi = 0.0;
    for (int i = 0; i < m_mesh->m_nodesNumEdges[currentNode]; i++)
    {
        volxi += 0.5 * (m_Divxi[currentTopology][i] * m_xisCache[i] + m_Diveta[currentTopology][i] * m_etasCache[i]);
    }
    if (volxi == 0.0)
    {
        volxi = 1.0;
    }

    for (int i = 0; i < m_mesh->m_nodesNumEdges[currentNode]; i++)
    {
        m_Divxi[currentTopology][i] = m_Divxi[currentTopology][i] / volxi;
        m_Diveta[currentTopology][i] = m_Diveta[currentTopology][i] / volxi;
    }

    //compute the node-to-node gradients
    for (int f = 0; f < m_numTopologyFaces[currentTopology]; f++)
    {
        // internal edge
        if (!m_mesh->IsEdgeOnBoundary(m_mesh->m_nodesEdges[currentNode][f]))
        {
            int rightNode = f - 1;
            if (rightNode < 0)
            {
                rightNode += m_mesh->m_nodesNumEdges[currentNode];
            }
            for (int i = 0; i < m_numTopologyNodes[currentTopology]; i++)
            {
                m_Jxi[currentTopology][i] += m_Divxi[currentTopology][f] * 0.5 * (m_Az[currentTopology][f][i] + m_Az[currentTopology][rightNode][i]);
                m_Jeta[currentTopology][i] += m_Diveta[currentTopology][f] * 0.5 * (m_Az[currentTopology][f][i] + m_Az[currentTopology][rightNode][i]);
            }
        }
        else
        {
            m_Jxi[currentTopology][0] += m_Divxi[currentTopology][f] * 0.5;
            m_Jxi[currentTopology][f + 1] += m_Divxi[currentTopology][f] * 0.5;
            m_Jeta[currentTopology][0] += m_Diveta[currentTopology][f] * 0.5;
            m_Jeta[currentTopology][f + 1] += m_Diveta[currentTopology][f] * 0.5;
        }
    }

    //compute the weights in the Laplacian smoother
    std::fill(m_ww2[currentTopology].begin(), m_ww2[currentTopology].end(), 0.0);
    for (int n = 0; n < m_mesh->m_nodesNumEdges[currentNode]; n++)
    {
        for (int i = 0; i < m_numTopologyNodes[currentTopology]; i++)
        {
            m_ww2[currentTopology][i] += m_Divxi[currentTopology][n] * m_Gxi[currentTopology][n][i] + m_Diveta[currentTopology][n] * m_Geta[currentTopology][n][i];
        }
    }
}

void meshkernel::Smoother::ComputeNodeXiEta(int currentNode,
                                            int numSharedFaces,
                                            int numConnectedNodes)
{
    // the angles for the squared nodes connected to the stencil nodes, first the ones directly connected, then the others
    std::vector<double> thetaSquare(numConnectedNodes, doubleMissingValue);
    // for each shared face, a boolean indicating if it is squared or not
    std::vector<bool> isSquareFace(numSharedFaces, false);

    int numNonStencilQuad = 0;

    //loop over the connected edges
    for (int f = 0; f < numSharedFaces; f++)
    {
        auto edgeIndex = m_mesh->m_nodesEdges[currentNode][f];
        auto nextNode = m_connectedNodesCache[f + 1]; // the first entry is always the stencil node
        int faceLeft = m_mesh->m_edgesFaces[edgeIndex][0];
        int faceRigth = faceLeft;

        if (!m_mesh->IsEdgeOnBoundary(edgeIndex))
        {
            faceRigth = m_mesh->m_edgesFaces[edgeIndex][1];
        }

        //check if it is a rectangular node (not currentNode itself)
        bool isSquare = true;
        for (int e = 0; e < m_mesh->m_nodesNumEdges[nextNode]; e++)
        {
            auto edge = m_mesh->m_nodesEdges[nextNode][e];
            for (int ff = 0; ff < m_mesh->m_edgesNumFaces[edge]; ff++)
            {
                auto face = m_mesh->m_edgesFaces[edge][ff];
                if (face != faceLeft && face != faceRigth)
                {
                    isSquare = isSquare && m_mesh->GetNumFaceEdges(face) == 4;
                }
            }
            if (!isSquare)
            {
                break;
            }
        }

        //Compute optimal angle thetaSquare based on the node type
        int leftFaceIndex = f - 1;
        if (leftFaceIndex < 0)
        {
            leftFaceIndex = leftFaceIndex + numSharedFaces;
        }

        if (isSquare)
        {
            if (m_mesh->m_nodesTypes[nextNode] == 1 || m_mesh->m_nodesTypes[nextNode] == 4)
            {
                // Inner node
                numNonStencilQuad = m_mesh->m_nodesNumEdges[nextNode] - 2;
                thetaSquare[f + 1] = (2.0 - double(numNonStencilQuad) * 0.5) * M_PI;
            }
            if (m_mesh->m_nodesTypes[nextNode] == 2)
            {
                // boundary node
                numNonStencilQuad = m_mesh->m_nodesNumEdges[nextNode] - 1 - m_mesh->m_edgesNumFaces[edgeIndex];
                thetaSquare[f + 1] = (1.0 - double(numNonStencilQuad) * 0.5) * M_PI;
            }
            if (m_mesh->m_nodesTypes[nextNode] == 3)
            {
                //corner node
                thetaSquare[f + 1] = 0.5 * M_PI;
            }

            if (m_sharedFacesCache[f] > 1 && m_mesh->GetNumFaceEdges(m_sharedFacesCache[f]) == 4)
            {
                numNonStencilQuad += 1;
            }
            if (m_sharedFacesCache[leftFaceIndex] > 1 && m_mesh->GetNumFaceEdges(m_sharedFacesCache[leftFaceIndex]) == 4)
            {
                numNonStencilQuad += 1;
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
        if (m_sharedFacesCache[f] < 0)
            continue;

        // non boundary face
        if (m_mesh->GetNumFaceEdges(m_sharedFacesCache[f]) == 4)
        {
            for (int n = 0; n < m_mesh->GetNumFaceEdges(m_sharedFacesCache[f]); n++)
            {
                if (m_faceNodeMappingCache[f][n] <= numSharedFaces)
                {
                    continue;
                }
                thetaSquare[m_faceNodeMappingCache[f][n]] = 0.5 * M_PI;
            }
        }
    }

    // Compute internal angle
    int numSquaredTriangles = 0;
    int numTriangles = 0;
    double phiSquaredTriangles = 0.0;
    double phiQuads = 0.0;
    double phiTriangles = 0.0;
    double phiTot = 0.0;
    numNonStencilQuad = 0;
    for (int f = 0; f < numSharedFaces; f++)
    {
        // boundary face
        if (m_sharedFacesCache[f] < 0)
        {
            continue;
        }

        auto numFaceNodes = m_mesh->GetNumFaceEdges(m_sharedFacesCache[f]);
        double phi = OptimalEdgeAngle(numFaceNodes);

        if (isSquareFace[f] || numFaceNodes == 4)
        {
            int nextNode = f + 2;
            if (nextNode > numSharedFaces)
            {
                nextNode = nextNode - numSharedFaces;
            }

            phi = OptimalEdgeAngle(numFaceNodes, thetaSquare[f + 1], thetaSquare[nextNode], m_mesh->IsEdgeOnBoundary(m_mesh->m_nodesEdges[currentNode][f]));
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
    if (m_mesh->m_nodesTypes[currentNode] == 2)
    {
        factor = 0.5;
    }

    if (m_mesh->m_nodesTypes[currentNode] == 3)
    {
        factor = 0.25;
    }

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
        //TODO: add cirr(xk(k0), yk(k0), ncolhl)
        m_nodeXErrors.emplace_back(m_mesh->m_nodes[currentNode].x);
        m_nodeXErrors.emplace_back(m_mesh->m_nodes[currentNode].y);
        throw AlgorithmError("Smoother::ComputeNodeXiEta: Fatal error (phiTot=0).");
    }

    double phi0 = 0.0;
    double dPhi0 = 0.0;
    double dPhi = 0.0;
    double dTheta = 0.0;
    for (int f = 0; f < numSharedFaces; f++)
    {
        phi0 = phi0 + 0.5 * dPhi;
        if (m_sharedFacesCache[f] < 0)
        {
            if (m_mesh->m_nodesTypes[currentNode] == 2)
            {
                dPhi = M_PI;
            }
            else if (m_mesh->m_nodesTypes[currentNode] == 3)
            {
                dPhi = 1.5 * M_PI;
            }
            else
            {
                //TODO: add cirr(xk(k0), yk(k0), ncolhl)
                m_nodeXErrors.emplace_back(m_mesh->m_nodes[currentNode].x);
                m_nodeXErrors.emplace_back(m_mesh->m_nodes[currentNode].y);
                throw AlgorithmError("Smoother::ComputeNodeXiEta: Inappropriate fictitious boundary cell.");
            }
            phi0 = phi0 + 0.5 * dPhi;
            continue;
        }

        int numFaceNodes = m_mesh->GetNumFaceEdges(m_sharedFacesCache[f]);
        if (numFaceNodes > maximumNumberOfEdgesPerNode)
        {
            throw AlgorithmError("Smoother::ComputeNodeXiEta: The number of face nodes is greater than the maximum number of edges per node.");
        }

        dPhi0 = OptimalEdgeAngle(numFaceNodes);
        if (isSquareFace[f])
        {
            int nextNode = f + 2;
            if (nextNode > numSharedFaces)
            {
                nextNode = nextNode - numSharedFaces;
            }

            dPhi0 = OptimalEdgeAngle(numFaceNodes, thetaSquare[f + 1], thetaSquare[nextNode], m_mesh->IsEdgeOnBoundary(m_mesh->m_nodesEdges[currentNode][f]));
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
        const auto nodeIndex = FindIndex(m_mesh->m_facesNodes[m_sharedFacesCache[f]], currentNode);

        // optimal angle
        dTheta = 2.0 * M_PI / double(numFaceNodes);

        // orientation of the face (necessary for folded cells)
        int previousNode = NextCircularForwardIndex(nodeIndex, numFaceNodes);
        int nextNode = NextCircularBackwardIndex(nodeIndex, numFaceNodes);

        if ((m_faceNodeMappingCache[f][nextNode] - m_faceNodeMappingCache[f][previousNode]) == -1 ||
            (m_faceNodeMappingCache[f][nextNode] - m_faceNodeMappingCache[f][previousNode]) == m_mesh->m_nodesNumEdges[currentNode])
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

            m_xiCache[m_faceNodeMappingCache[f][n]] = xip * std::cos(phi0) - aspectRatio * ethap * std::sin(phi0);
            m_etaCache[m_faceNodeMappingCache[f][n]] = xip * std::sin(phi0) + aspectRatio * ethap * std::cos(phi0);
        }
    }
}

void meshkernel::Smoother::NodeAdministration(int currentNode,
                                              int& numSharedFaces,
                                              int& numConnectedNodes)
{

    numSharedFaces = 0;
    numConnectedNodes = 0;
    if (m_mesh->m_nodesNumEdges[currentNode] < 2)
    {
        return;
    }

    // For the currentNode, find the shared faces
    int newFaceIndex = intMissingValue;
    for (int e = 0; e < m_mesh->m_nodesNumEdges[currentNode]; e++)
    {
        const auto firstEdge = m_mesh->m_nodesEdges[currentNode][e];

        int secondEdgeIndex = e + 1;
        if (secondEdgeIndex >= m_mesh->m_nodesNumEdges[currentNode])
        {
            secondEdgeIndex = 0;
        }

        const auto secondEdge = m_mesh->m_nodesEdges[currentNode][secondEdgeIndex];
        if (m_mesh->m_edgesNumFaces[firstEdge] < 1 || m_mesh->m_edgesNumFaces[secondEdge] < 1)
        {
            continue;
        }

        // find the face shared by the two edges
        const int firstFace = std::max(std::min(m_mesh->m_edgesNumFaces[firstEdge], int(2)), int(1)) - 1;
        const int secondFace = std::max(std::min(m_mesh->m_edgesNumFaces[secondEdge], int(2)), int(1)) - 1;

        if (m_mesh->m_edgesFaces[firstEdge][0] != newFaceIndex &&
            (m_mesh->m_edgesFaces[firstEdge][0] == m_mesh->m_edgesFaces[secondEdge][0] ||
             m_mesh->m_edgesFaces[firstEdge][0] == m_mesh->m_edgesFaces[secondEdge][secondFace]))
        {
            newFaceIndex = m_mesh->m_edgesFaces[firstEdge][0];
        }
        else if (m_mesh->m_edgesFaces[firstEdge][firstFace] != newFaceIndex &&
                 (m_mesh->m_edgesFaces[firstEdge][firstFace] == m_mesh->m_edgesFaces[secondEdge][0] ||
                  m_mesh->m_edgesFaces[firstEdge][firstFace] == m_mesh->m_edgesFaces[secondEdge][secondFace]))
        {
            newFaceIndex = m_mesh->m_edgesFaces[firstEdge][firstFace];
        }
        else
        {
            newFaceIndex = intMissingValue;
        }

        //corner face (already found in the first iteration)
        if (m_mesh->m_nodesNumEdges[currentNode] == 2 && e == 1 && m_mesh->m_nodesTypes[currentNode] == 3 && m_sharedFacesCache[0] == newFaceIndex)
        {
            newFaceIndex = intMissingValue;
        }
        m_sharedFacesCache[numSharedFaces] = newFaceIndex;
        numSharedFaces += 1;
    }

    // no shared face found
    if (numSharedFaces < 1)
    {
        return;
    }

    int connectedNodesIndex = 0;
    m_connectedNodesCache[connectedNodesIndex] = currentNode;

    // edge connected nodes
    for (int e = 0; e < m_mesh->m_nodesNumEdges[currentNode]; e++)
    {
        const auto edgeIndex = m_mesh->m_nodesEdges[currentNode][e];
        const auto node = m_mesh->m_edges[edgeIndex].first + m_mesh->m_edges[edgeIndex].second - currentNode;
        connectedNodesIndex++;
        m_connectedNodesCache[connectedNodesIndex] = node;
        if (int(m_connectedNodes[currentNode].size()) < connectedNodesIndex + 1)
        {
            m_connectedNodes[currentNode].resize(connectedNodesIndex);
        }
        m_connectedNodes[currentNode][connectedNodesIndex] = node;
    }

    // for each face store the positions of the its nodes in the connectedNodes (compressed array)
    if (m_faceNodeMappingCache.size() < numSharedFaces)
    {
        m_faceNodeMappingCache.resize(numSharedFaces, std::vector<size_t>(maximumNumberOfNodesPerFace, 0));
    }
    for (int f = 0; f < numSharedFaces; f++)
    {
        const auto faceIndex = m_sharedFacesCache[f];
        if (faceIndex < 0)
        {
            continue;
        }

        // find the stencil node position  in the current face
        int faceNodeIndex = 0;
        const auto numFaceNodes = m_mesh->GetNumFaceEdges(faceIndex);
        for (int n = 0; n < numFaceNodes; n++)
        {
            if (m_mesh->m_facesNodes[faceIndex][n] == currentNode)
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

            const auto node = m_mesh->m_facesNodes[faceIndex][faceNodeIndex];

            bool isNewNode = true;
            for (int i = 0; i < connectedNodesIndex + 1; i++)
            {
                if (node == m_connectedNodesCache[i])
                {
                    isNewNode = false;
                    m_faceNodeMappingCache[f][faceNodeIndex] = i;
                    break;
                }
            }

            if (isNewNode)
            {
                connectedNodesIndex++;
                m_connectedNodesCache[connectedNodesIndex] = node;
                m_faceNodeMappingCache[f][faceNodeIndex] = connectedNodesIndex;
                if (int(m_connectedNodes[currentNode].size()) < connectedNodesIndex + 1)
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

    //update connected nodes (kkc)
    m_numConnectedNodes[currentNode] = numConnectedNodes;
}

double meshkernel::Smoother::OptimalEdgeAngle(int numFaceNodes, double theta1, double theta2, bool isBoundaryEdge) const
{
    double angle = M_PI * (1.0 - 2.0 / double(numFaceNodes));

    if (std::abs(theta1 + 1.0) > std::numeric_limits<double>::epsilon() &&
        std::abs(theta2 + 1.0) > std::numeric_limits<double>::epsilon() &&
        numFaceNodes == 3)
    {
        angle = 0.25 * M_PI;
        if (std::abs(theta1 + theta2 - M_PI) <= std::numeric_limits<double>::epsilon() && !isBoundaryEdge)
        {
            angle = 0.5 * M_PI;
        }
    }
    return angle;
}

double meshkernel::Smoother::MatrixNorm(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& matCoefficents) const
{
    const double norm = (matCoefficents[0] * x[0] + matCoefficents[1] * x[1]) * y[0] + (matCoefficents[2] * x[0] + matCoefficents[3] * x[1]) * y[1];
    return norm;
}

void meshkernel::Smoother::Initialize()
{
    // local matrices caches
    m_numConnectedNodes.resize(m_mesh->GetNumNodes());
    std::fill(m_numConnectedNodes.begin(), m_numConnectedNodes.end(), 0.0);

    m_connectedNodes.resize(m_mesh->GetNumNodes());
    std::fill(m_connectedNodes.begin(), m_connectedNodes.end(), std::vector<size_t>(maximumNumberOfConnectedNodes, 0));

    m_sharedFacesCache.resize(maximumNumberOfEdgesPerNode, -1);
    std::fill(m_sharedFacesCache.begin(), m_sharedFacesCache.end(), -1);

    m_connectedNodesCache.resize(maximumNumberOfConnectedNodes, 0);
    std::fill(m_connectedNodesCache.begin(), m_connectedNodesCache.end(), 0);

    m_faceNodeMappingCache.resize(maximumNumberOfConnectedNodes);
    std::fill(m_faceNodeMappingCache.begin(), m_faceNodeMappingCache.end(), std::vector<size_t>(maximumNumberOfNodesPerFace, 0));

    m_xiCache.resize(maximumNumberOfConnectedNodes, 0.0);
    std::fill(m_xiCache.begin(), m_xiCache.end(), 0.0);

    m_etaCache.resize(maximumNumberOfConnectedNodes, 0.0);
    std::fill(m_etaCache.begin(), m_etaCache.end(), 0.0);

    // topology
    m_numTopologies = 0;

    m_nodeTopologyMapping.resize(m_mesh->GetNumNodes());
    std::fill(m_nodeTopologyMapping.begin(), m_nodeTopologyMapping.end(), -1);

    m_numTopologyNodes.resize(m_topologyInitialSize);
    std::fill(m_numTopologyNodes.begin(), m_numTopologyNodes.end(), -1);

    m_numTopologyFaces.resize(m_topologyInitialSize);
    std::fill(m_numTopologyFaces.begin(), m_numTopologyFaces.end(), -1);

    m_topologyXi.resize(m_topologyInitialSize);
    std::fill(m_topologyXi.begin(), m_topologyXi.end(), std::vector<double>(maximumNumberOfConnectedNodes, 0.0));

    m_topologyEta.resize(m_topologyInitialSize);
    std::fill(m_topologyEta.begin(), m_topologyEta.end(), std::vector<double>(maximumNumberOfConnectedNodes, 0.0));

    m_topologySharedFaces.resize(m_topologyInitialSize);
    std::fill(m_topologySharedFaces.begin(), m_topologySharedFaces.end(), std::vector<int>(maximumNumberOfConnectedNodes, -1));

    m_topologyConnectedNodes.resize(m_topologyInitialSize);
    std::fill(m_topologyConnectedNodes.begin(), m_topologyConnectedNodes.end(), std::vector<size_t>(maximumNumberOfConnectedNodes, 0));

    m_topologyFaceNodeMapping.resize(m_topologyInitialSize);
    std::fill(m_topologyFaceNodeMapping.begin(), m_topologyFaceNodeMapping.end(), std::vector<std::vector<size_t>>(maximumNumberOfConnectedNodes, std::vector<size_t>(maximumNumberOfConnectedNodes, 0)));
}

void meshkernel::Smoother::AllocateNodeOperators(int topologyIndex)
{
    const auto numSharedFaces = m_numTopologyFaces[topologyIndex];
    const auto numConnectedNodes = m_numTopologyNodes[topologyIndex];

    //// will reallocate only if necessary
    m_Az[topologyIndex].resize(numSharedFaces);
    std::fill(m_Az[topologyIndex].begin(), m_Az[topologyIndex].end(), std::vector<double>(numConnectedNodes, 0.0));

    m_Gxi[topologyIndex].resize(numSharedFaces);
    std::fill(m_Gxi[topologyIndex].begin(), m_Gxi[topologyIndex].end(), std::vector<double>(numConnectedNodes, 0.0));

    m_Geta[topologyIndex].resize(numSharedFaces);
    std::fill(m_Geta[topologyIndex].begin(), m_Geta[topologyIndex].end(), std::vector<double>(numConnectedNodes, 0.0));

    m_Divxi[topologyIndex].resize(numSharedFaces);
    std::fill(m_Divxi[topologyIndex].begin(), m_Divxi[topologyIndex].end(), 0.0);

    m_Diveta[topologyIndex].resize(numSharedFaces);
    std::fill(m_Diveta[topologyIndex].begin(), m_Diveta[topologyIndex].end(), 0.0);

    m_Jxi[topologyIndex].resize(numConnectedNodes);
    std::fill(m_Jxi[topologyIndex].begin(), m_Jxi[topologyIndex].end(), 0.0);

    m_Jeta[topologyIndex].resize(numConnectedNodes);
    std::fill(m_Jeta[topologyIndex].begin(), m_Jeta[topologyIndex].end(), 0.0);

    m_ww2[topologyIndex].resize(numConnectedNodes);
    std::fill(m_ww2[topologyIndex].begin(), m_ww2[topologyIndex].end(), 0.0);
}

void meshkernel::Smoother::SaveNodeTopologyIfNeeded(int currentNode,
                                                    int numSharedFaces,
                                                    int numConnectedNodes)
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
            double thetaLoc = std::atan2(m_etaCache[n], m_xiCache[n]);
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
            m_topologyConnectedNodes.resize(int(m_numTopologies * 1.5), std::vector<size_t>(maximumNumberOfConnectedNodes, 0));
            m_topologyFaceNodeMapping.resize(int(m_numTopologies * 1.5), std::vector<std::vector<size_t>>(maximumNumberOfConnectedNodes, std::vector<size_t>(maximumNumberOfConnectedNodes, 0)));
        }

        int topologyIndex = m_numTopologies - 1;
        m_numTopologyNodes[topologyIndex] = numConnectedNodes;
        m_topologyConnectedNodes[topologyIndex] = m_connectedNodesCache;
        m_numTopologyFaces[topologyIndex] = numSharedFaces;
        m_topologySharedFaces[topologyIndex] = m_sharedFacesCache;
        m_topologyXi[topologyIndex] = m_xiCache;
        m_topologyEta[topologyIndex] = m_etaCache;
        m_topologyFaceNodeMapping[topologyIndex] = m_faceNodeMappingCache;
        m_nodeTopologyMapping[currentNode] = topologyIndex;
    }
}

void meshkernel::Smoother::ComputeJacobian(int currentNode, std::vector<double>& J) const
{
    const auto currentTopology = m_nodeTopologyMapping[currentNode];
    const auto numNodes = m_numTopologyNodes[currentTopology];
    if (m_mesh->m_projection == Projections::cartesian)
    {
        J[0] = 0.0;
        J[1] = 0.0;
        J[2] = 0.0;
        J[3] = 0.0;
        for (int i = 0; i < numNodes; i++)
        {
            J[0] += m_Jxi[currentTopology][i] * m_mesh->m_nodes[m_topologyConnectedNodes[currentTopology][i]].x;
            J[1] += m_Jxi[currentTopology][i] * m_mesh->m_nodes[m_topologyConnectedNodes[currentTopology][i]].y;
            J[2] += m_Jeta[currentTopology][i] * m_mesh->m_nodes[m_topologyConnectedNodes[currentTopology][i]].x;
            J[3] += m_Jeta[currentTopology][i] * m_mesh->m_nodes[m_topologyConnectedNodes[currentTopology][i]].y;
        }
    }
    if (m_mesh->m_projection == Projections::spherical || m_mesh->m_projection == Projections::sphericalAccurate)
    {
        double cosFactor = std::cos(m_mesh->m_nodes[currentNode].y * degrad_hp);
        J[0] = 0.0;
        J[1] = 0.0;
        J[2] = 0.0;
        J[3] = 0.0;
        for (int i = 0; i < numNodes; i++)
        {
            J[0] += m_Jxi[currentTopology][i] * m_mesh->m_nodes[m_topologyConnectedNodes[currentTopology][i]].x * cosFactor;
            J[1] += m_Jxi[currentTopology][i] * m_mesh->m_nodes[m_topologyConnectedNodes[currentTopology][i]].y;
            J[2] += m_Jeta[currentTopology][i] * m_mesh->m_nodes[m_topologyConnectedNodes[currentTopology][i]].x * cosFactor;
            J[3] += m_Jeta[currentTopology][i] * m_mesh->m_nodes[m_topologyConnectedNodes[currentTopology][i]].y;
        }
    }
}
