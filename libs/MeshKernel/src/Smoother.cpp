//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2021.
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

#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Smoother.hpp>

using meshkernel::Mesh2D;
using meshkernel::Smoother;

Smoother::Smoother(std::shared_ptr<Mesh2D> mesh) : m_mesh(mesh)
{
}

void Smoother::Compute()
{
    // compute smoother topologies
    ComputeTopologies();

    // compute smoother operators
    ComputeOperators();

    // compute weights
    ComputeWeights();
}

void Smoother::ComputeTopologies()
{
    Initialize();

    for (Index n = 0; n < m_mesh->GetNumNodes(); n++)
    {
        NodeAdministration(n);

        ComputeNodeXiEta(n);

        SaveNodeTopologyIfNeeded(n);

        m_maximumNumConnectedNodes = std::max(m_maximumNumConnectedNodes, static_cast<Index>(m_connectedNodesCache.size()));
        m_maximumNumSharedFaces = std::max(m_maximumNumSharedFaces, static_cast<Index>(m_sharedFacesCache.size()));
    }
}

void Smoother::ComputeOperators()
{
    // allocate local operators for unique topologies
    m_Az.resize(m_topologyConnectedNodes.size());
    m_Gxi.resize(m_topologyConnectedNodes.size());
    m_Geta.resize(m_topologyConnectedNodes.size());
    m_Divxi.resize(m_topologyConnectedNodes.size());
    m_Diveta.resize(m_topologyConnectedNodes.size());
    m_Jxi.resize(m_topologyConnectedNodes.size());
    m_Jeta.resize(m_topologyConnectedNodes.size());
    m_ww2.resize(m_topologyConnectedNodes.size());

    // allocate caches
    m_boundaryEdgesCache.resize(2, constants::missing::sizetValue);
    m_leftXFaceCenterCache.resize(Mesh::m_maximumNumberOfEdgesPerNode, 0.0);
    m_leftYFaceCenterCache.resize(Mesh::m_maximumNumberOfEdgesPerNode, 0.0);
    m_rightXFaceCenterCache.resize(Mesh::m_maximumNumberOfEdgesPerNode, 0.0);
    m_rightYFaceCenterCache.resize(Mesh::m_maximumNumberOfEdgesPerNode, 0.0);
    m_xisCache.resize(Mesh::m_maximumNumberOfEdgesPerNode, 0.0);
    m_etasCache.resize(Mesh::m_maximumNumberOfEdgesPerNode, 0.0);

    std::vector<bool> isNewTopology(m_topologyConnectedNodes.size(), true);

    for (Index n = 0; n < m_mesh->GetNumNodes(); n++)
    {
        if (m_mesh->m_nodesTypes[n] != 1 && m_mesh->m_nodesTypes[n] != 2 && m_mesh->m_nodesTypes[n] != 3 && m_mesh->m_nodesTypes[n] != 4)
        {
            continue;
        }

        // for each node, the associated topology
        const auto currentTopology = m_nodeTopologyMapping[n];

        if (isNewTopology[currentTopology])
        {
            isNewTopology[currentTopology] = false;

            // Compute node operators

            AllocateNodeOperators(currentTopology);
            ComputeOperatorsNode(n);
        }
    }
}

void Smoother::ComputeWeights()
{
    std::vector<std::vector<double>> J(m_mesh->GetNumNodes(), std::vector<double>(4, 0));    // Jacobian
    std::vector<std::vector<double>> Ginv(m_mesh->GetNumNodes(), std::vector<double>(4, 0)); // Mesh2D monitor matrices

    for (Index n = 0; n < m_mesh->GetNumNodes(); n++)
    {
        if (m_mesh->m_nodesTypes[n] != 1 && m_mesh->m_nodesTypes[n] != 2 && m_mesh->m_nodesTypes[n] != 4)
        {
            continue;
        }
        ComputeJacobian(n, J[n]);
    }

    // TODO: Account for samples: call orthonet_comp_Ginv(u, ops, J, Ginv)
    for (Index n = 0; n < m_mesh->GetNumNodes(); n++)
    {
        Ginv[n][0] = 1.0;
        Ginv[n][1] = 0.0;
        Ginv[n][2] = 0.0;
        Ginv[n][3] = 1.0;
    }

    ResizeAndFill2DVector(m_weights, m_mesh->GetNumNodes(), m_maximumNumConnectedNodes, true, 0.0);
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
    for (Index n = 0; n < m_mesh->GetNumNodes(); n++)
    {

        if (m_mesh->m_nodesNumEdges[n] < 2)
            continue;

        // Internal nodes and boundary nodes
        if (m_mesh->m_nodesTypes[n] == 1 || m_mesh->m_nodesTypes[n] == 2)
        {
            const auto currentTopology = m_nodeTopologyMapping[n];

            ComputeJacobian(n, J[n]);

            // compute the contravariant base vectors
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
            for (Index i = 0; i < m_topologyConnectedNodes[currentTopology].size(); i++)
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
            for (Index i = 0; i < m_topologyConnectedNodes[currentTopology].size(); i++)
            {
                for (Index j = 0; j < m_Divxi[currentTopology].size(); j++)
                {
                    GxiByDivxi[i] += m_Gxi[currentTopology][j][i] * m_Divxi[currentTopology][j];
                    GxiByDiveta[i] += m_Gxi[currentTopology][j][i] * m_Diveta[currentTopology][j];
                    GetaByDivxi[i] += m_Geta[currentTopology][j][i] * m_Divxi[currentTopology][j];
                    GetaByDiveta[i] += m_Geta[currentTopology][j][i] * m_Diveta[currentTopology][j];
                }
            }

            for (Index i = 0; i < m_topologyConnectedNodes[currentTopology].size(); i++)
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
            for (Index i = 1; i < m_topologyConnectedNodes[currentTopology].size(); i++)
            {
                alpha = std::max(alpha, -m_weights[n][i]) / std::max(1.0, m_ww2[currentTopology][i]);
            }

            double sumValues = 0.0;
            for (Index i = 1; i < m_topologyConnectedNodes[currentTopology].size(); i++)
            {
                m_weights[n][i] = m_weights[n][i] + alpha * std::max(1.0, m_ww2[currentTopology][i]);
                sumValues += m_weights[n][i];
            }
            m_weights[n][0] = -sumValues;
            for (Index i = 0; i < m_topologyConnectedNodes[currentTopology].size(); i++)
            {
                m_weights[n][i] = -m_weights[n][i] / (-sumValues + 1e-8);
            }
        }
    }
}

void Smoother::ComputeOperatorsNode(Index currentNode)
{
    // the current topology index
    const auto currentTopology = m_nodeTopologyMapping[currentNode];

    for (Index f = 0; f < m_topologySharedFaces[currentTopology].size(); f++)
    {
        if (m_topologySharedFaces[currentTopology][f] == constants::missing::sizetValue || m_mesh->m_nodesTypes[currentNode] == 3)
        {
            continue;
        }

        Index edgeLeft = f + 1;
        Index edgeRight = edgeLeft + 1;
        if (edgeRight > m_topologySharedFaces[currentTopology].size())
        {
            edgeRight -= m_topologySharedFaces[currentTopology].size();
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
            const auto nodeIndex = FindIndex(m_mesh->m_facesNodes[m_topologySharedFaces[currentTopology][f]], static_cast<Index>(currentNode));
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
            for (Index i = 0; i < numFaceNodes; ++i)
            {
                const auto element = m_topologyFaceNodeMapping[currentTopology][f][i];
                m_Az[currentTopology][f][element] = 1.0 / static_cast<double>(numFaceNodes);
            }
        }
    }

    // Initialize caches
    std::fill(m_boundaryEdgesCache.begin(), m_boundaryEdgesCache.end(), constants::missing::sizetValue);
    std::fill(m_leftXFaceCenterCache.begin(), m_leftXFaceCenterCache.end(), 0.0);
    std::fill(m_leftYFaceCenterCache.begin(), m_leftYFaceCenterCache.end(), 0.0);
    std::fill(m_rightXFaceCenterCache.begin(), m_rightXFaceCenterCache.end(), 0.0);
    std::fill(m_rightYFaceCenterCache.begin(), m_rightYFaceCenterCache.end(), 0.0);
    std::fill(m_xisCache.begin(), m_xisCache.end(), 0.0);
    std::fill(m_etasCache.begin(), m_etasCache.end(), 0.0);

    Index faceRightIndex = 0;
    Index faceLeftIndex = 0;
    double xiBoundary = 0.0;
    double etaBoundary = 0.0;

    for (Index f = 0; f < m_topologySharedFaces[currentTopology].size(); f++)
    {
        auto edgeIndex = m_mesh->m_nodesEdges[currentNode][f];
        auto otherNode = OtherNodeOfEdge(m_mesh->m_edges[edgeIndex], currentNode);

        auto leftFace = m_mesh->m_edgesFaces[edgeIndex][0];
        faceLeftIndex = FindIndex(m_topologySharedFaces[currentTopology], leftFace);

        if (m_topologySharedFaces[currentTopology][faceLeftIndex] != leftFace)
        {
            throw std::invalid_argument("Smoother::ComputeOperatorsNode: Face could not be found, this happens when the face is outside of the polygon.");
        }

        // by construction
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
            if (m_boundaryEdgesCache[0] == constants::missing::sizetValue)
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
            for (Index i = 0; i < m_topologyConnectedNodes[currentTopology].size(); i++)
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
            faceRightIndex = NextCircularBackwardIndex(faceLeftIndex, m_topologySharedFaces[currentTopology].size());

            if (faceRightIndex == constants::missing::sizetValue)
            {
                continue;
            }

            const auto faceLeft = m_topologySharedFaces[currentTopology][faceLeftIndex];
            const auto faceRight = m_topologySharedFaces[currentTopology][faceRightIndex];

            if ((faceLeft != m_mesh->m_edgesFaces[edgeIndex][0] && faceLeft != m_mesh->m_edgesFaces[edgeIndex][1]) ||
                (faceRight != m_mesh->m_edgesFaces[edgeIndex][0] && faceRight != m_mesh->m_edgesFaces[edgeIndex][1]))
            {
                throw std::invalid_argument("Smoother::ComputeOperatorsNode: Invalid argument.");
            }

            for (Index i = 0; i < m_topologyConnectedNodes[currentTopology].size(); i++)
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

        // boundary link
        if (m_mesh->IsEdgeOnBoundary(edgeIndex))
        {
            facxi1 += -facxiL * 2.0 * alpha_x;
            facxi0 += -facxiL * 2.0 * (1.0 - alpha_x);
            facxiL += facxiL;
            // note that facxiR does not exist
            faceta1 += -facetaL * 2.0 * alpha_x;
            faceta0 += -facetaL * 2.0 * (1.0 - alpha_x);
            facetaL = 2.0 * facetaL;
        }

        Index node1 = f + 1;
        Index node0 = 0;
        for (Index i = 0; i < m_topologyConnectedNodes[currentTopology].size(); i++)
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

        // fill the node-based gradient matrix
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
    for (Index i = 0; i < m_mesh->m_nodesNumEdges[currentNode]; i++)
    {
        volxi += 0.5 * (m_Divxi[currentTopology][i] * m_xisCache[i] + m_Diveta[currentTopology][i] * m_etasCache[i]);
    }
    if (volxi == 0.0)
    {
        volxi = 1.0;
    }

    for (Index i = 0; i < m_mesh->m_nodesNumEdges[currentNode]; i++)
    {
        m_Divxi[currentTopology][i] = m_Divxi[currentTopology][i] / volxi;
        m_Diveta[currentTopology][i] = m_Diveta[currentTopology][i] / volxi;
    }

    // compute the node-to-node gradients
    for (Index f = 0; f < m_topologySharedFaces[currentTopology].size(); f++)
    {
        // internal edge
        if (!m_mesh->IsEdgeOnBoundary(m_mesh->m_nodesEdges[currentNode][f]))
        {
            Index rightNode;
            if (f == 0)
            {
                rightNode = f + m_mesh->m_nodesNumEdges[currentNode] - 1;
            }
            else
            {
                rightNode = f - 1;
            }

            for (Index i = 0; i < m_topologyConnectedNodes[currentTopology].size(); i++)
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

    // compute the weights in the Laplacian smoother
    std::fill(m_ww2[currentTopology].begin(), m_ww2[currentTopology].end(), 0.0);
    for (Index n = 0; n < m_mesh->m_nodesNumEdges[currentNode]; n++)
    {
        for (Index i = 0; i < m_topologyConnectedNodes[currentTopology].size(); i++)
        {
            m_ww2[currentTopology][i] += m_Divxi[currentTopology][n] * m_Gxi[currentTopology][n][i] + m_Diveta[currentTopology][n] * m_Geta[currentTopology][n][i];
        }
    }
}

void Smoother::ComputeNodeXiEta(Index currentNode)
{
    // set caches to 0
    std::fill(m_xiCache.begin(), m_xiCache.end(), 0.0);
    std::fill(m_etaCache.begin(), m_etaCache.end(), 0.0);

    const auto numSharedFaces = m_sharedFacesCache.size();
    // the angles for the squared nodes connected to the stencil nodes, first the ones directly connected, then the others
    std::vector<double> thetaSquare(m_connectedNodesCache.size(), constants::missing::doubleValue);
    // for each shared face, a boolean indicating if it is squared or not
    std::vector<bool> isSquareFace(numSharedFaces, false);

    Index numNonStencilQuad = 0;

    // loop over the connected edges
    for (Index f = 0; f < numSharedFaces; f++)
    {
        auto edgeIndex = m_mesh->m_nodesEdges[currentNode][f];
        auto nextNode = m_connectedNodesCache[f + 1]; // the first entry is always the stencil node
        auto faceLeft = m_mesh->m_edgesFaces[edgeIndex][0];
        auto faceRight = faceLeft;

        if (!m_mesh->IsEdgeOnBoundary(edgeIndex))
        {
            faceRight = m_mesh->m_edgesFaces[edgeIndex][1];
        }

        // check if it is a rectangular node (not currentNode itself)
        bool isSquare = true;
        for (Index e = 0; e < m_mesh->m_nodesNumEdges[nextNode]; e++)
        {
            auto edge = m_mesh->m_nodesEdges[nextNode][e];
            for (Index ff = 0; ff < m_mesh->m_edgesNumFaces[edge]; ff++)
            {
                auto face = m_mesh->m_edgesFaces[edge][ff];
                if (face != faceLeft && face != faceRight)
                {
                    isSquare = isSquare && m_mesh->GetNumFaceEdges(face) == 4;
                }
            }
            if (!isSquare)
            {
                break;
            }
        }

        // Compute optimal angle thetaSquare based on the node type
        Index leftFaceIndex;
        if (f == 0)
        {
            leftFaceIndex = f + numSharedFaces - static_cast<Index>(1);
        }
        else
        {
            leftFaceIndex = f - static_cast<Index>(1);
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
                // corner node
                thetaSquare[f + 1] = 0.5 * M_PI;
            }

            if (m_sharedFacesCache[f] != constants::missing::sizetValue && m_sharedFacesCache[f] > 1 && m_mesh->GetNumFaceEdges(m_sharedFacesCache[f]) == 4)
            {
                numNonStencilQuad += 1;
            }
            if (m_sharedFacesCache[leftFaceIndex] != constants::missing::sizetValue && m_sharedFacesCache[leftFaceIndex] > 1 && m_mesh->GetNumFaceEdges(m_sharedFacesCache[leftFaceIndex]) == 4)
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

    for (Index f = 0; f < numSharedFaces; f++)
    {
        // boundary face
        if (m_sharedFacesCache[f] == constants::missing::sizetValue)
            continue;

        // non boundary face
        if (m_mesh->GetNumFaceEdges(m_sharedFacesCache[f]) == 4)
        {
            for (Index n = 0; n < m_mesh->GetNumFaceEdges(m_sharedFacesCache[f]); n++)
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
    Index numSquaredTriangles = 0;
    Index numTriangles = 0;
    double phiSquaredTriangles = 0.0;
    double phiQuads = 0.0;
    double phiTriangles = 0.0;
    double phiTot = 0.0;
    numNonStencilQuad = 0;
    for (Index f = 0; f < numSharedFaces; f++)
    {
        // boundary face
        if (m_sharedFacesCache[f] == constants::missing::sizetValue)
        {
            continue;
        }

        auto numFaceNodes = m_mesh->GetNumFaceEdges(m_sharedFacesCache[f]);
        double phi = OptimalEdgeAngle(numFaceNodes);

        if (isSquareFace[f] || numFaceNodes == 4)
        {
            Index nextNode = static_cast<Index>(f) + static_cast<Index>(2);
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
        muSquaredTriangles = std::max(factor * 2.0 * M_PI - (phiTot - phiSquaredTriangles), static_cast<double>(numSquaredTriangles) * minPhi) / phiSquaredTriangles;
    }

    if (phiTot > 1e-18)
    {
        mu = factor * 2.0 * M_PI / (phiTot - (1.0 - muTriangles) * phiTriangles - (1.0 - muSquaredTriangles) * phiSquaredTriangles);
    }
    else if (numSharedFaces > 0)
    {
        throw MeshGeometryError("Smoother::ComputeNodeXiEta: Fatal error (phiTot=0)", currentNode, Mesh::Location::Nodes);
    }

    double phi0 = 0.0;
    double dPhi0 = 0.0;
    double dPhi = 0.0;
    double dTheta = 0.0;
    for (Index f = 0; f < numSharedFaces; f++)
    {
        phi0 = phi0 + 0.5 * dPhi;
        if (m_sharedFacesCache[f] == constants::missing::sizetValue)
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
                throw MeshGeometryError("Smoother::ComputeNodeXiEta: Inappropriate fictitious boundary face", currentNode, Mesh::Location::Nodes);
            }
            phi0 = phi0 + 0.5 * dPhi;
            continue;
        }

        const auto numFaceNodes = m_mesh->GetNumFaceEdges(m_sharedFacesCache[f]);
        if (numFaceNodes > Mesh::m_maximumNumberOfEdgesPerNode)
        {
            throw AlgorithmError("Smoother::ComputeNodeXiEta: The number of face nodes is greater than the maximum number of edges per node.");
        }

        dPhi0 = OptimalEdgeAngle(numFaceNodes);
        if (isSquareFace[f])
        {
            auto nextNode = f + static_cast<Index>(2);
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
        const Index nodeIndex = FindIndex(m_mesh->m_facesNodes[m_sharedFacesCache[f]], static_cast<Index>(currentNode));

        // optimal angle
        dTheta = 2.0 * M_PI / static_cast<double>(numFaceNodes);

        // orientation of the face (necessary for folded cells)
        const auto previousNode = NextCircularForwardIndex(nodeIndex, numFaceNodes);
        const auto nextNode = NextCircularBackwardIndex(nodeIndex, numFaceNodes);

        if (m_faceNodeMappingCache[f][nextNode] + 1 == m_faceNodeMappingCache[f][previousNode] ||
            m_faceNodeMappingCache[f][nextNode] - m_faceNodeMappingCache[f][previousNode] == m_mesh->m_nodesNumEdges[currentNode])
        {
            dTheta = -dTheta;
        }

        double aspectRatio = (1.0 - std::cos(dTheta)) / std::sin(std::abs(dTheta)) * std::tan(0.5 * dPhi);
        double radius = std::cos(0.5 * dPhi) / (1.0 - cos(dTheta));

        for (Index n = 0; n < numFaceNodes; n++)
        {
            double theta = dTheta * (static_cast<int>(n) - static_cast<int>(nodeIndex));
            double xip = radius - radius * std::cos(theta);
            double ethap = -radius * std::sin(theta);

            m_xiCache[m_faceNodeMappingCache[f][n]] = xip * std::cos(phi0) - aspectRatio * ethap * std::sin(phi0);
            m_etaCache[m_faceNodeMappingCache[f][n]] = xip * std::sin(phi0) + aspectRatio * ethap * std::cos(phi0);
        }
    }
}

void Smoother::NodeAdministration(Index currentNode)
{
    m_sharedFacesCache.clear();
    m_connectedNodesCache.clear();
    m_connectedNodes[currentNode].clear();

    if (m_mesh->m_nodesNumEdges[currentNode] < 2)
    {
        return;
    }

    // For the currentNode, find the shared faces
    Index newFaceIndex = constants::missing::sizetValue;
    for (Index e = 0; e < m_mesh->m_nodesNumEdges[currentNode]; e++)
    {
        const auto firstEdge = m_mesh->m_nodesEdges[currentNode][e];

        Index secondEdgeIndex = e + 1;
        if (secondEdgeIndex >= m_mesh->m_nodesNumEdges[currentNode])
        {
            secondEdgeIndex = 0;
        }

        const auto secondEdge = m_mesh->m_nodesEdges[currentNode][secondEdgeIndex];
        if (m_mesh->m_edgesNumFaces[firstEdge] == 0 || m_mesh->m_edgesNumFaces[secondEdge] == 0)
        {
            continue;
        }

        // find the face shared by the two edges
        const auto firstFace = std::max(std::min(m_mesh->m_edgesNumFaces[firstEdge], 2U), 1U) - 1U;
        const auto secondFace = std::max(std::min(m_mesh->m_edgesNumFaces[secondEdge], static_cast<Index>(2)), static_cast<Index>(1)) - 1;

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
            newFaceIndex = constants::missing::sizetValue;
        }

        // corner face (already found in the first iteration)
        if (m_mesh->m_nodesNumEdges[currentNode] == 2 &&
            e == 1 &&
            m_mesh->m_nodesTypes[currentNode] == 3 &&
            !m_sharedFacesCache.empty() &&
            m_sharedFacesCache[0] == newFaceIndex)
        {
            newFaceIndex = constants::missing::sizetValue;
        }
        m_sharedFacesCache.emplace_back(newFaceIndex);
    }

    // no shared face found
    if (m_sharedFacesCache.empty())
    {
        return;
    }

    m_connectedNodes[currentNode].emplace_back(currentNode);
    m_connectedNodesCache.emplace_back(currentNode);

    // edge connected nodes
    for (Index e = 0; e < m_mesh->m_nodesNumEdges[currentNode]; e++)
    {
        const auto edgeIndex = m_mesh->m_nodesEdges[currentNode][e];
        const auto node = OtherNodeOfEdge(m_mesh->m_edges[edgeIndex], currentNode);
        m_connectedNodes[currentNode].emplace_back(node);
        m_connectedNodesCache.emplace_back(node);
    }

    // for each face store the positions of the its nodes in the connectedNodes (compressed array)
    if (m_faceNodeMappingCache.size() < m_sharedFacesCache.size())
    {
        ResizeAndFill2DVector(m_faceNodeMappingCache, m_sharedFacesCache.size(), Mesh::m_maximumNumberOfNodesPerFace);
    }
    for (Index f = 0; f < m_sharedFacesCache.size(); f++)
    {
        const auto faceIndex = m_sharedFacesCache[f];
        if (faceIndex == constants::missing::sizetValue)
        {
            continue;
        }

        // find the stencil node position  in the current face
        Index faceNodeIndex = 0;
        const auto numFaceNodes = m_mesh->GetNumFaceEdges(faceIndex);
        for (Index n = 0; n < numFaceNodes; n++)
        {
            if (m_mesh->m_facesNodes[faceIndex][n] == currentNode)
            {
                faceNodeIndex = n;
                break;
            }
        }

        for (Index n = 0; n < numFaceNodes; n++)
        {

            if (faceNodeIndex >= numFaceNodes)
            {
                faceNodeIndex -= numFaceNodes;
            }

            const auto node = m_mesh->m_facesNodes[faceIndex][faceNodeIndex];

            bool isNewNode = true;
            for (Index i = 0; i < m_connectedNodesCache.size(); i++)
            {
                if (node == m_connectedNodesCache[i])
                {
                    isNewNode = false;
                    m_faceNodeMappingCache[f][faceNodeIndex] = static_cast<Index>(i);
                    break;
                }
            }

            if (isNewNode)
            {
                m_connectedNodesCache.emplace_back(node);
                m_faceNodeMappingCache[f][faceNodeIndex] = static_cast<Index>(m_connectedNodesCache.size() - 1);
                m_connectedNodes[currentNode].emplace_back(node);
            }

            // update node index
            faceNodeIndex += 1;
        }
    }

    // update connected nodes (kkc)
    m_numConnectedNodes[currentNode] = static_cast<Index>(m_connectedNodesCache.size());
}

double Smoother::OptimalEdgeAngle(Index numFaceNodes, double theta1, double theta2, bool isBoundaryEdge) const
{
    double angle = M_PI * (1.0 - 2.0 / static_cast<double>(numFaceNodes));

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

void Smoother::Initialize()
{
    // local matrices caches
    m_numConnectedNodes.resize(m_mesh->GetNumNodes());
    std::fill(m_numConnectedNodes.begin(), m_numConnectedNodes.end(), 0);

    m_connectedNodes.resize(m_mesh->GetNumNodes());
    std::fill(m_connectedNodes.begin(), m_connectedNodes.end(), std::vector<Index>(Mesh::m_maximumNumberOfConnectedNodes, 0));

    m_sharedFacesCache.reserve(Mesh::m_maximumNumberOfEdgesPerNode);

    m_connectedNodesCache.reserve(Mesh::m_maximumNumberOfConnectedNodes);

    m_faceNodeMappingCache.resize(Mesh::m_maximumNumberOfConnectedNodes);
    std::fill(m_faceNodeMappingCache.begin(), m_faceNodeMappingCache.end(), std::vector<Index>(Mesh::m_maximumNumberOfNodesPerFace, 0));

    m_xiCache.resize(Mesh::m_maximumNumberOfConnectedNodes, 0.0);
    std::fill(m_xiCache.begin(), m_xiCache.end(), 0.0);

    m_etaCache.resize(Mesh::m_maximumNumberOfConnectedNodes, 0.0);
    std::fill(m_etaCache.begin(), m_etaCache.end(), 0.0);

    m_nodeTopologyMapping.resize(m_mesh->GetNumNodes());
    std::fill(m_nodeTopologyMapping.begin(), m_nodeTopologyMapping.end(), constants::missing::sizetValue);
}

void Smoother::AllocateNodeOperators(Index topologyIndex)
{
    const auto numSharedFaces = m_topologySharedFaces[topologyIndex].size();
    const auto numConnectedNodes = m_topologyConnectedNodes[topologyIndex].size();

    // will reallocate only if necessary
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

void Smoother::SaveNodeTopologyIfNeeded(Index currentNode)
{
    bool isNewTopology = true;
    for (Index topo = 0; topo < m_topologyConnectedNodes.size(); topo++)
    {
        if (m_sharedFacesCache.size() != m_topologySharedFaces[topo].size() || m_connectedNodesCache.size() != m_topologyConnectedNodes[topo].size())
        {
            continue;
        }

        isNewTopology = false;
        for (Index n = 1; n < m_connectedNodesCache.size(); n++)
        {
            const double thetaLoc = std::atan2(m_etaCache[n], m_xiCache[n]);
            const double thetaTopology = std::atan2(m_topologyEta[topo][n], m_topologyXi[topo][n]);
            if (std::abs(thetaLoc - thetaTopology) > m_thetaTolerance)
            {
                isNewTopology = true;
                break;
            }
        }

        if (!isNewTopology)
        {
            m_nodeTopologyMapping[currentNode] = static_cast<Index>(topo);
            break;
        }
    }

    if (isNewTopology)
    {
        m_topologyConnectedNodes.emplace_back(m_connectedNodesCache);
        m_topologySharedFaces.emplace_back(m_sharedFacesCache);
        m_topologyXi.emplace_back(m_xiCache);
        m_topologyEta.emplace_back(m_etaCache);
        m_topologyFaceNodeMapping.emplace_back(m_faceNodeMappingCache);
        m_nodeTopologyMapping[currentNode] = static_cast<Index>(m_topologyConnectedNodes.size() - 1);
    }
}

void Smoother::ComputeJacobian(Index currentNode, std::vector<double>& J) const
{
    const auto currentTopology = m_nodeTopologyMapping[currentNode];
    const auto numNodes = m_topologyConnectedNodes[currentTopology].size();
    if (m_mesh->m_projection == Projection::cartesian)
    {
        J[0] = 0.0;
        J[1] = 0.0;
        J[2] = 0.0;
        J[3] = 0.0;
        for (Index i = 0; i < numNodes; i++)
        {
            J[0] += m_Jxi[currentTopology][i] * m_mesh->m_nodes[m_topologyConnectedNodes[currentTopology][i]].x;
            J[1] += m_Jxi[currentTopology][i] * m_mesh->m_nodes[m_topologyConnectedNodes[currentTopology][i]].y;
            J[2] += m_Jeta[currentTopology][i] * m_mesh->m_nodes[m_topologyConnectedNodes[currentTopology][i]].x;
            J[3] += m_Jeta[currentTopology][i] * m_mesh->m_nodes[m_topologyConnectedNodes[currentTopology][i]].y;
        }
    }
    if (m_mesh->m_projection == Projection::spherical || m_mesh->m_projection == Projection::sphericalAccurate)
    {
        const auto cosFactor = std::cos(m_mesh->m_nodes[currentNode].y * constants::conversion::degToRad);
        J[0] = 0.0;
        J[1] = 0.0;
        J[2] = 0.0;
        J[3] = 0.0;
        for (Index i = 0; i < numNodes; i++)
        {
            J[0] += m_Jxi[currentTopology][i] * m_mesh->m_nodes[m_topologyConnectedNodes[currentTopology][i]].x * cosFactor;
            J[1] += m_Jxi[currentTopology][i] * m_mesh->m_nodes[m_topologyConnectedNodes[currentTopology][i]].y;
            J[2] += m_Jeta[currentTopology][i] * m_mesh->m_nodes[m_topologyConnectedNodes[currentTopology][i]].x * cosFactor;
            J[3] += m_Jeta[currentTopology][i] * m_mesh->m_nodes[m_topologyConnectedNodes[currentTopology][i]].y;
        }
    }
}
