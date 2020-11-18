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
#include <memory>

namespace meshkernel
{
    class Mesh;

    /// <summary>
    /// Orthogonalizion (optimize the aspect ratios) and and mesh smoothing (optimize internal face angles or area).
    /// </summary>
    class Smoother
    {

    public:
        /// <summary>
        /// Mesh ctor
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        explicit Smoother(std::shared_ptr<Mesh> mesh);

        /// @brief Computes the smoother weights
        void Compute();

        /// <summary>
        /// Gets the weight for a certain node and connected node
        /// </summary>
        /// <param name="node"></param>
        /// <param name="connectedNode"></param>
        /// <returns></returns>
        [[nodiscard]] inline auto GetWeight(int node, int connectedNode)
        {
            return m_weights[node][connectedNode];
        }

        /// <summary>
        /// Get the index of the coonected node as assigned by the smoother administration
        /// </summary>
        /// <param name="node"></param>
        /// <param name="connectedNode"></param>
        /// <returns></returns>
        [[nodiscard]] inline auto GetCoonectedNodeIndex(int node, int connectedNode)
        {
            return m_connectedNodes[node][connectedNode];
        }

        /// <summary>
        /// Get number of connected nodes
        /// </summary>
        /// <param name="node"></param>
        /// <returns></returns>
        [[nodiscard]] inline auto GetNumConnectedNodes(int node)
        {
            return m_numConnectedNodes[node];
        }

    private:
        /// @brief Initialize smoother topologies. A topology is determined by how many nodes are connected to the current node.
        ///        There are at maximum mesh.m_numNodes topologies, most likely much less
        void Initialize();

        /// @brief Computes all topologies of the elliptic smoother
        void ComputeTopologies();

        /// @brief Computes all operators of the elliptic smoother
        void ComputeOperators();

        /// <summary>
        /// Compute nodes local coordinates, sice-effects only for sphericalAccurate projection (comp_local_coords)
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool ComputeCoordinates();

        /// @brief Computes the smoother weights from the operators (orthonet_compweights_smooth)
        void ComputeWeights();

        /// Computes operators of the elliptic smoother by node (orthonet_comp_operators)
        /// @param[in] currentNode
        void ComputeOperatorsNode(int currentNode);

        /// @brief Computes m_faceNodeMappingCache, m_sharedFacesCache, m_connectedNodes for the current node, required before computing xi and eta
        /// @param[in] currentNode
        /// @param[out] numSharedFaces
        /// @param[out] numConnectedNodes
        void NodeAdministration(int currentNode,
                                int& numSharedFaces,
                                int& numConnectedNodes);

        /// @brief Compute compute current node xi and eta (orthonet_assign_xieta)
        /// @param[in] currentNode
        /// @param[in] numSharedFaces
        /// @param[in] numConnectedNodes
        void ComputeNodeXiEta(int currentNode,
                              int numSharedFaces,
                              int numConnectedNodes);

        /// <summary>
        /// Compute optimal edge angle
        /// </summary>
        /// <param name="numFaceNodes"></param>
        /// <param name="theta1"></param>
        /// <param name="theta2"></param>
        /// <param name="isBoundaryEdge"></param>
        /// <returns></returns>
        [[nodiscard]] double OptimalEdgeAngle(int numFaceNodes,
                                              double theta1 = -1.0,
                                              double theta2 = -1.0,
                                              bool isBoundaryEdge = false) const;

        /// @brief Allocate smoother operators
        /// @param[in] topologyIndex
        void AllocateNodeOperators(int topologyIndex);

        /// @brief If it is a new topology, save it
        /// @param[in] currentNode
        /// @param[in] numSharedFaces
        /// @param[in] numConnectedNodes
        void SaveNodeTopologyIfNeeded(int currentNode,
                                      int numSharedFaces,
                                      int numConnectedNodes);

        /// @brief Computes local coordinates jacobian from the mapped jacobians m_Jxi and m_Jeta
        /// @param[in] currentNode
        /// @param[out] J
        void ComputeJacobian(int currentNode, std::vector<double>& J) const;

        /// <summary>
        /// Compute the matrix norm
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="matCoefficents"></param>
        /// <returns></returns>
        [[nodiscard]] double MatrixNorm(const std::vector<double>& x,
                                        const std::vector<double>& y,
                                        const std::vector<double>& matCoefficents) const;

        // The mesh to smooth
        std::shared_ptr<Mesh> m_mesh;

        // Smoother weights
        std::vector<std::vector<double>> m_weights;

        // Smoother operators
        std::vector<std::vector<std::vector<double>>> m_Gxi;  // Node to edge xi derivative
        std::vector<std::vector<std::vector<double>>> m_Geta; // Node to edge etha derivative
        std::vector<std::vector<double>> m_Divxi;             // Edge to node xi derivative
        std::vector<std::vector<double>> m_Diveta;            // Edge to node etha derivative
        std::vector<std::vector<std::vector<double>>> m_Az;   // Coefficients to estimate values at cell circumcenters
        std::vector<std::vector<double>> m_Jxi;               // Node to node xi derivative (Jacobian)
        std::vector<std::vector<double>> m_Jeta;              // Node to node eta derivative (Jacobian)
        std::vector<std::vector<double>> m_ww2;               // weights

        // Smoother local caches
        std::vector<int> m_sharedFacesCache;
        std::vector<size_t> m_connectedNodesCache;
        std::vector<std::vector<size_t>> m_faceNodeMappingCache;
        std::vector<double> m_xiCache;
        std::vector<double> m_etaCache;
        std::vector<int> m_boundaryEdgesCache;
        std::vector<double> m_leftXFaceCenterCache;
        std::vector<double> m_leftYFaceCenterCache;
        std::vector<double> m_rightXFaceCenterCache;
        std::vector<double> m_rightYFaceCenterCache;
        std::vector<double> m_xisCache;
        std::vector<double> m_etasCache;

        // Smoother topologies
        int m_numTopologies = 0;
        std::vector<int> m_nodeTopologyMapping;
        std::vector<int> m_numTopologyNodes;
        std::vector<int> m_numTopologyFaces;
        std::vector<std::vector<double>> m_topologyXi;
        std::vector<std::vector<double>> m_topologyEta;
        std::vector<std::vector<int>> m_topologySharedFaces;
        std::vector<std::vector<std::vector<size_t>>> m_topologyFaceNodeMapping;
        std::vector<std::vector<size_t>> m_topologyConnectedNodes;

        std::vector<int> m_numConnectedNodes;              // (nmk2)
        std::vector<std::vector<size_t>> m_connectedNodes; // (kk2)

        // Class variables
        int m_maximumNumConnectedNodes = 0;
        int m_maximumNumSharedFaces = 0;

        // nodes with errors
        std::vector<double> m_nodeXErrors;
        std::vector<double> m_nodeYErrors;
        std::vector<int> m_nodeErrorCode;

        static constexpr int m_topologyInitialSize = 10;
        static constexpr double m_thetaTolerance = 1e-4;
    };
} // namespace meshkernel
