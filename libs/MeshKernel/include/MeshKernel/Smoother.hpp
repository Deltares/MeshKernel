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

#pragma once

#include <vector>

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Definitions.hpp"

namespace meshkernel
{
    class Mesh2D;

    /// <summary>
    /// Orthogonalizion (optimize the aspect ratios) and mesh smoothing (optimize internal face angles or area).
    /// </summary>
    class Smoother
    {

    public:
        /// @brief Mesh2D constructor
        /// @param[in] mesh Mesh
        /// @param[in] nodeType Type of node in the mesh, may be different from the node type contained in the Mesh.
        explicit Smoother(const Mesh2D& mesh, const std::vector<MeshNodeType>& nodeType);

        /// @brief Computes the smoother weights
        void Compute();

        /// @brief Gets the weight for a certain node and connected node
        [[nodiscard]] auto GetWeight(UInt node, int connectedNode) const
        {
            return m_weights[node][connectedNode];
        }

        /// @brief Get the index of the connected node as assigned by the smoother administration
        /// @brief node
        /// @brief connectedNode
        /// @returns
        [[nodiscard]] auto GetConnectedNodeIndex(UInt node, int connectedNode) const
        {
            return m_connectedNodes[node][connectedNode];
        }

        /// @brief Get number of connected nodes
        /// @brief node
        /// @returns
        [[nodiscard]] auto GetNumConnectedNodes(UInt node) const
        {
            return m_numConnectedNodes[node];
        }

    private:
        /// @brief Contains internal edge angle information
        struct InternalAngleData
        {
            UInt numSquaredTriangles = 0;     ///< Number of  squared triangles
            UInt numTriangles = 0;            ///< Number of triangles sharing node
            double phiSquaredTriangles = 0.0; ///< Sum of optimal edge angles for squared triangles
            double phiQuads = 0.0;            ///< Sum of optimal edge angles for quadrilaterals
            double phiTriangles = 0.0;        ///< Sum of optimal edge angles for triangles
            double phiTot = 0.0;              ///< Sum of optimal edge angles for all element shapes
        };

        /// @brief Compute number and sum of all angles interior angles
        void ComputeInternalAngle(const UInt currentNode,
                                  const UInt numSharedFaces,
                                  const std::vector<double>& thetaSquare,
                                  const std::vector<bool>& isSquareFace,
                                  InternalAngleData& internalAngleData,
                                  UInt& numNonStencilQuad);

        /// @brief Update theta squared for interior faces.
        void UpdateThetaForInteriorFaces(const UInt numSharedFaces, std::vector<double>& thetaSquare);

        /// @brief Updata xi- and eta-caches for shared faces of a node.
        void UpdateXiEtaForSharedFace(const UInt currentNode,
                                      const UInt currentFace,
                                      const UInt numFaceNodes,
                                      const double dPhi,
                                      const double phi0);

        /// @brief Compute the optimal angle for all theshared faces attached to a node.
        void ComputeOptimalAngleForSharedFaces(const UInt currentNode,
                                               const UInt numSharedFaces,
                                               const MeshNodeType nodeType,
                                               const std::vector<double>& thetaSquare,
                                               const std::vector<bool>& isSquareFace,
                                               const double mu,
                                               const double muSquaredTriangles,
                                               const double muTriangles);

        /// @brief Initialize smoother topologies. A topology is determined by how many nodes are connected to the current node.
        ///        There are at maximum mesh.m_numNodes topologies, most likely much less
        void Initialize();

        /// @brief Computes all topologies of the elliptic smoother
        void ComputeTopologies();

        /// @brief Computes all operators of the elliptic smoother
        void ComputeOperators();

        /// @brief Computes the smoother weights from the operators (orthonet_compweights_smooth)
        void ComputeWeights();

        /// @brief Compute elliptic smoother operators coefficients for boundary nodes
        std::tuple<double, double> ComputeOperatorsForBoundaryNode(const UInt f, const UInt faceLeftIndex, const UInt currentTopology);

        /// @brief Compute elliptic smoother operators coefficients for interior nodes
        std::tuple<double, double, double, double> ComputeOperatorsForInteriorNode(const UInt f,
                                                                                   const UInt edgeIndex,
                                                                                   const UInt faceLeftIndex,
                                                                                   const UInt faceRightIndex,
                                                                                   const UInt currentTopology);

        /// @brief Compute the node to edge derivatives, gxi and geta
        void ComputeNodeEdgeDerivative(const UInt f,
                                       const UInt edgeIndex,
                                       const UInt currentTopology,
                                       const UInt faceLeftIndex,
                                       const UInt faceRightIndex,
                                       const double facxiL,
                                       const double facetaL,
                                       const double facxiR,
                                       const double facetaR);

        /// Computes operators of the elliptic smoother by node (orthonet_comp_operators)
        /// @param[in] currentNode
        /// @param[in] nodeType Node type of current node
        void ComputeOperatorsNode(UInt currentNode, const MeshNodeType nodeType);

        /// @brief Computes m_faceNodeMappingCache, m_sharedFacesCache, m_connectedNodes for the current node, required before computing xi and eta
        /// @param[in] currentNode
        void NodeAdministration(UInt currentNode);

        /// @brief Compute compute current node xi and eta (orthonet_assign_xieta)
        /// @param[in] currentNode
        void ComputeNodeXiEta(UInt currentNode);

        /// @brief Compute optimal edge angle
        /// @brief numFaceNodes
        /// @brief theta1
        /// @brief theta2
        /// @brief isBoundaryEdge
        /// @returns
        [[nodiscard]] double OptimalEdgeAngle(UInt numFaceNodes,
                                              double theta1 = -1.0,
                                              double theta2 = -1.0,
                                              bool isBoundaryEdge = false) const;

        /// @brief Allocate smoother operators
        /// @param[in] topologyIndex
        void AllocateNodeOperators(UInt topologyIndex);

        /// @brief If it is a new topology, save it
        /// @param[in] currentNode
        void SaveNodeTopologyIfNeeded(UInt currentNode);

        /// @brief Computes local coordinates jacobian from the mapped jacobians m_Jxi and m_Jeta
        /// @param[in] currentNode
        /// @param[out] J
        void ComputeJacobian(UInt currentNode, std::vector<double>& J) const;

        /// @brief Compute the coefficients to estimate values at cell circumcenters, filling m_Az.
        void ComputeCellCircumcentreCoefficients(const UInt currentNode, const UInt currentTopology, const MeshNodeType nodeType);

        /// @brief Compute the segment gradients
        void ComputeNodeToNodeGradients(const UInt currentNode, const UInt currentTopology);

        /// @brief Compute the weights for the Laplacian smoother
        void ComputeLaplacianSmootherWeights(const UInt currentNode, const UInt currentTopology);

        // The mesh to smooth
        const Mesh2D& m_mesh; ///< A reference to mesh

        const std::vector<MeshNodeType>& m_nodeType; ///< Node types

        // Smoother weights
        std::vector<std::vector<double>> m_weights; ///< Weights

        // Smoother operators
        std::vector<std::vector<std::vector<double>>> m_Gxi;  ///< Node to edge xi derivative
        std::vector<std::vector<std::vector<double>>> m_Geta; ///< Node to edge etha derivative
        std::vector<std::vector<double>> m_Divxi;             ///< Edge to node xi derivative
        std::vector<std::vector<double>> m_Diveta;            ///< Edge to node etha derivative
        std::vector<std::vector<std::vector<double>>> m_Az;   ///< Coefficients to estimate values at cell circumcenters
        std::vector<std::vector<double>> m_Jxi;               ///< Node to node xi derivative (Jacobian)
        std::vector<std::vector<double>> m_Jeta;              ///< Node to node eta derivative (Jacobian)
        std::vector<std::vector<double>> m_ww2;               ///< weights

        // Smoother local caches
        std::vector<UInt> m_sharedFacesCache;                  ///< Cache for shared faces
        std::vector<UInt> m_connectedNodesCache;               ///< Cache for connected nodes
        std::vector<std::vector<UInt>> m_faceNodeMappingCache; ///< Cache for face node mapping
        std::vector<double> m_xiCache;                         ///< Cache for xi
        std::vector<double> m_etaCache;                        ///< Cache for eta
        std::vector<double> m_leftXFaceCenterCache;            ///< Cache for left x face center
        std::vector<double> m_leftYFaceCenterCache;            ///< Cache for left y face center
        std::vector<double> m_rightXFaceCenterCache;           ///< Cache for right x face center
        std::vector<double> m_rightYFaceCenterCache;           ///< Cache for right y face center
        std::vector<double> m_xisCache;                        ///< Cache for xis
        std::vector<double> m_etasCache;                       ///< Cache for etas

        // Smoother topologies
        std::vector<UInt> m_nodeTopologyMapping;                               ///< Node topology mapping
        std::vector<std::vector<double>> m_topologyXi;                         ///< Topology xi
        std::vector<std::vector<double>> m_topologyEta;                        ///< Topology eta
        std::vector<std::vector<UInt>> m_topologySharedFaces;                  ///< Topology shared faces
        std::vector<std::vector<std::vector<UInt>>> m_topologyFaceNodeMapping; ///< Topology face node mapping
        std::vector<std::vector<UInt>> m_topologyConnectedNodes;               ///< Topology connected nodes

        std::vector<UInt> m_numConnectedNodes;           ///< Number of connected nodes (nmk2)
        std::vector<std::vector<UInt>> m_connectedNodes; ///< Connected nodes (kk2)

        // Class variables
        UInt m_maximumNumConnectedNodes = 0; ///< Maximum number of connected nodes
        UInt m_maximumNumSharedFaces = 0;    ///< Maximum number of shared faces

        static constexpr int m_topologyInitialSize = 10; ///< Initial size of topology vectors
        static constexpr double m_thetaTolerance = 1e-4; ///< Tolerance
    };
} // namespace meshkernel
