#pragma once
#include <vector>

namespace GridGeom
{
    class Mesh;

    /// <summary>
    /// Orthogonalizion (optimize the aspect ratios) and and mesh smoothing (optimize internal face angles or area).
    /// </summary>
    class Smoother
    {
        

    public:

        /// <summary>
        /// Ctor
        /// </summary>
        /// <returns></returns>
        Smoother();
        
        /// <summary>
        /// Computes the smoother weights
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool Compute(const Mesh& mesh);

        
        /// <summary>
        /// Gets the smoother weight for a certain node and connected node
        /// </summary>
        /// <param name="node"></param>
        /// <param name="connectedNode"></param>
        /// <returns></returns>
        inline double GetWeight(int node, int connectedNode)
        {
            return m_weights[node][connectedNode];
        }

        /// <summary>
        /// Get the index of the coonected node as assigned by the smoother administration
        /// </summary>
        /// <param name="node"></param>
        /// <param name="connectedNode"></param>
        /// <returns></returns>
        inline std::size_t GetCoonectedNodeIndex(int node, int connectedNode)
        {
            return m_connectedNodes[node][connectedNode];
        }

        /// <summary>
        /// Get number of connected nodes
        /// </summary>
        /// <param name="node"></param>
        /// <returns></returns>
        inline int GetNumConnectedNodes(int node)
        {
            return m_numConnectedNodes[node];
        }

    private:
        /// <summary>
        /// Initialize smoother topologies. A topology is determined by how many nodes are connected to the current node.
        /// There are at maximum mesh.m_numNodes topologies, most likeley much less
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool Initialize(const Mesh& mesh);

        /// <summary>
        /// Computes all topologies of the elliptic smoother 
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool ComputeTopologies(const Mesh& mesh);

        /// <summary>
        /// Computes all operators of the elliptic smoother 
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool ComputeOperators(const Mesh& mesh);

        /// <summary>
        /// Compute nodes local coordinates, sice-effects only for sphericalAccurate projection (comp_local_coords)
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool ComputeCoordinates(const Mesh& mesh);

        /// <summary>
        /// Computes the smoother weights from the operators (orthonet_compweights_smooth)
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool ComputeWeights(const Mesh& mesh);

        /// <summary>
        /// Computes operators of the elliptic smoother by node (orthonet_comp_operators)
        /// </summary>
        /// <param name="mesh"></param>
        /// <param name="currentNode"></param>
        /// <returns></returns>
        bool ComputeOperatorsNode(const Mesh& mesh,
                                          int currentNode);

        /// <summary>
        /// Computes m_faceNodeMappingCache, m_sharedFacesCache, m_connectedNodes for the current node, required before computing xi and eta
        /// </summary>
        /// <param name="mesh"></param>
        /// <param name="currentNode"></param>
        /// <param name="numSharedFaces"></param>
        /// <param name="numConnectedNodes"></param>
        /// <returns></returns>
        bool NodeAdministration(const Mesh& mesh,
            const int currentNode,
            int& numSharedFaces,
            int& numConnectedNodes);

        /// <summary>
        /// Compute compute current node xi and eta (orthonet_assign_xieta)
        /// </summary>
        /// <param name="mesh"></param>
        /// <param name="currentNode"></param>
        /// <param name="numSharedFaces"></param>
        /// <param name="numConnectedNodes"></param>
        /// <returns></returns>
        bool ComputeNodeXiEta(const Mesh& mesh, 
                                      int currentNode, 
                                      const int& numSharedFaces, 
                                      const int& numConnectedNodes);

        /// <summary>
        /// Compute optimal edge angle
        /// </summary>
        /// <param name="numFaceNodes"></param>
        /// <param name="theta1"></param>
        /// <param name="theta2"></param>
        /// <param name="isBoundaryEdge"></param>
        /// <returns></returns>
        double OptimalEdgeAngle(int numFaceNodes, 
                                double theta1 = -1.0, 
                                double theta2 = -1.0, 
                                bool isBoundaryEdge = false);

        /// <summary>
        /// Allocate smoother operators
        /// </summary>
        /// <param name="topologyIndex"></param>
        /// <returns></returns>
        bool AllocateNodeOperators(int topologyIndex);

        /// <summary>
        /// If it is a new topology, save it
        /// </summary>
        /// <param name="currentNode"></param>
        /// <param name="numSharedFaces"></param>
        /// <param name="numConnectedNodes"></param>
        /// <returns></returns>
        bool SaveNodeTopologyIfNeeded(int currentNode, 
                                      int numSharedFaces, 
                                      int numConnectedNodes);


        /// <summary>
        /// Computes local coordinates jacobian from the mapped jacobians m_Jxi and m_Jeta
        /// </summary>
        /// <param name="currentNode"></param>
        /// <param name="mesh"></param>
        /// <param name="J"></param>
        /// <returns></returns>
        bool ComputeJacobian(int currentNode, const Mesh& mesh, std::vector<double>& J) const;

        /// <summary>
        /// Compute the matrix norm
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="matCoefficents"></param>
        /// <returns></returns>
        double MatrixNorm(const std::vector<double>& x,
                          const std::vector<double>& y,
                          const std::vector<double>& matCoefficents);

        
        // Smoother weights
        std::vector<std::vector<double>>                   m_weights;

        // Smoother operators
        std::vector<std::vector<std::vector<double>>>      m_Gxi;                    // Node to edge xi derivative
        std::vector<std::vector<std::vector<double>>>      m_Geta;                   // Node to edge etha derivative
        std::vector<std::vector<double>>                   m_Divxi;                  // Edge to node xi derivative
        std::vector<std::vector<double>>                   m_Diveta;                 // Edge to node etha derivative
        std::vector<std::vector<std::vector<double>>>      m_Az;                     // Coefficents to estimate values at cell circumcenters                                                                         
        std::vector<std::vector<double>>                   m_Jxi;                    // Node to node xi derivative (Jacobian)
        std::vector<std::vector<double>>                   m_Jeta;                   // Node to node eta derivative (Jacobian)
        std::vector<std::vector<double>>                   m_ww2;                    // weights
           
        // Smoother local caches
        std::vector<int>                                   m_sharedFacesCache;
        std::vector<std::size_t>                           m_connectedNodesCache;
        std::vector<std::vector<std::size_t>>              m_faceNodeMappingCache;
        std::vector<double>                                m_xiCache;
        std::vector<double>                                m_etaCache;
        std::vector<int>                                   m_boundaryEdgesCache;
        std::vector<double>                                m_leftXFaceCenterCache;
        std::vector<double>                                m_leftYFaceCenterCache;
        std::vector<double>                                m_rightXFaceCenterCache;
        std::vector<double>                                m_rightYFaceCenterCache;
        std::vector<double>                                m_xisCache;
        std::vector<double>                                m_etasCache;

        // Smoother topologies
        int m_numTopologies = 0;                           
        std::vector<int>                                   m_nodeTopologyMapping;
        std::vector<int>                                   m_numTopologyNodes;
        std::vector<int>                                   m_numTopologyFaces;
        std::vector<std::vector<double>>                   m_topologyXi;
        std::vector<std::vector<double>>                   m_topologyEta;
        std::vector<std::vector<int>>                      m_topologySharedFaces;
        std::vector<std::vector<std::vector<std::size_t>>> m_topologyFaceNodeMapping;
        std::vector < std::vector<std::size_t>>            m_topologyConnectedNodes;

        std::vector<int>                                   m_numConnectedNodes;        // (nmk2)
        std::vector<std::vector<std::size_t>>              m_connectedNodes;           // (kk2)

        // Class variables
        int m_maximumNumConnectedNodes = 0;
        int m_maximumNumSharedFaces = 0;

        // nodes with errors
        std::vector<double>                                m_nodeXErrors;
        std::vector<double>                                m_nodeYErrors;
        std::vector<int>                                   m_nodeErrorCode;

        static constexpr int m_topologyInitialSize = 10;
        static constexpr double m_thetaTolerance = 1e-4;
    };
}