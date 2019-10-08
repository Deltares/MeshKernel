#pragma once

#define _USE_MATH_DEFINES
#include <vector>
#include "Entities.hpp"

namespace GridGeom
{
    template<typename Mesh>
    class Orthogonalization
    {
    public:
        
        bool initialize(const Mesh& mesh);
        bool iterate(Mesh& mesh);

    private:

        typedef typename Mesh::Operations Operations;

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
        bool projectOnBoundary(Mesh& mesh, std::vector<int>& nearestPoints, const std::vector<GridGeom::Point>& originalNodes);

        // orthonet_compweights_smooth
        // inverse - mapping elliptic smoother
        // computes weight ww in :
        // sum_kk ww(kk, k0)* x1(kk2(kk, k0)) = 0
        // sum_kk ww(kk, k0) * y1(kk2(kk, k0)) = 0
        bool computeWeightsSmoother(const Mesh& mesh);

        bool computeSmootherOperators(const Mesh& mesh);

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
            const std::vector<std::vector<size_t>>& faceNodeMapping);

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
            std::vector<double>& eta);

        bool computeFacesNumEdges(const Mesh& mesh);

        // computes the shared faces and the connected nodes of a stencil node and the faceNodeMapping in the connectedNodes array for each shared face.
        bool orthogonalizationAdministration(const Mesh& mesh, const int currentNode, std::vector<int>& sharedFaces, int& numSharedFaces, std::vector<size_t>& connectedNodes, int& numConnectedNodes, std::vector<std::vector<size_t>>& faceNodeMapping);

        double optimalEdgeAngle(const int numFaceNodes, const double theta1 = -1, const double theta2 = -1, bool isBoundaryEdge = false);

        // compute link - based aspect ratios
        // orthonet_compute_aspect
        bool aspectRatio(const Mesh& mesh);

        //MAKENETNODESCODING
        bool classifyNodes(const Mesh& mesh);

        //orthonet_compweights
        // compute weights wwand right - hand side rhs in orthogonizer :
        // sum_kk ww(kk, k0)* (x1(kk1(kk, k0)) - x1(k0)) = rhs(1, k0)
        // sum_kk ww(kk, k0) * (y1(kk1(kk, k0)) - y1(k0)) = rhs(2, k0)
        bool computeWeightsOrthogonalizer(const Mesh& mesh);

        double matrixNorm(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& matCoefficents);

        bool initializeTopologies(const Mesh& mesh);

        bool allocateNodeOperators(const int topologyIndex);

        // save only the unique topologies
        bool saveTopology(const int currentNode,
            const std::vector<int>& sharedFaces,
            const int numSharedFaces,
            const std::vector<size_t>& connectedNodes,
            const int numConnectedNodes,
            const std::vector<std::vector<size_t>>& faceNodeMapping,
            const std::vector<double>& xi,
            const std::vector<double>& eta);
    };
}