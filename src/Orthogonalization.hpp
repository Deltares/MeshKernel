#pragma once

#include <vector>
#include "Entities.hpp"
#include "Mesh.hpp"
#include "OrthogonalizationParametersNative.hpp"
#include "GeometryListNative.hpp"

namespace GridGeom
{
    class Orthogonalization
    {
    public:
        
        bool initialize(const Mesh& mesh,
            int& isTriangulationRequired,
            int& isAccountingForLandBoundariesRequired,
            int& projectToLandBoundaryOption,
            GridGeomApi::OrthogonalizationParametersNative& orthogonalizationParametersNative,
            GridGeomApi::GeometryListNative& geometryListNativePolygon,
            GridGeomApi::GeometryListNative& geometryListNativeLandBoundaries);

        bool iterate(Mesh& mesh);

        bool prapareOuterIteration(const Mesh& mesh);

        bool innerIteration(Mesh& mesh);

        bool finalizeOuterIteration(Mesh& mesh);

        /// @brief retrive orthogonality
        bool getOrthogonality(const Mesh& mesh, double* orthogonality);

        /// @brief retrive smoothness
        bool getSmoothness(const Mesh& mesh, double* smoothness);
    
    private:

        /// @brief orthonet_project_on_boundary: project boundary-nodes back to the boundary of an original net
        bool projectOnBoundary(Mesh& mesh);

        /// @brief orthonet_compweights_smooth: inverse - mapping elliptic smoother
        bool computeWeightsSmoother(const Mesh& mesh);

        bool computeSmootherOperators(const Mesh& mesh);

        /// @brief orthonet_comp_operators, compute coefficient matrix G of gradient at link, compute coefficientmatrix Div of gradient in node, compute coefficientmatrix Az of cell - center in cell
        bool computeOperatorsNode(const Mesh& mesh, const int currentNode, const std::size_t& numConnectedNodes, const std::vector<std::size_t>& connectedNodes, const std::size_t& numSharedFaces, const std::vector<int>& sharedFaces,
            const std::vector<double>& xi, const std::vector<double>& eta, const std::vector<std::vector<std::size_t>>& faceNodeMapping);

        /// @brief orthonet_assign_xieta: assign xiand eta to all nodes in the stencil
        bool computeXiEta(const Mesh& mesh, int currentNode, const std::vector<int>& sharedFaces, const int& numSharedFaces, const std::vector<std::size_t>& connectedNodes,
            const std::size_t& numConnectedNodes, const std::vector<std::vector<std::size_t>>& faceNodeMapping, std::vector<double>& xi, std::vector<double>& eta);

        bool computeFacesNumEdges(const Mesh& mesh);

        /// @brief  computes the shared faces and the connected nodes of a stencil node and the faceNodeMapping in the connectedNodes array for each shared face.
        bool orthogonalizationAdministration(const Mesh& mesh, const int currentNode, std::vector<int>& sharedFaces, int& numSharedFaces, std::vector<std::size_t>& connectedNodes, int& numConnectedNodes, std::vector<std::vector<std::size_t>>& faceNodeMapping);

        double optimalEdgeAngle(int numFaceNodes, double theta1 = -1.0, double theta2 = -1.0, bool isBoundaryEdge = false);

        /// @brief  orthonet_compute_aspect: compute link - based aspect ratios
        bool aspectRatio(const Mesh& mesh);

        /// @brief orthonet_compweights: compute weights wwand right - hand side rhs in orthogonalizer
        bool computeWeightsOrthogonalizer(const Mesh& mesh);

        double matrixNorm(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& matCoefficents);

        bool initializeTopologies(const Mesh& mesh);

        bool allocateNodeOperators(const int topologyIndex);

        /// @brief save only the unique topologies
        bool saveTopology(int currentNode, const std::vector<int>& sharedFaces, int numSharedFaces, const std::vector<std::size_t>& connectedNodes, int numConnectedNodes,
            const std::vector<std::vector<std::size_t>>& faceNodeMapping, const std::vector<double>& xi, const std::vector<double>& eta);

        bool computeIncrements(const Mesh& mesh);

        bool allocateCaches(const Mesh& mesh);

        bool deallocateCaches();

        enum class NodeTypes
        {
            internalNode,
            onRing,
            cornerNode,
            hangingNode,
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
        std::vector<std::vector<std::vector<std::size_t>>> m_topologyFaceNodeMapping;
        std::vector < std::vector<std::size_t>>  m_topologyConnectedNodes;

        std::vector<double> m_aspectRatios;
        std::vector<std::vector<double>> m_ww2Global;
        std::vector<std::size_t> m_numConnectedNodes;                            // nmk2, determined from local node administration
        std::vector<std::vector<std::size_t>> m_connectedNodes;                  // kk2, determined from local node administration
        std::vector<int> m_localCoordinates;                                     // iloc
                                                                                 
        // run-time options                                                      
        bool m_keepCircumcentersAndMassCenters = false;                          
        double m_orthogonalizationToSmoothingFactor = 0.975;                     // Factor(0. <= ATPF <= 1.) between grid smoothing and grid ortho resp.
        double m_orthogonalizationToSmoothingFactorBoundary = 1.0;               // minimum ATPF on the boundary
        double m_smoothorarea = 1.0;                                             // Factor between smoother(1.0) and area - homogenizer(0.0)
        int m_orthogonalizationOuterIterations = 2;
        int m_orthogonalizationBoundaryIterations = 25;
        int m_orthogonalizationInnerIterations = 25;

        static constexpr int m_topologyInitialSize = 10;
        static constexpr double m_thetaTolerance = 1e-4;

        std::vector<std::vector<double>>  m_weights;
        std::vector<std::vector<double>>  m_rightHandSide;
        std::vector<int> m_faceNumNodes;                                         //number of face nodes

        int m_maximumNumConnectedNodes = 0;
        int m_maximumNumSharedFaces = 0;
        std::size_t m_maxNumNeighbours;
        std::vector< std::vector<int>> m_nodesNodes;                              //node neighbours 

        //local caches (avoid re-allocation)
        std::vector<int> m_boundaryEdges;
        std::vector<double> m_leftXFaceCenter;
        std::vector<double> m_leftYFaceCenter;
        std::vector<double> m_rightXFaceCenter;
        std::vector<double> m_rightYFaceCenter;
        std::vector<double> m_xis;
        std::vector<double> m_etas;

        // orthogonalization iterations
        std::vector<std::vector<double>> m_ww2x;
        std::vector<std::vector<double>> m_ww2y;
        std::vector<Point> m_orthogonalCoordinates;
        std::vector<int> m_nearestPoints;
        std::vector<Point> m_originalNodes;

        std::vector<int> m_k1;
        std::vector<double> m_wwx;
        std::vector<double> m_wwy;
        std::vector<double> m_increments;
        std::vector<double> m_rightHandSideCache;
        std::vector<int> m_startCacheIndex;
        std::vector<int> m_endCacheIndex;
        int m_cacheSize = 0;
        
        double m_mumax;
        double m_mu;

        // nodes with errors
        std::vector<double> m_nodeXErrors;
        std::vector<double> m_nodeYErrors;
        std::vector<int>    m_nodeErrorCode;

    };
}