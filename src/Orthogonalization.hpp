#pragma once

#include <vector>
#include "Entities.hpp"
#include "Mesh.hpp"
#include "OrthogonalizationParametersNative.hpp"
#include "LandBoundaries.hpp"
#include "Polygons.hpp"

namespace GridGeom
{
    class Orthogonalization
    {
    public:
        
        bool Set(const Mesh& mesh,
            int& isTriangulationRequired,
            int& isAccountingForLandBoundariesRequired,
            int& projectToLandBoundaryOption,
            GridGeomApi::OrthogonalizationParametersNative& orthogonalizationParametersNative,
            std::vector<Point>& polygon,
            std::vector<Point>& landBoundaries);

        bool Iterate(Mesh& mesh);

        bool PrapareOuterIteration(const Mesh& mesh);

        bool InnerIteration(Mesh& mesh);

        bool FinalizeOuterIteration(Mesh& mesh);

        bool GetOrthogonality(const Mesh& mesh, double* orthogonality);

        bool GetSmoothness(const Mesh& mesh, double* smoothness);

    private:

        /// orthonet_project_on_boundary: project boundary-nodes back to the boundary of an original net
        bool ProjectOnBoundary(Mesh& mesh);

        /// snapping nodes to land boundaries
        bool SnapToLandBoundary(Mesh& mesh, const LandBoundaries& landBoundaries);

        /// orthonet_compweights_smooth: inverse - mapping elliptic smoother
        bool ComputeWeightsSmoother(const Mesh& mesh);

        bool ComputeSmootherOperators(const Mesh& mesh);

        /// orthonet_comp_operators, compute coefficient matrix G of gradient at link, compute coefficientmatrix Div of gradient in node, compute coefficientmatrix Az of cell - center in cell
        bool ComputeOperatorsNode(const Mesh& mesh, const int currentNode, const std::size_t& numConnectedNodes, const std::vector<std::size_t>& connectedNodes, const std::size_t& numSharedFaces, const std::vector<int>& sharedFaces,
            const std::vector<double>& xi, const std::vector<double>& eta, const std::vector<std::vector<std::size_t>>& faceNodeMapping);

        /// orthonet_assign_xieta: assign xiand eta to all nodes in the stencil
        bool ComputeXiEta(const Mesh& mesh, int currentNode, const std::vector<int>& sharedFaces, const int& numSharedFaces, const std::vector<std::size_t>& connectedNodes,
            const std::size_t& numConnectedNodes, const std::vector<std::vector<std::size_t>>& faceNodeMapping, std::vector<double>& xi, std::vector<double>& eta);

        bool ComputeFacesNumEdges(const Mesh& mesh);

        ///  computes the shared faces and the connected nodes of a stencil node and the faceNodeMapping in the connectedNodes array for each shared face.
        bool OrthogonalizationAdministration(const Mesh& mesh, const int currentNode, std::vector<int>& sharedFaces, int& numSharedFaces, std::vector<std::size_t>& connectedNodes, int& numConnectedNodes, std::vector<std::vector<std::size_t>>& faceNodeMapping);

        double OptimalEdgeAngle(int numFaceNodes, double theta1 = -1.0, double theta2 = -1.0, bool isBoundaryEdge = false);

        ///  orthonet_compute_aspect: compute link - based aspect ratios
        bool AspectRatio(const Mesh& mesh);

        /// orthonet_compweights: compute weights wwand right - hand side rhs in orthogonalizer
        bool ComputeWeightsOrthogonalizer(const Mesh& mesh);

        double MatrixNorm(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& matCoefficents);

        bool InitializeTopologies(const Mesh& mesh);

        bool AllocateNodeOperators(const int topologyIndex);

        /// save only the unique topologies
        bool SaveTopology(int currentNode, const std::vector<int>& sharedFaces, int numSharedFaces, const std::vector<std::size_t>& connectedNodes, int numConnectedNodes,
            const std::vector<std::vector<std::size_t>>& faceNodeMapping, const std::vector<double>& xi, const std::vector<double>& eta);

        bool ComputeIncrements(const Mesh& mesh);

        bool AllocateCaches(const Mesh& mesh);

        bool DeallocateCaches();

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
        std::vector<int> m_numConnectedNodes;                            // nmk2, determined from local node administration
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

        // the land boundary
        LandBoundaries m_landBoundaries;

        // the polygons
        Polygons m_polygons;

        int m_isTriangulationRequired;

        int m_isAccountingForLandBoundariesRequired;
        
        int m_projectToLandBoundaryOption;
    };
}