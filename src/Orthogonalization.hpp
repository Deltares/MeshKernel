#pragma once

#include <vector>
#include "LandBoundaries.hpp"
#include "Polygons.hpp"
#include "OrthogonalizationParametersNative.hpp"

namespace GridGeom
{
    struct Point;
    enum class Projections;
    class Mesh;

    /// <summary>
    /// Orthogonalizion (optimize the aspect ratios) and and mesh smoothing (optimize internal face angles or area).
    /// </summary>
    class Orthogonalization
    {

    public:
        
        /// <summary>
        /// Set algorithm parameters
        /// </summary>
        /// <param name="mesh"></param>
        /// <param name="isTriangulationRequired"></param>
        /// <param name="isAccountingForLandBoundariesRequired"></param>
        /// <param name="projectToLandBoundaryOption"></param>
        /// <param name="orthogonalizationParametersNative"></param>
        /// <param name="polygon"></param>
        /// <param name="landBoundaries"></param>
        /// <returns></returns>
        bool Set( Mesh& mesh,
                  int& isTriangulationRequired,
                  int& isAccountingForLandBoundariesRequired,
                  int& projectToLandBoundaryOption,
                  GridGeomApi::OrthogonalizationParametersNative& orthogonalizationParametersNative,
                  const Polygons& polygon,
                  std::vector<Point>& landBoundaries );

        /// <summary>
        /// Executes the algorithm
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool Compute(Mesh& mesh);

        /// <summary>
        /// Prepares the outer iteration, calculates orthogonalizer and smoother coefficents
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool PrapareOuterIteration(const Mesh& mesh);

        /// <summary>
        /// Performs an inner iteration (update of node positions)
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool InnerIteration(Mesh& mesh);

        /// <summary>
        /// Performs an outer iteration (re-computes the operators)
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool FinalizeOuterIteration(Mesh& mesh);

        /// <summary>
        /// Gets the orthogonality values (innerproduct edges and face circumcenter connecting segments)
        /// </summary>
        /// <param name="mesh"></param>
        /// <param name="orthogonality"></param>
        /// <returns></returns>
        bool GetOrthogonality(const Mesh& mesh, double* orthogonality);

        /// <summary>
        /// Gets the smoothness values (face area ratios)
        /// </summary>
        /// <param name="mesh"></param>
        /// <param name="smoothness"></param>
        /// <returns></returns>
        bool GetSmoothness(const Mesh& mesh, double* smoothness);

    private:

        /// <summary>
        /// Computes the aspect ratio of each edge (orthonet_compute_aspect)
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool AspectRatio(const Mesh& mesh);

        /// <summary>
        /// Computes orthogonalizer weights equation 3.10 of dflowfm technical reference manual (orthonet_compweights)
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool ComputeWeightsAndRhsOrthogonalizer(const Mesh& mesh);

        /// <summary>
        /// Initialize smoother topologies. A topology is determined by how many nodes are connected to the current node.
        /// There are at maximum mesh.m_numNodes topologies, most likeley much less
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool InitializeSmoother(const Mesh& mesh);

        /// <summary>
        /// Computes all topologies of the elliptic smoother 
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool ComputeSmootherTopologies(const Mesh& mesh);

        /// <summary>
        /// Computes all operators of the elliptic smoother 
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool ComputeSmootherOperators(const Mesh& mesh);

        /// <summary>
        /// Compute nodes local coordinates, sice-effects only for sphericalAccurate projection (comp_local_coords)
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool ComputeLocalCoordinates(const Mesh& mesh);

        /// <summary>
        /// Computes the smoother weights from the operators (orthonet_compweights_smooth)
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool ComputeSmootherWeights(const Mesh& mesh);

        /// <summary>
        /// Computes operators of the elliptic smoother by node (orthonet_comp_operators)
        /// </summary>
        /// <param name="mesh"></param>
        /// <param name="currentNode"></param>
        /// <returns></returns>
        bool ComputeSmootherOperatorsNode(const Mesh& mesh,
                                          int currentNode);

        /// <summary>
        /// Computes m_faceNodeMappingCache, m_sharedFacesCache, m_connectedNodes for the current node, required before computing xi and eta
        /// </summary>
        /// <param name="mesh"></param>
        /// <param name="currentNode"></param>
        /// <param name="numSharedFaces"></param>
        /// <param name="numConnectedNodes"></param>
        /// <returns></returns>
        bool SmootherNodeAdministration(const Mesh& mesh,
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
        bool SmootherComputeNodeXiEta(const Mesh& mesh, 
                                      int currentNode, 
                                      const int& numSharedFaces, 
                                      const int& numConnectedNodes);



        /// <summary>
        /// Project mesh nodes back to the boundary of an original mesh (orthonet_project_on_boundary)
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool ProjectOnOriginalMeshBoundary(Mesh& mesh);

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
        bool AllocateSmootherNodeOperators(int topologyIndex);

        /// <summary>
        /// If it is a new topology, save it
        /// </summary>
        /// <param name="currentNode"></param>
        /// <param name="numSharedFaces"></param>
        /// <param name="numConnectedNodes"></param>
        /// <returns></returns>
        bool SaveSmootherNodeTopologyIfNeeded(int currentNode, 
                                              int numSharedFaces, 
                                              int numConnectedNodes);

        /// <summary>
        /// Computes how much the coordinates have to be incremented every inner iteration.
        /// Assembles the contributions of smoother and orthogonalizer
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool ComputeLinearSystemTerms(const Mesh& mesh);

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

        /// <summary>
        /// 
        /// </summary>
        /// <param name="wwx"></param>
        /// <param name="wwy"></param>
        /// <param name="currentNode"></param>
        /// <param name="n"></param>
        /// <param name="mesh"></param>
        /// <param name="dx0"></param>
        /// <param name="dy0"></param>
        /// <param name="increments"></param>
        /// <returns></returns>
        bool ComputeLocalIncrements(double wwx, 
                                    double wwy, 
                                    int currentNode, 
                                    int n, 
                                    const Mesh& mesh, 
                                    double& dx0, 
                                    double& dy0, 
                                    double* increments);

        /// <summary>
        /// Compute orthogonal coordinates
        /// </summary>
        /// <param name="nodeIndex"></param>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool UpdateNodeCoordinates(int nodeIndex, const Mesh& mesh);

        /// <summary>
        /// Allocate linear system arrays
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool AllocateLinearSystem(const Mesh& mesh);

        /// <summary>
        /// Deallocate linear system arrays
        /// </summary>
        /// <returns></returns>
        bool DeallocateLinearSystem();

        // Land boundaries
        LandBoundaries m_landBoundaries;

        // Polygons
        Polygons m_polygons;

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
        std::vector<int>                                   m_compressedEndNodeIndex;
        std::vector<int>                                   m_compressedStartNodeIndex;
        int m_nodeCacheSize = 0;
                  
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

        std::vector<double>                                m_aspectRatios;
        std::vector<std::vector<double>>                   m_wSmoother;
        std::vector<int>                                   m_numConnectedNodes;        // (nmk2)
        std::vector<std::vector<std::size_t>>              m_connectedNodes;           // (kk2)
        std::vector<int>                                   m_localCoordinatesIndexes;  // (iloc)
        std::vector<Point>                                 m_localCoordinates;         // (xloc,yloc) 

        std::vector<std::vector<double>>                   m_wOrthogonalizer;
        std::vector<std::vector<double>>                   m_rhsOrthogonalizer;
        std::size_t                                        m_maxNumNeighbours;
        std::vector<std::vector<int>>                      m_nodesNodes;               // node neighbours 

        // orthogonalization iterations
        std::vector<Point>                                 m_orthogonalCoordinates;
        std::vector<int>                                   m_nearestPoints;
        std::vector<Point>                                 m_originalNodes;

        // Linear system terms
        std::vector<double>                                m_compressedWeightX;
        std::vector<double>                                m_compressedWeightY;
        std::vector<double>                                m_compressedRhs;
        std::vector<int>                                   m_compressedNodesNodes;

        // Class variables
        int m_maximumNumConnectedNodes = 0;
        int m_maximumNumSharedFaces = 0;
        double m_mumax;
        double m_mu;

        // nodes with errors
        std::vector<double>                                m_nodeXErrors;
        std::vector<double>                                m_nodeYErrors;
        std::vector<int>                                   m_nodeErrorCode;

        // run-time algorithm parameters                                                      
        bool m_keepCircumcentersAndMassCenters = false;                          
        double m_orthogonalizationToSmoothingFactor = 0.975;                          // Factor(0. <= ATPF <= 1.) between grid smoothing and grid ortho resp.
        double m_orthogonalizationToSmoothingFactorBoundary = 1.0;                    // ATPF_B minimum ATPF on the boundary
        double m_smoothorarea = 1.0;                                                  // Factor between smoother(1.0) and area - homogenizer(0.0)
        int m_orthogonalizationOuterIterations = 2;
        int m_orthogonalizationBoundaryIterations = 25;
        int m_orthogonalizationInnerIterations = 25;
        int m_isTriangulationRequired;
        int m_isAccountingForLandBoundariesRequired;
        int m_projectToLandBoundaryOption;

        static constexpr int m_topologyInitialSize = 10;
        static constexpr double m_thetaTolerance = 1e-4;
    };
}