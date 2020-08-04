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

    class Orthogonalization
    {

    public:
        
        bool Set(Mesh& mesh,
            int& isTriangulationRequired,
            int& isAccountingForLandBoundariesRequired,
            int& projectToLandBoundaryOption,
            GridGeomApi::OrthogonalizationParametersNative& orthogonalizationParametersNative,
            const Polygons& polygon,
            std::vector<Point>& landBoundaries);

        bool Iterate(Mesh& mesh);

        bool PrapareOuterIteration(const Mesh& mesh);

        bool InnerIteration(Mesh& mesh);

        bool FinalizeOuterIteration(Mesh& mesh);

        bool GetOrthogonality(const Mesh& mesh, double* orthogonality);

        bool GetSmoothness(const Mesh& mesh, double* smoothness);

    private:

        /// <summary>
        /// Project mesh nodes back to the boundary of an original mesh (orthonet_project_on_boundary)
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool ProjectOnOriginalMeshBoundary(Mesh& mesh);

        /// snapping nodes to land boundaries

        /// <summary>
        /// Project nodes on land boundaries 
        /// </summary>
        /// <param name="mesh"></param>
        /// <param name="landBoundaries"></param>
        /// <returns></returns>
        bool ProjectOnLandBoundary(Mesh& mesh, const LandBoundaries& landBoundaries);

        /// <summary>
        /// Inverse-mapping elliptic smoother (orthonet_compweights_smooth)
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool ComputeWeightsSmoother(const Mesh& mesh);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool ComputeSmootherOperators(const Mesh& mesh);

        /// comp_local_coords

        /// <summary>
        /// (comp_local_coords)
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool ComputeLocalCoordinates(const Mesh& mesh);

        /// <summary>
        /// Compute coefficient matrix G of gradient at edge, compute coefficientmatrix Div of gradient in node, compute coefficientmatrix Az of face center (orthonet_comp_operators)
        /// </summary>
        /// <param name="mesh"></param>
        /// <param name="currentNode"></param>
        /// <param name="numConnectedNodes"></param>
        /// <param name="connectedNodes"></param>
        /// <param name="numSharedFaces"></param>
        /// <param name="sharedFaces"></param>
        /// <param name="xi"></param>
        /// <param name="eta"></param>
        /// <param name="faceNodeMapping"></param>
        /// <returns></returns>
        bool ComputeOperatorsNode(const Mesh& mesh, 
                                  int currentNode, 
                                  const std::size_t& numConnectedNodes, 
                                  const std::vector<std::size_t>& connectedNodes, 
                                  const std::size_t& numSharedFaces, 
                                  const std::vector<int>& sharedFaces,
                                  const std::vector<double>& xi, 
                                  const std::vector<double>& eta, 
                                  const std::vector<std::vector<std::size_t>>& faceNodeMapping);

        /// <summary>
        /// Assign xi and eta to all nodes in the stencil (orthonet_assign_xieta)
        /// </summary>
        /// <param name="mesh"></param>
        /// <param name="currentNode"></param>
        /// <param name="sharedFaces"></param>
        /// <param name="numSharedFaces"></param>
        /// <param name="connectedNodes"></param>
        /// <param name="numConnectedNodes"></param>
        /// <param name="faceNodeMapping"></param>
        /// <param name="xi"></param>
        /// <param name="eta"></param>
        /// <returns></returns>
        bool ComputeXiEta(const Mesh& mesh, 
                          int currentNode, 
                          const std::vector<int>& sharedFaces, 
                          const int& numSharedFaces, 
                          const std::vector<std::size_t>& connectedNodes,
                          const std::size_t& numConnectedNodes, 
                          const std::vector<std::vector<std::size_t>>& faceNodeMapping, 
                          std::vector<double>& xi, 
                          std::vector<double>& eta);

        /// <summary>
        /// Computes the shared faces and the connected nodes of a stencil node and the faceNodeMapping in the connectedNodes array for each shared face
        /// </summary>
        /// <param name="mesh"></param>
        /// <param name="currentNode"></param>
        /// <param name="sharedFaces"></param>
        /// <param name="numSharedFaces"></param>
        /// <param name="connectedNodes"></param>
        /// <param name="numConnectedNodes"></param>
        /// <param name="faceNodeMapping"></param>
        /// <returns></returns>
        bool OrthogonalizationAdministration(const Mesh& mesh, 
                                             const int currentNode, 
                                             std::vector<int>& sharedFaces, 
                                             int& numSharedFaces, 
                                             std::vector<std::size_t>& connectedNodes, 
                                             int& numConnectedNodes, 
                                             std::vector<std::vector<std::size_t>>& faceNodeMapping);

        /// <summary>
        /// Compute optimal angle
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
        /// Computes the aspect ratio of each edge (orthonet_compute_aspect)
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool AspectRatio(const Mesh& mesh);

        /// <summary>
        /// compute weights ww and right - hand side rhs in orthogonalizer (orthonet_compweights)
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool ComputeWeightsOrthogonalizer(const Mesh& mesh);

        /// <summary>
        /// Initialize available mesh topologies
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool InitializeTopologies(const Mesh& mesh);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="topologyIndex"></param>
        /// <returns></returns>
        bool AllocateNodeOperators(int topologyIndex);

        /// save only the unique topologies

        /// <summary>
        /// Save unique topologies
        /// </summary>
        /// <param name="currentNode"></param>
        /// <param name="sharedFaces"></param>
        /// <param name="numSharedFaces"></param>
        /// <param name="connectedNodes"></param>
        /// <param name="numConnectedNodes"></param>
        /// <param name="faceNodeMapping"></param>
        /// <param name="xi"></param>
        /// <param name="eta"></param>
        /// <returns></returns>
        bool SaveTopology(int currentNode, 
                          const std::vector<int>& sharedFaces, 
                          int numSharedFaces, 
                          const std::vector<std::size_t>& connectedNodes, 
                          int numConnectedNodes,
                          const std::vector<std::vector<std::size_t>>& faceNodeMapping, 
                          const std::vector<double>& xi, 
                          const std::vector<double>& eta);

        bool ComputeIncrements(const Mesh& mesh);

        bool ComputeJacobian(int currentNode, const Mesh& mesh, std::vector<double>& J) const;

        /// <summary>
        /// Compute the matric norm
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="matCoefficents"></param>
        /// <returns></returns>
        double MatrixNorm(const std::vector<double>& x,
                          const std::vector<double>& y,
                          const std::vector<double>& matCoefficents);

        /// <summary>
        /// Compute local increments
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
        bool ComputeOrthogonalCoordinates(int nodeIndex, const Mesh& mesh);

        /// <summary>
        /// Allocate caches
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool AllocateCaches(const Mesh& mesh);

        /// <summary>
        /// Deallocate caches
        /// </summary>
        /// <returns></returns>
        bool DeallocateCaches();

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
        std::vector<int> m_numConnectedNodes;                                    // nmk2, determined from local node administration
        std::vector<std::vector<std::size_t>> m_connectedNodes;                  // kk2, determined from local node administration
        std::vector<int> m_localCoordinatesIndexes;                              // iloc
        std::vector<Point> m_localCoordinates;                                   // xloc,yloc 
                                                                                 
        // run-time options                                                      
        bool m_keepCircumcentersAndMassCenters = false;                          
        double m_orthogonalizationToSmoothingFactor = 0.975;                     // Factor(0. <= ATPF <= 1.) between grid smoothing and grid ortho resp.
        double m_orthogonalizationToSmoothingFactorBoundary = 1.0;               // ATPF_B minimum ATPF on the boundary
        double m_smoothorarea = 1.0;                                             // Factor between smoother(1.0) and area - homogenizer(0.0)
        int m_orthogonalizationOuterIterations = 2;
        int m_orthogonalizationBoundaryIterations = 25;
        int m_orthogonalizationInnerIterations = 25;

        static constexpr int m_topologyInitialSize = 10;
        static constexpr double m_thetaTolerance = 1e-4;

        std::vector<std::vector<double>>  m_weights;
        std::vector<std::vector<double>>  m_rightHandSide;

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

        // the land boundaries
        LandBoundaries m_landBoundaries;

        // polygons
        Polygons m_polygons;

        int m_isTriangulationRequired;

        int m_isAccountingForLandBoundariesRequired;
        
        int m_projectToLandBoundaryOption;
    };
}