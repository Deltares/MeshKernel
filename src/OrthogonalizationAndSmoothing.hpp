#pragma once

#include <vector>
#include "LandBoundaries.hpp"
#include "Polygons.hpp"
#include "Smoother.hpp"
#include "Orthogonalizer.hpp"
#include "OrthogonalizationParametersNative.hpp"

namespace GridGeom
{
    struct Point;
    enum class Projections;
    class Mesh;

    /// <summary>
    /// Orthogonalizion (optimize the aspect ratios) and and mesh smoothing (optimize internal face angles or area).
    /// </summary>
    class OrthogonalizationAndSmoothing
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
        bool PrapareOuterIteration(Mesh& mesh);

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
        /// Project mesh nodes back to the boundary of an original mesh (orthonet_project_on_boundary)
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool ProjectOnOriginalMeshBoundary(Mesh& mesh);

        /// <summary>
        /// Computes how much the coordinates have to be incremented every inner iteration.
        /// Assembles the contributions of smoother and orthogonalizer
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool ComputeLinearSystemTerms(const Mesh& mesh);

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

        /// <summary>
        /// Compute nodes local coordinates, sice-effects only for sphericalAccurate projection (comp_local_coords)
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool ComputeCoordinates(const Mesh& mesh);

        // Land boundaries
        LandBoundaries m_landBoundaries;

        // Polygons
        Polygons m_polygons;

        // Smoother
        Smoother m_smoother;

        // Orthogonalizer
        Orthogonalizer m_orthogonalizer;
        
        // Local coordinates for sphericalAccurate projection
        std::vector<int>                                   m_localCoordinatesIndexes;  // (iloc)
        std::vector<Point>                                 m_localCoordinates;         // (xloc,yloc) 

        // orthogonalization iterations
        std::vector<Point>                                 m_orthogonalCoordinates;
        std::vector<Point>                                 m_originalNodes;

        // Linear system terms
        int m_nodeCacheSize = 0;
        std::vector<int>                                   m_compressedEndNodeIndex;
        std::vector<int>                                   m_compressedStartNodeIndex;
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
    };
}