#include <vector>
#include <algorithm>
#include <numeric>
#include <string>

#include "Operations.cpp"
#include "Smoother.hpp"
#include "Orthogonalizer.hpp"
#include "OrthogonalizationAndSmoothing.hpp"
#include "Entities.hpp"
#include "Mesh.hpp"

bool GridGeom::OrthogonalizationAndSmoothing::Set(Mesh& mesh,
    int& isTriangulationRequired,
    int& isAccountingForLandBoundariesRequired,
    int& projectToLandBoundaryOption,
    GridGeomApi::OrthogonalizationParametersNative& orthogonalizationParametersNative,
    const Polygons& polygon,
    std::vector<Point>& landBoundaries)
{
    
    m_polygons = polygon;
    m_smoother = Smoother();
    m_orthogonalizer = Orthogonalizer();

    // Sets the node mask
    mesh.MaskNodesInPolygons(m_polygons, true);
    // Flag nodes outside the polygon as corner points
    for (auto n = 0; n < mesh.GetNumNodes(); n++)
    {
        if (mesh.m_nodeMask[n] == 0) 
        {
            mesh.m_nodesTypes[n] = 3;
        }
    }

    // TODO: calculate volume weights for areal smoother
    m_mumax = (1.0 - m_smoothorarea) * 0.5;
    m_mu = std::min(1e-2, m_mumax);
    m_orthogonalCoordinates.resize(mesh.GetNumNodes() );


    // back-up original nodes, for projection on original mesh boundary
    m_originalNodes = mesh.m_nodes;
    m_orthogonalCoordinates = mesh.m_nodes;

    // algorithm settings
    m_orthogonalizationToSmoothingFactor = orthogonalizationParametersNative.OrthogonalizationToSmoothingFactor;
    m_orthogonalizationToSmoothingFactorBoundary = orthogonalizationParametersNative.OrthogonalizationToSmoothingFactorBoundary;
    m_smoothorarea = orthogonalizationParametersNative.Smoothorarea;
    m_orthogonalizationOuterIterations = orthogonalizationParametersNative.OuterIterations;
    m_orthogonalizationBoundaryIterations = orthogonalizationParametersNative.BoundaryIterations;
    m_orthogonalizationInnerIterations = orthogonalizationParametersNative.InnerIterations;

    m_landBoundaries.Set(landBoundaries);

    m_isTriangulationRequired = isTriangulationRequired;

    m_isAccountingForLandBoundariesRequired = isAccountingForLandBoundariesRequired;

    m_projectToLandBoundaryOption = projectToLandBoundaryOption;

    // project on land boundary
    if (m_projectToLandBoundaryOption >= 1)
    {
        // account for enclosing polygon
        m_landBoundaries.Administrate(mesh, m_polygons);
        m_landBoundaries.FindNearestMeshBoundary(mesh, m_polygons, m_projectToLandBoundaryOption);
    }

    // for spherical accurate computations we need to call orthonet_comp_ops 
    if(mesh.m_projection == Projections::sphericalAccurate)
    {
        if(m_orthogonalizationToSmoothingFactor<1.0)
        {
            bool successful = PrapareOuterIteration(mesh);
            if (!successful)
            {
                return false;
            }
        }

        m_localCoordinatesIndexes.resize(mesh.GetNumNodes() + 1);
        m_localCoordinatesIndexes[0] = 1;
        for (int n = 0; n < mesh.GetNumNodes(); ++n)
        {
            m_localCoordinatesIndexes[n + 1] = m_localCoordinatesIndexes[n] + std::max(mesh.m_nodesNumEdges[n] + 1, m_smoother.GetNumConnectedNodes(n));
        }

        m_localCoordinates.resize(m_localCoordinatesIndexes.back() - 1, { doubleMissingValue, doubleMissingValue });
    }

    return true;
}

bool GridGeom::OrthogonalizationAndSmoothing::Compute(Mesh& mesh)
{
    bool successful = true;

    for (auto outerIter = 0; outerIter < m_orthogonalizationOuterIterations; outerIter++)
    {
        if (successful)
        {
            successful = PrapareOuterIteration(mesh);
        }
        for (auto boundaryIter = 0; boundaryIter < m_orthogonalizationBoundaryIterations; boundaryIter++)
        {
            for (auto innerIter = 0; innerIter < m_orthogonalizationInnerIterations; innerIter++)
            {
                if (successful)
                {
                    successful = InnerIteration(mesh);
                }        
            } // inner iteration
        } // boundary iter

        //update mu
        if (successful)
        {
            successful = FinalizeOuterIteration(mesh);
        }    
    }// outer iter

    DeallocateLinearSystem();

    return true;
}


bool GridGeom::OrthogonalizationAndSmoothing::PrapareOuterIteration(Mesh& mesh) 
{

    bool successful = true;

    // compute weights and rhs of orthogonalizer
    if (successful)
    {
        successful = m_orthogonalizer.Compute(mesh);
    }

    // computes the smoother weights
    if (successful)
    {
        successful = m_smoother.Compute(mesh);
    }
     
    // allocate linear system for smoother and orthogonalizer
    if (successful)
    {
        successful = AllocateLinearSystem(mesh);
    }

    // compute linear system terms for smoother and orthogonalizer
    if (successful)
    {
        successful = ComputeLinearSystemTerms(mesh);
    }
    return successful;
}

bool GridGeom::OrthogonalizationAndSmoothing::AllocateLinearSystem(const Mesh& mesh) 
{
    bool successful = true;
    // reallocate caches
    if (successful && m_nodeCacheSize == 0)
    {
        m_compressedRhs.resize(mesh.GetNumNodes() * 2);
        std::fill(m_compressedRhs.begin(), m_compressedRhs.end(), 0.0);

        m_compressedEndNodeIndex.resize(mesh.GetNumNodes());
        std::fill(m_compressedEndNodeIndex.begin(), m_compressedEndNodeIndex.end(), 0.0);

        m_compressedStartNodeIndex.resize(mesh.GetNumNodes() );
        std::fill(m_compressedStartNodeIndex.begin(), m_compressedStartNodeIndex.end(), 0.0);

        for (int n = 0; n < mesh.GetNumNodes() ; n++)
        {
            m_compressedEndNodeIndex[n] = m_nodeCacheSize;
            m_nodeCacheSize += std::max(mesh.m_nodesNumEdges[n] + 1, m_smoother.GetNumConnectedNodes(n));
            m_compressedStartNodeIndex[n] = m_nodeCacheSize;
        }

        m_compressedNodesNodes.resize(m_nodeCacheSize);
        m_compressedWeightX.resize(m_nodeCacheSize);
        m_compressedWeightY.resize(m_nodeCacheSize);
    }
    return successful;
}


bool GridGeom::OrthogonalizationAndSmoothing::DeallocateLinearSystem() 
{
    m_compressedRhs.resize(0);
    m_compressedEndNodeIndex.resize(0);
    m_compressedStartNodeIndex.resize(0);
    m_compressedNodesNodes.resize(0);
    m_compressedWeightX.resize(0);
    m_compressedWeightY.resize(0);
    m_nodeCacheSize = 0;

    return true;
}

bool GridGeom::OrthogonalizationAndSmoothing::FinalizeOuterIteration(Mesh& mesh)
{
    m_mu = std::min(2.0 * m_mu, m_mumax);

    //compute new faces circumcenters
    if (!m_keepCircumcentersAndMassCenters)
    {
        mesh.ComputeFaceCircumcentersMassCentersAndAreas();
    }

    return true;
}

bool GridGeom::OrthogonalizationAndSmoothing::ComputeLinearSystemTerms(const Mesh& mesh)
{
    double max_aptf = std::max(m_orthogonalizationToSmoothingFactorBoundary, m_orthogonalizationToSmoothingFactor);
#pragma omp parallel for
	for (int n = 0; n < mesh.GetNumNodes() ; n++)
    {    
        if ((mesh.m_nodesTypes[n] != 1 && mesh.m_nodesTypes[n] != 2) || mesh.m_nodesNumEdges[n] < 2)
        {
            continue;
        }
        if (m_keepCircumcentersAndMassCenters != false && (mesh.m_nodesNumEdges[n] != 3 || mesh.m_nodesNumEdges[n] != 1))
        {
            continue;
        }

        const double atpfLoc  = mesh.m_nodesTypes[n] == 2 ? max_aptf : m_orthogonalizationToSmoothingFactor;
        const double atpf1Loc = 1.0 - atpfLoc;
        double mumat    = m_mu;
        int maxnn = m_compressedStartNodeIndex[n] - m_compressedEndNodeIndex[n];
        for (int nn = 1, cacheIndex = m_compressedEndNodeIndex[n]; nn < maxnn; nn++, cacheIndex++)
        {
            double wwx = 0.0;
            double wwy = 0.0;

            // Smoother
            if (atpf1Loc > 0.0 && mesh.m_nodesTypes[n] == 1)
            {
                wwx = atpf1Loc * m_smoother.GetWeight(n, nn);
                wwy = atpf1Loc * m_smoother.GetWeight(n, nn);
            }
            
            // Orthogonalizer
            if (nn < mesh.m_nodesNumEdges[n] + 1)
            {
                wwx += atpfLoc * m_orthogonalizer.GetWeight(n,nn - 1);
                wwy += atpfLoc * m_orthogonalizer.GetWeight(n,nn - 1);
                m_compressedNodesNodes[cacheIndex] = mesh.m_nodesNodes[n][nn - 1];
            }
            else
            {
                m_compressedNodesNodes[cacheIndex] = m_smoother.GetCoonectedNodeIndex(n,nn);
            }

            m_compressedWeightX[cacheIndex] = wwx;
            m_compressedWeightY[cacheIndex] = wwy;
        }
        int firstCacheIndex = n * 2;
        m_compressedRhs[firstCacheIndex] = atpfLoc * m_orthogonalizer.GetRightHandSide(n, 0);
        m_compressedRhs[firstCacheIndex+1] = atpfLoc * m_orthogonalizer.GetRightHandSide(n, 1);
	}

	return true;
}


bool GridGeom::OrthogonalizationAndSmoothing::InnerIteration(Mesh& mesh)
{
#pragma omp parallel for
	for (int n = 0; n < mesh.GetNumNodes() ; n++)
    {
	    UpdateNodeCoordinates(n, mesh);   
    }
	
    // update mesh node coordinates
    mesh.m_nodes = m_orthogonalCoordinates;

    // project on the original net boundary
    ProjectOnOriginalMeshBoundary(mesh);

    // compute local coordinates
    ComputeCoordinates(mesh);

    // project on land boundary
    if (m_projectToLandBoundaryOption >= 1)
    {
        m_landBoundaries.SnapMeshToLandBoundaries(mesh);
    }

    return true;
}

bool GridGeom::OrthogonalizationAndSmoothing::ProjectOnOriginalMeshBoundary(Mesh& mesh)
{
    Point firstPoint;
    Point secondPoint;
    Point thirdPoint;
    Point normalSecondPoint;
    Point normalThirdPoint;

    // in this case the nearest point is the point itself
    std::vector<int>   nearestPoints(mesh.GetNumNodes(), 0);
    std::iota(nearestPoints.begin(), nearestPoints.end(), 0);

    for (auto n = 0; n < mesh.GetNumNodes() ; n++)
    {
        int nearestPointIndex = nearestPoints[n];
        if (mesh.m_nodesTypes[n] == 2 && mesh.m_nodesNumEdges[n] > 0 && mesh.m_nodesNumEdges[nearestPointIndex] > 0)
        {
            firstPoint = mesh.m_nodes[n];
            int numEdges = mesh.m_nodesNumEdges[nearestPointIndex];
            int numNodes = 0;
            int leftNode;
            int rightNode;
            for (auto nn = 0; nn < numEdges; nn++)
            {
                auto edgeIndex = mesh.m_nodesEdges[nearestPointIndex][nn];
                if (mesh.m_edgesNumFaces[edgeIndex] == 1)
                {
                    numNodes++;
                    if (numNodes == 1)
                    {
                        leftNode = mesh.m_nodesNodes[n][nn];
                        if (leftNode == intMissingValue)
                        {
                            return false;
                        }
                        secondPoint = m_originalNodes[leftNode];
                    }
                    else if (numNodes == 2)
                    {
                        rightNode = mesh.m_nodesNodes[n][nn];
                        if (rightNode == intMissingValue)
                        {
                            return false;
                        }
                        thirdPoint = m_originalNodes[rightNode];
                    }
                }
            }

            //Project the moved boundary point back onto the closest original edge (either between 0 and 2 or 0 and 3)
            double rl2;
            double dis2 = DistanceFromLine(firstPoint, m_originalNodes[nearestPointIndex], secondPoint, normalSecondPoint, rl2, mesh.m_projection);

            double rl3;
            double dis3 = DistanceFromLine(firstPoint, m_originalNodes[nearestPointIndex], thirdPoint, normalThirdPoint, rl3, mesh.m_projection);

            if (dis2 < dis3)
            {
                mesh.m_nodes[n] = normalSecondPoint;
                if (rl2 > 0.5 && mesh.m_nodesTypes[n] != 3)
                {
                    nearestPoints[n] = leftNode;
                }
            }
            else
            {
                mesh.m_nodes[n] = normalThirdPoint;
                if (rl3 > 0.5 && mesh.m_nodesTypes[n] != 3)
                {
                    nearestPoints[n] = rightNode;
                }
            }
        }
    }
    return true;
}


bool GridGeom::OrthogonalizationAndSmoothing::ComputeCoordinates(const Mesh& mesh)
{
    if(mesh.m_projection == Projections::sphericalAccurate)
    {
        //TODO : missing implementation
        return true;
    }

    return true;
}

bool GridGeom::OrthogonalizationAndSmoothing::GetOrthogonality(const Mesh& mesh, double* orthogonality)
{
    for(int e=0; e < mesh.GetNumEdges() ; e++)
    {
        orthogonality[e] = doubleMissingValue;
        int firstVertex = mesh.m_edges[e].first;
        int secondVertex = mesh.m_edges[e].second;

        if (firstVertex!=0 && secondVertex !=0)
        {
            if (e < mesh.GetNumEdges() && mesh.m_edgesNumFaces[e]==2 )
            {
                orthogonality[e] = NormalizedInnerProductTwoSegments(mesh.m_nodes[firstVertex], mesh.m_nodes[secondVertex],
                    mesh.m_facesCircumcenters[mesh.m_edgesFaces[e][0]], mesh.m_facesCircumcenters[mesh.m_edgesFaces[e][1]], mesh.m_projection);
                if (orthogonality[e] != doubleMissingValue)
                {
                    orthogonality[e] = std::abs(orthogonality[e]);

                }
            }
        }
    }
    return true;
}

bool GridGeom::OrthogonalizationAndSmoothing::GetSmoothness(const Mesh& mesh, double* smoothness)
{
    for (int e = 0; e < mesh.GetNumEdges(); e++)
    {
        smoothness[e] = doubleMissingValue;
        int firstVertex = mesh.m_edges[e].first;
        int secondVertex = mesh.m_edges[e].second;

        if (firstVertex != 0 && secondVertex != 0)
        {
            if (e < mesh.GetNumEdges() && mesh.m_edgesNumFaces[e] == 2)
            {
                const auto leftFace = mesh.m_edgesFaces[e][0];
                const auto rightFace = mesh.m_edgesFaces[e][1];
                const auto leftFaceArea = mesh.m_faceArea[leftFace];
                const auto rightFaceArea = mesh.m_faceArea[rightFace];

                if (leftFaceArea< minimumCellArea || rightFaceArea< minimumCellArea)
                {
                    smoothness[e] = rightFaceArea / leftFaceArea;
                }
                if (smoothness[e] < 1.0) 
                {
                    smoothness[e] = 1.0 / smoothness[e];
                }
            }
        }
    }
    return true;
}

bool GridGeom::OrthogonalizationAndSmoothing::UpdateNodeCoordinates(int nodeIndex, const Mesh& mesh)
{
    int numConnectedNodes = m_compressedStartNodeIndex[nodeIndex] - m_compressedEndNodeIndex[nodeIndex];    
    double dx0 = 0.0;
    double dy0 = 0.0;
    double increments[2]{ 0.0, 0.0 };
    for (int nn = 1, cacheIndex = m_compressedEndNodeIndex[nodeIndex]; nn < numConnectedNodes; nn++, cacheIndex++)
    {
        ComputeLocalIncrements(m_compressedWeightX[cacheIndex], m_compressedWeightY[cacheIndex], m_compressedNodesNodes[cacheIndex], nodeIndex, mesh, dx0, dy0, increments);
    }

    if (increments[0] <= 1e-8 || increments[1] <= 1e-8)
    {
        return true;
    }

    int firstCacheIndex = nodeIndex * 2;
    dx0 = (dx0 + m_compressedRhs[firstCacheIndex]) / increments[0];
    dy0 = (dy0 + m_compressedRhs[firstCacheIndex + 1]) / increments[1];
    constexpr double relaxationFactor = 0.75;
    if (mesh.m_projection == Projections::cartesian || mesh.m_projection == Projections::spherical)
    {
        double x0 = mesh.m_nodes[nodeIndex].x + dx0;
        double y0 = mesh.m_nodes[nodeIndex].y + dy0;
        static constexpr double relaxationFactorCoordinates = 1.0 - relaxationFactor;

        m_orthogonalCoordinates[nodeIndex].x = relaxationFactor * x0 + relaxationFactorCoordinates * mesh.m_nodes[nodeIndex].x;
        m_orthogonalCoordinates[nodeIndex].y = relaxationFactor * y0 + relaxationFactorCoordinates * mesh.m_nodes[nodeIndex].y;
    }
    if (mesh.m_projection == Projections::sphericalAccurate)
    {
        Point localPoint{ relaxationFactor * dx0, relaxationFactor * dy0 };

        double exxp[3];
        double eyyp[3];
        double ezzp[3];
        ComputeThreeBaseComponents(mesh.m_nodes[nodeIndex], exxp, eyyp, ezzp);

        //get 3D-coordinates in rotated frame
        Cartesian3DPoint cartesianLocalPoint;
        SphericalToCartesian(localPoint, cartesianLocalPoint);

        //project to fixed frame
        Cartesian3DPoint transformedCartesianLocalPoint;
        transformedCartesianLocalPoint.x = exxp[0] * cartesianLocalPoint.x + eyyp[0] * cartesianLocalPoint.y + ezzp[0] * cartesianLocalPoint.z;
        transformedCartesianLocalPoint.y = exxp[1] * cartesianLocalPoint.x + eyyp[1] * cartesianLocalPoint.y + ezzp[1] * cartesianLocalPoint.z;
        transformedCartesianLocalPoint.z = exxp[2] * cartesianLocalPoint.x + eyyp[2] * cartesianLocalPoint.y + ezzp[2] * cartesianLocalPoint.z;

        //tranform to spherical coordinates
        CartesianToSpherical(transformedCartesianLocalPoint, mesh.m_nodes[nodeIndex].x, m_orthogonalCoordinates[nodeIndex]);
    }
    return true;
}

bool GridGeom::OrthogonalizationAndSmoothing::ComputeLocalIncrements(double wwx, double wwy, int currentNode, int n, const Mesh& mesh, double& dx0, double& dy0, double* increments)
{

    double wwxTransformed;
    double wwyTransformed;
    if (mesh.m_projection == Projections::cartesian)
    {
        wwxTransformed = wwx;
        wwyTransformed = wwy;
        dx0 = dx0 + wwxTransformed * (mesh.m_nodes[currentNode].x - mesh.m_nodes[n].x);
        dy0 = dy0 + wwyTransformed * (mesh.m_nodes[currentNode].y - mesh.m_nodes[n].y);
    }

    if (mesh.m_projection == Projections::spherical)
    {
        wwxTransformed = wwx * earth_radius * degrad_hp *
            std::cos(0.5 * (mesh.m_nodes[n].y + mesh.m_nodes[currentNode].y) * degrad_hp);
        wwyTransformed = wwy * earth_radius * degrad_hp;

        dx0 = dx0 + wwxTransformed * (mesh.m_nodes[currentNode].x - mesh.m_nodes[n].x);
        dy0 = dy0 + wwyTransformed * (mesh.m_nodes[currentNode].y - mesh.m_nodes[n].y);

    }
    if (mesh.m_projection == Projections::sphericalAccurate)
    {
        wwxTransformed = wwx * earth_radius * degrad_hp;
        wwyTransformed = wwy * earth_radius * degrad_hp;

        dx0 = dx0 + wwxTransformed * m_localCoordinates[m_localCoordinatesIndexes[n] + currentNode - 1].x;
        dy0 = dy0 + wwyTransformed * m_localCoordinates[m_localCoordinatesIndexes[n] + currentNode - 1].y;
    }

    increments[0] += wwxTransformed;
    increments[1] += wwyTransformed;

    return true;
}
