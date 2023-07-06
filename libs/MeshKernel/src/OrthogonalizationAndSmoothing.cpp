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

#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/LandBoundaries.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/OrthogonalizationAndSmoothing.hpp>
#include <MeshKernel/Orthogonalizer.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/Smoother.hpp>

using meshkernel::Mesh2D;
using meshkernel::OrthogonalizationAndSmoothing;

OrthogonalizationAndSmoothing::OrthogonalizationAndSmoothing(std::shared_ptr<Mesh2D> mesh,
                                                             std::shared_ptr<Smoother> smoother,
                                                             std::shared_ptr<Orthogonalizer> orthogonalizer,
                                                             std::shared_ptr<Polygons> polygon,
                                                             std::shared_ptr<LandBoundaries> landBoundaries,
                                                             LandBoundaries::ProjectToLandBoundaryOption projectToLandBoundaryOption,
                                                             const OrthogonalizationParameters& orthogonalizationParameters)
    : m_mesh(mesh),
      m_smoother(smoother),
      m_orthogonalizer(orthogonalizer),
      m_polygons(polygon),
      m_landBoundaries(landBoundaries),
      m_projectToLandBoundaryOption(projectToLandBoundaryOption),
      m_orthogonalizationParameters(orthogonalizationParameters)
{
}

void OrthogonalizationAndSmoothing::Initialize()
{
    // Sets the node mask
    m_mesh->Administrate();
    const auto nodeMask = m_mesh->NodeMaskFromPolygon(*m_polygons, true);

    // Flag nodes outside the polygon as corner points
    for (size_t n = 0; n < nodeMask.size(); n++)
    {
        if (nodeMask[n] == 0)
        {
            m_mesh->m_nodesTypes[n] = 3;
        }
    }

    // TODO: calculate volume weights for areal smoother
    m_mumax = (1.0 - m_orthogonalizationParameters.areal_to_angle_smoothing_factor) * 0.5;
    m_mu = std::min(1e-2, m_mumax);

    // back-up original nodes, for projection on original mesh boundary
    m_originalNodes = m_mesh->m_nodes;
    m_orthogonalCoordinates = m_mesh->m_nodes;

    // account for enclosing polygon
    m_landBoundaries->FindNearestMeshBoundary(m_projectToLandBoundaryOption);

    // for spherical accurate computations we need to call PrapareOuterIteration (orthonet_comp_ops)
    if (m_mesh->m_projection == Projection::sphericalAccurate)
    {
        if (m_orthogonalizationParameters.orthogonalization_to_smoothing_factor < 1.0)
        {
            PrepareOuterIteration();
        }

        m_localCoordinatesIndices.resize(m_mesh->GetNumNodes() + 1);
        m_localCoordinatesIndices[0] = 1;
        for (size_t n = 0; n < m_mesh->GetNumNodes(); ++n)
        {
            m_localCoordinatesIndices[n + 1] = m_localCoordinatesIndices[n] + std::max(m_mesh->m_nodesNumEdges[n] + 1, m_smoother->GetNumConnectedNodes(n));
        }

        m_localCoordinates.resize(m_localCoordinatesIndices.back() - 1, {constants::missing::doubleValue, constants::missing::doubleValue});
    }
}

void OrthogonalizationAndSmoothing::Compute()
{
    for (auto outerIter = 0; outerIter < m_orthogonalizationParameters.outer_iterations; outerIter++)
    {
        PrepareOuterIteration();
        for (auto boundaryIter = 0; boundaryIter < m_orthogonalizationParameters.boundary_iterations; boundaryIter++)
        {
            for (auto innerIter = 0; innerIter < m_orthogonalizationParameters.inner_iterations; innerIter++)
            {
                Solve();

            } // inner iteration
        }     // boundary iter

        // update mu
        FinalizeOuterIteration();
    } // outer iter
}

void OrthogonalizationAndSmoothing::PrepareOuterIteration()
{
    // compute weights and rhs of orthogonalizer
    m_orthogonalizer->Compute();

    // computes the smoother weights
    m_smoother->Compute();

    // allocate linear system for smoother and orthogonalizer
    AllocateLinearSystem();

    // compute linear system terms for smoother and orthogonalizer
    ComputeLinearSystemTerms();
}

void OrthogonalizationAndSmoothing::AllocateLinearSystem()
{
    // reallocate caches
    if (m_nodeCacheSize == 0)
    {
        m_compressedRhs.resize(m_mesh->GetNumNodes() * 2);
        std::fill(m_compressedRhs.begin(), m_compressedRhs.end(), 0.0);

        m_compressedEndNodeIndex.resize(m_mesh->GetNumNodes());
        std::fill(m_compressedEndNodeIndex.begin(), m_compressedEndNodeIndex.end(), 0);

        m_compressedStartNodeIndex.resize(m_mesh->GetNumNodes());
        std::fill(m_compressedStartNodeIndex.begin(), m_compressedStartNodeIndex.end(), 0);

        for (size_t n = 0; n < m_mesh->GetNumNodes(); n++)
        {
            m_compressedEndNodeIndex[n] = m_nodeCacheSize;
            m_nodeCacheSize += std::max(m_mesh->m_nodesNumEdges[n] + 1, m_smoother->GetNumConnectedNodes(n));
            m_compressedStartNodeIndex[n] = m_nodeCacheSize;
        }

        m_compressedNodesNodes.resize(m_nodeCacheSize);
        m_compressedWeightX.resize(m_nodeCacheSize);
        m_compressedWeightY.resize(m_nodeCacheSize);
    }
}

void OrthogonalizationAndSmoothing::FinalizeOuterIteration()
{
    m_mu = std::min(2.0 * m_mu, m_mumax);
    m_mesh->ComputeCircumcentersMassCentersAndFaceAreas(true);
}

void OrthogonalizationAndSmoothing::ComputeLinearSystemTerms()
{
    const double max_aptf = std::max(m_orthogonalizationParameters.orthogonalization_to_smoothing_factor_at_boundary, m_orthogonalizationParameters.orthogonalization_to_smoothing_factor);
#pragma omp parallel for
    for (int n = 0; n < static_cast<int>(m_mesh->GetNumNodes()); n++)
    {
        if ((m_mesh->m_nodesTypes[n] != 1 && m_mesh->m_nodesTypes[n] != 2) || m_mesh->m_nodesNumEdges[n] < 2)
        {
            continue;
        }

        const double atpfLoc = m_mesh->m_nodesTypes[n] == 2 ? max_aptf : m_orthogonalizationParameters.orthogonalization_to_smoothing_factor;
        const double atpf1Loc = 1.0 - atpfLoc;
        const auto maxnn = m_compressedStartNodeIndex[n] - m_compressedEndNodeIndex[n];
        auto cacheIndex = m_compressedEndNodeIndex[n];
        for (int nn = 1; nn < static_cast<int>(maxnn); nn++)
        {
            double wwx = 0.0;
            double wwy = 0.0;

            // Smoother
            if (atpf1Loc > 0.0 && m_mesh->m_nodesTypes[n] == 1)
            {
                wwx = atpf1Loc * m_smoother->GetWeight(n, nn);
                wwy = atpf1Loc * m_smoother->GetWeight(n, nn);
            }

            // Orthogonalizer
            if (nn < static_cast<int>(m_mesh->m_nodesNumEdges[n]) + 1)
            {
                wwx += atpfLoc * m_orthogonalizer->GetWeight(n, nn - 1);
                wwy += atpfLoc * m_orthogonalizer->GetWeight(n, nn - 1);
                m_compressedNodesNodes[cacheIndex] = m_mesh->m_nodesNodes[n][nn - 1];
            }
            else
            {
                m_compressedNodesNodes[cacheIndex] = m_smoother->GetConnectedNodeIndex(n, nn);
            }

            m_compressedWeightX[cacheIndex] = wwx;
            m_compressedWeightY[cacheIndex] = wwy;
            cacheIndex++;
        }
        const Index firstCacheIndex = n * 2;
        m_compressedRhs[firstCacheIndex] = atpfLoc * m_orthogonalizer->GetRightHandSide(n, 0);
        m_compressedRhs[firstCacheIndex + 1] = atpfLoc * m_orthogonalizer->GetRightHandSide(n, 1);
    }
}

void OrthogonalizationAndSmoothing::Solve()
{
#pragma omp parallel for
    for (int n = 0; n < static_cast<int>(m_mesh->GetNumNodes()); n++)
    {
        UpdateNodeCoordinates(n);
    }

    // update mesh node coordinates
    m_mesh->m_nodes = m_orthogonalCoordinates;

    // project on the original net boundary
    SnapMeshToOriginalMeshBoundary();

    // compute local coordinates
    // TODO: Not implemented yet ComputeCoordinates();

    // project on land boundary
    m_landBoundaries->SnapMeshToLandBoundaries();
}

void OrthogonalizationAndSmoothing::SnapMeshToOriginalMeshBoundary()
{
    // in this case the nearest point is the point itself
    std::vector<Index> nearestPoints(m_mesh->GetNumNodes(), 0);
    std::iota(nearestPoints.begin(), nearestPoints.end(), 0);

    for (Index n = 0; n < m_mesh->GetNumNodes(); n++)
    {
        const auto nearestPointIndex = nearestPoints[n];
        if (m_mesh->m_nodesTypes[n] == 2 && m_mesh->m_nodesNumEdges[n] > 0 && m_mesh->m_nodesNumEdges[nearestPointIndex] > 0)
        {
            Point firstPoint = m_mesh->m_nodes[n];
            if (!firstPoint.IsValid())
            {
                continue;
            }

            const auto numEdges = m_mesh->m_nodesNumEdges[nearestPointIndex];
            Index numNodes = 0;
            Index leftNode = constants::missing::sizetValue;
            Index rightNode = constants::missing::sizetValue;
            Point secondPoint{constants::missing::doubleValue, constants::missing::doubleValue};
            Point thirdPoint{constants::missing::doubleValue, constants::missing::doubleValue};
            for (Index nn = 0; nn < numEdges; nn++)
            {
                const auto edgeIndex = m_mesh->m_nodesEdges[nearestPointIndex][nn];
                if (m_mesh->IsEdgeOnBoundary(edgeIndex))
                {
                    numNodes++;
                    if (numNodes == 1)
                    {
                        leftNode = m_mesh->m_nodesNodes[n][nn];
                        if (leftNode == constants::missing::sizetValue)
                        {
                            throw AlgorithmError("OrthogonalizationAndSmoothing::SnapMeshToOriginalMeshBoundary: The left node is invalid.");
                        }
                        secondPoint = m_originalNodes[leftNode];
                    }
                    else if (numNodes == 2)
                    {
                        rightNode = m_mesh->m_nodesNodes[n][nn];
                        if (rightNode == constants::missing::sizetValue)
                        {
                            throw AlgorithmError("OrthogonalizationAndSmoothing::SnapMeshToOriginalMeshBoundary: The right node is invalid.");
                        }
                        thirdPoint = m_originalNodes[rightNode];
                    }
                }
            }

            if (!secondPoint.IsValid() || !thirdPoint.IsValid())
            {
                continue;
            }

            // Project the moved boundary point back onto the closest original edge (either between 0 and 2 or 0 and 3)

            const auto [distanceSecondPoint, normalSecondPoint, ratioSecondPoint] = DistanceFromLine(firstPoint, m_originalNodes[nearestPointIndex], secondPoint, m_mesh->m_projection);

            const auto [distanceThirdPoint, normalThirdPoint, ratioThirdPoint] = DistanceFromLine(firstPoint, m_originalNodes[nearestPointIndex], thirdPoint, m_mesh->m_projection);

            if (distanceSecondPoint < distanceThirdPoint)
            {
                m_mesh->m_nodes[n] = normalSecondPoint;
                if (ratioSecondPoint > 0.5 && m_mesh->m_nodesTypes[n] != 3)
                {
                    nearestPoints[n] = leftNode;
                }
            }
            else
            {
                m_mesh->m_nodes[n] = normalThirdPoint;
                if (ratioThirdPoint > 0.5 && m_mesh->m_nodesTypes[n] != 3)
                {
                    nearestPoints[n] = rightNode;
                }
            }
        }
    }
}

void OrthogonalizationAndSmoothing::ComputeCoordinates() const
{
    // TODO :  implementation for m_mesh->m_projection == Projection::sphericalAccurate

    if (m_mesh->m_projection == Projection::sphericalAccurate)
    {
    }
    throw AlgorithmError("OrthogonalizationAndSmoothing::ComputeCoordinates: This functionality is not implemented yet.");
}

void OrthogonalizationAndSmoothing::UpdateNodeCoordinates(Index nodeIndex)
{

    double dx0 = 0.0;
    double dy0 = 0.0;
    std::array<double, 2> increments{0.0, 0.0};
    ComputeLocalIncrements(nodeIndex, dx0, dy0, increments);

    if (increments[0] <= 1e-8 || increments[1] <= 1e-8)
    {
        return;
    }

    const auto firstCacheIndex = nodeIndex * 2;
    dx0 = (dx0 + m_compressedRhs[firstCacheIndex]) / increments[0];
    dy0 = (dy0 + m_compressedRhs[firstCacheIndex + 1]) / increments[1];
    constexpr double relaxationFactor = 0.75;
    if (m_mesh->m_projection == Projection::cartesian || m_mesh->m_projection == Projection::spherical)
    {
        const double x0 = m_mesh->m_nodes[nodeIndex].x + dx0;
        const double y0 = m_mesh->m_nodes[nodeIndex].y + dy0;
        static constexpr double relaxationFactorCoordinates = 1.0 - relaxationFactor;

        m_orthogonalCoordinates[nodeIndex].x = relaxationFactor * x0 + relaxationFactorCoordinates * m_mesh->m_nodes[nodeIndex].x;
        m_orthogonalCoordinates[nodeIndex].y = relaxationFactor * y0 + relaxationFactorCoordinates * m_mesh->m_nodes[nodeIndex].y;
    }
    if (m_mesh->m_projection == Projection::sphericalAccurate)
    {
        const Point localPoint{relaxationFactor * dx0, relaxationFactor * dy0};

        std::array<double, 3> exxp{0.0, 0.0, 0.0};
        std::array<double, 3> eyyp{0.0, 0.0, 0.0};
        std::array<double, 3> ezzp{0.0, 0.0, 0.0};
        ComputeThreeBaseComponents(m_mesh->m_nodes[nodeIndex], exxp, eyyp, ezzp);

        // get 3D-coordinates in rotated frame
        const Cartesian3DPoint cartesianLocalPoint{SphericalToCartesian3D(localPoint)};

        // project to fixed frame
        Cartesian3DPoint transformedCartesianLocalPoint;
        transformedCartesianLocalPoint.x = exxp[0] * cartesianLocalPoint.x + eyyp[0] * cartesianLocalPoint.y + ezzp[0] * cartesianLocalPoint.z;
        transformedCartesianLocalPoint.y = exxp[1] * cartesianLocalPoint.x + eyyp[1] * cartesianLocalPoint.y + ezzp[1] * cartesianLocalPoint.z;
        transformedCartesianLocalPoint.z = exxp[2] * cartesianLocalPoint.x + eyyp[2] * cartesianLocalPoint.y + ezzp[2] * cartesianLocalPoint.z;

        // transform to spherical coordinates
        m_orthogonalCoordinates[nodeIndex] = Cartesian3DToSpherical(transformedCartesianLocalPoint, m_mesh->m_nodes[nodeIndex].x);
    }
}

void OrthogonalizationAndSmoothing::ComputeLocalIncrements(Index nodeIndex, double& dx0, double& dy0, std::array<double, 2>& weightsSum)
{
    const auto numConnectedNodes = m_compressedStartNodeIndex[nodeIndex] - m_compressedEndNodeIndex[nodeIndex];
    auto cacheIndex = m_compressedEndNodeIndex[nodeIndex];
    for (Index nn = 1; nn < numConnectedNodes; nn++)
    {
        const auto wwx = m_compressedWeightX[cacheIndex];
        const auto wwy = m_compressedWeightY[cacheIndex];
        const auto currentNode = m_compressedNodesNodes[cacheIndex];
        cacheIndex++;

        if (m_mesh->m_projection == Projection::cartesian)
        {
            const double wwxTransformed = wwx;
            const double wwyTransformed = wwy;
            dx0 = dx0 + wwxTransformed * (m_mesh->m_nodes[currentNode].x - m_mesh->m_nodes[nodeIndex].x);
            dy0 = dy0 + wwyTransformed * (m_mesh->m_nodes[currentNode].y - m_mesh->m_nodes[nodeIndex].y);
            weightsSum[0] += wwxTransformed;
            weightsSum[1] += wwyTransformed;
        }

        if (m_mesh->m_projection == Projection::spherical)
        {
            const double wwxTransformed = wwx * constants::geometric::earth_radius * constants::conversion::degToRad *
                                          std::cos(0.5 * (m_mesh->m_nodes[nodeIndex].y + m_mesh->m_nodes[currentNode].y) * constants::conversion::degToRad);
            const double wwyTransformed = wwy * constants::geometric::earth_radius * constants::conversion::degToRad;

            dx0 = dx0 + wwxTransformed * (m_mesh->m_nodes[currentNode].x - m_mesh->m_nodes[nodeIndex].x);
            dy0 = dy0 + wwyTransformed * (m_mesh->m_nodes[currentNode].y - m_mesh->m_nodes[nodeIndex].y);
            weightsSum[0] += wwxTransformed;
            weightsSum[1] += wwyTransformed;
        }
        if (m_mesh->m_projection == Projection::sphericalAccurate)
        {
            const double wwxTransformed = wwx * constants::geometric::earth_radius * constants::conversion::degToRad;
            const double wwyTransformed = wwy * constants::geometric::earth_radius * constants::conversion::degToRad;

            dx0 = dx0 + wwxTransformed * m_localCoordinates[m_localCoordinatesIndices[nodeIndex] + currentNode - 1].x;
            dy0 = dy0 + wwyTransformed * m_localCoordinates[m_localCoordinatesIndices[nodeIndex] + currentNode - 1].y;
            weightsSum[0] += wwxTransformed;
            weightsSum[1] += wwyTransformed;
        }
    }
}
