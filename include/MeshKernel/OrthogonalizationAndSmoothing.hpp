//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2020.
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
#include <memory>
#include <MeshKernel/OrthogonalizationParametersNative.hpp>

namespace meshkernel
{
    // Forward declare everything to reduce compile time dependency
    struct Point;
    class Mesh;
    class Smoother;
    class Orthogonalizer;
    class LandBoundaries;
    class Polygons;
    enum class Projections;

    /// <summary>
    /// Orthogonalizion (optimize the aspect ratios) and and mesh smoothing (optimize internal face angles or area).
    /// </summary>
    class OrthogonalizationAndSmoothing
    {

    public:
        /// Set the parameters
        /// @param[in] mesh The mesh to orthogonalize
        /// @param[in] smoother The mesh to smoother
        /// @param[in] orthogonalizer The mesh to orthogonalizer
        /// @param[in] isTriangulationRequired Not used
        /// @param[in] isAccountingForLandBoundariesRequired Not used
        /// @param[in] projectToLandBoundaryOption Snap to land boundaries (1) or not (0)
        /// @param[in] orthogonalizationParametersNative The orthogonalization parameters
        /// @param[in] polygon The polygon where orthogonalization should occur
        /// @param[in] landBoundaries The land boundaries
        OrthogonalizationAndSmoothing(std::shared_ptr<Mesh> mesh,
                                      std::shared_ptr<Smoother> smoother,
                                      std::shared_ptr<Orthogonalizer> orthogonalizer,
                                      std::shared_ptr<Polygons> polygon,
                                      std::shared_ptr<LandBoundaries> landBoundaries,
                                      int projectToLandBoundaryOption,
                                      const meshkernelapi::OrthogonalizationParametersNative& orthogonalizationParametersNative);

        void Initialize();

        /// @brief Executes the entire algorithm
        void Compute();

        /// @brief Prepares the outer iteration, calculates orthogonalizer and smoother coefficients and assable the linear system
        void PrepareOuterIteration();

        /// @brief Performs an inner iteration, update the mesh node positions
        void InnerIteration();

        /// @brief Finalize the outer iteration, computes new mu and face areas, masscenters, circumcenters
        void FinalizeOuterIteration();

    private:
        /// @brief Project mesh nodes back to the original mesh boundary (orthonet_project_on_boundary)
        void ProjectOnOriginalMeshBoundary();

        /// @brief Assembles the contributions of smoother and orthogonalizer
        void ComputeLinearSystemTerms();

        /// Computes how much the coordinates of a node need to be incremented at each inner iteration.
        /// @param[in] nodeIndex The node index
        /// @param[out] dx0 The computed x increment
        /// @param[out] dy0 The computed y increment
        /// @param[out] weightsSum The sum of the weights in x and y
        void ComputeLocalIncrements(int nodeIndex,
                                    double& dx0,
                                    double& dy0,
                                    std::array<double, 2>& weightsSum);

        /// @brief Update the nodal coordinates based on the increments
        /// @param[in] nodeIndex
        void UpdateNodeCoordinates(int nodeIndex);

        /// @brief Allocate linear system vectors
        void AllocateLinearSystem();

        /// @brief Compute nodes local coordinates, sice-effects only for sphericalAccurate projection (comp_local_coords)
        void ComputeCoordinates() const;

        std::shared_ptr<Mesh> m_mesh;                                                         // A pointer to mesh
        std::shared_ptr<Smoother> m_smoother;                                                 // A pointer to the smoother
        std::shared_ptr<Orthogonalizer> m_orthogonalizer;                                     // A pointer to the orthogonalizer
        std::shared_ptr<Polygons> m_polygons;                                                 // The polygon where to perform the orthogonalization
        std::shared_ptr<LandBoundaries> m_landBoundaries;                                     // The land boundaries
        int m_projectToLandBoundaryOption;                                                    // The project to land boundary option
        meshkernelapi::OrthogonalizationParametersNative m_orthogonalizationParametersNative; // The orthogonalization parameters

        std::vector<int> m_localCoordinatesIndexes; // Used in sphericalAccurate projection (iloc)
        std::vector<Point> m_localCoordinates;      // Used in sphericalAccurate projection (xloc,yloc)
        std::vector<Point> m_orthogonalCoordinates; // A copy of the mesh node, orthogonalized
        std::vector<Point> m_originalNodes;         // The original mesh

        // Linear system terms
        int m_nodeCacheSize = 0;
        std::vector<int> m_compressedEndNodeIndex;   // Start index in m_compressedWeightX
        std::vector<int> m_compressedStartNodeIndex; // End index in m_compressedWeightY
        std::vector<double> m_compressedWeightX;     // The computed weights X
        std::vector<double> m_compressedWeightY;     // The computed weights Y
        std::vector<double> m_compressedRhs;         // The right hand side
        std::vector<int> m_compressedNodesNodes;     // The indices of the neighbouring nodes
        std::vector<int> m_nodeErrorCode;            // nodes with errors

        // run-time parameters
        double m_mumax;
        double m_mu;
        bool m_keepCircumcentersAndMassCenters = false;
    };
} // namespace meshkernel
