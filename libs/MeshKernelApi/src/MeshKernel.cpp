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

#include "MeshKernel/CurvilinearGrid/CurvilinearGridDeleteExterior.hpp"
#include "MeshKernel/Mesh2DIntersections.hpp"

#include "MeshKernel/CurvilinearGrid/CurvilinearGridDeleteInterior.hpp"
#include <MeshKernel/AveragingInterpolation.hpp>
#include <MeshKernel/BilinearInterpolationOnGriddedSamples.hpp>
#include <MeshKernel/ConnectMeshes.hpp>
#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Contacts.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridAlgorithm.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridDeRefinement.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridFromPolygon.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridFromSplines.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridFromSplinesTransfinite.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLineAttractionRepulsion.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLineMirror.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridOrthogonalization.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridRectangular.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridRefinement.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridSmoothing.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/FlipEdges.hpp>
#include <MeshKernel/LandBoundaries.hpp>
#include <MeshKernel/LandBoundary.hpp>
#include <MeshKernel/Mesh.hpp>
#include <MeshKernel/Mesh1D.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/MeshRefinement.hpp>
#include <MeshKernel/MeshTransformation.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/OrthogonalizationAndSmoothing.hpp>
#include <MeshKernel/Orthogonalizer.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/RemoveDisconnectedRegions.hpp>
#include <MeshKernel/Smoother.hpp>
#include <MeshKernel/SplineAlgorithms.hpp>
#include <MeshKernel/Splines.hpp>
#include <MeshKernel/TriangulationInterpolation.hpp>
#include <MeshKernel/Utilities/LinearAlgebra.hpp>

#include <MeshKernelApi/MeshKernel.hpp>
#include <MeshKernelApi/State.hpp>
#include <MeshKernelApi/Utils.hpp>

#include <Version/Version.hpp>

#include <cstring>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace meshkernelapi
{
    // The state held by MeshKernel
    static std::unordered_map<int, MeshKernelState> meshKernelState;
    static int meshKernelStateCounter = 0;

    // Error state
    static size_t constexpr bufferSize = 512;
    static size_t constexpr maxCharsToCopy = bufferSize - 1; // make sure destination string is null-terminated when strncpy is used
    static char exceptionMessage[bufferSize] = "";
    static meshkernel::ExitCode lastExitCode = meshkernel::ExitCode::Success;
    static meshkernel::UInt invalidMeshIndex{0};
    static meshkernel::Mesh::Location invalidMeshLocation{meshkernel::Mesh::Location::Unknown};

    static meshkernel::ExitCode HandleException(std::exception_ptr exception_ptr = std::current_exception())
    {
        try
        {
            std::rethrow_exception(exception_ptr);
        }
        catch (meshkernel::MeshGeometryError const& e)
        {
            std::strncpy(exceptionMessage, e.what(), maxCharsToCopy);
            invalidMeshIndex = e.MeshIndex();
            invalidMeshLocation = e.MeshLocation();
            return e.Code();
        }
        catch (meshkernel::MeshKernelError const& e)
        {
            std::strncpy(exceptionMessage, e.what(), maxCharsToCopy);
            return e.Code();
        }
        catch (std::exception const& e)
        {
            std::strncpy(exceptionMessage, e.what(), maxCharsToCopy);
            return meshkernel::ExitCode::StdLibExceptionCode;
        }
        catch (...)
        {
            std::strncpy(exceptionMessage, "Unknown exception", maxCharsToCopy);
            return meshkernel::ExitCode::UnknownExceptionCode;
        }
    }

    MKERNEL_API int mkernel_allocate_state(int projectionType, int& meshKernelId)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        meshKernelId = meshKernelStateCounter++;
        auto const projection = static_cast<meshkernel::Projection>(projectionType);
        meshKernelState.insert({meshKernelId, MeshKernelState(projection)});
        return lastExitCode;
    }

    MKERNEL_API int mkernel_deallocate_state(int meshKernelId)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }
            meshKernelState.erase(meshKernelId);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_delete(int meshKernelId, const GeometryList& polygon, int deletionOption, int invertDeletion)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }
            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() <= 0)
            {
                throw meshkernel::ConstraintError("The 2d mesh contains no nodes.");
            }

            const auto polygonPoints = ConvertGeometryListToPointVector(polygon);

            const bool invertDeletionBool = invertDeletion == 1 ? true : false;
            const meshkernel::Polygons meshKernelPolygon(polygonPoints, meshKernelState[meshKernelId].m_mesh2d->m_projection);
            meshKernelState[meshKernelId].m_mesh2d->DeleteMesh(meshKernelPolygon, deletionOption, invertDeletionBool);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_set(int meshKernelId, const Mesh2D& mesh2d)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            // convert raw arrays to containers
            const auto edges2d = meshkernel::ConvertToEdgeNodesVector(mesh2d.num_edges,
                                                                      mesh2d.edge_nodes);

            const auto nodes2d = meshkernel::ConvertToNodesVector(mesh2d.num_nodes,
                                                                  mesh2d.node_x,
                                                                  mesh2d.node_y);

            if (mesh2d.num_faces > 0 && mesh2d.face_nodes != nullptr && mesh2d.nodes_per_face != nullptr)
            {

                const auto face_nodes = meshkernel::ConvertToFaceNodesVector(mesh2d.num_faces, mesh2d.face_nodes, mesh2d.nodes_per_face);

                std::vector<meshkernel::UInt> num_face_nodes;
                num_face_nodes.reserve(mesh2d.num_faces);
                for (auto n = 0; n < mesh2d.num_faces; n++)
                {
                    num_face_nodes.emplace_back(static_cast<meshkernel::UInt>(mesh2d.nodes_per_face[n]));
                }

                // Do not change the pointer, just the object it is pointing to
                *meshKernelState[meshKernelId].m_mesh2d += meshkernel::Mesh2D(edges2d,
                                                                              nodes2d,
                                                                              face_nodes,
                                                                              num_face_nodes,
                                                                              meshKernelState[meshKernelId].m_projection);
            }
            else
            {
                // Do not change the pointer, just the object it is pointing to
                // Compute the faces
                *meshKernelState[meshKernelId].m_mesh2d += meshkernel::Mesh2D(edges2d, nodes2d, meshKernelState[meshKernelId].m_projection);
            }
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh1d_set(int meshKernelId,
                                       const Mesh1D& mesh1d)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }
            // convert raw arrays to containers
            const auto edges1d = meshkernel::ConvertToEdgeNodesVector(mesh1d.num_edges,
                                                                      mesh1d.edge_nodes);
            const auto nodes1d = meshkernel::ConvertToNodesVector(mesh1d.num_nodes,
                                                                  mesh1d.node_x,
                                                                  mesh1d.node_y);
            // Do not change the pointer, just the object it is pointing to
            *meshKernelState[meshKernelId].m_mesh1d = meshkernel::Mesh1D(edges1d, nodes1d, meshKernelState[meshKernelId].m_projection);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_network1d_set(int meshKernelId, const GeometryList& polylines)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            auto const localPolylines = ConvertGeometryListToVectorOfPointVectors(polylines);

            // Do not change the pointer, just the object it is pointing to
            *meshKernelState[meshKernelId].m_network1d = meshkernel::Network1D(localPolylines, meshKernelState[meshKernelId].m_projection);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_network1d_compute_fixed_chainages(int meshKernelId, double* fixedChainages, int sizeFixedChainages, double minFaceSize, double fixedChainagesOffset)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            std::vector<double> localFixedChainages(sizeFixedChainages);
            for (auto i = 0; i < sizeFixedChainages; ++i)
            {
                localFixedChainages[i] = fixedChainages[i];
            }
            const auto fixedChainagesByPolyline = ConvertVectorToVectorOfVectors(localFixedChainages, mkernel_get_separator());

            // Do not change the pointer, just the object it is pointing to
            meshKernelState[meshKernelId].m_network1d->ComputeFixedChainages(fixedChainagesByPolyline, minFaceSize, fixedChainagesOffset);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_network1d_compute_offsetted_chainages(int meshKernelId, double offset)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            // Do not change the pointer, just the object it is pointing to
            meshKernelState[meshKernelId].m_network1d->ComputeOffsettedChainages(offset);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_network1d_to_mesh1d(int meshKernelId, double minFaceSize)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            // Do not change the pointer, just the object it is pointing to (add to the existing mesh1d stored in the instance)
            *meshKernelState[meshKernelId].m_mesh1d += meshkernel::Mesh1D(*meshKernelState[meshKernelId].m_network1d, minFaceSize);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_dimensions(int meshKernelId, Mesh2D& mesh2d)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }
            meshKernelState[meshKernelId].m_mesh2d->Administrate();
            SetMesh2dApiDimensions(meshKernelState[meshKernelId].m_mesh2d, mesh2d);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_data(int meshKernelId, Mesh2D& mesh2d)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            SetMesh2dApiData(meshKernelState[meshKernelId].m_mesh2d, mesh2d);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh1d_get_dimensions(int meshKernelId, Mesh1D& mesh1d)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            mesh1d.num_nodes = static_cast<int>(meshKernelState[meshKernelId].m_mesh1d->GetNumNodes());
            mesh1d.num_edges = static_cast<int>(meshKernelState[meshKernelId].m_mesh1d->GetNumEdges());
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh1d_get_data(int meshKernelId,
                                            Mesh1D& mesh1d)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            SetMesh1dApiData(meshKernelState[meshKernelId].m_mesh1d, mesh1d);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_get_dimensions(int meshKernelId, CurvilinearGrid& curvilinearGrid)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }
            curvilinearGrid.num_n = static_cast<int>(meshKernelState[meshKernelId].m_curvilinearGrid->m_numN);
            curvilinearGrid.num_m = static_cast<int>(meshKernelState[meshKernelId].m_curvilinearGrid->m_numM);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_get_data(int meshKernelId,
                                                 CurvilinearGrid& curvilinearGrid)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }
            meshKernelState[meshKernelId].m_curvilinearGrid->SetFlatCopies();
            SetCurvilinearGridApiData(meshKernelState[meshKernelId].m_curvilinearGrid, curvilinearGrid);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_contacts_get_dimensions(int meshKernelId, Contacts& contacts)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }
            contacts.num_contacts = static_cast<int>(meshKernelState[meshKernelId].m_contacts->Mesh2dIndices().size());
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_contacts_get_data(int meshKernelId, Contacts& contacts)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            auto const& mesh1dIndices = meshKernelState[meshKernelId].m_contacts->Mesh1dIndices();
            auto const& mesh2dIndices = meshKernelState[meshKernelId].m_contacts->Mesh2dIndices();
            for (auto i = 0; i < contacts.num_contacts; ++i)
            {
                contacts.mesh1d_indices[i] = static_cast<int>(mesh1dIndices[i]);
                contacts.mesh2d_indices[i] = static_cast<int>(mesh2dIndices[i]);
            }
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_count_hanging_edges(int meshKernelId, int& numHangingEdges)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }
            meshKernelState[meshKernelId].m_mesh2d->Administrate();
            const auto hangingEdges = meshKernelState[meshKernelId].m_mesh2d->GetHangingEdges();
            numHangingEdges = static_cast<int>(hangingEdges.size());
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_hanging_edges(int meshKernelId, int* edges)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }
            const auto hangingEdges = meshKernelState[meshKernelId].m_mesh2d->GetHangingEdges();
            for (size_t i = 0; i < hangingEdges.size(); ++i)
            {
                edges[i] = static_cast<int>(hangingEdges[i]);
            }
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_delete_hanging_edges(int meshKernelId)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }
            meshKernelState[meshKernelId].m_mesh2d->DeleteHangingEdges();
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_compute_orthogonalization(int meshKernelId,
                                                             int projectToLandBoundaryOption,
                                                             const meshkernel::OrthogonalizationParameters& orthogonalizationParameters,
                                                             const GeometryList& selectingPolygon,
                                                             const GeometryList& landBoundaries)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }
            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() <= 0)
            {
                throw meshkernel::MeshKernelError("The 2d mesh contains no nodes.");
            }

            // build the selecting polygon
            auto const polygonNodes = ConvertGeometryListToPointVector(selectingPolygon);

            // build the land boundary
            auto const landBoundariesPoints = ConvertGeometryListToPointVector(landBoundaries);

            // Construct all dependencies
            auto const smoother = std::make_shared<meshkernel::Smoother>(meshKernelState[meshKernelId].m_mesh2d);
            auto const orthogonalizer = std::make_shared<meshkernel::Orthogonalizer>(meshKernelState[meshKernelId].m_mesh2d);
            auto const polygon = std::make_shared<meshkernel::Polygons>(polygonNodes, meshKernelState[meshKernelId].m_mesh2d->m_projection);
            auto const landBoundary = std::make_shared<meshkernel::LandBoundaries>(landBoundariesPoints, meshKernelState[meshKernelId].m_mesh2d, polygon);

            meshkernel::OrthogonalizationAndSmoothing ortogonalization(meshKernelState[meshKernelId].m_mesh2d,
                                                                       smoother,
                                                                       orthogonalizer,
                                                                       polygon,
                                                                       landBoundary,
                                                                       static_cast<meshkernel::LandBoundaries::ProjectToLandBoundaryOption>(projectToLandBoundaryOption),
                                                                       orthogonalizationParameters);
            ortogonalization.Initialize();
            ortogonalization.Compute();
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_initialize_orthogonalization(int meshKernelId,
                                                                int projectToLandBoundaryOption,
                                                                meshkernel::OrthogonalizationParameters& orthogonalizationParameters,
                                                                const GeometryList& selectingPolygon,
                                                                const GeometryList& landBoundaries)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() <= 0)
            {
                return lastExitCode;
            }

            // build the selecting polygon
            auto const polygonNodesVector = ConvertGeometryListToPointVector(selectingPolygon);

            // build the land boundary
            auto const landBoundariesNodeVector = ConvertGeometryListToPointVector(landBoundaries);

            // Construct all dependencies
            auto const smoother = std::make_shared<meshkernel::Smoother>(meshKernelState[meshKernelId].m_mesh2d);
            auto const orthogonalizer = std::make_shared<meshkernel::Orthogonalizer>(meshKernelState[meshKernelId].m_mesh2d);
            auto const polygon = std::make_shared<meshkernel::Polygons>(polygonNodesVector, meshKernelState[meshKernelId].m_mesh2d->m_projection);
            auto const landBoundary = std::make_shared<meshkernel::LandBoundaries>(landBoundariesNodeVector, meshKernelState[meshKernelId].m_mesh2d, polygon);

            meshKernelState[meshKernelId].m_meshOrthogonalization = std::make_shared<meshkernel::OrthogonalizationAndSmoothing>(meshKernelState[meshKernelId].m_mesh2d,
                                                                                                                                smoother,
                                                                                                                                orthogonalizer,
                                                                                                                                polygon,
                                                                                                                                landBoundary,
                                                                                                                                static_cast<meshkernel::LandBoundaries::ProjectToLandBoundaryOption>(projectToLandBoundaryOption),
                                                                                                                                orthogonalizationParameters);
            meshKernelState[meshKernelId].m_meshOrthogonalization->Initialize();
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_prepare_outer_iteration_orthogonalization(int meshKernelId)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() <= 0)
            {
                return lastExitCode;
            }

            meshKernelState[meshKernelId].m_meshOrthogonalization->PrepareOuterIteration();
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_compute_inner_ortogonalization_iteration(int meshKernelId)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() <= 0)
            {
                return lastExitCode;
            }

            meshKernelState[meshKernelId].m_meshOrthogonalization->Solve();
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_finalize_inner_ortogonalization_iteration(int meshKernelId)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() <= 0)
            {
                return lastExitCode;
            }

            meshKernelState[meshKernelId].m_meshOrthogonalization->FinalizeOuterIteration();
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_delete_orthogonalization(int meshKernelId)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() <= 0)
            {
                return lastExitCode;
            }

            meshKernelState[meshKernelId].m_meshOrthogonalization.reset();
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_orthogonality(int meshKernelId, GeometryList& geometryList)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() <= 0)
            {
                return lastExitCode;
            }

            const auto result = meshKernelState[meshKernelId].m_mesh2d->GetOrthogonality();

            if (static_cast<size_t>(geometryList.num_coordinates) != result.size())
            {
                throw meshkernel::MeshKernelError("The value array has not the same size of the result array storing the orthogonality values at the edges");
            }

            for (auto i = 0; i < geometryList.num_coordinates; ++i)
            {
                geometryList.values[i] = result[i];
            }
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_smoothness(int meshKernelId, GeometryList& geometryList)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() <= 0)
            {
                return lastExitCode;
            }

            const auto result = meshKernelState[meshKernelId].m_mesh2d->GetSmoothness();

            for (auto i = 0; i < geometryList.num_coordinates; ++i)
            {
                geometryList.values[i] = result[i];
            }
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_get_splines(const GeometryList& geometryListIn,
                                        GeometryList& geometryListOut,
                                        int numberOfPointsBetweenNodes)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (geometryListIn.num_coordinates == 0)
            {
                throw meshkernel::MeshKernelError("The number of coordinates of the given geometry is zero.");
            }

            std::vector<meshkernel::Point> splines(geometryListIn.num_coordinates);
            for (auto i = 0; i < geometryListIn.num_coordinates; ++i)
            {
                splines[i].x = geometryListIn.coordinates_x[i];
                splines[i].y = geometryListIn.coordinates_y[i];
            }

            const auto indices = FindIndices(splines, 0, static_cast<meshkernel::UInt>(splines.size()), meshkernel::constants::missing::doubleValue);
            const auto numSplines = static_cast<meshkernel::UInt>(indices.size());

            int index = 0;
            for (meshkernel::UInt s = 0; s < numSplines; s++)
            {
                const auto& [startIndex, endIndex] = indices[s];
                std::vector<meshkernel::Point> coordinates(splines.begin() + startIndex, splines.begin() + static_cast<int>(endIndex) + 1);
                const int numNodes = static_cast<int>(endIndex) - static_cast<int>(startIndex) + 1;
                const auto coordinatesDerivatives = meshkernel::SplineAlgorithms::SecondOrderDerivative(coordinates, 0, static_cast<meshkernel::UInt>(coordinates.size()) - 1);

                for (auto n = 0; n < numNodes - 1; n++)
                {
                    // Add the first point
                    geometryListOut.coordinates_x[index] = coordinates[n].x;
                    geometryListOut.coordinates_y[index] = coordinates[n].y;
                    index++;
                    for (auto p = 1; p <= numberOfPointsBetweenNodes; p++)
                    {
                        const double pointAdimensionalCoordinate = n + static_cast<double>(p) / static_cast<double>(numberOfPointsBetweenNodes + 1);
                        const auto pointCoordinate = ComputePointOnSplineAtAdimensionalDistance(coordinates, coordinatesDerivatives, pointAdimensionalCoordinate);
                        if (!pointCoordinate.IsValid())
                        {
                            break;
                        }

                        geometryListOut.coordinates_x[index] = pointCoordinate.x;
                        geometryListOut.coordinates_y[index] = pointCoordinate.y;
                        geometryListOut.values[index] = meshkernel::constants::missing::doubleValue;
                        index++;
                    }
                }

                geometryListOut.coordinates_x[index] = coordinates.back().x;
                geometryListOut.coordinates_y[index] = coordinates.back().y;
                geometryListOut.values[index] = meshkernel::constants::missing::doubleValue;
                index++;

                if (s != numSplines - 1)
                {
                    geometryListOut.coordinates_x[index] = meshkernel::constants::missing::doubleValue;
                    geometryListOut.coordinates_y[index] = meshkernel::constants::missing::doubleValue;
                    geometryListOut.values[index] = meshkernel::constants::missing::doubleValue;
                    index++;
                }
            }
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_make_triangular_mesh_from_polygon(int meshKernelId, const GeometryList& polygonPoints)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            auto const polygonPointsVector = ConvertGeometryListToPointVector(polygonPoints);

            const meshkernel::Polygons polygon(polygonPointsVector, meshKernelState[meshKernelId].m_mesh2d->m_projection);

            // generate samples in all polygons
            auto const generatedPoints = polygon.ComputePointsInPolygons();

            const meshkernel::Mesh2D mesh(generatedPoints[0], polygon, meshKernelState[meshKernelId].m_mesh2d->m_projection);
            *meshKernelState[meshKernelId].m_mesh2d += mesh;
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_make_triangular_mesh_from_samples(int meshKernelId, const GeometryList& samples)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }
            auto sampleVector = ConvertGeometryListToPointVector(samples);

            meshkernel::Polygons polygon;
            const meshkernel::Mesh2D mesh(sampleVector, polygon, meshKernelState[meshKernelId].m_mesh2d->m_projection);
            *meshKernelState[meshKernelId].m_mesh2d += mesh;
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_make_rectangular_mesh(int meshKernelId,
                                                         const meshkernel::MakeGridParameters& makeGridParameters)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            const auto projection = meshKernelState[meshKernelId].m_projection;
            const auto curvilinearGrid = CreateRectangularCurvilinearGrid(makeGridParameters, projection);

            auto const [nodes, edges, gridIndices] = curvilinearGrid.ConvertCurvilinearToNodesAndEdges();
            *meshKernelState[meshKernelId].m_mesh2d += meshkernel::Mesh2D(edges, nodes, projection);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_make_rectangular_mesh_from_polygon(int meshKernelId,
                                                                      const meshkernel::MakeGridParameters& makeGridParameters,
                                                                      const GeometryList& geometryList)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            const auto projection = meshKernelState[meshKernelId].m_projection;
            const auto curvilinearGrid = CreateRectangularCurvilinearGridFromPolygons(makeGridParameters, geometryList, projection);

            auto const [nodes, edges, gridIndices] = curvilinearGrid.ConvertCurvilinearToNodesAndEdges();
            *meshKernelState[meshKernelId].m_mesh2d += meshkernel::Mesh2D(edges, nodes, projection);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_make_rectangular_mesh_on_extension(int meshKernelId,
                                                                      const meshkernel::MakeGridParameters& makeGridParameters)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            const auto projection = meshKernelState[meshKernelId].m_projection;
            auto const curvilinearGrid = CreateRectangularCurvilinearGridOnExtension(makeGridParameters, projection);

            auto const [nodes, edges, gridIndices] = curvilinearGrid.ConvertCurvilinearToNodesAndEdges();
            *meshKernelState[meshKernelId].m_mesh2d += meshkernel::Mesh2D(edges, nodes, meshKernelState[meshKernelId].m_curvilinearGrid->m_projection);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_mesh_boundaries_as_polygons(int meshKernelId, GeometryList& boundaryPolygons)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            const std::vector<meshkernel::Point> polygonNodes;
            const auto meshBoundaryPolygon = meshKernelState[meshKernelId].m_mesh2d->MeshBoundaryToPolygon(polygonNodes);

            ConvertPointVectorToGeometryList(meshBoundaryPolygon, boundaryPolygons);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_count_mesh_boundaries_as_polygons(int meshKernelId, int& numberOfPolygonNodes)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            const std::vector<meshkernel::Point> polygonNodes;
            const auto meshBoundaryPolygon = meshKernelState[meshKernelId].m_mesh2d->MeshBoundaryToPolygon(polygonNodes);
            numberOfPolygonNodes = static_cast<int>(meshBoundaryPolygon.size()); // last value is a separator
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_intersections_from_polygon(int meshKernelId,
                                                              const GeometryList& boundaryPolygon,
                                                              int* edgeNodes,
                                                              int* edgeIndex,
                                                              double* edgeDistances,
                                                              double* segmentDistances,
                                                              int* segmentIndexes,
                                                              int* faceIndexes,
                                                              int* faceNumEdges,
                                                              int* faceEdgeIndex)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            auto const boundaryPolygonPoints = ConvertGeometryListToPointVector(boundaryPolygon);

            const meshkernel::Polygons polygons(boundaryPolygonPoints, meshKernelState[meshKernelId].m_projection);

            meshkernel::Mesh2DIntersections mesh2DIntersections(*meshKernelState[meshKernelId].m_mesh2d);
            mesh2DIntersections.Compute(polygons);
            auto edgeIntersections = mesh2DIntersections.EdgeIntersections();
            auto faceIntersections = mesh2DIntersections.FaceIntersections();

            meshkernel::Mesh2DIntersections::sortAndEraseIntersections(edgeIntersections);
            meshkernel::Mesh2DIntersections::sortAndEraseIntersections(faceIntersections);

            int edgeNodesCount = 0;
            int edgeCount = 0;
            for (size_t i = 0; i < edgeIntersections.size(); ++i)
            {
                const auto& edgeIntersection = edgeIntersections[i];

                // edge information must be stored only once
                edgeNodes[edgeNodesCount] = static_cast<int>(edgeIntersection.edgeFirstNode);
                edgeNodesCount++;
                edgeNodes[edgeNodesCount] = static_cast<int>(edgeIntersection.edgeSecondNode);
                edgeNodesCount++;

                // the edge count
                edgeDistances[edgeCount] = edgeIntersection.edgeDistance;
                segmentIndexes[edgeCount] = edgeIntersection.polylineSegmentIndex;
                segmentDistances[edgeCount] = edgeIntersection.adimensionalPolylineSegmentDistance;
                edgeIndex[edgeCount] = static_cast<int>(edgeIntersection.edgeIndex);
                edgeCount++;
            }

            int faceEdgesCount = 0;
            int faceCount = 0;
            for (const auto& intersection : faceIntersections)
            {
                faceNumEdges[faceCount] = static_cast<int>(intersection.edgeIndices.size());
                faceCount++;
                for (size_t i = 0; i < intersection.edgeIndices.size(); ++i)
                {
                    faceIndexes[faceEdgesCount] = static_cast<int>(intersection.faceIndex);
                    faceEdgeIndex[faceEdgesCount] = static_cast<int>(intersection.edgeIndices[i]);
                    faceEdgesCount++;
                }
            }
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_polygon_refine(int meshKernelId, const GeometryList& polygonToRefine, int firstNodeIndex, int secondNodeIndex, double targetEdgeLength, GeometryList& refinedPolygon)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }
            auto const polygonVector = ConvertGeometryListToPointVector(polygonToRefine);

            const meshkernel::Polygons polygon(polygonVector, meshKernelState[meshKernelId].m_mesh2d->m_projection);
            auto const refinementResult = polygon.RefineFirstPolygon(firstNodeIndex, secondNodeIndex, targetEdgeLength);

            ConvertPointVectorToGeometryList(refinementResult, refinedPolygon);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_polygon_count_refine(int meshKernelId,
                                                 const GeometryList& polygonToRefine,
                                                 int firstIndex,
                                                 int secondIndex,
                                                 double distance,
                                                 int& numberOfPolygonNodes)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            auto const polygonVector = ConvertGeometryListToPointVector(polygonToRefine);

            const meshkernel::Polygons polygon(polygonVector, meshKernelState[meshKernelId].m_mesh2d->m_projection);

            const auto refinedPolygon = polygon.RefineFirstPolygon(firstIndex, secondIndex, distance);

            numberOfPolygonNodes = int(refinedPolygon.size());
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_merge_nodes(int meshKernelId, const GeometryList& geometryListIn)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            auto const polygonVector = ConvertGeometryListToPointVector(geometryListIn);

            const meshkernel::Polygons polygon(polygonVector, meshKernelState[meshKernelId].m_mesh2d->m_projection);
            meshKernelState[meshKernelId].m_mesh2d->ComputeEdgesLengths();

            const auto minEdgeLength = meshKernelState[meshKernelId].m_mesh2d->ComputeMinEdgeLength(polygon);
            const auto searchRadius = std::max(1e-6, minEdgeLength * 0.1);
            meshKernelState[meshKernelId].m_mesh2d->MergeNodesInPolygon(polygon, searchRadius);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_merge_nodes_with_merging_distance(int meshKernelId, const GeometryList& geometryListIn, double mergingDistance)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            auto const polygonVector = ConvertGeometryListToPointVector(geometryListIn);

            const meshkernel::Polygons polygon(polygonVector, meshKernelState[meshKernelId].m_mesh2d->m_projection);

            meshKernelState[meshKernelId].m_mesh2d->MergeNodesInPolygon(polygon, mergingDistance);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_merge_two_nodes(int meshKernelId, int firstNode, int secondNode)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }
            meshKernelState[meshKernelId].m_mesh2d->MergeTwoNodes(firstNode, secondNode);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_nodes_in_polygons(int meshKernelId,
                                                         const GeometryList& geometryListIn,
                                                         int inside,
                                                         int* selectedNodes)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            auto const polygonVector = ConvertGeometryListToPointVector(geometryListIn);

            const meshkernel::Polygons polygon(polygonVector, meshKernelState[meshKernelId].m_mesh2d->m_projection);

            const bool selectInside = inside == 1 ? true : false;
            const auto nodeMask = meshKernelState[meshKernelId].m_mesh2d->NodeMaskFromPolygon(polygon, selectInside);

            int index = 0;
            for (size_t i = 0; i < meshKernelState[meshKernelId].m_mesh2d->GetNumNodes(); ++i)
            {
                if (nodeMask[i] > 0)
                {
                    selectedNodes[index] = static_cast<int>(i);
                    index++;
                }
            }
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_count_nodes_in_polygons(int meshKernelId,
                                                           const GeometryList& geometryListIn,
                                                           int inside,
                                                           int& numberOfMeshNodes)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }
            auto const polygonVector = ConvertGeometryListToPointVector(geometryListIn);

            const meshkernel::Polygons polygon(polygonVector, meshKernelState[meshKernelId].m_mesh2d->m_projection);

            const bool selectInside = inside == 1 ? true : false;
            const auto nodeMask = meshKernelState[meshKernelId].m_mesh2d->NodeMaskFromPolygon(polygon, selectInside);

            numberOfMeshNodes = 0;
            for (size_t i = 0; i < meshKernelState[meshKernelId].m_mesh2d->GetNumNodes(); ++i)
            {
                if (nodeMask[i] > 0)
                {
                    numberOfMeshNodes++;
                }
            }
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_insert_edge(int meshKernelId, int startNode, int endNode, int& new_edge_index)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            new_edge_index = static_cast<int>(meshKernelState[meshKernelId].m_mesh2d->ConnectNodes(startNode, endNode));
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_insert_node(int meshKernelId, double xCoordinate, double yCoordinate, int& nodeIndex)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            meshkernel::Point const nodeCoordinateVector{xCoordinate, yCoordinate};

            nodeIndex = static_cast<int>(meshKernelState[meshKernelId].m_mesh2d->InsertNode(nodeCoordinateVector));
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_delete_node(int meshKernelId, int nodeIndex)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            meshKernelState[meshKernelId].m_mesh2d->DeleteNode(nodeIndex);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_move_node(int meshKernelId, double xCoordinate, double yCoordinate, int nodeIndex)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            meshkernel::Point newPosition{xCoordinate, yCoordinate};

            meshKernelState[meshKernelId].m_mesh2d->MoveNode(newPosition, nodeIndex);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_delete_edge(int meshKernelId,
                                               double xCoordinate,
                                               double yCoordinate,
                                               double xLowerLeftBoundingBox,
                                               double yLowerLeftBoundingBox,
                                               double xUpperRightBoundingBox,
                                               double yUpperRightBoundingBox)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            meshkernel::Point const point{xCoordinate, yCoordinate};

            meshkernel::BoundingBox boundingBox{{xLowerLeftBoundingBox, yLowerLeftBoundingBox}, {xUpperRightBoundingBox, yUpperRightBoundingBox}};

            meshKernelState[meshKernelId].m_mesh2d->BuildTree(meshkernel::Mesh::Location::Edges, boundingBox);
            const auto edgeIndex = meshKernelState[meshKernelId].m_mesh2d->FindEdgeCloseToAPoint(point);

            meshKernelState[meshKernelId].m_mesh2d->DeleteEdge(edgeIndex);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_edge(int meshKernelId,
                                            double xCoordinate,
                                            double yCoordinate,
                                            double xLowerLeftBoundingBox,
                                            double yLowerLeftBoundingBox,
                                            double xUpperRightBoundingBox,
                                            double yUpperRightBoundingBox,
                                            int& edgeIndex)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            meshkernel::Point const point{xCoordinate, yCoordinate};

            meshkernel::BoundingBox boundingBox{{xLowerLeftBoundingBox, yLowerLeftBoundingBox}, {xUpperRightBoundingBox, yUpperRightBoundingBox}};

            meshKernelState[meshKernelId].m_mesh2d->BuildTree(meshkernel::Mesh::Location::Edges, boundingBox);

            edgeIndex = static_cast<int>(meshKernelState[meshKernelId].m_mesh2d->FindEdgeCloseToAPoint(point));
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_polygon_get_offset(int meshKernelId, const GeometryList& geometryListIn, int inWard, double distance, GeometryList& geometryListOut)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            auto const polygonVector = ConvertGeometryListToPointVector(geometryListIn);

            const meshkernel::Polygons polygon(polygonVector, meshKernelState[meshKernelId].m_mesh2d->m_projection);

            const bool inWardBool = inWard == 1;
            const auto newPolygon = polygon.OffsetCopy(distance, inWardBool);

            ConvertPointVectorToGeometryList(newPolygon.GatherAllEnclosureNodes(), geometryListOut);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_polygon_count_offset(int meshKernelId, const GeometryList& geometryListIn, int innerPolygon, double distance, int& numberOfPolygonNodes)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }
            auto polygonPoints = ConvertGeometryListToPointVector(geometryListIn);

            const bool innerPolygonBool = innerPolygon == 1;
            const meshkernel::Polygons polygon(polygonPoints, meshKernelState[meshKernelId].m_mesh2d->m_projection);
            const auto newPolygon = polygon.OffsetCopy(distance, innerPolygonBool);

            numberOfPolygonNodes = static_cast<int>(newPolygon.GetNumNodes());
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_refine_based_on_samples(int meshKernelId,
                                                           const GeometryList& samples,
                                                           double relativeSearchRadius,
                                                           int minimumNumSamples,
                                                           const meshkernel::MeshRefinementParameters& meshRefinementParameters)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }
            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() <= 0)
            {
                throw meshkernel::ConstraintError("The selected mesh has no nodes.");
            }

            auto samplesVector = ConvertGeometryListToSampleVector(samples);

            meshkernel::AveragingInterpolation::Method averagingMethod;
            if (meshRefinementParameters.refinement_type == 1)
            {
                averagingMethod = meshkernel::AveragingInterpolation::Method::MinAbsValue;
            }
            if (meshRefinementParameters.refinement_type == 2)
            {
                averagingMethod = meshkernel::AveragingInterpolation::Method::Max;
            }

            const bool refineOutsideFace = meshRefinementParameters.account_for_samples_outside == 1 ? true : false;
            const bool transformSamples = meshRefinementParameters.refinement_type == 2 ? true : false;

            const auto averaging = std::make_shared<meshkernel::AveragingInterpolation>(*meshKernelState[meshKernelId].m_mesh2d,
                                                                                        samplesVector,
                                                                                        averagingMethod,
                                                                                        meshkernel::Mesh::Location::Faces,
                                                                                        relativeSearchRadius,
                                                                                        refineOutsideFace,
                                                                                        transformSamples,
                                                                                        static_cast<meshkernel::UInt>(minimumNumSamples));

            meshkernel::MeshRefinement meshRefinement(meshKernelState[meshKernelId].m_mesh2d, averaging, meshRefinementParameters);
            meshRefinement.Compute();
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_refine_based_on_gridded_samples(int meshKernelId,
                                                                   const GriddedSamples& griddedSamples,
                                                                   const meshkernel::MeshRefinementParameters& meshRefinementParameters,
                                                                   bool useNodalRefinement)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }
            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() <= 0)
            {
                throw meshkernel::ConstraintError("The selected mesh has no nodes.");
            }

            std::vector values(griddedSamples.num_x * griddedSamples.num_y, 0.0);
            for (size_t i = 0; i < values.size(); ++i)
            {
                values[i] = griddedSamples.values[i];
            }

            std::shared_ptr<meshkernel::MeshInterpolation> interpolant;
            if (griddedSamples.x_coordinates == nullptr && griddedSamples.y_coordinates == nullptr)
            {
                meshkernel::Point origin{griddedSamples.x_origin, griddedSamples.y_origin};
                interpolant = std::make_shared<meshkernel::BilinearInterpolationOnGriddedSamples>(*meshKernelState[meshKernelId].m_mesh2d,
                                                                                                  griddedSamples.num_x,
                                                                                                  griddedSamples.num_y,
                                                                                                  origin,
                                                                                                  griddedSamples.cell_size,
                                                                                                  values);
            }
            else
            {
                if (griddedSamples.x_coordinates == nullptr)
                {
                    throw meshkernel::MeshKernelError("griddedSamples.x_coordinates is nullptr");
                }

                if (griddedSamples.y_coordinates == nullptr)
                {
                    throw meshkernel::MeshKernelError("griddedSamples.y_coordinates is nullptr");
                }

                std::vector<double> xCoordinates(griddedSamples.num_x);
                for (size_t i = 0; i < xCoordinates.size(); ++i)
                {
                    xCoordinates[i] = griddedSamples.x_coordinates[i];
                }
                std::vector<double> yCoordinates(griddedSamples.num_y);
                for (size_t i = 0; i < yCoordinates.size(); ++i)
                {
                    yCoordinates[i] = griddedSamples.y_coordinates[i];
                }

                interpolant = std::make_shared<meshkernel::BilinearInterpolationOnGriddedSamples>(*meshKernelState[meshKernelId].m_mesh2d,
                                                                                                  xCoordinates,
                                                                                                  yCoordinates,
                                                                                                  values);
            }

            meshkernel::MeshRefinement meshRefinement(meshKernelState[meshKernelId].m_mesh2d, interpolant, meshRefinementParameters, useNodalRefinement);
            meshRefinement.Compute();
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_refine_based_on_polygon(int meshKernelId,
                                                           const GeometryList& geometryList,
                                                           const meshkernel::MeshRefinementParameters& meshRefinementParameters)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }
            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() <= 0)
            {
                throw meshkernel::ConstraintError("The selected mesh has no nodes.");
            }

            auto points = ConvertGeometryListToPointVector(geometryList);

            const meshkernel::Polygons polygon(points, meshKernelState[meshKernelId].m_mesh2d->m_projection);

            meshkernel::MeshRefinement meshRefinement(meshKernelState[meshKernelId].m_mesh2d, polygon, meshRefinementParameters);
            meshRefinement.Compute();
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_remove_disconnected_regions(int meshKernelId)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }
            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() <= 0)
            {
                throw meshkernel::ConstraintError("The selected mesh has no nodes.");
            }

            meshkernel::RemoveDisconnectedRegions removeDisconnectedRegions;
            removeDisconnectedRegions.Compute(*meshKernelState[meshKernelId].m_mesh2d);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_node_index(int meshKernelId,
                                                  double xCoordinate,
                                                  double yCoordinate,
                                                  double searchRadius,
                                                  double xLowerLeftBoundingBox,
                                                  double yLowerLeftBoundingBox,
                                                  double xUpperRightBoundingBox,
                                                  double yUpperRightBoundingBox,
                                                  int& nodeIndex)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }
            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() <= 0)
            {
                throw meshkernel::ConstraintError("The selected mesh has no nodes.");
            }

            meshkernel::Point const point{xCoordinate, yCoordinate};
            meshkernel::BoundingBox boundingBox{{xLowerLeftBoundingBox, yLowerLeftBoundingBox}, {xUpperRightBoundingBox, yUpperRightBoundingBox}};
            meshKernelState[meshKernelId].m_mesh2d->BuildTree(meshkernel::Mesh::Location::Nodes, boundingBox);
            nodeIndex = static_cast<int>(meshKernelState[meshKernelId].m_mesh2d->FindNodeCloseToAPoint(point, searchRadius));
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_closest_node(int meshKernelId,
                                                    double xCoordinateIn,
                                                    double yCoordinateIn,
                                                    double searchRadius,
                                                    double xLowerLeftBoundingBox,
                                                    double yLowerLeftBoundingBox,
                                                    double xUpperRightBoundingBox,
                                                    double yUpperRightBoundingBox,
                                                    double& xCoordinateOut,
                                                    double& yCoordinateOut)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            int nodeIndex;
            lastExitCode = static_cast<meshkernel::ExitCode>(mkernel_mesh2d_get_node_index(meshKernelId,
                                                                                           xCoordinateIn,
                                                                                           yCoordinateIn,
                                                                                           searchRadius,
                                                                                           xLowerLeftBoundingBox,
                                                                                           yLowerLeftBoundingBox,
                                                                                           xUpperRightBoundingBox,
                                                                                           yUpperRightBoundingBox,
                                                                                           nodeIndex));

            // Set the node coordinate
            const auto foundNode = meshKernelState[meshKernelId].m_mesh2d->m_nodes[nodeIndex];
            xCoordinateOut = foundNode.x;
            yCoordinateOut = foundNode.y;
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_polygon_get_included_points(int meshKernelId, const GeometryList& selectingPolygon, const GeometryList& polygonToSelect, GeometryList& selectionResults)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }
            auto const polygonVector = ConvertGeometryListToPointVector(selectingPolygon);

            auto const points = ConvertGeometryListToPointVector(polygonToSelect);

            const meshkernel::Polygons localPolygon(polygonVector, meshKernelState[meshKernelId].m_mesh2d->m_projection);

            for (size_t i = 0; i < points.size(); ++i)
            {
                selectionResults.values[i] = localPolygon.IsPointInPolygon(points[i], 0) ? 1.0 : 0.0;
            }
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_rotate(int meshKernelId, double centreX, double centreY, double theta)
    {
        lastExitCode = meshkernel::ExitCode::Success;

        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            meshkernel::RigidBodyTransformation transformation;
            meshkernel::Translation translation;

            // Translate centre of rotation to origin
            transformation.compose(meshkernel::Translation(meshkernel::Vector(-centreX, -centreY)));
            // Add rotation
            transformation.compose(meshkernel::Rotation(theta));
            // Translate origin back to centre of rotation
            transformation.compose(meshkernel::Translation(meshkernel::Vector(centreX, centreY)));

            meshkernel::MeshTransformation::Compute(*meshKernelState[meshKernelId].m_mesh2d, transformation);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_translate(int meshKernelId, double translationX, double translationY)
    {
        lastExitCode = meshkernel::ExitCode::Success;

        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            meshkernel::Translation translation(meshkernel::Vector(translationX, translationY));
            meshkernel::MeshTransformation::Compute(*meshKernelState[meshKernelId].m_mesh2d, translation);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_splines_snap_to_landboundary(int meshKernelId,
                                                         const GeometryList& land,
                                                         GeometryList& splines,
                                                         int startSplineIndex,
                                                         int endSplineIndex)
    {
        lastExitCode = meshkernel::ExitCode::Success;

        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (startSplineIndex < 0)
            {
                throw meshkernel::ConstraintError("Start spline index is less than zero: {}", startSplineIndex);
            }

            if (startSplineIndex > endSplineIndex)
            {
                throw meshkernel::ConstraintError("Invalid spline range: {} > {}",
                                                  startSplineIndex,
                                                  endSplineIndex);
            }

            if (land.num_coordinates == 0)
            {
                throw meshkernel::MeshKernelError("Land boundary has no point values.");
            }

            if (land.coordinates_x == nullptr || land.coordinates_y == nullptr)
            {
                throw meshkernel::MeshKernelError("Land boundary data is null.");
            }

            if (splines.num_coordinates == 0)
            {
                throw meshkernel::MeshKernelError("Spline has no point values.");
            }

            if (splines.coordinates_x == nullptr || splines.coordinates_y == nullptr)
            {
                throw meshkernel::MeshKernelError("Spline data is null.");
            }

            if (startSplineIndex > splines.num_coordinates)
            {
                throw meshkernel::ConstraintError("Invalid spline range: start greater than number of spline coordinates {} > {}",
                                                  startSplineIndex,
                                                  splines.num_coordinates);
            }

            if (endSplineIndex >= splines.num_coordinates)
            {
                throw meshkernel::ConstraintError("Invalid spline range: end greater than number of spline coordinates {} >= {}",
                                                  endSplineIndex,
                                                  splines.num_coordinates);
            }

            std::vector<meshkernel::Point> landBoundaryPoints(ConvertGeometryListToPointVector(land));
            std::vector<meshkernel::Point> splinePoints(ConvertGeometryListToPointVector(splines));

            meshkernel::LandBoundary landBoundary(landBoundaryPoints);
            meshkernel::Splines splineValues(meshKernelState[meshKernelId].m_mesh2d->m_projection);

            splineValues.AddSpline(splinePoints, startSplineIndex, static_cast<meshkernel::UInt>(splinePoints.size()));

            //--------------------------------
            // Snap specified splines to the land boundary
            splineValues.SnapSpline(0, landBoundary);

            //--------------------------------
            // Now copy back the snapped spline values

            int splinePointIndex = startSplineIndex;

            // Now copy back to spline (geometry-list)
            for (int i = startSplineIndex; i <= endSplineIndex; ++i)
            {
                const meshkernel::Point& splinePoint = splineValues.m_splineNodes[0][i];
                splines.coordinates_x[splinePointIndex] = splinePoint.x;
                splines.coordinates_y[splinePointIndex] = splinePoint.y;
                ++splinePointIndex;
            }
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_polygon_snap_to_landboundary(int meshKernelId,
                                                         const GeometryList& land,
                                                         GeometryList& polygon,
                                                         int startIndex,
                                                         int endIndex)
    {
        lastExitCode = meshkernel::ExitCode::Success;

        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (land.num_coordinates == 0)
            {
                throw meshkernel::MeshKernelError("Land boundary has no point values.");
            }

            if (land.coordinates_x == nullptr || land.coordinates_y == nullptr)
            {
                throw meshkernel::MeshKernelError("Land boundary data is null.");
            }

            if (polygon.num_coordinates == 0)
            {
                throw meshkernel::MeshKernelError("Polygon has no point values.");
            }

            if (polygon.coordinates_x == nullptr || polygon.coordinates_y == nullptr)
            {
                throw meshkernel::MeshKernelError("Polygon data is null.");
            }

            if (startIndex < 0 || endIndex < 0)
            {
                throw meshkernel::ConstraintError("Invalid polygon points range: startIndex and/or endIndex {} < 0 and/or {} < 0",
                                                  startIndex,
                                                  endIndex);
            }

            if (startIndex > endIndex)
            {
                throw meshkernel::ConstraintError("Invalid polygon points range: startIndex greater than endIndex {} > {}",
                                                  startIndex,
                                                  endIndex);
            }

            if (endIndex >= polygon.num_coordinates)
            {
                throw meshkernel::ConstraintError("Invalid polygon points range: endIndex greater than number of polygon coordinates {} >= {}",
                                                  endIndex,
                                                  polygon.num_coordinates);
            }

            std::vector<meshkernel::Point> landBoundaryPoints(ConvertGeometryListToPointVector(land));
            std::vector<meshkernel::Point> polygonPoints(ConvertGeometryListToPointVector(polygon));

            meshkernel::LandBoundary landBoundary(landBoundaryPoints);
            meshkernel::Polygons polygons(polygonPoints, meshKernelState[meshKernelId].m_mesh2d->m_projection);

            polygons.SnapToLandBoundary(landBoundary, startIndex, endIndex);
            const auto [enclosureIndex, polygonStartNode, polygonEndNode] = polygons.PolygonIndex(startIndex, endIndex);

            //--------------------------------
            // Now copy back the polygon values

            const std::vector<meshkernel::Point>& snappedPolygonPoints = polygons.Enclosure(enclosureIndex).Outer().Nodes();

            for (int i = startIndex; i <= endIndex; ++i)
            {
                polygon.coordinates_x[i] = snappedPolygonPoints[i].x;
                polygon.coordinates_y[i] = snappedPolygonPoints[i].y;
            }
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_flip_edges(int meshKernelId,
                                              int isTriangulationRequired,
                                              int projectToLandBoundaryRequired,
                                              const GeometryList& selectingPolygon,
                                              const GeometryList& landBoundaries)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            // build the selecting polygon
            auto const polygonNodesVector = ConvertGeometryListToPointVector(selectingPolygon);

            // build the land boundary
            auto const landBoundariesNodeVector = ConvertGeometryListToPointVector(landBoundaries);

            // construct all dependencies
            auto const polygon = std::make_shared<meshkernel::Polygons>(polygonNodesVector, meshKernelState[meshKernelId].m_mesh2d->m_projection);
            auto const landBoundary = std::make_shared<meshkernel::LandBoundaries>(landBoundariesNodeVector, meshKernelState[meshKernelId].m_mesh2d, polygon);
            bool const triangulateFaces = isTriangulationRequired == 0 ? false : true;
            bool const projectToLandBoundary = projectToLandBoundaryRequired == 0 ? false : true;

            const meshkernel::FlipEdges flipEdges(meshKernelState[meshKernelId].m_mesh2d, landBoundary, triangulateFaces, projectToLandBoundary);

            flipEdges.Compute();
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_count_small_flow_edge_centers(int meshKernelId, double smallFlowEdgesLengthThreshold, int& numSmallFlowEdges)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }
            const auto edgesCrossingSmallFlowEdges = meshKernelState[meshKernelId].m_mesh2d->GetEdgesCrossingSmallFlowEdges(smallFlowEdgesLengthThreshold);
            const auto smallFlowEdgeCenters = meshKernelState[meshKernelId].m_mesh2d->GetFlowEdgesCenters(edgesCrossingSmallFlowEdges);

            numSmallFlowEdges = static_cast<int>(smallFlowEdgeCenters.size());
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_small_flow_edge_centers(int meshKernelId, double smallFlowEdgesThreshold, GeometryList& result)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            const auto edgesCrossingSmallFlowEdges = meshKernelState[meshKernelId].m_mesh2d->GetEdgesCrossingSmallFlowEdges(smallFlowEdgesThreshold);
            const auto smallFlowEdgeCenters = meshKernelState[meshKernelId].m_mesh2d->GetFlowEdgesCenters(edgesCrossingSmallFlowEdges);

            ConvertPointVectorToGeometryList(smallFlowEdgeCenters, result);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_get_error(char* errorMessage)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        std::memcpy(errorMessage, exceptionMessage, sizeof exceptionMessage);
        return lastExitCode;
    }

    MKERNEL_API int mkernel_get_exit_code_success(int& exitCode)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        exitCode = meshkernel::ExitCode::Success;
        return lastExitCode;
    }

    MKERNEL_API int mkernel_get_exit_code_meshkernel_error(int& exitCode)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        exitCode = meshkernel::ExitCode::MeshKernelErrorCode;
        return lastExitCode;
    }

    MKERNEL_API int mkernel_get_exit_code_not_implemented_error(int& exitCode)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        exitCode = meshkernel::ExitCode::NotImplementedErrorCode;
        return lastExitCode;
    }

    MKERNEL_API int mkernel_get_exit_code_algorithm_error(int& exitCode)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        exitCode = meshkernel::ExitCode::AlgorithmErrorCode;
        return lastExitCode;
    }

    MKERNEL_API int mkernel_get_exit_code_constraint_error(int& exitCode)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        exitCode = meshkernel::ExitCode::ConstraintErrorCode;
        return lastExitCode;
    }

    MKERNEL_API int mkernel_get_exit_code_mesh_geometry_error(int& exitCode)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        exitCode = meshkernel::ExitCode::MeshGeometryErrorCode;
        return lastExitCode;
    }

    MKERNEL_API int mkernel_get_exit_code_linear_algebra_error(int& exitCode)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        exitCode = meshkernel::ExitCode::LinearAlgebraErrorCode;
        return lastExitCode;
    }

    MKERNEL_API int mkernel_get_exit_code_range_error(int& exitCode)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        exitCode = meshkernel::ExitCode::RangeErrorCode;
        return lastExitCode;
    }

    MKERNEL_API int mkernel_get_exit_code_stdlib_exception(int& exitCode)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        exitCode = meshkernel::ExitCode::StdLibExceptionCode;
        return lastExitCode;
    }

    MKERNEL_API int mkernel_get_exit_code_unknown_exception(int& exitCode)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        exitCode = meshkernel::ExitCode::UnknownExceptionCode;
        return lastExitCode;
    }

    MKERNEL_API int mkernel_get_version(char* version)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        std::memcpy(version, versionString, sizeof versionString);
        return lastExitCode;
    }

    MKERNEL_API int mkernel_get_geometry_error(int& meshIndex, int& meshLocation)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        meshIndex = static_cast<int>(invalidMeshIndex);
        meshLocation = static_cast<int>(invalidMeshLocation);
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_count_obtuse_triangles(int meshKernelId, int& numObtuseTriangles)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            const auto obtuseTriangles = meshKernelState[meshKernelId].m_mesh2d->GetObtuseTrianglesCenters();

            numObtuseTriangles = static_cast<int>(obtuseTriangles.size());
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_obtuse_triangles_mass_centers(int meshKernelId, GeometryList& result)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            const auto obtuseTriangles = meshKernelState[meshKernelId].m_mesh2d->GetObtuseTrianglesCenters();

            ConvertPointVectorToGeometryList(obtuseTriangles, result);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_delete_small_flow_edges_and_small_triangles(int meshKernelId, double smallFlowEdgesThreshold, double minFractionalAreaTriangles)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            meshKernelState[meshKernelId].m_mesh2d->DeleteSmallFlowEdges(smallFlowEdgesThreshold);
            meshKernelState[meshKernelId].m_mesh2d->DeleteSmallTrianglesAtBoundaries(minFractionalAreaTriangles);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_contacts_compute_single(int meshKernelId,
                                                    const int* oneDNodeMask,
                                                    const GeometryList& polygons,
                                                    double projectionFactor)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            // Convert 1D node mask from int** to vector<bool>
            auto num1DNodes = meshKernelState[meshKernelId].m_mesh1d->GetNumNodes();
            auto meshKernel1DNodeMask = ConvertIntegerArrayToBoolVector(oneDNodeMask,
                                                                        num1DNodes);

            // Convert polygon date from GeometryList to Polygons
            auto polygonPoints = ConvertGeometryListToPointVector(polygons);
            const meshkernel::Polygons meshKernelPolygons(polygonPoints,
                                                          meshKernelState[meshKernelId].m_mesh2d->m_projection);

            // Execute
            meshKernelState[meshKernelId].m_contacts->ComputeSingleContacts(meshKernel1DNodeMask, meshKernelPolygons, projectionFactor);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_contacts_compute_multiple(int meshKernelId,
                                                      const int* oneDNodeMask)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            // Convert 1D node mask from int** to vector<bool>
            auto num1DNodes = meshKernelState[meshKernelId].m_mesh1d->GetNumNodes();
            auto meshKernel1DNodeMask = ConvertIntegerArrayToBoolVector(oneDNodeMask,
                                                                        num1DNodes);

            // Execute
            meshKernelState[meshKernelId].m_contacts->ComputeMultipleContacts(meshKernel1DNodeMask);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_contacts_compute_with_polygons(int meshKernelId,
                                                           const int* oneDNodeMask,
                                                           const GeometryList& polygons)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            // Convert 1D node mask from int* to vector<bool>
            auto num1DNodes = meshKernelState[meshKernelId].m_mesh1d->GetNumNodes();
            auto meshKernel1DNodeMask = ConvertIntegerArrayToBoolVector(oneDNodeMask,
                                                                        num1DNodes);

            // Convert polygon date from GeometryList to Polygons
            auto polygonPoints = ConvertGeometryListToPointVector(polygons);
            const meshkernel::Polygons meshKernelPolygons(polygonPoints,
                                                          meshKernelState[meshKernelId].m_mesh2d->m_projection);

            // Execute
            meshKernelState[meshKernelId].m_contacts->ComputeContactsWithPolygons(meshKernel1DNodeMask,
                                                                                  meshKernelPolygons);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }
    MKERNEL_API int mkernel_contacts_compute_with_points(int meshKernelId,
                                                         const int* oneDNodeMask,
                                                         const GeometryList& points)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            // Convert 1D node mask from int* to vector<bool>
            auto num1DNodes = meshKernelState[meshKernelId].m_mesh1d->GetNumNodes();
            auto meshKernel1DNodeMask = ConvertIntegerArrayToBoolVector(oneDNodeMask,
                                                                        num1DNodes);

            // Convert polygon date from GeometryList to Point vector
            auto meshKernelPoints = ConvertGeometryListToPointVector(points);
            // Execute
            meshKernelState[meshKernelId].m_contacts->ComputeContactsWithPoints(meshKernel1DNodeMask, meshKernelPoints);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_connect_meshes(int meshKernelId, const Mesh2D& mesh2d, double searchFraction)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            // convert raw arrays to containers
            const auto edges2d = meshkernel::ConvertToEdgeNodesVector(mesh2d.num_edges,
                                                                      mesh2d.edge_nodes);

            const auto nodes2d = meshkernel::ConvertToNodesVector(mesh2d.num_nodes,
                                                                  mesh2d.node_x,
                                                                  mesh2d.node_y);

            meshkernel::Mesh2D meshToConnect;

            if (mesh2d.num_faces > 0 && mesh2d.face_nodes != nullptr && mesh2d.nodes_per_face != nullptr)
            {

                const auto face_nodes = meshkernel::ConvertToFaceNodesVector(mesh2d.num_faces, mesh2d.face_nodes, mesh2d.nodes_per_face);

                std::vector<meshkernel::UInt> num_face_nodes;
                num_face_nodes.reserve(mesh2d.num_faces);

                for (auto n = 0; n < mesh2d.num_faces; n++)
                {
                    num_face_nodes.emplace_back(static_cast<meshkernel::UInt>(mesh2d.nodes_per_face[n]));
                }

                meshToConnect = meshkernel::Mesh2D(edges2d,
                                                   nodes2d,
                                                   face_nodes,
                                                   num_face_nodes,
                                                   meshKernelState[meshKernelId].m_projection);
            }
            else
            {
                // Do not change the pointer, just the object it is pointing to
                // Compute the faces
                meshToConnect = meshkernel::Mesh2D(edges2d, nodes2d, meshKernelState[meshKernelId].m_projection);
            }

            meshkernel::Mesh2D mergedMeshes = meshkernel::Mesh2D::Merge(*meshKernelState[meshKernelId].m_mesh2d, meshToConnect);
            meshkernel::ConnectMeshes connectMeshes;
            connectMeshes.Compute(mergedMeshes, searchFraction);
            *meshKernelState[meshKernelId].m_mesh2d = mergedMeshes;
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_contacts_compute_boundary(int meshKernelId,
                                                      const int* oneDNodeMask,
                                                      const GeometryList& polygons,
                                                      double searchRadius)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            // Convert 1D node mask from int* to vector<bool>
            auto const num1DNodes = meshKernelState[meshKernelId].m_mesh1d->GetNumNodes();
            auto const meshKernel1DNodeMask = ConvertIntegerArrayToBoolVector(oneDNodeMask, num1DNodes);

            // Convert polygon date from GeometryList to Polygons
            auto const polygonPoints = ConvertGeometryListToPointVector(polygons);
            const meshkernel::Polygons meshKernelPolygons(polygonPoints, meshKernelState[meshKernelId].m_mesh2d->m_projection);

            // Execute
            meshKernelState[meshKernelId].m_contacts->ComputeBoundaryContacts(meshKernel1DNodeMask, meshKernelPolygons, searchRadius);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_refine(int meshKernelId,
                                               double xLowerLeftCorner,
                                               double yLowerLeftCorner,
                                               double xUpperRightCorner,
                                               double yUpperRightCorner,
                                               int refinement)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }
            meshkernel::Point const firstPoint{xLowerLeftCorner, yLowerLeftCorner};
            meshkernel::Point const secondPoint{xUpperRightCorner, yUpperRightCorner};

            // Execute
            meshkernel::CurvilinearGridRefinement curvilinearGridRefinement(*meshKernelState[meshKernelId].m_curvilinearGrid, refinement);
            curvilinearGridRefinement.SetBlock(firstPoint, secondPoint);
            curvilinearGridRefinement.Compute();
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_derefine(int meshKernelId,
                                                 double xLowerLeftCorner,
                                                 double yLowerLeftCorner,
                                                 double xUpperRightCorner,
                                                 double yUpperRightCorner)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            meshkernel::Point const firstPoint{xLowerLeftCorner, yLowerLeftCorner};
            meshkernel::Point const secondPoint{xUpperRightCorner, yUpperRightCorner};

            // Execute
            meshkernel::CurvilinearGridDeRefinement curvilinearGridDeRefinement(*meshKernelState[meshKernelId].m_curvilinearGrid);

            curvilinearGridDeRefinement.SetBlock(firstPoint, secondPoint);
            curvilinearGridDeRefinement.Compute();
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_compute_transfinite_from_splines(int meshKernelId,
                                                                         const GeometryList& splines,
                                                                         const meshkernel::CurvilinearParameters& curvilinearParameters)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            // Use the default constructor, no instance present
            const auto meshKernelSplines = std::make_shared<meshkernel::Splines>(meshKernelState[meshKernelId].m_projection);
            SetSplines(splines, *meshKernelSplines);

            // Create algorithm and set the splines
            meshkernel::CurvilinearGridFromSplinesTransfinite curvilinearGridFromSplinesTransfinite(meshKernelSplines, curvilinearParameters);

            // Compute the curvilinear grid
            const auto curvilinearGrid = curvilinearGridFromSplinesTransfinite.Compute();

            // Set the state
            meshKernelState[meshKernelId].m_curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(curvilinearGrid);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_compute_transfinite_from_polygon(int meshKernelId,
                                                                         const GeometryList& polygons,
                                                                         int firstNode,
                                                                         int secondNode,
                                                                         int thirdNode,
                                                                         int useFourthSide)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            auto polygonPoints = ConvertGeometryListToPointVector(polygons);

            meshkernel::Polygon localPolygon(polygonPoints, meshKernelState[meshKernelId].m_projection);

            const meshkernel::CurvilinearGridFromPolygon curvilinearGridFromPolygon(localPolygon);

            const bool useFourthSideBool = useFourthSide == 1 ? true : false;
            const auto curvilinearGrid = curvilinearGridFromPolygon.Compute(firstNode, secondNode, thirdNode, useFourthSideBool);

            // set the curvilinear state
            meshKernelState[meshKernelId].m_curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(curvilinearGrid);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_compute_transfinite_from_triangle(int meshKernelId,
                                                                          const GeometryList& polygon,
                                                                          int firstNode,
                                                                          int secondNode,
                                                                          int thirdNode)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            auto polygonPoints = ConvertGeometryListToPointVector(polygon);

            meshkernel::Polygon localPolygon(polygonPoints, meshKernelState[meshKernelId].m_projection);

            const meshkernel::CurvilinearGridFromPolygon curvilinearGridFromPolygon(localPolygon);

            const auto curvilinearGrid = curvilinearGridFromPolygon.Compute(firstNode, secondNode, thirdNode);

            // set the curvilinear state
            meshKernelState[meshKernelId].m_curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(curvilinearGrid);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_compute_orthogonal_grid_from_splines(int meshKernelId,
                                                                             const GeometryList& geometryListIn,
                                                                             const meshkernel::CurvilinearParameters& curvilinearParameters,
                                                                             const meshkernel::SplinesToCurvilinearParameters& splinesToCurvilinearParameters)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            // use the default constructor, no instance present
            const auto spline = std::make_shared<meshkernel::Splines>(meshKernelState[meshKernelId].m_projection);
            SetSplines(geometryListIn, *spline);

            meshkernel::CurvilinearGridFromSplines curvilinearGridFromSplines(spline, curvilinearParameters, splinesToCurvilinearParameters);

            // set the curvilinear state
            meshKernelState[meshKernelId].m_curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(curvilinearGridFromSplines.Compute());
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_initialize_orthogonal_grid_from_splines(int meshKernelId,
                                                                                const GeometryList& geometryList,
                                                                                const meshkernel::CurvilinearParameters& curvilinearParameters,
                                                                                const meshkernel::SplinesToCurvilinearParameters& splinesToCurvilinearParameters)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            auto spline = std::make_shared<meshkernel::Splines>(meshKernelState[meshKernelId].m_projection);
            SetSplines(geometryList, *spline);

            meshKernelState[meshKernelId].m_curvilinearGridFromSplines = std::make_shared<meshkernel::CurvilinearGridFromSplines>(spline, curvilinearParameters, splinesToCurvilinearParameters);

            meshKernelState[meshKernelId].m_curvilinearGridFromSplines->Initialize();
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_iterate_orthogonal_grid_from_splines(int meshKernelId, int layer)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }
            if (meshKernelState[meshKernelId].m_curvilinearGridFromSplines == nullptr)
            {
                throw meshkernel::MeshKernelError("CurvilinearGridFromSplines not instantiated.");
            }

            meshKernelState[meshKernelId].m_curvilinearGridFromSplines->Iterate(layer);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_refresh_orthogonal_grid_from_splines(int meshKernelId)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGridFromSplines == nullptr)
            {
                throw meshkernel::MeshKernelError("CurvilinearGridFromSplines not instantiated.");
            }

            const auto curvilinearGrid = meshKernelState[meshKernelId].m_curvilinearGridFromSplines->ComputeCurvilinearGridFromGridPoints();

            meshKernelState[meshKernelId].m_curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(curvilinearGrid);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_delete_orthogonal_grid_from_splines(int meshKernelId)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGridFromSplines == nullptr)
            {
                throw meshkernel::MeshKernelError("CurvilinearGridFromSplines not instantiated.");
            }

            meshKernelState[meshKernelId].m_curvilinearGridFromSplines.reset();
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_compute_rectangular_grid(int meshKernelId,
                                                                 const meshkernel::MakeGridParameters& makeGridParameters)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            *meshKernelState[meshKernelId].m_curvilinearGrid = CreateRectangularCurvilinearGrid(makeGridParameters,
                                                                                                meshKernelState[meshKernelId].m_projection);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_compute_rectangular_grid_from_polygon(int meshKernelId,
                                                                              const meshkernel::MakeGridParameters& makeGridParameters,
                                                                              const GeometryList& geometryList)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            *meshKernelState[meshKernelId].m_curvilinearGrid = CreateRectangularCurvilinearGridFromPolygons(makeGridParameters,
                                                                                                            geometryList,
                                                                                                            meshKernelState[meshKernelId].m_projection);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_compute_rectangular_grid_on_extension(int meshKernelId,
                                                                              const meshkernel::MakeGridParameters& makeGridParameters)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            *meshKernelState[meshKernelId].m_curvilinearGrid = CreateRectangularCurvilinearGridOnExtension(makeGridParameters,
                                                                                                           meshKernelState[meshKernelId].m_projection);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_initialize_orthogonalize(int meshKernelId,
                                                                 const meshkernel::OrthogonalizationParameters& orthogonalizationParameters)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            meshKernelState[meshKernelId].m_curvilinearGridOrthogonalization = std::make_shared<meshkernel::CurvilinearGridOrthogonalization>(*meshKernelState[meshKernelId].m_curvilinearGrid,
                                                                                                                                              orthogonalizationParameters);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_set_block_orthogonalize(int meshKernelId,
                                                                double xLowerLeftCorner,
                                                                double yLowerLeftCorner,
                                                                double xUpperRightCorner,
                                                                double yUpperRightCorner)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel state does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGridOrthogonalization == nullptr)
            {
                throw meshkernel::MeshKernelError("CurvilinearGridOrthogonalization not instantiated.");
            }

            meshkernel::Point firstPoint{xLowerLeftCorner, yLowerLeftCorner};
            meshkernel::Point secondPoint{xUpperRightCorner, yUpperRightCorner};

            // Execute
            meshKernelState[meshKernelId].m_curvilinearGridOrthogonalization->SetBlock(firstPoint, secondPoint);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_set_frozen_lines_orthogonalize(int meshKernelId,
                                                                       double xFirstGridLineNode,
                                                                       double yFirstGridLineNode,
                                                                       double xSecondGridLineNode,
                                                                       double ySecondGridLineNode)

    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel state does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGridOrthogonalization == nullptr)
            {
                throw meshkernel::MeshKernelError("CurvilinearGridOrthogonalization not instantiated.");
            }

            meshkernel::Point const firstPoint{xFirstGridLineNode, yFirstGridLineNode};
            meshkernel::Point const secondPoint{xSecondGridLineNode, ySecondGridLineNode};

            // Execute
            meshKernelState[meshKernelId].m_curvilinearGridOrthogonalization->SetLine(firstPoint, secondPoint);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_orthogonalize(int meshKernelId)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel state does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGridOrthogonalization == nullptr)
            {
                throw meshkernel::MeshKernelError("CurvilinearGridOrthogonalization not instantiated.");
            }

            // Execute
            meshKernelState[meshKernelId].m_curvilinearGridOrthogonalization->Compute();
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_finalize_orthogonalize(int meshKernelId)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel state does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGridOrthogonalization == nullptr)
            {
                throw meshkernel::MeshKernelError("CurvilinearGridOrthogonalization not instantiated.");
            }

            meshKernelState[meshKernelId].m_curvilinearGridOrthogonalization.reset();
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_smoothing(int meshKernelId,
                                                  int smoothingIterations,
                                                  double xLowerLeftCorner,
                                                  double yLowerLeftCorner,
                                                  double xUpperRightCorner,
                                                  double yUpperRightCorner)

    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel state does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGrid == nullptr)
            {
                throw meshkernel::MeshKernelError("Not a valid curvilinear grid instance.");
            }

            if (!meshKernelState[meshKernelId].m_curvilinearGrid->IsValid())
            {
                throw meshkernel::MeshKernelError("Not valid curvilinear grid.");
            }

            const meshkernel::Point firstPoint{xLowerLeftCorner, yLowerLeftCorner};
            const meshkernel::Point secondPoint{xUpperRightCorner, yUpperRightCorner};

            // Execute
            meshkernel::CurvilinearGridSmoothing curvilinearGridSmoothing(*meshKernelState[meshKernelId].m_curvilinearGrid,
                                                                          static_cast<meshkernel::UInt>(smoothingIterations));

            curvilinearGridSmoothing.SetBlock(firstPoint, secondPoint);
            curvilinearGridSmoothing.Compute();
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_smoothing_directional(int meshKernelId,
                                                              int smoothingIterations,
                                                              double xFirstGridlineNode,
                                                              double yFirstGridlineNode,
                                                              double xSecondGridLineNode,
                                                              double ySecondGridLineNode,
                                                              double xLowerLeftCornerSmoothingArea,
                                                              double yLowerLeftCornerSmoothingArea,
                                                              double xUpperRightCornerSmootingArea,
                                                              double yUpperRightCornerSmootingArea)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel state does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGrid == nullptr)
            {
                throw meshkernel::MeshKernelError("Not a valid curvilinear grid instance.");
            }

            if (!meshKernelState[meshKernelId].m_curvilinearGrid->IsValid())
            {
                throw meshkernel::MeshKernelError("Not valid curvilinear grid.");
            }

            meshkernel::Point const firstNode{xFirstGridlineNode, yFirstGridlineNode};
            meshkernel::Point const secondNode{xSecondGridLineNode, ySecondGridLineNode};
            meshkernel::Point const lowerLeft{xLowerLeftCornerSmoothingArea, yLowerLeftCornerSmoothingArea};
            meshkernel::Point const upperRight{xUpperRightCornerSmootingArea, yUpperRightCornerSmootingArea};

            // Execute
            meshkernel::CurvilinearGridSmoothing curvilinearGridSmoothing(*meshKernelState[meshKernelId].m_curvilinearGrid, smoothingIterations);

            curvilinearGridSmoothing.SetLine(firstNode, secondNode);
            curvilinearGridSmoothing.SetBlock(lowerLeft, upperRight);

            *meshKernelState[meshKernelId].m_curvilinearGrid = curvilinearGridSmoothing.ComputeDirectional();
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_initialize_line_shift(int meshKernelId)
    {

        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel state does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGrid == nullptr)
            {
                throw meshkernel::MeshKernelError("Not a valid curvilinear grid instance.");
            }

            if (!meshKernelState[meshKernelId].m_curvilinearGrid->IsValid())
            {
                throw meshkernel::MeshKernelError("Not valid curvilinear grid.");
            }

            meshKernelState[meshKernelId].m_curvilinearGridLineShift = std::make_shared<meshkernel::CurvilinearGridLineShift>(*meshKernelState[meshKernelId].m_curvilinearGrid);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_set_line_line_shift(int meshKernelId,
                                                            double xFirstGridLineNode,
                                                            double yFirstGridLineNode,
                                                            double xSecondGridLineNode,
                                                            double ySecondGridLineNode)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel state does not exist.");
            }

            meshkernel::Point const firstNode{xFirstGridLineNode, yFirstGridLineNode};
            meshkernel::Point const secondNode{xSecondGridLineNode, ySecondGridLineNode};

            meshKernelState[meshKernelId].m_curvilinearGridLineShift->SetLine(firstNode, secondNode);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_set(int meshKernelId, const CurvilinearGrid& grid)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel state does not exist.");
            }

            lin_alg::Matrix<meshkernel::Point> curviGridPoints(grid.num_m, grid.num_n);
            int nodeIndex = 0;
            for (int i = 0; i < grid.num_m; ++i)
            {
                for (int j = 0; j < grid.num_n; ++j)
                {

                    curviGridPoints(i, j) = meshkernel::Point(grid.node_x[nodeIndex], grid.node_y[nodeIndex]);
                    nodeIndex++;
                }
            }

            const auto& projection = meshKernelState[meshKernelId].m_projection;
            meshKernelState[meshKernelId].m_curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(curviGridPoints, projection);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_set_block_line_shift(int meshKernelId,
                                                             double xLowerLeftCorner,
                                                             double yLowerLeftCorner,
                                                             double xUpperRightCorner,
                                                             double yUpperRightCorner)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel state does not exist.");
            }

            meshkernel::Point const lowerLeftPoint{xLowerLeftCorner, yLowerLeftCorner};
            meshkernel::Point const upperRightPoint{xUpperRightCorner, yUpperRightCorner};

            meshKernelState[meshKernelId].m_curvilinearGridLineShift->SetBlock(lowerLeftPoint, upperRightPoint);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_move_node_line_shift(int meshKernelId,
                                                             double xFromCoordinate,
                                                             double yFromCoordinate,
                                                             double xToCoordinate,
                                                             double yToCoordinate)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel state does not exist.");
            }
            meshkernel::Point const fromPoint{xFromCoordinate, yFromCoordinate};
            meshkernel::Point const toPoint{xToCoordinate, yToCoordinate};
            meshKernelState[meshKernelId].m_curvilinearGridLineShift->MoveNode(fromPoint, toPoint);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_line_shift(int meshKernelId)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel state does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGridLineShift == nullptr)
            {
                throw meshkernel::MeshKernelError("Curvilinear grid line shift algorithm instance is null.");
            }

            meshKernelState[meshKernelId].m_curvilinearGridLineShift->Compute();
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_finalize_line_shift(int meshKernelId)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel state does not exist.");
            }

            meshKernelState[meshKernelId].m_curvilinearGridLineShift.reset();
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_insert_face(int meshKernelId, double xCoordinate, double yCoordinate)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel state does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGrid == nullptr)
            {
                throw meshkernel::MeshKernelError("Empty curvilinear grid");
            }

            if (!meshKernelState[meshKernelId].m_curvilinearGrid->IsValid())
            {
                throw meshkernel::MeshKernelError("Not valid curvilinear grid.");
            }

            meshkernel::Point const point{xCoordinate, yCoordinate};

            meshKernelState[meshKernelId].m_curvilinearGrid->InsertFace(point);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_convert_to_mesh2d(int meshKernelId)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() > 0 && meshKernelState[meshKernelId].m_curvilinearGrid->m_projection != meshKernelState[meshKernelId].m_mesh2d->m_projection)
            {
                throw meshkernel::MeshKernelError("The existing mesh2d projection is not equal to the curvilinear grid projection");
            }

            const auto [nodes, edges, gridIndices] = meshKernelState[meshKernelId].m_curvilinearGrid->ConvertCurvilinearToNodesAndEdges();

            const auto curviTemp = meshKernelState[meshKernelId].m_curvilinearGrid;

            *meshKernelState[meshKernelId].m_mesh2d += meshkernel::Mesh2D(edges, nodes, meshKernelState[meshKernelId].m_curvilinearGrid->m_projection);

            // curvilinear grid must be re-setted
            *meshKernelState[meshKernelId].m_curvilinearGrid = meshkernel::CurvilinearGrid();
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_delete_exterior(int meshKernelId,
                                                        double xFirstPointCoordinate,
                                                        double yFirstPointCoordinate,
                                                        double xSecondPointCoordinate,
                                                        double ySecondPointCoordinate)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGrid == nullptr)
            {
                throw meshkernel::MeshKernelError("Not a valid curvilinear grid instance.");
            }

            if (!meshKernelState[meshKernelId].m_curvilinearGrid->IsValid())
            {
                throw meshkernel::MeshKernelError("Not valid curvilinear grid.");
            }
            meshkernel::CurvilinearGridDeleteExterior curvilinearDeleteExterior(*meshKernelState[meshKernelId].m_curvilinearGrid);

            curvilinearDeleteExterior.SetBlock({xFirstPointCoordinate, yFirstPointCoordinate},
                                               {xSecondPointCoordinate, ySecondPointCoordinate});

            curvilinearDeleteExterior.Compute();
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_delete_interior(int meshKernelId,
                                                        double xFirstPointCoordinate,
                                                        double yFirstPointCoordinate,
                                                        double xSecondPointCoordinate,
                                                        double ySecondPointCoordinate)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGrid == nullptr)
            {
                throw meshkernel::MeshKernelError("Not a valid curvilinear grid instance.");
            }

            if (!meshKernelState[meshKernelId].m_curvilinearGrid->IsValid())
            {
                throw meshkernel::MeshKernelError("Not valid curvilinear grid.");
            }
            meshkernel::CurvilinearGridDeleteInterior curvilinearDeleteInterior(*meshKernelState[meshKernelId].m_curvilinearGrid);

            curvilinearDeleteInterior.SetBlock({xFirstPointCoordinate, yFirstPointCoordinate},
                                               {xSecondPointCoordinate, ySecondPointCoordinate});

            curvilinearDeleteInterior.Compute();
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_line_attraction_repulsion(int meshKernelId,
                                                                  double repulsionParameter,
                                                                  double xFirstNodeOnTheLine,
                                                                  double yFirstNodeOnTheLine,
                                                                  double xSecondNodeOnTheLine,
                                                                  double ySecondNodeOnTheLine,
                                                                  double xLowerLeftCorner,
                                                                  double yLowerLeftCorner,
                                                                  double xUpperRightCorner,
                                                                  double yUpperRightCorner)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel state does not exist.");
            }

            meshkernel::CurvilinearGridLineAttractionRepulsion curvilinearLineAttractionRepulsion(*meshKernelState[meshKernelId].m_curvilinearGrid, repulsionParameter);

            meshkernel::Point const lineFrom{xFirstNodeOnTheLine, yFirstNodeOnTheLine};
            meshkernel::Point const lineTo{xSecondNodeOnTheLine, ySecondNodeOnTheLine};
            curvilinearLineAttractionRepulsion.SetLine(lineFrom, lineTo);

            meshkernel::Point const lowerLeft{xLowerLeftCorner, yLowerLeftCorner};
            meshkernel::Point const upperRight{xUpperRightCorner, yUpperRightCorner};
            curvilinearLineAttractionRepulsion.SetBlock(lowerLeft, upperRight);

            curvilinearLineAttractionRepulsion.Compute();
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_line_mirror(int meshKernelId,
                                                    double mirroringFactor,
                                                    double xFirstGridLineNode,
                                                    double yFirstGridLineNode,
                                                    double xSecondGridLineNode,
                                                    double ySecondGridLineNode)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGrid == nullptr)
            {
                throw meshkernel::MeshKernelError("Not a valid curvilinear grid instance.");
            }

            if (!meshKernelState[meshKernelId].m_curvilinearGrid->IsValid())
            {
                throw meshkernel::MeshKernelError("Not valid curvilinear grid.");
            }

            auto curvilinearGridLineMirror = meshkernel::CurvilinearGridLineMirror(*meshKernelState[meshKernelId].m_curvilinearGrid, mirroringFactor);

            curvilinearGridLineMirror.SetLine({xFirstGridLineNode, yFirstGridLineNode}, {xSecondGridLineNode, ySecondGridLineNode});

            curvilinearGridLineMirror.Compute();
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_delete_node(int meshKernelId,
                                                    double xPointCoordinate,
                                                    double yPointCoordinate)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGrid == nullptr)
            {
                throw meshkernel::MeshKernelError("Not a valid curvilinear grid instance.");
            }

            if (!meshKernelState[meshKernelId].m_curvilinearGrid->IsValid())
            {
                throw meshkernel::MeshKernelError("Not valid curvilinear grid.");
            }

            meshKernelState[meshKernelId].m_curvilinearGrid->DeleteNode({xPointCoordinate, yPointCoordinate});
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_move_node(int meshKernelId,
                                                  double xFromPoint,
                                                  double yFromPoint,
                                                  double xToPoint,
                                                  double yToPoint)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGrid == nullptr)
            {
                throw meshkernel::MeshKernelError("Not a valid curvilinear grid instance.");
            }

            meshkernel::Point const fromPoint{xFromPoint, yFromPoint};
            meshkernel::Point const toPoint{xToPoint, yToPoint};

            meshKernelState[meshKernelId].m_curvilinearGrid->MoveNode(fromPoint, toPoint);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API double mkernel_get_separator()
    {
        return meshkernel::constants::missing::doubleValue;
    }

    MKERNEL_API double mkernel_get_inner_outer_separator()
    {
        return meshkernel::constants::missing::innerOuterSeparator;
    }

    MKERNEL_API int mkernel_mesh2d_averaging_interpolation(int meshKernelId,
                                                           const GeometryList& samples,
                                                           int locationType,
                                                           int averagingMethodType,
                                                           double relativeSearchSize,
                                                           size_t minNumSamples,
                                                           GeometryList& results)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() == 0)
            {
                throw meshkernel::MeshKernelError("The mesh is empty.");
            }

            auto sampleValues = ConvertGeometryListToSampleVector(samples);
            auto const meshLocation = static_cast<meshkernel::Mesh::Location>(locationType);
            auto const averagingMethod = static_cast<meshkernel::AveragingInterpolation::Method>(averagingMethodType);

            meshkernel::AveragingInterpolation averaging(*meshKernelState[meshKernelId].m_mesh2d,
                                                         sampleValues,
                                                         averagingMethod,
                                                         meshLocation,
                                                         relativeSearchSize,
                                                         false,
                                                         false,
                                                         static_cast<meshkernel::UInt>(minNumSamples));

            averaging.Compute();

            // Get the results
            std::vector<double> interpolationResults;
            if (meshLocation == meshkernel::Mesh::Location::Nodes)
            {
                interpolationResults = averaging.GetNodeResults();
            }
            else if (meshLocation == meshkernel::Mesh::Location::Edges)
            {
                interpolationResults = averaging.GetEdgeResults();
            }
            else if (meshLocation == meshkernel::Mesh::Location::Faces)
            {
                interpolationResults = averaging.GetFaceResults();
            }

            auto const locations = meshKernelState[meshKernelId].m_mesh2d->ComputeLocations(meshLocation);
            ConvertSampleVectorToGeometryList(locations, interpolationResults, results);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_triangulation_interpolation(int meshKernelId,
                                                               const GeometryList& samples,
                                                               int locationType,
                                                               GeometryList& results)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() == 0)
            {
                throw meshkernel::MeshKernelError("The mesh is empty.");
            }

            // Locations
            auto const sampleValues = ConvertGeometryListToSampleVector(samples);
            auto const meshLocation = static_cast<meshkernel::Mesh::Location>(locationType);
            auto const locations = meshKernelState[meshKernelId].m_mesh2d->ComputeLocations(meshLocation);

            // Execute triangulation
            meshkernel::TriangulationInterpolation triangulationInterpolation(locations, sampleValues, meshKernelState[meshKernelId].m_mesh2d->m_projection);
            triangulationInterpolation.Compute();

            // Get the results and copy them back to the results vector
            auto const& interpolationResults = triangulationInterpolation.GetResults();
            ConvertSampleVectorToGeometryList(locations, interpolationResults, results);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_get_edges_location_type(int& type)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        type = static_cast<int>(meshkernel::Mesh::Location::Edges);
        return lastExitCode;
    }
    MKERNEL_API int mkernel_get_nodes_location_type(int& type)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        type = static_cast<int>(meshkernel::Mesh::Location::Nodes);
        return lastExitCode;
    }
    MKERNEL_API int mkernel_get_faces_location_type(int& type)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        type = static_cast<int>(meshkernel::Mesh::Location::Faces);
        return lastExitCode;
    }

    MKERNEL_API int mkernel_get_averaging_method_simple_averaging(int& method)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        method = static_cast<int>(meshkernel::AveragingInterpolation::Method::SimpleAveraging);
        return lastExitCode;
    }

    MKERNEL_API int mkernel_get_averaging_method_closest_point(int& method)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        method = static_cast<int>(meshkernel::AveragingInterpolation::Method::Closest);
        return lastExitCode;
    }

    MKERNEL_API int mkernel_get_averaging_method_max(int& method)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        method = static_cast<int>(meshkernel::AveragingInterpolation::Method::Max);
        return lastExitCode;
    }

    MKERNEL_API int mkernel_get_averaging_method_min(int& method)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        method = static_cast<int>(meshkernel::AveragingInterpolation::Method::Min);
        return lastExitCode;
    }

    MKERNEL_API int mkernel_get_averaging_method_inverse_distance_weighting(int& method)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        method = static_cast<int>(meshkernel::AveragingInterpolation::Method::InverseWeightedDistance);
        return lastExitCode;
    }

    MKERNEL_API int mkernel_get_averaging_method_min_absolute_value(int& method)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        method = static_cast<int>(meshkernel::AveragingInterpolation::Method::MinAbsValue);
        return lastExitCode;
    }

    MKERNEL_API int mkernel_get_projection_cartesian(int& projection)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        projection = static_cast<int>(meshkernel::Projection::cartesian);
        return lastExitCode;
    }

    MKERNEL_API int mkernel_get_projection_spherical(int& projection)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        projection = static_cast<int>(meshkernel::Projection::spherical);
        return lastExitCode;
    }

    MKERNEL_API int mkernel_get_projection_spherical_accurate(int& projection)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        projection = static_cast<int>(meshkernel::Projection::sphericalAccurate);
        return lastExitCode;
    }

    MKERNEL_API int mkernel_get_projection(int meshKernelId, int& projection)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        projection = static_cast<int>(meshKernelState[meshKernelId].m_projection);
        return lastExitCode;
    }

} // namespace meshkernelapi
