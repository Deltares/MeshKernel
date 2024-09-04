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
#include "MeshKernel/Mesh2DToCurvilinear.hpp"
#include "MeshKernel/SamplesHessianCalculator.hpp"

#include "MeshKernel/CurvilinearGrid/CurvilinearGridDeleteInterior.hpp"
#include "MeshKernelApi/BoundingBox.hpp"

#include <MeshKernel/AveragingInterpolation.hpp>
#include <MeshKernel/BilinearInterpolationOnGriddedSamples.hpp>
#include <MeshKernel/CasulliDeRefinement.hpp>
#include <MeshKernel/CasulliRefinement.hpp>
#include <MeshKernel/ConnectMeshes.hpp>
#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Contacts.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridAlgorithm.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridCurvature.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridDeRefinement.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridFromPolygon.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridFromSplines.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridFromSplinesTransfinite.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridFullRefinement.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridGenerateCircularGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLineAttractionRepulsion.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLineMirror.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridOrthogonalization.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridRefinement.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridSmoothing.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridSmoothness.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridSnapGridToLandBoundary.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridSnapGridToSpline.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridSplineToGrid.hpp>
#include <MeshKernel/Definitions.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/FlipEdges.hpp>
#include <MeshKernel/LandBoundaries.hpp>
#include <MeshKernel/LandBoundary.hpp>
#include <MeshKernel/Mesh.hpp>
#include <MeshKernel/Mesh1D.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Mesh2DGenerateGlobal.hpp>
#include <MeshKernel/MeshConversion.hpp>
#include <MeshKernel/MeshRefinement.hpp>
#include <MeshKernel/MeshTransformation.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/OrthogonalizationAndSmoothing.hpp>
#include <MeshKernel/Orthogonalizer.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/ProjectionConversions.hpp>
#include <MeshKernel/RangeCheck.hpp>
#include <MeshKernel/RemoveDisconnectedRegions.hpp>
#include <MeshKernel/Smoother.hpp>
#include <MeshKernel/SplineAlgorithms.hpp>
#include <MeshKernel/Splines.hpp>
#include <MeshKernel/SplitRowColumnOfMesh.hpp>
#include <MeshKernel/TriangulationInterpolation.hpp>
#include <MeshKernel/UndoActions/CompoundUndoAction.hpp>
#include <MeshKernel/UndoActions/NoActionUndo.hpp>
#include <MeshKernel/UndoActions/UndoAction.hpp>
#include <MeshKernel/UndoActions/UndoActionStack.hpp>
#include <MeshKernel/Utilities/LinearAlgebra.hpp>

#include <MeshKernelApi/MKStateUndoAction.hpp>
#include <MeshKernelApi/MeshKernel.hpp>
#include <MeshKernelApi/State.hpp>
#include <MeshKernelApi/Utils.hpp>

#include <Version/Version.hpp>

#include <cstring>
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
    static meshkernel::Location invalidMeshLocation{meshkernel::Location::Unknown};

    /// @brief Stack of undo actions
    static meshkernel::UndoActionStack meshKernelUndoStack;

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
        try
        {
            meshKernelId = meshKernelStateCounter++;
            meshkernel::range_check::CheckOneOf<int>(projectionType, meshkernel::GetValidProjections(), "Projection");
            auto const projection = static_cast<meshkernel::Projection>(projectionType);
            meshKernelState.insert({meshKernelId, MeshKernelState(projection)});
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_is_valid_state(int meshKernelId, bool& isValid)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        isValid = false;

        try
        {
            isValid = meshKernelState.contains(meshKernelId);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
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

            auto undoAction = MKStateUndoAction::Create(meshKernelState[meshKernelId]);

            MeshKernelState& mkState = meshKernelState[meshKernelId];

            mkState.m_mesh1d = std::make_shared<meshkernel::Mesh1D>(mkState.m_projection);
            mkState.m_mesh2d = std::make_shared<meshkernel::Mesh2D>(mkState.m_projection);
            mkState.m_network1d = std::make_shared<meshkernel::Network1D>(mkState.m_projection);
            mkState.m_contacts = std::make_shared<meshkernel::Contacts>(*mkState.m_mesh1d, *mkState.m_mesh2d);
            mkState.m_curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(mkState.m_projection);

            mkState.m_meshOrthogonalization.reset();
            mkState.m_curvilinearGridFromSplines.reset();
            mkState.m_curvilinearGridOrthogonalization.reset();
            mkState.m_curvilinearGridLineShift.reset();

            meshKernelUndoStack.Add(std::move(undoAction), meshKernelId);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_expunge_state(int meshKernelId)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            meshKernelState.erase(meshKernelId);
            meshKernelUndoStack.Remove(meshKernelId);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_set_undo_size(int undoStackSize)
    {
        lastExitCode = meshkernel::ExitCode::Success;

        try
        {
            if (undoStackSize < 0)
            {
                throw meshkernel::MeshKernelError("Incorrect undo stack size: {}", undoStackSize);
            }

            meshKernelUndoStack.SetMaximumSize(static_cast<meshkernel::UInt>(undoStackSize));
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_undo_state(bool& undone, int& meshKernelId)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        undone = false;
        meshKernelId = meshkernel::constants::missing::intValue;

        try
        {
            if (auto undoOption = meshKernelUndoStack.Undo())
            {
                undone = true;
                meshKernelId = *undoOption;
            }
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_undo_state_count(int& committedCount, int& restoredCount)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        committedCount = 0;
        restoredCount = 0;

        try
        {
            committedCount = static_cast<int>(meshKernelUndoStack.CommittedSize());
            restoredCount = static_cast<int>(meshKernelUndoStack.RestoredSize());
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_undo_state_count_for_id(int meshKernelId, int& committedCount, int& restoredCount)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        committedCount = 0;
        restoredCount = 0;

        try
        {
            if (meshKernelState.contains(meshKernelId))
            {
                committedCount = static_cast<int>(meshKernelUndoStack.CommittedSize(meshKernelId));
                restoredCount = static_cast<int>(meshKernelUndoStack.RestoredSize(meshKernelId));
            }
            else
            {
                committedCount = 0;
                restoredCount = 0;
            }
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_redo_state(bool& redone, int& meshKernelId)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        redone = false;
        meshKernelId = meshkernel::constants::missing::intValue;

        try
        {
            if (auto redoOption = meshKernelUndoStack.Commit())
            {
                redone = true;
                meshKernelId = *redoOption;
            }
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_clear_state()
    {
        lastExitCode = meshkernel::ExitCode::Success;

        try
        {
            meshKernelUndoStack.Clear();
            meshKernelState.clear();
            meshKernelStateCounter = 0;
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_clear_undo_state()
    {
        lastExitCode = meshkernel::ExitCode::Success;

        try
        {
            meshKernelUndoStack.Clear();
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_clear_undo_state_for_id(int meshKernelId)
    {
        lastExitCode = meshkernel::ExitCode::Success;

        try
        {
            meshKernelUndoStack.Remove(meshKernelId);
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

            meshkernel::range_check::CheckOneOf<int>(deletionOption, meshkernel::GetValidDeletionOptions(), "Deletion");

            const std::vector<meshkernel::Point> polygonPoints = ConvertGeometryListToPointVector(polygon);

            const bool invertDeletionBool = invertDeletion == 1;
            const meshkernel::Polygons meshKernelPolygon(polygonPoints, meshKernelState[meshKernelId].m_mesh2d->m_projection);
            const auto deletionOptionEnum = static_cast<meshkernel::Mesh2D::DeleteMeshOptions>(deletionOption);

            meshKernelUndoStack.Add(meshKernelState[meshKernelId].m_mesh2d->DeleteMesh(meshKernelPolygon, deletionOptionEnum, invertDeletionBool), meshKernelId);
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

            auto undoAction = MKStateUndoAction::Create(meshKernelState[meshKernelId]);

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
                meshKernelState[meshKernelId].m_mesh2d = std::make_unique<meshkernel::Mesh2D>(edges2d,
                                                                                              nodes2d,
                                                                                              face_nodes,
                                                                                              num_face_nodes,
                                                                                              meshKernelState[meshKernelId].m_projection);
            }
            else
            {

                // Do not change the pointer, just the object it is pointing to
                // Compute the faces
                meshKernelState[meshKernelId].m_mesh2d = std::make_unique<meshkernel::Mesh2D>(edges2d,
                                                                                              nodes2d,
                                                                                              meshKernelState[meshKernelId].m_projection);
            }

            meshKernelUndoStack.Add(std::move(undoAction), meshKernelId);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_snap_to_landboundary(int meshKernelId, const GeometryList& selectingPolygon, const GeometryList& landBoundaries)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_mesh2d == nullptr)
            {
                throw meshkernel::MeshKernelError("The selected mesh not exist.");
            }

            auto const landBoundariesPoints = ConvertGeometryListToPointVector(landBoundaries);
            auto const polygonNodes = ConvertGeometryListToPointVector(selectingPolygon);

            // Construct all dependencies
            const auto polygon = meshkernel::Polygons(polygonNodes, meshKernelState[meshKernelId].m_projection);
            auto landBoundary = meshkernel::LandBoundaries(landBoundariesPoints, *meshKernelState[meshKernelId].m_mesh2d, polygon);
            landBoundary.FindNearestMeshBoundary(meshkernel::LandBoundaries::ProjectToLandBoundaryOption::InnerAndOuterMeshBoundaryToLandBoundary);

            // Execute algorithm
            meshKernelUndoStack.Add(landBoundary.SnapMeshToLandBoundaries());
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_split_row(int meshKernelId,
                                             int firstNode,
                                             int secondNode)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_mesh2d == nullptr)
            {
                throw meshkernel::MeshKernelError("The selected mesh not exist.");
            }

            meshkernel::UInt edgeId = meshKernelState[meshKernelId].m_mesh2d->FindEdge(static_cast<meshkernel::UInt>(firstNode),
                                                                                       static_cast<meshkernel::UInt>(secondNode));

            meshkernel::SplitRowColumnOfMesh splitAlongRow;
            meshKernelUndoStack.Add(splitAlongRow.Compute(*meshKernelState[meshKernelId].m_mesh2d, edgeId), meshKernelId);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_add(int meshKernelId, const Mesh2D& mesh2d)
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

            std::unique_ptr<meshkernel::UndoAction> undoAction;

            if (mesh2d.num_faces > 0 && mesh2d.face_nodes != nullptr && mesh2d.nodes_per_face != nullptr)
            {
                const auto face_nodes = meshkernel::ConvertToFaceNodesVector(mesh2d.num_faces, mesh2d.face_nodes, mesh2d.nodes_per_face);

                std::vector<meshkernel::UInt> num_face_nodes;
                num_face_nodes.reserve(mesh2d.num_faces);
                for (auto n = 0; n < mesh2d.num_faces; n++)
                {
                    num_face_nodes.emplace_back(static_cast<meshkernel::UInt>(mesh2d.nodes_per_face[n]));
                }

                undoAction = meshKernelState[meshKernelId].m_mesh2d->Join(meshkernel::Mesh2D(edges2d, nodes2d, face_nodes, num_face_nodes,
                                                                                             meshKernelState[meshKernelId].m_projection));
            }
            else
            {
                // Compute the faces
                undoAction = meshKernelState[meshKernelId].m_mesh2d->Join(meshkernel::Mesh2D(edges2d, nodes2d, meshKernelState[meshKernelId].m_projection));
            }

            meshKernelUndoStack.Add(std::move(undoAction), meshKernelId);
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

            auto undoAction = MKStateUndoAction::Create(meshKernelState[meshKernelId]);

            // Do not change the pointer, just the object it is pointing to
            meshKernelState[meshKernelId].m_mesh1d = std::make_unique<meshkernel::Mesh1D>(edges1d, nodes1d, meshKernelState[meshKernelId].m_projection);

            meshKernelUndoStack.Add(std::move(undoAction), meshKernelId);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh1d_add(int meshKernelId,
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

            meshKernelUndoStack.Add(meshKernelState[meshKernelId].m_mesh1d->Join(meshkernel::Mesh1D(edges1d, nodes1d, meshKernelState[meshKernelId].m_projection)),
                                    meshKernelId);
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

            auto undoAction = MKStateUndoAction::Create(meshKernelState[meshKernelId]);

            // Do not change the pointer, just the object it is pointing to
            *meshKernelState[meshKernelId].m_network1d = meshkernel::Network1D(localPolylines, meshKernelState[meshKernelId].m_projection);

            meshKernelUndoStack.Add(std::move(undoAction), meshKernelId);
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

            meshKernelUndoStack.Add(meshKernelState[meshKernelId].m_mesh1d->Join(meshkernel::Mesh1D(*meshKernelState[meshKernelId].m_network1d, minFaceSize)),
                                    meshKernelId);
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
            SetMesh2dApiDimensions(*meshKernelState[meshKernelId].m_mesh2d, mesh2d);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_node_edge_data(int meshKernelId, Mesh2D& mesh2d)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (mesh2d.node_x == nullptr || mesh2d.node_y == nullptr || mesh2d.edge_nodes == nullptr)
            {
                throw meshkernel::MeshKernelError("The meshkernel api Mesh2D has not been initialised correctly.");
            }

            SetMesh2dApiNodeEdgeData(*meshKernelState[meshKernelId].m_mesh2d, mesh2d);
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

            SetMesh2dApiData(*meshKernelState[meshKernelId].m_mesh2d, mesh2d);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_orthogonality_property_type(int& type)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        type = static_cast<int>(meshkernel::Mesh2D::Property::Orthogonality);
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

            SetMesh1dApiDimension(*meshKernelState[meshKernelId].m_mesh1d, mesh1d);
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

            SetMesh1dApiData(*meshKernelState[meshKernelId].m_mesh1d, mesh1d);
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

            auto mk = meshKernelState[meshKernelId];

            if (meshKernelState[meshKernelId].m_curvilinearGrid->IsValid())
            {
                curvilinearGrid.num_n = static_cast<int>(meshKernelState[meshKernelId].m_curvilinearGrid->NumN());
                curvilinearGrid.num_m = static_cast<int>(meshKernelState[meshKernelId].m_curvilinearGrid->NumM());
            }
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_get_location_index(int meshKernelId,
                                                           double xCoordinate,
                                                           double yCoordinate,
                                                           int locationType,
                                                           const BoundingBox& boundingBox,
                                                           int& locationIndex)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGrid->GetNumNodes() <= 0)
            {
                throw meshkernel::ConstraintError("The selected curvilinear grid has no nodes.");
            }

            auto const meshLocation = static_cast<meshkernel::Location>(locationType);

            meshkernel::Point const point{xCoordinate, yCoordinate};
            meshkernel::BoundingBox box{{boundingBox.xLowerLeft, boundingBox.yLowerLeft}, {boundingBox.xUpperRight, boundingBox.yUpperRight}};

            locationIndex = static_cast<int>(meshKernelState[meshKernelId].m_curvilinearGrid->FindLocationIndex(point, meshLocation, box));
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
            SetCurvilinearGridApiData(*meshKernelState[meshKernelId].m_curvilinearGrid, curvilinearGrid);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_get_boundaries_as_polygons(int meshKernelId, int lowerLeftN, int lowerLeftM, int upperRightN, int upperRightM, GeometryList& boundaryPolygons)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (!meshKernelState[meshKernelId].m_curvilinearGrid->IsValid())
            {
                throw meshkernel::MeshKernelError("Invalid curvilinear grid");
            }

            auto lowerLeftNUnsigned = static_cast<meshkernel::UInt>(lowerLeftN);
            auto lowerLeftMUnsigned = static_cast<meshkernel::UInt>(lowerLeftM);
            auto upperRightNUnsigned = static_cast<meshkernel::UInt>(upperRightN);
            auto upperRightMUnsigned = static_cast<meshkernel::UInt>(upperRightM);

            const auto minN = std::min(lowerLeftNUnsigned, upperRightNUnsigned);
            const auto maxN = std::max(lowerLeftNUnsigned, upperRightNUnsigned);
            const auto minM = std::min(lowerLeftMUnsigned, upperRightMUnsigned);
            const auto maxM = std::max(lowerLeftMUnsigned, upperRightMUnsigned);

            const auto boundaryPolygon = meshKernelState[meshKernelId].m_curvilinearGrid->ComputeBoundaryPolygons({minN, minM},
                                                                                                                  {maxN, maxM});
            ConvertPointVectorToGeometryList(boundaryPolygon, boundaryPolygons);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_count_boundaries_as_polygons(int meshKernelId, int lowerLeftN, int lowerLeftM, int upperRightN, int upperRightM, int& numberOfPolygonNodes)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (!meshKernelState[meshKernelId].m_curvilinearGrid->IsValid())
            {
                throw meshkernel::MeshKernelError("Invalid curvilinear grid");
            }

            const auto lowerLeftNUnsigned = static_cast<meshkernel::UInt>(lowerLeftN);
            const auto lowerLeftMUnsigned = static_cast<meshkernel::UInt>(lowerLeftM);
            const auto upperRightNUnsigned = static_cast<meshkernel::UInt>(upperRightN);
            const auto upperRightMUnsigned = static_cast<meshkernel::UInt>(upperRightM);

            const auto minN = std::min(lowerLeftNUnsigned, upperRightNUnsigned);
            const auto maxN = std::max(lowerLeftNUnsigned, upperRightNUnsigned);
            const auto minM = std::min(lowerLeftMUnsigned, upperRightMUnsigned);
            const auto maxM = std::max(lowerLeftMUnsigned, upperRightMUnsigned);

            const auto boundaryPolygon = meshKernelState[meshKernelId].m_curvilinearGrid->ComputeBoundaryPolygons({minN, minM},
                                                                                                                  {maxN, maxM});
            numberOfPolygonNodes = static_cast<int>(boundaryPolygon.size());
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

    MKERNEL_API int mkernel_contacts_set(int meshKernelId, const Contacts& contacts)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            std::vector<meshkernel::UInt> mesh1dIndices(contacts.num_contacts);
            std::vector<meshkernel::UInt> mesh2dIndices(contacts.num_contacts);

            for (int i = 0; i < contacts.num_contacts; ++i)
            {
                mesh1dIndices[i] = contacts.mesh1d_indices[i];
                mesh2dIndices[i] = contacts.mesh2d_indices[i];
            }

            meshKernelState[meshKernelId].m_contacts->SetIndices(mesh1dIndices, mesh2dIndices);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_convert_projection(int meshKernelId, int projectionType, const char* const zoneString)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            const meshkernel::Projection targetProjection = meshkernel::GetProjectionValue(projectionType);

            const meshkernel::Projection& sourceProjection = meshKernelState[meshKernelId].m_mesh2d->m_projection;

            if (sourceProjection != targetProjection)
            {
                if (sourceProjection == meshkernel::Projection::cartesian)
                {
                    meshkernel::ConvertCartesianToSpherical conversion(zoneString);
                    meshKernelUndoStack.Add(meshkernel::MeshConversion::Compute(*meshKernelState[meshKernelId].m_mesh2d, conversion), meshKernelId);
                    meshKernelState[meshKernelId].m_projection = conversion.TargetProjection();
                }
                else if (sourceProjection == meshkernel::Projection::spherical)
                {
                    meshkernel::ConvertSphericalToCartesian conversion(zoneString);
                    meshKernelUndoStack.Add(meshkernel::MeshConversion::Compute(*meshKernelState[meshKernelId].m_mesh2d, conversion), meshKernelId);
                    meshKernelState[meshKernelId].m_projection = conversion.TargetProjection();
                }
                else
                {
                    throw meshkernel::MeshKernelError("Mesh conversion between projection {} and {} has not been implemented.",
                                                      meshkernel::ProjectionToString(sourceProjection),
                                                      meshkernel::ProjectionToString(targetProjection));
                }
            }
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_convert_to_curvilinear(int meshKernelId, double xPointCoordinate, double yPointCoordinate)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            auto undoAction = MKStateUndoAction::Create(meshKernelState[meshKernelId]);
            meshkernel::Mesh2DToCurvilinear mesh2DToCurvilinear(*meshKernelState[meshKernelId].m_mesh2d);

            meshKernelState[meshKernelId].m_curvilinearGrid = mesh2DToCurvilinear.Compute({xPointCoordinate, yPointCoordinate});
            meshKernelState[meshKernelId].m_mesh2d = std::make_unique<meshkernel::Mesh2D>(meshKernelState[meshKernelId].m_projection);
            meshKernelUndoStack.Add(std::move(undoAction), meshKernelId);
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

            meshKernelUndoStack.Add(meshKernelState[meshKernelId].m_mesh2d->DeleteHangingEdges(), meshKernelId);
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
            auto smoother = std::make_unique<meshkernel::Smoother>(*meshKernelState[meshKernelId].m_mesh2d);
            auto orthogonalizer = std::make_unique<meshkernel::Orthogonalizer>(*meshKernelState[meshKernelId].m_mesh2d);
            auto polygon = std::make_unique<meshkernel::Polygons>(polygonNodes, meshKernelState[meshKernelId].m_mesh2d->m_projection);
            auto landBoundary = std::make_unique<meshkernel::LandBoundaries>(landBoundariesPoints, *meshKernelState[meshKernelId].m_mesh2d, *polygon);

            meshkernel::OrthogonalizationAndSmoothing ortogonalization(*meshKernelState[meshKernelId].m_mesh2d,
                                                                       std::move(smoother),
                                                                       std::move(orthogonalizer),
                                                                       std::move(polygon),
                                                                       std::move(landBoundary),
                                                                       static_cast<meshkernel::LandBoundaries::ProjectToLandBoundaryOption>(projectToLandBoundaryOption),
                                                                       orthogonalizationParameters);
            meshKernelUndoStack.Add(ortogonalization.Initialize(), meshKernelId);
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
            auto smoother = std::make_unique<meshkernel::Smoother>(*meshKernelState[meshKernelId].m_mesh2d);
            auto orthogonalizer = std::make_unique<meshkernel::Orthogonalizer>(*meshKernelState[meshKernelId].m_mesh2d);
            auto polygon = std::make_unique<meshkernel::Polygons>(polygonNodesVector, meshKernelState[meshKernelId].m_mesh2d->m_projection);
            auto landBoundary = std::make_unique<meshkernel::LandBoundaries>(landBoundariesNodeVector, *meshKernelState[meshKernelId].m_mesh2d, *polygon);

            meshKernelState[meshKernelId].m_meshOrthogonalization = std::make_unique<meshkernel::OrthogonalizationAndSmoothing>(*meshKernelState[meshKernelId].m_mesh2d,
                                                                                                                                std::move(smoother),
                                                                                                                                std::move(orthogonalizer),
                                                                                                                                std::move(polygon),
                                                                                                                                std::move(landBoundary),
                                                                                                                                static_cast<meshkernel::LandBoundaries::ProjectToLandBoundaryOption>(projectToLandBoundaryOption),
                                                                                                                                orthogonalizationParameters);
            meshKernelUndoStack.Add(meshKernelState[meshKernelId].m_meshOrthogonalization->Initialize(), meshKernelId);
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

    MKERNEL_API int mkernel_mesh2d_get_location_index(int meshKernelId,
                                                      double xCoordinate,
                                                      double yCoordinate,
                                                      int locationType,
                                                      const BoundingBox& boundingBox,
                                                      int& locationIndex)
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
                throw meshkernel::ConstraintError("The selected mesh grid has no nodes.");
            }

            auto const meshLocation = static_cast<meshkernel::Location>(locationType);

            meshkernel::Point const point{xCoordinate, yCoordinate};
            meshkernel::BoundingBox box{{boundingBox.xLowerLeft, boundingBox.yLowerLeft}, {boundingBox.xUpperRight, boundingBox.yUpperRight}};

            locationIndex = static_cast<int>(meshKernelState[meshKernelId].m_mesh2d->FindLocationIndex(point, meshLocation, {}, box));
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

    MKERNEL_API int mkernel_mesh2d_get_property(int meshKernelId, int propertyValue, const GeometryList& geometryList)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            const auto& mesh2d = meshKernelState.at(meshKernelId).m_mesh2d;

            if (!mesh2d || mesh2d->GetNumNodes() <= 0)
            {
                return lastExitCode;
            }

            const auto propertyValueEnum = static_cast<meshkernel::Mesh2D::Property>(propertyValue);
            switch (propertyValueEnum)
            {
            case meshkernel::Mesh2D::Property::Orthogonality:
            {
                std::vector<double> values = mesh2d->GetOrthogonality();
                if (static_cast<size_t>(geometryList.num_coordinates) < values.size())
                {
                    throw meshkernel::MeshKernelError("GeometryList with wrong dimensions");
                }
                std::copy(values.begin(), values.end(), geometryList.values);
            }
            break;
            case meshkernel::Mesh2D::Property::EdgeLength:
            {
                mesh2d->ComputeEdgesLengths();
                std::vector<double> values = mesh2d->m_edgeLengths;
                if (static_cast<size_t>(geometryList.num_coordinates) < values.size())
                {
                    throw meshkernel::MeshKernelError("GeometryList with wrong dimensions");
                }
                std::copy(values.begin(), values.end(), geometryList.values);
            }
            break;
            default:
                throw meshkernel::MeshKernelError("Property not supported");
            }
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_property_dimension(int meshKernelId, int propertyValue, int& dimension)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            const auto& mesh2d = meshKernelState.at(meshKernelId).m_mesh2d;

            if (!mesh2d || mesh2d->GetNumNodes() <= 0)
            {
                return lastExitCode;
            }

            const auto propertyValueEnum = static_cast<meshkernel::Mesh2D::Property>(propertyValue);
            dimension = -1;
            switch (propertyValueEnum)
            {
            case meshkernel::Mesh2D::Property::Orthogonality:
                dimension = static_cast<int>(mesh2d->GetOrthogonality().size());
                break;

            case meshkernel::Mesh2D::Property::EdgeLength:
                mesh2d->ComputeEdgesLengths();
                dimension = static_cast<int>(mesh2d->m_edgeLengths.size());
                break;
            default:
                throw meshkernel::MeshKernelError("Property not supported");
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

    MKERNEL_API int mkernel_mesh2d_make_global(int meshKernelId, int numLongitudeNodes, int numLatitudeNodes)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (numLongitudeNodes == 0)
            {
                throw meshkernel::MeshKernelError("The number of longitude nodes cannot be 0");
            }

            if (numLatitudeNodes == 0)
            {
                throw meshkernel::MeshKernelError("The number of latitude nodes cannot be 0");
            }

            const auto mesh = meshkernel::Mesh2DGenerateGlobal::Compute(numLongitudeNodes, numLatitudeNodes, meshKernelState[meshKernelId].m_projection);
            meshKernelUndoStack.Add(meshKernelState[meshKernelId].m_mesh2d->Join(*mesh), meshKernelId);
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
            meshKernelUndoStack.Add(meshKernelState[meshKernelId].m_mesh2d->Join(mesh), meshKernelId);
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
            meshKernelUndoStack.Add(meshKernelState[meshKernelId].m_mesh2d->Join(mesh), meshKernelId);
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

            auto const nodes = curvilinearGrid->ComputeNodes();
            auto const edges = curvilinearGrid->ComputeEdges();
            meshKernelUndoStack.Add(meshKernelState[meshKernelId].m_mesh2d->Join(meshkernel::Mesh2D(edges, nodes, projection)), meshKernelId);
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

            const auto edges = curvilinearGrid->ComputeEdges();
            const auto nodes = curvilinearGrid->ComputeNodes();

            meshKernelUndoStack.Add(meshKernelState[meshKernelId].m_mesh2d->Join(meshkernel::Mesh2D(edges, nodes, projection)), meshKernelId);
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

            const auto edges = curvilinearGrid->ComputeEdges();
            const auto nodes = curvilinearGrid->ComputeNodes();

            meshKernelUndoStack.Add(meshKernelState[meshKernelId].m_mesh2d->Join(meshkernel::Mesh2D(edges, nodes, meshKernelState[meshKernelId].m_curvilinearGrid->projection())),
                                    meshKernelId);
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
            const auto meshBoundaryPolygon = meshKernelState[meshKernelId].m_mesh2d->ComputeBoundaryPolygons(polygonNodes);

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
            const auto meshBoundaryPolygon = meshKernelState[meshKernelId].m_mesh2d->ComputeBoundaryPolygons(polygonNodes);
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

            const meshkernel::Polygons polygon(polygonVector, meshKernelState[meshKernelId].m_projection);
            auto const refinementResult = polygon.RefineFirstPolygon(firstNodeIndex, secondNodeIndex, targetEdgeLength);

            ConvertPointVectorToGeometryList(refinementResult, refinedPolygon);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_polygon_linear_refine(int meshKernelId, const GeometryList& polygonToRefine, int firstNodeIndex, int secondNodeIndex, GeometryList& refinedPolygon)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            auto const polygonVector = ConvertGeometryListToPointVector(polygonToRefine);

            const meshkernel::Polygons polygon(polygonVector, meshKernelState[meshKernelId].m_projection);
            const auto firstNodeIndexUnsigned = static_cast<meshkernel::UInt>(firstNodeIndex);
            const auto secondNodeUnsigned = static_cast<meshkernel::UInt>(secondNodeIndex);
            const auto refinementResult = polygon.LinearRefinePolygon(0, firstNodeIndexUnsigned, secondNodeUnsigned);

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

            const meshkernel::Polygons polygon(polygonVector, meshKernelState[meshKernelId].m_projection);

            const auto refinedPolygon = polygon.RefineFirstPolygon(firstIndex, secondIndex, distance);

            numberOfPolygonNodes = static_cast<int>(refinedPolygon.size());
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_polygon_count_linear_refine(int meshKernelId, const GeometryList& polygonToRefine, int firstNodeIndex, int secondNodeIndex, int& numberOfPolygonNodes)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            auto const polygonVector = ConvertGeometryListToPointVector(polygonToRefine);

            const meshkernel::Polygons polygon(polygonVector, meshKernelState[meshKernelId].m_projection);
            const auto firstNodeIndexUnsigned = static_cast<meshkernel::UInt>(firstNodeIndex);
            const auto secondNodeUnsigned = static_cast<meshkernel::UInt>(secondNodeIndex);
            const auto refinementResult = polygon.LinearRefinePolygon(0, firstNodeIndexUnsigned, secondNodeUnsigned);

            numberOfPolygonNodes = static_cast<int>(refinementResult.size());
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
            meshKernelUndoStack.Add(meshKernelState[meshKernelId].m_mesh2d->MergeNodesInPolygon(polygon, searchRadius), meshKernelId);
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

            meshKernelUndoStack.Add(meshKernelState[meshKernelId].m_mesh2d->MergeNodesInPolygon(polygon, mergingDistance), meshKernelId);
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

            meshKernelUndoStack.Add(meshKernelState[meshKernelId].m_mesh2d->MergeTwoNodes(firstNode, secondNode), meshKernelId);
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

            auto [edgeId, undoAction] = meshKernelState[meshKernelId].m_mesh2d->ConnectNodes(startNode, endNode);
            meshKernelUndoStack.Add(std::move(undoAction), meshKernelId);

            new_edge_index = static_cast<int>(edgeId);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_insert_edge_from_coordinates(int meshKernelId,
                                                                double firstNodeX,
                                                                double firstNodeY,
                                                                double secondNodeX,
                                                                double secondNodeY,
                                                                int& firstNodeIndex,
                                                                int& secondNodeIndex,
                                                                int& edgeIndex)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            std::unique_ptr<meshkernel::CompoundUndoAction> compoundUndoAction = meshkernel::CompoundUndoAction::Create();
            meshkernel::Point const firstNodeCoordinates{firstNodeX, firstNodeY};
            meshkernel::Point const secondNodeCoordinates{secondNodeX, secondNodeY};

            const auto newEdgeLength = ComputeDistance(firstNodeCoordinates, secondNodeCoordinates, meshKernelState[meshKernelId].m_projection);
            const auto& edgeLengths = meshKernelState[meshKernelId].m_mesh2d->m_edgeLengths;
            constexpr auto lengthFraction = 0.01;

            const auto minMeshEdgeLength = edgeLengths.empty() ? newEdgeLength : *std::ranges::min_element(edgeLengths);
            const auto searchRadius = std::min(newEdgeLength * lengthFraction, minMeshEdgeLength * lengthFraction);

            if (searchRadius <= 0.0)
            {
                throw meshkernel::MeshKernelError("The first and the second node are coinciding.");
            }

            meshKernelState[meshKernelId].m_mesh2d->BuildTree(meshkernel::Location::Nodes);

            auto firstNodeId = meshKernelState[meshKernelId].m_mesh2d->FindNodeCloseToAPoint(firstNodeCoordinates, searchRadius);
            if (firstNodeId == meshkernel::constants::missing::uintValue)
            {
                auto result = meshKernelState[meshKernelId].m_mesh2d->InsertNode(firstNodeCoordinates);
                std::tie(firstNodeId, std::ignore) = result;
                compoundUndoAction->Add(std::move(std::get<1>(result)));
                meshKernelState[meshKernelId].m_mesh2d->BuildTree(meshkernel::Location::Nodes);
            }
            firstNodeIndex = static_cast<int>(firstNodeId);

            auto secondNodeId = meshKernelState[meshKernelId].m_mesh2d->FindNodeCloseToAPoint(secondNodeCoordinates, searchRadius);
            if (secondNodeId == meshkernel::constants::missing::uintValue)
            {
                auto result = meshKernelState[meshKernelId].m_mesh2d->InsertNode(secondNodeCoordinates);
                std::tie(secondNodeId, std::ignore) = result;
                compoundUndoAction->Add(std::move(std::get<1>(result)));
            }
            secondNodeIndex = static_cast<int>(secondNodeId);

            meshKernelState[meshKernelId].m_mesh2d->BuildTree(meshkernel::Location::Edges);

            auto [edgeId, action] = meshKernelState[meshKernelId].m_mesh2d->ConnectNodes(firstNodeId, secondNodeId);
            if (edgeId != meshkernel::constants::missing::uintValue)
            {
                compoundUndoAction->Add(std::move(action));
                meshKernelUndoStack.Add(std::move(compoundUndoAction), meshKernelId);
            }
            edgeIndex = static_cast<int>(edgeId);
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

            auto [nodeId, undoAction] = meshKernelState[meshKernelId].m_mesh2d->InsertNode(nodeCoordinateVector);
            meshKernelUndoStack.Add(std::move(undoAction), meshKernelId);

            nodeIndex = static_cast<int>(nodeId);
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

            meshKernelUndoStack.Add(meshKernelState[meshKernelId].m_mesh2d->DeleteNode(nodeIndex), meshKernelId);
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

            meshKernelUndoStack.Add(meshKernelState[meshKernelId].m_mesh2d->MoveNode(newPosition, nodeIndex), meshKernelId);
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

            const auto edgeIndex = meshKernelState[meshKernelId].m_mesh2d->FindLocationIndex(point, meshkernel::Location::Edges, {}, boundingBox);
            meshKernelUndoStack.Add(meshKernelState[meshKernelId].m_mesh2d->DeleteEdge(edgeIndex), meshKernelId);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_delete_edge_by_index(int meshKernelId, int edgeIndex)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            meshKernelUndoStack.Add(meshKernelState[meshKernelId].m_mesh2d->DeleteEdge(edgeIndex), meshKernelId);
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

            edgeIndex = static_cast<int>(meshKernelState[meshKernelId].m_mesh2d->FindLocationIndex(point, meshkernel::Location::Edges, {}, boundingBox));
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_face_polygons(int meshKernelId, int numEdges, const GeometryList& facePolygons)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }
            meshKernelState[meshKernelId].m_mesh2d->Administrate();
            const auto numFaces = meshKernelState[meshKernelId].m_mesh2d->GetNumFaces();
            std::vector<bool> validFace(numFaces, false);
            for (meshkernel::UInt f = 0; f < numFaces; ++f)
            {
                const auto& faceNodes = meshKernelState[meshKernelId].m_mesh2d->m_facesNodes[f];
                const auto faceNumEdges = static_cast<int>(faceNodes.size());
                if (faceNumEdges == numEdges)
                {
                    validFace[f] = true;
                }
            }
            FillFacePolygons(meshKernelState[meshKernelId].m_mesh2d, validFace, facePolygons);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_face_polygons_dimension(int meshKernelId, int numEdges, int& geometryListDimension)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }
            meshKernelState[meshKernelId].m_mesh2d->Administrate();
            const auto numFaces = meshKernelState[meshKernelId].m_mesh2d->GetNumFaces();
            int numMatchingFaces = 0;
            for (meshkernel::UInt f = 0; f < numFaces; ++f)
            {
                const auto faceNumEdges = static_cast<int>(meshKernelState[meshKernelId].m_mesh2d->m_facesNodes[f].size());
                if (faceNumEdges != numEdges)
                {
                    continue;
                }
                numMatchingFaces += 1;
            }
            if (numMatchingFaces > 0)
            {
                geometryListDimension = (numEdges + 2) * (numMatchingFaces - 1) + numEdges + 1;
            }
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_filtered_face_polygons_dimension(int meshKernelId,
                                                                        int propertyValue,
                                                                        double minValue,
                                                                        double maxValue,
                                                                        int& geometryListDimension)
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
            const auto filterEnum = static_cast<meshkernel::Mesh2D::Property>(propertyValue);
            const auto filterMask = meshKernelState[meshKernelId].m_mesh2d->FilterBasedOnMetric(meshkernel::Location::Faces,
                                                                                                filterEnum,
                                                                                                minValue,
                                                                                                maxValue);
            geometryListDimension = 0;
            for (meshkernel::UInt f = 0; f < filterMask.size(); ++f)
            {
                if (!filterMask[f])
                {
                    continue;
                }
                const auto faceNumEdges = static_cast<int>(meshKernelState[meshKernelId].m_mesh2d->m_facesNodes[f].size());
                geometryListDimension += faceNumEdges + 2;
            }
            if (geometryListDimension > 0)
            {
                geometryListDimension -= 1;
            }
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_filtered_face_polygons(int meshKernelId,
                                                              int propertyValue,
                                                              double minValue,
                                                              double maxValue,
                                                              const GeometryList& facePolygons)
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

            const auto filterEnum = static_cast<meshkernel::Mesh2D::Property>(propertyValue);
            const auto filterMask = meshKernelState[meshKernelId].m_mesh2d->FilterBasedOnMetric(meshkernel::Location::Faces,
                                                                                                filterEnum,
                                                                                                minValue,
                                                                                                maxValue);
            FillFacePolygons(meshKernelState[meshKernelId].m_mesh2d, filterMask, facePolygons);
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

            // averagingMethod may be used uninitialised;
            meshkernel::AveragingInterpolation::Method averagingMethod;

            if (meshRefinementParameters.refinement_type == static_cast<int>(meshkernel::MeshRefinement::RefinementType::WaveCourant))
            {
                averagingMethod = meshkernel::AveragingInterpolation::Method::MinAbsValue;
            }
            else if (meshRefinementParameters.refinement_type == static_cast<int>(meshkernel::MeshRefinement::RefinementType::RefinementLevels))
            {
                averagingMethod = meshkernel::AveragingInterpolation::Method::Max;
            }
            else
            {
                throw meshkernel::MeshKernelError("Invalid mesh refinement type.");
            }

            const bool refineOutsideFace = meshRefinementParameters.account_for_samples_outside == 1 ? true : false;
            const bool transformSamples = meshRefinementParameters.refinement_type == 2 ? true : false;

            auto averaging = std::make_unique<meshkernel::AveragingInterpolation>(*meshKernelState[meshKernelId].m_mesh2d,
                                                                                  samplesVector,
                                                                                  averagingMethod,
                                                                                  meshkernel::Location::Faces,
                                                                                  relativeSearchRadius,
                                                                                  refineOutsideFace,
                                                                                  transformSamples,
                                                                                  static_cast<meshkernel::UInt>(minimumNumSamples));

            meshkernel::MeshRefinement meshRefinement(*meshKernelState[meshKernelId].m_mesh2d,
                                                      std::move(averaging),
                                                      meshRefinementParameters);
            meshKernelUndoStack.Add(meshRefinement.Compute(), meshKernelId);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_refine_ridges_based_on_gridded_samples(int meshKernelId,
                                                                          const GriddedSamples& samples,
                                                                          double relativeSearchRadius,
                                                                          int minimumNumSamples,
                                                                          int numberOfSmoothingIterations,
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
            if (meshRefinementParameters.refinement_type != static_cast<int>(meshkernel::MeshRefinement::RefinementType::RidgeDetection))
            {
                throw meshkernel::MeshKernelError("The mesh refinement type in MeshRefinementParameters must be set equal to ridge refinement.");
            }

            const auto samplesVector = ConvertGriddedData(samples);

            auto samplesHessian = meshkernel::SamplesHessianCalculator::ComputeSamplesHessian(samplesVector,
                                                                                              meshKernelState[meshKernelId].m_projection,
                                                                                              numberOfSmoothingIterations,
                                                                                              samples.num_x,
                                                                                              samples.num_y);

            auto averaging = std::make_unique<meshkernel::AveragingInterpolation>(*meshKernelState[meshKernelId].m_mesh2d,
                                                                                  samplesHessian,
                                                                                  meshkernel::AveragingInterpolation::Method::Max,
                                                                                  meshkernel::Location::Faces,
                                                                                  relativeSearchRadius,
                                                                                  false,
                                                                                  false,
                                                                                  static_cast<meshkernel::UInt>(minimumNumSamples));

            meshkernel::MeshRefinement meshRefinement(*meshKernelState[meshKernelId].m_mesh2d,
                                                      std::move(averaging),
                                                      meshRefinementParameters);
            meshKernelUndoStack.Add(meshRefinement.Compute(), meshKernelId);
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

            auto interpolant = CreateBilinearInterpolatorBasedOnType(griddedSamples, *meshKernelState[meshKernelId].m_mesh2d);

            meshkernel::MeshRefinement meshRefinement(*meshKernelState[meshKernelId].m_mesh2d,
                                                      std::move(interpolant),
                                                      meshRefinementParameters,
                                                      useNodalRefinement);
            meshKernelUndoStack.Add(meshRefinement.Compute(), meshKernelId);
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

            const auto points = ConvertGeometryListToPointVector(geometryList);

            const auto polygon = meshkernel::Polygons(points, meshKernelState[meshKernelId].m_mesh2d->m_projection);

            meshkernel::MeshRefinement meshRefinement(*meshKernelState[meshKernelId].m_mesh2d, polygon, meshRefinementParameters);
            meshKernelUndoStack.Add(meshRefinement.Compute(), meshKernelId);
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
            meshKernelUndoStack.Add(removeDisconnectedRegions.Compute(*meshKernelState[meshKernelId].m_mesh2d), meshKernelId);
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

            meshKernelState[meshKernelId].m_mesh2d->BuildTree(meshkernel::Location::Nodes);

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
            const auto foundNode = meshKernelState[meshKernelId].m_mesh2d->Node(nodeIndex);
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

            meshKernelUndoStack.Add(meshkernel::MeshTransformation::Compute(*meshKernelState[meshKernelId].m_mesh2d, transformation), meshKernelId);
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
            meshKernelUndoStack.Add(meshkernel::MeshTransformation::Compute(*meshKernelState[meshKernelId].m_mesh2d, translation), meshKernelId);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_casulli_refinement(int meshKernelId)
    {
        lastExitCode = meshkernel::ExitCode::Success;

        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            meshKernelUndoStack.Add(meshkernel::CasulliRefinement::Compute(*meshKernelState[meshKernelId].m_mesh2d), meshKernelId);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_casulli_refinement_on_polygon(int meshKernelId, const GeometryList& polygons)
    {
        lastExitCode = meshkernel::ExitCode::Success;

        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            // Convert polygon date from GeometryList to Polygons
            auto polygonPoints = ConvertGeometryListToPointVector(polygons);
            const meshkernel::Polygons meshKernelPolygons(polygonPoints,
                                                          meshKernelState[meshKernelId].m_mesh2d->m_projection);

            meshKernelUndoStack.Add(meshkernel::CasulliRefinement::Compute(*meshKernelState[meshKernelId].m_mesh2d,
                                                                           meshKernelPolygons),
                                    meshKernelId);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_casulli_derefinement_elements(int meshKernelId, GeometryList& elements)
    {
        lastExitCode = meshkernel::ExitCode::Success;

        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (elements.coordinates_x == nullptr || elements.coordinates_y == nullptr)
            {
                throw meshkernel::MeshKernelError("element coordinate list is null.");
            }

            meshkernel::Polygons polygons; // empty polygonal region, i.e. the entire mesh

            std::vector<meshkernel::Point> elementCentres(meshkernel::CasulliDeRefinement::ElementsToDelete(*meshKernelState[meshKernelId].m_mesh2d, polygons));

            elements.num_coordinates = static_cast<int>(elementCentres.size());

            for (size_t i = 0; i < elementCentres.size(); ++i)
            {
                elements.coordinates_x[i] = elementCentres[i].x;
                elements.coordinates_y[i] = elementCentres[i].y;
            }
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_casulli_derefinement_elements_on_polygon(int meshKernelId, const GeometryList& polygonGeometry, GeometryList& elements)
    {
        lastExitCode = meshkernel::ExitCode::Success;

        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (elements.coordinates_x == nullptr || elements.coordinates_y == nullptr)
            {
                throw meshkernel::MeshKernelError("element coordinate list is null.");
            }

            // Convert polygon date from GeometryList to Polygons
            auto polygonPoints = ConvertGeometryListToPointVector(polygonGeometry);
            const meshkernel::Polygons polygons(polygonPoints,
                                                meshKernelState[meshKernelId].m_mesh2d->m_projection);

            std::vector<meshkernel::Point> elementCentres(meshkernel::CasulliDeRefinement::ElementsToDelete(*meshKernelState[meshKernelId].m_mesh2d, polygons));

            elements.num_coordinates = static_cast<int>(elementCentres.size());

            for (size_t i = 0; i < elementCentres.size(); ++i)
            {
                elements.coordinates_x[i] = elementCentres[i].x;
                elements.coordinates_y[i] = elementCentres[i].y;
            }
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_casulli_derefinement(int meshKernelId)
    {
        lastExitCode = meshkernel::ExitCode::Success;

        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            meshKernelUndoStack.Add(meshkernel::CasulliDeRefinement::Compute(*meshKernelState[meshKernelId].m_mesh2d), meshKernelId);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_mesh2d_casulli_derefinement_on_polygon(int meshKernelId, const GeometryList& polygons)
    {
        lastExitCode = meshkernel::ExitCode::Success;

        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            // Convert polygon date from GeometryList to Polygons
            auto polygonPoints = ConvertGeometryListToPointVector(polygons);
            const meshkernel::Polygons meshKernelPolygons(polygonPoints,
                                                          meshKernelState[meshKernelId].m_mesh2d->m_projection);
            meshKernelUndoStack.Add(meshkernel::CasulliDeRefinement::Compute(*meshKernelState[meshKernelId].m_mesh2d,
                                                                             meshKernelPolygons),
                                    meshKernelId);
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
            auto const polygon = meshkernel::Polygons(polygonNodesVector, meshKernelState[meshKernelId].m_mesh2d->m_projection);
            auto landBoundary = meshkernel::LandBoundaries(landBoundariesNodeVector, *meshKernelState[meshKernelId].m_mesh2d, polygon);
            bool const triangulateFaces = isTriangulationRequired == 0 ? false : true;
            bool const projectToLandBoundary = projectToLandBoundaryRequired == 0 ? false : true;

            const meshkernel::FlipEdges flipEdges(*meshKernelState[meshKernelId].m_mesh2d, landBoundary, triangulateFaces, projectToLandBoundary);

            meshKernelUndoStack.Add(flipEdges.Compute(), meshKernelId);
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

            std::unique_ptr<meshkernel::CompoundUndoAction> undoSmallFlowEdges = meshkernel::CompoundUndoAction::Create();
            undoSmallFlowEdges->Add(meshKernelState[meshKernelId].m_mesh2d->DeleteSmallFlowEdges(smallFlowEdgesThreshold));
            undoSmallFlowEdges->Add(meshKernelState[meshKernelId].m_mesh2d->DeleteSmallTrianglesAtBoundaries(minFractionalAreaTriangles));
            meshKernelUndoStack.Add(std::move(undoSmallFlowEdges), meshKernelId);
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

            auto undoAction = MKStateUndoAction::Create(meshKernelState[meshKernelId]);

            // Execute
            meshKernelState[meshKernelId].m_contacts = std::make_unique<meshkernel::Contacts>(*meshKernelState[meshKernelId].m_mesh1d, *meshKernelState[meshKernelId].m_mesh2d);
            meshKernelState[meshKernelId].m_contacts->ComputeSingleContacts(meshKernel1DNodeMask, meshKernelPolygons, projectionFactor);

            meshKernelUndoStack.Add(std::move(undoAction), meshKernelId);
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

            auto undoAction = MKStateUndoAction::Create(meshKernelState[meshKernelId]);

            // Execute
            meshKernelState[meshKernelId].m_contacts = std::make_unique<meshkernel::Contacts>(*meshKernelState[meshKernelId].m_mesh1d, *meshKernelState[meshKernelId].m_mesh2d);
            meshKernelState[meshKernelId].m_contacts->ComputeMultipleContacts(meshKernel1DNodeMask);

            meshKernelUndoStack.Add(std::move(undoAction), meshKernelId);
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

            auto undoAction = MKStateUndoAction::Create(meshKernelState[meshKernelId]);

            // Execute
            meshKernelState[meshKernelId].m_contacts = std::make_unique<meshkernel::Contacts>(*meshKernelState[meshKernelId].m_mesh1d, *meshKernelState[meshKernelId].m_mesh2d);
            meshKernelState[meshKernelId].m_contacts->ComputeContactsWithPolygons(meshKernel1DNodeMask,
                                                                                  meshKernelPolygons);

            meshKernelUndoStack.Add(std::move(undoAction), meshKernelId);
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

            auto undoAction = MKStateUndoAction::Create(meshKernelState[meshKernelId]);

            // Execute
            meshKernelState[meshKernelId].m_contacts = std::make_unique<meshkernel::Contacts>(*meshKernelState[meshKernelId].m_mesh1d, *meshKernelState[meshKernelId].m_mesh2d);
            meshKernelState[meshKernelId].m_contacts->ComputeContactsWithPoints(meshKernel1DNodeMask, meshKernelPoints);

            meshKernelUndoStack.Add(std::move(undoAction), meshKernelId);
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

            // Check mesh has been initialised

            // convert raw arrays to containers
            const auto edges2d = meshkernel::ConvertToEdgeNodesVector(mesh2d.num_edges,
                                                                      mesh2d.edge_nodes);

            const auto nodes2d = meshkernel::ConvertToNodesVector(mesh2d.num_nodes,
                                                                  mesh2d.node_x,
                                                                  mesh2d.node_y);

            std::unique_ptr<meshkernel::Mesh2D> meshToConnect;

            if (mesh2d.num_faces > 0 && mesh2d.face_nodes != nullptr && mesh2d.nodes_per_face != nullptr)
            {

                const auto face_nodes = meshkernel::ConvertToFaceNodesVector(mesh2d.num_faces, mesh2d.face_nodes, mesh2d.nodes_per_face);

                std::vector<meshkernel::UInt> num_face_nodes;
                num_face_nodes.reserve(mesh2d.num_faces);

                for (auto n = 0; n < mesh2d.num_faces; n++)
                {
                    num_face_nodes.emplace_back(static_cast<meshkernel::UInt>(mesh2d.nodes_per_face[n]));
                }

                meshToConnect = std::make_unique<meshkernel::Mesh2D>(edges2d,
                                                                     nodes2d,
                                                                     face_nodes,
                                                                     num_face_nodes,
                                                                     meshKernelState[meshKernelId].m_projection);
            }
            else
            {
                meshToConnect = std::make_unique<meshkernel::Mesh2D>(edges2d, nodes2d, meshKernelState[meshKernelId].m_projection);
            }

            const auto mergedMeshes = meshkernel::Mesh2D::Merge(*meshKernelState[meshKernelId].m_mesh2d, *meshToConnect);
            // Keep existing mesh to restore with undo
            auto undoAction = meshkernel::FullUnstructuredGridUndo::Create(*meshKernelState[meshKernelId].m_mesh2d);
            // The undo information collected from the ConnectMeshes::Compute is not needed here.
            [[maybe_unused]] auto undo = meshkernel::ConnectMeshes::Compute(*mergedMeshes, searchFraction);
            meshKernelState[meshKernelId].m_mesh2d->SetNodes(mergedMeshes->Nodes());
            meshKernelState[meshKernelId].m_mesh2d->SetEdges(mergedMeshes->Edges());
            meshKernelState[meshKernelId].m_mesh2d->Administrate();
            meshKernelUndoStack.Add(std::move(undoAction), meshKernelId);
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

            auto undoAction = MKStateUndoAction::Create(meshKernelState[meshKernelId]);

            // Execute
            meshKernelState[meshKernelId].m_contacts = std::make_unique<meshkernel::Contacts>(*meshKernelState[meshKernelId].m_mesh1d, *meshKernelState[meshKernelId].m_mesh2d);
            meshKernelState[meshKernelId].m_contacts->ComputeBoundaryContacts(meshKernel1DNodeMask, meshKernelPolygons, searchRadius);

            meshKernelUndoStack.Add(std::move(undoAction), meshKernelId);
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
            meshKernelUndoStack.Add(curvilinearGridRefinement.Compute(), meshKernelId);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_full_refine(int meshKernelId,
                                                    int mRefinement,
                                                    int nRefinement)
    {
        lastExitCode = meshkernel::ExitCode::Success;

        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (!meshKernelState[meshKernelId].m_curvilinearGrid->IsValid())
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id is valid, but the expected curvilinear grid is not a valid grid");
            }

            if (mRefinement <= 0 || nRefinement <= 0)
            {
                throw meshkernel::MeshKernelError("Invalid mesh refinement factors: m-refinement {}, n-refinement {} ",
                                                  mRefinement, nRefinement);
            }

            meshkernel::CurvilinearGridFullRefinement gridRefinement;
            meshKernelUndoStack.Add(gridRefinement.Compute(*meshKernelState[meshKernelId].m_curvilinearGrid,
                                                           static_cast<meshkernel::UInt>(mRefinement),
                                                           static_cast<meshkernel::UInt>(nRefinement)),
                                    meshKernelId);
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
            meshKernelUndoStack.Add(curvilinearGridDeRefinement.Compute(), meshKernelId);
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

            auto undoAction = MKStateUndoAction::Create(meshKernelState[meshKernelId]);

            // Set the state
            meshKernelState[meshKernelId].m_curvilinearGrid = curvilinearGridFromSplinesTransfinite.Compute();

            meshKernelUndoStack.Add(std::move(undoAction), meshKernelId);
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

            auto undoAction = MKStateUndoAction::Create(meshKernelState[meshKernelId]);

            // set the curvilinear state
            meshKernelState[meshKernelId].m_curvilinearGrid = curvilinearGridFromPolygon.Compute(firstNode, secondNode, thirdNode, useFourthSideBool);

            meshKernelUndoStack.Add(std::move(undoAction), meshKernelId);
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

            auto undoAction = MKStateUndoAction::Create(meshKernelState[meshKernelId]);

            // set the curvilinear state
            meshKernelState[meshKernelId].m_curvilinearGrid = curvilinearGridFromPolygon.Compute(firstNode, secondNode, thirdNode);

            meshKernelUndoStack.Add(std::move(undoAction), meshKernelId);
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

            auto undoAction = MKStateUndoAction::Create(meshKernelState[meshKernelId]);

            // set the curvilinear state
            meshKernelState[meshKernelId].m_curvilinearGrid = curvilinearGridFromSplines.Compute();

            meshKernelUndoStack.Add(std::move(undoAction), meshKernelId);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_compute_grid_from_splines(int meshKernelId,
                                                                  const GeometryList& geometryListIn,
                                                                  const meshkernel::CurvilinearParameters& curvilinearParameters)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            meshkernel::Splines splines(meshKernelState[meshKernelId].m_projection);
            SetSplines(geometryListIn, splines);

            meshkernel::CurvilinearGridSplineToGrid splineToGrid;

            auto undoAction = MKStateUndoAction::Create(meshKernelState[meshKernelId]);

            //  the curvilinear grid
            meshKernelState[meshKernelId].m_curvilinearGrid = std::make_unique<meshkernel::CurvilinearGrid>(splineToGrid.Compute(splines, curvilinearParameters));

            meshKernelUndoStack.Add(std::move(undoAction), meshKernelId);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_compute_curvature(int meshKernelId, int direction, double* curvature)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (curvature == nullptr)
            {
                throw meshkernel::ConstraintError("The curvautre array is null");
            }

            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id, {}, does not exist.", meshKernelId);
            }

            if (meshKernelState[meshKernelId].m_curvilinearGrid == nullptr)
            {
                throw meshkernel::MeshKernelError("The curvilinear grid id, {}, does not exist.", meshKernelId);
            }

            meshkernel::CurvilinearDirection directionEnum = meshkernel::GetCurvilinearDirectionValue(direction);
            const meshkernel::CurvilinearGrid& grid = *meshKernelState[meshKernelId].m_curvilinearGrid;
            lin_alg::Matrix<double> curvatureMatrix;

            meshkernel::CurvilinearGridCurvature::Compute(grid, directionEnum, curvatureMatrix);
            Eigen::Map<lin_alg::Matrix<double>>(curvature, curvatureMatrix.rows(), curvatureMatrix.cols()) = curvatureMatrix;
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_compute_smoothness(int meshKernelId, int direction, double* smoothness)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (smoothness == nullptr)
            {
                throw meshkernel::ConstraintError("The smoothness array is null");
            }

            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id, {}, does not exist.", meshKernelId);
            }

            if (meshKernelState[meshKernelId].m_curvilinearGrid == nullptr)
            {
                throw meshkernel::MeshKernelError("The curvilinear grid id, {}, does not exist.", meshKernelId);
            }

            meshkernel::CurvilinearDirection directionEnum = meshkernel::GetCurvilinearDirectionValue(direction);
            const meshkernel::CurvilinearGrid& grid = *meshKernelState[meshKernelId].m_curvilinearGrid;
            lin_alg::Matrix<double> smoothnessMatrix;

            meshkernel::CurvilinearGridSmoothness::Compute(grid, directionEnum, smoothnessMatrix);
            Eigen::Map<lin_alg::Matrix<double>>(smoothness, smoothnessMatrix.rows(), smoothnessMatrix.cols()) = smoothnessMatrix;
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

            meshKernelState[meshKernelId].m_curvilinearGridFromSplines = std::make_unique<meshkernel::CurvilinearGridFromSplines>(spline, curvilinearParameters, splinesToCurvilinearParameters);

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

            auto undoAction = MKStateUndoAction::Create(meshKernelState[meshKernelId]);

            meshKernelState[meshKernelId].m_curvilinearGrid = meshKernelState[meshKernelId].m_curvilinearGridFromSplines->ComputeCurvilinearGridFromGridPoints();

            meshKernelUndoStack.Add(std::move(undoAction), meshKernelId);
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

            auto undoAction = MKStateUndoAction::Create(meshKernelState[meshKernelId]);

            meshKernelState[meshKernelId].m_curvilinearGrid = CreateRectangularCurvilinearGrid(makeGridParameters, meshKernelState[meshKernelId].m_projection);

            meshKernelUndoStack.Add(std::move(undoAction), meshKernelId);
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

            auto undoAction = MKStateUndoAction::Create(meshKernelState[meshKernelId]);

            meshKernelState[meshKernelId].m_curvilinearGrid = CreateRectangularCurvilinearGridFromPolygons(makeGridParameters,
                                                                                                           geometryList,
                                                                                                           meshKernelState[meshKernelId].m_projection);

            meshKernelUndoStack.Add(std::move(undoAction), meshKernelId);
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

            auto undoAction = MKStateUndoAction::Create(meshKernelState[meshKernelId]);

            meshKernelState[meshKernelId].m_curvilinearGrid = CreateRectangularCurvilinearGridOnExtension(makeGridParameters,
                                                                                                          meshKernelState[meshKernelId].m_projection);

            meshKernelUndoStack.Add(std::move(undoAction), meshKernelId);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_compute_circular_grid(int meshKernelId,
                                                              const meshkernel::MakeGridParameters& parameters)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGrid->IsValid())
            {
                throw meshkernel::MeshKernelError("The selected mesh already contains a valid grid.");
            }

            auto undoAction = MKStateUndoAction::Create(meshKernelState[meshKernelId]);
            *meshKernelState[meshKernelId].m_curvilinearGrid = meshkernel::CurvilinearGridGenerateCircularGrid::Compute(parameters,
                                                                                                                        meshKernelState[meshKernelId].m_projection);

            meshKernelUndoStack.Add(std::move(undoAction), meshKernelId);
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

            meshKernelState[meshKernelId].m_curvilinearGridOrthogonalization = std::make_unique<meshkernel::CurvilinearGridOrthogonalization>(*meshKernelState[meshKernelId].m_curvilinearGrid,
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
            meshKernelUndoStack.Add(meshKernelState[meshKernelId].m_curvilinearGridOrthogonalization->Compute(), meshKernelId);
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
            meshKernelUndoStack.Add(curvilinearGridSmoothing.Compute(), meshKernelId);
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

            auto undoAction = MKStateUndoAction::Create(meshKernelState[meshKernelId]);

            meshKernelState[meshKernelId].m_curvilinearGrid = curvilinearGridSmoothing.ComputeDirectional();

            meshKernelUndoStack.Add(std::move(undoAction), meshKernelId);
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

            meshKernelState[meshKernelId].m_curvilinearGridLineShift = std::make_unique<meshkernel::CurvilinearGridLineShift>(*meshKernelState[meshKernelId].m_curvilinearGrid);
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

            lin_alg::Matrix<meshkernel::Point> curviGridPoints(grid.num_n, grid.num_m);
            int nodeIndex = 0;
            for (int i = 0; i < grid.num_n; ++i)
            {
                for (int j = 0; j < grid.num_m; ++j)
                {

                    curviGridPoints(i, j) = meshkernel::Point(grid.node_x[nodeIndex], grid.node_y[nodeIndex]);
                    nodeIndex++;
                }
            }

            auto undoAction = MKStateUndoAction::Create(meshKernelState[meshKernelId]);

            const auto& projection = meshKernelState[meshKernelId].m_projection;
            meshKernelState[meshKernelId].m_curvilinearGrid = std::make_unique<meshkernel::CurvilinearGrid>(curviGridPoints, projection);

            meshKernelUndoStack.Add(std::move(undoAction), meshKernelId);
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
            meshKernelUndoStack.Add(meshKernelState[meshKernelId].m_curvilinearGridLineShift->MoveNode(fromPoint, toPoint), meshKernelId);
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

            meshKernelUndoStack.Add(meshKernelState[meshKernelId].m_curvilinearGridLineShift->Compute(), meshKernelId);
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

            meshKernelUndoStack.Add(meshKernelState[meshKernelId].m_curvilinearGrid->InsertFace(point), meshKernelId);
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

            if (!meshKernelState[meshKernelId].m_curvilinearGrid->IsValid())
            {
                throw meshkernel::MeshKernelError("Invalid curvilinear grid");
            }

            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() > 0 &&
                meshKernelState[meshKernelId].m_curvilinearGrid->projection() != meshKernelState[meshKernelId].m_mesh2d->m_projection)
            {
                throw meshkernel::MeshKernelError("The existing mesh2d projection is not equal to the curvilinear grid projection");
            }

            const auto edges = meshKernelState[meshKernelId].m_curvilinearGrid->ComputeEdges();
            const auto nodes = meshKernelState[meshKernelId].m_curvilinearGrid->ComputeNodes();

            // The undo action for conversion of clg to m2d is made in two steps
            std::unique_ptr<meshkernel::CompoundUndoAction> undoAction = meshkernel::CompoundUndoAction::Create();

            // 1. Keep the mesh kernel state to be able to restore (manily) the curvilinear grid and (secondly) other members.
            undoAction->Add(MKStateUndoAction::Create(meshKernelState[meshKernelId]));

            // 2. Keep track of the undo required to restore the mesh2d to its pre-converted state.
            undoAction->Add(meshKernelState[meshKernelId].m_mesh2d->Join(meshkernel::Mesh2D(edges, nodes, meshKernelState[meshKernelId].m_curvilinearGrid->projection())));

            // curvilinear grid must be reset to an empty curvilinear grid
            meshKernelState[meshKernelId].m_curvilinearGrid = std::make_unique<meshkernel::CurvilinearGrid>();
            meshKernelUndoStack.Add(std::move(undoAction), meshKernelId);
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

            meshKernelUndoStack.Add(curvilinearDeleteExterior.Compute(), meshKernelId);
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

            meshKernelUndoStack.Add(curvilinearDeleteInterior.Compute(), meshKernelId);
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

            meshKernelUndoStack.Add(curvilinearLineAttractionRepulsion.Compute(), meshKernelId);
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

            meshKernelUndoStack.Add(curvilinearGridLineMirror.Compute(), meshKernelId);
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

            meshKernelUndoStack.Add(meshKernelState[meshKernelId].m_curvilinearGrid->DeleteNode({xPointCoordinate, yPointCoordinate}), meshKernelId);
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

            meshKernelUndoStack.Add(meshKernelState[meshKernelId].m_curvilinearGrid->MoveNode(fromPoint, toPoint), meshKernelId);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_snap_to_landboundary(int meshKernelId,
                                                             const GeometryList& land,
                                                             double sectionControlPoint1x,
                                                             double sectionControlPoint1y,
                                                             double sectionControlPoint2x,
                                                             double sectionControlPoint2y,
                                                             double regionControlPointX,
                                                             double regionControlPointY)
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

            std::vector<meshkernel::Point> landBoundaryPoints(ConvertGeometryListToPointVector(land));
            std::vector<meshkernel::Point> controlPoints(regionControlPointX == meshkernel::constants::missing::doubleValue ? 2 : 3);

            controlPoints[0] = {sectionControlPoint1x, sectionControlPoint1y};
            controlPoints[1] = {sectionControlPoint2x, sectionControlPoint2y};

            if (regionControlPointX != meshkernel::constants::missing::doubleValue)
            {
                controlPoints[2] = {regionControlPointX, regionControlPointY};
            }

            meshkernel::LandBoundary landBoundary(landBoundaryPoints);

            //--------------------------------
            // Snap curvilinear grid to the land boundary
            meshkernel::CurvilinearGridSnapGridToLandBoundary gridSnapping(*meshKernelState[meshKernelId].m_curvilinearGrid, landBoundary, controlPoints);
            meshKernelUndoStack.Add(gridSnapping.Compute(), meshKernelId);
        }
        catch (...)
        {
            lastExitCode = HandleException();
        }
        return lastExitCode;
    }

    MKERNEL_API int mkernel_curvilinear_snap_to_spline(int meshKernelId,
                                                       const GeometryList& spline,
                                                       double sectionControlPoint1x,
                                                       double sectionControlPoint1y,
                                                       double sectionControlPoint2x,
                                                       double sectionControlPoint2y,
                                                       double regionControlPointX,
                                                       double regionControlPointY)
    {
        lastExitCode = meshkernel::ExitCode::Success;

        try
        {
            if (!meshKernelState.contains(meshKernelId))
            {
                throw meshkernel::MeshKernelError("The selected mesh kernel id does not exist.");
            }

            if (spline.num_coordinates == 0)
            {
                throw meshkernel::MeshKernelError("Spline has no point values.");
            }

            if (spline.coordinates_x == nullptr || spline.coordinates_y == nullptr)
            {
                throw meshkernel::MeshKernelError("Spline data is null.");
            }

            std::vector<meshkernel::Point> splinePoints(ConvertGeometryListToPointVector(spline));
            std::vector<meshkernel::Point> controlPoints(regionControlPointX == meshkernel::constants::missing::doubleValue ? 2 : 3);

            controlPoints[0] = {sectionControlPoint1x, sectionControlPoint1y};
            controlPoints[1] = {sectionControlPoint2x, sectionControlPoint2y};

            if (regionControlPointX != meshkernel::constants::missing::doubleValue)
            {
                controlPoints[2] = {regionControlPointX, regionControlPointY};
            }

            meshkernel::Splines mkSpline(meshKernelState[meshKernelId].m_curvilinearGrid->projection());

            mkSpline.AddSpline(splinePoints);

            //--------------------------------
            // Snap curvilinear grid to the spline

            meshkernel::CurvilinearGridSnapGridToSpline gridSnapping(*meshKernelState[meshKernelId].m_curvilinearGrid, mkSpline, controlPoints);
            meshKernelUndoStack.Add(gridSnapping.Compute(), meshKernelId);
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
            auto const meshLocation = static_cast<meshkernel::Location>(locationType);
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
            if (meshLocation == meshkernel::Location::Nodes)
            {
                interpolationResults = averaging.GetNodeResults();
            }
            else if (meshLocation == meshkernel::Location::Edges)
            {
                interpolationResults = averaging.GetEdgeResults();
            }
            else if (meshLocation == meshkernel::Location::Faces)
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
            auto const meshLocation = static_cast<meshkernel::Location>(locationType);
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
        type = static_cast<int>(meshkernel::Location::Edges);
        return lastExitCode;
    }
    MKERNEL_API int mkernel_get_nodes_location_type(int& type)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        type = static_cast<int>(meshkernel::Location::Nodes);
        return lastExitCode;
    }
    MKERNEL_API int mkernel_get_faces_location_type(int& type)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        type = static_cast<int>(meshkernel::Location::Faces);
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

    MKERNEL_API int mkernel_get_interpolation_type_short(int& type)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        type = static_cast<int>(meshkernel::InterpolationDataTypes::Short);
        return lastExitCode;
    }

    MKERNEL_API int mkernel_get_interpolation_type_float(int& type)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        type = static_cast<int>(meshkernel::InterpolationDataTypes::Float);
        return lastExitCode;
    }

    MKERNEL_API int mkernel_get_interpolation_type_int(int& type)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        type = static_cast<int>(meshkernel::InterpolationDataTypes::Int);
        return lastExitCode;
    }

    MKERNEL_API int mkernel_get_interpolation_type_double(int& type)
    {
        lastExitCode = meshkernel::ExitCode::Success;
        type = static_cast<int>(meshkernel::InterpolationDataTypes::Double);
        return lastExitCode;
    }

} // namespace meshkernelapi
