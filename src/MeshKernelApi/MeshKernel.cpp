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

#include <map>
#include <stdexcept>
#include <string.h>
#include <vector>

#include <MeshKernel/AveragingInterpolation.hpp>
#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Contacts.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridAlgorithm.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridCreateUniform.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridDeRefinement.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridFromPolygon.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridFromSplines.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridFromSplinesTransfinite.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLineAttractionRepulsion.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLineMirror.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridOrthogonalization.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridRefinement.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridSmoothing.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/FlipEdges.hpp>
#include <MeshKernel/LandBoundaries.hpp>
#include <MeshKernel/Mesh1D.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/MeshRefinement.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/OrthogonalizationAndSmoothing.hpp>
#include <MeshKernel/Orthogonalizer.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/Smoother.hpp>
#include <MeshKernel/Splines.hpp>
#include <MeshKernel/TriangulationInterpolation.hpp>

#include <MeshKernelApi/CurvilinearParameters.hpp>
#include <MeshKernelApi/MeshKernel.hpp>
#include <MeshKernelApi/SplinesToCurvilinearParameters.hpp>
#include <MeshKernelApi/State.hpp>
#include <MeshKernelApi/Utils.hpp>

namespace meshkernelapi
{
    // The state held by MeshKernel
    static std::map<int, MeshKernelState> meshKernelState;
    static int meshKernelStateCounter = 0;

    // Error state
    static char exceptionMessage[512] = "";
    static meshkernel::MeshGeometryError meshGeometryError = meshkernel::MeshGeometryError();

    int HandleExceptions(const std::exception_ptr exceptionPtr)
    {
        try
        {
            std::rethrow_exception(exceptionPtr);
        }
        catch (const meshkernel::MeshGeometryError& e)
        {
            meshGeometryError = e;
            strcpy(exceptionMessage, e.what());
            return InvalidGeometry;
        }
        catch (const std::exception& e)
        {
            strcpy(exceptionMessage, e.what());
            return Exception;
        }
    }

    MKERNEL_API int mkernel_allocate_state(int isGeographic, int& meshKernelId)
    {
        meshKernelId = meshKernelStateCounter++;
        auto const projection = static_cast<meshkernel::Projection>(isGeographic);
        meshKernelState.insert({meshKernelId, MeshKernelState(projection)});
        return Success;
    };

    MKERNEL_API int mkernel_deallocate_state(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            meshKernelState.erase(meshKernelId);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_delete(int meshKernelId, const GeometryList& polygon, int deletionOption, int invertDeletion)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() <= 0)
            {
                throw std::invalid_argument("MeshKernel: The 2d mesh contains no nodes.");
            }

            auto polygonPoints = ConvertGeometryListToPointVector(polygon);

            const bool invertDeletionBool = invertDeletion == 1 ? true : false;
            const meshkernel::Polygons meshKernelPolygon(polygonPoints, meshKernelState[meshKernelId].m_mesh2d->m_projection);
            meshKernelState[meshKernelId].m_mesh2d->DeleteMesh(meshKernelPolygon, deletionOption, invertDeletionBool);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_set(int meshKernelId,
                                       const Mesh2D& mesh2d)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            // convert raw arrays to containers
            const auto edges2d = meshkernel::ConvertToEdgeNodesVector(mesh2d.num_edges,
                                                                      mesh2d.edge_nodes);
            const auto nodes2d = meshkernel::ConvertToNodesVector(mesh2d.num_nodes,
                                                                  mesh2d.node_x,
                                                                  mesh2d.node_y);

            // Do not change the pointer, just the object it is pointing to
            *meshKernelState[meshKernelId].m_mesh2d = meshkernel::Mesh2D(edges2d, nodes2d, meshKernelState[meshKernelId].m_projection);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh1d_set(int meshKernelId,
                                       const Mesh1D& mesh1d)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
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
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_dimensions(int meshKernelId,
                                                  Mesh2D& mesh2d)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            meshKernelState[meshKernelId].m_mesh2d->Administrate(meshkernel::Mesh2D::AdministrationOption::AdministrateMeshEdgesAndFaces);
            SetMesh2dDimensions(meshKernelState[meshKernelId].m_mesh2d, mesh2d);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_data(int meshKernelId, Mesh2D& mesh2d)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            SetMesh2d(meshKernelState[meshKernelId].m_mesh2d, mesh2d);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh1d_get_dimensions(int meshKernelId,
                                                  Mesh1D& mesh1d)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            mesh1d.num_nodes = static_cast<int>(meshKernelState[meshKernelId].m_mesh1d->GetNumNodes());
            mesh1d.num_edges = static_cast<int>(meshKernelState[meshKernelId].m_mesh1d->GetNumEdges());
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh1d_get_data(int meshKernelId,
                                            Mesh1D& mesh1d)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            SetMesh1d(meshKernelState[meshKernelId].m_mesh1d, mesh1d);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_get_dimensions(int meshKernelId,
                                                       CurvilinearGrid& curvilinearGrid)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            curvilinearGrid.num_nodes = static_cast<int>(meshKernelState[meshKernelId].m_curvilinearGrid->GetNumNodes());
            curvilinearGrid.num_edges = static_cast<int>(meshKernelState[meshKernelId].m_curvilinearGrid->GetNumEdges());
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_get_data(int meshKernelId,
                                                 CurvilinearGrid& curvilinearGrid)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            meshKernelState[meshKernelId].m_curvilinearGrid->SetFlatCopies();
            SetCurvilinear(meshKernelState[meshKernelId].m_curvilinearGrid, curvilinearGrid);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_contacts_get_dimensions(int meshKernelId,
                                                    Contacts& contacts)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            contacts.num_contacts = static_cast<int>(meshKernelState[meshKernelId].m_contacts->m_mesh2dIndices.size());
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_contacts_get_data(int meshKernelId,
                                              Contacts& contacts)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            for (auto i = 0; i < contacts.num_contacts; ++i)
            {
                contacts.mesh1d_indices[i] = static_cast<int>(meshKernelState[meshKernelId].m_contacts->m_mesh1dIndices[i]);
                contacts.mesh2d_indices[i] = static_cast<int>(meshKernelState[meshKernelId].m_contacts->m_mesh2dIndices[i]);
            }
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_count_hanging_edges(int meshKernelId, int& numHangingEdges)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            meshKernelState[meshKernelId].m_mesh2d->Administrate(meshkernel::Mesh2D::AdministrationOption::AdministrateMeshEdges);
            const auto hangingEdges = meshKernelState[meshKernelId].m_mesh2d->GetHangingEdges();
            numHangingEdges = static_cast<int>(hangingEdges.size());
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_hanging_edges(int meshKernelId, int* edges)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            const auto hangingEdges = meshKernelState[meshKernelId].m_mesh2d->GetHangingEdges();
            for (auto i = 0; i < hangingEdges.size(); ++i)
            {
                edges[i] = static_cast<int>(hangingEdges[i]);
            }
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_delete_hanging_edges(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            meshKernelState[meshKernelId].m_mesh2d->DeleteHangingEdges();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_compute_orthogonalization(int meshKernelId,
                                                             int projectToLandBoundaryOption,
                                                             const OrthogonalizationParameters& orthogonalizationParameters,
                                                             const GeometryList& selectingPolygon,
                                                             const GeometryList& landBoundaries)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() <= 0)
            {
                throw std::invalid_argument("MeshKernel: The 2d mesh contains no nodes.");
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
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_initialize_orthogonalization(int meshKernelId,
                                                                int projectToLandBoundaryOption,
                                                                OrthogonalizationParameters& orthogonalizationParameters,
                                                                const GeometryList& selectingPolygon,
                                                                const GeometryList& landBoundaries)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() <= 0)
            {
                return exitCode;
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
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_prepare_outer_iteration_orthogonalization(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() <= 0)
            {
                return exitCode;
            }

            meshKernelState[meshKernelId].m_meshOrthogonalization->PrepareOuterIteration();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_compute_inner_ortogonalization_iteration(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() <= 0)
            {
                return exitCode;
            }

            meshKernelState[meshKernelId].m_meshOrthogonalization->Solve();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_finalize_inner_ortogonalization_iteration(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() <= 0)
            {
                return exitCode;
            }

            meshKernelState[meshKernelId].m_meshOrthogonalization->FinalizeOuterIteration();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_delete_orthogonalization(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() <= 0)
            {
                return exitCode;
            }

            meshKernelState[meshKernelId].m_meshOrthogonalization.reset();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_orthogonality(int meshKernelId, GeometryList& geometryList)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() <= 0)
            {
                return exitCode;
            }

            const auto result = meshKernelState[meshKernelId].m_mesh2d->GetOrthogonality();

            if (geometryList.num_coordinates != result.size())
            {
                throw std::invalid_argument("MeshKernel: The value array has not the same size of the result array storing the orthogonality values at the edges");
            }

            for (auto i = 0; i < geometryList.num_coordinates; ++i)
            {
                geometryList.values[i] = result[i];
            }
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_smoothness(int meshKernelId, GeometryList& geometryList)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() <= 0)
            {
                return exitCode;
            }

            const auto result = meshKernelState[meshKernelId].m_mesh2d->GetSmoothness();

            for (auto i = 0; i < geometryList.num_coordinates; ++i)
            {
                geometryList.values[i] = result[i];
            }
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_splines(const GeometryList& geometryListIn,
                                        GeometryList& geometryListOut,
                                        int numberOfPointsBetweenNodes)
    {
        int exitCode = Success;
        try
        {
            if (geometryListIn.num_coordinates == 0)
            {
                throw std::invalid_argument("MeshKernel: The number of coordinates of the given geometry is zero.");
            }

            std::vector<meshkernel::Point> splines(geometryListIn.num_coordinates);
            for (auto i = 0; i < geometryListIn.num_coordinates; i++)
            {
                splines[i].x = geometryListIn.coordinates_x[i];
                splines[i].y = geometryListIn.coordinates_y[i];
            }

            const auto indices = FindIndices(splines, 0, splines.size(), meshkernel::doubleMissingValue);
            const auto numSplines = indices.size();

            int index = 0;
            for (auto s = 0; s < numSplines; s++)
            {
                std::vector<meshkernel::Point> coordinates(splines.begin() + indices[s][0], splines.begin() + static_cast<int>(indices[s][1]) + 1);
                const int numNodes = static_cast<int>(indices[s][1]) - static_cast<int>(indices[s][0]) + 1;
                const auto coordinatesDerivatives = meshkernel::Splines::SecondOrderDerivative(coordinates, 0, coordinates.size() - 1);

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
                        geometryListOut.values[index] = meshkernel::doubleMissingValue;
                        index++;
                    }
                }

                geometryListOut.coordinates_x[index] = coordinates.back().x;
                geometryListOut.coordinates_y[index] = coordinates.back().y;
                geometryListOut.values[index] = meshkernel::doubleMissingValue;
                index++;

                geometryListOut.coordinates_x[index] = meshkernel::doubleMissingValue;
                geometryListOut.coordinates_y[index] = meshkernel::doubleMissingValue;
                geometryListOut.values[index] = meshkernel::doubleMissingValue;
                index++;
            }
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_make_mesh_from_polygon(int meshKernelId, const GeometryList& polygonPoints)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
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
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_make_mesh_from_samples(int meshKernelId, const GeometryList& samples)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            auto sampleVector = ConvertGeometryListToPointVector(samples);

            meshkernel::Polygons polygon;
            const meshkernel::Mesh2D mesh(sampleVector, polygon, meshKernelState[meshKernelId].m_mesh2d->m_projection);
            *meshKernelState[meshKernelId].m_mesh2d += mesh;
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_mesh_boundaries_as_polygons(int meshKernelId, GeometryList& boundaryPolygons)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            const std::vector<meshkernel::Point> polygonNodes;
            const auto meshBoundaryPolygon = meshKernelState[meshKernelId].m_mesh2d->MeshBoundaryToPolygon(polygonNodes);

            ConvertPointVectorToGeometryList(meshBoundaryPolygon, boundaryPolygons);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_count_mesh_boundaries_as_polygons(int meshKernelId, int& numberOfPolygonNodes)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            const std::vector<meshkernel::Point> polygonNodes;
            const auto meshBoundaryPolygon = meshKernelState[meshKernelId].m_mesh2d->MeshBoundaryToPolygon(polygonNodes);
            numberOfPolygonNodes = static_cast<int>(meshBoundaryPolygon.size() - 1); // last value is a separator
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_polygon_refine(int meshKernelId, const GeometryList& polygonToRefine, int firstNodeIndex, int secondNodeIndex, double targetEdgeLength, GeometryList& refinedPolygon)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            auto const polygonVector = ConvertGeometryListToPointVector(polygonToRefine);

            const meshkernel::Polygons polygon(polygonVector, meshKernelState[meshKernelId].m_mesh2d->m_projection);
            auto const refinementResult = polygon.RefineFirstPolygon(firstNodeIndex, secondNodeIndex, targetEdgeLength);

            ConvertPointVectorToGeometryList(refinementResult, refinedPolygon);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_polygon_count_refine(int meshKernelId,
                                                 const GeometryList& polygonToRefine,
                                                 int firstIndex,
                                                 int secondIndex,
                                                 double distance,
                                                 int& numberOfPolygonNodes)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            auto const polygonVector = ConvertGeometryListToPointVector(polygonToRefine);

            const meshkernel::Polygons polygon(polygonVector, meshKernelState[meshKernelId].m_mesh2d->m_projection);

            const auto refinedPolygon = polygon.RefineFirstPolygon(firstIndex, secondIndex, distance);

            numberOfPolygonNodes = int(refinedPolygon.size());
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_merge_nodes(int meshKernelId, const GeometryList& geometryListIn, double mergingDistance)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            auto const polygonVector = ConvertGeometryListToPointVector(geometryListIn);

            const meshkernel::Polygons polygon(polygonVector, meshKernelState[meshKernelId].m_mesh2d->m_projection);

            meshKernelState[meshKernelId].m_mesh2d->MergeNodesInPolygon(polygon, mergingDistance);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_merge_two_nodes(int meshKernelId, int firstNode, int secondNode)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            meshKernelState[meshKernelId].m_mesh2d->MergeTwoNodes(firstNode, secondNode);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_nodes_in_polygons(int meshKernelId,
                                                         const GeometryList& geometryListIn,
                                                         int inside,
                                                         int* selectedNodes)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            auto const polygonVector = ConvertGeometryListToPointVector(geometryListIn);

            const meshkernel::Polygons polygon(polygonVector, meshKernelState[meshKernelId].m_mesh2d->m_projection);

            const bool selectInside = inside == 1 ? true : false;
            meshKernelState[meshKernelId].m_mesh2d->MaskNodesInPolygons(polygon, selectInside);

            int index = 0;
            for (auto i = 0; i < meshKernelState[meshKernelId].m_mesh2d->GetNumNodes(); ++i)
            {
                if (meshKernelState[meshKernelId].m_mesh2d->m_nodeMask[i] > 0)
                {
                    selectedNodes[index] = i;
                    index++;
                }
            }
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_count_nodes_in_polygons(int meshKernelId,
                                                           const GeometryList& geometryListIn,
                                                           int inside,
                                                           int& numberOfMeshNodes)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            auto const polygonVector = ConvertGeometryListToPointVector(geometryListIn);

            const meshkernel::Polygons polygon(polygonVector, meshKernelState[meshKernelId].m_mesh2d->m_projection);

            const bool selectInside = inside == 1 ? true : false;
            meshKernelState[meshKernelId].m_mesh2d->MaskNodesInPolygons(polygon, selectInside);

            numberOfMeshNodes = 0;
            for (auto i = 0; i < meshKernelState[meshKernelId].m_mesh2d->GetNumNodes(); ++i)
            {
                if (meshKernelState[meshKernelId].m_mesh2d->m_nodeMask[i] > 0)
                {
                    numberOfMeshNodes++;
                }
            }
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_insert_edge(int meshKernelId, int startNode, int endNode, int& new_edge_index)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            new_edge_index = static_cast<int>(meshKernelState[meshKernelId].m_mesh2d->ConnectNodes(startNode, endNode));
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_insert_node(int meshKernelId, double xCoordinate, double yCoordinate, int& nodeIndex)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            meshkernel::Point const nodeCoordinateVector{xCoordinate, yCoordinate};

            nodeIndex = static_cast<int>(meshKernelState[meshKernelId].m_mesh2d->InsertNode(nodeCoordinateVector));
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_delete_node(int meshKernelId, int nodeIndex)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            meshKernelState[meshKernelId].m_mesh2d->DeleteNode(nodeIndex);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_move_node(int meshKernelId, double xCoordinate, double yCoordinate, int nodeIndex)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            meshkernel::Point newPosition{xCoordinate, yCoordinate};

            meshKernelState[meshKernelId].m_mesh2d->MoveNode(newPosition, nodeIndex);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_delete_edge(int meshKernelId, double xCoordinate, double yCoordinate)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            meshkernel::Point const point{xCoordinate, yCoordinate};

            const auto edgeIndex = meshKernelState[meshKernelId].m_mesh2d->FindEdgeCloseToAPoint(point);

            meshKernelState[meshKernelId].m_mesh2d->DeleteEdge(edgeIndex);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_edge(int meshKernelId, double xCoordinate, double yCoordinate, int& edgeIndex)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            meshkernel::Point const point{xCoordinate, yCoordinate};

            edgeIndex = static_cast<int>(meshKernelState[meshKernelId].m_mesh2d->FindEdgeCloseToAPoint(point));
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_polygon_get_offset(int meshKernelId, const GeometryList& geometryListIn, int inWard, double distance, GeometryList& geometryListOut)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            auto const polygonVector = ConvertGeometryListToPointVector(geometryListIn);

            const meshkernel::Polygons polygon(polygonVector, meshKernelState[meshKernelId].m_mesh2d->m_projection);

            const bool inWardBool = inWard == 1 ? true : false;
            const auto newPolygon = polygon.OffsetCopy(distance, inWardBool);

            ConvertPointVectorToGeometryList(newPolygon.m_nodes, geometryListOut);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_polygon_count_offset(int meshKernelId, const GeometryList& geometryListIn, int innerPolygon, double distance, int& numberOfPolygonNodes)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            auto polygonPoints = ConvertGeometryListToPointVector(geometryListIn);

            const bool innerPolygonBool = innerPolygon == 1 ? true : false;
            const meshkernel::Polygons polygon(polygonPoints, meshKernelState[meshKernelId].m_mesh2d->m_projection);
            const auto newPolygon = polygon.OffsetCopy(distance, innerPolygonBool);

            numberOfPolygonNodes = static_cast<int>(newPolygon.GetNumNodes());
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_refine_based_on_samples(int meshKernelId,
                                                           const GeometryList& samples,
                                                           double relativeSearchRadius,
                                                           int minimumNumSamples,
                                                           const MeshRefinementParameters& meshRefinementParameters)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() <= 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh has no nodes.");
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

            const auto averaging = std::make_shared<meshkernel::AveragingInterpolation>(meshKernelState[meshKernelId].m_mesh2d,
                                                                                        samplesVector,
                                                                                        averagingMethod,
                                                                                        meshkernel::MeshLocations::Faces,
                                                                                        relativeSearchRadius,
                                                                                        refineOutsideFace,
                                                                                        transformSamples,
                                                                                        static_cast<size_t>(minimumNumSamples));

            meshkernel::MeshRefinement meshRefinement(meshKernelState[meshKernelId].m_mesh2d, averaging, meshRefinementParameters);
            meshRefinement.Compute();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_refine_based_on_polygon(int meshKernelId, const GeometryList& geometryList, const MeshRefinementParameters& meshRefinementParameters)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() <= 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh has no nodes.");
            }

            auto points = ConvertGeometryListToPointVector(geometryList);

            const meshkernel::Polygons polygon(points, meshKernelState[meshKernelId].m_mesh2d->m_projection);

            meshkernel::MeshRefinement meshRefinement(meshKernelState[meshKernelId].m_mesh2d, polygon, meshRefinementParameters);
            meshRefinement.Compute();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_node_index(int meshKernelId, double xCoordinate, double yCoordinate, double searchRadius, int& nodeIndex)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() <= 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh has no nodes.");
            }

            meshkernel::Point const point{xCoordinate, yCoordinate};

            nodeIndex = static_cast<int>(meshKernelState[meshKernelId].m_mesh2d->FindNodeCloseToAPoint(point, searchRadius));
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_closest_node(int meshKernelId,
                                                    double xCoordinateIn,
                                                    double yCoordinateIn,
                                                    double searchRadius,
                                                    double& xCoordinateOut,
                                                    double& yCoordinateOut)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() <= 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh has no nodes.");
            }

            meshkernel::Point const point{xCoordinateIn, yCoordinateIn};

            const auto nodeIndex = meshKernelState[meshKernelId].m_mesh2d->FindNodeCloseToAPoint(point, searchRadius);

            // Set the node coordinate
            auto foundNode = meshKernelState[meshKernelId].m_mesh2d->m_nodes[nodeIndex];
            xCoordinateOut = foundNode.x;
            yCoordinateOut = foundNode.y;
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_polygon_get_included_points(int meshKernelId, const GeometryList& selectingPolygon, const GeometryList& polygonToSelect, GeometryList& selectionResults)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            auto const polygonVector = ConvertGeometryListToPointVector(selectingPolygon);

            auto const points = ConvertGeometryListToPointVector(polygonToSelect);

            const meshkernel::Polygons localPolygon(polygonVector, meshKernelState[meshKernelId].m_mesh2d->m_projection);

            for (auto i = 0; i < points.size(); i++)
            {
                selectionResults.values[i] = localPolygon.IsPointInPolygon(points[i], 0) ? 1.0 : 0.0;
            }
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_flip_edges(int meshKernelId,
                                              int isTriangulationRequired,
                                              int projectToLandBoundaryRequired,
                                              const GeometryList& selectingPolygon,
                                              const GeometryList& landBoundaries)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
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
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_count_small_flow_edge_centers(int meshKernelId, double smallFlowEdgesLengthThreshold, int& numSmallFlowEdges)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            const auto edgesCrossingSmallFlowEdges = meshKernelState[meshKernelId].m_mesh2d->GetEdgesCrossingSmallFlowEdges(smallFlowEdgesLengthThreshold);
            const auto smallFlowEdgeCenters = meshKernelState[meshKernelId].m_mesh2d->GetFlowEdgesCenters(edgesCrossingSmallFlowEdges);

            numSmallFlowEdges = static_cast<int>(smallFlowEdgeCenters.size());
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_small_flow_edge_centers(int meshKernelId, double smallFlowEdgesThreshold, GeometryList& result)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            const auto edgesCrossingSmallFlowEdges = meshKernelState[meshKernelId].m_mesh2d->GetEdgesCrossingSmallFlowEdges(smallFlowEdgesThreshold);
            const auto smallFlowEdgeCenters = meshKernelState[meshKernelId].m_mesh2d->GetFlowEdgesCenters(edgesCrossingSmallFlowEdges);

            ConvertPointVectorToGeometryList(smallFlowEdgeCenters, result);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_error(const char*& error_message)
    {
        error_message = exceptionMessage;
        return Success;
    }

    MKERNEL_API int mkernel_get_geometry_error(int& invalidIndex, int& type)
    {
        invalidIndex = static_cast<int>(meshGeometryError.m_invalidIndex);
        type = static_cast<int>(meshGeometryError.m_location);
        return Success;
    }

    MKERNEL_API int mkernel_mesh2d_count_obtuse_triangles(int meshKernelId, int& numObtuseTriangles)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            const auto obtuseTriangles = meshKernelState[meshKernelId].m_mesh2d->GetObtuseTrianglesCenters();

            numObtuseTriangles = static_cast<int>(obtuseTriangles.size());
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_get_obtuse_triangles_mass_centers(int meshKernelId, GeometryList& result)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            const auto obtuseTriangles = meshKernelState[meshKernelId].m_mesh2d->GetObtuseTrianglesCenters();

            ConvertPointVectorToGeometryList(obtuseTriangles, result);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_delete_small_flow_edges_and_small_triangles(int meshKernelId, double smallFlowEdgesThreshold, double minFractionalAreaTriangles)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            meshKernelState[meshKernelId].m_mesh2d->DeleteSmallFlowEdges(smallFlowEdgesThreshold);
            meshKernelState[meshKernelId].m_mesh2d->DeleteSmallTrianglesAtBoundaries(minFractionalAreaTriangles);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_contacts_compute_single(int meshKernelId,
                                                    const int* oneDNodeMask,
                                                    const GeometryList& polygons)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
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
            meshKernelState[meshKernelId].m_contacts->ComputeSingleContacts(meshKernel1DNodeMask,
                                                                            meshKernelPolygons);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_contacts_compute_multiple(int meshKernelId,
                                                      const int* oneDNodeMask)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
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
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_contacts_compute_with_polygons(int meshKernelId,
                                                           const int* oneDNodeMask,
                                                           const GeometryList& polygons)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
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
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }
    MKERNEL_API int mkernel_contacts_compute_with_points(int meshKernelId,
                                                         const int* oneDNodeMask,
                                                         const GeometryList& points)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
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
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_contacts_compute_boundary(int meshKernelId,
                                                      const int* oneDNodeMask,
                                                      const GeometryList& polygons,
                                                      double searchRadius)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
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
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_refine(int meshKernelId, double xLowerLeftCorner, double yLowerLeftCorner, double xUpperRightCorner, double yUpperRightCorner, int refinement)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            meshkernel::Point const firstPoint{xLowerLeftCorner, yLowerLeftCorner};
            meshkernel::Point const secondPoint{xUpperRightCorner, yUpperRightCorner};

            // Execute
            meshkernel::CurvilinearGridRefinement curvilinearGridRefinement(meshKernelState[meshKernelId].m_curvilinearGrid, refinement);
            curvilinearGridRefinement.SetBlock(firstPoint, secondPoint);
            meshKernelState[meshKernelId].m_curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(curvilinearGridRefinement.Compute());
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_derefine(int meshKernelId,
                                                 double xLowerLeftCorner,
                                                 double yLowerLeftCorner,
                                                 double xUpperRightCorner,
                                                 double yUpperRightCorner)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            meshkernel::Point const firstPoint{xLowerLeftCorner, yLowerLeftCorner};
            meshkernel::Point const secondPoint{xUpperRightCorner, yUpperRightCorner};

            // Execute
            meshkernel::CurvilinearGridDeRefinement curvilinearGridDeRefinement(meshKernelState[meshKernelId].m_curvilinearGrid);

            curvilinearGridDeRefinement.SetBlock(firstPoint, secondPoint);

            meshKernelState[meshKernelId].m_curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(curvilinearGridDeRefinement.Compute());
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_compute_transfinite_from_splines(int meshKernelId,
                                                                         const GeometryList& splines,
                                                                         const CurvilinearParameters& curvilinearParameters)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
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
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_compute_transfinite_from_polygon(int meshKernelId,
                                                                         const GeometryList& polygons,
                                                                         int firstNode,
                                                                         int secondNode,
                                                                         int thirdNode,
                                                                         int useFourthSide)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            auto polygonPoints = ConvertGeometryListToPointVector(polygons);

            const auto localPolygon = std::make_shared<meshkernel::Polygons>(polygonPoints, meshKernelState[meshKernelId].m_projection);

            const meshkernel::CurvilinearGridFromPolygon curvilinearGridFromPolygon(localPolygon);

            const bool useFourthSideBool = useFourthSide == 1 ? true : false;
            const auto curvilinearGrid = curvilinearGridFromPolygon.Compute(firstNode, secondNode, thirdNode, useFourthSideBool);

            // set the curvilinear state
            meshKernelState[meshKernelId].m_curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(curvilinearGrid);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_compute_transfinite_from_triangle(int meshKernelId,
                                                                          const GeometryList& polygon,
                                                                          int firstNode,
                                                                          int secondNode,
                                                                          int thirdNode)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            auto polygonPoints = ConvertGeometryListToPointVector(polygon);

            const auto localPolygon = std::make_shared<meshkernel::Polygons>(polygonPoints, meshKernelState[meshKernelId].m_projection);

            const meshkernel::CurvilinearGridFromPolygon curvilinearGridFromPolygon(localPolygon);

            const auto curvilinearGrid = curvilinearGridFromPolygon.Compute(firstNode, secondNode, thirdNode);

            // set the curvilinear state
            meshKernelState[meshKernelId].m_curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(curvilinearGrid);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_compute_orthogonal_grid_from_splines(int meshKernelId,
                                                                             const GeometryList& geometryListIn,
                                                                             const CurvilinearParameters& curvilinearParameters,
                                                                             const SplinesToCurvilinearParameters& splinesToCurvilinearParameters)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
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
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_initialize_orthogonal_grid_from_splines(int meshKernelId,
                                                                                const GeometryList& geometryList,
                                                                                const CurvilinearParameters& curvilinearParameters,
                                                                                const SplinesToCurvilinearParameters& splinesToCurvilinearParameters)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            auto spline = std::make_shared<meshkernel::Splines>(meshKernelState[meshKernelId].m_projection);
            SetSplines(geometryList, *spline);

            meshKernelState[meshKernelId].m_curvilinearGridFromSplines = std::make_shared<meshkernel::CurvilinearGridFromSplines>(spline, curvilinearParameters, splinesToCurvilinearParameters);

            meshKernelState[meshKernelId].m_curvilinearGridFromSplines->Initialize();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_iterate_orthogonal_grid_from_splines(int meshKernelId, int layer)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            if (meshKernelState[meshKernelId].m_curvilinearGridFromSplines == nullptr)
            {
                throw std::invalid_argument("MeshKernel: CurvilinearGridFromSplines not instantiated.");
            }

            meshKernelState[meshKernelId].m_curvilinearGridFromSplines->Iterate(layer);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_refresh_orthogonal_grid_from_splines(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGridFromSplines == nullptr)
            {
                throw std::invalid_argument("MeshKernel: CurvilinearGridFromSplines not instantiated.");
            }

            const auto curvilinearGrid = meshKernelState[meshKernelId].m_curvilinearGridFromSplines->ComputeCurvilinearGridFromGridPoints();

            meshKernelState[meshKernelId].m_curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(curvilinearGrid);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_delete_orthogonal_grid_from_splines(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGridFromSplines == nullptr)
            {
                throw std::invalid_argument("MeshKernel: CurvilinearGridFromSplines not instantiated.");
            }

            meshKernelState[meshKernelId].m_curvilinearGridFromSplines.reset();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_make_uniform(int meshKernelId,
                                                     const MakeMeshParameters& makeGridParameters,
                                                     const GeometryList& geometryList)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            auto polygonNodes = ConvertGeometryListToPointVector(geometryList);

            const auto polygon = std::make_shared<meshkernel::Polygons>(polygonNodes, meshKernelState[meshKernelId].m_projection);

            meshkernel::CurvilinearGridCreateUniform curvilinearGridCreateUniform(makeGridParameters, polygon);

            meshKernelState[meshKernelId].m_curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(curvilinearGridCreateUniform.Compute());
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_initialize_orthogonalize(int meshKernelId,
                                                                 const OrthogonalizationParameters& orthogonalizationParameters)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            meshKernelState[meshKernelId].m_curvilinearGridOrthogonalization = std::make_shared<meshkernel::CurvilinearGridOrthogonalization>(meshKernelState[meshKernelId].m_curvilinearGrid,
                                                                                                                                              orthogonalizationParameters);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_set_block_orthogonalize(int meshKernelId,
                                                                double xLowerLeftCorner,
                                                                double yLowerLeftCorner,
                                                                double xUpperRightCorner,
                                                                double yUpperRightCorner)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel state does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGridOrthogonalization == nullptr)
            {
                throw std::invalid_argument("MeshKernel: CurvilinearGridOrthogonalization not instantiated.");
            }

            meshkernel::Point firstPoint{xLowerLeftCorner, yLowerLeftCorner};
            meshkernel::Point secondPoint{xUpperRightCorner, yUpperRightCorner};

            // Execute
            meshKernelState[meshKernelId].m_curvilinearGridOrthogonalization->SetBlock(firstPoint, secondPoint);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_set_frozen_lines_orthogonalize(int meshKernelId,
                                                                       double xFirstGridLineNode,
                                                                       double yFirstGridLineNode,
                                                                       double xSecondGridLineNode,
                                                                       double ySecondGridLineNode)

    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel state does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGridOrthogonalization == nullptr)
            {
                throw std::invalid_argument("MeshKernel: CurvilinearGridOrthogonalization not instantiated.");
            }

            meshkernel::Point const firstPoint{xFirstGridLineNode, yFirstGridLineNode};
            meshkernel::Point const secondPoint{xSecondGridLineNode, ySecondGridLineNode};

            // Execute
            meshKernelState[meshKernelId].m_curvilinearGridOrthogonalization->SetLine(firstPoint, secondPoint);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_orthogonalize(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel state does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGridOrthogonalization == nullptr)
            {
                throw std::invalid_argument("MeshKernel: CurvilinearGridOrthogonalization not instantiated.");
            }

            // Execute
            meshKernelState[meshKernelId].m_curvilinearGridOrthogonalization->Compute();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_finalize_orthogonalize(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel state does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGridOrthogonalization == nullptr)
            {
                throw std::invalid_argument("MeshKernel: CurvilinearGridOrthogonalization not instantiated.");
            }

            meshKernelState[meshKernelId].m_curvilinearGridOrthogonalization.reset();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_smoothing(int meshKernelId,
                                                  int smoothingIterations,
                                                  double xLowerLeftCorner,
                                                  double yLowerLeftCorner,
                                                  double xUpperRightCorner,
                                                  double yUpperRightCorner)

    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel state does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGrid == nullptr)
            {
                throw std::invalid_argument("MeshKernel: Not a valid curvilinear grid instance.");
            }

            if (!meshKernelState[meshKernelId].m_curvilinearGrid->IsValid())
            {
                throw std::invalid_argument("MeshKernel: Not valid curvilinear grid.");
            }

            const meshkernel::Point firstPoint{xLowerLeftCorner, yLowerLeftCorner};
            const meshkernel::Point secondPoint{xUpperRightCorner, yUpperRightCorner};

            // Execute
            meshkernel::CurvilinearGridSmoothing curvilinearGridSmoothing(meshKernelState[meshKernelId].m_curvilinearGrid,
                                                                          static_cast<size_t>(smoothingIterations));

            curvilinearGridSmoothing.SetBlock(firstPoint, secondPoint);
            curvilinearGridSmoothing.Compute();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
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
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel state does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGrid == nullptr)
            {
                throw std::invalid_argument("MeshKernel: Not a valid curvilinear grid instance.");
            }

            if (!meshKernelState[meshKernelId].m_curvilinearGrid->IsValid())
            {
                throw std::invalid_argument("MeshKernel: Not valid curvilinear grid.");
            }

            meshkernel::Point const firstNode{xFirstGridlineNode, yFirstGridlineNode};
            meshkernel::Point const secondNode{xSecondGridLineNode, ySecondGridLineNode};
            meshkernel::Point const lowerLeft{xLowerLeftCornerSmoothingArea, yLowerLeftCornerSmoothingArea};
            meshkernel::Point const upperRight{xUpperRightCornerSmootingArea, yUpperRightCornerSmootingArea};

            // Execute
            meshkernel::CurvilinearGridSmoothing curvilinearGridSmoothing(meshKernelState[meshKernelId].m_curvilinearGrid, smoothingIterations);

            curvilinearGridSmoothing.SetLine(firstNode, secondNode);
            curvilinearGridSmoothing.SetBlock(lowerLeft, upperRight);

            curvilinearGridSmoothing.ComputeDirectional();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_initialize_line_shift(int meshKernelId)
    {

        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel state does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGrid == nullptr)
            {
                throw std::invalid_argument("MeshKernel: Not a valid curvilinear grid instance.");
            }

            if (!meshKernelState[meshKernelId].m_curvilinearGrid->IsValid())
            {
                throw std::invalid_argument("MeshKernel: Not valid curvilinear grid.");
            }

            meshKernelState[meshKernelId].m_curvilinearGridLineShift = std::make_shared<meshkernel::CurvilinearGridLineShift>(meshKernelState[meshKernelId].m_curvilinearGrid);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_set_line_line_shift(int meshKernelId,
                                                            double xFirstGridLineNode,
                                                            double yFirstGridLineNode,
                                                            double xSecondGridLineNode,
                                                            double ySecondGridLineNode)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel state does not exist.");
            }

            meshkernel::Point const firstNode{xFirstGridLineNode, yFirstGridLineNode};
            meshkernel::Point const secondNode{xSecondGridLineNode, ySecondGridLineNode};

            meshKernelState[meshKernelId].m_curvilinearGridLineShift->SetLine(firstNode, secondNode);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_set_block_line_shift(int meshKernelId,
                                                             double xLowerLeftCorner,
                                                             double yLowerLeftCorner,
                                                             double xUpperRightCorner,
                                                             double yUpperRightCorner)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel state does not exist.");
            }

            meshkernel::Point const lowerLeftPoint{xLowerLeftCorner, yLowerLeftCorner};
            meshkernel::Point const upperRightPoint{xUpperRightCorner, yUpperRightCorner};

            meshKernelState[meshKernelId].m_curvilinearGridLineShift->SetBlock(lowerLeftPoint, upperRightPoint);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_move_node_line_shift(int meshKernelId,
                                                             double xFromCoordinate,
                                                             double yFromCoordinate,
                                                             double xToCoordinate,
                                                             double yToCoordinate)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel state does not exist.");
            }
            meshkernel::Point const fromPoint{xFromCoordinate, yFromCoordinate};
            meshkernel::Point const toPoint{xToCoordinate, yToCoordinate};
            meshKernelState[meshKernelId].m_curvilinearGridLineShift->MoveNode(fromPoint, toPoint);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_line_shift(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel state does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGridLineShift == nullptr)
            {
                throw std::invalid_argument("MeshKernel: Curvilinear grid line shift algorithm instance is null.");
            }

            meshKernelState[meshKernelId].m_curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(meshKernelState[meshKernelId].m_curvilinearGridLineShift->Compute());
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_finalize_line_shift(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel state does not exist.");
            }

            meshKernelState[meshKernelId].m_curvilinearGridLineShift.reset();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_insert_face(int meshKernelId, double xCoordinate, double yCoordinate)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel state does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGrid == nullptr)
            {
                throw std::invalid_argument("MeshKernel: Empty curvilinear grid");
            }

            if (!meshKernelState[meshKernelId].m_curvilinearGrid->IsValid())
            {
                throw std::invalid_argument("MeshKernel: Not valid curvilinear grid.");
            }

            meshkernel::Point const point{xCoordinate, yCoordinate};

            meshKernelState[meshKernelId].m_curvilinearGrid->InsertFace(point);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_convert_to_mesh2d(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() > 0 && meshKernelState[meshKernelId].m_curvilinearGrid->m_projection != meshKernelState[meshKernelId].m_mesh2d->m_projection)
            {
                throw std::invalid_argument("MeshKernel: The existing mesh2d projection is not equal to the curvilinear grid projection");
            }

            const auto [nodes, edges, gridIndices] = meshKernelState[meshKernelId].m_curvilinearGrid->ConvertCurvilinearToNodesAndEdges();

            *meshKernelState[meshKernelId].m_mesh2d += meshkernel::Mesh2D(edges, nodes, meshKernelState[meshKernelId].m_curvilinearGrid->m_projection);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
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
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel state does not exist.");
            }

            meshkernel::CurvilinearGridLineAttractionRepulsion curvilinearLineAttractionRepulsion(meshKernelState[meshKernelId].m_curvilinearGrid, repulsionParameter);

            meshkernel::Point const lineFrom{xFirstNodeOnTheLine, yFirstNodeOnTheLine};
            meshkernel::Point const lineTo{xSecondNodeOnTheLine, ySecondNodeOnTheLine};
            curvilinearLineAttractionRepulsion.SetLine(lineFrom, lineTo);

            meshkernel::Point const lowerLeft{xLowerLeftCorner, yLowerLeftCorner};
            meshkernel::Point const upperRight{xUpperRightCorner, yUpperRightCorner};
            curvilinearLineAttractionRepulsion.SetBlock(lowerLeft, upperRight);

            *meshKernelState[meshKernelId].m_curvilinearGrid = curvilinearLineAttractionRepulsion.Compute();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_line_mirror(int meshKernelId,
                                                    double mirroring,
                                                    double xFirstGridLineNode,
                                                    double yFirstGridLineNode,
                                                    double xSecondGridLineNode,
                                                    double ySecondGridLineNode)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_curvilinearGrid == nullptr)
            {
                throw std::invalid_argument("MeshKernel: Not a valid curvilinear grid instance.");
            }

            if (!meshKernelState[meshKernelId].m_curvilinearGrid->IsValid())
            {
                throw std::invalid_argument("MeshKernel: Not valid curvilinear grid.");
            }

            auto curvilinearGridLineMirror = meshkernel::CurvilinearGridLineMirror(meshKernelState[meshKernelId].m_curvilinearGrid, mirroring);

            curvilinearGridLineMirror.SetLine({xFirstGridLineNode, yFirstGridLineNode}, {xSecondGridLineNode, ySecondGridLineNode});

            *meshKernelState[meshKernelId].m_curvilinearGrid = curvilinearGridLineMirror.Compute();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API double mkernel_get_separator()
    {
        return meshkernel::doubleMissingValue;
    }

    MKERNEL_API double mkernel_get_inner_outer_separator()
    {
        return meshkernel::innerOuterSeparator;
    }

    MKERNEL_API int mkernel_mesh2d_averaging_interpolation(int meshKernelId,
                                                           const GeometryList& samples,
                                                           int locationType,
                                                           int averagingMethodType,
                                                           double relativeSearchSize,
                                                           size_t minNumSamples,
                                                           GeometryList& results)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() == 0)
            {
                throw std::invalid_argument("MeshKernel: The mesh is empty.");
            }

            auto sampleValues = ConvertGeometryListToSampleVector(samples);
            auto const meshLocation = static_cast<meshkernel::MeshLocations>(locationType);
            auto const averagingMethod = static_cast<meshkernel::AveragingInterpolation::Method>(averagingMethodType);

            meshkernel::AveragingInterpolation averaging(meshKernelState[meshKernelId].m_mesh2d,
                                                         sampleValues,
                                                         averagingMethod,
                                                         meshLocation,
                                                         relativeSearchSize,
                                                         false,
                                                         false,
                                                         minNumSamples);
            // Execute averaging
            averaging.Compute();

            // Get the results and copy them to the result vector
            auto const& interpolationResults = averaging.GetResults();
            auto const locations = meshKernelState[meshKernelId].m_mesh2d->ComputeLocations(meshLocation);

            ConvertSampleVectorToGeometryList(locations, interpolationResults, results);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_mesh2d_triangulation_interpolation(int meshKernelId,
                                                               const GeometryList& samples,
                                                               int locationType,
                                                               GeometryList& results)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() == 0)
            {
                throw std::invalid_argument("MeshKernel: The mesh is empty.");
            }

            // Locations
            auto const sampleValues = ConvertGeometryListToSampleVector(samples);
            auto const meshLocation = static_cast<meshkernel::MeshLocations>(locationType);
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
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

} // namespace meshkernelapi
