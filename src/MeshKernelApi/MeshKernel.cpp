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
#include <MeshKernel/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGridCreateUniform.hpp>
#include <MeshKernel/CurvilinearGridDeRefinement.hpp>
#include <MeshKernel/CurvilinearGridFromPolygon.hpp>
#include <MeshKernel/CurvilinearGridFromSplines.hpp>
#include <MeshKernel/CurvilinearGridFromSplinesTransfinite.hpp>
#include <MeshKernel/CurvilinearGridOrthogonalization.hpp>
#include <MeshKernel/CurvilinearGridRefinement.hpp>
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
        meshkernel::Projection projection = isGeographic == 1 ? meshkernel::Projection::spherical : meshkernel::Projection::cartesian;
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

    MKERNEL_API int mkernel_delete_mesh2d(int meshKernelId, const GeometryList& polygon, int deletionOption, bool invertDeletion)
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

            const meshkernel::Polygons meshKernelPolygon(polygonPoints, meshKernelState[meshKernelId].m_mesh2d->m_projection);
            meshKernelState[meshKernelId].m_mesh2d->DeleteMesh(meshKernelPolygon, deletionOption, invertDeletion);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_set_mesh2d(int meshKernelId,
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

            *(meshKernelState[meshKernelId].m_mesh2d) = meshkernel::Mesh2D(edges2d,
                                                                           nodes2d,
                                                                           meshKernelState[meshKernelId].m_projection);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_set_mesh1d(int meshKernelId,
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
                                                                  mesh1d.nodex,
                                                                  mesh1d.nodey);

            *meshKernelState[meshKernelId].m_mesh1d = meshkernel::Mesh1D(edges1d,
                                                                         nodes1d,
                                                                         meshKernelState[meshKernelId].m_projection);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_mesh2d(int meshKernelId,
                                       Mesh2D& mesh2d)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            meshKernelState[meshKernelId].m_mesh2d->Administrate(meshkernel::Mesh2D::AdministrationOptions::AdministrateMeshEdges);
            meshKernelState[meshKernelId].m_mesh2d->SetFlatCopies();
            SetMesh(meshKernelState[meshKernelId].m_mesh2d, mesh2d);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_curvilinear(int meshKernelId,
                                            Mesh2D& mesh2d)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            meshKernelState[meshKernelId].m_curvilinearGrid->SetFlatCopies();

            // cast the curvilinear grid to mesh, because an unstructured mesh is communicated
            const auto mesh = std::dynamic_pointer_cast<meshkernel::Mesh>(meshKernelState[meshKernelId].m_curvilinearGrid);

            SetMesh(mesh, mesh2d);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_count_contacts(int meshKernelId,
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

    MKERNEL_API int mkernel_get_contacts(int meshKernelId,
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

    MKERNEL_API int mkernel_get_mesh1d(int meshKernelId,
                                       Mesh1D& mesh1d)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            meshKernelState[meshKernelId].m_mesh1d->SetFlatCopies();
            SetMesh1D(meshKernelState[meshKernelId].m_mesh1d, mesh1d);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_faces_mesh2d(int meshKernelId, Mesh2D& mesh2d)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            meshKernelState[meshKernelId].m_mesh2d->Administrate(meshkernel::Mesh2D::AdministrationOptions::AdministrateMeshEdgesAndFaces);
            meshKernelState[meshKernelId].m_mesh2d->SetFlatCopies();

            SetMesh(meshKernelState[meshKernelId].m_mesh2d, mesh2d);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_count_hanging_edges_mesh2d(int meshKernelId, int& numHangingEdges)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            const auto hangingEdges = meshKernelState[meshKernelId].m_mesh2d->GetHangingEdges();
            numHangingEdges = static_cast<int>(hangingEdges.size());
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_hanging_edges_mesh2d(int meshKernelId, int** hangingEdgesIndices)
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
                *(hangingEdgesIndices)[i] = static_cast<int>(hangingEdges[i]);
            }
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_delete_hanging_edges_mesh2d(int meshKernelId)
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

    MKERNEL_API int mkernel_compute_orthogonalization_mesh2d(int meshKernelId,
                                                             int projectToLandBoundaryOption,
                                                             const OrthogonalizationParameters& orthogonalizationParameters,
                                                             const GeometryList& polygons,
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

            // build enclosing polygon
            std::vector<meshkernel::Point> nodes(polygons.numberOfCoordinates);
            for (auto i = 0; i < polygons.numberOfCoordinates; i++)
            {
                nodes[i].x = polygons.xCoordinates[i];
                nodes[i].y = polygons.yCoordinates[i];
            }

            auto meshKernelPolygons = std::make_shared<meshkernel::Polygons>(nodes, meshKernelState[meshKernelId].m_mesh2d->m_projection);

            // build land boundary
            std::vector<meshkernel::Point> landBoundariesPoints(landBoundaries.numberOfCoordinates);
            for (auto i = 0; i < landBoundaries.numberOfCoordinates; i++)
            {
                landBoundariesPoints[i].x = landBoundaries.xCoordinates[i];
                landBoundariesPoints[i].y = landBoundaries.yCoordinates[i];
            }

            const auto orthogonalizer = std::make_shared<meshkernel::Orthogonalizer>(meshKernelState[meshKernelId].m_mesh2d);
            const auto smoother = std::make_shared<meshkernel::Smoother>(meshKernelState[meshKernelId].m_mesh2d);
            const auto meshKernelLandBoundary = std::make_shared<meshkernel::LandBoundaries>(landBoundariesPoints, meshKernelState[meshKernelId].m_mesh2d, meshKernelPolygons);

            meshkernel::OrthogonalizationAndSmoothing ortogonalization(meshKernelState[meshKernelId].m_mesh2d,
                                                                       smoother,
                                                                       orthogonalizer,
                                                                       meshKernelPolygons,
                                                                       meshKernelLandBoundary,
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

    MKERNEL_API int mkernel_initialize_orthogonalization_mesh2d(int meshKernelId,
                                                                int projectToLandBoundaryOption,
                                                                OrthogonalizationParameters& orthogonalizationParameters,
                                                                const GeometryList& geometryListPolygon,
                                                                const GeometryList& geometryListLandBoundaries)
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

            // build enclosing polygon
            std::vector<meshkernel::Point> nodes(geometryListPolygon.numberOfCoordinates);
            for (auto i = 0; i < geometryListPolygon.numberOfCoordinates; i++)
            {
                nodes[i].x = geometryListPolygon.xCoordinates[i];
                nodes[i].y = geometryListPolygon.yCoordinates[i];
            }

            // build land boundary
            std::vector<meshkernel::Point> landBoundaries(geometryListLandBoundaries.numberOfCoordinates);
            for (auto i = 0; i < geometryListLandBoundaries.numberOfCoordinates; i++)
            {
                landBoundaries[i].x = geometryListLandBoundaries.xCoordinates[i];
                landBoundaries[i].y = geometryListLandBoundaries.yCoordinates[i];
            }

            auto orthogonalizer = std::make_shared<meshkernel::Orthogonalizer>(meshKernelState[meshKernelId].m_mesh2d);
            auto smoother = std::make_shared<meshkernel::Smoother>(meshKernelState[meshKernelId].m_mesh2d);
            auto polygon = std::make_shared<meshkernel::Polygons>(nodes, meshKernelState[meshKernelId].m_mesh2d->m_projection);
            auto landBoundary = std::make_shared<meshkernel::LandBoundaries>(landBoundaries, meshKernelState[meshKernelId].m_mesh2d, polygon);

            meshKernelState[meshKernelId].m_orthogonalization = std::make_shared<meshkernel::OrthogonalizationAndSmoothing>(meshKernelState[meshKernelId].m_mesh2d,
                                                                                                                            smoother,
                                                                                                                            orthogonalizer,
                                                                                                                            polygon,
                                                                                                                            landBoundary,
                                                                                                                            static_cast<meshkernel::LandBoundaries::ProjectToLandBoundaryOption>(projectToLandBoundaryOption),
                                                                                                                            orthogonalizationParameters);
            meshKernelState[meshKernelId].m_orthogonalization->Initialize();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_prepare_outer_iteration_orthogonalization_mesh2d(int meshKernelId)
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

            meshKernelState[meshKernelId].m_orthogonalization->PrepareOuterIteration();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_compute_inner_ortogonalization_iteration_mesh2d(int meshKernelId)
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

            meshKernelState[meshKernelId].m_orthogonalization->InnerIteration();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_finalize_inner_ortogonalization_iteration_mesh2d(int meshKernelId)
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

            meshKernelState[meshKernelId].m_orthogonalization->FinalizeOuterIteration();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_delete_orthogonalization_mesh2d(int meshKernelId)
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

            meshKernelState[meshKernelId].m_orthogonalization.reset();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_orthogonality_mesh2d(int meshKernelId, GeometryList& geometryList)
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

            for (auto i = 0; i < geometryList.numberOfCoordinates; ++i)
            {
                geometryList.zCoordinates[i] = result[i];
            }
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_smoothness_mesh2d(int meshKernelId, GeometryList& geometryList)
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

            for (auto i = 0; i < geometryList.numberOfCoordinates; ++i)
            {
                geometryList.zCoordinates[i] = result[i];
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
            if (geometryListIn.numberOfCoordinates == 0)
            {
                throw std::invalid_argument("MeshKernel: The number of coordinates of the given geometry is zero.");
            }

            std::vector<meshkernel::Point> splines(geometryListIn.numberOfCoordinates);
            for (auto i = 0; i < geometryListIn.numberOfCoordinates; i++)
            {
                splines[i].x = geometryListIn.xCoordinates[i];
                splines[i].y = geometryListIn.yCoordinates[i];
            }

            const auto indices = FindIndices(splines, 0, splines.size(), meshkernel::doubleMissingValue);
            const auto numSplines = indices.size();

            int index = 0;
            for (auto s = 0; s < numSplines; s++)
            {
                std::vector<meshkernel::Point> coordinates(splines.begin() + indices[s][0], splines.begin() + int(indices[s][1]) + 1);
                const int numNodes = int(indices[s][1]) - int(indices[s][0]) + 1;
                const auto coordinatesDerivatives = meshkernel::Splines::SecondOrderDerivative(coordinates, 0, coordinates.size() - 1);

                for (auto n = 0; n < numNodes - 1; n++)
                {
                    for (auto p = 0; p <= numberOfPointsBetweenNodes; p++)
                    {
                        const double pointAdimensionalCoordinate = n + double(p) / double(numberOfPointsBetweenNodes);
                        const auto pointCoordinate = ComputePointOnSplineAtAdimensionalDistance(coordinates, coordinatesDerivatives, pointAdimensionalCoordinate);
                        if (!pointCoordinate.IsValid())
                        {
                            break;
                        }

                        geometryListOut.xCoordinates[index] = pointCoordinate.x;
                        geometryListOut.yCoordinates[index] = pointCoordinate.y;
                        geometryListOut.zCoordinates[index] = meshkernel::doubleMissingValue;
                        index++;
                    }
                }

                geometryListOut.xCoordinates[index] = meshkernel::doubleMissingValue;
                geometryListOut.yCoordinates[index] = meshkernel::doubleMissingValue;
                geometryListOut.zCoordinates[index] = meshkernel::doubleMissingValue;
                index++;
            }

            geometryListOut.numberOfCoordinates = index - 1;
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_make_mesh_from_polygon_mesh2d(int meshKernelId, const GeometryList& disposableGeometryListIn)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            auto result = ConvertGeometryListToPointVector(disposableGeometryListIn);

            const meshkernel::Polygons polygon(result, meshKernelState[meshKernelId].m_mesh2d->m_projection);

            // generate samples in all polygons
            const auto generatedPoints = polygon.ComputePointsInPolygons();

            const meshkernel::Mesh2D mesh(generatedPoints[0], polygon, meshKernelState[meshKernelId].m_mesh2d->m_projection);
            *meshKernelState[meshKernelId].m_mesh2d += mesh;
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_make_mesh_from_samples_mesh2d(int meshKernelId, const GeometryList& geometryList)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            auto samplePoints = ConvertGeometryListToPointVector(geometryList);

            meshkernel::Polygons polygon;
            const meshkernel::Mesh2D mesh(samplePoints, polygon, meshKernelState[meshKernelId].m_mesh2d->m_projection);
            *meshKernelState[meshKernelId].m_mesh2d += mesh;
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_mesh_boundaries_to_polygon_mesh2d(int meshKernelId, GeometryList& geometryList)
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

            ConvertPointVectorToGeometryList(meshBoundaryPolygon, geometryList);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_count_mesh_boundaries_to_polygon_mesh2d(int meshKernelId, int& numberOfPolygonNodes)
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

    MKERNEL_API int mkernel_refine_polygon(int meshKernelId, const GeometryList& geometryListIn, int firstIndex, int secondIndex, double distance, GeometryList& geometryListOut)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            auto polygonPoints = ConvertGeometryListToPointVector(geometryListIn);

            const meshkernel::Polygons polygon(polygonPoints, meshKernelState[meshKernelId].m_mesh2d->m_projection);
            const auto refinedPolygon = polygon.RefineFirstPolygon(firstIndex, secondIndex, distance);

            ConvertPointVectorToGeometryList(refinedPolygon, geometryListOut);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_count_refine_polygon(int meshKernelId,
                                                 const GeometryList& geometryListIn,
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

            auto polygonPoints = ConvertGeometryListToPointVector(geometryListIn);

            const meshkernel::Polygons polygon(polygonPoints, meshKernelState[meshKernelId].m_mesh2d->m_projection);

            const auto refinedPolygon = polygon.RefineFirstPolygon(firstIndex, secondIndex, distance);

            numberOfPolygonNodes = int(refinedPolygon.size());
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_merge_nodes_mesh2d(int meshKernelId, const GeometryList& geometryListIn)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            auto polygonPoints = ConvertGeometryListToPointVector(geometryListIn);

            const meshkernel::Polygons polygon(polygonPoints, meshKernelState[meshKernelId].m_mesh2d->m_projection);

            meshKernelState[meshKernelId].m_mesh2d->MergeNodesInPolygon(polygon);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_merge_two_nodes_mesh2d(int meshKernelId, int startNode, int endNode)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            meshKernelState[meshKernelId].m_mesh2d->MergeTwoNodes(startNode, endNode);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_nodes_in_polygons(int meshKernelId,
                                              const GeometryList& geometryListIn,
                                              int inside,
                                              int** selectedNodes)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            auto polygonPoints = ConvertGeometryListToPointVector(geometryListIn);

            const meshkernel::Polygons polygon(polygonPoints, meshKernelState[meshKernelId].m_mesh2d->m_projection);

            const bool selectInside = inside == 1 ? true : false;
            meshKernelState[meshKernelId].m_mesh2d->MaskNodesInPolygons(polygon, selectInside);

            int index = 0;
            for (auto i = 0; i < meshKernelState[meshKernelId].m_mesh2d->GetNumNodes(); ++i)
            {
                if (meshKernelState[meshKernelId].m_mesh2d->m_nodeMask[i] > 0)
                {
                    (*selectedNodes)[index] = i;
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

    MKERNEL_API int mkernel_count_nodes_in_polygons(int meshKernelId,
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
            auto polygonPoints = ConvertGeometryListToPointVector(geometryListIn);

            const meshkernel::Polygons polygon(polygonPoints, meshKernelState[meshKernelId].m_mesh2d->m_projection);

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

    MKERNEL_API int mkernel_insert_edge_mesh2d(int meshKernelId, int startNode, int endNode, int& new_edge_index)
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

    MKERNEL_API int mkernel_insert_node_mesh2d(int meshKernelId, double xCoordinate, double yCoordinate, int& nodeIndex)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                //create a valid instance, by default cartesian
                *meshKernelState[meshKernelId].m_mesh2d = meshkernel::Mesh2D();
                meshKernelState[meshKernelId].m_mesh2d->m_projection = meshkernel::Projection::cartesian;
            }

            const meshkernel::Point newNode{xCoordinate, yCoordinate};
            nodeIndex = static_cast<int>(meshKernelState[meshKernelId].m_mesh2d->InsertNode(newNode));
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_delete_node_mesh2d(int meshKernelId, int nodeIndex)
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

    MKERNEL_API int mkernel_move_node_mesh2d(int meshKernelId, const GeometryList& geometryListIn, int nodeIndex)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            auto newPoint = ConvertGeometryListToPointVector(geometryListIn);

            meshKernelState[meshKernelId].m_mesh2d->MoveNode(newPoint[0], nodeIndex);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_delete_edge_mesh2d(int meshKernelId, const GeometryList& geometryListIn)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            auto newPoint = ConvertGeometryListToPointVector(geometryListIn);

            const auto edgeIndex = meshKernelState[meshKernelId].m_mesh2d->FindEdgeCloseToAPoint(newPoint[0]);

            meshKernelState[meshKernelId].m_mesh2d->DeleteEdge(edgeIndex);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_find_edge_mesh2d(int meshKernelId, const GeometryList& geometryListIn, int& edgeIndex)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            auto newPoint = ConvertGeometryListToPointVector(geometryListIn);

            edgeIndex = static_cast<int>(meshKernelState[meshKernelId].m_mesh2d->FindEdgeCloseToAPoint(newPoint[0]));
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_offsetted_polygon(int meshKernelId, const GeometryList& geometryListIn, bool innerPolygon, double distance, GeometryList& geometryListOut)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            auto polygonPoints = ConvertGeometryListToPointVector(geometryListIn);

            const meshkernel::Polygons polygon(polygonPoints, meshKernelState[meshKernelId].m_mesh2d->m_projection);

            const auto newPolygon = polygon.OffsetCopy(distance, innerPolygon);

            ConvertPointVectorToGeometryList(newPolygon.m_nodes, geometryListOut);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_count_offsetted_polygon(int meshKernelId, const GeometryList& geometryListIn, bool innerPolygon, double distance, int& numberOfPolygonNodes)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            auto polygonPoints = ConvertGeometryListToPointVector(geometryListIn);

            const meshkernel::Polygons polygon(polygonPoints, meshKernelState[meshKernelId].m_mesh2d->m_projection);
            const auto newPolygon = polygon.OffsetCopy(distance, innerPolygon);

            numberOfPolygonNodes = static_cast<int>(newPolygon.GetNumNodes());
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_refine_based_on_samples_mesh2d(int meshKernelId,
                                                           const GeometryList& geometryListIn,
                                                           const InterpolationParameters& interpolationParameters,
                                                           const SampleRefineParameters& sampleRefineParameters)
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

            auto samples = ConvertGeometryListToSampleVector(geometryListIn);

            meshkernel::AveragingInterpolation::Method averagingMethod;
            if (sampleRefineParameters.RefinementType == 2)
            {
                averagingMethod = meshkernel::AveragingInterpolation::Method::MinAbsValue;
            }
            if (sampleRefineParameters.RefinementType == 3)
            {

                averagingMethod = meshkernel::AveragingInterpolation::Method::Max;
            }

            const bool refineOutsideFace = sampleRefineParameters.AccountForSamplesOutside == 1 ? true : false;
            const bool transformSamples = sampleRefineParameters.RefinementType == 3 ? true : false;

            const auto averaging = std::make_shared<meshkernel::AveragingInterpolation>(meshKernelState[meshKernelId].m_mesh2d,
                                                                                        samples,
                                                                                        averagingMethod,
                                                                                        meshkernel::MeshLocations::Faces,
                                                                                        1.0,
                                                                                        refineOutsideFace,
                                                                                        transformSamples);

            meshkernel::MeshRefinement meshRefinement(meshKernelState[meshKernelId].m_mesh2d, averaging, sampleRefineParameters, interpolationParameters);
            meshRefinement.Compute();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_refine_based_on_polygon_mesh2d(int meshKernelId, const GeometryList& geometryList, const InterpolationParameters& interpolationParameters)
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

            meshkernel::MeshRefinement meshRefinement(meshKernelState[meshKernelId].m_mesh2d, polygon, interpolationParameters);
            meshRefinement.Compute();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_node_index_mesh2d(int meshKernelId, const GeometryList& geometryListIn, double searchRadius, int& nodeIndex)
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

            auto polygonPoints = ConvertGeometryListToPointVector(geometryListIn);

            nodeIndex = static_cast<int>(meshKernelState[meshKernelId].m_mesh2d->FindNodeCloseToAPoint(polygonPoints[0], searchRadius));
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_closest_node_mesh2d(int meshKernelId, const GeometryList& geometryListIn, double searchRadius, GeometryList& geometryListOut)
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

            if (geometryListOut.numberOfCoordinates <= 0)
            {
                throw std::invalid_argument("MeshKernel: The output-geometry has no coordinates.");
            }

            auto polygonPoints = ConvertGeometryListToPointVector(geometryListIn);

            const auto nodeIndex = meshKernelState[meshKernelId].m_mesh2d->FindNodeCloseToAPoint(polygonPoints[0], searchRadius);

            // Set the node coordinate
            const auto node = meshKernelState[meshKernelId].m_mesh2d->m_nodes[nodeIndex];
            std::vector<meshkernel::Point> pointVector;
            pointVector.emplace_back(node);
            ConvertPointVectorToGeometryList(pointVector, geometryListOut);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_points_in_polygon(int meshKernelId, const GeometryList& polygon, const GeometryList& pointsNative, GeometryList& selectedPointsNative)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            auto polygonNodes = ConvertGeometryListToPointVector(polygon);

            auto points = ConvertGeometryListToPointVector(pointsNative);
            const meshkernel::Polygons localPolygon(polygonNodes, meshKernelState[meshKernelId].m_mesh2d->m_projection);

            for (auto i = 0; i < points.size(); i++)
            {
                selectedPointsNative.zCoordinates[i] = localPolygon.IsPointInPolygon(points[i], 0) ? 1.0 : 0.0;
            }
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_flip_edges_mesh2d(int meshKernelId,
                                              int isTriangulationRequired,
                                              int projectToLandBoundaryRequired)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            //set landboundaries
            auto polygon = std::make_shared<meshkernel::Polygons>();

            std::vector<meshkernel::Point> landBoundary;
            const auto landBoundaries = std::make_shared<meshkernel::LandBoundaries>(landBoundary, meshKernelState[meshKernelId].m_mesh2d, polygon);

            const bool triangulateFaces = isTriangulationRequired == 0 ? false : true;
            const bool projectToLandBoundary = projectToLandBoundaryRequired == 0 ? false : true;
            const meshkernel::FlipEdges flipEdges(meshKernelState[meshKernelId].m_mesh2d, landBoundaries, triangulateFaces, projectToLandBoundary);

            flipEdges.Compute();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_count_small_flow_edge_centers_mesh2d(int meshKernelId, double smallFlowEdgesThreshold, int& numSmallFlowEdges)
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

            numSmallFlowEdges = static_cast<int>(smallFlowEdgeCenters.size());
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_small_flow_edge_centers_mesh2d(int meshKernelId, double smallFlowEdgesThreshold, GeometryList& result)
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

    MKERNEL_API int mkernel_count_obtuse_triangles_mesh2d(int meshKernelId, int& numObtuseTriangles)
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

    MKERNEL_API int mkernel_get_obtuse_triangles_mass_centers_mesh2d(int meshKernelId, GeometryList& result)
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

    MKERNEL_API int mkernel_delete_small_flow_edges_mesh2d(int meshKernelId, double smallFlowEdgesThreshold, double minFractionalAreaTriangles)
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

    MKERNEL_API int mkernel_compute_single_contacts(int meshKernelId,
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

    MKERNEL_API int mkernel_compute_multiple_contacts(int meshKernelId,
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

    MKERNEL_API int mkernel_compute_with_polygons_contacts(int meshKernelId,
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
            meshKernelState[meshKernelId].m_contacts->ComputeContactsWithPolygons(meshKernel1DNodeMask,
                                                                                  meshKernelPolygons);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }
    MKERNEL_API int mkernel_compute_with_points_contacts(int meshKernelId,
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

            // Convert 1D node mask from int** to vector<bool>
            auto num1DNodes = meshKernelState[meshKernelId].m_mesh1d->GetNumNodes();
            auto meshKernel1DNodeMask = ConvertIntegerArrayToBoolVector(oneDNodeMask,
                                                                        num1DNodes);

            // Convert polygon date from GeometryList to Point vector
            auto meshKernelPoints = ConvertGeometryListToPointVector(points);
            // Execute
            meshKernelState[meshKernelId].m_contacts->ComputeContactsWithPoints(meshKernel1DNodeMask,
                                                                                meshKernelPoints);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_compute_boundary_contacts(int meshKernelId,
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

            // Convert 1D node mask from int** to vector<bool>
            auto num1DNodes = meshKernelState[meshKernelId].m_mesh1d->GetNumNodes();
            auto meshKernel1DNodeMask = ConvertIntegerArrayToBoolVector(oneDNodeMask,
                                                                        num1DNodes);

            // Convert polygon date from GeometryList to Polygons
            auto polygonPoints = ConvertGeometryListToPointVector(polygons);
            const meshkernel::Polygons meshKernelPolygons(polygonPoints,
                                                          meshKernelState[meshKernelId].m_mesh2d->m_projection);

            // Execute
            meshKernelState[meshKernelId].m_contacts->ComputeBoundaryContacts(meshKernel1DNodeMask,
                                                                              meshKernelPolygons,
                                                                              searchRadius);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_refine_curvilinear(int meshKernelId, const GeometryList& geometryListFirstPoint, const GeometryList& geometryListSecondPoint, int refinement)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            const auto firstPoint = ConvertGeometryListToPointVector(geometryListFirstPoint);

            if (firstPoint.empty())
            {
                throw std::invalid_argument("MeshKernel: No first node of the segment defining the refinement zone has been provided.");
            }

            const auto secondPoint = ConvertGeometryListToPointVector(geometryListSecondPoint);
            if (secondPoint.empty())
            {
                throw std::invalid_argument("MeshKernel: No second node of the segment defining the refinement zone has been provided.");
            }

            // Execute
            meshkernel::CurvilinearGridRefinement curvilinearGridRefinement(meshKernelState[meshKernelId].m_curvilinearGrid, firstPoint[0], secondPoint[0], refinement);
            curvilinearGridRefinement.Compute();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_derefine_curvilinear(int meshKernelId,
                                                 const GeometryList& geometryListFirstPoint,
                                                 const GeometryList& geometryListSecondPoint)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            const auto firstPoint = ConvertGeometryListToPointVector(geometryListFirstPoint);

            if (firstPoint.empty())
            {
                throw std::invalid_argument("MeshKernel: No first node of the segment defining the refinement zone has been provided.");
            }

            const auto secondPoint = ConvertGeometryListToPointVector(geometryListSecondPoint);
            if (secondPoint.empty())
            {
                throw std::invalid_argument("MeshKernel: No second node of the segment defining the refinement zone has been provided.");
            }

            // Execute
            meshkernel::CurvilinearGridDeRefinement curvilinearGridDeRefinement(meshKernelState[meshKernelId].m_curvilinearGrid, firstPoint[0], secondPoint[0]);

            meshKernelState[meshKernelId].m_curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(curvilinearGridDeRefinement.Compute());
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_compute_transfinite_from_splines_curvilinear(int meshKernelId,
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

    MKERNEL_API int mkernel_compute_transfinite_from_polygon_curvilinear(int meshKernelId,
                                                                         const GeometryList& polygons,
                                                                         int firstNode,
                                                                         int secondNode,
                                                                         int thirdNode,
                                                                         bool useFourthSide)
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

            const auto curvilinearGrid = curvilinearGridFromPolygon.Compute(firstNode, secondNode, thirdNode, useFourthSide);

            // set the curvilinear state
            meshKernelState[meshKernelId].m_curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(curvilinearGrid);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_compute_transfinite_from_triangle_curvilinear(int meshKernelId,
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

    MKERNEL_API int mkernel_compute_orthogonal_curvilinear(int meshKernelId,
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

    MKERNEL_API int mkernel_initialize_orthogonal_curvilinear(int meshKernelId,
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

    MKERNEL_API int mkernel_iterate_orthogonal_curvilinear(int meshKernelId, int layer)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }

            meshKernelState[meshKernelId].m_curvilinearGridFromSplines->Iterate(layer);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_refresh_orthogonal_curvilinear(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
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

    MKERNEL_API int mkernel_delete_orthogonal_curvilinear(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel id does not exist.");
            }
            meshKernelState[meshKernelId].m_curvilinearGridFromSplines.reset();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_make_uniform_curvilinear(int meshKernelId,
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

    MKERNEL_API int mkernel_orthogonalize_curvilinear(int meshKernelId,
                                                      const OrthogonalizationParameters& orthogonalizationParameters,
                                                      const GeometryList& geometryListFirstPoint,
                                                      const GeometryList& geometryListSecondPoint)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh kernel state does not exist.");
            }
            const auto firstPoint = ConvertGeometryListToPointVector(geometryListFirstPoint);

            if (firstPoint.empty())
            {
                throw std::invalid_argument("MeshKernel: No first node of the segment defining the refinement zone has been provided.");
            }

            const auto secondPoint = ConvertGeometryListToPointVector(geometryListSecondPoint);
            if (secondPoint.empty())
            {
                throw std::invalid_argument("MeshKernel: No second node of the segment defining the refinement zone has been provided.");
            }

            // Execute
            meshkernel::CurvilinearGridOrthogonalization curvilinearGridOrthogonalization(meshKernelState[meshKernelId].m_curvilinearGrid,
                                                                                          orthogonalizationParameters,
                                                                                          firstPoint[0],
                                                                                          secondPoint[0]);

            curvilinearGridOrthogonalization.Compute();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_convert_curvilinear_to_mesh2d(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelState.count(meshKernelId) == 0)
            {
                throw std::invalid_argument("mkernel_convert_curvilinear_to_mesh2d: The selected mesh kernel id does not exist.");
            }

            if (meshKernelState[meshKernelId].m_mesh2d->GetNumNodes() > 0 && meshKernelState[meshKernelId].m_curvilinearGrid->m_projection != meshKernelState[meshKernelId].m_mesh2d->m_projection)
            {
                throw std::invalid_argument("mkernel_convert_curvilinear_to_mesh2d: The existing mesh2d projection is not equal to the curvilinear grid projection");
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

    MKERNEL_API double mkernel_get_separator()
    {
        return meshkernel::doubleMissingValue;
    }

    MKERNEL_API double mkernel_get_inner_outer_separator()
    {
        return meshkernel::innerOuterSeparator;
    }

    MKERNEL_API int averaging(const Mesh2D& mesh2d,
                              const int& startIndex,
                              const double** samplesXCoordinate,
                              const double** samplesYCoordinate,
                              const double** samplesValue,
                              const int& numSamples,
                              double** results,
                              const int& locationType,
                              const double& Wu1Duni,
                              const int& averagingMethod,
                              const int& minNumberOfSamples,
                              const double& relativeSearchSize,
                              const int& spherical,
                              const int& sphericalAccurate)
    {
        int exitCode = Success;
        try
        {
            // Projection
            auto projection = meshkernel::Projection::cartesian;
            if (spherical == 1)
            {
                projection = meshkernel::Projection::spherical;
            }
            if (sphericalAccurate == 1)
            {
                projection = meshkernel::Projection::sphericalAccurate;
            }

            // Set the mesh
            const auto edges = meshkernel::ConvertToEdgeNodesVector(mesh2d.num_edges, mesh2d.edge_nodes);
            const auto nodes = meshkernel::ConvertToNodesVector(mesh2d.num_nodes, mesh2d.node_x, mesh2d.node_y);
            const auto mesh = std::make_shared<meshkernel::Mesh2D>(edges, nodes, projection);

            // Build the samples
            std::vector<meshkernel::Sample> samples(numSamples);
            for (auto i = 0; i < samples.size(); ++i)
            {
                samples[i].x = (*samplesXCoordinate)[i];
                samples[i].y = (*samplesYCoordinate)[i];
                samples[i].value = (*samplesValue)[i];
            }

            // Execute averaging
            meshkernel::AveragingInterpolation averaging(mesh,
                                                         samples,
                                                         static_cast<meshkernel::AveragingInterpolation::Method>(averagingMethod),
                                                         static_cast<meshkernel::MeshLocations>(locationType),
                                                         relativeSearchSize,
                                                         false,
                                                         false);
            averaging.Compute();

            // Get the results and copy them to the result vector
            auto interpolationResults = averaging.GetResults();
            for (auto i = 0; i < interpolationResults.size(); ++i)
            {
                (*results)[i] = interpolationResults[i];
            }
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    // ec_module dll (stateless)
    MKERNEL_API int triangulation(const Mesh2D& mesh2d,
                                  const double** samplesXCoordinate,
                                  const double** samplesYCoordinate,
                                  const double** samplesValue,
                                  const int& numSamples,
                                  const int& locationType,
                                  const int& spherical,
                                  const int& sphericalAccurate,
                                  double** results)
    {
        int exitCode = Success;
        try
        {
            // Projection
            auto projection = meshkernel::Projection::cartesian;
            if (spherical == 1)
            {
                projection = meshkernel::Projection::spherical;
            }
            if (sphericalAccurate == 1)
            {
                projection = meshkernel::Projection::sphericalAccurate;
            }

            // Locations
            const auto location = static_cast<meshkernel::MeshLocations>(locationType);
            const auto locations = ComputeLocations(mesh2d, location);

            // Build the samples
            const auto samples = meshkernel::Sample::ConvertToSamples(numSamples, samplesXCoordinate, samplesYCoordinate, samplesValue);

            // Execute triangulation
            meshkernel::TriangulationInterpolation triangulationInterpolation(locations, samples, projection);
            triangulationInterpolation.Compute();

            // Get the results and copy them back to the results vector
            auto interpolationResults = triangulationInterpolation.GetResults();
            for (auto i = 0; i < interpolationResults.size(); ++i)
            {
                (*results)[i] = interpolationResults[i];
            }
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

} // namespace meshkernelapi
