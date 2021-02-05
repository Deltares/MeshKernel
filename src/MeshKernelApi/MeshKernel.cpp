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
#include <vector>

#include <MeshKernel/AveragingInterpolation.hpp>
#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Contacts.hpp>
#include <MeshKernel/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGridFromPolygon.hpp>
#include <MeshKernel/CurvilinearGridFromSplines.hpp>
#include <MeshKernel/CurvilinearGridFromSplinesTransfinite.hpp>
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
#include <MeshKernelApi/Utils.hpp>

namespace meshkernelapi
{
    // The vectors containing the mesh instances
    static std::vector<std::shared_ptr<meshkernel::Mesh2D>> mesh2dInstances;
    static std::vector<std::shared_ptr<meshkernel::Mesh1D>> mesh1dInstances;
    static std::vector<std::shared_ptr<meshkernel::Contacts>> contactsInstances;

    // For interactivity
    static std::map<int, std::shared_ptr<meshkernel::OrthogonalizationAndSmoothing>> orthogonalizationInstances;
    static std::map<int, std::shared_ptr<meshkernel::CurvilinearGridFromSplines>> curvilinearGridFromSplinesInstances;

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
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            return InvalidGeometry;
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            return Exception;
        }
    }

    MKERNEL_API int mkernel_new_mesh(int& meshKernelId)
    {
        meshKernelId = int(mesh2dInstances.size());
        mesh2dInstances.emplace_back(std::make_shared<meshkernel::Mesh2D>());
        mesh1dInstances.emplace_back(std::make_shared<meshkernel::Mesh1D>());
        return Success;
    };

    MKERNEL_API int mkernel_deallocate_state(int meshKernelId)
    {
        if (meshKernelId < 0 && meshKernelId >= mesh2dInstances.size())
        {
            return -1;
        }

        mesh2dInstances.erase(mesh2dInstances.begin() + meshKernelId);
        mesh1dInstances.erase(mesh1dInstances.begin() + meshKernelId);
        return 0;
    }

    MKERNEL_API int mkernel_delete_mesh(int meshKernelId, const GeometryList& polygon, int deletionOption, bool invertDeletion)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            if (mesh2dInstances[meshKernelId]->GetNumNodes() <= 0)
            {
                return exitCode;
            }

            auto polygonPoints = ConvertGeometryListToPointVector(polygon);

            const meshkernel::Polygons polygon(polygonPoints, mesh2dInstances[meshKernelId]->m_projection);
            mesh2dInstances[meshKernelId]->DeleteMesh(polygon, deletionOption, invertDeletion);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_set_state(int meshKernelId,
                                      const MeshGeometryDimensions& meshGeometryDimensions,
                                      const MeshGeometry& meshGeometry,
                                      const Mesh1D& mesh1d,
                                      bool isGeographic)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            // convert raw arrays to containers
            const auto edges1d = meshkernel::ConvertToEdgeNodesVector(mesh1d.num_edges,
                                                                      mesh1d.edge_nodes);
            const auto nodes1d = meshkernel::ConvertToNodesVector(mesh1d.num_nodes,
                                                                  mesh1d.nodex,
                                                                  mesh1d.nodey);

            const auto edges2d = meshkernel::ConvertToEdgeNodesVector(meshGeometryDimensions.numedge,
                                                                      meshGeometry.edge_nodes);
            const auto nodes2d = meshkernel::ConvertToNodesVector(meshGeometryDimensions.numnode,
                                                                  meshGeometry.nodex,
                                                                  meshGeometry.nodey);

            // spherical or cartesian
            meshkernel::Projection projection = isGeographic ? meshkernel::Projection::spherical : meshkernel::Projection::cartesian;

            mesh1dInstances[meshKernelId] = std::make_shared<meshkernel::Mesh1D>(edges1d,
                                                                                 nodes1d,
                                                                                 projection);

            mesh2dInstances[meshKernelId] = std::make_shared<meshkernel::Mesh2D>(edges2d,
                                                                                 nodes2d,
                                                                                 projection);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_mesh(int meshKernelId,
                                     MeshGeometryDimensions& meshGeometryDimensions,
                                     MeshGeometry& meshGeometry,
                                     Mesh1D& mesh1d)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            mesh2dInstances[meshKernelId]->SetFlatCopies(meshkernel::Mesh2D::AdministrationOptions::AdministrateMeshEdges);
            mesh1dInstances[meshKernelId]->SetFlatCopies();

            SetMesh2DGeometry(mesh2dInstances, meshKernelId, meshGeometryDimensions, meshGeometry);
            SetMesh1DGeometry(mesh1dInstances, meshKernelId, mesh1d);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_find_faces(int meshKernelId, MeshGeometryDimensions& meshGeometryDimensions, MeshGeometry& meshGeometry)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            mesh2dInstances[meshKernelId]->SetFlatCopies(meshkernel::Mesh2D::AdministrationOptions::AdministrateMeshEdgesAndFaces);

            SetMesh2DGeometry(mesh2dInstances, meshKernelId, meshGeometryDimensions, meshGeometry);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_count_hanging_edges(int meshKernelId, int& numHangingEdges)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            const auto hangingEdges = mesh2dInstances[meshKernelId]->GetHangingEdges();
            numHangingEdges = hangingEdges.size();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_hanging_edges(int meshKernelId, int** hangingEdgesIndices)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            const auto hangingEdges = mesh2dInstances[meshKernelId]->GetHangingEdges();
            for (auto i = 0; i < hangingEdges.size(); ++i)
            {
                *(hangingEdgesIndices)[i] = hangingEdges[i];
            }
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_delete_hanging_edges(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            mesh2dInstances[meshKernelId]->DeleteHangingEdges();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_orthogonalize(int meshKernelId,
                                          int projectToLandBoundaryOption,
                                          const OrthogonalizationParameters& orthogonalizationParameters,
                                          const GeometryList& polygons,
                                          const GeometryList& landBoundaries)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            if (mesh2dInstances[meshKernelId]->GetNumNodes() <= 0)
            {
                return exitCode;
            }

            // build enclosing polygon
            std::vector<meshkernel::Point> nodes(polygons.numberOfCoordinates);
            for (auto i = 0; i < polygons.numberOfCoordinates; i++)
            {
                nodes[i].x = polygons.xCoordinates[i];
                nodes[i].y = polygons.yCoordinates[i];
            }

            auto meshKernelPolygons = std::make_shared<meshkernel::Polygons>(nodes, mesh2dInstances[meshKernelId]->m_projection);

            // build land boundary
            std::vector<meshkernel::Point> landBoundariesPoints(landBoundaries.numberOfCoordinates);
            for (auto i = 0; i < landBoundaries.numberOfCoordinates; i++)
            {
                landBoundariesPoints[i].x = landBoundaries.xCoordinates[i];
                landBoundariesPoints[i].y = landBoundaries.yCoordinates[i];
            }

            const auto orthogonalizer = std::make_shared<meshkernel::Orthogonalizer>(mesh2dInstances[meshKernelId]);
            const auto smoother = std::make_shared<meshkernel::Smoother>(mesh2dInstances[meshKernelId]);
            const auto meshKernelLandBoundary = std::make_shared<meshkernel::LandBoundaries>(landBoundariesPoints, mesh2dInstances[meshKernelId], meshKernelPolygons);

            meshkernel::OrthogonalizationAndSmoothing ortogonalization(mesh2dInstances[meshKernelId],
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

    MKERNEL_API int mkernel_orthogonalize_initialize(int meshKernelId,
                                                     int projectToLandBoundaryOption,
                                                     OrthogonalizationParameters& orthogonalizationParameters,
                                                     const GeometryList& geometryListPolygon,
                                                     const GeometryList& geometryListLandBoundaries)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            if (mesh2dInstances[meshKernelId]->GetNumNodes() <= 0)
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

            auto orthogonalizer = std::make_shared<meshkernel::Orthogonalizer>(mesh2dInstances[meshKernelId]);
            auto smoother = std::make_shared<meshkernel::Smoother>(mesh2dInstances[meshKernelId]);
            auto polygon = std::make_shared<meshkernel::Polygons>(nodes, mesh2dInstances[meshKernelId]->m_projection);
            auto landBoundary = std::make_shared<meshkernel::LandBoundaries>(landBoundaries, mesh2dInstances[meshKernelId], polygon);

            auto orthogonalizationInstance = std::make_shared<meshkernel::OrthogonalizationAndSmoothing>(mesh2dInstances[meshKernelId],
                                                                                                         smoother,
                                                                                                         orthogonalizer,
                                                                                                         polygon,
                                                                                                         landBoundary,
                                                                                                         static_cast<meshkernel::LandBoundaries::ProjectToLandBoundaryOption>(projectToLandBoundaryOption),
                                                                                                         orthogonalizationParameters);
            orthogonalizationInstance->Initialize();

            orthogonalizationInstances.insert({meshKernelId, orthogonalizationInstance});
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_orthogonalize_prepare_outer_iteration(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            if (mesh2dInstances[meshKernelId]->GetNumNodes() <= 0)
            {
                return exitCode;
            }

            orthogonalizationInstances[meshKernelId]->PrepareOuterIteration();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_orthogonalize_inner_iteration(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            if (mesh2dInstances[meshKernelId]->GetNumNodes() <= 0)
            {
                return exitCode;
            }

            orthogonalizationInstances[meshKernelId]->InnerIteration();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_orthogonalize_finalize_outer_iteration(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            if (mesh2dInstances[meshKernelId]->GetNumNodes() <= 0)
            {
                return exitCode;
            }

            orthogonalizationInstances[meshKernelId]->FinalizeOuterIteration();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_orthogonalize_delete(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            if (mesh2dInstances[meshKernelId]->GetNumNodes() <= 0)
            {
                return exitCode;
            }

            orthogonalizationInstances.erase(meshKernelId);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_orthogonality(int meshKernelId, GeometryList& geometryList)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            if (mesh2dInstances[meshKernelId]->GetNumNodes() <= 0)
            {
                return exitCode;
            }

            const auto result = mesh2dInstances[meshKernelId]->GetOrthogonality();

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

    MKERNEL_API int mkernel_get_smoothness(int meshKernelId, GeometryList& geometryList)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            if (mesh2dInstances[meshKernelId]->GetNumNodes() <= 0)
            {
                return exitCode;
            }

            const auto result = mesh2dInstances[meshKernelId]->GetSmoothness();

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
            std::vector<meshkernel::Point> coordinatesDerivatives(geometryListIn.numberOfCoordinates);

            int index = 0;
            for (auto s = 0; s < numSplines; s++)
            {
                std::vector<meshkernel::Point> coordinates(splines.begin() + indices[s][0], splines.begin() + int(indices[s][1]) + 1);
                const int numNodes = int(indices[s][1]) - int(indices[s][0]) + 1;
                meshkernel::Splines::SecondOrderDerivative(coordinates, numNodes, coordinatesDerivatives);

                for (auto n = 0; n < numNodes - 1; n++)
                {
                    for (auto p = 0; p <= numberOfPointsBetweenNodes; p++)
                    {
                        const double pointAdimensionalCoordinate = n + double(p) / double(numberOfPointsBetweenNodes);
                        meshkernel::Point pointCoordinate{meshkernel::doubleMissingValue, meshkernel::doubleMissingValue};
                        const bool successful = InterpolateSplinePoint(coordinates, coordinatesDerivatives, pointAdimensionalCoordinate, pointCoordinate);
                        if (!successful)
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

    MKERNEL_API int mkernel_make_mesh(int meshKernelId, const MakeMeshParameters& makeGridParameters, const GeometryList& geometryList)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            auto result = ConvertGeometryListToPointVector(geometryList);

            const meshkernel::Polygons polygon(result, mesh2dInstances[meshKernelId]->m_projection);

            meshkernel::Mesh2D mesh;
            mesh.MakeMesh(makeGridParameters, polygon);

            *mesh2dInstances[meshKernelId] += mesh;
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_make_mesh_from_polygon(int meshKernelId, const GeometryList& disposableGeometryListIn)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            auto result = ConvertGeometryListToPointVector(disposableGeometryListIn);

            const meshkernel::Polygons polygon(result, mesh2dInstances[meshKernelId]->m_projection);

            // generate samples in all polygons
            const auto generatedPoints = polygon.ComputePointsInPolygons();

            const meshkernel::Mesh2D mesh(generatedPoints[0], polygon, mesh2dInstances[meshKernelId]->m_projection);
            *mesh2dInstances[meshKernelId] += mesh;
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_make_mesh_from_samples(int meshKernelId, GeometryList& geometryList)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            auto samplePoints = ConvertGeometryListToPointVector(geometryList);

            meshkernel::Polygons polygon;
            const meshkernel::Mesh2D mesh(samplePoints, polygon, mesh2dInstances[meshKernelId]->m_projection);
            *mesh2dInstances[meshKernelId] += mesh;
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_copy_mesh_boundaries_to_polygon(int meshKernelId, GeometryList& geometryList)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            const std::vector<meshkernel::Point> polygonNodes;
            const auto meshBoundaryPolygon = mesh2dInstances[meshKernelId]->MeshBoundaryToPolygon(polygonNodes);

            ConvertPointVectorToGeometryList(meshBoundaryPolygon, geometryList);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_copy_mesh_boundaries_to_polygon_count_nodes(int meshKernelId, int& numberOfPolygonNodes)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            const std::vector<meshkernel::Point> polygonNodes;
            const auto meshBoundaryPolygon = mesh2dInstances[meshKernelId]->MeshBoundaryToPolygon(polygonNodes);
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
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            auto polygonPoints = ConvertGeometryListToPointVector(geometryListIn);

            const meshkernel::Polygons polygon(polygonPoints, mesh2dInstances[meshKernelId]->m_projection);
            const auto refinedPolygon = polygon.RefineFirstPolygon(firstIndex, secondIndex, distance);

            ConvertPointVectorToGeometryList(refinedPolygon, geometryListOut);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_refine_polygon_count(int meshKernelId, GeometryList& geometryListIn, int firstIndex, int secondIndex, double distance, int& numberOfPolygonNodes)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            auto polygonPoints = ConvertGeometryListToPointVector(geometryListIn);

            const meshkernel::Polygons polygon(polygonPoints, mesh2dInstances[meshKernelId]->m_projection);

            const auto refinedPolygon = polygon.RefineFirstPolygon(firstIndex, secondIndex, distance);

            numberOfPolygonNodes = int(refinedPolygon.size());
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_merge_nodes(int meshKernelId, const GeometryList& geometryListIn)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            auto polygonPoints = ConvertGeometryListToPointVector(geometryListIn);

            const meshkernel::Polygons polygon(polygonPoints, mesh2dInstances[meshKernelId]->m_projection);

            mesh2dInstances[meshKernelId]->MergeNodesInPolygon(polygon);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_merge_two_nodes(int meshKernelId, int startNode, int endNode)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            mesh2dInstances[meshKernelId]->MergeTwoNodes(startNode, endNode);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_nodes_in_polygons(int meshKernelId, GeometryList& geometryListIn, int inside, int numberOfMeshNodes, int** selectedNodes)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            auto polygonPoints = ConvertGeometryListToPointVector(geometryListIn);

            const meshkernel::Polygons polygon(polygonPoints, mesh2dInstances[meshKernelId]->m_projection);

            const bool selectInside = inside == 1 ? true : false;
            mesh2dInstances[meshKernelId]->MaskNodesInPolygons(polygon, selectInside);

            int index = 0;
            for (auto i = 0; i < mesh2dInstances[meshKernelId]->GetNumNodes(); ++i)
            {
                if (mesh2dInstances[meshKernelId]->m_nodeMask[i] > 0)
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

    MKERNEL_API int mkernel_count_nodes_in_polygons(int meshKernelId, GeometryList& geometryListIn, int inside, int& numberOfMeshNodes)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            auto polygonPoints = ConvertGeometryListToPointVector(geometryListIn);

            const meshkernel::Polygons polygon(polygonPoints, mesh2dInstances[meshKernelId]->m_projection);

            const bool selectInside = inside == 1 ? true : false;
            mesh2dInstances[meshKernelId]->MaskNodesInPolygons(polygon, selectInside);

            numberOfMeshNodes = 0;
            for (auto i = 0; i < mesh2dInstances[meshKernelId]->GetNumNodes(); ++i)
            {
                if (mesh2dInstances[meshKernelId]->m_nodeMask[i] > 0)
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

    MKERNEL_API int mkernel_insert_edge(int meshKernelId, int startNode, int endNode, int& new_edge_index)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            new_edge_index = mesh2dInstances[meshKernelId]->ConnectNodes(startNode, endNode);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_insert_node(int meshKernelId, double xCoordinate, double yCoordinate, int& nodeIndex)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                //create a valid instance, by default cartesian
                *mesh2dInstances[meshKernelId] = meshkernel::Mesh2D();
                mesh2dInstances[meshKernelId]->m_projection = meshkernel::Projection::cartesian;
            }

            const meshkernel::Point newNode{xCoordinate, yCoordinate};
            nodeIndex = mesh2dInstances[meshKernelId]->InsertNode(newNode);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_delete_node(int meshKernelId, int nodeIndex)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            mesh2dInstances[meshKernelId]->DeleteNode(nodeIndex);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_move_node(int meshKernelId, const GeometryList& geometryListIn, int nodeIndex)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            auto newPoint = ConvertGeometryListToPointVector(geometryListIn);

            mesh2dInstances[meshKernelId]->MoveNode(newPoint[0], nodeIndex);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_delete_edge(int meshKernelId, const GeometryList& geometryListIn)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            auto newPoint = ConvertGeometryListToPointVector(geometryListIn);

            const auto edgeIndex = mesh2dInstances[meshKernelId]->FindEdgeCloseToAPoint(newPoint[0]);

            mesh2dInstances[meshKernelId]->DeleteEdge(edgeIndex);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_find_edge(int meshKernelId, const GeometryList& geometryListIn, int& edgeIndex)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            auto newPoint = ConvertGeometryListToPointVector(geometryListIn);

            edgeIndex = static_cast<int>(mesh2dInstances[meshKernelId]->FindEdgeCloseToAPoint(newPoint[0]));
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_offsetted_polygon(int meshKernelId, const GeometryList& geometryListIn, bool innerPolygon, double distance, GeometryList& geometryListOut)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            auto polygonPoints = ConvertGeometryListToPointVector(geometryListIn);

            const meshkernel::Polygons polygon(polygonPoints, mesh2dInstances[meshKernelId]->m_projection);

            const auto newPolygon = polygon.OffsetCopy(distance, innerPolygon);

            ConvertPointVectorToGeometryList(newPolygon.m_nodes, geometryListOut);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_offsetted_polygon_count(int meshKernelId, const GeometryList& geometryListIn, bool innerPolygon, double distance, int& numberOfPolygonNodes)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            auto polygonPoints = ConvertGeometryListToPointVector(geometryListIn);

            const meshkernel::Polygons polygon(polygonPoints, mesh2dInstances[meshKernelId]->m_projection);
            const auto newPolygon = polygon.OffsetCopy(distance, innerPolygon);

            numberOfPolygonNodes = newPolygon.GetNumNodes();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_refine_mesh_based_on_samples(int meshKernelId,
                                                         const GeometryList& geometryListIn,
                                                         const InterpolationParameters& interpolationParameters,
                                                         const SampleRefineParameters& sampleRefineParameters)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            if (mesh2dInstances[meshKernelId]->GetNumNodes() <= 0)
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

            const auto averaging = std::make_shared<meshkernel::AveragingInterpolation>(mesh2dInstances[meshKernelId],
                                                                                        samples,
                                                                                        averagingMethod,
                                                                                        meshkernel::MeshLocations::Faces,
                                                                                        1.0,
                                                                                        refineOutsideFace,
                                                                                        transformSamples);

            meshkernel::MeshRefinement meshRefinement(mesh2dInstances[meshKernelId], averaging, sampleRefineParameters, interpolationParameters);
            meshRefinement.Compute();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_refine_mesh_based_on_polygon(int meshKernelId, const GeometryList& geometryList, const InterpolationParameters& interpolationParameters)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            if (mesh2dInstances[meshKernelId]->GetNumNodes() <= 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh has no nodes.");
            }

            auto points = ConvertGeometryListToPointVector(geometryList);

            const meshkernel::Polygons polygon(points, mesh2dInstances[meshKernelId]->m_projection);

            meshkernel::MeshRefinement meshRefinement(mesh2dInstances[meshKernelId], polygon, interpolationParameters);
            meshRefinement.Compute();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_node_index(int meshKernelId, const GeometryList& geometryListIn, double searchRadius, int& nodeIndex)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            if (mesh2dInstances[meshKernelId]->GetNumNodes() <= 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh has no nodes.");
            }

            auto polygonPoints = ConvertGeometryListToPointVector(geometryListIn);

            nodeIndex = static_cast<int>(mesh2dInstances[meshKernelId]->FindNodeCloseToAPoint(polygonPoints[0], searchRadius));
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_node_coordinate(int meshKernelId, const GeometryList& geometryListIn, double searchRadius, GeometryList& geometryListOut)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            if (mesh2dInstances[meshKernelId]->GetNumNodes() <= 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh has no nodes.");
            }

            if (geometryListOut.numberOfCoordinates <= 0)
            {
                throw std::invalid_argument("MeshKernel: The output-geometry has no coordinates.");
            }

            auto polygonPoints = ConvertGeometryListToPointVector(geometryListIn);

            const auto nodeIndex = mesh2dInstances[meshKernelId]->FindNodeCloseToAPoint(polygonPoints[0], searchRadius);

            // Set the node coordinate
            const auto node = mesh2dInstances[meshKernelId]->m_nodes[nodeIndex];
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

    MKERNEL_API int mkernel_curvilinear_mesh_from_splines_ortho(int meshKernelId,
                                                                const GeometryList& geometryListIn,
                                                                const CurvilinearParameters& curvilinearParameters,
                                                                const SplinesToCurvilinearParameters& splinesToCurvilinearParameters)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            // use the default constructor, no instance present
            const auto spline = std::make_shared<meshkernel::Splines>(mesh2dInstances[meshKernelId]->m_projection);
            SetSplines(geometryListIn, *spline);

            meshkernel::CurvilinearGridFromSplines curvilinearGridFromSplines(spline, curvilinearParameters, splinesToCurvilinearParameters);

            meshkernel::CurvilinearGrid curvilinearGrid;
            curvilinearGridFromSplines.Compute(curvilinearGrid);
            *mesh2dInstances[meshKernelId] += meshkernel::Mesh2D(curvilinearGrid, mesh2dInstances[meshKernelId]->m_projection);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_mesh_from_splines_ortho_initialize(int meshKernelId, const GeometryList& geometryList, const CurvilinearParameters& curvilinearParameters, const SplinesToCurvilinearParameters& splinesToCurvilinearParameters)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            auto spline = std::make_shared<meshkernel::Splines>(mesh2dInstances[meshKernelId]->m_projection);
            SetSplines(geometryList, *spline);

            auto curvilinearGridFromSplines = std::make_shared<meshkernel::CurvilinearGridFromSplines>(spline, curvilinearParameters, splinesToCurvilinearParameters);

            curvilinearGridFromSplinesInstances.insert({meshKernelId, curvilinearGridFromSplines});

            curvilinearGridFromSplinesInstances[meshKernelId]->Initialize();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_mesh_from_splines_ortho_iteration(int meshKernelId, int layer)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            curvilinearGridFromSplinesInstances[meshKernelId]->Iterate(layer);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_mesh_from_splines_ortho_refresh_mesh(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            meshkernel::CurvilinearGrid curvilinearGrid;
            curvilinearGridFromSplinesInstances[meshKernelId]->ComputeCurvilinearGrid(curvilinearGrid);

            *mesh2dInstances[meshKernelId] += meshkernel::Mesh2D(curvilinearGrid, mesh2dInstances[meshKernelId]->m_projection);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_mesh_from_splines_ortho_delete(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            curvilinearGridFromSplinesInstances.erase(meshKernelId);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_points_in_polygon(int meshKernelId, const GeometryList& polygon, const GeometryList& pointsNative, GeometryList& selectedPointsNative)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            auto polygonNodes = ConvertGeometryListToPointVector(polygon);

            auto points = ConvertGeometryListToPointVector(pointsNative);
            const meshkernel::Polygons localPolygon(polygonNodes, mesh2dInstances[meshKernelId]->m_projection);

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

    MKERNEL_API int mkernel_flip_edges(int meshKernelId,
                                       int isTriangulationRequired,
                                       int projectToLandBoundaryRequired)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            //set landboundaries
            auto polygon = std::make_shared<meshkernel::Polygons>();

            std::vector<meshkernel::Point> landBoundary;
            const auto landBoundaries = std::make_shared<meshkernel::LandBoundaries>(landBoundary, mesh2dInstances[meshKernelId], polygon);

            const bool triangulateFaces = isTriangulationRequired == 0 ? false : true;
            const bool projectToLandBoundary = projectToLandBoundaryRequired == 0 ? false : true;
            const meshkernel::FlipEdges flipEdges(mesh2dInstances[meshKernelId], landBoundaries, triangulateFaces, projectToLandBoundary);

            flipEdges.Compute();
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_mesh_from_splines(int meshKernelId,
                                                          const GeometryList& splines,
                                                          const CurvilinearParameters& curvilinearParameters)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            // Use the default constructor, no instance present
            const auto meshKernelSplines = std::make_shared<meshkernel::Splines>(mesh2dInstances[meshKernelId]->m_projection);
            SetSplines(splines, *meshKernelSplines);

            // Create algorithm and set the splines
            meshkernel::CurvilinearGridFromSplinesTransfinite curvilinearGridFromSplinesTransfinite(meshKernelSplines, curvilinearParameters);

            // Compute the curvilinear grid
            meshkernel::CurvilinearGrid curvilinearGrid;
            curvilinearGridFromSplinesTransfinite.Compute(curvilinearGrid);

            // Transform and set mesh pointer
            *mesh2dInstances[meshKernelId] += meshkernel::Mesh2D(curvilinearGrid, mesh2dInstances[meshKernelId]->m_projection);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_from_polygon(int meshKernelId,
                                                     const GeometryList& polygons,
                                                     int firstNode,
                                                     int secondNode,
                                                     int thirdNode,
                                                     bool useFourthSide)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            auto polygonPoints = ConvertGeometryListToPointVector(polygons);

            const auto localPolygon = std::make_shared<meshkernel::Polygons>(polygonPoints, mesh2dInstances[meshKernelId]->m_projection);

            meshkernel::CurvilinearGrid curvilinearGrid;
            const meshkernel::CurvilinearGridFromPolygon curvilinearGridFromPolygon(localPolygon);
            curvilinearGridFromPolygon.Compute(firstNode, secondNode, thirdNode, useFourthSide, curvilinearGrid);

            // convert to curvilinear grid and add it to the current mesh
            *mesh2dInstances[meshKernelId] += meshkernel::Mesh2D(curvilinearGrid, mesh2dInstances[meshKernelId]->m_projection);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_from_triangle(int meshKernelId,
                                                      const GeometryList& polygon,
                                                      int firstNode,
                                                      int secondNode,
                                                      int thirdNode)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            auto polygonPoints = ConvertGeometryListToPointVector(polygon);

            const auto localPolygon = std::make_shared<meshkernel::Polygons>(polygonPoints, mesh2dInstances[meshKernelId]->m_projection);

            meshkernel::CurvilinearGrid curvilinearGrid;
            const meshkernel::CurvilinearGridFromPolygon curvilinearGridFromPolygon(localPolygon);
            curvilinearGridFromPolygon.Compute(firstNode, secondNode, thirdNode, curvilinearGrid);

            // convert to curvilinear grid and add it to the current mesh
            *mesh2dInstances[meshKernelId] += meshkernel::Mesh2D(curvilinearGrid, mesh2dInstances[meshKernelId]->m_projection);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_small_flow_edge_centers_count(int meshKernelId, double smallFlowEdgesThreshold, int& numSmallFlowEdges)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            const auto edgesCrossingSmallFlowEdges = mesh2dInstances[meshKernelId]->GetEdgesCrossingSmallFlowEdges(smallFlowEdgesThreshold);
            const auto smallFlowEdgeCenters = mesh2dInstances[meshKernelId]->GetFlowEdgesCenters(edgesCrossingSmallFlowEdges);

            numSmallFlowEdges = static_cast<int>(smallFlowEdgeCenters.size());
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_small_flow_edge_centers(int meshKernelId, double smallFlowEdgesThreshold, GeometryList& result)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            const auto edgesCrossingSmallFlowEdges = mesh2dInstances[meshKernelId]->GetEdgesCrossingSmallFlowEdges(smallFlowEdgesThreshold);
            const auto smallFlowEdgeCenters = mesh2dInstances[meshKernelId]->GetFlowEdgesCenters(edgesCrossingSmallFlowEdges);

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
        invalidIndex = meshGeometryError.m_invalidIndex;
        type = static_cast<int>(meshGeometryError.m_location);
        return Success;
    }

    MKERNEL_API int mkernel_get_obtuse_triangles_count(int meshKernelId, int& numObtuseTriangles)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            const auto obtuseTriangles = mesh2dInstances[meshKernelId]->GetObtuseTrianglesCenters();

            numObtuseTriangles = static_cast<int>(obtuseTriangles.size());
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_obtuse_triangles(int meshKernelId, GeometryList& result)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            const auto obtuseTriangles = mesh2dInstances[meshKernelId]->GetObtuseTrianglesCenters();

            ConvertPointVectorToGeometryList(obtuseTriangles, result);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_delete_small_flow_edges(int meshKernelId, double smallFlowEdgesThreshold, double minFractionalAreaTriangles)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            mesh2dInstances[meshKernelId]->DeleteSmallFlowEdges(smallFlowEdgesThreshold);
            mesh2dInstances[meshKernelId]->DeleteSmallTrianglesAtBoundaries(minFractionalAreaTriangles);
        }
        catch (...)
        {
            exitCode = HandleExceptions(std::current_exception());
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_compute_single_contacts(int meshKernelId,
                                                    int** oneDNodeMask,
                                                    const GeometryList& polygons,
                                                    Contacts& contacts)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= mesh2dInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            // Convert 1D node mask from int** to vector<bool>
            auto num1DNodes = mesh1dInstances[meshKernelId]->GetNumNodes();
            std::vector<bool> meshKernel1DNodeMask(num1DNodes);
            for (auto i = 0; i < num1DNodes; ++i)
            {
                switch ((*oneDNodeMask)[i])
                {
                case 0:
                    meshKernel1DNodeMask[i] = false;
                    break;
                case 1:
                    meshKernel1DNodeMask[i] = true;
                    break;
                default:
                    throw std::invalid_argument("MeshKernel: Invalid 1D mask.");
                    break;
                }
            }

            // Init Contacts instance
            meshkernel::Contacts meshKernelContacts(mesh1dInstances[meshKernelId],
                                                    mesh2dInstances[meshKernelId]);

            // Convert polygon date from GeometryList to Polygons
            auto polygonPoints = ConvertGeometryListToPointVector(polygons);
            const meshkernel::Polygons meshKernelPolygons(polygonPoints,
                                                          mesh2dInstances[meshKernelId]->m_projection);

            // Execute
            meshKernelContacts.ComputeSingleContacts(meshKernelPolygons,
                                                     meshKernel1DNodeMask);
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

    MKERNEL_API int averaging(const MeshGeometryDimensions& meshGeometryDimensions,
                              const MeshGeometry& meshGeometry,
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
            const auto edges = meshkernel::ConvertToEdgeNodesVector(meshGeometryDimensions.numedge, meshGeometry.edge_nodes);
            const auto nodes = meshkernel::ConvertToNodesVector(meshGeometryDimensions.numnode, meshGeometry.nodex, meshGeometry.nodey);
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
    MKERNEL_API int triangulation(const MeshGeometryDimensions& meshGeometryDimensions,
                                  const MeshGeometry& meshGeometry,
                                  int& startIndex,
                                  const double** samplesXCoordinate,
                                  const double** samplesYCoordinate,
                                  const double** samplesValue,
                                  int& numSamples,
                                  double** results,
                                  int& locationType,
                                  int& spherical,
                                  int& sphericalAccurate)
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
            const auto locations = ComputeLocations(meshGeometryDimensions, meshGeometry, location);

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
