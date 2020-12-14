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

#include <map>
#include <stdexcept>
#include <vector>

#include <MeshKernel/AveragingInterpolation.hpp>
#include <MeshKernel/Constants.hpp>
#include <MeshKernel/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGridFromPolygon.hpp>
#include <MeshKernel/CurvilinearGridFromSplines.hpp>
#include <MeshKernel/CurvilinearGridFromSplinesTransfinite.hpp>
#include <MeshKernel/CurvilinearParameters.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/FlipEdges.hpp>
#include <MeshKernel/LandBoundaries.hpp>
#include <MeshKernel/Mesh.hpp>
#include <MeshKernel/MeshKernel.hpp>
#include <MeshKernel/MeshRefinement.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/OrthogonalizationAndSmoothing.hpp>
#include <MeshKernel/Orthogonalizer.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/Smoother.hpp>
#include <MeshKernel/Splines.hpp>
#include <MeshKernel/SplinesToCurvilinearParameters.hpp>
#include <MeshKernel/TriangulationInterpolation.hpp>

namespace meshkernelapi
{
    // The vector containing the mesh instances
    static std::vector<std::shared_ptr<meshkernel::Mesh>> meshInstances;

    // For interactivity
    static std::map<int, std::shared_ptr<meshkernel::OrthogonalizationAndSmoothing>> orthogonalizationInstances;
    static std::map<int, std::shared_ptr<meshkernel::CurvilinearGridFromSplines>> curvilinearGridFromSplinesInstances;

    enum error
    {
        Success = 0,         // 0b0000
        Exception = 1,       // 0b0001
        InvalidGeometry = 2, // 0b0010
    };

    static char exceptionMessage[512] = "";

    // TODO: Return result instead of relying on second input parameter
    static void ConvertGeometryListToPointVector(const GeometryList& geometryListIn, std::vector<meshkernel::Point>& result)
    {
        if (geometryListIn.numberOfCoordinates == 0)
        {
            return;
        }
        result.resize(geometryListIn.numberOfCoordinates);

        for (int i = 0; i < geometryListIn.numberOfCoordinates; i++)
        {
            result[i] = {geometryListIn.xCoordinates[i], geometryListIn.yCoordinates[i]};
        }
    }

    // TODO: Return result instead of relying on second input parameter
    static void ConvertGeometryListToSampleVector(const GeometryList& geometryListIn, std::vector<meshkernel::Sample>& result)
    {
        if (geometryListIn.numberOfCoordinates == 0)
        {
            throw std::invalid_argument("MeshKernel: The samples are empty.");
        }
        result.resize(geometryListIn.numberOfCoordinates);

        for (int i = 0; i < geometryListIn.numberOfCoordinates; i++)
        {
            result[i] = {geometryListIn.xCoordinates[i], geometryListIn.yCoordinates[i], geometryListIn.zCoordinates[i]};
        }
    }

    // TODO: Return result instead of relying on second input parameter
    static void ConvertPointVectorToGeometryList(std::vector<meshkernel::Point> pointVector, GeometryList& result)
    {
        if (pointVector.size() < result.numberOfCoordinates)
        {
            throw std::invalid_argument("MeshKernel: Invalid memory allocation, the point-vector size is smaller than the number of coordinates.");
        }

        for (int i = 0; i < result.numberOfCoordinates; i++)
        {
            result.xCoordinates[i] = pointVector[i].x;
            result.yCoordinates[i] = pointVector[i].y;
        }
    }

    static bool SetSplines(const GeometryList& geometryListIn, meshkernel::Splines& spline)
    {
        if (geometryListIn.numberOfCoordinates == 0)
        {
            return false;
        }

        std::vector<meshkernel::Point> splineCornerPoints;
        ConvertGeometryListToPointVector(geometryListIn, splineCornerPoints);

        const auto indexes = FindIndices(splineCornerPoints, 0, splineCornerPoints.size(), meshkernel::doubleMissingValue);

        for (const auto& index : indexes)
        {
            const auto size = index[1] - index[0] + 1;
            if (size > 0)
            {
                spline.AddSpline(splineCornerPoints, index[0], size);
            }
        }

        return true;
    }

    static bool SetMeshGeometry(int meshKernelId, MeshGeometryDimensions& meshGeometryDimensions, MeshGeometry& meshGeometry)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return false;
        }

        meshGeometry.nodex = &(meshInstances[meshKernelId]->m_nodex[0]);
        meshGeometry.nodey = &(meshInstances[meshKernelId]->m_nodey[0]);
        meshGeometry.nodez = &(meshInstances[meshKernelId]->m_nodez[0]);
        meshGeometry.edge_nodes = &(meshInstances[meshKernelId]->m_edgeNodes[0]);

        meshGeometryDimensions.maxnumfacenodes = meshkernel::maximumNumberOfNodesPerFace;
        meshGeometryDimensions.numface = meshInstances[meshKernelId]->GetNumFaces();
        if (meshGeometryDimensions.numface > 0)
        {
            meshGeometry.face_nodes = &(meshInstances[meshKernelId]->m_faceNodes[0]);
            meshGeometry.facex = &(meshInstances[meshKernelId]->m_facesCircumcentersx[0]);
            meshGeometry.facey = &(meshInstances[meshKernelId]->m_facesCircumcentersy[0]);
            meshGeometry.facez = &(meshInstances[meshKernelId]->m_facesCircumcentersz[0]);
        }

        if (meshInstances[meshKernelId]->GetNumNodes() == 1)
        {
            meshGeometryDimensions.numnode = 0;
            meshGeometryDimensions.numedge = 0;
        }
        else
        {
            meshGeometryDimensions.numnode = meshInstances[meshKernelId]->GetNumNodes();
            meshGeometryDimensions.numedge = meshInstances[meshKernelId]->GetNumEdges();
        }

        return true;
    }

    static std::vector<meshkernel::Point> ComputeLocations(const MeshGeometryDimensions& meshGeometryDimensions, const MeshGeometry& meshGeometry, meshkernel::InterpolationLocation interpolationLocation)
    {
        std::vector<meshkernel::Point> locations;
        if (interpolationLocation == meshkernel::InterpolationLocation::Nodes)
        {
            locations = meshkernel::ConvertToNodesVector(meshGeometryDimensions.numnode, meshGeometry.nodex, meshGeometry.nodey);
        }
        if (interpolationLocation == meshkernel::InterpolationLocation::Edges)
        {
            const auto edges = meshkernel::ConvertToEdgeNodesVector(meshGeometryDimensions.numedge, meshGeometry.edge_nodes);
            const auto nodes = meshkernel::ConvertToNodesVector(meshGeometryDimensions.numnode, meshGeometry.nodex, meshGeometry.nodey);
            locations = ComputeEdgeCenters(nodes, edges);
        }
        if (interpolationLocation == meshkernel::InterpolationLocation::Faces)
        {
            locations = meshkernel::ConvertToFaceCentersVector(meshGeometryDimensions.numface, meshGeometry.facex, meshGeometry.facey);
        }

        return locations;
    }

    MKERNEL_API int mkernel_new_mesh(int& meshKernelId)
    {
        meshKernelId = int(meshInstances.size());
        meshInstances.emplace_back(std::make_shared<meshkernel::Mesh>());
        return 0;
    };

    MKERNEL_API int mkernel_deallocate_state(int meshKernelId)
    {
        if (meshKernelId < 0 && meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        meshInstances.erase(meshInstances.begin() + meshKernelId);
        return 0;
    }

    MKERNEL_API int mkernel_delete_mesh(int meshKernelId, const GeometryList& geometryListIn, int deletionOption, bool invertDeletion)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            if (meshInstances[meshKernelId]->GetNumNodes() <= 0)
            {
                return exitCode;
            }
            std::vector<meshkernel::Point> polygonPoints;
            ConvertGeometryListToPointVector(geometryListIn, polygonPoints);

            meshkernel::Polygons polygon(polygonPoints, meshInstances[meshKernelId]->m_projection);
            meshInstances[meshKernelId]->DeleteMesh(polygon, deletionOption, invertDeletion);
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_set_state(int meshKernelId, const MeshGeometryDimensions& meshGeometryDimensions, const MeshGeometry& meshGeometry, bool isGeographic)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            const auto edges = meshkernel::ConvertToEdgeNodesVector(meshGeometryDimensions.numedge, meshGeometry.edge_nodes);
            const auto nodes = meshkernel::ConvertToNodesVector(meshGeometryDimensions.numnode, meshGeometry.nodex, meshGeometry.nodey);

            // spherical or cartesian
            if (isGeographic)
            {
                meshInstances[meshKernelId] = std::make_shared<meshkernel::Mesh>(edges, nodes, meshkernel::Projection::spherical);
            }
            else
            {
                meshInstances[meshKernelId] = std::make_shared<meshkernel::Mesh>(edges, nodes, meshkernel::Projection::cartesian);
            }
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_mesh(int meshKernelId, MeshGeometryDimensions& meshGeometryDimensions, MeshGeometry& meshGeometry)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            meshInstances[meshKernelId]->SetFlatCopies(meshkernel::Mesh::AdministrationOptions::AdministrateMeshEdges);

            SetMeshGeometry(meshKernelId, meshGeometryDimensions, meshGeometry);
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_find_faces(int meshKernelId, MeshGeometryDimensions& meshGeometryDimensions, MeshGeometry& meshGeometry)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            meshInstances[meshKernelId]->SetFlatCopies(meshkernel::Mesh::AdministrationOptions::AdministrateMeshEdgesAndFaces);

            SetMeshGeometry(meshKernelId, meshGeometryDimensions, meshGeometry);
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_count_hanging_edges(int meshKernelId, int& numHangingEdges)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            const auto hangingEdges = meshInstances[meshKernelId]->GetHangingEdges();
            numHangingEdges = hangingEdges.size();
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_hanging_edges(int meshKernelId, int** hangingEdgesIndices)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            const auto hangingEdges = meshInstances[meshKernelId]->GetHangingEdges();
            for (auto i = 0; i < hangingEdges.size(); ++i)
            {
                *(hangingEdgesIndices)[i] = hangingEdges[i];
            }
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_delete_hanging_edges(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            meshInstances[meshKernelId]->DeleteHangingEdges();
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_orthogonalize(int meshKernelId,
                                          int projectToLandBoundaryOption,
                                          const OrthogonalizationParameters& orthogonalizationParameters,
                                          const GeometryList& geometryListPolygon,
                                          const GeometryList& geometryListLandBoundaries)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            if (meshInstances[meshKernelId]->GetNumNodes() <= 0)
            {
                return exitCode;
            }

            // build enclosing polygon
            std::vector<meshkernel::Point> nodes(geometryListPolygon.numberOfCoordinates);
            for (int i = 0; i < geometryListPolygon.numberOfCoordinates; i++)
            {
                nodes[i].x = geometryListPolygon.xCoordinates[i];
                nodes[i].y = geometryListPolygon.yCoordinates[i];
            }

            auto polygon = std::make_shared<meshkernel::Polygons>(nodes, meshInstances[meshKernelId]->m_projection);

            // build land boundary
            std::vector<meshkernel::Point> landBoundaries(geometryListLandBoundaries.numberOfCoordinates);
            for (int i = 0; i < geometryListLandBoundaries.numberOfCoordinates; i++)
            {
                landBoundaries[i].x = geometryListLandBoundaries.xCoordinates[i];
                landBoundaries[i].y = geometryListLandBoundaries.yCoordinates[i];
            }

            auto orthogonalizer = std::make_shared<meshkernel::Orthogonalizer>(meshInstances[meshKernelId]);
            auto smoother = std::make_shared<meshkernel::Smoother>(meshInstances[meshKernelId]);
            auto landBoundary = std::make_shared<meshkernel::LandBoundaries>(landBoundaries, meshInstances[meshKernelId], polygon);

            meshkernel::OrthogonalizationAndSmoothing ortogonalization(meshInstances[meshKernelId],
                                                                       smoother,
                                                                       orthogonalizer,
                                                                       polygon,
                                                                       landBoundary,
                                                                       static_cast<meshkernel::LandBoundaries::ProjectToLandBoundaryOption>(projectToLandBoundaryOption),
                                                                       orthogonalizationParameters);
            ortogonalization.Initialize();
            ortogonalization.Compute();
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
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
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            if (meshInstances[meshKernelId]->GetNumNodes() <= 0)
            {
                return exitCode;
            }

            // build enclosing polygon
            std::vector<meshkernel::Point> nodes(geometryListPolygon.numberOfCoordinates);
            for (int i = 0; i < geometryListPolygon.numberOfCoordinates; i++)
            {
                nodes[i].x = geometryListPolygon.xCoordinates[i];
                nodes[i].y = geometryListPolygon.yCoordinates[i];
            }

            // build land boundary
            std::vector<meshkernel::Point> landBoundaries(geometryListLandBoundaries.numberOfCoordinates);
            for (int i = 0; i < geometryListLandBoundaries.numberOfCoordinates; i++)
            {
                landBoundaries[i].x = geometryListLandBoundaries.xCoordinates[i];
                landBoundaries[i].y = geometryListLandBoundaries.yCoordinates[i];
            }

            auto orthogonalizer = std::make_shared<meshkernel::Orthogonalizer>(meshInstances[meshKernelId]);
            auto smoother = std::make_shared<meshkernel::Smoother>(meshInstances[meshKernelId]);
            auto polygon = std::make_shared<meshkernel::Polygons>(nodes, meshInstances[meshKernelId]->m_projection);
            auto landBoundary = std::make_shared<meshkernel::LandBoundaries>(landBoundaries, meshInstances[meshKernelId], polygon);

            auto orthogonalizationInstance = std::make_shared<meshkernel::OrthogonalizationAndSmoothing>(meshInstances[meshKernelId],
                                                                                                         smoother,
                                                                                                         orthogonalizer,
                                                                                                         polygon,
                                                                                                         landBoundary,
                                                                                                         static_cast<meshkernel::LandBoundaries::ProjectToLandBoundaryOption>(projectToLandBoundaryOption),
                                                                                                         orthogonalizationParameters);
            orthogonalizationInstance->Initialize();

            orthogonalizationInstances.insert({meshKernelId, orthogonalizationInstance});
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_orthogonalize_prepare_outer_iteration(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            if (meshInstances[meshKernelId]->GetNumNodes() <= 0)
            {
                return exitCode;
            }

            orthogonalizationInstances[meshKernelId]->PrepareOuterIteration();
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_orthogonalize_inner_iteration(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            if (meshInstances[meshKernelId]->GetNumNodes() <= 0)
            {
                return exitCode;
            }

            orthogonalizationInstances[meshKernelId]->InnerIteration();
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_orthogonalize_finalize_outer_iteration(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            if (meshInstances[meshKernelId]->GetNumNodes() <= 0)
            {
                return exitCode;
            }

            orthogonalizationInstances[meshKernelId]->FinalizeOuterIteration();
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_orthogonalize_delete(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            if (meshInstances[meshKernelId]->GetNumNodes() <= 0)
            {
                return exitCode;
            }

            orthogonalizationInstances.erase(meshKernelId);
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_orthogonality(int meshKernelId, GeometryList& geometryList)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            if (meshInstances[meshKernelId]->GetNumNodes() <= 0)
            {
                return exitCode;
            }

            const auto result = meshInstances[meshKernelId]->GetOrthogonality();

            for (auto i = 0; i < geometryList.numberOfCoordinates; ++i)
            {
                geometryList.zCoordinates[i] = result[i];
            }
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_smoothness(int meshKernelId, GeometryList& geometryList)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            if (meshInstances[meshKernelId]->GetNumNodes() <= 0)
            {
                return exitCode;
            }

            const auto result = meshInstances[meshKernelId]->GetSmoothness();

            for (auto i = 0; i < geometryList.numberOfCoordinates; ++i)
            {
                geometryList.zCoordinates[i] = result[i];
            }
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
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
            for (int i = 0; i < geometryListIn.numberOfCoordinates; i++)
            {
                splines[i].x = geometryListIn.xCoordinates[i];
                splines[i].y = geometryListIn.yCoordinates[i];
            }

            const auto indexes = FindIndices(splines, 0, splines.size(), meshkernel::doubleMissingValue);
            const auto numSplines = indexes.size();
            std::vector<meshkernel::Point> coordinatesDerivatives(geometryListIn.numberOfCoordinates);

            int index = 0;
            for (int s = 0; s < numSplines; s++)
            {
                std::vector<meshkernel::Point> coordinates(splines.begin() + indexes[s][0], splines.begin() + int(indexes[s][1]) + 1);
                int numNodes = int(indexes[s][1]) - int(indexes[s][0]) + 1;
                meshkernel::Splines::SecondOrderDerivative(coordinates, numNodes, coordinatesDerivatives);

                for (int n = 0; n < numNodes - 1; n++)
                {
                    for (int p = 0; p <= numberOfPointsBetweenNodes; p++)
                    {

                        double pointAdimensionalCoordinate = n + double(p) / double(numberOfPointsBetweenNodes);
                        meshkernel::Point pointCoordinate{meshkernel::doubleMissingValue, meshkernel::doubleMissingValue};
                        bool successful = InterpolateSplinePoint(coordinates, coordinatesDerivatives, pointAdimensionalCoordinate, pointCoordinate);
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
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_make_mesh(int meshKernelId, const MakeMeshParameters& makeGridParameters, const GeometryList& geometryList)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            std::vector<meshkernel::Point> result;
            ConvertGeometryListToPointVector(geometryList, result);

            meshkernel::Polygons polygon(result, meshInstances[meshKernelId]->m_projection);

            meshkernel::Mesh mesh;
            mesh.MakeMesh(makeGridParameters, polygon);

            *meshInstances[meshKernelId] += mesh;
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_make_mesh_from_polygon(int meshKernelId, const GeometryList& disposableGeometryListIn)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            std::vector<meshkernel::Point> result;
            ConvertGeometryListToPointVector(disposableGeometryListIn, result);

            meshkernel::Polygons polygon(result, meshInstances[meshKernelId]->m_projection);

            // generate samples in all polygons
            const auto generatedPoints = polygon.ComputePointsInPolygons();

            meshkernel::Mesh mesh(generatedPoints[0], polygon, meshInstances[meshKernelId]->m_projection);
            *meshInstances[meshKernelId] += mesh;
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_make_mesh_from_samples(int meshKernelId, GeometryList& geometryList)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            std::vector<meshkernel::Point> samplePoints;
            ConvertGeometryListToPointVector(geometryList, samplePoints);

            meshkernel::Polygons polygon;
            meshkernel::Mesh mesh(samplePoints, polygon, meshInstances[meshKernelId]->m_projection);
            *meshInstances[meshKernelId] += mesh;
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_copy_mesh_boundaries_to_polygon(int meshKernelId, GeometryList& geometryList)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            std::vector<meshkernel::Point> polygonNodes;
            const auto meshBoundaryPolygon = meshInstances[meshKernelId]->MeshBoundaryToPolygon(polygonNodes);

            ConvertPointVectorToGeometryList(meshBoundaryPolygon, geometryList);
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_copy_mesh_boundaries_to_polygon_count_nodes(int meshKernelId, int& numberOfPolygonNodes)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            std::vector<meshkernel::Point> polygonNodes;
            const auto meshBoundaryPolygon = meshInstances[meshKernelId]->MeshBoundaryToPolygon(polygonNodes);
            numberOfPolygonNodes = static_cast<int>(meshBoundaryPolygon.size() - 1); // last value is a separator
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_refine_polygon(int meshKernelId, const GeometryList& geometryListIn, int firstIndex, int secondIndex, double distance, GeometryList& geometryListOut)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            std::vector<meshkernel::Point> polygonPoints;
            ConvertGeometryListToPointVector(geometryListIn, polygonPoints);

            meshkernel::Polygons polygon(polygonPoints, meshInstances[meshKernelId]->m_projection);
            const auto refinedPolygon = polygon.RefineFirstPolygon(firstIndex, secondIndex, distance);

            ConvertPointVectorToGeometryList(refinedPolygon, geometryListOut);
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_refine_polygon_count(int meshKernelId, GeometryList& geometryListIn, int firstIndex, int secondIndex, double distance, int& numberOfPolygonNodes)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            std::vector<meshkernel::Point> polygonPoints;
            ConvertGeometryListToPointVector(geometryListIn, polygonPoints);

            meshkernel::Polygons polygon(polygonPoints, meshInstances[meshKernelId]->m_projection);

            const auto refinedPolygon = polygon.RefineFirstPolygon(firstIndex, secondIndex, distance);

            numberOfPolygonNodes = int(refinedPolygon.size());
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_merge_nodes(int meshKernelId, const GeometryList& geometryListIn)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            std::vector<meshkernel::Point> polygonPoints;
            ConvertGeometryListToPointVector(geometryListIn, polygonPoints);

            meshkernel::Polygons polygon(polygonPoints, meshInstances[meshKernelId]->m_projection);

            meshInstances[meshKernelId]->MergeNodesInPolygon(polygon);
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_merge_two_nodes(int meshKernelId, int startNode, int endNode)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            meshInstances[meshKernelId]->MergeTwoNodes(startNode, endNode);
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_nodes_in_polygons(int meshKernelId, GeometryList& geometryListIn, int inside, int numberOfMeshNodes, int** selectedNodes)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            std::vector<meshkernel::Point> polygonPoints;
            ConvertGeometryListToPointVector(geometryListIn, polygonPoints);

            meshkernel::Polygons polygon(polygonPoints, meshInstances[meshKernelId]->m_projection);

            bool selectInside = inside == 1 ? true : false;
            meshInstances[meshKernelId]->MaskNodesInPolygons(polygon, selectInside);

            int index = 0;
            for (int i = 0; i < meshInstances[meshKernelId]->GetNumNodes(); ++i)
            {
                if (meshInstances[meshKernelId]->m_nodeMask[i] > 0)
                {
                    (*selectedNodes)[index] = i;
                    index++;
                }
            }
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_count_nodes_in_polygons(int meshKernelId, GeometryList& geometryListIn, int inside, int& numberOfMeshNodes)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            std::vector<meshkernel::Point> polygonPoints;
            ConvertGeometryListToPointVector(geometryListIn, polygonPoints);

            meshkernel::Polygons polygon(polygonPoints, meshInstances[meshKernelId]->m_projection);

            bool selectInside = inside == 1 ? true : false;
            meshInstances[meshKernelId]->MaskNodesInPolygons(polygon, selectInside);

            numberOfMeshNodes = 0;
            for (auto i = 0; i < meshInstances[meshKernelId]->GetNumNodes(); ++i)
            {
                if (meshInstances[meshKernelId]->m_nodeMask[i] > 0)
                {
                    numberOfMeshNodes++;
                }
            }
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_insert_edge(int meshKernelId, int startNode, int endNode, int& new_edge_index)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            new_edge_index = meshInstances[meshKernelId]->ConnectNodes(startNode, endNode);
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_insert_node(int meshKernelId, double xCoordinate, double yCoordinate, int& nodeIndex)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                //create a valid instance, by default cartesian
                *meshInstances[meshKernelId] = meshkernel::Mesh();
                meshInstances[meshKernelId]->m_projection = meshkernel::Projection::cartesian;
            }

            meshkernel::Point newNode{xCoordinate, yCoordinate};
            nodeIndex = meshInstances[meshKernelId]->InsertNode(newNode);
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_delete_node(int meshKernelId, int nodeIndex)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            meshInstances[meshKernelId]->DeleteNode(nodeIndex);
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_move_node(int meshKernelId, const GeometryList& geometryListIn, int nodeIndex)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            std::vector<meshkernel::Point> newPoint;
            ConvertGeometryListToPointVector(geometryListIn, newPoint);

            meshInstances[meshKernelId]->MoveNode(newPoint[0], nodeIndex);
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_delete_edge(int meshKernelId, const GeometryList& geometryListIn)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            std::vector<meshkernel::Point> newPoint;
            ConvertGeometryListToPointVector(geometryListIn, newPoint);

            int edgeIndex = meshInstances[meshKernelId]->FindEdgeCloseToAPoint(newPoint[0]);

            meshInstances[meshKernelId]->DeleteEdge(edgeIndex);
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_find_edge(int meshKernelId, const GeometryList& geometryListIn, int& edgeIndex)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            std::vector<meshkernel::Point> newPoint;
            ConvertGeometryListToPointVector(geometryListIn, newPoint);

            edgeIndex = meshInstances[meshKernelId]->FindEdgeCloseToAPoint(newPoint[0]);
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_offsetted_polygon(int meshKernelId, const GeometryList& geometryListIn, bool innerPolygon, double distance, GeometryList& geometryListOut)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            std::vector<meshkernel::Point> polygonPoints;
            ConvertGeometryListToPointVector(geometryListIn, polygonPoints);

            meshkernel::Polygons polygon(polygonPoints, meshInstances[meshKernelId]->m_projection);

            const auto newPolygon = polygon.OffsetCopy(distance, innerPolygon);

            ConvertPointVectorToGeometryList(newPolygon.m_nodes, geometryListOut);
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_offsetted_polygon_count(int meshKernelId, const GeometryList& geometryListIn, bool innerPolygon, double distance, int& numberOfPolygonNodes)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            std::vector<meshkernel::Point> polygonPoints;
            ConvertGeometryListToPointVector(geometryListIn, polygonPoints);

            meshkernel::Polygons polygon(polygonPoints, meshInstances[meshKernelId]->m_projection);
            const auto newPolygon = polygon.OffsetCopy(distance, innerPolygon);

            numberOfPolygonNodes = newPolygon.GetNumNodes();
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
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
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            if (meshInstances[meshKernelId]->GetNumNodes() <= 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh has no nodes.");
            }
            std::vector<meshkernel::Sample> samples;
            ConvertGeometryListToSampleVector(geometryListIn, samples);

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

            const auto averaging = std::make_shared<meshkernel::AveragingInterpolation>(meshInstances[meshKernelId],
                                                                                        samples,
                                                                                        averagingMethod,
                                                                                        meshkernel::InterpolationLocation::Faces,
                                                                                        1.0,
                                                                                        refineOutsideFace,
                                                                                        transformSamples);

            meshkernel::MeshRefinement meshRefinement(meshInstances[meshKernelId], averaging, sampleRefineParameters, interpolationParameters);
            meshRefinement.Compute();
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_refine_mesh_based_on_polygon(int meshKernelId, const GeometryList& geometryList, const InterpolationParameters& interpolationParameters)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            if (meshInstances[meshKernelId]->GetNumNodes() <= 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh has no nodes.");
            }

            std::vector<meshkernel::Point> points;
            ConvertGeometryListToPointVector(geometryList, points);

            meshkernel::Polygons polygon(points, meshInstances[meshKernelId]->m_projection);

            meshkernel::MeshRefinement meshRefinement(meshInstances[meshKernelId], polygon, interpolationParameters);
            meshRefinement.Compute();
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_node_index(int meshKernelId, const GeometryList& geometryListIn, double searchRadius, int& nodeIndex)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            if (meshInstances[meshKernelId]->GetNumNodes() <= 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh has no nodes.");
            }

            std::vector<meshkernel::Point> polygonPoints;
            ConvertGeometryListToPointVector(geometryListIn, polygonPoints);

            nodeIndex = meshInstances[meshKernelId]->GetNodeIndex(polygonPoints[0], searchRadius);
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_node_coordinate(int meshKernelId, const GeometryList& geometryListIn, double searchRadius, GeometryList& geometryListOut)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            if (meshInstances[meshKernelId]->GetNumNodes() <= 0)
            {
                throw std::invalid_argument("MeshKernel: The selected mesh has no nodes.");
            }

            if (geometryListOut.numberOfCoordinates <= 0)
            {
                throw std::invalid_argument("MeshKernel: The output-geometry has no coordinates.");
            }

            std::vector<meshkernel::Point> polygonPoints;
            ConvertGeometryListToPointVector(geometryListIn, polygonPoints);

            int nodeIndex = meshInstances[meshKernelId]->GetNodeIndex(polygonPoints[0], searchRadius);

            // Set the node coordinate
            auto node = meshInstances[meshKernelId]->m_nodes[nodeIndex];
            std::vector<meshkernel::Point> pointVector;
            pointVector.emplace_back(node);
            ConvertPointVectorToGeometryList(pointVector, geometryListOut);
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
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
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            // use the default constructor, no instance present
            auto spline = std::make_shared<meshkernel::Splines>(meshInstances[meshKernelId]->m_projection);
            SetSplines(geometryListIn, *spline);

            meshkernel::CurvilinearGridFromSplines curvilinearGridFromSplines(spline, curvilinearParameters, splinesToCurvilinearParameters);

            meshkernel::CurvilinearGrid curvilinearGrid;
            curvilinearGridFromSplines.Compute(curvilinearGrid);
            *meshInstances[meshKernelId] += meshkernel::Mesh(curvilinearGrid, meshInstances[meshKernelId]->m_projection);
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_mesh_from_splines_ortho_initialize(int meshKernelId, const GeometryList& geometryList, const CurvilinearParameters& curvilinearParameters, const SplinesToCurvilinearParameters& splinesToCurvilinearParameters)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            auto spline = std::make_shared<meshkernel::Splines>(meshInstances[meshKernelId]->m_projection);
            SetSplines(geometryList, *spline);

            auto curvilinearGridFromSplines = std::make_shared<meshkernel::CurvilinearGridFromSplines>(spline, curvilinearParameters, splinesToCurvilinearParameters);

            curvilinearGridFromSplinesInstances.insert({meshKernelId, curvilinearGridFromSplines});

            curvilinearGridFromSplinesInstances[meshKernelId]->Initialize();
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_mesh_from_splines_iteration(int meshKernelId, int layer)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            curvilinearGridFromSplinesInstances[meshKernelId]->Iterate(layer);
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_mesh_from_splines_ortho_refresh_mesh(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            meshkernel::CurvilinearGrid curvilinearGrid;
            curvilinearGridFromSplinesInstances[meshKernelId]->ComputeCurvilinearGrid(curvilinearGrid);

            *meshInstances[meshKernelId] += meshkernel::Mesh(curvilinearGrid, meshInstances[meshKernelId]->m_projection);
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_mesh_from_splines_ortho_delete(int meshKernelId)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            curvilinearGridFromSplinesInstances.erase(meshKernelId);
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_points_in_polygon(int meshKernelId, const GeometryList& polygon, const GeometryList& pointsNative, GeometryList& selectedPointsNative)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            std::vector<meshkernel::Point> polygonNodes;
            ConvertGeometryListToPointVector(polygon, polygonNodes);

            std::vector<meshkernel::Point> points;
            ConvertGeometryListToPointVector(pointsNative, points);
            meshkernel::Polygons localPolygon(polygonNodes, meshInstances[meshKernelId]->m_projection);

            for (int i = 0; i < points.size(); i++)
            {
                selectedPointsNative.zCoordinates[i] = localPolygon.IsPointInPolygon(points[i], 0) ? 1.0 : 0.0;
            }
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
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
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            //set landboundaries
            auto polygon = std::make_shared<meshkernel::Polygons>();

            std::vector<meshkernel::Point> landBoundary;
            auto landBoundaries = std::make_shared<meshkernel::LandBoundaries>(landBoundary, meshInstances[meshKernelId], polygon);

            bool triangulateFaces = isTriangulationRequired == 0 ? false : true;
            bool projectToLandBoundary = projectToLandBoundaryRequired == 0 ? false : true;
            meshkernel::FlipEdges flipEdges(meshInstances[meshKernelId], landBoundaries, triangulateFaces, projectToLandBoundary);

            flipEdges.Compute();
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_mesh_from_splines(int meshKernelId,
                                                          const GeometryList& geometryListIn,
                                                          const CurvilinearParameters& curvilinearParameters)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            // Use the default constructor, no instance present
            auto spline = std::make_shared<meshkernel::Splines>(meshInstances[meshKernelId]->m_projection);
            SetSplines(geometryListIn, *spline);

            // Create algorithm and set the splines
            meshkernel::CurvilinearGridFromSplinesTransfinite curvilinearGridFromSplinesTransfinite(spline, curvilinearParameters);

            // Compute the curvilinear grid
            meshkernel::CurvilinearGrid curvilinearGrid;
            curvilinearGridFromSplinesTransfinite.Compute(curvilinearGrid);

            // Transform and set mesh pointer
            *meshInstances[meshKernelId] += meshkernel::Mesh(curvilinearGrid, meshInstances[meshKernelId]->m_projection);
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_curvilinear_from_polygon(int meshKernelId,
                                                     const GeometryList& polygon,
                                                     int firstNode,
                                                     int secondNode,
                                                     int thirdNode,
                                                     bool useFourthSide)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            std::vector<meshkernel::Point> polygonPoints;
            ConvertGeometryListToPointVector(polygon, polygonPoints);

            auto localPolygon = std::make_shared<meshkernel::Polygons>(polygonPoints, meshInstances[meshKernelId]->m_projection);

            meshkernel::CurvilinearGrid curvilinearGrid;
            meshkernel::CurvilinearGridFromPolygon curvilinearGridFromPolygon(localPolygon);
            curvilinearGridFromPolygon.Compute(firstNode, secondNode, thirdNode, useFourthSide, curvilinearGrid);

            // convert to curvilinear grid and add it to the current mesh
            *meshInstances[meshKernelId] += meshkernel::Mesh(curvilinearGrid, meshInstances[meshKernelId]->m_projection);
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
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
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            std::vector<meshkernel::Point> polygonPoints;
            ConvertGeometryListToPointVector(polygon, polygonPoints);

            auto localPolygon = std::make_shared<meshkernel::Polygons>(polygonPoints, meshInstances[meshKernelId]->m_projection);

            meshkernel::CurvilinearGrid curvilinearGrid;
            meshkernel::CurvilinearGridFromPolygon curvilinearGridFromPolygon(localPolygon);
            curvilinearGridFromPolygon.Compute(firstNode, secondNode, thirdNode, curvilinearGrid);

            // convert to curvilinear grid and add it to the current mesh
            *meshInstances[meshKernelId] += meshkernel::Mesh(curvilinearGrid, meshInstances[meshKernelId]->m_projection);
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_small_flow_edge_centers_count(int meshKernelId, double smallFlowEdgesThreshold, int& numSmallFlowEdges)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }
            const auto edgesCrossingSmallFlowEdges = meshInstances[meshKernelId]->GetEdgesCrossingSmallFlowEdges(smallFlowEdgesThreshold);
            const auto smallFlowEdgeCenters = meshInstances[meshKernelId]->GetFlowEdgesCenters(edgesCrossingSmallFlowEdges);

            numSmallFlowEdges = static_cast<int>(smallFlowEdgeCenters.size());
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_small_flow_edge_centers(int meshKernelId, double smallFlowEdgesThreshold, GeometryList& result)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            const auto edgesCrossingSmallFlowEdges = meshInstances[meshKernelId]->GetEdgesCrossingSmallFlowEdges(smallFlowEdgesThreshold);
            const auto smallFlowEdgeCenters = meshInstances[meshKernelId]->GetFlowEdgesCenters(edgesCrossingSmallFlowEdges);

            ConvertPointVectorToGeometryList(smallFlowEdgeCenters, result);
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_error(const char*& error_message)
    {
        int exitCode = Success;
        error_message = exceptionMessage;
        return exitCode;
    }

    MKERNEL_API int mkernel_get_obtuse_triangles_count(int meshKernelId, int& numObtuseTriangles)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            const auto obtuseTriangles = meshInstances[meshKernelId]->GetObtuseTrianglesCenters();

            numObtuseTriangles = static_cast<int>(obtuseTriangles.size());
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_get_obtuse_triangles(int meshKernelId, GeometryList& result)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            const auto obtuseTriangles = meshInstances[meshKernelId]->GetObtuseTrianglesCenters();

            ConvertPointVectorToGeometryList(obtuseTriangles, result);
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
        }
        return exitCode;
    }

    MKERNEL_API int mkernel_delete_small_flow_edges(int meshKernelId, double smallFlowEdgesThreshold, double minFractionalAreaTriangles)
    {
        int exitCode = Success;
        try
        {
            if (meshKernelId >= meshInstances.size())
            {
                throw std::invalid_argument("MeshKernel: The selected mesh does not exist.");
            }

            meshInstances[meshKernelId]->DeleteSmallFlowEdges(smallFlowEdgesThreshold);
            meshInstances[meshKernelId]->DeleteSmallTrianglesAtBoundaries(minFractionalAreaTriangles);
        }
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
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
            const auto mesh = std::make_shared<meshkernel::Mesh>(edges, nodes, projection);

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
                                                         static_cast<meshkernel::InterpolationLocation>(locationType),
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
        catch (const std::exception& e)
        {
            strcpy_s(exceptionMessage, sizeof exceptionMessage, e.what());
            exitCode |= Exception;
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
        const auto location = static_cast<meshkernel::InterpolationLocation>(locationType);
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

        return 0;
    }

} // namespace meshkernelapi
