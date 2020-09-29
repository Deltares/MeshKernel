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

#include "MeshKernel.hpp"
#include "Mesh.hpp"
#include "Smoother.hpp"
#include "Orthogonalizer.hpp"
#include "OrthogonalizationAndSmoothing.hpp"
#include "CurvilinearGridFromSplines.hpp"
#include "CurvilinearGridFromSplinesTransfinite.hpp"
#include "CurvilinearGridFromPolygon.hpp"
#include "FlipEdges.hpp"
#include "CurvilinearGrid.hpp"
#include "Splines.hpp"
#include "Entities.hpp"
#include "MeshRefinement.hpp"
#include "LandBoundaries.hpp"
#include "Constants.cpp"

// The vector containing the mesh instances 
static std::vector<std::shared_ptr<meshkernel::Mesh>> meshInstances;

// For interactivity
static std::map<int, std::shared_ptr<meshkernel::OrthogonalizationAndSmoothing>> orthogonalizationInstances;
static std::map<int, std::shared_ptr<meshkernel::CurvilinearGridFromSplines>> curvilinearGridFromSplinesInstances;

namespace meshkernelapi
{
    static bool ConvertGeometryListNativeToPointVector(const GeometryListNative& geometryListIn, std::vector<meshkernel::Point>& result)
    {
        if (geometryListIn.numberOfCoordinates == 0)
        {
            // empty polygon
            return true;
        }
        result.resize(geometryListIn.numberOfCoordinates);

        for (int i = 0; i < geometryListIn.numberOfCoordinates; i++)
        {
            result[i] = { geometryListIn.xCoordinates[i] , geometryListIn.yCoordinates[i] };
        }
        return true;
    }

    static bool ConvertGeometryListNativeToSampleVector(const GeometryListNative& geometryListIn, std::vector<meshkernel::Sample>& result)
    {
        if (geometryListIn.numberOfCoordinates == 0)
        {
            // empty samples
            return true;
        }
        result.resize(geometryListIn.numberOfCoordinates);

        for (int i = 0; i < geometryListIn.numberOfCoordinates; i++)
        {
            result[i] = { geometryListIn.xCoordinates[i] , geometryListIn.yCoordinates[i], geometryListIn.zCoordinates[i] };
        }
        return true;
    }

    static bool ConvertPointVectorToGeometryListNative(std::vector<meshkernel::Point> pointVector, GeometryListNative& result)
    {
        // invalid memory allocation
        if (pointVector.size() < result.numberOfCoordinates)
        {
            return false;
        }

        for (int i = 0; i < result.numberOfCoordinates; i++)
        {
            result.xCoordinates[i] = pointVector[i].x;
            result.yCoordinates[i] = pointVector[i].y;
        }
        return true;
    }

    static bool SetSplines(const GeometryListNative& geometryListIn, meshkernel::Splines& spline)
    {
        if (geometryListIn.numberOfCoordinates == 0)
        {
            return false;
        }

        std::vector<meshkernel::Point> splineCornerPoints;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, splineCornerPoints);
        if (!successful)
        {
            return false;
        }

        std::vector<std::vector<size_t>> indexes(splineCornerPoints.size(), std::vector<size_t>(2,0));
        int pos = FindIndexes(splineCornerPoints, 0, splineCornerPoints.size(), meshkernel::doubleMissingValue, indexes);
        indexes.resize(pos);

        for (auto i = 0; i < indexes.size(); i++)
        {
            int size = indexes[i][1] - indexes[i][0] + 1;
            if (size > 0)
            {
                spline.AddSpline(splineCornerPoints, indexes[i][0], size);
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
        if(meshGeometryDimensions.numface > 0)
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

    MKERNEL_API int mkernel_new_mesh(int& meshKernelId)
    {   
        meshKernelId = int(meshInstances.size());
        meshInstances.push_back(std::make_shared<meshkernel::Mesh>());
        return 0;
    };

    MKERNEL_API int mkernel_deallocate_state(int meshKernelId)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        meshInstances.pop_back();
        return 0;
    }

    MKERNEL_API int mkernel_delete_mesh(int meshKernelId, GeometryListNative& geometryListIn, int deletionOption, bool invertDeletion)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        if (meshInstances[meshKernelId]->GetNumNodes() <= 0)
        {
            return 0;
        }

        std::vector<meshkernel::Point> polygonPoints;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, polygonPoints);
        if (!successful)
        {
            return -1;
        }

        meshkernel::Polygons polygon(polygonPoints, meshInstances[meshKernelId]->m_projection);
        meshInstances[meshKernelId]->DeleteMesh(polygon, deletionOption, invertDeletion);

        return 0;
    }


    MKERNEL_API int mkernel_set_state(int meshKernelId, const MeshGeometryDimensions& meshGeometryDimensions, const MeshGeometry& meshGeometry, bool isGeographic)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        std::vector<meshkernel::Edge> edges(meshGeometryDimensions.numedge);
        int ei = 0;
        for (int e = 0; e < edges.size(); e++)
        {
            edges[e].first = meshGeometry.edge_nodes[ei];
            ei++;
            edges[e].second = meshGeometry.edge_nodes[ei];
            ei++;
        }

        std::vector<meshkernel::Point> nodes(meshGeometryDimensions.numnode);
        for (int n = 0; n < nodes.size(); n++)
        {
            nodes[n].x = meshGeometry.nodex[n];
            nodes[n].y = meshGeometry.nodey[n];
        }

        // spherical or cartesian
        if (isGeographic)
        {
            meshInstances[meshKernelId]->Set(edges, nodes, meshkernel::Projections::spherical);
        }
        else
        {
            meshInstances[meshKernelId]->Set(edges, nodes, meshkernel::Projections::cartesian);
        }

        return 0;
    }

    MKERNEL_API int mkernel_get_mesh(int meshKernelId, MeshGeometryDimensions& meshGeometryDimensions, MeshGeometry& meshGeometry)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        meshInstances[meshKernelId]->SetFlatCopies(meshkernel::Mesh::AdministrationOptions::AdministrateMeshEdges);

        bool successful = SetMeshGeometry(meshKernelId, meshGeometryDimensions, meshGeometry);
        if (!successful)
        {
            return -1;
        }

        return 0;
    }

    MKERNEL_API int mkernel_find_faces(int meshKernelId, MeshGeometryDimensions& meshGeometryDimensions, MeshGeometry& meshGeometry)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        meshInstances[meshKernelId]->SetFlatCopies(meshkernel::Mesh::AdministrationOptions::AdministrateMeshEdgesAndFaces);

        bool successful = SetMeshGeometry(meshKernelId,meshGeometryDimensions, meshGeometry);
        if(!successful)
        {
            return -1;
        }

        return 0;
    }

    MKERNEL_API int mkernel_orthogonalize( int meshKernelId, 
                                           int isTriangulationRequired, 
                                           int isAccountingForLandBoundariesRequired, 
                                           int projectToLandBoundaryOption,
                                           const OrthogonalizationParametersNative& orthogonalizationParametersNative, 
                                           const GeometryListNative& geometryListNativePolygon,
                                           const GeometryListNative& geometryListNativeLandBoundaries)
    {

        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        if (meshInstances[meshKernelId]->GetNumNodes() <= 0)
        {
            return 0;
        }
        
        // build enclosing polygon
        std::vector<meshkernel::Point> nodes(geometryListNativePolygon.numberOfCoordinates);
        for (int i = 0; i < geometryListNativePolygon.numberOfCoordinates; i++)
        {
            nodes[i].x = geometryListNativePolygon.xCoordinates[i];
            nodes[i].y = geometryListNativePolygon.yCoordinates[i];
        }

        auto polygon = std::make_shared<meshkernel::Polygons>(nodes, meshInstances[meshKernelId]->m_projection);

        // build land boundary
        std::vector<meshkernel::Point> landBoundaries(geometryListNativeLandBoundaries.numberOfCoordinates);
        for (int i = 0; i < geometryListNativeLandBoundaries.numberOfCoordinates; i++)
        {
            landBoundaries[i].x = geometryListNativeLandBoundaries.xCoordinates[i];
            landBoundaries[i].y = geometryListNativeLandBoundaries.yCoordinates[i];
        }

        auto orthogonalizer = std::make_shared<meshkernel::Orthogonalizer>(meshInstances[meshKernelId]);
        auto smoother = std::make_shared<meshkernel::Smoother>(meshInstances[meshKernelId]);
        auto landBoundary = std::make_shared < meshkernel::LandBoundaries>(landBoundaries, meshInstances[meshKernelId], polygon);

        meshkernel::OrthogonalizationAndSmoothing ortogonalization(meshInstances[meshKernelId],
                                                                   smoother,
                                                                   orthogonalizer,
                                                                   polygon,
                                                                   landBoundary,
                                                                   isTriangulationRequired,
                                                                   isAccountingForLandBoundariesRequired,
                                                                   projectToLandBoundaryOption,
                                                                   orthogonalizationParametersNative);
        bool successful = ortogonalization.Initialize();
        if (!successful)
        {
            return -1;
        }

        successful = ortogonalization.Compute();
        if (!successful)
        {
            return -1;
        }
    }

    MKERNEL_API int mkernel_orthogonalize_initialize(int meshKernelId,
        int isTriangulationRequired,
        int isAccountingForLandBoundariesRequired,
        int projectToLandBoundaryOption,
        OrthogonalizationParametersNative& orthogonalizationParametersNative,
        GeometryListNative& geometryListNativePolygon,
        GeometryListNative& geometryListNativeLandBoundaries)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        if (meshInstances[meshKernelId]->GetNumNodes() <= 0)
        {
            return 0;
        }

        // build enclosing polygon
        std::vector<meshkernel::Point> nodes(geometryListNativePolygon.numberOfCoordinates);
        for (int i = 0; i < geometryListNativePolygon.numberOfCoordinates; i++)
        {
            nodes[i].x = geometryListNativePolygon.xCoordinates[i];
            nodes[i].y = geometryListNativePolygon.yCoordinates[i];
        }

        // build land boundary
        std::vector<meshkernel::Point> landBoundaries(geometryListNativeLandBoundaries.numberOfCoordinates);
        for (int i = 0; i < geometryListNativeLandBoundaries.numberOfCoordinates; i++)
        {
            landBoundaries[i].x = geometryListNativeLandBoundaries.xCoordinates[i];
            landBoundaries[i].y = geometryListNativeLandBoundaries.yCoordinates[i];
        }


        auto orthogonalizer = std::make_shared<meshkernel::Orthogonalizer>(meshInstances[meshKernelId]);
        auto smoother = std::make_shared<meshkernel::Smoother>(meshInstances[meshKernelId]);
        auto polygon = std::make_shared<meshkernel::Polygons>(nodes, meshInstances[meshKernelId]->m_projection);
        auto landBoundary = std::make_shared<meshkernel::LandBoundaries>(landBoundaries, meshInstances[meshKernelId], polygon);

        auto orthogonalizationInstance = std::make_shared<meshkernel::OrthogonalizationAndSmoothing>( meshInstances[meshKernelId],
                                                                                                      smoother,
                                                                                                      orthogonalizer,
                                                                                                      polygon,
                                                                                                      landBoundary,
                                                                                                      isTriangulationRequired,
                                                                                                      isAccountingForLandBoundariesRequired,
                                                                                                      projectToLandBoundaryOption,
                                                                                                      orthogonalizationParametersNative );
        bool successful = orthogonalizationInstance->Initialize();
        if (!successful)
        {
            return -1;
        }

        orthogonalizationInstances.insert({ meshKernelId, orthogonalizationInstance });

        return 0;
    }

    MKERNEL_API int mkernel_orthogonalize_prepare_outer_iteration(int meshKernelId)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        if (meshInstances[meshKernelId]->GetNumNodes() <= 0)
        {
            return 0;
        }

        bool status = orthogonalizationInstances[meshKernelId]->PrapareOuterIteration();
        return status == true ? 0 : 1;
    }

    MKERNEL_API int mkernel_orthogonalize_inner_iteration(int meshKernelId)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        if (meshInstances[meshKernelId]->GetNumNodes() <= 0)
        {
            return 0;
        }
        const bool status = orthogonalizationInstances[meshKernelId]->InnerIteration();
        return status == true ? 0 : 1;
    }

    MKERNEL_API int mkernel_orthogonalize_finalize_outer_iteration(int meshKernelId)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        if (meshInstances[meshKernelId]->GetNumNodes() <= 0)
        {
            return 0;
        }
        const bool status = orthogonalizationInstances[meshKernelId]->FinalizeOuterIteration();
        return status == true ? 0 : 1;
    }

    MKERNEL_API int mkernel_orthogonalize_delete(int meshKernelId)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        if (meshInstances[meshKernelId]->GetNumNodes() <= 0)
        {
            return 0;
        }

        const auto returnValue = orthogonalizationInstances.erase(meshKernelId);
        return returnValue == 1 ? 0 : 1;
    }

    MKERNEL_API int mkernel_get_orthogonality(int meshKernelId, GeometryListNative& geometryList)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        if (meshInstances[meshKernelId]->GetNumNodes() <= 0)
        {
            return 0;
        }
        const bool status = meshInstances[meshKernelId]->GetOrthogonality(geometryList.zCoordinates);
        return status == true ? 0 : 1;
    }

    MKERNEL_API int mkernel_get_smoothness(int meshKernelId, GeometryListNative& geometryList)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        if (meshInstances[meshKernelId]->GetNumNodes() <= 0)
        {
            return 0;
        }

        const bool status = meshInstances[meshKernelId]->GetSmoothness(geometryList.zCoordinates);
        return status == true ? 0 : 1;
    }

    MKERNEL_API int mkernel_get_splines( const GeometryListNative& geometryListIn, 
                                         GeometryListNative& geometry_list_out, 
                                         int number_of_points_between_vertices )
    {

        if (geometryListIn.numberOfCoordinates == 0)
        {
            return -1;
        }

        std::vector<meshkernel::Point> splines(geometryListIn.numberOfCoordinates);
        for (int i = 0; i < geometryListIn.numberOfCoordinates; i++)
        {
            splines[i].x = geometryListIn.xCoordinates[i];
            splines[i].y = geometryListIn.yCoordinates[i];
        }

        std::vector<std::vector<size_t>> indexes(geometryListIn.numberOfCoordinates, std::vector<size_t>(2));
        int numSplines = FindIndexes(splines, 0, splines.size(), meshkernel::doubleMissingValue, indexes);
        std::vector<meshkernel::Point> coordinatesDerivatives(geometryListIn.numberOfCoordinates);

        int index = 0;
        for (int s = 0; s < numSplines; s++)
        {            
            std::vector<meshkernel::Point> coordinates(splines.begin() + indexes[s][0], splines.begin() + indexes[s][1] + 1);
            int numNodes = indexes[s][1] - indexes[s][0] + 1;
            meshkernel::Splines::SecondOrderDerivative(coordinates, numNodes, coordinatesDerivatives);

            for (int n = 0; n < numNodes - 1; n++)
            {
                for (int p = 0; p <= number_of_points_between_vertices; p++)
                {

                    double pointAdimensionalCoordinate = n + double(p)/ double(number_of_points_between_vertices);
                    meshkernel::Point pointCoordinate{meshkernel::doubleMissingValue,meshkernel::doubleMissingValue};
                    bool successful = InterpolateSplinePoint(coordinates, coordinatesDerivatives, pointAdimensionalCoordinate, pointCoordinate);
                    if (!successful) 
                    {
                        break;
                    }

                    geometry_list_out.xCoordinates[index] = pointCoordinate.x;
                    geometry_list_out.yCoordinates[index] = pointCoordinate.y;
                    geometry_list_out.zCoordinates[index] = meshkernel::doubleMissingValue;
                    index++;
                }
            }

            geometry_list_out.xCoordinates[index] = meshkernel::doubleMissingValue;
            geometry_list_out.yCoordinates[index] = meshkernel::doubleMissingValue;
            geometry_list_out.zCoordinates[index] = meshkernel::doubleMissingValue;
            index++;
        }

        geometry_list_out.numberOfCoordinates = index - 1;
        return 0;
    }

    MKERNEL_API int mkernel_make_mesh(int meshKernelId, const MakeGridParametersNative& makeGridParameters, const GeometryListNative& geometryListNative)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        std::vector<meshkernel::Point> result;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListNative, result);
        if (!successful)
        {
            return -1;
        }

        meshkernel::Polygons polygon(result, meshInstances[meshKernelId]->m_projection);

        meshkernel::Mesh mesh;
        successful = mesh.MakeMesh(makeGridParameters, polygon);
        if (!successful)
        {
            return -1;
        }

        *meshInstances[meshKernelId] += mesh;

        return 0;
    }

    MKERNEL_API int mkernel_make_mesh_from_polygon(int meshKernelId, const GeometryListNative& disposableGeometryListIn)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        std::vector<meshkernel::Point> result;
        bool successful = ConvertGeometryListNativeToPointVector(disposableGeometryListIn, result);
        if (!successful)
        {
            return -1;
        }

        meshkernel::Polygons polygon(result, meshInstances[meshKernelId]->m_projection);

        std::vector<std::vector<meshkernel::Point>> generatedPoints;
        successful = polygon.CreatePointsInPolygons(generatedPoints);
        if (!successful)
        {
            return -1;
        }

        meshkernel::Mesh mesh(generatedPoints[0], polygon, meshInstances[meshKernelId]->m_projection);
        *meshInstances[meshKernelId] += mesh;
        return 0;
    }

    MKERNEL_API int mkernel_make_mesh_from_samples(int meshKernelId, GeometryListNative& geometryListNative)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        std::vector<meshkernel::Point> samplePoints;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListNative, samplePoints);
        if (!successful)
        {
            return -1;
        }
        meshkernel::Polygons polygon;
        meshkernel::Mesh mesh(samplePoints, polygon, meshInstances[meshKernelId]->m_projection);
        *meshInstances[meshKernelId] += mesh;
        return 0;
    }

    MKERNEL_API int mkernel_copy_mesh_boundaries_to_polygon(int meshKernelId, GeometryListNative& geometryListNative)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        meshkernel::Polygons polygon;
        std::vector<meshkernel::Point> meshBoundaryPolygon;

        int numNodesBoundaryPolygons;
        polygon.MeshBoundaryToPolygon(*meshInstances[meshKernelId], meshBoundaryPolygon, numNodesBoundaryPolygons);

        bool successful = ConvertPointVectorToGeometryListNative(meshBoundaryPolygon, geometryListNative);
        if (!successful)
        {
            return -1;
        }

        return 0;
    }

    MKERNEL_API int mkernel_copy_mesh_boundaries_to_polygon_count_vertices(int meshKernelId, int& numberOfPolygonVertices)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        meshkernel::Polygons polygon;
        std::vector<meshkernel::Point> meshBoundaryPolygon;

        polygon.MeshBoundaryToPolygon(*meshInstances[meshKernelId], meshBoundaryPolygon, numberOfPolygonVertices);

        return 0;
    }

    MKERNEL_API int mkernel_refine_polygon(int meshKernelId, const GeometryListNative& geometryListIn, int firstIndex, int secondIndex, double distance, GeometryListNative& geometryListOut)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        std::vector<meshkernel::Point> polygonPoints;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, polygonPoints);
        if (!successful)
        {
            return -1;
        }

        meshkernel::Polygons polygon(polygonPoints, meshInstances[meshKernelId]->m_projection);

        std::vector<meshkernel::Point> refinedPolygon;
        successful = polygon.RefinePolygonPart(firstIndex, secondIndex, distance, refinedPolygon);
        if (!successful)
        {
            return -1;
        }

        successful = ConvertPointVectorToGeometryListNative(refinedPolygon, geometryListOut);
        if (!successful)
        {
            return -1;
        }

        return 0;
    }

    MKERNEL_API int mkernel_refine_polygon_count(int meshKernelId, GeometryListNative& geometryListIn, int firstIndex, int secondIndex, double distance, int& numberOfPolygonVertices)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        std::vector<meshkernel::Point> polygonPoints;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, polygonPoints);
        if (!successful)
        {
            return -1;
        }

        meshkernel::Polygons polygon(polygonPoints, meshInstances[meshKernelId]->m_projection);

        std::vector<meshkernel::Point> refinedPolygon;
        successful = polygon.RefinePolygonPart(firstIndex, secondIndex, distance, refinedPolygon);
        if (!successful)
        {
            return -1;
        }

        numberOfPolygonVertices = int(refinedPolygon.size());

        return 0;
    }


    MKERNEL_API int mkernel_merge_nodes(int meshKernelId, GeometryListNative& geometryListIn) 
    {

        if (meshKernelId >= meshInstances.size()) 
        {
            return 0;
        }

        std::vector<meshkernel::Point> polygonPoints;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, polygonPoints);
        if (!successful)
        {
            return -1;
        }

        meshkernel::Polygons polygon(polygonPoints, meshInstances[meshKernelId]->m_projection);

        successful = meshInstances[meshKernelId]->MergeNodesInPolygon(polygon);
        if (!successful)
        {
            return -1;
        }
        return 0;
    }

    MKERNEL_API int mkernel_merge_two_nodes(int meshKernelId, int startNode, int endNode)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return 0;
        }

        bool successful = meshInstances[meshKernelId]->MergeTwoNodes(startNode, endNode);
        if (!successful)
        {
            return -1;
        }

        return 0;
    }
    // TO-DO: remove numberOfMeshVertices
    MKERNEL_API int mkernel_nodes_in_polygons(int meshKernelId, GeometryListNative& geometryListIn, int inside, int numberOfMeshVertices, int** selectedVertices)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return 0;
        }

        std::vector<meshkernel::Point> polygonPoints;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, polygonPoints);
        if (!successful)
        {
            return -1;
        }

        meshkernel::Polygons polygon(polygonPoints, meshInstances[meshKernelId]->m_projection);

        bool selectInside = inside == 1 ? true : false;
        successful = meshInstances[meshKernelId]->MaskNodesInPolygons(polygon, selectInside);
        if (!successful)
        {
            return -1;
        }

        int index = 0;
        for (int i = 0; i < meshInstances[meshKernelId]->GetNumNodes() ; ++i)
        {
            if (meshInstances[meshKernelId]->m_nodeMask[i]>0) 
            {
                (*selectedVertices)[index] = i;
                index++;
            }
        }

        return 0;
    }

    MKERNEL_API int mkernel_count_nodes_in_polygons(int meshKernelId, GeometryListNative& geometryListIn, int inside, int& numberOfMeshVertices)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return 0;
        }

        std::vector<meshkernel::Point> polygonPoints;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, polygonPoints);
        if (!successful)
        {
            return -1;
        }

        meshkernel::Polygons polygon(polygonPoints, meshInstances[meshKernelId]->m_projection);

        bool selectInside = inside == 1 ? true : false;
        successful = meshInstances[meshKernelId]->MaskNodesInPolygons(polygon, selectInside);
        if (!successful)
        {
            return -1;
        }

        numberOfMeshVertices = 0;
        for (auto i = 0; i < meshInstances[meshKernelId]->GetNumNodes(); ++i)
        {
            if (meshInstances[meshKernelId]->m_nodeMask[i] > 0)
            {
                numberOfMeshVertices++;
            }
        }

        return 0;
    }

    MKERNEL_API int mkernel_insert_edge(int meshKernelId, int startNode, int endNode, int& new_edge_index)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return 0;
        }

        bool successful = meshInstances[meshKernelId]->ConnectNodes(startNode, endNode, new_edge_index);
        if (!successful)
        {
            return -1;
        }


        return 0;
    }

    MKERNEL_API int mkernel_insert_node(int meshKernelId, double xCoordinate, double yCoordinate, double zCoordinate, int& nodeIndex)
    {
        if (meshKernelId >= meshInstances.size())
        {
            //create a valid instance, by default cartesian
            *meshInstances[meshKernelId] = meshkernel::Mesh();
            meshInstances[meshKernelId]->m_projection = meshkernel::Projections::cartesian;
        }

        meshkernel::Point newNode{ xCoordinate, yCoordinate };
        bool successful = meshInstances[meshKernelId]->InsertNode(newNode, nodeIndex);
        if (!successful)
        {
            return -1;
        }

        return 0;
    }


    MKERNEL_API int mkernel_delete_node(int meshKernelId, int nodeIndex)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return 0;
        }

        bool successful = meshInstances[meshKernelId]->DeleteNode(nodeIndex);
        if (!successful)
        {
            return -1;
        }

        return 0;
    }

    MKERNEL_API int mkernel_move_node(int meshKernelId, GeometryListNative& geometryListIn, int nodeIndex)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return 0;
        }

        std::vector<meshkernel::Point> newPoint;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, newPoint);
        if (!successful || newPoint.size() != 1)
        {
            return -1;
        }

        successful = meshInstances[meshKernelId]->MoveNode(newPoint[0],nodeIndex);
        if (!successful)
        {
            return -1;
        }

        return 0;
    }

    MKERNEL_API int mkernel_delete_edge(int meshKernelId, GeometryListNative& geometryListIn, double searchRadius)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return 0;
        }

        std::vector<meshkernel::Point> newPoint;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, newPoint);
        if (!successful || newPoint.size() != 1)
        {
            return -1;
        }

        auto edgeIndex = meshInstances[meshKernelId]->FindEdgeCloseToAPoint(newPoint[0], searchRadius);
        successful = meshInstances[meshKernelId]->DeleteEdge(edgeIndex);
        if (!successful)
        {
            return -1;
        }

        return 0;
    }

    MKERNEL_API int mkernel_find_edge(int meshKernelId, GeometryListNative& geometryListIn, double searchRadius, int& edgeIndex)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return 0;
        }

        std::vector<meshkernel::Point> newPoint;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, newPoint);
        if (!successful || newPoint.size() != 1)
        {
            return -1;
        }

        edgeIndex = meshInstances[meshKernelId]->FindEdgeCloseToAPoint(newPoint[0], searchRadius);

        return 0;
    }

    MKERNEL_API int mkernel_offsetted_polygon(int meshKernelId, GeometryListNative& geometryListIn, bool innerAndOuter, double distance, GeometryListNative& geometryListOut)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return 0;
        }

        std::vector<meshkernel::Point> polygonPoints;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, polygonPoints);
        if (!successful)
        {
            return -1;
        }

        meshkernel::Polygons polygon(polygonPoints, meshInstances[meshKernelId]->m_projection);

        meshkernel::Polygons newPolygon;
        successful = polygon.OffsetCopy(distance, innerAndOuter, newPolygon);
        if (!successful)
        {
            return -1;
        }

        successful = ConvertPointVectorToGeometryListNative(newPolygon.m_nodes, geometryListOut);
        if (!successful)
        {
            return -1;
        }

        return 0;
    }

    MKERNEL_API int mkernel_offsetted_polygon_count(int meshKernelId, GeometryListNative& geometryListIn, bool innerAndOuter, double distance, int& numberOfPolygonVertices)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return 0;
        }

        std::vector<meshkernel::Point> polygonPoints;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, polygonPoints);
        if (!successful)
        {
            return -1;
        }

        meshkernel::Polygons polygon(polygonPoints, meshInstances[meshKernelId]->m_projection);

        meshkernel::Polygons newPolygon;
        successful = polygon.OffsetCopy(distance, innerAndOuter, newPolygon);
        if (!successful)
        {
            return -1;
        }

        numberOfPolygonVertices = newPolygon.GetNumNodes();
    
        return 0;
    }

    MKERNEL_API int mkernel_refine_mesh_based_on_samples(int meshKernelId, GeometryListNative& geometryListIn, InterpolationParametersNative& interpolationParametersNative, SampleRefineParametersNative& sampleRefineParametersNative)
    {
        if (meshKernelId >= meshInstances.size() || meshInstances[meshKernelId]->GetNumNodes() <= 0)
        {
            return 0;
        }

        std::vector<meshkernel::Sample> samples;
        bool successful = ConvertGeometryListNativeToSampleVector(geometryListIn, samples);
        if (!successful)
        {
            return -1;
        }

        if (samples.empty())
        {
            return 0;
        }

        meshkernel::MeshRefinement meshRefinement(meshInstances[meshKernelId]);

        // polygon could be passed as api parameter
        meshkernel::Polygons polygon;

        successful = meshRefinement.Refine(samples, polygon, sampleRefineParametersNative, interpolationParametersNative);
        if (!successful)
        {
            return -1;
        }

        return 0;
    }

    MKERNEL_API int mkernel_refine_mesh_based_on_polygon(int meshKernelId, GeometryListNative& geometryListNative, InterpolationParametersNative& interpolationParametersNative)
    {
        
        if (meshKernelId >= meshInstances.size() || meshInstances[meshKernelId]->GetNumNodes()<=0)
        {
            return 0;
        }


        std::vector<meshkernel::Point> points;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListNative, points);
        if (!successful)
        {
            return -1;
        }

        if (points.empty())
        {
            return 0;
        }


        meshkernel::Polygons polygon(points, meshInstances[meshKernelId]->m_projection);

        meshkernel::MeshRefinement meshRefinement(meshInstances[meshKernelId]);

        // polygon could be passed as api parameter
        std::vector<meshkernel::Sample> samples; 
        SampleRefineParametersNative sampleRefineParametersNative;
        successful = meshRefinement.Refine(samples, polygon, sampleRefineParametersNative, interpolationParametersNative);
        if (!successful)
        {
            return -1;
        }

        return 0;
    }

    MKERNEL_API int mkernel_get_node_index(int meshKernelId, GeometryListNative& geometryListIn, double searchRadius, int& nodeIndex)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        if (meshInstances[meshKernelId]->GetNumNodes() <= 0)
        {
            return -1;
        }

        std::vector<meshkernel::Point> polygonPoints;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, polygonPoints);
        if (!successful || polygonPoints.empty())
        {
            return -1;
        }

        successful = meshInstances[meshKernelId]->GetNodeIndex(polygonPoints[0], searchRadius, nodeIndex);
        if (!successful)
        {
            return -1;
        }

        return 0;
    }

    MKERNEL_API int mkernel_get_node_coordinate(int meshKernelId, GeometryListNative& geometryListIn, double searchRadius, GeometryListNative& geometryListOut)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        if (meshInstances[meshKernelId]->GetNumNodes() <= 0)
        {
            return -1;
        }

        if (geometryListOut.numberOfCoordinates <= 0)
        {
            return -1;
        }

        std::vector<meshkernel::Point> polygonPoints;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, polygonPoints);
        if (!successful || polygonPoints.empty())
        {
            return -1;
        }

        int nodeIndex = -1;
        successful = meshInstances[meshKernelId]->GetNodeIndex(polygonPoints[0], searchRadius, nodeIndex);
        if (!successful || nodeIndex < 0)
        {
            return -1;
        }

        // Set the node coordinate
        auto node = meshInstances[meshKernelId]->m_nodes[nodeIndex];
        std::vector<meshkernel::Point> pointVector;
        pointVector.push_back(node);
        ConvertPointVectorToGeometryListNative(pointVector, geometryListOut);

        if (!successful)
        {
            return -1;
        }

        return 0;
    }


    MKERNEL_API int mkernel_curvilinear_mesh_from_splines_ortho(int meshKernelId,
                                                                const GeometryListNative& geometryListIn,
                                                                const CurvilinearParametersNative& curvilinearParameters,
                                                                const SplinesToCurvilinearParametersNative& splineToCurvilinearParameters)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        // use the default constructor, no instance present
        auto spline = std::make_shared<meshkernel::Splines>(meshInstances[meshKernelId]->m_projection);
        bool successful = SetSplines(geometryListIn, *spline);
        if(!successful)
        {
            return -1;
        }

        meshkernel::CurvilinearGridFromSplines curvilinearGridFromSplines(spline, curvilinearParameters, splineToCurvilinearParameters);

        meshkernel::CurvilinearGrid curvilinearGrid;
        curvilinearGridFromSplines.Compute(curvilinearGrid);
        *meshInstances[meshKernelId] += meshkernel::Mesh(curvilinearGrid, meshInstances[meshKernelId]->m_projection);

        return 0;
    }


    MKERNEL_API int mkernel_curvilinear_mesh_from_splines_ortho_initialize(int meshKernelId, const GeometryListNative& geometryListNative, const CurvilinearParametersNative& curvilinearParametersNative, const SplinesToCurvilinearParametersNative& splinesToCurvilinearParametersNative)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return 0;
        }

        auto spline = std::make_shared<meshkernel::Splines>(meshInstances[meshKernelId]->m_projection);
        bool successful = SetSplines(geometryListNative, *spline);
        if (!successful)
        {
            return -1;
        }
        auto curvilinearGridFromSplines =std::make_shared<meshkernel::CurvilinearGridFromSplines>(spline, curvilinearParametersNative, splinesToCurvilinearParametersNative);

        curvilinearGridFromSplinesInstances.insert({ meshKernelId, curvilinearGridFromSplines });
 
        successful = curvilinearGridFromSplinesInstances[meshKernelId]->Initialize();
        if (!successful)
        {
            return -1;
        }
        return 0;
    }

    MKERNEL_API int mkernel_curvilinear_mesh_from_splines_iteration(int meshKernelId, int layer)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return 0;
        }

        const bool successful = curvilinearGridFromSplinesInstances[meshKernelId]->Iterate(layer);
        return successful ? 0 : 1;
    }


    MKERNEL_API int mkernel_curvilinear_mesh_from_splines_ortho_refresh_mesh(int meshKernelId)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return 0;
        }

        meshkernel::CurvilinearGrid curvilinearGrid;
        const bool successful = curvilinearGridFromSplinesInstances[meshKernelId]->ComputeCurvilinearGrid(curvilinearGrid);

        *meshInstances[meshKernelId] += meshkernel::Mesh(curvilinearGrid, meshInstances[meshKernelId]->m_projection);

        return successful ? 0 : 1;
    }

    MKERNEL_API int mkernel_curvilinear_mesh_from_splines_ortho_delete(int meshKernelId)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return 0;
        }
        const auto returnValue = curvilinearGridFromSplinesInstances.erase(meshKernelId);
        return returnValue == 1 ? 0 : 1;
    }

    MKERNEL_API int mkernel_points_in_polygon(int meshKernelId, GeometryListNative& polygonNative, GeometryListNative& pointsNative, GeometryListNative& selectedPointsNative)
    {
    
        std::vector<meshkernel::Point> polygonNodes;
        bool successful = ConvertGeometryListNativeToPointVector(polygonNative, polygonNodes);
        if (!successful)
        {
            return -1;
        }

        std::vector<meshkernel::Point> points;
        successful = ConvertGeometryListNativeToPointVector(pointsNative, points);
        if (!successful || points.empty() || points.size()!= selectedPointsNative.numberOfCoordinates)
        {
            return -1;
        }

        
        meshkernel::Polygons polygon(polygonNodes, meshInstances[meshKernelId]->m_projection);

        for (int i = 0; i < points.size(); i++)
        {
            selectedPointsNative.zCoordinates[i]= polygon.IsPointInPolygon(points[i], 0)? 1.0 : 0.0;
        }
        
        return successful ? 0 : 1;
    }

    MKERNEL_API int mkernel_flip_edges(int meshKernelId, 
                                       int isTriangulationRequired, 
                                       int isAccountingForLandBoundariesRequired, 
                                       int projectToLandBoundaryRequired)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        //set landboundaries
        auto polygon = std::make_shared<meshkernel::Polygons>();

        std::vector<meshkernel::Point> landBoundary;
        auto landBoundaries = std::make_shared<meshkernel::LandBoundaries>(landBoundary, meshInstances[meshKernelId], polygon);

        bool triangulateFaces = isTriangulationRequired == 0 ? false : true;
        bool projectToLandBoundary = projectToLandBoundaryRequired == 0 ? false : true;
        meshkernel::FlipEdges flipEdges(meshInstances[meshKernelId], landBoundaries, triangulateFaces, projectToLandBoundary);

        // compute Flip edges
        auto successful = flipEdges.Compute();

        return successful ? 0 : 1;
    }

    MKERNEL_API int mkernel_curvilinear_mesh_from_splines(int meshKernelId, GeometryListNative& geometryListNativeIn, CurvilinearParametersNative& curvilinearParametersNative)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        // Use the default constructor, no instance present
        auto spline = std::make_shared<meshkernel::Splines>(meshInstances[meshKernelId]->m_projection);
        bool successful = SetSplines(geometryListNativeIn, *spline);
        if (!successful)
        {
            return -1;
        }

        // Create algorithm and set the splines
        meshkernel::CurvilinearGridFromSplinesTransfinite curvilinearGridFromSplinesTransfinite(spline, curvilinearParametersNative);

        
        // Compute the curvilinear grid
        meshkernel::CurvilinearGrid curvilinearGrid;
        successful = curvilinearGridFromSplinesTransfinite.Compute(curvilinearGrid);
        if (!successful)
        {
            return -1;
        }

        // Transform and set mesh pointer
        *meshInstances[meshKernelId] += meshkernel::Mesh(curvilinearGrid, meshInstances[meshKernelId]->m_projection);

        return  0;
    }

    MKERNEL_API int mkernel_curvilinear_from_polygon( int meshKernelId,
                                                       GeometryListNative& polygonNative,
                                                       int firstNode,
                                                       int secondNode,
                                                       int thirdNode,
                                                       bool useFourthSide )
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        std::vector<meshkernel::Point> polygonPoints;
        bool successful = ConvertGeometryListNativeToPointVector(polygonNative, polygonPoints);
        if (!successful)
        {
            return -1;
        }

        auto polygon = std::make_shared<meshkernel::Polygons>(polygonPoints, meshInstances[meshKernelId]->m_projection);
        if (!successful)
        {
            return -1;
        }

        meshkernel::CurvilinearGrid curvilinearGrid;
        meshkernel::CurvilinearGridFromPolygon curvilinearGridFromPolygon(polygon);
        successful = curvilinearGridFromPolygon.Compute(firstNode, secondNode, thirdNode, useFourthSide, curvilinearGrid);

        if (!successful)
        {
            return -1;
        }

        // convert to curvilinear grid and add it to the current mesh
        *meshInstances[meshKernelId] += meshkernel::Mesh(curvilinearGrid, meshInstances[meshKernelId]->m_projection);

        return successful ? 0 : 1;
    }

    MKERNEL_API int mkernel_curvilinear_from_triangle( int meshKernelId,
                                                        GeometryListNative& polygonNative,
                                                        int firstNode,
                                                        int secondNode,
                                                        int thirdNode )
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }


        std::vector<meshkernel::Point> polygonPoints;
        bool successful = ConvertGeometryListNativeToPointVector(polygonNative, polygonPoints);
        if (!successful)
        {
            return -1;
        }

        auto polygon = std::make_shared<meshkernel::Polygons>(polygonPoints, meshInstances[meshKernelId]->m_projection);
        if (!successful)
        {
            return -1;
        }

        meshkernel::CurvilinearGrid curvilinearGrid;
        meshkernel::CurvilinearGridFromPolygon curvilinearGridFromPolygon(polygon);
        successful = curvilinearGridFromPolygon.Compute(firstNode, secondNode, thirdNode, curvilinearGrid);

        if (!successful)
        {
            return -1;
        }

        // convert to curvilinear grid and add it to the current mesh
        *meshInstances[meshKernelId] += meshkernel::Mesh(curvilinearGrid, meshInstances[meshKernelId]->m_projection);

        return successful ? 0 : 1;
    }

}
