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
#include "FlipEdges.hpp"
#include "CurvilinearGrid.hpp"
#include "Splines.hpp"
#include "Entities.hpp"
#include "MeshRefinement.hpp"
#include "LandBoundaries.hpp"
#include "Constants.cpp"

// The vector containing the mesh instances 
static std::vector<std::shared_ptr<MeshKernel::Mesh>> meshInstances;

// For supporting interactivity for orthogonalization and orthogonal curvilinear grid from splines, we need to save some instances
static std::map<int, MeshKernel::OrthogonalizationAndSmoothing> orthogonalizationInstances;
static std::map<int, MeshKernel::CurvilinearGridFromSplines> splineInstances;

namespace MeshKernelApi
{
    static bool ConvertGeometryListNativeToPointVector(GeometryListNative& geometryListIn, std::vector<MeshKernel::Point>& result)
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

    static bool ConvertGeometryListNativeToSampleVector(GeometryListNative& geometryListIn, std::vector<MeshKernel::Sample>& result)
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

    static bool ConvertPointVectorToGeometryListNative(std::vector<MeshKernel::Point> pointVector, GeometryListNative& result)
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

    static bool SetSplines(GeometryListNative& geometryListIn, MeshKernel::Splines& spline)
    {
        if (geometryListIn.numberOfCoordinates == 0)
        {
            return false;
        }

        std::vector<MeshKernel::Point> splineCornerPoints;
        bool success = ConvertGeometryListNativeToPointVector(geometryListIn, splineCornerPoints);
        if (!success)
        {
            return false;
        }

        std::vector<std::vector<int>> indexes(splineCornerPoints.size(), std::vector<int>(2,0));
        int pos = FindIndexes(splineCornerPoints, 0, splineCornerPoints.size(), MeshKernel::doubleMissingValue, indexes);
        indexes.resize(pos);

        for (int i = 0; i < indexes.size(); i++)
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

        meshGeometryDimensions.maxnumfacenodes = MeshKernel::maximumNumberOfNodesPerFace;
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

    MESHKERNEL_API int ggeo_new_mesh(int& meshKernelId)
    {
        int instanceSize = meshInstances.size();
        meshInstances.resize(instanceSize + 1);
        meshKernelId = instanceSize;
        meshInstances[meshKernelId] = std::make_shared<MeshKernel::Mesh>();
        return 0;
    };

    MESHKERNEL_API int ggeo_deallocate_state(int meshKernelId)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        meshInstances.pop_back();
        return 0;
    }

    MESHKERNEL_API int ggeo_delete_mesh(int meshKernelId, GeometryListNative& geometryListIn, int deletionOption, bool invertDeletion)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        if (meshInstances[meshKernelId]->GetNumNodes() <= 0)
        {
            return 0;
        }

        std::vector<MeshKernel::Point> polygonPoints;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, polygonPoints);
        if (!successful)
        {
            return -1;
        }

        MeshKernel::Polygons polygon;
        successful = polygon.Set(polygonPoints, meshInstances[meshKernelId]->m_projection);
        if (!successful)
        {
            return -1;
        }

        meshInstances[meshKernelId]->DeleteMesh(polygon, deletionOption, invertDeletion);

        return 0;
    }


    MESHKERNEL_API int ggeo_set_state(int meshKernelId, MeshGeometryDimensions& meshGeometryDimensions, MeshGeometry& meshGeometry, bool isGeographic)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        std::vector<MeshKernel::Edge> edges(meshGeometryDimensions.numedge);
        int ei = 0;
        for (int e = 0; e < edges.size(); e++)
        {
            edges[e].first = meshGeometry.edge_nodes[ei];
            ei++;
            edges[e].second = meshGeometry.edge_nodes[ei];
            ei++;
        }

        std::vector<MeshKernel::Point> nodes(meshGeometryDimensions.numnode);
        for (int n = 0; n < nodes.size(); n++)
        {
            nodes[n].x = meshGeometry.nodex[n];
            nodes[n].y = meshGeometry.nodey[n];
        }

        // spherical or cartesian
        if (isGeographic)
        {
            meshInstances[meshKernelId]->Set(edges, nodes, MeshKernel::Projections::spherical);
        }
        else
        {
            meshInstances[meshKernelId]->Set(edges, nodes, MeshKernel::Projections::cartesian);
        }

        return 0;
    }

    MESHKERNEL_API int ggeo_get_mesh(int meshKernelId, MeshGeometryDimensions& meshGeometryDimensions, MeshGeometry& meshGeometry)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        meshInstances[meshKernelId]->SetFlatCopies(MeshKernel::Mesh::AdministrationOptions::AdministrateMeshEdges);

        bool successful = SetMeshGeometry(meshKernelId, meshGeometryDimensions, meshGeometry);
        if (!successful)
        {
            return -1;
        }

        return 0;
    }

    MESHKERNEL_API int ggeo_find_faces(int meshKernelId, MeshGeometryDimensions& meshGeometryDimensions, MeshGeometry& meshGeometry)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        meshInstances[meshKernelId]->SetFlatCopies(MeshKernel::Mesh::AdministrationOptions::AdministrateMeshEdgesAndFaces);

        bool successful = SetMeshGeometry(meshKernelId,meshGeometryDimensions, meshGeometry);
        if(!successful)
        {
            return -1;
        }

        return 0;
    }

    MESHKERNEL_API int ggeo_orthogonalize(int meshKernelId, int isTriangulationRequired, int isAccountingForLandBoundariesRequired, int projectToLandBoundaryOption,
        OrthogonalizationParametersNative& orthogonalizationParametersNative, GeometryListNative& geometryListNativePolygon, GeometryListNative& geometryListNativeLandBoundaries)
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
        std::vector<MeshKernel::Point> nodes(geometryListNativePolygon.numberOfCoordinates);
        for (int i = 0; i < geometryListNativePolygon.numberOfCoordinates; i++)
        {
            nodes[i].x = geometryListNativePolygon.xCoordinates[i];
            nodes[i].y = geometryListNativePolygon.yCoordinates[i];
        }

        auto polygon = std::make_shared<MeshKernel::Polygons>();
        polygon->Set(nodes, meshInstances[meshKernelId]->m_projection);

        // build land boundary
        std::vector<MeshKernel::Point> landBoundaries(geometryListNativeLandBoundaries.numberOfCoordinates);
        for (int i = 0; i < geometryListNativeLandBoundaries.numberOfCoordinates; i++)
        {
            landBoundaries[i].x = geometryListNativeLandBoundaries.xCoordinates[i];
            landBoundaries[i].y = geometryListNativeLandBoundaries.yCoordinates[i];
        }

        auto orthogonalizer = std::make_shared<MeshKernel::Orthogonalizer>(meshInstances[meshKernelId]);
        auto smoother = std::make_shared<MeshKernel::Smoother>(meshInstances[meshKernelId]);
        auto landBoundary = std::make_shared < MeshKernel::LandBoundaries>();
        landBoundary->Set(landBoundaries, meshInstances[meshKernelId], polygon);

        MeshKernel::OrthogonalizationAndSmoothing ortogonalization;
        ortogonalization.Set(meshInstances[meshKernelId],
                             smoother,
                             orthogonalizer,
                             polygon,
                             landBoundary,
                             isTriangulationRequired,
                             isAccountingForLandBoundariesRequired,
                             projectToLandBoundaryOption,
                             orthogonalizationParametersNative);
        ortogonalization.Compute();
        return 0;
    }

    MESHKERNEL_API int ggeo_orthogonalize_initialize(int meshKernelId,
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
        std::vector<MeshKernel::Point> nodes(geometryListNativePolygon.numberOfCoordinates);
        for (int i = 0; i < geometryListNativePolygon.numberOfCoordinates; i++)
        {
            nodes[i].x = geometryListNativePolygon.xCoordinates[i];
            nodes[i].y = geometryListNativePolygon.yCoordinates[i];
        }

        // build land boundary
        std::vector<MeshKernel::Point> landBoundaries(geometryListNativeLandBoundaries.numberOfCoordinates);
        for (int i = 0; i < geometryListNativeLandBoundaries.numberOfCoordinates; i++)
        {
            landBoundaries[i].x = geometryListNativeLandBoundaries.xCoordinates[i];
            landBoundaries[i].y = geometryListNativeLandBoundaries.yCoordinates[i];
        }


        auto orthogonalizer = std::make_shared<MeshKernel::Orthogonalizer>(meshInstances[meshKernelId]);
        auto smoother = std::make_shared<MeshKernel::Smoother>(meshInstances[meshKernelId]);
        auto polygon = std::make_shared<MeshKernel::Polygons>();
        polygon->Set(nodes, meshInstances[meshKernelId]->m_projection);
        auto landBoundary = std::make_shared<MeshKernel::LandBoundaries>();
        landBoundary->Set(landBoundaries, meshInstances[meshKernelId], polygon);

        MeshKernel::OrthogonalizationAndSmoothing orthogonalizationInstance;
        orthogonalizationInstance.Set(meshInstances[meshKernelId],
                                      smoother,
                                      orthogonalizer,
                                      polygon,
                                      landBoundary,
                                      isTriangulationRequired,
                                      isAccountingForLandBoundariesRequired,
                                      projectToLandBoundaryOption,
                                      orthogonalizationParametersNative);

        orthogonalizationInstances.insert({ meshKernelId, orthogonalizationInstance });

        return 0;
    }

    MESHKERNEL_API int ggeo_orthogonalize_prepare_outer_iteration(int meshKernelId)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        if (meshInstances[meshKernelId]->GetNumNodes() <= 0)
        {
            return 0;
        }

        bool status = orthogonalizationInstances[meshKernelId].PrapareOuterIteration();
        return status == true ? 0 : 1;
    }

    MESHKERNEL_API int ggeo_orthogonalize_inner_iteration(int meshKernelId)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        if (meshInstances[meshKernelId]->GetNumNodes() <= 0)
        {
            return 0;
        }
        const bool status = orthogonalizationInstances[meshKernelId].InnerIteration();
        return status == true ? 0 : 1;
    }

    MESHKERNEL_API int ggeo_orthogonalize_finalize_outer_iteration(int meshKernelId)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        if (meshInstances[meshKernelId]->GetNumNodes() <= 0)
        {
            return 0;
        }
        const bool status = orthogonalizationInstances[meshKernelId].FinalizeOuterIteration();
        return status == true ? 0 : 1;
    }

    MESHKERNEL_API int ggeo_orthogonalize_delete(int meshKernelId)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        if (meshInstances[meshKernelId]->GetNumNodes() <= 0)
        {
            return 0;
        }

        const int returnValue = orthogonalizationInstances.erase(meshKernelId);
        return returnValue == 1 ? 0 : 1;
    }

    MESHKERNEL_API int ggeo_get_orthogonality(int meshKernelId, GeometryListNative& geometryList)
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

    MESHKERNEL_API int ggeo_get_smoothness(int meshKernelId, GeometryListNative& geometryList)
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

    MESHKERNEL_API int ggeo_get_splines(GeometryListNative& geometryListIn, GeometryListNative& geometry_list_out, int number_of_points_between_vertices)
    {

        if (geometryListIn.numberOfCoordinates == 0)
        {
            return -1;
        }

        std::vector<MeshKernel::Point> splines(geometryListIn.numberOfCoordinates);
        for (int i = 0; i < geometryListIn.numberOfCoordinates; i++)
        {
            splines[i].x = geometryListIn.xCoordinates[i];
            splines[i].y = geometryListIn.yCoordinates[i];
        }

        std::vector<std::vector<int>> indexes(geometryListIn.numberOfCoordinates, std::vector<int>(2));
        int numSplines = FindIndexes(splines, 0, splines.size(), MeshKernel::doubleMissingValue, indexes);
        std::vector<MeshKernel::Point> coordinatesDerivatives(geometryListIn.numberOfCoordinates);

        int index = 0;
        for (int s = 0; s < numSplines; s++)
        {            
            std::vector<MeshKernel::Point> coordinates(splines.begin() + indexes[s][0], splines.begin() + indexes[s][1] + 1);
            int numNodes = indexes[s][1] - indexes[s][0] + 1;
            MeshKernel::Splines::SecondOrderDerivative(coordinates, numNodes, coordinatesDerivatives);

            for (int n = 0; n < numNodes - 1; n++)
            {
                for (int p = 0; p <= number_of_points_between_vertices; p++)
                {

                    double pointAdimensionalCoordinate = n + double(p)/ double(number_of_points_between_vertices);
                    MeshKernel::Point pointCoordinate;
                    InterpolateSplinePoint(coordinates, coordinatesDerivatives, pointAdimensionalCoordinate, pointCoordinate);
                    geometry_list_out.xCoordinates[index] = pointCoordinate.x;
                    geometry_list_out.yCoordinates[index] = pointCoordinate.y;
                    geometry_list_out.zCoordinates[index] = MeshKernel::doubleMissingValue;
                    index++;
                }
            }

            geometry_list_out.xCoordinates[index] = MeshKernel::doubleMissingValue;
            geometry_list_out.yCoordinates[index] = MeshKernel::doubleMissingValue;
            geometry_list_out.zCoordinates[index] = MeshKernel::doubleMissingValue;
            index++;
        }

        geometry_list_out.numberOfCoordinates = index - 1;
        return 0;
    }

    MESHKERNEL_API int ggeo_make_mesh(int meshKernelId, MakeGridParametersNative& makeGridParameters, GeometryListNative& geometryListNative)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        MeshKernel::Polygons polygon;

        std::vector<MeshKernel::Point> result;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListNative, result);
        if (!successful)
        {
            return -1;
        }

        successful = polygon.Set(result, meshInstances[meshKernelId]->m_projection);
        if(!successful)
        {
            return -1;
        }

        MeshKernel::Mesh mesh;
        successful = mesh.MakeMesh(makeGridParameters, polygon);
        if (!successful)
        {
            return -1;
        }

        *meshInstances[meshKernelId] += mesh;

        return 0;
    }

    MESHKERNEL_API int ggeo_make_mesh_from_polygon(int meshKernelId, GeometryListNative& disposableGeometryListIn)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        MeshKernel::Polygons polygon;

        std::vector<MeshKernel::Point> result;
        bool successful = ConvertGeometryListNativeToPointVector(disposableGeometryListIn, result);
        if (!successful)
        {
            return -1;
        }

        successful = polygon.Set(result, meshInstances[meshKernelId]->m_projection);
        if (!successful)
        {
            return -1;
        }

        std::vector<std::vector<MeshKernel::Point>> generatedPoints;
        successful = polygon.CreatePointsInPolygons(generatedPoints);
        if (!successful)
        {
            return -1;
        }

        MeshKernel::Mesh mesh(generatedPoints[0], polygon, meshInstances[meshKernelId]->m_projection);
        *meshInstances[meshKernelId] += mesh;
        return 0;
    }

    MESHKERNEL_API int ggeo_make_mesh_from_samples(int meshKernelId, GeometryListNative& geometryListNative)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        std::vector<MeshKernel::Point> samplePoints;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListNative, samplePoints);
        if (!successful)
        {
            return -1;
        }
        MeshKernel::Polygons polygon;
        MeshKernel::Mesh mesh(samplePoints, polygon, meshInstances[meshKernelId]->m_projection);
        *meshInstances[meshKernelId] += mesh;
        return 0;
    }

    MESHKERNEL_API int ggeo_copy_mesh_boundaries_to_polygon(int meshKernelId, GeometryListNative& geometryListNative)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        MeshKernel::Polygons polygon;
        std::vector<MeshKernel::Point> meshBoundaryPolygon;

        int counterClockWise = 0;
        int numNodesBoundaryPolygons;
        polygon.MeshBoundaryToPolygon(*meshInstances[meshKernelId], counterClockWise, meshBoundaryPolygon, numNodesBoundaryPolygons);

        bool successful = ConvertPointVectorToGeometryListNative(meshBoundaryPolygon, geometryListNative);
        if (!successful)
        {
            return -1;
        }

        return 0;
    }

    MESHKERNEL_API int ggeo_copy_mesh_boundaries_to_polygon_count_vertices(int meshKernelId, int& numberOfPolygonVertices)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        MeshKernel::Polygons polygon;
        std::vector<MeshKernel::Point> meshBoundaryPolygon;

        int counterClockWise =0;
        polygon.MeshBoundaryToPolygon(*meshInstances[meshKernelId], counterClockWise, meshBoundaryPolygon, numberOfPolygonVertices);

        return 0;
    }

    MESHKERNEL_API int ggeo_refine_polygon(int meshKernelId, GeometryListNative& geometryListIn, int& firstIndex, int& secondIndex, double& distance, GeometryListNative& geometryListOut)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        std::vector<MeshKernel::Point> polygonPoints;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, polygonPoints);
        if (!successful)
        {
            return -1;
        }

        MeshKernel::Polygons polygon;
        polygon.Set(polygonPoints, meshInstances[meshKernelId]->m_projection);

        std::vector<MeshKernel::Point> refinedPolygon;
        successful = polygon.RefinePart(firstIndex, secondIndex, distance, refinedPolygon);
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

    MESHKERNEL_API int ggeo_refine_polygon_count(int meshKernelId, GeometryListNative& geometryListIn, int& firstIndex, int& secondIndex, double& distance, int& numberOfPolygonVertices)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        std::vector<MeshKernel::Point> polygonPoints;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, polygonPoints);
        if (!successful)
        {
            return -1;
        }

        MeshKernel::Polygons polygon;
        polygon.Set(polygonPoints, meshInstances[meshKernelId]->m_projection);

        std::vector<MeshKernel::Point> refinedPolygon;
        successful = polygon.RefinePart(firstIndex, secondIndex, distance, refinedPolygon);
        if (!successful)
        {
            return -1;
        }

        numberOfPolygonVertices = refinedPolygon.size();

        return 0;
    }


    MESHKERNEL_API int ggeo_merge_nodes(int meshKernelId, GeometryListNative& geometryListIn) 
    {

        if (meshKernelId >= meshInstances.size()) 
        {
            return 0;
        }

        std::vector<MeshKernel::Point> polygonPoints;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, polygonPoints);
        if (!successful)
        {
            return -1;
        }

        MeshKernel::Polygons polygon;
        polygon.Set(polygonPoints, meshInstances[meshKernelId]->m_projection);

        successful = meshInstances[meshKernelId]->MergeNodesInPolygon(polygon);
        if (!successful)
        {
            return -1;
        }
        return 0;
    }

    MESHKERNEL_API int ggeo_merge_two_nodes(int meshKernelId, int startNode, int endNode)
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

    MESHKERNEL_API int ggeo_nodes_in_polygons(int meshKernelId, GeometryListNative& geometryListIn, int inside, int numberOfMeshVertices, int** selectedVertices)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return 0;
        }

        std::vector<MeshKernel::Point> polygonPoints;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, polygonPoints);
        if (!successful)
        {
            return -1;
        }

        MeshKernel::Polygons polygon;
        polygon.Set(polygonPoints, meshInstances[meshKernelId]->m_projection);

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

    MESHKERNEL_API int ggeo_count_nodes_in_polygons(int meshKernelId, GeometryListNative& geometryListIn, int inside, int& numberOfMeshVertices)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return 0;
        }

        std::vector<MeshKernel::Point> polygonPoints;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, polygonPoints);
        if (!successful)
        {
            return -1;
        }

        MeshKernel::Polygons polygon;
        polygon.Set(polygonPoints, meshInstances[meshKernelId]->m_projection);

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

    MESHKERNEL_API int ggeo_insert_edge(int meshKernelId, int startNode, int endNode, int& new_edge_index)
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

    MESHKERNEL_API int ggeo_insert_node(int meshKernelId, double xCoordinate, double yCoordinate, double zCoordinate, int& vertexIndex)
    {
        if (meshKernelId >= meshInstances.size())
        {
            //create a valid instance, by default cartesian
            *meshInstances[meshKernelId] = MeshKernel::Mesh();
            meshInstances[meshKernelId]->m_projection = MeshKernel::Projections::cartesian;
        }

        MeshKernel::Point newNode{ xCoordinate, yCoordinate };
        bool successful = meshInstances[meshKernelId]->InsertNode(newNode, vertexIndex);
        if (!successful)
        {
            return -1;
        }

        return 0;
    }


    MESHKERNEL_API int ggeo_delete_node(int meshKernelId, int nodeIndex)
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

    MESHKERNEL_API int ggeo_move_node(int meshKernelId, GeometryListNative& geometryListIn, int nodeIndex)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return 0;
        }

        std::vector<MeshKernel::Point> newPoint;
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

    MESHKERNEL_API int ggeo_delete_edge(int meshKernelId, GeometryListNative& geometryListIn, double searchRadius)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return 0;
        }

        std::vector<MeshKernel::Point> newPoint;
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

    MESHKERNEL_API int ggeo_find_edge(int meshKernelId, GeometryListNative& geometryListIn, double searchRadius, int& edgeIndex)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return 0;
        }

        std::vector<MeshKernel::Point> newPoint;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, newPoint);
        if (!successful || newPoint.size() != 1)
        {
            return -1;
        }

        edgeIndex = meshInstances[meshKernelId]->FindEdgeCloseToAPoint(newPoint[0], searchRadius);

        return 0;
    }

    MESHKERNEL_API int ggeo_offsetted_polygon(int meshKernelId, GeometryListNative& geometryListIn, bool innerAndOuter, double distance, GeometryListNative& geometryListOut)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return 0;
        }

        std::vector<MeshKernel::Point> polygonPoints;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, polygonPoints);
        if (!successful)
        {
            return -1;
        }

        MeshKernel::Polygons polygon;
        successful = polygon.Set(polygonPoints, meshInstances[meshKernelId]->m_projection);
        if (!successful)
        {
            return -1;
        }

        MeshKernel::Polygons newPolygon;
        successful = polygon.OffsetCopy(0, distance, innerAndOuter, newPolygon);
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

    MESHKERNEL_API int ggeo_offsetted_polygon_count(int meshKernelId, GeometryListNative& geometryListIn, bool innerAndOuter, double distance, int& numberOfPolygonVertices)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return 0;
        }

        std::vector<MeshKernel::Point> polygonPoints;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, polygonPoints);
        if (!successful)
        {
            return -1;
        }

        MeshKernel::Polygons polygon;
        successful = polygon.Set(polygonPoints, meshInstances[meshKernelId]->m_projection);
        if (!successful)
        {
            return -1;
        }

        MeshKernel::Polygons newPolygon;
        successful = polygon.OffsetCopy(0, distance, innerAndOuter, newPolygon);
        if (!successful)
        {
            return -1;
        }

        numberOfPolygonVertices = newPolygon.GetNumNodes();
    
        return 0;
    }

    MESHKERNEL_API int ggeo_refine_mesh_based_on_samples(int meshKernelId, GeometryListNative& geometryListIn, InterpolationParametersNative& interpolationParametersNative, SampleRefineParametersNative& sampleRefineParametersNative)
    {
        if (meshKernelId >= meshInstances.size() || meshInstances[meshKernelId]->GetNumNodes() <= 0)
        {
            return 0;
        }

        std::vector<MeshKernel::Sample> samples;
        bool successful = ConvertGeometryListNativeToSampleVector(geometryListIn, samples);
        if (!successful)
        {
            return -1;
        }

        if (samples.empty())
        {
            return 0;
        }

        MeshKernel::MeshRefinement meshRefinement(meshInstances[meshKernelId]);

        // polygon could be passed as api parameter
        MeshKernel::Polygons polygon;

        successful = meshRefinement.Refine(samples, polygon, sampleRefineParametersNative, interpolationParametersNative);
        if (!successful)
        {
            return -1;
        }

        return 0;
    }

    MESHKERNEL_API int ggeo_refine_mesh_based_on_polygon(int meshKernelId, GeometryListNative& geometryListNative, InterpolationParametersNative& interpolationParametersNative)
    {
        
        if (meshKernelId >= meshInstances.size() || meshInstances[meshKernelId]->GetNumNodes()<=0)
        {
            return 0;
        }


        std::vector<MeshKernel::Point> points;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListNative, points);
        if (!successful)
        {
            return -1;
        }

        if (points.empty())
        {
            return 0;
        }


        MeshKernel::Polygons polygon;
        successful = polygon.Set(points, meshInstances[meshKernelId]->m_projection);
        if (!successful)
        {
            return -1;
        }

        MeshKernel::MeshRefinement meshRefinement(meshInstances[meshKernelId]);

        // polygon could be passed as api parameter
        std::vector<MeshKernel::Sample> samples; 
        SampleRefineParametersNative sampleRefineParametersNative;
        successful = meshRefinement.Refine(samples, polygon, sampleRefineParametersNative, interpolationParametersNative);
        if (!successful)
        {
            return -1;
        }

        return 0;
    }

    MESHKERNEL_API int ggeo_get_node_index(int meshKernelId, GeometryListNative& geometryListIn, double searchRadius, int& vertexIndex)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return 0;
        }

        std::vector<MeshKernel::Point> polygonPoints;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, polygonPoints);
        if (!successful || polygonPoints.empty())
        {
            return -1;
        }

        successful = meshInstances[meshKernelId]->GetNodeIndex(polygonPoints[0], searchRadius, vertexIndex);
        if (!successful)
        {
            return -1;
        }

        return 0;
    }

    MESHKERNEL_API int ggeo_curvilinear_mesh_from_splines_ortho(int meshKernelId,
        GeometryListNative& geometryListIn,
        CurvilinearParametersNative& curvilinearParameters,
        SplinesToCurvilinearParametersNative& splineToCurvilinearParameters)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        // use the default constructor, no instance present
        auto spline = std::make_shared<MeshKernel::Splines>(meshInstances[meshKernelId]->m_projection);
        bool successful = SetSplines(geometryListIn, *spline);
        if(!successful)
        {
            return -1;
        }

        MeshKernel::CurvilinearGridFromSplines curvilinearGridFromSplines(spline);
        curvilinearGridFromSplines.SetParameters(curvilinearParameters, splineToCurvilinearParameters);
        MeshKernel::CurvilinearGrid curvilinearGrid;
        curvilinearGridFromSplines.Compute(curvilinearGrid);
        *meshInstances[meshKernelId] += MeshKernel::Mesh(curvilinearGrid, meshInstances[meshKernelId]->m_projection);

        return 0;
    }


    MESHKERNEL_API int ggeo_curvilinear_mesh_from_splines_ortho_initialize(int meshKernelId, GeometryListNative& geometryListNative, CurvilinearParametersNative& curvilinearParametersNative, SplinesToCurvilinearParametersNative& splinesToCurvilinearParametersNative)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return 0;
        }

        auto spline = std::make_shared<MeshKernel::Splines>(meshInstances[meshKernelId]->m_projection);
        bool successful = SetSplines(geometryListNative, *spline);
        if (!successful)
        {
            return -1;
        }
        MeshKernel::CurvilinearGridFromSplines curvilinearGridFromSplines(spline);

        successful = curvilinearGridFromSplines.SetParameters(curvilinearParametersNative, splinesToCurvilinearParametersNative);
        if (!successful)
        {
            return -1;
        }

        splineInstances.insert({ meshKernelId, curvilinearGridFromSplines });
 
        successful = splineInstances[meshKernelId].Initialize();
        if (!successful)
        {
            return -1;
        }
        return 0;
    }

    MESHKERNEL_API int ggeo_curvilinear_mesh_from_splines_iteration(int meshKernelId, int layer)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return 0;
        }

        const bool successful = splineInstances[meshKernelId].Iterate(layer);
        return successful ? 0 : 1;
    }


    MESHKERNEL_API int ggeo_curvilinear_mesh_from_splines_ortho_refresh_mesh(int meshKernelId)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return 0;
        }

        MeshKernel::CurvilinearGrid curvilinearGrid;
        const bool successful = splineInstances[meshKernelId].ComputeCurvilinearGrid(curvilinearGrid);

        *meshInstances[meshKernelId] += MeshKernel::Mesh(curvilinearGrid, meshInstances[meshKernelId]->m_projection);

        return successful ? 0 : 1;
    }

    MESHKERNEL_API int ggeo_curvilinear_mesh_from_splines_ortho_delete(int meshKernelId)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return 0;
        }
        const int successful = splineInstances.erase(meshKernelId);
        return successful;
    }

    MESHKERNEL_API int ggeo_points_in_polygon(int meshKernelId, GeometryListNative& polygonNative, GeometryListNative& pointsNative, GeometryListNative& selectedPointsNative)
    {
    
        std::vector<MeshKernel::Point> polygonNodes;
        bool successful = ConvertGeometryListNativeToPointVector(polygonNative, polygonNodes);
        if (!successful)
        {
            return -1;
        }

        std::vector<MeshKernel::Point> points;
        successful = ConvertGeometryListNativeToPointVector(pointsNative, points);
        if (!successful || points.empty() || points.size()!= selectedPointsNative.numberOfCoordinates)
        {
            return -1;
        }

        
        MeshKernel::Polygons polygon;
        polygon.Set(polygonNodes, meshInstances[meshKernelId]->m_projection);

        for (int i = 0; i < points.size(); i++)
        {
            selectedPointsNative.zCoordinates[i]= polygon.IsPointInPolygon(points[i], 0)? 1.0 : 0.0;
        }
        
        return successful ? 0 : 1;
    }

    MESHKERNEL_API int ggeo_flip_edges(int meshKernelId, 
                                     int isTriangulationRequired, 
                                     int isAccountingForLandBoundariesRequired, 
                                     int projectToLandBoundaryRequired)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        std::vector<MeshKernel::Point> landBoundary;
        auto landBoundaries = std::make_shared<MeshKernel::LandBoundaries>();

        bool triangulateFaces = isTriangulationRequired == 0 ? false : true;
        bool projectToLandBoundary = projectToLandBoundaryRequired == 0 ? false : true;
        MeshKernel::FlipEdges flipEdges(meshInstances[meshKernelId], landBoundaries, triangulateFaces, projectToLandBoundary);

        // compute Flip edges
        auto successful = flipEdges.Compute();

        return successful;
    }

    MESHKERNEL_API int ggeo_curvilinear_mesh_from_splines(int meshKernelId, GeometryListNative& geometryListNativeIn, CurvilinearParametersNative& curvilinearParametersNative)
    {
        if (meshKernelId >= meshInstances.size())
        {
            return -1;
        }

        // Use the default constructor, no instance present
        auto spline = std::make_shared<MeshKernel::Splines>(meshInstances[meshKernelId]->m_projection);
        bool successful = SetSplines(geometryListNativeIn, *spline);
        if (!successful)
        {
            return -1;
        }

        // Create algorithm and set the splines
        MeshKernel::CurvilinearGridFromSplinesTransfinite curvilinearGridFromSplinesTransfinite(spline);
        curvilinearGridFromSplinesTransfinite.Set(curvilinearParametersNative);
        
        // Compute the curvilinear grid
        MeshKernel::CurvilinearGrid curvilinearGrid;
        bool success = curvilinearGridFromSplinesTransfinite.Compute(curvilinearGrid);

        // Transform and set mesh pointer
        meshInstances[meshKernelId] = std::make_unique<MeshKernel::Mesh>(curvilinearGrid, meshInstances[meshKernelId]->m_projection);

        return 0;
    }

}