#include <map>

#include "Gridgeom.hpp"
#include "Mesh.hpp"
#include "Orthogonalization.hpp"
#include "CurvilinearGrid.hpp"
#include "Splines.hpp"
#include "Entities.hpp"
#include "MeshRefinement.hpp"

static std::vector<GridGeom::Mesh> meshInstances;
static std::map<int, GridGeom::Orthogonalization> orthogonalizationInstances;
static std::map<int, GridGeom::Polygons> polygonInstances;

namespace GridGeomApi
{
    static bool ConvertGeometryListNativeToPointVector(GeometryListNative& geometryListIn, std::vector<GridGeom::Point>& result)
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

    static bool ConvertGeometryListNativeToSampleVector(GeometryListNative& geometryListIn, std::vector<GridGeom::Sample>& result)
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

    static bool ConvertPointVectorToGeometryListNative(std::vector<GridGeom::Point> pointVector, GeometryListNative& result)
    {
        // invalid memory allocation
        if (pointVector.size()!= result.numberOfCoordinates)
        {
            return -1;
        }

        for (int i = 0; i < pointVector.size(); i++)
        {
            result.xCoordinates[i] = pointVector[i].x;
            result.yCoordinates[i] = pointVector[i].y;
        }
        return true;
    }

    static int SetSplines(GeometryListNative& geometryListIn, GridGeom::Splines& spline)
    {
        if (geometryListIn.numberOfCoordinates == 0)
        {
            return -1;
        }

        std::vector<GridGeom::Point> splineCornerPoints;
        bool success = ConvertGeometryListNativeToPointVector(geometryListIn, splineCornerPoints);
        if (!success)
        {
            return -1;
        }

        std::vector<std::vector<int>> indexes(splineCornerPoints.size(), std::vector<int>(2,0));
        int pos = FindIndexes(splineCornerPoints, 0, splineCornerPoints.size(), GridGeom::doubleMissingValue, indexes);
        indexes.resize(pos);

        for (int i = 0; i < indexes.size(); i++)
        {
            int size = indexes[i][1] - indexes[i][0] + 1;
            if (size > 0)
            {
                spline.AddSpline(splineCornerPoints, indexes[i][0], size);
            }
        }

        return 0;
    }

    static bool IsValidInstance(int gridStateId)
    {
        if (gridStateId >= meshInstances.size() || !meshInstances[gridStateId].IsSet())
        {
            return false;
        }
        return true;
    }

    GRIDGEOM_API int ggeo_new_grid(int& gridStateId)
    {
        int instanceSize = meshInstances.size();
        meshInstances.resize(instanceSize + 1);
        gridStateId = instanceSize;
        return 0;
    };

    GRIDGEOM_API int ggeo_deallocate_state(int& gridStateId)
    {
        meshInstances[gridStateId].DeleteFlatCopies();
        return 0;
    }

    GRIDGEOM_API int ggeo_delete_mesh(int& gridStateId, GeometryListNative& geometryListNativePolygon, int& deletionOption)
    {
        meshInstances[gridStateId].DeleteFlatCopies();
        return 0;
    }

    GRIDGEOM_API int ggeo_set_state(int& gridStateId, MeshGeometryDimensions& meshGeometryDimensions, MeshGeometry& meshGeometry, bool IsGeographic)
    {

        std::vector<GridGeom::Edge> edges(meshGeometryDimensions.numedge);
        int ei = 0;
        for (int e = 0; e < edges.size(); e++)
        {
            edges[e].first = meshGeometry.edge_nodes[ei];
            ei++;
            edges[e].second = meshGeometry.edge_nodes[ei];
            ei++;
        }

        std::vector<GridGeom::Point> nodes(meshGeometryDimensions.numnode);
        for (int n = 0; n < nodes.size(); n++)
        {
            nodes[n].x = meshGeometry.nodex[n];
            nodes[n].y = meshGeometry.nodey[n];
        }

        // TODO: re-enable switch
        //if (IsGeographic)
        //{
        meshInstances[gridStateId].Set(edges, nodes, GridGeom::Projections::cartesian);
        //}
        //else
        //{
        //    auto instance = std::make_unique<GridGeom::Mesh<GridGeom::OperationTypes::sphericalOperations>>();
        //    instance->Set(edges, nodes);
        //}

        return 0;
    }

    GRIDGEOM_API int ggeo_get_mesh(int& gridStateId, MeshGeometryDimensions& meshGeometryDimensions, MeshGeometry& meshGeometry)
    {

        meshInstances[gridStateId].SetFlatCopies();

        meshGeometry.nodex = &meshInstances[gridStateId].m_nodex[0];
        meshGeometry.nodey = &meshInstances[gridStateId].m_nodey[0];
        meshGeometry.nodez = &meshInstances[gridStateId].m_nodez[0];
        meshGeometry.edge_nodes = &meshInstances[gridStateId].m_edgeNodes[0];

        if (meshInstances[gridStateId].GetNumNodes() == 1)
        {
            meshGeometryDimensions.numnode = 0;
            meshGeometryDimensions.numedge = 0;
            meshGeometryDimensions.numface = 0;
            meshGeometryDimensions.maxnumfacenodes = 4;
        }
        else
        {
            meshGeometryDimensions.numnode = meshInstances[gridStateId].GetNumNodes();
            meshGeometryDimensions.numedge = meshInstances[gridStateId].GetNumEdges();
            meshGeometryDimensions.numface = meshInstances[gridStateId].GetNumFaces();
            meshGeometryDimensions.maxnumfacenodes = 4;
        }

        return 0;
    }

    GRIDGEOM_API int ggeo_orthogonalize(int& gridStateId, int& isTriangulationRequired, int& isAccountingForLandBoundariesRequired, int& projectToLandBoundaryOption,
        OrthogonalizationParametersNative& orthogonalizationParametersNative, GeometryListNative& geometryListNativePolygon, GeometryListNative& geometryListNativeLandBoundaries)
    {
        GridGeom::Orthogonalization ortogonalization;

        // build enclosing polygon
        std::vector<GridGeom::Point> polygon(geometryListNativePolygon.numberOfCoordinates);
        for (int i = 0; i < geometryListNativePolygon.numberOfCoordinates; i++)
        {
            polygon[i].x = geometryListNativePolygon.xCoordinates[i];
            polygon[i].y = geometryListNativePolygon.yCoordinates[i];
        }

        // build land boundary
        std::vector<GridGeom::Point> landBoundaries(geometryListNativeLandBoundaries.numberOfCoordinates);
        for (int i = 0; i < geometryListNativeLandBoundaries.numberOfCoordinates; i++)
        {
            landBoundaries[i].x = geometryListNativeLandBoundaries.xCoordinates[i];
            landBoundaries[i].y = geometryListNativeLandBoundaries.yCoordinates[i];
        }

        ortogonalization.Set(meshInstances[gridStateId],
            isTriangulationRequired,
            isAccountingForLandBoundariesRequired,
            projectToLandBoundaryOption,
            orthogonalizationParametersNative,
            polygon,
            landBoundaries);
        ortogonalization.Iterate(meshInstances[gridStateId]);
        return 0;
    }

    GRIDGEOM_API int ggeo_orthogonalize_initialize(int& gridStateId,
        int& isTriangulationRequired,
        int& isAccountingForLandBoundariesRequired,
        int& projectToLandBoundaryOption,
        OrthogonalizationParametersNative& orthogonalizationParametersNative,
        GeometryListNative& geometryListNativePolygon,
        GeometryListNative& geometryListNativeLandBoundaries)
    {
        // build enclosing polygon
        std::vector<GridGeom::Point> polygon(geometryListNativePolygon.numberOfCoordinates);
        for (int i = 0; i < geometryListNativePolygon.numberOfCoordinates; i++)
        {
            polygon[i].x = geometryListNativePolygon.xCoordinates[i];
            polygon[i].y = geometryListNativePolygon.yCoordinates[i];
        }

        // build land boundary
        std::vector<GridGeom::Point> landBoundaries(geometryListNativeLandBoundaries.numberOfCoordinates);
        for (int i = 0; i < geometryListNativeLandBoundaries.numberOfCoordinates; i++)
        {
            landBoundaries[i].x = geometryListNativeLandBoundaries.xCoordinates[i];
            landBoundaries[i].y = geometryListNativeLandBoundaries.yCoordinates[i];
        }

        orthogonalizationInstances[gridStateId].Set(meshInstances[gridStateId],
            isTriangulationRequired,
            isAccountingForLandBoundariesRequired,
            projectToLandBoundaryOption,
            orthogonalizationParametersNative,
            polygon,
            landBoundaries);
        return 0;
    }

    GRIDGEOM_API int ggeo_orthogonalize_prepare_outer_iteration(int& gridStateId)
    {
        bool status = orthogonalizationInstances[gridStateId].PrapareOuterIteration(meshInstances[gridStateId]);
        return status == true ? 0 : 1;
    }

    GRIDGEOM_API int ggeo_orthogonalize_inner_iteration(int& gridStateId)
    {
        const bool status = orthogonalizationInstances[gridStateId].InnerIteration(meshInstances[gridStateId]);
        return status == true ? 0 : 1;
    }

    GRIDGEOM_API int ggeo_orthogonalize_finalize_outer_iteration(int& gridStateId)
    {
        const bool status = orthogonalizationInstances[gridStateId].FinalizeOuterIteration(meshInstances[gridStateId]);
        return status == true ? 0 : 1;
    }

    GRIDGEOM_API int ggeo_orthogonalize_delete(int& gridStateId)
    {
        const int returnValue = orthogonalizationInstances.erase(gridStateId);
        return returnValue == 1 ? 0 : 1;
    }

    GRIDGEOM_API int ggeo_get_orthogonality(int& gridStateId, GeometryListNative& geometryList)
    {
        const bool status = orthogonalizationInstances[gridStateId].GetOrthogonality(meshInstances[gridStateId], geometryList.zCoordinates);
        return status == true ? 0 : 1;
    }

    GRIDGEOM_API int ggeo_get_smoothness(int& gridStateId, GeometryListNative& geometryList)
    {
        const bool status = orthogonalizationInstances[gridStateId].GetSmoothness(meshInstances[gridStateId], geometryList.zCoordinates);
        return status == true ? 0 : 1;
    }

    GRIDGEOM_API int ggeo_get_splines(GeometryListNative& geometryListIn, GeometryListNative& geometry_list_out, int& number_of_points_between_vertices)
    {

        if (geometryListIn.numberOfCoordinates == 0)
        {
            return -1;
        }

        std::vector<GridGeom::Point> splines(geometryListIn.numberOfCoordinates);
        for (int i = 0; i < geometryListIn.numberOfCoordinates; i++)
        {
            splines[i].x = geometryListIn.xCoordinates[i];
            splines[i].y = geometryListIn.yCoordinates[i];
        }

        std::vector<std::vector<int>> indexes(geometryListIn.numberOfCoordinates, std::vector<int>(2));
        int numSplines = FindIndexes(splines, 0, splines.size(), GridGeom::doubleMissingValue, indexes);
        std::vector<GridGeom::Point> coordinatesDerivatives(geometryListIn.numberOfCoordinates);

        int index = 0;
        for (int s = 0; s < numSplines; s++)
        {            
            std::vector<GridGeom::Point> coordinates(splines.begin() + indexes[s][0], splines.begin() + indexes[s][1] + 1);
            int numNodes = indexes[s][1] - indexes[s][0] + 1;
            GridGeom::Splines::SecondOrderDerivative(coordinates, numNodes, coordinatesDerivatives);

            for (int n = 0; n < numNodes - 1; n++)
            {
                for (int p = 0; p <= number_of_points_between_vertices; p++)
                {

                    double pointAdimensionalCoordinate = n + double(p)/ double(number_of_points_between_vertices);
                    GridGeom::Point pointCoordinate;
                    GridGeom::Splines::Interpolate(coordinates, coordinatesDerivatives, pointAdimensionalCoordinate, pointCoordinate);
                    geometry_list_out.xCoordinates[index] = pointCoordinate.x;
                    geometry_list_out.yCoordinates[index] = pointCoordinate.y;
                    geometry_list_out.zCoordinates[index] = GridGeom::doubleMissingValue;
                    index++;
                }
            }

            geometry_list_out.xCoordinates[index] = GridGeom::doubleMissingValue;
            geometry_list_out.yCoordinates[index] = GridGeom::doubleMissingValue;
            geometry_list_out.zCoordinates[index] = GridGeom::doubleMissingValue;
            index++;
        }

        geometry_list_out.numberOfCoordinates = index - 1;
        return 0;
    }

    GRIDGEOM_API int ggeo_curvilinear_mesh_from_splines_ortho(int& gridStateId, 
        GeometryListNative& geometryListIn, 
        CurvilinearParametersNative& curvilinearParameters, 
        SplinesToCurvilinearParametersNative& splineToCurvilinearParameters)
    {
        // use the default constructor, no instance present
        GridGeom::Polygons polygon;
        GridGeom::Splines spline(meshInstances[gridStateId].m_projection, polygon);

        int success = SetSplines(geometryListIn, spline);
        spline.SetParameters(curvilinearParameters, splineToCurvilinearParameters);
        GridGeom::CurvilinearGrid curvilinearGrid;
        spline.OrthogonalCurvilinearGridFromSplines(curvilinearGrid);
        meshInstances[gridStateId] = GridGeom::Mesh(curvilinearGrid, meshInstances[gridStateId].m_projection);

        return 0;
    }


    GRIDGEOM_API int ggeo_make_net(int& gridStateId, MakeGridParametersNative& makeGridParameters, GeometryListNative& disposableGeometryListIn)
    {
        GridGeom::Polygons polygon;

        std::vector<GridGeom::Point> result;
        bool successful = ConvertGeometryListNativeToPointVector(disposableGeometryListIn, result);
        if (!successful)
        {
            return -1;
        }

        successful = polygon.Set(result, meshInstances[gridStateId].m_projection);
        if(!successful)
        {
            return -1;
        }

        successful = meshInstances[gridStateId].MakeMesh(makeGridParameters, polygon);
        if (!successful)
        {
            return -1;
        }
        return 0;
    }

    GRIDGEOM_API int ggeo_mesh_from_polygon(int& gridStateId, GeometryListNative& disposableGeometryListIn)
    {
        GridGeom::Polygons polygon;

        std::vector<GridGeom::Point> result;
        bool successful = ConvertGeometryListNativeToPointVector(disposableGeometryListIn, result);
        if (!successful)
        {
            return -1;
        }

        successful = polygon.Set(result, meshInstances[gridStateId].m_projection);
        if (!successful)
        {
            return -1;
        }

        std::vector<std::vector<GridGeom::Point>> generatedPoints;
        successful = polygon.CreatePointsInPolygons(generatedPoints);
        if (!successful)
        {
            return -1;
        }

        GridGeom::Mesh mesh(generatedPoints[0], polygon, meshInstances[gridStateId].m_projection);
        meshInstances[gridStateId] = mesh;
        return 0;
    }

    GRIDGEOM_API int ggeo_mesh_from_samples(int& gridStateId, GeometryListNative& disposableGeometryListIn)
    {
        std::vector<GridGeom::Point> samplePoints;
        bool successful = ConvertGeometryListNativeToPointVector(disposableGeometryListIn, samplePoints);
        if (!successful)
        {
            return -1;
        }
        GridGeom::Polygons polygon;
        GridGeom::Mesh mesh(samplePoints, polygon, meshInstances[gridStateId].m_projection);
        meshInstances[gridStateId] = mesh;
        return 0;
    }

    GRIDGEOM_API int ggeo_copy_mesh_boundaries_to_polygon_count_edges(int& gridStateId, int& numberOfPolygonVertices)
    {
        GridGeom::Polygons polygon;
        std::vector<GridGeom::Point> meshBoundaryPolygon;

        int counterClockWise =0;
        int setMeshState = 0;
        polygon.MeshBoundaryToPolygon(meshInstances[gridStateId], counterClockWise, setMeshState, meshBoundaryPolygon, numberOfPolygonVertices);

        return 0;
    }

    GRIDGEOM_API int ggeo_copy_mesh_boundaries_to_polygon(int& gridStateId, GeometryListNative& disposableGeometryListInOut)
    {
        GridGeom::Polygons polygon;
        std::vector<GridGeom::Point> meshBoundaryPolygon;

        int counterClockWise = 0;
        int setMeshState = 0;
        int numNodesBoundaryPolygons;
        polygon.MeshBoundaryToPolygon(meshInstances[gridStateId], counterClockWise, setMeshState, meshBoundaryPolygon, numNodesBoundaryPolygons);

        bool successful = ConvertPointVectorToGeometryListNative(meshBoundaryPolygon, disposableGeometryListInOut);
        if (!successful)
        {
            return -1;
        }

        return 0;
    }

    GRIDGEOM_API int ggeo_refine_polygon_count(int& gridStateId, GeometryListNative& geometryListIn, int& firstIndex, int& secondIndex, double& distance, int& numberOfPolygonVertices)
    {
        std::vector<GridGeom::Point> polygonPoints;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, polygonPoints);
        if (!successful)
        {
            return -1;
        }

        GridGeom::Polygons polygon;
        polygon.Set(polygonPoints, meshInstances[gridStateId].m_projection);

        std::vector<GridGeom::Point> refinedPolygon;
        successful = polygon.RefinePart(firstIndex, secondIndex, distance, refinedPolygon);
        if (!successful)
        {
            return -1;
        }

        numberOfPolygonVertices = refinedPolygon.size();

        return 0;
    }

    GRIDGEOM_API int ggeo_refine_polygon(int& gridStateId, GeometryListNative& geometryListIn, int& firstIndex, int& secondIndex, double& distance, GeometryListNative& geometryListOut)
    {
        std::vector<GridGeom::Point> polygonPoints;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, polygonPoints);
        if (!successful)
        {
            return -1;
        }

        GridGeom::Polygons polygon;
        polygon.Set(polygonPoints, meshInstances[gridStateId].m_projection);

        std::vector<GridGeom::Point> refinedPolygon;
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


    GRIDGEOM_API int ggeo_merge_nodes(int& gridStateId, GeometryListNative& geometryListIn) 
    {

        if (!IsValidInstance(gridStateId)) 
        {
            return 0;
        }

        if (!meshInstances[gridStateId].IsSet())
        {
            return 0;
        }

        std::vector<GridGeom::Point> polygonPoints;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, polygonPoints);
        if (!successful)
        {
            return -1;
        }

        GridGeom::Polygons polygon;
        polygon.Set(polygonPoints, meshInstances[gridStateId].m_projection);

        successful = meshInstances[gridStateId].MergeNodesInPolygon(polygon);
        if (!successful)
        {
            return -1;
        }
        return 0;
    }


    GRIDGEOM_API int ggeo_count_vertices_in_polygons(int& gridStateId, GeometryListNative& geometryListIn, int& inside, int& numberOfMeshVertices) 
    {
        if (!IsValidInstance(gridStateId))
        {
            return 0;
        }

        std::vector<GridGeom::Point> polygonPoints;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, polygonPoints);
        if (!successful)
        {
            return -1;
        }

        GridGeom::Polygons polygon;
        polygon.Set(polygonPoints, meshInstances[gridStateId].m_projection);

        successful = meshInstances[gridStateId].SelectNodesInPolygon(polygon, inside);
        if (!successful)
        {
            return -1;
        }
        
        numberOfMeshVertices = 0;
        for (auto i = 0; i < meshInstances[gridStateId].GetNumNodes() ; ++i) 
        {
            if (meshInstances[gridStateId].m_nodeMask[i])
            {
                numberOfMeshVertices++;
            }
        }

        return 0;
    }

    GRIDGEOM_API int ggeo_vertices_in_polygons(int& gridStateId, GeometryListNative& geometryListIn, int& inside, int& numberOfMeshVertices, int* selectedVertices)
    {
        if (!IsValidInstance(gridStateId))
        {
            return 0;
        }

        std::vector<GridGeom::Point> polygonPoints;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, polygonPoints);
        if (!successful)
        {
            return -1;
        }

        GridGeom::Polygons polygon;
        polygon.Set(polygonPoints, meshInstances[gridStateId].m_projection);

        successful = meshInstances[gridStateId].SelectNodesInPolygon(polygon,inside);
        if (!successful)
        {
            return -1;
        }

        numberOfMeshVertices = 0;
        for (int i = 0; i < meshInstances[gridStateId].GetNumNodes() ; ++i)
        {
            if (meshInstances[gridStateId].m_nodeMask[i]) 
            {
                selectedVertices[numberOfMeshVertices] = i;
                numberOfMeshVertices++;
            }
        }
        return 0;
    }

    GRIDGEOM_API int ggeo_insert_edge(int& gridStateId, int& start_node, int& end_node, int& new_edge_index)
    {
        if (!IsValidInstance(gridStateId))
        {
            return 0;
        }

        bool successful = meshInstances[gridStateId].ConnectNodes(start_node, end_node, new_edge_index);
        if (!successful)
        {
            return -1;
        }


        return 0;
    }

    GRIDGEOM_API int ggeo_insert_node(int& gridStateId, double& xCoordinate, double& yCoordinate, double& zCoordinate, int& vertexIndex)
    {
        if (!IsValidInstance(gridStateId))
        {
            return 0;
        }

        GridGeom::Point newNode{ xCoordinate, yCoordinate };
        bool successful = meshInstances[gridStateId].InsertNode(newNode, vertexIndex);
        if (!successful)
        {
            return -1;
        }

        return 0;
    }


    GRIDGEOM_API int ggeo_delete_node(int& gridStateId, int& nodeIndex)
    {
        if (!IsValidInstance(gridStateId))
        {
            return 0;
        }

        bool successful = meshInstances[gridStateId].DeleteNode(nodeIndex);
        if (!successful)
        {
            return -1;
        }

        return 0;
    }

    GRIDGEOM_API int ggeo_offsetted_polygon_count(int& gridStateId, GeometryListNative& geometryListIn, bool& innerAndOuter, double& distance, int& numberOfPolygonVertices)
    {
        if (!IsValidInstance(gridStateId))
        {
            return 0;
        }

        std::vector<GridGeom::Point> polygonPoints;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, polygonPoints);
        if (!successful)
        {
            return -1;
        }

        GridGeom::Polygons polygon;
        successful = polygon.Set(polygonPoints, meshInstances[gridStateId].m_projection);
        if (!successful)
        {
            return -1;
        }

        GridGeom::Polygons newPolygon;
        successful = polygon.OffsetCopy(0, distance, innerAndOuter, newPolygon);
        if (!successful)
        {
            return -1;
        }

        numberOfPolygonVertices = newPolygon.GetNumNodes();
    
        return 0;
    }

    GRIDGEOM_API int ggeo_offsetted_polygon(int& gridStateId, GeometryListNative& geometryListIn, bool& innerAndOuter, double& distance, GeometryListNative& geometryListOut)
    {
        if (!IsValidInstance(gridStateId))
        {
            return 0;
        }

        std::vector<GridGeom::Point> polygonPoints;
        bool successful = ConvertGeometryListNativeToPointVector(geometryListIn, polygonPoints);
        if (!successful)
        {
            return -1;
        }

        GridGeom::Polygons polygon;
        successful = polygon.Set(polygonPoints, meshInstances[gridStateId].m_projection);
        if (!successful)
        {
            return -1;
        }

        GridGeom::Polygons newPolygon;
        successful = polygon.OffsetCopy(0, distance, innerAndOuter, newPolygon);
        if (!successful)
        {
            return -1;
        }
    
        return 0;
    }

    GRIDGEOM_API int ggeo_refine_mesh_based_on_samples(int& gridStateId, GeometryListNative& geometryListIn, InterpolationParametersNative& interpolationParametersNative, SampleRefineParametersNative& sampleRefineParametersNative)
    {
        if (!IsValidInstance(gridStateId))
        {
            return 0;
        }

        std::vector<GridGeom::Sample> samples;
        bool successful = ConvertGeometryListNativeToSampleVector(geometryListIn, samples);
        if (!successful)
        {
            return -1;
        }

        if (samples.empty())
        {
            return 0;
        }

        GridGeom::MeshRefinement meshRefinement(meshInstances[gridStateId]);

        // polygon could be passed as api parameter
        GridGeom::Polygons polygon;

        successful = meshRefinement.RefineMeshBasedOnSamples(samples, polygon, sampleRefineParametersNative, interpolationParametersNative);
        if (!successful)
        {
            return -1;
        }

        return 0;
    }

}