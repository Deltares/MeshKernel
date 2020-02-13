#include <map>

#include "Gridgeom.hpp"
#include "Mesh.hpp"
#include "Orthogonalization.hpp"
#include "CurvilinearGrid.hpp"
#include "Splines.hpp"

static std::vector<GridGeom::Mesh> meshInstances;
static std::map<int, GridGeom::Orthogonalization> orthogonalizationInstances;
static std::map<int, GridGeom::Splines> splinesInstances;
static std::map<int, GridGeom::Polygons> polygonInstances;

namespace GridGeomApi
{
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

        if (meshInstances[gridStateId].m_nodex.size() == 1)
        {
            meshGeometryDimensions.numnode = 0;
            meshGeometryDimensions.numedge = 0;
            meshGeometryDimensions.numface = 0;
            meshGeometryDimensions.maxnumfacenodes = 4;
        }
        else
        {
            meshGeometryDimensions.numnode = meshInstances[gridStateId].m_nodex.size();
            meshGeometryDimensions.numedge = meshInstances[gridStateId].m_edgeNodes.size() / 2;
            meshGeometryDimensions.numface = meshInstances[gridStateId].m_numFaces;
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

    GRIDGEOM_API int ggeo_set_splines(int& gridStateId, GeometryListNative& geometryListIn)
    {
        if (geometryListIn.numberOfCoordinates == 0)
        {
            return -1;
        }

        // use the default constructor, no instance present
        if (splinesInstances.count(gridStateId) == 0)
        {
            GridGeom::Polygons polygon;
            GridGeom::Splines splines(meshInstances[gridStateId].m_projection, polygon);
            splinesInstances[gridStateId] = splines;
        }

        std::vector<GridGeom::Point> spline(geometryListIn.numberOfCoordinates);
        int ind = 0;
        for (int i = 0; i < geometryListIn.numberOfCoordinates; i++)
        {
            if (geometryListIn.xCoordinates[i] == GridGeom::doubleMissingValue)
            {
                splinesInstances[gridStateId].AddSpline(spline, 0, ind + 1);
                ind = 0;
                continue;
            }
            if (i == geometryListIn.numberOfCoordinates - 1) 
            {
                spline[ind] = { geometryListIn.xCoordinates[i] , geometryListIn.yCoordinates[i] };
                splinesInstances[gridStateId].AddSpline(spline, 0, ind + 1);
            }
            else 
            {
                spline[ind] = { geometryListIn.xCoordinates[i] , geometryListIn.yCoordinates[i] };
                ind++;
            }
        }

        return 0;
    }

    GRIDGEOM_API int ggeo_curvilinear_mesh_from_splines_ortho(int& gridStateId, 
        GeometryListNative& geometryListIn, 
        CurvilinearParametersNative& curvilinearParameters, 
        SplinesToCurvilinearParametersNative& splineToCurvilinearParameters)
    {
        // use the default constructor, no instance present
        if (splinesInstances.count(gridStateId) == 0)
        {
            GridGeom::Polygons polygon;
            GridGeom::Splines spline(meshInstances[gridStateId].m_projection, polygon);
            splinesInstances[gridStateId] = spline;
        }

        // use the default constructor, no instance present
        if (polygonInstances.count(gridStateId) == 0)
        {
            polygonInstances[gridStateId] = GridGeom::Polygons();
        }

        int success = ggeo_set_splines(gridStateId, geometryListIn);
        GridGeom::CurvilinearGrid curvilinearGrid;
        splinesInstances[gridStateId].SetParameters(curvilinearParameters, splineToCurvilinearParameters);
        splinesInstances[gridStateId].OrthogonalCurvilinearGridFromSplines(curvilinearGrid);
        meshInstances[gridStateId] = GridGeom::Mesh(curvilinearGrid);
        
        return 0;
    }

}