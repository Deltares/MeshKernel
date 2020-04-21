#pragma once
#include "MeshGeometryDimensions.hpp"
#include "MeshGeometry.hpp"
#include "GeometryListNative.hpp"
#include "OrthogonalizationParametersNative.hpp"
#include "CurvilinearParametersNative.hpp"
#include "SplinesToCurvilinearParametersNative.hpp"
#include "MakeGridParametersNative.hpp"
#include "InterpolationParametersNative.hpp"
#include "SampleRefineParametersNative.hpp"

#if defined(_WIN32) 
#if !defined(GRIDGEOM_API)
#define GRIDGEOM_API __declspec(dllexport)
#endif
#else  
#define GRIDGEOM_API __attribute__((visibility("default")))
#endif 

// contains all mesh instances
namespace GridGeomApi
{
#ifdef __cplusplus
    extern "C"
    {
#endif
        GRIDGEOM_API int ggeo_new_grid(int& gridStateId);

        GRIDGEOM_API int ggeo_deallocate_state(int& gridStateId);

        GRIDGEOM_API int ggeo_delete_mesh(int& gridStateId, GeometryListNative& geometryListNativePolygon, int& deletionOption);

        GRIDGEOM_API int ggeo_set_state(int& gridStateId, MeshGeometryDimensions& meshGeometryDimensions, MeshGeometry& meshGeometry, bool IsGeographic);

        GRIDGEOM_API int ggeo_orthogonalize(int& gridStateId, int& isTriangulationRequired, int& isAccountingForLandBoundariesRequired, int& projectToLandBoundaryOption,
            OrthogonalizationParametersNative& orthogonalizationParametersNative, GeometryListNative& geometryListNativePolygon, GeometryListNative& geometryListNativeLandBoundaries);

        GRIDGEOM_API int ggeo_get_orthogonality(int& gridStateId, GeometryListNative& geometryList);

        GRIDGEOM_API int ggeo_get_smoothness(int& gridStateId, GeometryListNative& geometryList);

        GRIDGEOM_API int ggeo_get_mesh(int& gridStateId, MeshGeometryDimensions& meshGeometryDimensions, MeshGeometry& meshGeometry);
        
        GRIDGEOM_API int ggeo_orthogonalize_initialize(int& gridStateId,
            int& isTriangulationRequired, 
            int& isAccountingForLandBoundariesRequired, 
            int& projectToLandBoundaryOption,
            OrthogonalizationParametersNative& orthogonalizationParametersNative, 
            GeometryListNative& geometryListNativePolygon, 
            GeometryListNative& geometryListNativeLandBoundaries);

        GRIDGEOM_API int ggeo_orthogonalize_prepare_outer_iteration(int& gridStateId);

        GRIDGEOM_API int ggeo_orthogonalize_inner_iteration(int& gridStateId);

        GRIDGEOM_API int ggeo_orthogonalize_finalize_outer_iteration(int& gridStateId);

        GRIDGEOM_API int ggeo_orthogonalize_delete(int& gridStateId);

        GRIDGEOM_API int ggeo_get_splines(GeometryListNative& geometryListIn, GeometryListNative& geometry_list_out, int& number_of_points_between_vertices);

        GRIDGEOM_API int ggeo_set_splines(int& gridStateId, GeometryListNative& geometryListIn);

        GRIDGEOM_API int ggeo_curvilinear_mesh_from_splines_ortho(int& gridStateId, GeometryListNative& geometryListIn, CurvilinearParametersNative& curvilinearParameters, SplinesToCurvilinearParametersNative& splineToCurvilinearParameters);

        GRIDGEOM_API int ggeo_make_net(int& gridStateId, MakeGridParametersNative& makeGridParameters, GeometryListNative& disposableGeometryListIn);

        GRIDGEOM_API int ggeo_mesh_from_polygon(int& gridStateId, GeometryListNative& disposableGeometryListIn);

        GRIDGEOM_API int ggeo_mesh_from_samples(int& gridStateId, GeometryListNative& disposableGeometryListIn);

        GRIDGEOM_API int ggeo_copy_mesh_boundaries_to_polygon_count_edges(int& gridStateId, int& numberOfPolygonVertices);

        GRIDGEOM_API int ggeo_copy_mesh_boundaries_to_polygon(int& gridStateId, GeometryListNative& disposableGeometryListInOut);

        GRIDGEOM_API int ggeo_refine_polygon_count(int& gridStateId, GeometryListNative& geometryListIn, int& firstIndex, int& secondIndex, double& distance, int& numberOfPolygonVertices);

        GRIDGEOM_API int ggeo_refine_polygon(int& gridStateId, GeometryListNative& geometryListIn, int& firstIndex, int& secondIndex, double& distance, GeometryListNative& geometryListOut);

        GRIDGEOM_API int ggeo_merge_nodes(int& gridStateId, GeometryListNative& geometryListIn);

        GRIDGEOM_API int ggeo_count_vertices_in_polygons(int& gridStateId, GeometryListNative& geometryListIn, int& inside, int& numberOfMeshVertices);

        GRIDGEOM_API int ggeo_vertices_in_polygons(int& gridStateId, GeometryListNative& geometryListIn, int& inside, int& numberOfMeshVertices, int* selectedVertices);
 
        GRIDGEOM_API int ggeo_insert_edge(int& gridStateId, int& start_node, int& end_node, int& new_edge_index);

        GRIDGEOM_API int ggeo_insert_node(int& gridStateId, double& xCoordinate, double& yCoordinate, double& zCoordinate, int& vertexIndex);

        GRIDGEOM_API int ggeo_delete_node(int& gridStateId, int& nodeIndex);

        GRIDGEOM_API int ggeo_offsetted_polygon_count(int& gridStateId, GeometryListNative& geometryListIn, bool& innerAndOuter, double& distance, int& numberOfPolygonVertices);

        GRIDGEOM_API int ggeo_offsetted_polygon(int& gridStateId, GeometryListNative& geometryListIn, bool& innerAndOuter, double& distance, GeometryListNative& geometryListOut);

        GRIDGEOM_API int ggeo_refine_mesh_based_on_samples(int& gridStateId, GeometryListNative& geometryListNative, InterpolationParametersNative& interpolationParametersNative, SampleRefineParametersNative& sampleRefineParametersNative);

#ifdef __cplusplus
    }
#endif
}