#pragma once
#include "MeshGeometryDimensions.hpp"
#include "MeshGeometry.hpp"
#include "GeometryListNative.hpp"
#include "OrthogonalizationParametersNative.hpp"

#ifdef GRIDGEOM_EXPORTS 
#define GRIDGEOM_API __declspec(dllexport)
#else  
#define GRIDGEOM_API __declspec(dllimport)   
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

        GRIDGEOM_API int ggeo_set_state(int& gridStateId, MeshGeometryDimensions& meshGeometryDimensions, MeshGeometry& meshGeometry, bool IsGeographic);

        GRIDGEOM_API int ggeo_orthogonalize(int& gridStateId, int& isTriangulationRequired, int& isAccountingForLandBoundariesRequired, int& projectToLandBoundaryOption,
            OrthogonalizationParametersNative& orthogonalizationParametersNative, GeometryListNative& geometryListNativePolygon, GeometryListNative& geometryListNativeLandBoundaries);

        GRIDGEOM_API int ggeo_get_mesh(int& gridStateId, MeshGeometryDimensions& meshGeometryDimensions, MeshGeometry& meshGeometry);
        
        GRIDGEOM_API int ggeo_orthogonalize_initialize(int& gridStateId);

        GRIDGEOM_API int ggeo_orthogonalize_prepare_outer_iteration(int& gridStateId);

        GRIDGEOM_API int ggeo_orthogonalize_inner_iteration(int& gridStateId);

        GRIDGEOM_API int ggeo_orthogonalize_finalize_outer_iteration(int& gridStateId);

        GRIDGEOM_API int ggeo_orthogonalize_delete(int& gridStateId);

#ifdef __cplusplus
    }
#endif
}