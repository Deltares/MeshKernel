#pragma once
#include <memory>
#include "Mesh.hpp"
#include "MeshGeometryDimensions.hpp"
#include "MeshGeometry.hpp"
#include "GeometryListNative.hpp"
#include "OrthogonalizationParametersNative.hpp"

#ifdef GRIDGEOM_EXPORTS 
#define GRIDGEOM_API __declspec(dllexport)
#else  
#define GRIDGEOM_API __declspec(dllimport)   
#endif 


namespace GridGeomApi
{
    // contains all mesh instances
    static std::vector<std::unique_ptr<GridGeom::MeshBase>> meshInstances;

#ifdef __cplusplus
    extern "C"
    {
#endif

        GRIDGEOM_API int ggeo_new_grid(int& gridStateId);

        GRIDGEOM_API int ggeo_deallocate_state(int& gridStateId);

        GRIDGEOM_API int ggeo_set_state(int gridStateId, GridGeomApi::MeshGeometryDimensions& meshGeometryDimensions, GridGeomApi::MeshGeometry& meshGeometry, bool IsGeographic);

        GRIDGEOM_API int ggeo_orthogonalize(int gridStateId, int isTriangulationRequired, int isAccountingForLandBoundariesRequired, int projectToLandBoundaryOption,
                               GridGeomApi::OrthogonalizationParametersNative& orthogonalizationParametersNative, GridGeomApi::GeometryListNative& geometryListNativePolygon, GridGeomApi::GeometryListNative& geometryListNativeLandBoundaries);

        GRIDGEOM_API int ggeo_get_mesh(int gridStateId, GridGeomApi::MeshGeometryDimensions& meshGeometryDimensions, GridGeomApi::MeshGeometry& meshGeometry);

#ifdef __cplusplus
    }
#endif
}

