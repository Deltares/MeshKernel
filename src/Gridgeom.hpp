#pragma once
#include <memory>
#include "Mesh.hpp"
#include "MeshGeometryDimensions.hpp"
#include "MeshGeometry.hpp"
#include "GeometryListNative.hpp"
#include "OrthogonalizationParametersNative.hpp"

namespace GridGeomApi
{
    // contains all mesh instances
    static std::vector<std::unique_ptr<GridGeom::MeshBase>> meshInstances;

#ifdef __cplusplus
    extern "C"
    {
#endif

        int ggeo_new_grid(int& gridStateId);

        int ggeo_deallocate_state(int& gridStateId);

        int ggeo_set_state(int gridStateId, MeshGeometryDimensions& meshGeometryDimensions, MeshGeometry& meshGeometry, bool IsGeographic);

        int ggeo_orthogonalize(int gridStateId, int isTriangulationRequired, int isAccountingForLandBoundariesRequired, int projectToLandBoundaryOption,
            OrthogonalizationParametersNative& orthogonalizationParametersNative, GeometryListNative& geometryListNativePolygon, GeometryListNative& geometryListNativeLandBoundaries);

        int ggeo_get_mesh(int gridStateId, MeshGeometryDimensions& meshGeometryDimensions, MeshGeometry& meshGeometry);

#ifdef __cplusplus
    }
#endif
}

