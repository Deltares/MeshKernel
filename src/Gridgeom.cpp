#pragma once

#include "Gridgeom.hpp"
#include "Mesh.hpp"
#include "Orthogonalization.hpp"
#include "OperationsCartesian.cpp"

static std::vector<GridGeom::Mesh> meshInstances;

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
        meshInstances.erase(meshInstances.begin() + gridStateId);
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

        GridGeom::OperationsCartesian operationsCartesian;
        meshInstances[gridStateId].m_operations = &operationsCartesian;
        meshInstances[gridStateId].setMesh(edges, nodes);

        //}
        //else
        //{
        //    auto instance = std::make_unique<GridGeom::Mesh<GridGeom::OperationTypes::sphericalOperations>>();
        //    instance->setMesh(edges, nodes);
        //}

        return 0;
    }

    GRIDGEOM_API int ggeo_get_mesh(int& gridStateId, MeshGeometryDimensions& meshGeometryDimensions, MeshGeometry& meshGeometry)
    {
        
        meshInstances[gridStateId].setState();
                
        meshGeometry.nodex = &meshInstances[gridStateId].m_nodex[0];
        meshGeometry.nodey = &meshInstances[gridStateId].m_nodey[0];
        meshGeometry.nodez = &meshInstances[gridStateId].m_nodez[0];
        meshGeometry.edge_nodes = &meshInstances[gridStateId].m_edgeNodes[0];

        meshGeometryDimensions.numnode = meshInstances[gridStateId].m_nodex.size();
        meshGeometryDimensions.numedge = meshInstances[gridStateId].m_edgeNodes.size() / 2;
        meshGeometryDimensions.numface = meshInstances[gridStateId].m_numFaces;
        meshGeometryDimensions.maxnumfacenodes = 4;

        return 0;
    }

    GRIDGEOM_API int ggeo_orthogonalize(int& gridStateId, int& isTriangulationRequired, int& isAccountingForLandBoundariesRequired, int& projectToLandBoundaryOption,
                           OrthogonalizationParametersNative& orthogonalizationParametersNative, GeometryListNative& geometryListNativePolygon, GeometryListNative& geometryListNativeLandBoundaries)
    {

            GridGeom::Orthogonalization ortogonalization;
            ortogonalization.initialize(meshInstances[gridStateId]);
            ortogonalization.iterate(meshInstances[gridStateId]);
            return 0;


        return -1;
    }
}