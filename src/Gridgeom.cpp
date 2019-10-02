#pragma once

#include "Gridgeom.hpp"
#include "Mesh.cpp"
#include "Orthogonalization.cpp"

namespace GridGeomApi
{
    int ggeo_new_grid(int& gridStateId)
    {
        gridStateId = meshInstances.size() + 1;
        meshInstances.reserve(gridStateId);
        return 0;
    };

    int ggeo_deallocate_state(int& gridStateId)
    {
        meshInstances.erase(meshInstances.begin() + gridStateId);
        return 0;
    }

    int ggeo_set_state(int gridStateId, MeshGeometryDimensions& meshGeometryDimensions, MeshGeometry& meshGeometry, bool IsGeographic)
    {

        std::vector<GridGeom::Edge> edges(meshGeometryDimensions.numedge);
        for (int e = 0; e < edges.size(); e++)
        {
            edges[e].first = 0;
            edges[e].second = 0;
        }

        std::vector<GridGeom::Point> nodes(meshGeometryDimensions.numnode);
        for (int n = 0; n < edges.size(); n++)
        {
            nodes[n].x = meshGeometry.nodex[n];
            nodes[n].y = meshGeometry.nodey[n];
        }

        if (IsGeographic)
        {
            auto instance = std::make_unique<GridGeom::Mesh<GridGeom::OperationTypes::cartesianOperations>>();
            instance->setMesh(edges, nodes);
            meshInstances[gridStateId] = std::move(instance);
        }
        else
        {
            auto instance = std::make_unique<GridGeom::Mesh<GridGeom::OperationTypes::sphericalOperations>>();
            instance->setMesh(edges, nodes);
            meshInstances[gridStateId] = std::move(instance);
        }

        return 0;
    }


    int ggeo_get_mesh(int gridStateId, MeshGeometryDimensions& meshGeometryDimensions, MeshGeometry& meshGeometry)
    {

        auto nodes = meshInstances[gridStateId]->getNodes();
        auto edges = meshInstances[gridStateId]->getEdges();
        int numFaces = meshInstances[gridStateId]->getNumFaces();

        meshGeometryDimensions.numnode = nodes.size();
        meshGeometryDimensions.numedge = nodes.size();
        meshGeometryDimensions.numface = numFaces;

        for (int n = 0; n < nodes.size(); n++)
        {
            meshGeometry.nodex[n] = nodes[n].x;
            meshGeometry.nodey[n] = nodes[n].y;
        }

        int ei = 0;
        for (int e = 0; e < edges.size(); e++)
        {
            meshGeometry.edge_nodes[ei] = edges[e].first;
            ei++;
            meshGeometry.edge_nodes[ei] = edges[e].second;
            ei++;
        }

        return 0;
    }

    int ggeo_orthogonalize(int gridStateId, int isTriangulationRequired, int isAccountingForLandBoundariesRequired, int projectToLandBoundaryOption,
        OrthogonalizationParametersNative& orthogonalizationParametersNative, GeometryListNative& geometryListNativePolygon, GeometryListNative& geometryListNativeLandBoundaries)
    {
        GridGeom::Mesh<GridGeom::cartesianOperations>* cartesianMeshPtr = dynamic_cast<GridGeom::Mesh<GridGeom::cartesianOperations>*>(meshInstances[gridStateId].get());
        if(cartesianMeshPtr !=nullptr)
        {
            GridGeom::Orthogonalization<GridGeom::Mesh<GridGeom::cartesianOperations>> ortogonalization;
            ortogonalization.initialize(*cartesianMeshPtr);
            ortogonalization.iterate(*cartesianMeshPtr);
            return 0;
        }

        GridGeom::Mesh<GridGeom::sphericalOperations>* sphericalMeshPtr = dynamic_cast<GridGeom::Mesh<GridGeom::sphericalOperations>*>(meshInstances[gridStateId].get());
        if (sphericalMeshPtr != nullptr)
        {
            GridGeom::Orthogonalization<GridGeom::Mesh<GridGeom::sphericalOperations>> ortogonalization;
            ortogonalization.initialize(*sphericalMeshPtr);
            ortogonalization.iterate(*sphericalMeshPtr);
            return 0;
        }

        return -1;
    }
}